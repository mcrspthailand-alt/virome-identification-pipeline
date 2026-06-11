#!/usr/bin/env python3
import argparse
import os
import re
import shlex
import subprocess
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd


NAME_CLASS_ALLOW = {
    "scientific name", "synonym", "equivalent name", "genbank synonym",
    "includes", "genbank common name", "common name", "acronym",
    "genbank acronym", "blast name",
}
RANK_TO_COL = {
    "superkingdom": "kingdom",
    "kingdom": "kingdom",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}
RANK_COLS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_config(path: str):
    cmd = f"set -a; source {shlex.quote(path)}; env"
    proc = subprocess.run(["bash", "-lc", cmd], check=True, text=True, capture_output=True)
    cfg = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k] = v
    return cfg


def normalize_name(s):
    s = "" if pd.isna(s) else str(s).strip()
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\s+\b(isolate|strain|variant|segment|protein|gene|partial|putative)\b.*$", "", s, flags=re.I)
    s = re.sub(r"\s*\(.*?\)\s*$", "", s)
    return re.sub(r"\s+", " ", s).strip(" ;,|")


def parse_organism(stitle):
    s = "" if pd.isna(stitle) else str(stitle).strip()
    if not s:
        return ""
    m = re.search(r"\bOS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)", s)
    if m:
        return normalize_name(m.group(1))
    m = re.search(r"\[([^\[\]]+)\]\s*$", s)
    if m:
        return normalize_name(m.group(1))
    s = re.split(r"\b(?:GN=|PE=|SV=|OX=|OS=)\b", s)[0].strip()
    if re.match(r"^[a-z]{2}\|[^|]+\|[^ ]+\s+", s, flags=re.I):
        s = re.sub(r"^[a-z]{2}\|[^|]+\|[^ ]+\s+", "", s, flags=re.I)
    m = re.search(r"([A-Z][A-Za-z0-9_.\- ]{0,160}?(?:virus|phage)[A-Za-z0-9_.\- ]*)", s, flags=re.I)
    if m:
        return normalize_name(m.group(1))
    return normalize_name(s)


def ensure_taxdump(tax_dir: Path, url: str):
    tax_dir.mkdir(parents=True, exist_ok=True)
    nodes = tax_dir / "nodes.dmp"
    names = tax_dir / "names.dmp"
    tar_path = tax_dir / "taxdump.tar.gz"
    if nodes.exists() and names.exists():
        return nodes, names
    if not tar_path.exists():
        print(f"[INFO] downloading NCBI taxdump -> {tar_path}")
        urllib.request.urlretrieve(url, tar_path)
    print(f"[INFO] extracting {tar_path}")
    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf.getmembers():
            if os.path.basename(member.name) in {"nodes.dmp", "names.dmp"}:
                tf.extract(member, tax_dir)
    if not nodes.exists() or not names.exists():
        raise SystemExit("[ERROR] taxdump extraction failed")
    return nodes, names


def iter_dmp(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            yield [x.strip() for x in line.split("|")]


def build_taxonomy(names_path, nodes_path, wanted_norm):
    name_to_taxid = {}
    all_hits = {}
    for parts in iter_dmp(names_path):
        if len(parts) < 4:
            continue
        taxid, name_txt, _, name_class = parts[:4]
        if name_class not in NAME_CLASS_ALLOW:
            continue
        norm = normalize_name(name_txt).lower()
        if norm in wanted_norm:
            all_hits.setdefault(norm, set()).add(taxid)
            if norm not in name_to_taxid or name_class == "scientific name":
                name_to_taxid[norm] = taxid

    parent, rank = {}, {}
    for parts in iter_dmp(nodes_path):
        if len(parts) < 3:
            continue
        taxid, parent_taxid, rank_name = parts[:3]
        parent[taxid] = parent_taxid
        rank[taxid] = rank_name

    needed = set(name_to_taxid.values())
    lineage_taxids = set(needed)
    for taxid in list(needed):
        cur, seen = taxid, set()
        while cur and cur not in seen and cur in parent:
            seen.add(cur)
            lineage_taxids.add(cur)
            if cur == parent.get(cur):
                break
            cur = parent.get(cur)

    scientific = {}
    for parts in iter_dmp(names_path):
        if len(parts) < 4:
            continue
        taxid, name_txt, _, name_class = parts[:4]
        if taxid in lineage_taxids and name_class == "scientific name" and taxid not in scientific:
            scientific[taxid] = name_txt.strip()

    return name_to_taxid, all_hits, parent, rank, scientific


def rank_row_for_taxid(taxid, parent, rank, scientific):
    out = {"taxid": taxid, **{c: "" for c in RANK_COLS}}
    cur, seen = taxid, set()
    while cur and cur not in seen and cur in parent:
        seen.add(cur)
        r = rank.get(cur, "")
        if r in RANK_TO_COL:
            col = RANK_TO_COL[r]
            if not out[col]:
                out[col] = scientific.get(cur, "")
        if cur == parent.get(cur):
            break
        cur = parent.get(cur)
    return out


def best_rank(row):
    for col in ["family", "order", "class", "phylum", "kingdom"]:
        val = str(row.get(col, "") or "").strip()
        if val:
            return val, col
    return "Unassigned", "none"


def lineage(row):
    return ";".join(str(row.get(c, "") or "").strip() for c in RANK_COLS if str(row.get(c, "") or "").strip())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)

    out_base = Path(cfg["OUT_BASE"])
    merged_dir = out_base / "merged"
    tax_dir = Path(cfg.get("TAXONOMY_DIR", str(out_base / "taxonomy")))
    taxdump_url = cfg.get("TAXDUMP_URL", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")

    evidence_path = merged_dir / "all_viral_evidence.tsv"
    evidence = pd.read_csv(evidence_path, sep="\t", dtype=str, low_memory=False)
    if evidence.empty:
        raise SystemExit(f"[ERROR] empty evidence table: {evidence_path}")

    stitle = evidence.get("blastx_stitle", pd.Series([""] * len(evidence))).fillna("").astype(str)
    fallback = evidence.get("blastp_stitle", pd.Series([""] * len(evidence))).fillna("").astype(str)
    stitle = stitle.mask(stitle.str.strip().eq(""), fallback)
    evidence["virus_name"] = stitle.map(parse_organism)
    evidence["virus_name_norm"] = evidence["virus_name"].map(lambda x: normalize_name(x).lower())
    evidence.to_csv(merged_dir / "all_viral_evidence.parsed_organisms.tsv", sep="\t", index=False)

    names = sorted({x for x in evidence["virus_name"].fillna("").astype(str) if x.strip()})
    tax_dir.mkdir(parents=True, exist_ok=True)
    with open(tax_dir / "virus_names_cleaned.txt", "w", encoding="utf-8") as handle:
        for name in names:
            handle.write(name + "\n")

    wanted = {normalize_name(x).lower() for x in names if normalize_name(x)}
    nodes_path, names_path = ensure_taxdump(tax_dir, taxdump_url)
    name_to_taxid, all_hits, parent, rank, scientific = build_taxonomy(names_path, nodes_path, wanted)

    name_rows = []
    for name in names:
        norm = normalize_name(name).lower()
        hits = sorted(all_hits.get(norm, set()))
        name_rows.append({
            "virus_name": name,
            "virus_name_norm": norm,
            "taxid": name_to_taxid.get(norm, ""),
            "n_taxid_hits": len(hits),
            "all_taxid_hits": ";".join(hits),
        })
    names_df = pd.DataFrame(name_rows)
    names_df.to_csv(merged_dir / "virus_names_with_taxid.tsv", sep="\t", index=False)

    mapped_taxids = sorted({x for x in names_df["taxid"].fillna("").astype(str) if x.strip()})
    ranks = pd.DataFrame([rank_row_for_taxid(t, parent, rank, scientific) for t in mapped_taxids])
    if ranks.empty:
        ranks = pd.DataFrame(columns=["taxid"] + RANK_COLS)
    ranks.to_csv(merged_dir / "taxid_ranks.tsv", sep="\t", index=False)

    tax_map = names_df[["virus_name", "taxid"]].merge(ranks, on="taxid", how="left")
    for c in RANK_COLS:
        if c not in tax_map.columns:
            tax_map[c] = ""
        tax_map[c] = tax_map[c].fillna("")
    br = tax_map.apply(lambda row: pd.Series(best_rank(row)), axis=1)
    tax_map["family_or_best"] = br[0]
    tax_map["best_rank"] = br[1]
    tax_map["lineage_ranks"] = tax_map.apply(lineage, axis=1)
    tax_map.loc[tax_map["taxid"].fillna("").astype(str).str.strip().eq(""), ["family_or_best", "best_rank", "lineage_ranks"]] = ["Unassigned", "none", ""]
    tax_map = tax_map[["virus_name", "taxid", "family_or_best", "best_rank", "lineage_ranks"] + RANK_COLS].drop_duplicates()
    tax_map.to_csv(merged_dir / "virus_taxonomy_map.tsv", sep="\t", index=False)

    reads = evidence.copy()
    reads = reads[reads["virus_name"].fillna("").astype(str).str.strip() != ""].copy()
    reads["mapped_reads"] = pd.to_numeric(reads.get("mapped_reads", 0), errors="coerce").fillna(0)
    reads["sample"] = reads["sample"].fillna("").astype(str).str.strip()
    if reads.empty:
        matrix = pd.DataFrame(columns=["virus_name"])
    else:
        g = reads.groupby(["virus_name", "sample"], as_index=False)["mapped_reads"].sum()
        matrix = g.pivot(index="virus_name", columns="sample", values="mapped_reads").fillna(0).astype(int).reset_index()
    matrix.to_csv(merged_dir / "virus_reads_matrix.tsv", sep="\t", index=False)

    sample_cols = [c for c in matrix.columns if c != "virus_name"]
    mat_tax = matrix.merge(tax_map[["virus_name", "taxid", "family_or_best", "best_rank", "lineage_ranks"]], on="virus_name", how="left")
    mat_tax.to_csv(merged_dir / "virus_reads_matrix_with_taxonomy.tsv", sep="\t", index=False)
    pa = mat_tax.copy()
    for c in sample_cols:
        pa[c] = (pd.to_numeric(pa[c], errors="coerce").fillna(0) > 0).astype(int)
    pa.to_csv(merged_dir / "virus_presence_absence_matrix.tsv", sep="\t", index=False)
    if sample_cols:
        long = mat_tax.melt(id_vars=[c for c in ["virus_name", "taxid", "family_or_best", "best_rank", "lineage_ranks"] if c in mat_tax.columns], value_vars=sample_cols, var_name="sample", value_name="reads")
    else:
        long = pd.DataFrame(columns=["virus_name", "taxid", "family_or_best", "best_rank", "lineage_ranks", "sample", "reads"])
    long.to_csv(merged_dir / "virus_reads_by_sample_long.tsv", sep="\t", index=False)

    report = merged_dir / "taxonomy_assignment_report.txt"
    mapped = int((tax_map["taxid"].fillna("").astype(str).str.strip() != "").sum())
    total = len(tax_map)
    with open(report, "w", encoding="utf-8") as handle:
        handle.write(f"unique_virus_names\t{total}\n")
        handle.write(f"mapped_taxid\t{mapped}\n")
        handle.write(f"mapping_rate_percent\t{(100.0 * mapped / total if total else 0):.2f}\n")
        handle.write("best_rank_counts\n")
        handle.write(tax_map["best_rank"].fillna("none").value_counts(dropna=False).to_string())
        handle.write("\n")

    print(f"[OK] wrote {merged_dir / 'all_viral_evidence.parsed_organisms.tsv'}")
    print(f"[OK] wrote {merged_dir / 'virus_taxonomy_map.tsv'}")
    print(f"[OK] wrote {merged_dir / 'virus_reads_matrix.tsv'}")
    print(f"[OK] wrote {merged_dir / 'virus_reads_matrix_with_taxonomy.tsv'}")
    print(f"[OK] wrote {merged_dir / 'virus_presence_absence_matrix.tsv'}")
    print(f"[OK] wrote {merged_dir / 'virus_reads_by_sample_long.tsv'}")


if __name__ == "__main__":
    main()
