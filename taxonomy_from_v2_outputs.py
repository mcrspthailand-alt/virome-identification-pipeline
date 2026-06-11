#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request
from collections import defaultdict
from functools import lru_cache
from typing import Dict, List, Set, Tuple

import pandas as pd

DEFAULT_PROJECT_ROOT = "/home/panda/workspace/projects/nextflow_metagenomic/eDNA/SED_WF"
DEFAULT_OUT_BASE = os.path.join(DEFAULT_PROJECT_ROOT, "viral_id_pipeline")
DEFAULT_MERGED_DIR = os.path.join(DEFAULT_OUT_BASE, "merged")
DEFAULT_MASTER = os.path.join(DEFAULT_MERGED_DIR, "ALL_samples.master_option1.combined.fixed.with_mapped_reads.tsv")
DEFAULT_READS_MATRIX = os.path.join(DEFAULT_MERGED_DIR, "tier1_reads_matrix.combined.tsv")
DEFAULT_TAX_DIR = os.path.join(DEFAULT_OUT_BASE, "taxonomy")
DEFAULT_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

NAME_CLASS_ALLOW = {
    "scientific name",
    "synonym",
    "equivalent name",
    "genbank synonym",
    "includes",
    "genbank common name",
    "common name",
    "acronym",
    "genbank acronym",
    "anamorph",
    "teleomorph",
    "blast name",
}

RANK_ORDER = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_TO_OUTCOL = {
    "superkingdom": "kingdom",
    "kingdom": "kingdom",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}

def eprint(*args):
    print(*args, file=sys.stderr)

def normalize_name(s: str) -> str:
    s = (s or "").strip()
    if not s:
        return ""
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\s+\b(isolate|strain|variant|segment|protein|gene|partial|putative)\b.*$", "", s, flags=re.I)
    s = re.sub(r"\s*\(.*?\)\s*$", "", s)
    s = re.sub(r"\s+", " ", s).strip(" ;,|")
    return s

def parse_organism_from_stitle(stitle: str) -> str:
    s = (stitle or "").strip()
    if not s:
        return ""

    # UniProt/SwissProt style: ... OS=Organism OX=...
    m = re.search(r"\bOS=([^=]+?)\s+OX=\d+", s)
    if m:
        return normalize_name(m.group(1))

    # Bracketed organism at end
    m = re.search(r"\[([^\[\]]+)\]\s*$", s)
    if m:
        return normalize_name(m.group(1))

    # NCBI-like "protein name, partial [Organism]" already handled above.
    # Try "virus/phage ..." free text fallback.
    # Keep a clean tail after first known metadata token.
    s2 = re.split(r"\b(?:GN=|PE=|SV=|OX=|OS=)\b", s)[0].strip()

    # Remove database accession prefix like sp|ID|ENTRY
    if re.match(r"^[a-z]{2}\|[^|]+\|[^ ]+\s+", s2, flags=re.I):
        s2 = re.sub(r"^[a-z]{2}\|[^|]+\|[^ ]+\s+", "", s2, flags=re.I)

    # Heuristic: look for virus/phage-containing chunk
    m = re.search(r"([A-Z][A-Za-z0-9_.\- ]{0,120}?(?:virus|phage)[A-Za-z0-9_.\- ]*)", s2, flags=re.I)
    if m:
        return normalize_name(m.group(1))

    return normalize_name(s2)

def ensure_taxdump(tax_dir: str, taxdump_url: str) -> Tuple[str, str]:
    os.makedirs(tax_dir, exist_ok=True)
    nodes = os.path.join(tax_dir, "nodes.dmp")
    names = os.path.join(tax_dir, "names.dmp")
    tar_path = os.path.join(tax_dir, "taxdump.tar.gz")

    if os.path.exists(nodes) and os.path.exists(names):
        return nodes, names

    if not os.path.exists(tar_path):
        eprint(f"[INFO] downloading taxdump -> {tar_path}")
        urllib.request.urlretrieve(taxdump_url, tar_path)

    eprint(f"[INFO] extracting {tar_path}")
    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf.getmembers():
            base = os.path.basename(member.name)
            if base in {"nodes.dmp", "names.dmp"}:
                tf.extract(member, tax_dir)

    if not (os.path.exists(nodes) and os.path.exists(names)):
        raise SystemExit("[ERROR] taxdump extraction failed; nodes.dmp or names.dmp missing")
    return nodes, names

def iter_dmp_rows(path: str):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [x.strip() for x in line.split("|")]
            yield parts

def build_name_to_taxid(names_path: str, wanted_names_norm: Set[str]) -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    mapping = {}
    all_hits = defaultdict(set)
    wanted = set(wanted_names_norm)

    for parts in iter_dmp_rows(names_path):
        if len(parts) < 4:
            continue
        taxid, name_txt, _, name_class = parts[:4]
        if name_class not in NAME_CLASS_ALLOW:
            continue
        n = normalize_name(name_txt)
        if not n:
            continue
        nl = n.lower()
        if nl in wanted:
            all_hits[nl].add(taxid)
            # prefer scientific name over others, otherwise first
            if nl not in mapping or name_class == "scientific name":
                mapping[nl] = taxid
    return mapping, all_hits

def load_nodes(nodes_path: str):
    parent = {}
    rank = {}
    for parts in iter_dmp_rows(nodes_path):
        if len(parts) < 3:
            continue
        taxid, parent_taxid, r = parts[:3]
        parent[taxid] = parent_taxid
        rank[taxid] = r
    return parent, rank

@lru_cache(maxsize=None)
def lineage_dict_for_taxid_cached(taxid: str, parent_tuple, rank_tuple):
    parent = dict(parent_tuple)
    rank = dict(rank_tuple)
    out = {v: "" for v in ["kingdom","phylum","class","order","family","genus","species"]}
    cur = taxid
    seen = set()
    while cur and cur not in seen and cur in parent:
        seen.add(cur)
        r = rank.get(cur, "")
        if r in RANK_TO_OUTCOL and not out[RANK_TO_OUTCOL[r]]:
            out[RANK_TO_OUTCOL[r]] = cur
        if cur == parent.get(cur):
            break
        cur = parent.get(cur)
    return out

def build_taxid_to_scientific_name(names_path: str, needed_taxids: Set[str]) -> Dict[str, str]:
    sci = {}
    need = set(needed_taxids)
    for parts in iter_dmp_rows(names_path):
        if len(parts) < 4:
            continue
        taxid, name_txt, _, name_class = parts[:4]
        if taxid in need and name_class == "scientific name" and taxid not in sci:
            sci[taxid] = name_txt.strip()
    return sci

def make_rank_table(taxids: List[str], parent: Dict[str, str], rank: Dict[str, str], sci_name: Dict[str, str]) -> pd.DataFrame:
    rows = []
    parent_tuple = tuple(parent.items())
    rank_tuple = tuple(rank.items())
    for taxid in taxids:
        outtax = {c: "" for c in ["kingdom","phylum","class","order","family","genus","species"]}
        cur = taxid
        seen = set()
        while cur and cur not in seen and cur in parent:
            seen.add(cur)
            r = rank.get(cur, "")
            if r in RANK_TO_OUTCOL:
                col = RANK_TO_OUTCOL[r]
                if not outtax[col]:
                    outtax[col] = sci_name.get(cur, "")
            if cur == parent.get(cur):
                break
            cur = parent.get(cur)
        rows.append({
            "taxid": taxid,
            **outtax,
        })
    return pd.DataFrame(rows)

def pick_best_rank(row):
    for rank in ["family","order","class","phylum","kingdom"]:
        v = str(row.get(rank, "") or "").strip()
        if v and v.lower() not in {"na", "n/a", "none"}:
            return v, rank
    return "Unassigned", "none"

def mk_lineage_ranks(row):
    parts = []
    for k in ["kingdom","phylum","class","order","family","genus","species"]:
        v = str(row.get(k, "") or "").strip()
        if v and v.lower() not in {"na","n/a","none"}:
            parts.append(v)
    return ";".join(parts)

def write_name_parsing_outputs(master_df: pd.DataFrame, tax_dir: str):
    unique_names = sorted({x for x in master_df["blastx_organism"].fillna("").astype(str) if x.strip()})
    with open(os.path.join(tax_dir, "tier1_virus_names.cleaned.txt"), "w", encoding="utf-8") as fh:
        for n in unique_names:
            fh.write(n + "\n")

def enrich_from_existing_outputs(args):
    master_path = args.master_tsv
    reads_matrix_path = args.reads_matrix_tsv
    out_dir = args.out_dir
    tax_dir = args.taxonomy_dir
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tax_dir, exist_ok=True)

    if not os.path.exists(master_path):
        raise SystemExit(f"[ERROR] missing master TSV: {master_path}")
    if not os.path.exists(reads_matrix_path):
        raise SystemExit(f"[ERROR] missing reads matrix TSV: {reads_matrix_path}")

    eprint(f"[INFO] reading {master_path}")
    m = pd.read_csv(master_path, sep="\t", dtype=str, low_memory=False)
    if len(m) == 0:
        raise SystemExit("[ERROR] master TSV is empty")

    if "stitle" not in m.columns:
        raise SystemExit("[ERROR] master TSV missing stitle column")

    m["blastx_organism_raw"] = m["stitle"].fillna("").astype(str).str.strip()
    m["blastx_organism"] = m["blastx_organism_raw"].map(parse_organism_from_stitle)
    m["blastx_organism_norm"] = m["blastx_organism"].map(lambda x: normalize_name(x).lower())

    parsed_master = os.path.join(out_dir, "ALL_samples.master_option1.combined.fixed.with_mapped_reads.parsed.tsv")
    m.to_csv(parsed_master, sep="\t", index=False)
    eprint(f"[OK] wrote {parsed_master}")

    write_name_parsing_outputs(m, tax_dir)

    wanted_norm = {x for x in m["blastx_organism_norm"] if isinstance(x, str) and x.strip()}
    eprint(f"[INFO] unique normalized organism names: {len(wanted_norm)}")

    nodes_path, names_path = ensure_taxdump(tax_dir, args.taxdump_url)
    eprint("[INFO] building name->taxid map from NCBI taxdump")
    exact_map, all_hits = build_name_to_taxid(names_path, wanted_norm)

    name_map_rows = []
    for raw_name, clean_name, norm_name in m[["blastx_organism_raw","blastx_organism","blastx_organism_norm"]].drop_duplicates().itertuples(index=False):
        taxid = exact_map.get(norm_name, "")
        hits = sorted(all_hits.get(norm_name, set()))
        name_map_rows.append({
            "blastx_organism_raw": raw_name,
            "blastx_organism": clean_name,
            "blastx_organism_norm": norm_name,
            "taxid": taxid,
            "n_taxid_hits": len(hits),
            "all_taxid_hits": ";".join(hits),
        })
    with_taxid = pd.DataFrame(name_map_rows)
    with_taxid_path = os.path.join(out_dir, "tier1_virus_names.with_taxid.tsv")
    with_taxid.to_csv(with_taxid_path, sep="\t", index=False)
    eprint(f"[OK] wrote {with_taxid_path}")

    parent, rank = load_nodes(nodes_path)
    mapped_taxids = sorted({x for x in with_taxid["taxid"].fillna("").astype(str) if x.strip() and x != "0"})
    sci = build_taxid_to_scientific_name(names_path, set(mapped_taxids))
    ranks_df = make_rank_table(mapped_taxids, parent, rank, sci)
    ranks_path = os.path.join(out_dir, "taxid.ranks.tsv")
    ranks_df.to_csv(ranks_path, sep="\t", header=False, index=False)
    eprint(f"[OK] wrote {ranks_path}")

    vt = with_taxid[["blastx_organism", "taxid"]].copy()
    map_df = vt.merge(ranks_df, on="taxid", how="left")
    best = map_df.apply(lambda r: pd.Series(pick_best_rank(r)), axis=1)
    map_df["family_or_best"] = best[0]
    map_df["best_rank"] = best[1]
    map_df["lineage_ranks"] = map_df.apply(mk_lineage_ranks, axis=1)
    map_df.loc[map_df["taxid"].fillna("").astype(str).str.strip().isin(["", "0"]), ["family_or_best","best_rank","lineage_ranks"]] = ["Unassigned","none",""]
    virus_map = map_df[["blastx_organism","taxid","family_or_best","best_rank","lineage_ranks"]].drop_duplicates().copy()
    map_path = os.path.join(out_dir, "virus_taxonomy_map.tsv")
    virus_map.to_csv(map_path, sep="\t", index=False)
    eprint(f"[OK] wrote {map_path}")

    reads_mat = pd.read_csv(reads_matrix_path, sep="\t", dtype=str, low_memory=False)
    if len(reads_mat) == 0:
        raise SystemExit("[ERROR] reads matrix is empty")

    # rebuild reads matrix using parsed names from master so keys match
    if "sample" not in m.columns or "mapped_reads" not in m.columns:
        raise SystemExit("[ERROR] parsed master missing sample or mapped_reads")
    m["sample"] = m["sample"].fillna("").astype(str).str.strip()
    m["mapped_reads"] = pd.to_numeric(m["mapped_reads"], errors="coerce").fillna(0)
    g = m[m["blastx_organism"].fillna("").astype(str).str.strip() != ""].groupby(["blastx_organism","sample"], as_index=False)["mapped_reads"].sum()
    rebuilt = g.pivot(index="blastx_organism", columns="sample", values="mapped_reads").fillna(0).astype(int).reset_index()
    rebuilt_reads = os.path.join(out_dir, "tier1_reads_matrix.cleaned.tsv")
    rebuilt.to_csv(rebuilt_reads, sep="\t", index=False)
    eprint(f"[OK] wrote {rebuilt_reads}")

    sample_cols = [c for c in rebuilt.columns if c != "blastx_organism"]
    df = rebuilt.merge(virus_map, on="blastx_organism", how="left")
    for c in ["taxid","family_or_best","best_rank","lineage_ranks"]:
        if c not in df.columns:
            df[c] = ""
    df["taxid"] = df["taxid"].fillna("").astype(str)
    df["family_or_best"] = df["family_or_best"].fillna("Unassigned")
    df["best_rank"] = df["best_rank"].fillna("none")
    df["lineage_ranks"] = df["lineage_ranks"].fillna("")
    meta_cols = ["blastx_organism","taxid","family_or_best","best_rank","lineage_ranks"]
    df = df[meta_cols + sample_cols]

    out_reads = os.path.join(out_dir, "tier1_class_reads_matrix.tsv")
    out_pa = os.path.join(out_dir, "tier1_class_presence_absence.tsv")
    out_long = os.path.join(out_dir, "tier1_class_by_sample.tsv")
    df.to_csv(out_reads, sep="\t", index=False)
    pa = df.copy()
    pa[sample_cols] = (pa[sample_cols] > 0).astype(int)
    pa.to_csv(out_pa, sep="\t", index=False)
    long = df.melt(id_vars=meta_cols, value_vars=sample_cols, var_name="sample", value_name="reads")
    long.to_csv(out_long, sep="\t", index=False)

    report_path = os.path.join(out_dir, "taxonomy_assignment_report.txt")
    mapped = (virus_map["taxid"].fillna("").astype(str).str.strip() != "").sum()
    total = len(virus_map)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(f"unique_organism_names\t{total}\n")
        fh.write(f"mapped_taxid\t{mapped}\n")
        fh.write(f"mapping_rate_percent\t{(mapped * 100.0 / total if total else 0):.2f}\n")
        fh.write("best_rank_counts\n")
        fh.write(virus_map["best_rank"].fillna("none").value_counts().to_string())
        fh.write("\n")
    eprint(f"[OK] wrote {out_reads}")
    eprint(f"[OK] wrote {out_pa}")
    eprint(f"[OK] wrote {out_long}")
    eprint(f"[OK] wrote {report_path}")

def main():
    ap = argparse.ArgumentParser(description="Post-process existing v2 sediment pipeline outputs and assign taxonomy from NCBI taxdump.")
    ap.add_argument("--master-tsv", default=DEFAULT_MASTER, help="Existing merged master TSV from v2 pipeline.")
    ap.add_argument("--reads-matrix-tsv", default=DEFAULT_READS_MATRIX, help="Existing reads matrix TSV from v2 pipeline.")
    ap.add_argument("--out-dir", default=DEFAULT_MERGED_DIR, help="Directory to write taxonomy-enriched outputs.")
    ap.add_argument("--taxonomy-dir", default=DEFAULT_TAX_DIR, help="Directory to store taxdump and taxonomy intermediates.")
    ap.add_argument("--taxdump-url", default=DEFAULT_TAXDUMP_URL, help="NCBI taxdump tarball URL.")
    args = ap.parse_args()
    enrich_from_existing_outputs(args)

if __name__ == "__main__":
    main()
