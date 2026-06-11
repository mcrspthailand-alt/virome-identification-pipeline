#!/usr/bin/env python3
import argparse
import glob
import os
import re
import shlex
import subprocess
from pathlib import Path

import pandas as pd


def load_config(path: str):
    cmd = f"set -a; source {shlex.quote(path)}; env"
    proc = subprocess.run(["bash", "-lc", cmd], check=True, text=True, capture_output=True)
    cfg = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k] = v
    return cfg


def read_tsv(path, **kwargs):
    try:
        return pd.read_csv(path, sep="\t", dtype=str, low_memory=False, **kwargs)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def clean_col(name):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(name).strip()).strip("_").lower()


def first_col(df, names):
    norm = {clean_col(c): c for c in df.columns}
    for name in names:
        key = clean_col(name)
        if key in norm:
            return norm[key]
    return None


def concat_tables(paths, sample_from_path=True, len_bin_from_path=True):
    frames = []
    for path in paths:
        df = read_tsv(path)
        if df.empty:
            continue
        parts = Path(path).parts
        if sample_from_path and "samples" in parts:
            df["sample"] = parts[parts.index("samples") + 1]
        if len_bin_from_path:
            if "07_diamond_ge1000" in parts:
                df["len_bin"] = "ge1000"
            elif "07_diamond_len500_999" in parts:
                df["len_bin"] = "len500_999"
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def collect_tool_annotations(samples_dir, subdir, preferred_patterns, fallback_keywords, prefix):
    rows = []
    for sample_dir in sorted(Path(samples_dir).glob("*")):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        tool_dir = sample_dir / subdir
        paths = []
        for pat in preferred_patterns:
            paths.extend(glob.glob(str(tool_dir / "**" / pat), recursive=True))
        if not paths:
            for path in glob.glob(str(tool_dir / "**" / "*.tsv"), recursive=True):
                low = os.path.basename(path).lower()
                if any(k in low for k in fallback_keywords):
                    paths.append(path)
        for path in sorted(set(paths)):
            df = read_tsv(path)
            if df.empty:
                continue
            contig_col = first_col(df, ["seqname", "seq_name", "contig", "contig_id", "sequence", "sequence_id", "name"])
            if contig_col is None:
                continue
            rel = os.path.relpath(path, tool_dir)
            for _, rec in df.iterrows():
                contig = str(rec.get(contig_col, "")).strip()
                if not contig:
                    continue
                row = {"sample": sample, "contig": contig, f"{prefix}_hit": True, f"{prefix}_source": rel}
                for col in df.columns:
                    if col == contig_col:
                        continue
                    row[f"{prefix}_{clean_col(col)}"] = rec.get(col, "")
                rows.append(row)
    if not rows:
        return pd.DataFrame(columns=["sample", "contig", f"{prefix}_hit", f"{prefix}_source"])

    out = pd.DataFrame(rows)
    for col in out.columns:
        out[col] = out[col].fillna("").astype(str)
    out[f"{prefix}_hit"] = True

    def first_nonempty(values):
        vals = [str(v).strip() for v in values if str(v).strip() and str(v).strip().lower() != "nan"]
        return vals[0] if vals else ""

    agg = {}
    for col in out.columns:
        if col in {"sample", "contig"}:
            continue
        if col == f"{prefix}_hit":
            agg[col] = lambda _: True
        elif col == f"{prefix}_source":
            agg[col] = lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))
        else:
            agg[col] = first_nonempty
    return out.groupby(["sample", "contig"], as_index=False).agg(agg)


def collect_checkv(samples_dir):
    frames = []
    for sample_dir in sorted(Path(samples_dir).glob("*")):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        for path in glob.glob(str(sample_dir / "06_checkv" / "**" / "quality_summary.tsv"), recursive=True):
            df = read_tsv(path)
            if df.empty:
                continue
            contig_col = first_col(df, ["contig_id", "contig", "seq_name", "seqname"])
            if contig_col is None:
                continue
            df = df.rename(columns={contig_col: "contig"})
            df["sample"] = sample
            rename = {}
            for col in df.columns:
                if col not in {"sample", "contig"}:
                    rename[col] = "checkv_" + clean_col(col)
            frames.append(df.rename(columns=rename))
    if not frames:
        return pd.DataFrame(columns=["sample", "contig"])
    return pd.concat(frames, ignore_index=True).drop_duplicates(subset=["sample", "contig"], keep="first")


def collect_prodigal(samples_dir):
    rows = []
    for faa in glob.glob(str(Path(samples_dir) / "*" / "07_diamond_*" / "prodigal" / "*.proteins.faa")):
        parts = Path(faa).parts
        sample = parts[parts.index("samples") + 1] if "samples" in parts else ""
        len_bin = "ge1000" if "07_diamond_ge1000" in parts else "len500_999" if "07_diamond_len500_999" in parts else ""
        counts = {}
        with open(faa, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                q = line[1:].split()[0]
                contig = re.sub(r"_[0-9]+$", "", q)
                counts[contig] = counts.get(contig, 0) + 1
        gff = re.sub(r"\.proteins\.faa$", ".gff", faa)
        for contig, n in counts.items():
            rows.append({
                "sample": sample,
                "len_bin": len_bin,
                "contig": contig,
                "prodigal_total_orfs": n,
                "prodigal_proteins_faa": faa,
                "prodigal_gff": gff if os.path.exists(gff) else "",
            })
    if not rows:
        return pd.DataFrame(columns=["sample", "len_bin", "contig", "prodigal_total_orfs", "prodigal_proteins_faa", "prodigal_gff"])
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)

    out_base = Path(cfg["OUT_BASE"])
    samples_dir = out_base / "samples"
    merged_dir = out_base / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)

    readcounts = read_tsv(merged_dir / "all_contig_read_counts.tsv")
    if readcounts.empty:
        raise SystemExit(f"[ERROR] missing or empty {merged_dir / 'all_contig_read_counts.tsv'}")
    for col in ["sample", "len_bin", "contig"]:
        readcounts[col] = readcounts[col].fillna("").astype(str).str.strip()

    blastx = concat_tables(glob.glob(str(samples_dir / "*" / "07_diamond_*" / "blastx_contigs" / "diamond_blastx_best.tsv")))
    blastp = concat_tables(glob.glob(str(samples_dir / "*" / "07_diamond_*" / "blastp_orfs" / "diamond_blastp_best.tsv")))

    if blastx.empty:
        blastx = pd.DataFrame(columns=["sample", "len_bin", "contig"])
    if blastp.empty:
        blastp = pd.DataFrame(columns=["sample", "len_bin", "contig"])
    for df in [blastx, blastp]:
        for col in ["sample", "len_bin", "contig"]:
            if col not in df.columns:
                df[col] = ""
            df[col] = df[col].fillna("").astype(str).str.strip()

    blastx.to_csv(merged_dir / "diamond_blastx_best_hits.tsv", sep="\t", index=False)
    blastp.to_csv(merged_dir / "diamond_blastp_best_hits.tsv", sep="\t", index=False)

    out = readcounts.merge(blastx, on=["sample", "len_bin", "contig"], how="outer")
    out = out.merge(blastp, on=["sample", "len_bin", "contig"], how="outer")

    vs2 = collect_tool_annotations(samples_dir, "04_virsorter2", ["final-viral-score.tsv", "final_viral_score.tsv", "final-viral-combined.fa.tsv"], ["viral", "score"], "vs2")
    genomad = collect_tool_annotations(samples_dir, "05_genomad", ["*_virus_summary.tsv", "virus_summary.tsv", "*summary.tsv"], ["virus", "summary"], "genomad")
    checkv = collect_checkv(samples_dir)
    prodigal = collect_prodigal(samples_dir)

    out = out.merge(vs2, on=["sample", "contig"], how="left")
    out = out.merge(genomad, on=["sample", "contig"], how="left")
    out = out.merge(checkv, on=["sample", "contig"], how="left")
    out = out.merge(prodigal, on=["sample", "len_bin", "contig"], how="left")

    if "blastx_sseqid" not in out.columns:
        out["blastx_sseqid"] = ""
    if "blastp_sseqid" not in out.columns:
        out["blastp_sseqid"] = ""
    out["blastx_hit"] = out["blastx_sseqid"].fillna("").astype(str).str.strip() != ""
    out["blastp_hit"] = out["blastp_sseqid"].fillna("").astype(str).str.strip() != ""
    for col in ["vs2_hit", "genomad_hit"]:
        if col not in out.columns:
            out[col] = False
        out[col] = out[col].fillna(False).astype(bool)

    evidence_cols = ["vs2_hit", "genomad_hit", "blastx_hit", "blastp_hit"]
    out["viral_evidence_count"] = out[evidence_cols].astype(int).sum(axis=1)
    out["putative_viral"] = out["viral_evidence_count"] >= 1
    out["viral_confidence"] = out["viral_evidence_count"].map(lambda n: "High" if n >= 2 else "Medium" if n == 1 else "None")
    out["viral_evidence"] = out.apply(
        lambda row: ";".join([
            name for name, col in [
                ("VirSorter2", "vs2_hit"),
                ("geNomad", "genomad_hit"),
                ("BLASTx", "blastx_hit"),
                ("BLASTp", "blastp_hit"),
            ] if bool(row.get(col, False))
        ]),
        axis=1,
    )

    out.to_csv(merged_dir / "all_viral_evidence.tsv", sep="\t", index=False)
    print(f"[OK] wrote {merged_dir / 'diamond_blastx_best_hits.tsv'} rows={len(blastx)}")
    print(f"[OK] wrote {merged_dir / 'diamond_blastp_best_hits.tsv'} rows={len(blastp)}")
    print(f"[OK] wrote {merged_dir / 'all_viral_evidence.tsv'} rows={len(out)}")


if __name__ == "__main__":
    main()
