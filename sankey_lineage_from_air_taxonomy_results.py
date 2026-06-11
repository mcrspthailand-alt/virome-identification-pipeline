#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd

RANK_COLS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

def clean_text(s):
    if pd.isna(s):
        return ""
    return str(s).strip()

def ensure_lineage_cols(df):
    for c in RANK_COLS:
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].fillna("").astype(str).str.strip()
    return df

def lineage_from_row(row, sep=";"):
    parts = []
    for c in RANK_COLS:
        v = clean_text(row.get(c, ""))
        if v and v.lower() not in {"na", "n/a", "none", "unassigned"}:
            parts.append(v)
    return sep.join(parts)

def best_available_lineage(row):
    if "lineage_ranks" in row and clean_text(row["lineage_ranks"]):
        return clean_text(row["lineage_ranks"])
    return lineage_from_row(row)

def load_rank_table(ranks_path):
    rk = read_tsv(ranks_path)
    # support either headered or headerless 8-column file
    if list(rk.columns[:8]) != ["taxid"] + RANK_COLS:
        if rk.shape[1] >= 8:
            rk = rk.iloc[:, :8].copy()
            rk.columns = ["taxid"] + RANK_COLS
    rk["taxid"] = rk["taxid"].fillna("").astype(str).str.strip()
    return ensure_lineage_cols(rk)

def main():
    ap = argparse.ArgumentParser(description="Create lineage-rich and Sankey-ready outputs from air taxonomy results")
    ap.add_argument("--merged-dir", default="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/air/non_hybrid/viral_id_pipeline/merged")
    ap.add_argument("--taxonomy-map", default=None, help="Path to virus_taxonomy_map.tsv")
    ap.add_argument("--ranks-tsv", default=None, help="Path to taxid.ranks.tsv")
    ap.add_argument("--by-sample", default=None, help="Path to tier1_class_by_sample.tsv")
    ap.add_argument("--out-prefix", default="air_taxonomy_lineage")
    args = ap.parse_args()

    merged_dir = args.merged_dir
    map_path = args.taxonomy_map or os.path.join(merged_dir, "virus_taxonomy_map.tsv")
    ranks_path = args.ranks_tsv or os.path.join(merged_dir, "taxid.ranks.tsv")
    by_sample_path = args.by_sample or os.path.join(merged_dir, "tier1_class_by_sample.tsv")

    if not os.path.exists(map_path):
        raise SystemExit(f"[ERROR] missing taxonomy map: {map_path}")
    if not os.path.exists(ranks_path):
        raise SystemExit(f"[ERROR] missing taxid ranks: {ranks_path}")
    if not os.path.exists(by_sample_path):
        raise SystemExit(f"[ERROR] missing by-sample table: {by_sample_path}")

    mp = read_tsv(map_path)
    rk = load_rank_table(ranks_path)
    bs = read_tsv(by_sample_path)

    if "taxid" not in mp.columns:
        raise SystemExit("[ERROR] taxonomy map has no taxid column")
    if "blastx_organism" not in mp.columns:
        raise SystemExit("[ERROR] taxonomy map has no blastx_organism column")

    mp["taxid"] = mp["taxid"].fillna("").astype(str).str.strip()
    mp["blastx_organism"] = mp["blastx_organism"].fillna("").astype(str).str.strip()

    # merge full rank columns into taxonomy map
    map_full = mp.merge(rk, on="taxid", how="left", suffixes=("", "_rk"))
    map_full = ensure_lineage_cols(map_full)
    if "family_or_best" not in map_full.columns:
        map_full["family_or_best"] = ""
    if "best_rank" not in map_full.columns:
        map_full["best_rank"] = ""
    if "lineage_ranks" not in map_full.columns:
        map_full["lineage_ranks"] = ""
    map_full["lineage_full"] = map_full.apply(best_available_lineage, axis=1)
    map_full["lineage_pipe"] = map_full["lineage_full"].str.replace(";", "|", regex=False)

    # write taxonomy map with explicit lineage ranks
    out_map = os.path.join(merged_dir, f"{args.out_prefix}.virus_taxonomy_map.with_lineage.tsv")
    cols = []
    for c in ["blastx_organism", "taxid", "family_or_best", "best_rank", "lineage_ranks"]:
        if c in map_full.columns:
            cols.append(c)
    cols += RANK_COLS + ["lineage_full", "lineage_pipe"]
    cols = [c for i, c in enumerate(cols) if c in map_full.columns and c not in cols[:i]]
    map_full[cols].to_csv(out_map, sep="\t", index=False)

    # enrich by-sample long table
    if "blastx_organism" not in bs.columns:
        raise SystemExit("[ERROR] by-sample table has no blastx_organism column")
    if "sample" not in bs.columns or "reads" not in bs.columns:
        raise SystemExit("[ERROR] by-sample table must contain sample and reads columns")

    bs["blastx_organism"] = bs["blastx_organism"].fillna("").astype(str).str.strip()
    bs["sample"] = bs["sample"].fillna("").astype(str).str.strip()
    bs["reads"] = pd.to_numeric(bs["reads"], errors="coerce").fillna(0)

    map_for_merge = map_full[[c for c in ["blastx_organism", "taxid", "family_or_best", "best_rank"] + RANK_COLS + ["lineage_full", "lineage_pipe"] if c in map_full.columns]].drop_duplicates()
    long_lineage = bs.merge(map_for_merge, on="blastx_organism", how="left")
    for c in ["taxid", "family_or_best", "best_rank"] + RANK_COLS + ["lineage_full", "lineage_pipe"]:
        if c not in long_lineage.columns:
            long_lineage[c] = ""
        long_lineage[c] = long_lineage[c].fillna("").astype(str)
    out_long = os.path.join(merged_dir, f"{args.out_prefix}.by_sample.with_lineage.tsv")
    long_lineage.to_csv(out_long, sep="\t", index=False)

    # Sankey-ready: sample -> each rank
    rank_edges = []
    use_cols = ["kingdom", "phylum", "class", "order", "family"]
    for _, row in long_lineage.iterrows():
        reads = float(row["reads"])
        if reads <= 0:
            continue
        sample = clean_text(row["sample"])
        prev = sample
        for rank in use_cols:
            cur = clean_text(row.get(rank, ""))
            if not cur:
                continue
            rank_edges.append((sample, prev, cur, rank, reads, clean_text(row["blastx_organism"])))
            prev = cur

    if rank_edges:
        edges = pd.DataFrame(rank_edges, columns=["sample", "source", "target", "target_rank", "reads", "blastx_organism"])
        sankey_edges = edges.groupby(["sample", "source", "target", "target_rank"], as_index=False)["reads"].sum()
    else:
        sankey_edges = pd.DataFrame(columns=["sample", "source", "target", "target_rank", "reads"])
    out_edges = os.path.join(merged_dir, f"{args.out_prefix}.sankey_edges.sample_to_family.tsv")
    sankey_edges.to_csv(out_edges, sep="\t", index=False)

    # Sankey-ready full leaf: sample -> lineage_full
    leaf = long_lineage.copy()
    leaf = leaf[leaf["reads"] > 0].copy()
    leaf["lineage_label"] = leaf["lineage_full"].where(leaf["lineage_full"].astype(str).str.strip() != "", "Unassigned")
    leaf_sankey = leaf.groupby(["sample", "lineage_label"], as_index=False)["reads"].sum()
    out_leaf = os.path.join(merged_dir, f"{args.out_prefix}.sankey_sample_to_lineage.tsv")
    leaf_sankey.to_csv(out_leaf, sep="\t", index=False)

    report_path = os.path.join(merged_dir, f"{args.out_prefix}.report.txt")
    n_total = len(map_full)
    n_lineage = int((map_full["lineage_full"].astype(str).str.strip() != "").sum())
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write("Lineage export report\n")
        fh.write(f"taxonomy_map_rows\t{n_total}\n")
        fh.write(f"rows_with_lineage\t{n_lineage}\n")
        fh.write(f"percent_with_lineage\t{(100.0*n_lineage/n_total if n_total else 0):.2f}\n")
        fh.write(f"output_map_with_lineage\t{out_map}\n")
        fh.write(f"output_by_sample_with_lineage\t{out_long}\n")
        fh.write(f"output_sankey_edges\t{out_edges}\n")
        fh.write(f"output_sankey_sample_to_lineage\t{out_leaf}\n")

    print(f"[OK] wrote {out_map}")
    print(f"[OK] wrote {out_long}")
    print(f"[OK] wrote {out_edges}")
    print(f"[OK] wrote {out_leaf}")
    print(f"[OK] wrote {report_path}")

if __name__ == "__main__":
    main()
