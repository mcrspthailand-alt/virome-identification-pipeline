#!/usr/bin/env python3
import argparse
import re
import shlex
import subprocess
from pathlib import Path

import pandas as pd


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


def read_tsv(path):
    try:
        return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def classify_virus_group(df):
    for col in ["virus_name", "family_or_best", "lineage_ranks", "blastx_stitle", "blastp_stitle"] + RANK_COLS:
        if col not in df.columns:
            df[col] = ""
    text = (
        df[["virus_name", "family_or_best", "lineage_ranks", "blastx_stitle", "blastp_stitle"] + RANK_COLS]
        .fillna("")
        .astype(str)
        .agg(" ; ".join, axis=1)
        .str.lower()
    )
    phage_pat = r"(bacteriophage|cyanophage|phage\b|caudovir|caudoviricetes|caudovirales|myoviridae|siphoviridae|podoviridae|drexlerviridae|schitoviridae|microviridae|autographiviridae|herelleviridae|inoviridae|tectiviridae|leviviridae)"
    euk_pat = r"(eukaryot|herpesvir|poxvir|adenovir|papillomavir|polyomavir|orthomyxo|orthomyxoviridae|retrovir|partitiviridae|phycodnaviridae|baculoviridae|circoviridae|parvoviridae|totiviridae|mitoviridae|endornaviridae|polymycoviridae|iridoviridae|nimaviridae|flaviviridae|hepeviridae|tombusviridae|prasinovirus|mimiviridae|marseilleviridae|pithovirus|pandoravirus)"
    arch_pat = r"(archaeal|archaea|haloarchae|rudiviridae|fuselloviridae|bicaudaviridae|ampullaviridae|clavaviridae|guttaviridae|hafunaviridae|lipothrixviridae|tristromaviridae)"

    df["virus_group"] = "Unidentified"
    df.loc[text.str.contains(euk_pat, regex=True, na=False), "virus_group"] = "Eukaryotic viruses"
    df.loc[text.str.contains(arch_pat, regex=True, na=False), "virus_group"] = "Archaeal viruses"
    df.loc[text.str.contains(phage_pat, regex=True, na=False), "virus_group"] = "Bacteriophage"
    return df


def reorder(df):
    front = [
        "sample", "contig", "len_bin", "putative_viral", "viral_confidence",
        "viral_evidence_count", "viral_evidence", "virus_group",
        "contig_length", "mapped_reads", "unmapped_reads",
        "vs2_hit", "genomad_hit", "blastx_hit", "blastp_hit",
        "virus_name", "taxid", "family_or_best", "best_rank",
        "kingdom", "phylum", "class", "order", "family", "genus", "species", "lineage_ranks",
        "blastx_sseqid", "blastx_pident", "blastx_aln_len", "blastx_evalue", "blastx_bitscore", "blastx_qcovhsp", "blastx_stitle",
        "blastp_sseqid", "blastp_pident", "blastp_aln_len", "blastp_evalue", "blastp_bitscore", "blastp_qcovhsp", "blastp_stitle", "blastp_n_orfs_strict",
        "prodigal_total_orfs",
    ]
    front = [c for c in front if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)

    merged_dir = Path(cfg["OUT_BASE"]) / "merged"
    parsed_path = merged_dir / "all_viral_evidence.parsed_organisms.tsv"
    evidence_path = parsed_path if parsed_path.exists() else merged_dir / "all_viral_evidence.tsv"
    evidence = read_tsv(evidence_path)
    if evidence.empty:
        raise SystemExit(f"[ERROR] missing evidence table: {evidence_path}")

    tax_map = read_tsv(merged_dir / "virus_taxonomy_map.tsv")
    if "virus_name" not in evidence.columns:
        evidence["virus_name"] = ""
    if tax_map.empty:
        tax_map = pd.DataFrame(columns=["virus_name", "taxid", "family_or_best", "best_rank", "lineage_ranks"] + RANK_COLS)

    for df in [evidence, tax_map]:
        df["virus_name"] = df["virus_name"].fillna("").astype(str).str.strip()
    keep_cols = [c for c in ["virus_name", "taxid", "family_or_best", "best_rank", "lineage_ranks"] + RANK_COLS if c in tax_map.columns]
    tax_map = tax_map[keep_cols].drop_duplicates(subset=["virus_name"], keep="first")
    out = evidence.merge(tax_map, on="virus_name", how="left")

    for col in ["putative_viral", "vs2_hit", "genomad_hit", "blastx_hit", "blastp_hit"]:
        if col not in out.columns:
            out[col] = False
        out[col] = out[col].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})

    out = classify_virus_group(out)
    out = reorder(out)
    viral = out[out["putative_viral"]].copy()

    all_annotated = merged_dir / "all_contigs_with_annotations.tsv"
    final_all = merged_dir / "final_virus_contigs.tsv"
    final_phage = merged_dir / "final_virus_contigs_bacteriophage.tsv"
    final_euk = merged_dir / "final_virus_contigs_eukaryotic.tsv"
    final_arch = merged_dir / "final_virus_contigs_archaeal.tsv"
    final_unid = merged_dir / "final_virus_contigs_unidentified.tsv"
    report = merged_dir / "final_virus_summary.txt"

    out.to_csv(all_annotated, sep="\t", index=False)
    viral.to_csv(final_all, sep="\t", index=False)
    viral[viral["virus_group"] == "Bacteriophage"].to_csv(final_phage, sep="\t", index=False)
    viral[viral["virus_group"] == "Eukaryotic viruses"].to_csv(final_euk, sep="\t", index=False)
    viral[viral["virus_group"] == "Archaeal viruses"].to_csv(final_arch, sep="\t", index=False)
    viral[viral["virus_group"] == "Unidentified"].to_csv(final_unid, sep="\t", index=False)

    with open(report, "w", encoding="utf-8") as handle:
        handle.write(f"all_contigs\t{len(out)}\n")
        handle.write(f"putative_viral_contigs\t{len(viral)}\n")
        handle.write("viral_confidence_counts\n")
        handle.write(viral.get("viral_confidence", pd.Series(dtype=str)).value_counts(dropna=False).to_string())
        handle.write("\nvirus_group_counts\n")
        handle.write(viral.get("virus_group", pd.Series(dtype=str)).value_counts(dropna=False).to_string())
        handle.write("\n")

    print(f"[OK] wrote {all_annotated}")
    print(f"[OK] wrote {final_all}")
    print(f"[OK] wrote {final_phage}")
    print(f"[OK] wrote {final_euk}")
    print(f"[OK] wrote {final_arch}")
    print(f"[OK] wrote {final_unid}")
    print(f"[OK] wrote {report}")


if __name__ == "__main__":
    main()
