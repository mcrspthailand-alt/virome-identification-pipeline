#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Sediment sd01-sd24: rerun DIAMOND with full alignment fields, then merge
# readcounts/length + taxonomy/lineage + virus grouping in ONE script.
#
# What this script does:
#   1) Re-run ONLY DIAMOND blastp from existing prodigal protein FASTA files
#      under viral_id_pipeline/samples for sd01-sd24
#   2) Keep richer alignment fields (qstart/qend/sstart/send/... etc.)
#   3) Rebuild merged master table
#   4) Merge readcounts / contig_length / mapped_reads / unmapped_reads / bam
#   5) Parse blastx_organism from stitle (OS=... / [organism] fallback)
#   6) Merge taxonomy from existing virus_taxonomy_map.tsv
#   7) Split lineage into kingdom/phylum/class/order/family/genus/species
#   8) Classify into:
#        - Bacteriophage
#        - Eukaryotic viruses
#        - Archaeal viruses
#        - Unidentified
#   9) Export:
#        - full merged table with blast+taxonomy+group
#        - phage/euk/archaeal/unidentified subsets
#        - group report
#
# IMPORTANT:
# - Designed for sediment samples sd01-sd24
# - Uses existing outputs from your SED_WF viral_id_pipeline
# - Assumes virus_taxonomy_map.tsv already exists and has lineage_ranks
###############################################################################

PROJECT_ROOT="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/SED_WF"
OUT_BASE="${PROJECT_ROOT}/viral_id_pipeline"
SAMPLES_DIR="${OUT_BASE}/samples"
MERGED_DIR="${OUT_BASE}/merged"
THREADS=32

DIAMOND_DB="/home/panda/workspace/databases/diamond/virus_nr/virus_nr.dmnd"

# sample list: sediment sd01-sd24
SAMPLES=(sd01 sd02 sd03 sd04 sd05 sd06 sd07 sd08 sd09 sd10 sd11 sd12 sd13 sd14 sd15 sd16 sd17 sd18 sd19 sd20 sd21 sd22 sd23 sd24)
LEN_BINS=(ge1000 len500_999)

LOG="${MERGED_DIR}/rerun_fullblast_and_finalize.log"
mkdir -p "${MERGED_DIR}"

echo "=== started: $(date) ===" | tee "${LOG}"
echo "[INFO] PROJECT_ROOT=${PROJECT_ROOT}" | tee -a "${LOG}"
echo "[INFO] SAMPLES_DIR=${SAMPLES_DIR}" | tee -a "${LOG}"
echo "[INFO] MERGED_DIR=${MERGED_DIR}" | tee -a "${LOG}"

[[ -d "${SAMPLES_DIR}" ]] || { echo "[ERROR] missing SAMPLES_DIR: ${SAMPLES_DIR}" | tee -a "${LOG}"; exit 1; }
[[ -s "${DIAMOND_DB}" ]] || { echo "[ERROR] missing DIAMOND_DB: ${DIAMOND_DB}" | tee -a "${LOG}"; exit 1; }
[[ -s "${MERGED_DIR}/ALL_samples.readcounts.by_lenbin.tsv" ]] || { echo "[ERROR] missing readcounts merged file" | tee -a "${LOG}"; exit 1; }
[[ -s "${MERGED_DIR}/virus_taxonomy_map.tsv" ]] || { echo "[ERROR] missing taxonomy map: ${MERGED_DIR}/virus_taxonomy_map.tsv" | tee -a "${LOG}"; exit 1; }

for exe in diamond python3; do
  command -v "${exe}" >/dev/null 2>&1 || { echo "[ERROR] ${exe} not found" | tee -a "${LOG}"; exit 1; }
done

###############################################################################
# 1) rerun diamond with full fields from existing prodigal protein fasta
###############################################################################

OUTFMT_FIELDS=(
  qseqid
  sseqid
  pident
  length
  mismatch
  gapopen
  qstart
  qend
  qlen
  sstart
  send
  slen
  evalue
  bitscore
  qcovhsp
  stitle
)

run_one_bin() {
  local sample="$1"
  local len_bin="$2"
  local d07="${SAMPLES_DIR}/${sample}/07_diamond_option1_all_depth_ge2_${len_bin}"
  local prod="${d07}/prodigal_out"
  local dout="${d07}/diamond_out_full"
  mkdir -p "${dout}"

  local prot_faa="${prod}/${sample}.${len_bin}.depth_ge2.proteins.faa"
  local raw="${dout}/viral_protein_hits.full.tsv"
  local best="${dout}/contig_best_hit.full.tsv"
  local best_filt="${dout}/contig_best_hit.full.filtered.tsv"

  if [[ ! -s "${prot_faa}" ]]; then
    echo "[WARN] ${sample} ${len_bin}: missing proteins faa -> ${prot_faa}" | tee -a "${LOG}"
    return 0
  fi

  echo "[INFO] diamond ${sample} ${len_bin}" | tee -a "${LOG}"
  diamond blastp \
    -d "${DIAMOND_DB}" \
    -q "${prot_faa}" \
    -o "${raw}" \
    -f 6 "${OUTFMT_FIELDS[@]}" \
    --evalue 1e-10 \
    --max-target-seqs 5 \
    --threads "${THREADS}"

  python3 - "${raw}" "${best}" "${best_filt}" "${sample}" "${len_bin}" <<'PY'
import sys
import pandas as pd

raw, best, best_filt, sample, len_bin = sys.argv[1:6]

cols = [
    "qseqid","sseqid","pident","aln_len","mismatch","gapopen",
    "qstart","qend","qlen","sstart","send","slen",
    "evalue","bitscore","qcovhsp","stitle"
]

try:
    df = pd.read_csv(raw, sep="\t", names=cols, dtype=str, low_memory=False)
except pd.errors.EmptyDataError:
    df = pd.DataFrame(columns=cols)

if len(df) == 0:
    out = pd.DataFrame(columns=["contig"] + cols[1:] + ["n_orfs_strict","sample","len_bin"])
    out.to_csv(best, sep="\t", index=False)
    out.to_csv(best_filt, sep="\t", index=False)
    raise SystemExit(0)

for c in ["pident","aln_len","mismatch","gapopen","qstart","qend","qlen","sstart","send","slen","evalue","bitscore","qcovhsp"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df["contig"] = df["qseqid"].astype(str).str.replace(r"_[0-9]+$", "", regex=True)

strict = (
    df[
        (df["evalue"] <= 1e-10) &
        (df["bitscore"] >= 100) &
        (df["aln_len"] >= 60)
    ][["contig","qseqid"]]
    .drop_duplicates()
    .groupby("contig")
    .size()
    .rename("n_orfs_strict")
)

best_df = (
    df.sort_values(["contig","bitscore","evalue"], ascending=[True,False,True])
      .groupby("contig", as_index=False)
      .head(1)
      .copy()
)

best_df = best_df.merge(strict, on="contig", how="left")
best_df["n_orfs_strict"] = best_df["n_orfs_strict"].fillna(0).astype(int)
best_df["sample"] = sample
best_df["len_bin"] = len_bin

for c in ["pident","aln_len","mismatch","gapopen","qstart","qend","qlen","sstart","send","slen","evalue","bitscore","qcovhsp"]:
    best_df[c] = best_df[c].astype(str)

best_df.to_csv(best, sep="\t", index=False)

filt = best_df[
    (pd.to_numeric(best_df["evalue"], errors="coerce") <= 1e-10) &
    (pd.to_numeric(best_df["bitscore"], errors="coerce") >= 100) &
    (pd.to_numeric(best_df["aln_len"], errors="coerce") >= 60)
].copy()

filt.to_csv(best_filt, sep="\t", index=False)
PY
}

for sample in "${SAMPLES[@]}"; do
  for len_bin in "${LEN_BINS[@]}"; do
    run_one_bin "${sample}" "${len_bin}"
  done
done

###############################################################################
# 2) rebuild merged full blast table
###############################################################################

MASTER_FULL="${MERGED_DIR}/ALL_samples.master_option1.combined.fullblast.tsv"

python3 - "${SAMPLES_DIR}" "${MASTER_FULL}" <<'PY'
import os, sys, glob, re
import pandas as pd

samples_dir, out = sys.argv[1:3]
frames = []

def extract_org(s):
    s = "" if pd.isna(s) else str(s).strip()
    if not s:
        return ""
    m = re.search(r"\bOS=([^=]+?)\s+OX=\d+\b", s)
    if m:
        return m.group(1).strip()
    m = re.search(r"\[([^\[\]]+)\]\s*$", s)
    if m:
        return m.group(1).strip()
    return s

for sample_dir in sorted(glob.glob(os.path.join(samples_dir, "*"))):
    sample = os.path.basename(sample_dir)
    for len_bin in ("ge1000", "len500_999"):
        p = os.path.join(
            sample_dir,
            f"07_diamond_option1_all_depth_ge2_{len_bin}",
            "diamond_out_full",
            "contig_best_hit.full.filtered.tsv"
        )
        if os.path.exists(p):
            df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
            if len(df) == 0:
                continue
            if "sample" not in df.columns:
                df["sample"] = sample
            if "len_bin" not in df.columns:
                df["len_bin"] = len_bin
            df["blastx_organism"] = df["stitle"].map(extract_org)
            df["blastx_organism_raw"] = df["stitle"].fillna("").astype(str).str.strip()
            frames.append(df)

if frames:
    outdf = pd.concat(frames, ignore_index=True)
else:
    outdf = pd.DataFrame()

outdf.to_csv(out, sep="\t", index=False)
print("[OK] wrote", out, "rows=", len(outdf))
PY

###############################################################################
# 3) merge readcounts / contig_length / mapped_reads / bam
###############################################################################

READS="${MERGED_DIR}/ALL_samples.readcounts.by_lenbin.tsv"
MASTER_FULL_WITH_READS="${MERGED_DIR}/ALL_samples.master_option1.combined.fullblast.with_reads.tsv"

python3 - "${MASTER_FULL}" "${READS}" "${MASTER_FULL_WITH_READS}" <<'PY'
import sys
import pandas as pd

master, reads, out = sys.argv[1:4]
m = pd.read_csv(master, sep="\t", dtype=str, low_memory=False)
r = pd.read_csv(reads, sep="\t", dtype=str, low_memory=False)

for c in ["sample","contig","len_bin"]:
    m[c] = m[c].astype(str).str.strip()
    r[c] = r[c].astype(str).str.strip()

outdf = m.merge(
    r[["sample","len_bin","contig","contig_length","mapped_reads","unmapped_reads","bam"]],
    on=["sample","len_bin","contig"],
    how="left"
)

if "length" not in outdf.columns:
    outdf["length"] = outdf["contig_length"]
else:
    outdf["length"] = outdf["length"].fillna("").astype(str)
    fill_mask = outdf["length"].eq("")
    outdf.loc[fill_mask, "length"] = outdf.loc[fill_mask, "contig_length"]

outdf.to_csv(out, sep="\t", index=False)
print("[OK] wrote", out, "rows=", len(outdf))
PY

###############################################################################
# 4) merge taxonomy + split lineage + classify virus groups + export
###############################################################################

FINAL_ALL="${MERGED_DIR}/final_annotated_with_taxonomy_fullblast.tsv"
FINAL_PHAGE="${MERGED_DIR}/final_phage_only_fullblast.tsv"
FINAL_EUK="${MERGED_DIR}/final_eukaryotic_only_fullblast.tsv"
FINAL_ARCH="${MERGED_DIR}/final_archaeal_only_fullblast.tsv"
FINAL_UNID="${MERGED_DIR}/final_unidentified_only_fullblast.tsv"
FINAL_REPORT="${MERGED_DIR}/final_group_report_fullblast.txt"

python3 - "${MASTER_FULL_WITH_READS}" "${MERGED_DIR}/virus_taxonomy_map.tsv" "${FINAL_ALL}" "${FINAL_PHAGE}" "${FINAL_EUK}" "${FINAL_ARCH}" "${FINAL_UNID}" "${FINAL_REPORT}" <<'PY'
import sys, re
import pandas as pd

master_path, map_path, out_all, out_phage, out_euk, out_arch, out_unid, out_report = sys.argv[1:9]

m = pd.read_csv(master_path, sep="\t", dtype=str, low_memory=False)
tax = pd.read_csv(map_path, sep="\t", dtype=str, low_memory=False)

def extract_org_series(series):
    s = series.fillna("").astype(str).str.strip()
    out = s.copy()
    os_match = s.str.extract(r"\bOS=([^=]+?)\s+OX=\d+\b", expand=False)
    out = os_match.fillna(out)
    bracket_match = out.str.extract(r"\[([^\[\]]+)\]\s*$", expand=False)
    out = bracket_match.fillna(out)
    return out.str.strip()

def normalize_series(series):
    s = series.fillna("").astype(str).str.strip()
    s = s.str.replace(r"\s+\b(isolate|strain|variant|segment|protein|gene)\b.*$", "", regex=True, flags=re.I)
    s = s.str.replace(r"\s*\(.*?\)\s*$", "", regex=True)
    s = s.str.replace(r"\s+", " ", regex=True).str.strip()
    return s

# Always trust stitle for parsing organism; do not trust previous blastx_organism blindly.
if "stitle" in m.columns:
    m["blastx_organism_raw"] = extract_org_series(m["stitle"])
else:
    m["blastx_organism_raw"] = m.get("blastx_organism", pd.Series([""] * len(m)))

m["blastx_organism"] = normalize_series(m["blastx_organism_raw"])

if "blastx_organism" not in tax.columns:
    raise SystemExit("[ERROR] taxonomy map missing blastx_organism")

tax["blastx_organism"] = normalize_series(tax["blastx_organism"])

# deduplicate taxonomy map by normalized organism name
keep_cols = [c for c in ["blastx_organism","taxid","family_or_best","best_rank","lineage_ranks"] if c in tax.columns]
tax = tax[keep_cols].drop_duplicates(subset=["blastx_organism"], keep="first").copy()

df = m.merge(tax, on="blastx_organism", how="left")

for c in ["taxid","family_or_best","best_rank","lineage_ranks"]:
    if c not in df.columns:
        df[c] = ""

df["lineage_ranks"] = df["lineage_ranks"].fillna("").astype(str)

parts = df["lineage_ranks"].str.split(";", expand=True)
rank_cols = ["kingdom","phylum","class","order","family","genus","species"]
for i, col in enumerate(rank_cols):
    if parts is not None and i in parts.columns:
        df[col] = parts[i].fillna("").astype(str).str.strip()
    else:
        df[col] = ""

df["lineage_full"] = df["lineage_ranks"]

# length fix
if "contig_length" in df.columns:
    df["contig_length"] = df["contig_length"].fillna("").astype(str)
else:
    df["contig_length"] = ""

if "length" not in df.columns:
    df["length"] = df["contig_length"]
else:
    df["length"] = df["length"].fillna("").astype(str)
    fill_mask = df["length"].eq("")
    df.loc[fill_mask, "length"] = df.loc[fill_mask, "contig_length"]

text = (
    df[["blastx_organism","family_or_best","lineage_full","kingdom","phylum","class","order","family","genus","species","stitle"]]
    .fillna("")
    .astype(str)
    .agg(" ; ".join, axis=1)
    .str.lower()
)

# taxonomy-driven grouping
phage_pat = r"(bacteriophage|cyanophage|phage\b|caudovir|caudoviricetes|caudovirales|myoviridae|siphoviridae|podoviridae|peduoviridae|drexlerviridae|schitoviridae|microviridae|autographiviridae|herelleviridae|inoviridae|tectiviridae|leviviridae)"
euk_pat = r"(eukaryot|herpesvir|poxvir|adenovir|papillomavir|polyomavir|orthomyxo|orthomyxoviridae|retrovir|partitiviridae|phycodnaviridae|baculoviridae|circoviridae|parvoviridae|totiviridae|mitoviridae|endornaviridae|polymycoviridae|iridoviridae|nimaviridae|flaviviridae|hepeviridae|tombusviridae|prasinovirus|mimiviridae|marseilleviridae|pithovirus|pandoravirus)"
arch_pat = r"(archaeal|archaea|haloarchae|rudiviridae|fuselloviridae|bicaudaviridae|ampullaviridae|clavaviridae|guttaviridae|hafunaviridae|lipothrixviridae|tristromaviridae)"

is_phage = text.str.contains(phage_pat, regex=True, na=False)
is_arch = text.str.contains(arch_pat, regex=True, na=False)
is_euk = text.str.contains(euk_pat, regex=True, na=False)

df["virus_group"] = "Unidentified"
df.loc[is_euk, "virus_group"] = "Eukaryotic viruses"
df.loc[is_arch, "virus_group"] = "Archaeal viruses"
df.loc[is_phage, "virus_group"] = "Bacteriophage"

# reorder: put important columns first
front = [
    "sample","contig","len_bin","length","contig_length","mapped_reads","unmapped_reads","bam",
    "pident","aln_len","mismatch","gapopen","qstart","qend","qlen","sstart","send","slen","evalue","bitscore","qcovhsp","n_orfs_strict",
    "blastx_organism","blastx_organism_raw","stitle","sseqid",
    "taxid","family_or_best","best_rank","kingdom","phylum","class","order","family","genus","species","lineage_ranks","lineage_full","virus_group"
]
front = [c for c in front if c in df.columns]
rest = [c for c in df.columns if c not in front]
df = df[front + rest]

df.to_csv(out_all, sep="\t", index=False)
df[df["virus_group"] == "Bacteriophage"].to_csv(out_phage, sep="\t", index=False)
df[df["virus_group"] == "Eukaryotic viruses"].to_csv(out_euk, sep="\t", index=False)
df[df["virus_group"] == "Archaeal viruses"].to_csv(out_arch, sep="\t", index=False)
df[df["virus_group"] == "Unidentified"].to_csv(out_unid, sep="\t", index=False)

with open(out_report, "w") as fh:
    fh.write(f"rows\t{len(df)}\n")
    fh.write("group_counts\n")
    fh.write(df["virus_group"].value_counts(dropna=False).to_string())
    fh.write("\n")
    fh.write(f"nonempty_length\t{int((df['length'].fillna('').astype(str).str.strip() != '').sum())}\n")
    fh.write(f"nonempty_contig_length\t{int((df['contig_length'].fillna('').astype(str).str.strip() != '').sum())}\n")
    fh.write(f"nonempty_lineage\t{int((df['lineage_full'].fillna('').astype(str).str.strip() != '').sum())}\n")
    fh.write(f"nonempty_qstart\t{int((df['qstart'].fillna('').astype(str).str.strip() != '').sum()) if 'qstart' in df.columns else 0}\n")
    fh.write(f"nonempty_qend\t{int((df['qend'].fillna('').astype(str).str.strip() != '').sum()) if 'qend' in df.columns else 0}\n")

print("[OK] wrote", out_all)
print("[OK] wrote", out_phage)
print("[OK] wrote", out_euk)
print("[OK] wrote", out_arch)
print("[OK] wrote", out_unid)
print("[OK] wrote", out_report)
PY

echo "=== finished: $(date) ===" | tee -a "${LOG}"
echo "[DONE] outputs:" | tee -a "${LOG}"
echo "  ${MASTER_FULL}" | tee -a "${LOG}"
echo "  ${MASTER_FULL_WITH_READS}" | tee -a "${LOG}"
echo "  ${FINAL_ALL}" | tee -a "${LOG}"
echo "  ${FINAL_PHAGE}" | tee -a "${LOG}"
echo "  ${FINAL_EUK}" | tee -a "${LOG}"
echo "  ${FINAL_ARCH}" | tee -a "${LOG}"
echo "  ${FINAL_UNID}" | tee -a "${LOG}"
echo "  ${FINAL_REPORT}" | tee -a "${LOG}"
