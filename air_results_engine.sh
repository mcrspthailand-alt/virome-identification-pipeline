#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Air non-hybrid results engine.
# AIR_STAGE=evidence: run DIAMOND full fields and build the evidence table.
# AIR_STAGE=finalize: merge taxonomy and write final result tables.
#
# What this script does:
#   1) Re-run ONLY DIAMOND blastp from existing prodigal protein FASTA files
#      under viral_id_pipeline/samples for all air non-hybrid samples
#   2) Keep richer alignment fields (qstart/qend/sstart/send/... etc.)
#   3) Rebuild merged master table
#   4) Merge readcounts / contig_length / mapped_reads / unmapped_reads / bam
#   5) Merge evidence from VirSorter2, geNomad, BLASTx, and CheckV
#   6) Parse blastx_organism from stitle (OS=... / [organism] fallback)
#   7) Optionally merge taxonomy from virus_taxonomy_map.tsv if available
#   8) Split lineage into kingdom/phylum/class/order/family/genus/species
#   9) Classify into:
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
# - Designed for air samples air13, air14, air15, air17, air19, air20, air21, air23, air24
# - Uses existing outputs from your air/non_hybrid viral_id_pipeline
# - CheckV is reported in the final table but is NOT used as a retention filter
###############################################################################

PROJECT_ROOT="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/air"
OUT_BASE="${PROJECT_ROOT}/non_hybrid/viral_id_pipeline"
SAMPLES_DIR="${OUT_BASE}/samples"
MERGED_DIR="${OUT_BASE}/merged"
THREADS=32

DIAMOND_DB="/home/panda/workspace/databases/diamond/virus_nr/virus_nr.dmnd"

# sample list: air non-hybrid cohort
SAMPLES=(air13 air14 air15 air17 air19 air20 air21 air23 air24)
LEN_BINS=(ge1000 len500_999)
AIR_STAGE="${AIR_STAGE:-all}"  # all, evidence, or finalize

LOG="${MERGED_DIR}/rerun_fullblast_and_finalize.log"
mkdir -p "${MERGED_DIR}"

MASTER_FULL="${MERGED_DIR}/diamond_best_hits_full.tsv"
MASTER_FULL_WITH_READS="${MERGED_DIR}/all_viral_evidence.tsv"
FINAL_ALL="${MERGED_DIR}/final_virus_contigs.tsv"
FINAL_PHAGE="${MERGED_DIR}/final_virus_contigs_bacteriophage.tsv"
FINAL_EUK="${MERGED_DIR}/final_virus_contigs_eukaryotic.tsv"
FINAL_ARCH="${MERGED_DIR}/final_virus_contigs_archaeal.tsv"
FINAL_UNID="${MERGED_DIR}/final_virus_contigs_unidentified.tsv"
FINAL_REPORT="${MERGED_DIR}/final_virus_summary.txt"

case "$AIR_STAGE" in
  all|evidence|finalize) ;;
  *) echo "[ERROR] AIR_STAGE must be one of: all, evidence, finalize" >&2; exit 1 ;;
esac

echo "=== started: $(date) ===" | tee "${LOG}"
echo "[INFO] PROJECT_ROOT=${PROJECT_ROOT}" | tee -a "${LOG}"
echo "[INFO] SAMPLES_DIR=${SAMPLES_DIR}" | tee -a "${LOG}"
echo "[INFO] MERGED_DIR=${MERGED_DIR}" | tee -a "${LOG}"
echo "[INFO] AIR_STAGE=${AIR_STAGE}" | tee -a "${LOG}"

[[ -d "${SAMPLES_DIR}" ]] || { echo "[ERROR] missing SAMPLES_DIR: ${SAMPLES_DIR}" | tee -a "${LOG}"; exit 1; }
[[ -s "${DIAMOND_DB}" ]] || { echo "[ERROR] missing DIAMOND_DB: ${DIAMOND_DB}" | tee -a "${LOG}"; exit 1; }
[[ -s "${MERGED_DIR}/all_contig_read_counts.tsv" ]] || { echo "[ERROR] missing readcounts merged file" | tee -a "${LOG}"; exit 1; }

for exe in diamond python3; do
  command -v "${exe}" >/dev/null 2>&1 || { echo "[ERROR] ${exe} not found" | tee -a "${LOG}"; exit 1; }
done

###############################################################################
# 1) rerun diamond with full fields from existing prodigal protein fasta
###############################################################################

if [[ "$AIR_STAGE" != "finalize" ]]; then

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
  local d07="${SAMPLES_DIR}/${sample}/07_diamond_${len_bin}"
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

  if [[ -s "${best_filt}" ]]; then
    echo "[INFO] skip diamond ${sample} ${len_bin}: existing ${best_filt}" | tee -a "${LOG}"
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
            f"07_diamond_{len_bin}",
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

READS="${MERGED_DIR}/all_contig_read_counts.tsv"

python3 - "${MASTER_FULL}" "${READS}" "${SAMPLES_DIR}" "${MASTER_FULL_WITH_READS}" <<'PY'
import os, sys, glob, re
import pandas as pd

master, reads, samples_dir, out = sys.argv[1:5]
m = pd.read_csv(master, sep="\t", dtype=str, low_memory=False)
r = pd.read_csv(reads, sep="\t", dtype=str, low_memory=False)

for c in ["sample", "contig", "len_bin"]:
    if c not in m.columns:
        m[c] = ""
    m[c] = m[c].fillna("").astype(str).str.strip()
    r[c] = r[c].astype(str).str.strip()

base = r[["sample", "len_bin", "contig", "contig_length", "mapped_reads", "unmapped_reads", "bam"]].copy()
blast = m.copy()
blast["blastx_hit"] = blast["contig"].astype(str).str.strip() != ""
outdf = base.merge(blast, on=["sample", "len_bin", "contig"], how="outer", suffixes=("", "_blast"))

for c in ["contig_length", "mapped_reads", "unmapped_reads", "bam"]:
    bc = f"{c}_blast"
    if bc in outdf.columns:
        outdf[c] = outdf[c].fillna(outdf[bc])
        outdf = outdf.drop(columns=[bc])

if "length" not in outdf.columns:
    outdf["length"] = outdf["contig_length"]
else:
    outdf["length"] = outdf["length"].fillna("").astype(str)
    fill_mask = outdf["length"].eq("")
    outdf.loc[fill_mask, "length"] = outdf.loc[fill_mask, "contig_length"]

def read_table(path):
    try:
        return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    except Exception:
        return pd.DataFrame()

def first_col(df, names):
    norm = {c.lower().replace(" ", "_").replace("-", "_"): c for c in df.columns}
    for n in names:
        key = n.lower().replace(" ", "_").replace("-", "_")
        if key in norm:
            return norm[key]
    return None

def clean_col(name):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(name).strip()).strip("_").lower()

def collect_tool_annotations(sample, tool_dir, preferred_patterns, fallback_keywords, prefix):
    rows = []
    paths = []
    for pat in preferred_patterns:
        paths.extend(glob.glob(os.path.join(tool_dir, "**", pat), recursive=True))
    if not paths:
        for p in glob.glob(os.path.join(tool_dir, "**", "*.tsv"), recursive=True):
            low = os.path.basename(p).lower()
            if any(k in low for k in fallback_keywords):
                paths.append(p)
    for p in sorted(set(paths)):
        df = read_table(p)
        if df.empty:
            continue
        c = first_col(df, ["seqname", "seq_name", "contig", "contig_id", "sequence", "sequence_id", "name"])
        if c is None:
            continue
        df = df.copy()
        df[c] = df[c].fillna("").astype(str).str.strip()
        source = os.path.relpath(p, tool_dir)
        for _, rec in df[df[c] != ""].iterrows():
            row = {"sample": sample, "contig": rec[c], f"{prefix}_source": source, f"{prefix}_hit": True}
            for col in df.columns:
                if col == c:
                    continue
                row[f"{prefix}_{clean_col(col)}"] = rec.get(col, "")
            rows.append(row)
    if not rows:
        return pd.DataFrame(columns=["sample", "contig", f"{prefix}_hit", f"{prefix}_source"])

    ann = pd.DataFrame(rows)
    for col in ann.columns:
        ann[col] = ann[col].fillna("").astype(str)
    ann[f"{prefix}_hit"] = True

    def first_nonempty(values):
        vals = [str(v).strip() for v in values if str(v).strip() and str(v).strip().lower() != "nan"]
        return vals[0] if vals else ""

    agg = {}
    for col in ann.columns:
        if col in {"sample", "contig"}:
            continue
        if col == f"{prefix}_hit":
            agg[col] = lambda x: True
        elif col == f"{prefix}_source":
            agg[col] = lambda x: ";".join(sorted(set(str(v).strip() for v in x if str(v).strip())))
        else:
            agg[col] = first_nonempty
    return ann.groupby(["sample", "contig"], as_index=False).agg(agg)

def count_prodigal_orfs(samples_dir):
    rows = []
    for faa in glob.glob(os.path.join(samples_dir, "*", "07_diamond_*", "prodigal_out", "*.proteins.faa"), recursive=True):
        sample = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(faa))))
        len_bin = "len500_999" if "len500_999" in faa else "ge1000" if "ge1000" in faa else ""
        counts = {}
        try:
            with open(faa, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    q = line[1:].split()[0]
                    contig = re.sub(r"_[0-9]+$", "", q)
                    counts[contig] = counts.get(contig, 0) + 1
        except OSError:
            continue
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

vs2_ann_frames = []
genomad_ann_frames = []
checkv_frames = []

for sample_dir in sorted(glob.glob(os.path.join(samples_dir, "*"))):
    sample = os.path.basename(sample_dir)
    vs2_ann_frames.append(collect_tool_annotations(
        sample,
        os.path.join(sample_dir, "04_virsorter2"),
        ["final-viral-score.tsv", "final_viral_score.tsv", "final-viral-combined.fa.tsv"],
        ["viral", "score"],
        "vs2",
    ))
    genomad_ann_frames.append(collect_tool_annotations(
        sample,
        os.path.join(sample_dir, "05_genomad"),
        ["*_virus_summary.tsv", "virus_summary.tsv", "*summary.tsv"],
        ["virus", "summary"],
        "genomad",
    ))

    for q in glob.glob(os.path.join(sample_dir, "06_checkv", "**", "quality_summary.tsv"), recursive=True):
        qdf = read_table(q)
        if qdf.empty:
            continue
        qdf["sample"] = sample
        checkv_frames.append(qdf)

def rows_to_df(rows, flag_name, source_name):
    if not rows:
        return pd.DataFrame(columns=["sample", "contig", flag_name, source_name])
    df = pd.DataFrame(rows, columns=["sample", "contig", source_name])
    df[flag_name] = True
    return (
        df.groupby(["sample", "contig"], as_index=False)
          .agg({flag_name: "max", source_name: lambda x: ";".join(sorted(set(map(str, x))))})
    )

vs2 = pd.concat([x for x in vs2_ann_frames if not x.empty], ignore_index=True) if any(not x.empty for x in vs2_ann_frames) else pd.DataFrame(columns=["sample", "contig", "vs2_hit", "vs2_source"])
gnm = pd.concat([x for x in genomad_ann_frames if not x.empty], ignore_index=True) if any(not x.empty for x in genomad_ann_frames) else pd.DataFrame(columns=["sample", "contig", "genomad_hit", "genomad_source"])
prod = count_prodigal_orfs(samples_dir)

if not vs2.empty:
    vs2 = vs2.drop_duplicates(subset=["sample", "contig"], keep="first")
if not gnm.empty:
    gnm = gnm.drop_duplicates(subset=["sample", "contig"], keep="first")

outdf = outdf.merge(vs2, on=["sample", "contig"], how="left")
outdf = outdf.merge(gnm, on=["sample", "contig"], how="left")
outdf = outdf.merge(prod, on=["sample", "len_bin", "contig"], how="left")

for old, new in [("virsorter2_hit", "vs2_hit"), ("virsorter2_source", "vs2_source")]:
    if old in outdf.columns and new not in outdf.columns:
        outdf[new] = outdf[old]

for c in ["vs2_hit", "genomad_hit", "blastx_hit"]:
    if c not in outdf.columns:
        outdf[c] = False
    outdf[c] = outdf[c].fillna(False).astype(bool)
for c in ["vs2_source", "genomad_source", "prodigal_proteins_faa", "prodigal_gff"]:
    if c not in outdf.columns:
        outdf[c] = ""
    outdf[c] = outdf[c].fillna("").astype(str)
if "prodigal_total_orfs" not in outdf.columns:
    outdf["prodigal_total_orfs"] = 0
outdf["prodigal_total_orfs"] = pd.to_numeric(outdf["prodigal_total_orfs"], errors="coerce").fillna(0).astype(int)

if checkv_frames:
    cv = pd.concat(checkv_frames, ignore_index=True)
    contig_col = first_col(cv, ["contig_id", "contig", "seq_name", "seqname"])
    if contig_col:
        cv = cv.rename(columns={contig_col: "contig"})
        cv["sample"] = cv["sample"].fillna("").astype(str).str.strip()
        cv["contig"] = cv["contig"].fillna("").astype(str).str.strip()
        rename = {}
        for c in cv.columns:
            if c not in {"sample", "contig"}:
                rename[c] = "checkv_" + re.sub(r"[^A-Za-z0-9_]+", "_", c.strip()).strip("_").lower()
        cv = cv.rename(columns=rename).drop_duplicates(subset=["sample", "contig"], keep="first")
        outdf = outdf.merge(cv, on=["sample", "contig"], how="left")

outdf["putative_viral"] = outdf[["vs2_hit", "genomad_hit", "blastx_hit"]].any(axis=1)
outdf["viral_evidence"] = outdf.apply(
    lambda row: ";".join([
        name for name, flag in [
            ("VirSorter2", row.get("vs2_hit", False)),
            ("geNomad", row.get("genomad_hit", False)),
            ("BLASTx", row.get("blastx_hit", False)),
        ] if bool(flag)
    ]),
    axis=1,
)

outdf.to_csv(out, sep="\t", index=False)
print("[OK] wrote", out, "rows=", len(outdf))
PY

fi

###############################################################################
# 4) merge taxonomy + split lineage + classify virus groups + export
###############################################################################

if [[ "$AIR_STAGE" != "evidence" ]]; then

[[ -s "${MASTER_FULL_WITH_READS}" ]] || { echo "[ERROR] missing evidence master table: ${MASTER_FULL_WITH_READS}" | tee -a "${LOG}"; exit 1; }

if [[ -s "${MERGED_DIR}/virus_taxonomy_map.tsv" ]]; then
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
else
  echo "[WARN] taxonomy map not found yet: ${MERGED_DIR}/virus_taxonomy_map.tsv" | tee -a "${LOG}"
  echo "[WARN] skipped taxonomy/group finalization; run taxonomy_from_air_outputs.py, then rerun this script or final merge." | tee -a "${LOG}"
fi

fi

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
