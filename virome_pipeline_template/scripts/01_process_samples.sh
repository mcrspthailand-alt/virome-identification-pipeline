#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config.env}"
[[ -s "$CONFIG" ]] || { echo "[ERROR] missing config: $CONFIG" >&2; exit 1; }
CONFIG_DIR="$(cd "$(dirname "$CONFIG")" && pwd)"
# shellcheck source=/dev/null
source "$CONFIG"

MAP_THREADS="${MAP_THREADS:-$THREADS}"
SORT_THREADS="${SORT_THREADS:-8}"
SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-1G}"

if [[ "$SAMPLES_TSV" != /* ]]; then
  SAMPLES_TSV="${CONFIG_DIR}/${SAMPLES_TSV}"
fi

OUT_BASE="${OUT_BASE%/}"
SAMPLES_DIR="${OUT_BASE}/samples"
MERGED_DIR="${OUT_BASE}/merged"
LOG_DIR="${OUT_BASE}/logs"
mkdir -p "$SAMPLES_DIR" "$MERGED_DIR" "$LOG_DIR"

need() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] missing command: $1" >&2; exit 1; }
}

for exe in awk gzip samtools minimap2 mosdepth prodigal diamond seqkit bbduk.sh python3 virsorter genomad checkv; do
  need "$exe"
done

[[ -s "$SAMPLES_TSV" ]] || { echo "[ERROR] missing sample sheet: $SAMPLES_TSV" >&2; exit 1; }
[[ -s "$DIAMOND_DB" ]] || { echo "[ERROR] missing DIAMOND_DB: $DIAMOND_DB" >&2; exit 1; }
[[ -d "$GENOMAD_DB" ]] || { echo "[ERROR] missing GENOMAD_DB: $GENOMAD_DB" >&2; exit 1; }
[[ -d "$VS2_DB_DIR" ]] || { echo "[ERROR] missing VS2_DB_DIR: $VS2_DB_DIR" >&2; exit 1; }
[[ -d "$CHECKV_DB" ]] || { echo "[ERROR] missing CHECKV_DB: $CHECKV_DB" >&2; exit 1; }

STATUS_TSV="${LOG_DIR}/sample_status.tsv"
echo -e "sample\tstatus\tmessage" > "$STATUS_TSV"

prepare_contigs() {
  local assembly="$1"
  local out_fa="$2"

  if [[ "$assembly" == *.gz ]]; then
    gzip -dc "$assembly" > "$out_fa"
  else
    cp -f "$assembly" "$out_fa"
  fi
}

write_empty_blastx() {
  local out="$1"
  echo -e "contig\tblastx_sseqid\tblastx_slen\tblastx_pident\tblastx_aln_len\tblastx_mismatch\tblastx_gapopen\tblastx_qstart\tblastx_qend\tblastx_qlen\tblastx_sstart\tblastx_send\tblastx_evalue\tblastx_bitscore\tblastx_qcovhsp\tblastx_stitle" > "$out"
}

write_empty_blastp() {
  local out="$1"
  echo -e "contig\tblastp_sseqid\tblastp_slen\tblastp_pident\tblastp_aln_len\tblastp_mismatch\tblastp_gapopen\tblastp_qstart\tblastp_qend\tblastp_qlen\tblastp_sstart\tblastp_send\tblastp_evalue\tblastp_bitscore\tblastp_qcovhsp\tblastp_stitle\tblastp_n_orfs_strict" > "$out"
}

log_msg() {
  echo "[$(date)] $*"
}

best_blastx_table() {
  local raw="$1"
  local best="$2"
  python3 - "$raw" "$best" "$DIAMOND_EVALUE" "$DIAMOND_BITSCORE" "$DIAMOND_ALNLEN" <<'PY'
import sys
import pandas as pd

raw, best, evalue, bitscore, alnlen = sys.argv[1:6]
evalue = float(evalue)
bitscore = float(bitscore)
alnlen = float(alnlen)

cols = [
    "contig", "blastx_qlen", "blastx_sseqid", "blastx_slen", "blastx_pident",
    "blastx_aln_len", "blastx_mismatch", "blastx_gapopen", "blastx_qstart",
    "blastx_qend", "blastx_sstart", "blastx_send", "blastx_evalue",
    "blastx_bitscore", "blastx_qcovhsp", "blastx_stitle",
]
out_cols = [
    "contig", "blastx_sseqid", "blastx_slen", "blastx_pident", "blastx_aln_len",
    "blastx_mismatch", "blastx_gapopen", "blastx_qstart", "blastx_qend",
    "blastx_qlen", "blastx_sstart", "blastx_send", "blastx_evalue",
    "blastx_bitscore", "blastx_qcovhsp", "blastx_stitle",
]
try:
    df = pd.read_csv(raw, sep="\t", names=cols, dtype=str, low_memory=False)
except pd.errors.EmptyDataError:
    pd.DataFrame(columns=out_cols).to_csv(best, sep="\t", index=False)
    raise SystemExit(0)

if df.empty:
    pd.DataFrame(columns=out_cols).to_csv(best, sep="\t", index=False)
    raise SystemExit(0)

for c in ["blastx_pident", "blastx_aln_len", "blastx_evalue", "blastx_bitscore", "blastx_qcovhsp"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df[(df["blastx_evalue"] <= evalue) & (df["blastx_bitscore"] >= bitscore) & (df["blastx_aln_len"] >= alnlen)].copy()
if df.empty:
    pd.DataFrame(columns=out_cols).to_csv(best, sep="\t", index=False)
    raise SystemExit(0)
best_df = (
    df.sort_values(["contig", "blastx_bitscore", "blastx_evalue"], ascending=[True, False, True])
      .groupby("contig", as_index=False)
      .head(1)
)
best_df[out_cols].to_csv(best, sep="\t", index=False)
PY
}

best_blastp_table() {
  local raw="$1"
  local best="$2"
  python3 - "$raw" "$best" "$DIAMOND_EVALUE" "$DIAMOND_BITSCORE" "$DIAMOND_ALNLEN" <<'PY'
import re
import sys
import pandas as pd

raw, best, evalue, bitscore, alnlen = sys.argv[1:6]
evalue = float(evalue)
bitscore = float(bitscore)
alnlen = float(alnlen)

cols = [
    "qseqid", "blastp_sseqid", "blastp_pident", "blastp_aln_len",
    "blastp_mismatch", "blastp_gapopen", "blastp_qstart", "blastp_qend",
    "blastp_qlen", "blastp_sstart", "blastp_send", "blastp_slen",
    "blastp_evalue", "blastp_bitscore", "blastp_qcovhsp", "blastp_stitle",
]
out_cols = [
    "contig", "blastp_sseqid", "blastp_slen", "blastp_pident", "blastp_aln_len",
    "blastp_mismatch", "blastp_gapopen", "blastp_qstart", "blastp_qend",
    "blastp_qlen", "blastp_sstart", "blastp_send", "blastp_evalue",
    "blastp_bitscore", "blastp_qcovhsp", "blastp_stitle", "blastp_n_orfs_strict",
]
try:
    df = pd.read_csv(raw, sep="\t", names=cols, dtype=str, low_memory=False)
except pd.errors.EmptyDataError:
    pd.DataFrame(columns=out_cols).to_csv(best, sep="\t", index=False)
    raise SystemExit(0)

if df.empty:
    pd.DataFrame(columns=out_cols).to_csv(best, sep="\t", index=False)
    raise SystemExit(0)

df["contig"] = df["qseqid"].astype(str).map(lambda x: re.sub(r"_[0-9]+$", "", x))
for c in ["blastp_pident", "blastp_aln_len", "blastp_evalue", "blastp_bitscore", "blastp_qcovhsp"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

strict = (
    df[(df["blastp_evalue"] <= evalue) & (df["blastp_bitscore"] >= bitscore) & (df["blastp_aln_len"] >= alnlen)]
    [["contig", "qseqid"]]
    .drop_duplicates()
    .groupby("contig")
    .size()
    .rename("blastp_n_orfs_strict")
)
df = df[(df["blastp_evalue"] <= evalue) & (df["blastp_bitscore"] >= bitscore) & (df["blastp_aln_len"] >= alnlen)].copy()
if df.empty:
    pd.DataFrame(columns=out_cols).to_csv(best, sep="\t", index=False)
    raise SystemExit(0)
best_df = (
    df.sort_values(["contig", "blastp_bitscore", "blastp_evalue"], ascending=[True, False, True])
      .groupby("contig", as_index=False)
      .head(1)
)
best_df = best_df.merge(strict, on="contig", how="left")
best_df["blastp_n_orfs_strict"] = best_df["blastp_n_orfs_strict"].fillna(0).astype(int)
best_df[out_cols].to_csv(best, sep="\t", index=False)
PY
}

run_diamond_bin() {
  local sample="$1"
  local len_bin="$2"
  local query_fa="$3"
  local out_dir="$4"

  local blastx_dir="${out_dir}/blastx_contigs"
  local prodigal_dir="${out_dir}/prodigal"
  local blastp_dir="${out_dir}/blastp_orfs"
  mkdir -p "$blastx_dir" "$prodigal_dir" "$blastp_dir"

  local blastx_raw="${blastx_dir}/diamond_blastx_raw.tsv"
  local blastx_best="${blastx_dir}/diamond_blastx_best.tsv"
  if [[ -s "$query_fa" ]]; then
    if [[ ! -s "$blastx_best" ]]; then
      diamond blastx \
        -d "$DIAMOND_DB" \
        -q "$query_fa" \
        -o "$blastx_raw" \
        -f 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle \
        --evalue "$DIAMOND_EVALUE" \
        --max-target-seqs "$DIAMOND_MAX_TARGET_SEQS" \
        --threads "$THREADS"
      best_blastx_table "$blastx_raw" "$blastx_best"
    fi
  else
    write_empty_blastx "$blastx_best"
  fi

  local prot_faa="${prodigal_dir}/${sample}.${len_bin}.proteins.faa"
  local gff="${prodigal_dir}/${sample}.${len_bin}.gff"
  local blastp_raw="${blastp_dir}/diamond_blastp_raw.tsv"
  local blastp_best="${blastp_dir}/diamond_blastp_best.tsv"
  if [[ -s "$query_fa" ]]; then
    [[ -s "$prot_faa" ]] || prodigal -i "$query_fa" -a "$prot_faa" -o "$gff" -p meta
    if [[ -s "$prot_faa" && ! -s "$blastp_best" ]]; then
      diamond blastp \
        -d "$DIAMOND_DB" \
        -q "$prot_faa" \
        -o "$blastp_raw" \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp stitle \
        --evalue "$DIAMOND_EVALUE" \
        --max-target-seqs "$DIAMOND_MAX_TARGET_SEQS" \
        --threads "$THREADS"
      best_blastp_table "$blastp_raw" "$blastp_best"
    elif [[ ! -s "$prot_faa" ]]; then
      write_empty_blastp "$blastp_best"
    fi
  else
    write_empty_blastp "$blastp_best"
  fi
}

run_one_sample() {
  local sample="$1"
  local assembly="$2"
  local r1="$3"
  local r2="$4"

  [[ -s "$assembly" ]] || { echo -e "${sample}\tFAILED\tmissing assembly: ${assembly}" >> "$STATUS_TSV"; return 1; }
  [[ -s "$r1" ]] || { echo -e "${sample}\tFAILED\tmissing R1: ${r1}" >> "$STATUS_TSV"; return 1; }
  [[ -s "$r2" ]] || { echo -e "${sample}\tFAILED\tmissing R2: ${r2}" >> "$STATUS_TSV"; return 1; }

  local sample_dir="${SAMPLES_DIR}/${sample}"
  local d00="${sample_dir}/00_input"
  local d01="${sample_dir}/01_contigs"
  local d02="${sample_dir}/02_mapping"
  local d03="${sample_dir}/03_depth_filtered_contigs"
  local d04="${sample_dir}/04_virsorter2"
  local d05="${sample_dir}/05_genomad"
  local d06="${sample_dir}/06_checkv"
  local d07x="${sample_dir}/07_diamond_ge1000"
  local d07s="${sample_dir}/07_diamond_len500_999"
  local d08="${sample_dir}/08_coverage"
  local log="${LOG_DIR}/${sample}.process.log"
  mkdir -p "$d00" "$d01" "$d02" "$d03" "$d04" "$d05" "$d06" "$d07x" "$d07s" "$d08"

  {
    echo "[INFO] sample=$sample"
    echo "[INFO] assembly=$assembly"
    echo "[INFO] r1=$r1"
    echo "[INFO] r2=$r2"

    ln -sfn "$assembly" "${d00}/${sample}.assembly"
    ln -sfn "$r1" "${d00}/${sample}_R1.fastq.gz"
    ln -sfn "$r2" "${d00}/${sample}_R2.fastq.gz"

    local raw_contigs="${d01}/contigs.raw.fa"
    local filtered="${d01}/contigs.filtered.ge${MIN_CONTIG_LEN}.fa"
    local ge1000="${d01}/contigs.filtered.ge1000.fa"
    local len500_999="${d01}/contigs.filtered.len500_999.fa"
    [[ -s "$raw_contigs" ]] || prepare_contigs "$assembly" "$raw_contigs"
    [[ -s "$filtered" ]] || bbduk.sh in="$raw_contigs" out="$filtered" entropy=0.70 entropywindow=50 entropyk=5 overwrite=t
    [[ -s "$ge1000" ]] || seqkit seq -m 1000 "$filtered" > "$ge1000"
    [[ -s "$len500_999" ]] || seqkit seq -m 500 "$filtered" | seqkit seq -M 999 > "$len500_999"
    seqkit stats "$filtered" "$ge1000" "$len500_999" > "${d01}/seqkit_stats.tsv"

    local bam="${d02}/${sample}.sorted.bam"
    local bam_tmp="${d02}/${sample}.sorted.tmp.bam"
    local sort_tmp_dir="${d02}/samtools_sort_tmp"
    mkdir -p "$sort_tmp_dir"
    if [[ -s "$bam" ]]; then
      if samtools quickcheck "$bam"; then
        log_msg "[SKIP] valid BAM exists: $bam"
      else
        log_msg "[WARN] existing BAM failed samtools quickcheck; rebuilding: $bam"
        rm -f "$bam" "${bam}.bai"
      fi
    fi
    if [[ ! -s "$bam" ]]; then
      log_msg "[START] mapping and sorting reads: minimap2_threads=${MAP_THREADS}, sort_threads=${SORT_THREADS}, sort_mem=${SAMTOOLS_SORT_MEM}"
      rm -f "$bam_tmp" "${bam_tmp}.bai"
      minimap2 -ax sr -t "$MAP_THREADS" "$filtered" "$r1" "$r2" \
        | samtools view -@ "$SORT_THREADS" -bS - \
        | samtools sort -@ "$SORT_THREADS" -m "$SAMTOOLS_SORT_MEM" -T "${sort_tmp_dir}/${sample}.sort" -o "$bam_tmp" -
      samtools quickcheck "$bam_tmp"
      mv -f "$bam_tmp" "$bam"
      log_msg "[DONE] mapping and sorting reads: $bam"
    fi
    if [[ ! -s "${bam}.bai" ]]; then
      log_msg "[START] samtools index: $bam"
      samtools index "$bam"
      log_msg "[DONE] samtools index: ${bam}.bai"
    fi

    local mos_prefix="${d08}/${sample}"
    local mos_summary="${mos_prefix}.mosdepth.summary.txt"
    if [[ ! -s "$mos_summary" ]]; then
      log_msg "[START] mosdepth"
      mosdepth -t "$THREADS" -n "$mos_prefix" "$bam"
      log_msg "[DONE] mosdepth"
    fi
    awk 'BEGIN{OFS="\t"; print "contig","length","mean_depth"}
         NR>1 && $1!="total" { if ($4 ~ /^[0-9.eE+-]+$/) print $1,$2,$4 }' "$mos_summary" > "${d08}/contig_depth.tsv"

    samtools idxstats "$bam" \
      | awk -v S="$sample" 'BEGIN{OFS="\t"; print "sample","len_bin","contig","contig_length","mapped_reads","unmapped_reads","bam"}
         $1!="*" {
           len=($2+0);
           bin=(len>=1000?"ge1000":(len>=500 && len<=999?"len500_999":"lt500"));
           print S,bin,$1,$2,$3,$4,"'"$bam"'"
         }' > "${d08}/contig_read_counts.tsv"

    local ids_ge1000="${d03}/ge1000.depth_ge${DEPTH_CUTOFF}.ids"
    local ids_len500="${d03}/len500_999.depth_ge${DEPTH_CUTOFF}.ids"
    local fa_ge1000="${d03}/ge1000.depth_ge${DEPTH_CUTOFF}.fa"
    local fa_len500="${d03}/len500_999.depth_ge${DEPTH_CUTOFF}.fa"
    awk -F'\t' -v d="$DEPTH_CUTOFF" 'NR>1 && ($2+0)>=1000 && ($3+0)>=d {print $1}' "${d08}/contig_depth.tsv" | sort -u > "$ids_ge1000"
    awk -F'\t' -v d="$DEPTH_CUTOFF" 'NR>1 && ($2+0)>=500 && ($2+0)<=999 && ($3+0)>=d {print $1}' "${d08}/contig_depth.tsv" | sort -u > "$ids_len500"
    if [[ -s "$ids_ge1000" ]]; then seqkit grep -f "$ids_ge1000" "$ge1000" > "$fa_ge1000" || true; else : > "$fa_ge1000"; fi
    if [[ -s "$ids_len500" ]]; then seqkit grep -f "$ids_len500" "$len500_999" > "$fa_len500" || true; else : > "$fa_len500"; fi

    if [[ ! -f "${d04}/.done" ]]; then
      log_msg "[START] VirSorter2"
      virsorter run -w "$d04" -i "$filtered" --db-dir "$VS2_DB_DIR" --min-length "$VS2_MIN_LEN" -j "$THREADS" all
      touch "${d04}/.done"
      log_msg "[DONE] VirSorter2"
    fi
    if [[ ! -f "${d05}/.done" ]]; then
      log_msg "[START] geNomad"
      genomad end-to-end --threads "$THREADS" "$filtered" "$d05" "$GENOMAD_DB"
      touch "${d05}/.done"
      log_msg "[DONE] geNomad"
    fi
    if [[ ! -f "${d06}/.done" ]]; then
      log_msg "[START] CheckV"
      checkv end_to_end "$filtered" "$d06" -t "$THREADS" -d "$CHECKV_DB"
      touch "${d06}/.done"
      log_msg "[DONE] CheckV"
    fi

    log_msg "[START] DIAMOND bin ge1000"
    run_diamond_bin "$sample" "ge1000" "$fa_ge1000" "$d07x"
    log_msg "[DONE] DIAMOND bin ge1000"
    log_msg "[START] DIAMOND bin len500_999"
    run_diamond_bin "$sample" "len500_999" "$fa_len500" "$d07s"
    log_msg "[DONE] DIAMOND bin len500_999"

    {
      echo -e "sample\t${sample}"
      echo -e "assembly\t${assembly}"
      echo -e "reads_r1\t${r1}"
      echo -e "reads_r2\t${r2}"
      echo -e "bam\t${bam}"
    } > "${sample_dir}/run_summary.tsv"
  } > "$log" 2>&1

  echo -e "${sample}\tDONE\tOK" >> "$STATUS_TSV"
}

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample assembly r1 r2 extra; do
  [[ -n "${sample:-}" ]] || continue
  if ! run_one_sample "$sample" "$assembly" "$r1" "$r2"; then
    echo "[WARN] sample failed: $sample" >&2
  fi
done

READCOUNT_MERGED="${MERGED_DIR}/all_contig_read_counts.tsv"
echo -e "sample\tlen_bin\tcontig\tcontig_length\tmapped_reads\tunmapped_reads\tbam" > "$READCOUNT_MERGED"
find "$SAMPLES_DIR" -path '*/08_coverage/contig_read_counts.tsv' -type f | sort | while read -r f; do
  tail -n +2 "$f" >> "$READCOUNT_MERGED"
done

echo "[OK] wrote $READCOUNT_MERGED"
