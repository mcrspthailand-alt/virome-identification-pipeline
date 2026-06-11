#!/usr/bin/env bash
set -euo pipefail

############################################
# GLOBAL PARAMETERS
############################################
THREADS=32
MIN_CONTIG_LEN=200

############################################
# INPUTS (from nf-core)
############################################
FASTP_BASE="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/SED_EF/QC_shortreads/fastp"

############################################
# OUTPUT BASE
############################################
OUT_BASE="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/SED_EF/overnight_megahit_blast_pipeline"
mkdir -p "$OUT_BASE"

############################################
# DATABASES
############################################
BLASTN_REFSEQ_DB="/home/panda/workspace/databases/blast/virus_refseq/ref_viruses_rep_genomes"
DIAMOND_DB="/home/panda/workspace/databases/diamond/RVDB/RVDB_U_virus_nr.dmnd"

############################################
# SAMPLES
############################################
SAMPLES=(sd01 sd05)

############################################
# TOOL CHECK
############################################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }; }

need megahit
need seqkit
need blastn
need diamond
need bowtie2
need bowtie2-build
need samtools
need awk
need sort
need gzip

############################################
# helper: find R1/R2 robustly
############################################
find_pair_reads () {
  local sdir="$1"
  local sample="$2"

  # common patterns from nf-core/mag fastp outputs
  local r1_candidates=(
    "${sdir}/${sample}_run0_fastp_1.fastp.fastq.gz"
    "${sdir}/${sample}_run0_fastp_1.fastq.gz"
    "${sdir}/${sample}*fastp*1*.fastq.gz"
    "${sdir}/${sample}*R1*.fastq.gz"
    "${sdir}/*1*.fastq.gz"
  )
  local r2_candidates=(
    "${sdir}/${sample}_run0_fastp_2.fastp.fastq.gz"
    "${sdir}/${sample}_run0_fastp_2.fastq.gz"
    "${sdir}/${sample}*fastp*2*.fastq.gz"
    "${sdir}/${sample}*R2*.fastq.gz"
    "${sdir}/*2*.fastq.gz"
  )

  local R1=""
  local R2=""

  for p in "${r1_candidates[@]}"; do
    for f in $p; do
      [[ -s "$f" ]] && { R1="$f"; break 2; }
    done
  done

  for p in "${r2_candidates[@]}"; do
    for f in $p; do
      [[ -s "$f" ]] && { R2="$f"; break 2; }
    done
  done

  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: cannot find paired fastq.gz for sample=${sample} in ${sdir}" >&2
    echo "  Found files:" >&2
    ls -lh "$sdir" >&2 || true
    exit 1
  fi

  echo "$R1" "$R2"
}

############################################
# MAIN LOOP
############################################
for SAMPLE in "${SAMPLES[@]}"; do
  echo "=================================================="
  echo ">>> SAMPLE: $SAMPLE"
  echo "=================================================="

  ##########################################
  # INPUT READS
  ##########################################
  SDIR="${FASTP_BASE}/${SAMPLE}"
  [[ -d "$SDIR" ]] || { echo "ERROR: missing dir $SDIR" >&2; exit 1; }

  read -r R1 R2 < <(find_pair_reads "$SDIR" "$SAMPLE")
  echo "[${SAMPLE}] R1: $R1"
  echo "[${SAMPLE}] R2: $R2"

  ##########################################
  # OUTPUT STRUCTURE
  ##########################################
  OD="${OUT_BASE}/${SAMPLE}"
  AS="${OD}/assembly"
  CTG="${OD}/contigs"
  BLN="${OD}/blastn_refseq"
  DIA="${OD}/diamond"
  MAP="${OD}/mapping"
  SUM="${OD}/summary"
  LOG="${OD}/logs"

  mkdir -p "$AS" "$CTG" "$BLN" "$DIA" "$MAP" "$SUM" "$LOG"

  ##########################################
  # 1) MEGAHIT ASSEMBLY
  ##########################################
  ASDIR="${AS}/assembly_${SAMPLE}"
  CONTIGS="${ASDIR}/final.contigs.fa"

  if [[ ! -s "$CONTIGS" ]]; then
    echo "[${SAMPLE}] MEGAHIT assembly (min contig ${MIN_CONTIG_LEN})"
    megahit \
      -1 "$R1" \
      -2 "$R2" \
      -o "$ASDIR" \
      -t "$THREADS" \
      --min-contig-len "$MIN_CONTIG_LEN" \
      --presets meta-sensitive \
      1> "${LOG}/${SAMPLE}.megahit.stdout.log" \
      2> "${LOG}/${SAMPLE}.megahit.stderr.log"
  fi

  [[ -s "$CONTIGS" ]] || { echo "ERROR: assembly failed for $SAMPLE" >&2; exit 1; }

  ##########################################
  # 2) FILTER CONTIGS >= MIN_CONTIG_LEN
  ##########################################
  QFA="${CTG}/${SAMPLE}_contigs_ge${MIN_CONTIG_LEN}.fa"
  if [[ ! -s "$QFA" ]]; then
    echo "[${SAMPLE}] Filter contigs >= ${MIN_CONTIG_LEN}"
    seqkit seq -m "$MIN_CONTIG_LEN" "$CONTIGS" > "$QFA"
  fi

  ##########################################
  # 3) BLASTn vs RefSeq viruses
  ##########################################
  BLASTN_OUT="${BLN}/${SAMPLE}_vs_refseqvir.m8"
  if [[ ! -s "$BLASTN_OUT" ]]; then
    echo "[${SAMPLE}] BLASTn vs RefSeq viruses"
    blastn \
      -task blastn \
      -query "$QFA" \
      -db "$BLASTN_REFSEQ_DB" \
      -num_threads "$THREADS" \
      -max_target_seqs 10 \
      -evalue 1e-10 \
      -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" \
      > "$BLASTN_OUT"
  fi

  BESTN="${SUM}/${SAMPLE}.blastn_refseq.besthit.tsv"
  if [[ ! -s "$BESTN" ]]; then
    sort -k1,1 -k13,13gr "$BLASTN_OUT" | awk '!seen[$1]++' > "$BESTN"
  fi

  ##########################################
  # 4) DIAMOND blastx vs RVDB
  ##########################################
  DIAMOND_OUT="${DIA}/${SAMPLE}_vs_rvdb.m8"
  if [[ ! -s "$DIAMOND_OUT" ]]; then
    echo "[${SAMPLE}] DIAMOND blastx vs RVDB"
    diamond blastx \
      -d "$DIAMOND_DB" \
      -q "$QFA" \
      -o "$DIAMOND_OUT" \
      -p "$THREADS" \
      --sensitive \
      -e 1e-5 \
      --max-target-seqs 10 \
      --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle
  fi

  ##########################################
  # 5) MAP READS BACK TO CONTIGS
  ##########################################
  IDX="${MAP}/${SAMPLE}_contigs_index"
  BAM="${MAP}/${SAMPLE}_vs_contigs.sorted.bam"

  if [[ ! -s "${IDX}.1.bt2" && ! -s "${IDX}.1.bt2l" ]]; then
    bowtie2-build "$QFA" "$IDX" > "${LOG}/${SAMPLE}_bowtie2build.log" 2>&1
  fi

  if [[ ! -s "$BAM" ]]; then
    bowtie2 -p "$THREADS" -x "$IDX" -1 "$R1" -2 "$R2" \
      2> "${LOG}/${SAMPLE}_map.log" \
    | samtools view -bS - \
    | samtools sort -@ 8 -o "$BAM" -
  fi

  samtools index "$BAM"
  samtools idxstats "$BAM" > "${MAP}/${SAMPLE}.idxstats.tsv"

  samtools depth -aa "$BAM" \
    | awk '{cov[$1]+=$3; n[$1]++} END{for (c in n) printf "%s\t%d\t%.6f\n", c, n[c], cov[c]/n[c]}' \
    | sort -k1,1 \
    > "${MAP}/${SAMPLE}.mean_depth.tsv"

  ##########################################
  # 6) MERGE FINAL TABLE
  ##########################################
  sort -k1,1 "$BESTN" > "${SUM}/blastn.sorted.tsv"
  sort -k1,1 "${MAP}/${SAMPLE}.idxstats.tsv" | awk '{print $1,$2,$3}' OFS="\t" > "${SUM}/idx.sorted.tsv"
  sort -k1,1 "${MAP}/${SAMPLE}.mean_depth.tsv" | awk '{print $1,$3}' OFS="\t" > "${SUM}/depth.sorted.tsv"

  sort -k1,1 "$DIAMOND_OUT" | awk '!seen[$1]++' > "${SUM}/blastx.best.tsv"

  join -t $'\t' -a1 -e NA -o auto "${SUM}/blastn.sorted.tsv" "${SUM}/idx.sorted.tsv" \
    | join -t $'\t' -a1 -e NA -o auto - "${SUM}/depth.sorted.tsv" \
    | join -t $'\t' -a1 -e NA -o auto - "${SUM}/blastx.best.tsv" \
    > "${SUM}/${SAMPLE}.final.noheader.tsv"

  {
    echo -e "contig_id\tblastn_ref_qlen\tblastn_ref_sseqid\tblastn_ref_pident\tblastn_ref_align_len\tblastn_ref_mismatch\tblastn_ref_gapopen\tblastn_ref_qstart\tblastn_ref_qend\tblastn_ref_sstart\tblastn_ref_send\tblastn_ref_evalue\tblastn_ref_bitscore\tblastn_ref_qcov\tblastn_ref_stitle\tcontig_length\tmapped_reads\tmean_depth\tblastx_qlen\tblastx_sseqid\tblastx_slen\tblastx_pident\tblastx_align_len\tblastx_mismatch\tblastx_gapopen\tblastx_qstart\tblastx_qend\tblastx_sstart\tblastx_send\tblastx_evalue\tblastx_bitscore\tblastx_qcov\tblastx_stitle"
    cat "${SUM}/${SAMPLE}.final.noheader.tsv"
  } > "${SUM}/${SAMPLE}.final_table.tsv"

  echo "[${SAMPLE}] DONE"
  echo
done

echo "?? ALL SAMPLES FINISHED. You can wake up and drink coffee ?"