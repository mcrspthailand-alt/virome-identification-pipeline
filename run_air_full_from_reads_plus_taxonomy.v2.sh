#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Unified air virome pipeline
# Combines the logic of:
#   1) run_all_bat_samples_full_lenbin.sh
#   2) run_one_bat_sample_full_lenbin.sh
#   3) taxonomy_enrich_tier1.py / make_taxonomy_enriched_matrices.py
#
# IMPORTANT
# - This version starts from assembly FASTA(.gz) + paired reads FASTQ(.gz),
#   because the original run_one_bat_sample_full_lenbin.sh builds the BAM
#   by mapping reads back to contigs.
# - Existing BAM files are optional only for later helper steps, not the main
#   starting input in the original flow.
###############################################################################

# =========================
# USER SETTINGS
# =========================
PROJECT_ROOT="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/air"
RUN_DIR_13_21="${PROJECT_ROOT}/results_13_21"
RUN_DIR_NON_HYBRID="${PROJECT_ROOT}/result_air_non_hybrid_use"
OUT_ROOT="${PROJECT_ROOT}/non_hybrid"
OUT_BASE="${OUT_ROOT}/viral_id_pipeline"
SCRIPT_DIR="${OUT_BASE}/scripts"
LOG_DIR="${OUT_BASE}/logs"
SAMPLES_DIR="${OUT_BASE}/samples"
MERGED_DIR="${OUT_BASE}/merged"

THREADS=32
DEPTH_CUTOFF="2"
VS2_MIN_LEN="500"

# Databases / tools
DIAMOND_DB="/home/panda/workspace/databases/diamond/virus_nr/virus_nr.dmnd"
GENOMAD_DB="/home/panda/workspace/databases/genomad/genomad_db"
VS2_DB_DIR="/home/panda/workspace/databases/virsorter2_db"
CHECKV_DB="/home/panda/workspace/databases/checkv/checkv-db-v1.5"

# Optional taxonomy inputs
WITH_TAXID=""   # e.g. /path/to/tier1_virus_names.with_taxid.tsv
RANKS_TSV=""    # e.g. /path/to/taxid.ranks.tsv
MAP_TSV=""      # e.g. /path/to/virus_taxonomy_map.tsv

# Sample list
SAMPLES=(air13 air14 air15 air17 air19 air20 air21 air23 air24)

sample_run_dir() {
  case "$1" in
    air13|air21) echo "$RUN_DIR_13_21" ;;
    air14|air15|air17|air19|air20|air23|air24) echo "$RUN_DIR_NON_HYBRID" ;;
    *) return 1 ;;
  esac
}

assembly_dir_for_sample() {
  local run_dir
  run_dir="$(sample_run_dir "$1")" || return 1
  echo "${run_dir}/Assembly/MEGAHIT"
}

reads_dir_for_sample() {
  local run_dir
  run_dir="$(sample_run_dir "$1")" || return 1
  echo "${run_dir}/QC_shortreads/remove_host"
}

find_assembly_for_sample() {
  local sample="$1"
  local adir
  adir="$(assembly_dir_for_sample "$sample")" || return 1

  local -a cands=(
    "${adir}/MEGAHIT-${sample}.contigs.fa.gz"
    "${adir}/MEGAHIT-${sample}.contigs.fa"
  )

  local f
  for f in "${cands[@]}"; do
    if [[ -s "$f" ]]; then
      echo "$f"
      return 0
    fi
  done
  return 1
}

# FASTQ naming patterns to try for each sample.
# The first existing pair will be used.
find_reads_for_sample() {
  local sample="$1"
  local READS_DIR
  READS_DIR="$(reads_dir_for_sample "$sample")" || return 1
  local r1="" r2=""

  local -a cands=(
    "${READS_DIR}/${sample}_run0_host_removed.unmapped_1.fastq.gz|${READS_DIR}/${sample}_run0_host_removed.unmapped_2.fastq.gz"
    "${READS_DIR}/${sample}_run0_host_removed.unmapped_1.fq.gz|${READS_DIR}/${sample}_run0_host_removed.unmapped_2.fq.gz"
    "${READS_DIR}/${sample}_R1.fastq.gz|${READS_DIR}/${sample}_R2.fastq.gz"
    "${READS_DIR}/${sample}_R1.fq.gz|${READS_DIR}/${sample}_R2.fq.gz"
    "${READS_DIR}/${sample}.R1.fastq.gz|${READS_DIR}/${sample}.R2.fastq.gz"
    "${READS_DIR}/${sample}.R1.fq.gz|${READS_DIR}/${sample}.R2.fq.gz"
    "${READS_DIR}/${sample}_1.fastq.gz|${READS_DIR}/${sample}_2.fastq.gz"
    "${READS_DIR}/${sample}_1.fq.gz|${READS_DIR}/${sample}_2.fq.gz"
  )

  local pair
  for pair in "${cands[@]}"; do
    r1="${pair%%|*}"
    r2="${pair##*|}"
    if [[ -s "$r1" && -s "$r2" ]]; then
      echo "$r1|$r2"
      return 0
    fi
  done
  return 1
}

# =========================
# TOOL CHECKS
# =========================
for exe in bash awk zcat samtools minimap2 mosdepth prodigal diamond seqkit bbduk.sh python3 checkv; do
  command -v "$exe" >/dev/null 2>&1 || { echo "[ERROR] $exe not found"; exit 1; }
done
command -v genomad >/dev/null 2>&1 || { echo "[ERROR] genomad not found"; exit 1; }
command -v virsorter >/dev/null 2>&1 || { echo "[ERROR] virsorter not found"; exit 1; }

[[ -s "$DIAMOND_DB" ]] || { echo "[ERROR] DIAMOND DB not found: $DIAMOND_DB"; exit 1; }
[[ -d "$GENOMAD_DB" ]] || { echo "[ERROR] geNomad DB dir not found: $GENOMAD_DB"; exit 1; }
[[ -d "$VS2_DB_DIR" ]] || { echo "[ERROR] VirSorter2 DB dir not found: $VS2_DB_DIR"; exit 1; }
[[ -d "$CHECKV_DB" ]] || { echo "[ERROR] CheckV DB dir not found: $CHECKV_DB"; exit 1; }

mkdir -p "$SCRIPT_DIR" "$LOG_DIR" "$SAMPLES_DIR" "$MERGED_DIR"

MASTER_LOG="${LOG_DIR}/air_full.master.log"
STATUS_TSV="${LOG_DIR}/air_full.status.tsv"
echo -e "sample\tstatus\tmessage" > "$STATUS_TSV"
echo "=== air full batch started: $(date) ===" | tee "$MASTER_LOG"

###############################################################################
# Per-sample pipeline function
###############################################################################
run_one_sample() {
  local SAMPLE="$1"
  local SAMPLE_OUT="${SAMPLES_DIR}/${SAMPLE}"
  local SAMPLE_LOG="${LOG_DIR}/${SAMPLE}.full.log"
  mkdir -p "$SAMPLE_OUT"

  local ASM
  if ! ASM="$(find_assembly_for_sample "$SAMPLE")"; then
    echo -e "${SAMPLE}\tFAILED\tmissing assembly in $(assembly_dir_for_sample "$SAMPLE" 2>/dev/null || echo NA)" >> "$STATUS_TSV"
    return 1
  fi

  local pair
  if ! pair="$(find_reads_for_sample "$SAMPLE")"; then
    echo -e "${SAMPLE}\tFAILED\tmissing paired FASTQ in $(reads_dir_for_sample "$SAMPLE" 2>/dev/null || echo NA)" >> "$STATUS_TSV"
    return 1
  fi
  local R1="${pair%%|*}"
  local R2="${pair##*|}"

  if [[ -s "${SAMPLE_OUT}/07_diamond_option1_all_depth_ge${DEPTH_CUTOFF}_ge1000/diamond_out/contig_best_hit.filtered.tsv" ]] \
     && [[ -s "${SAMPLE_OUT}/07_diamond_option1_all_depth_ge${DEPTH_CUTOFF}_len500_999/diamond_out/contig_best_hit.filtered.tsv" ]] \
     && [[ -f "${SAMPLE_OUT}/04_virsorter2/.done" ]] \
     && [[ -f "${SAMPLE_OUT}/05_genomad/.done" ]] \
     && [[ -f "${SAMPLE_OUT}/06_checkv/.done" ]] \
     && [[ -s "${SAMPLE_OUT}/08_coverage/cov.tsv" ]]; then
    echo -e "${SAMPLE}\tSKIP\talready_done" >> "$STATUS_TSV"
    return 0
  fi

  {
    echo "[INFO] sample=$SAMPLE"
    echo "[INFO] R1=$R1"
    echo "[INFO] R2=$R2"
    echo "[INFO] ASM=$ASM"

    local D00="${SAMPLE_OUT}/00_input"
    local D01="${SAMPLE_OUT}/01_filter_lenbin"
    local D02="${SAMPLE_OUT}/02_mapping"
    local D03="${SAMPLE_OUT}/03_depth_tiers"
    local D04="${SAMPLE_OUT}/04_virsorter2"
    local D05="${SAMPLE_OUT}/05_genomad"
    local D06="${SAMPLE_OUT}/06_checkv"
    local D07_500="${SAMPLE_OUT}/07_diamond_option1_all_depth_ge${DEPTH_CUTOFF}_len500_999"
    local D07_1000="${SAMPLE_OUT}/07_diamond_option1_all_depth_ge${DEPTH_CUTOFF}_ge1000"
    local D08="${SAMPLE_OUT}/08_coverage"
    mkdir -p "$D00" "$D01" "$D02" "$D03" "$D04" "$D05" "$D06" "$D07_500" "$D07_1000" "$D08"

    ln -sfn "$R1" "${D00}/${SAMPLE}_R1.fastq.gz"
    ln -sfn "$R2" "${D00}/${SAMPLE}_R2.fastq.gz"
    ln -sfn "$ASM" "${D00}/${SAMPLE}.assembly.contigs.fa${ASM##*.fa}"

    local CONTIGS="${D01}/contigs.fa"
    local FILTERED="${D01}/contigs.filtered.fa"
    local GE1000="${D01}/contigs.filtered.ge1000.fa"
    local LEN500_999="${D01}/contigs.filtered.500_999.fa"
    local SEQKIT_STATS="${D01}/seqkit_stats.txt"

    if [[ ! -s "$CONTIGS" ]]; then
      if [[ "$ASM" == *.gz ]]; then
        zcat "$ASM" > "$CONTIGS"
      else
        cp -f "$ASM" "$CONTIGS"
      fi
    fi
    [[ -s "$CONTIGS" ]] || { echo "[ERROR] empty contigs.fa"; exit 1; }

    if [[ ! -s "$FILTERED" ]]; then
      bbduk.sh in="$CONTIGS" out="$FILTERED" entropy=0.70 entropywindow=50 entropyk=5 overwrite=t
    fi
    if [[ ! -s "$GE1000" ]]; then
      seqkit seq -m 1000 "$FILTERED" > "$GE1000"
    fi
    if [[ ! -s "$LEN500_999" ]]; then
      seqkit seq -m 500 "$FILTERED" | seqkit seq -M 999 > "$LEN500_999"
    fi
    seqkit stats "$FILTERED" "$GE1000" "$LEN500_999" > "$SEQKIT_STATS"

    local BAM="${D02}/${SAMPLE}.filtered.sorted.bam"
    local BAI="${BAM}.bai"
    if [[ ! -s "$BAM" ]]; then
      minimap2 -ax sr -t "$THREADS" "$FILTERED" "$R1" "$R2" \
        | samtools view -@ "$THREADS" -bS - \
        | samtools sort -@ "$THREADS" -o "$BAM" -
    fi
    [[ -s "$BAM" ]] || { echo "[ERROR] BAM not created"; exit 1; }
    [[ -s "$BAI" ]] || samtools index "$BAM"

    local MOS_PREFIX="${D08}/${SAMPLE}"
    local MOS_SUMMARY="${MOS_PREFIX}.mosdepth.summary.txt"
    local COV_TSV="${D08}/cov.tsv"
    if [[ ! -s "$MOS_SUMMARY" ]]; then
      mosdepth -t "$THREADS" -n "$MOS_PREFIX" "$BAM"
    fi
    [[ -s "$MOS_SUMMARY" ]] || { echo "[ERROR] mosdepth summary missing"; exit 1; }
    awk 'BEGIN{OFS="\t"; print "contig","length","mean_depth"}
         NR>1 && $1!="total" { if ($4 ~ /^[0-9.eE+-]+$/) print $1,$2,$4 }' "$MOS_SUMMARY" > "$COV_TSV"

    # idxstats readcounts by lenbin
    samtools idxstats "$BAM" \
      | awk -v S="$SAMPLE" 'BEGIN{OFS="\t"; print "sample","len_bin","contig","contig_length","mapped_reads","unmapped_reads","bam"}
         $1!="*" {
           len=($2+0);
           bin=(len>=1000?"ge1000":(len>=500 && len<=999?"len500_999":"lt500"));
           print S,bin,$1,$2,$3,$4,"'"$BAM"'"
         }' > "${D08}/readcounts.by_lenbin.tsv"

    local DEPTH_GE2_IDS_GE1000="${D03}/contigs.ge1000.depth_ge${DEPTH_CUTOFF}.ids"
    local DEPTH_GE2_FA_GE1000="${D03}/contigs.ge1000.depth_ge${DEPTH_CUTOFF}.fa"
    local DEPTH_GE2_IDS_500="${D03}/contigs.len500_999.depth_ge${DEPTH_CUTOFF}.ids"
    local DEPTH_GE2_FA_500="${D03}/contigs.len500_999.depth_ge${DEPTH_CUTOFF}.fa"

    awk -F'\t' -v d="$DEPTH_CUTOFF" 'NR>1 && ($2+0)>=1000 && ($3+0)>=d {print $1}' "$COV_TSV" | sort -u > "$DEPTH_GE2_IDS_GE1000"
    awk -F'\t' -v d="$DEPTH_CUTOFF" 'NR>1 && ($2+0)>=500 && ($2+0)<=999 && ($3+0)>=d {print $1}' "$COV_TSV" | sort -u > "$DEPTH_GE2_IDS_500"

    if [[ -s "$DEPTH_GE2_IDS_GE1000" ]]; then seqkit grep -f "$DEPTH_GE2_IDS_GE1000" "$GE1000" > "$DEPTH_GE2_FA_GE1000" || true; else : > "$DEPTH_GE2_FA_GE1000"; fi
    if [[ -s "$DEPTH_GE2_IDS_500" ]]; then seqkit grep -f "$DEPTH_GE2_IDS_500" "$LEN500_999" > "$DEPTH_GE2_FA_500" || true; else : > "$DEPTH_GE2_FA_500"; fi

    if [[ ! -f "${D04}/.done" ]]; then
      virsorter run -w "$D04" -i "$FILTERED" --db-dir "$VS2_DB_DIR" --min-length "$VS2_MIN_LEN" -j "$THREADS" all
      touch "${D04}/.done"
    fi

    if [[ ! -f "${D05}/.done" ]]; then
      genomad end-to-end --threads "$THREADS" "$FILTERED" "$D05" "$GENOMAD_DB"
      touch "${D05}/.done"
    fi

    if [[ ! -f "${D06}/.done" ]]; then
      checkv end_to_end "$FILTERED" "$D06" -t "$THREADS" -d "$CHECKV_DB"
      touch "${D06}/.done"
    fi

    run_option1_diamond() {
      local LEN_BIN="$1"
      local DEPTH_FA="${D03}/contigs.${LEN_BIN}.depth_ge${DEPTH_CUTOFF}.fa"
      local D07="${SAMPLE_OUT}/07_diamond_option1_all_depth_ge${DEPTH_CUTOFF}_${LEN_BIN}"
      local PROD="${D07}/prodigal_out"
      local DOUT="${D07}/diamond_out"
      mkdir -p "$PROD" "$DOUT"

      local DIAMOND_BEST="${DOUT}/contig_best_hit.tsv"
      local DIAMOND_BEST_FILT="${DOUT}/contig_best_hit.filtered.tsv"
      if [[ ! -s "$DEPTH_FA" ]]; then
        echo "[WARN] ${SAMPLE}: empty ${DEPTH_FA} -> skip DIAMOND ${LEN_BIN}"
        echo -e "contig\tsseqid\tpident\taln_len\tevalue\tbitscore\tstitle\tn_orfs_strict" > "$DIAMOND_BEST"
        cp -f "$DIAMOND_BEST" "$DIAMOND_BEST_FILT"
        return 0
      fi

      local PROT_FAA="${PROD}/${SAMPLE}.${LEN_BIN}.depth_ge${DEPTH_CUTOFF}.proteins.faa"
      local GFF="${PROD}/${SAMPLE}.${LEN_BIN}.depth_ge${DEPTH_CUTOFF}.gff"
      local DIAMOND_RAW="${DOUT}/viral_protein_hits.tsv"
      local EVALUE="1e-10"
      local BITSCORE="100"
      local ALNLEN="60"

      [[ -s "$PROT_FAA" ]] || prodigal -i "$DEPTH_FA" -a "$PROT_FAA" -o "$GFF" -p meta
      if [[ ! -s "$PROT_FAA" ]]; then
        echo -e "contig\tsseqid\tpident\taln_len\tevalue\tbitscore\tstitle\tn_orfs_strict" > "$DIAMOND_BEST"
        cp -f "$DIAMOND_BEST" "$DIAMOND_BEST_FILT"
        return 0
      fi

      [[ -s "$DIAMOND_RAW" ]] || diamond blastp -d "$DIAMOND_DB" -q "$PROT_FAA" -o "$DIAMOND_RAW" -f 6 qseqid sseqid pident length qlen slen evalue bitscore stitle --evalue "$EVALUE" --max-target-seqs 5 --threads "$THREADS"
      if [[ ! -s "$DIAMOND_RAW" ]]; then
        echo -e "contig\tsseqid\tpident\taln_len\tevalue\tbitscore\tstitle\tn_orfs_strict" > "$DIAMOND_BEST"
        cp -f "$DIAMOND_BEST" "$DIAMOND_BEST_FILT"
        return 0
      fi

      awk -F'\t' 'BEGIN{ OFS="\t"; print "contig","sseqid","pident","aln_len","evalue","bitscore","stitle","n_orfs_strict" }
      {
        q=$1; contig=q; sub(/_[0-9]+$/, "", contig)
        e=$7+0; b=$8+0; al=$4+0; st=$9; for (i=10; i<=NF; i++) st = st " " $i
        if (e <= '"${EVALUE}"' && b >= '"${BITSCORE}"' && al >= '"${ALNLEN}"') {
          key = contig SUBSEP q; if (!(key in seen_orf)) { seen_orf[key]=1; strict_cnt[contig]++ }
        }
        if (!(contig in best_b) || b > best_b[contig] || (b==best_b[contig] && e < best_e[contig])) {
          best_b[contig]=b; best_e[contig]=e; best_line[contig]=contig OFS $2 OFS $3 OFS $4 OFS $7 OFS $8 OFS st
        }
      }
      END{ for (c in best_line) { n=(c in strict_cnt)?strict_cnt[c]:0; print best_line[c], n } }' "$DIAMOND_RAW" | sort > "$DIAMOND_BEST"

      awk -F'\t' 'NR==1 || ($5 <= '"${EVALUE}"' && $6 >= '"${BITSCORE}"' && $4 >= '"${ALNLEN}"')' "$DIAMOND_BEST" > "$DIAMOND_BEST_FILT"
    }

    run_option1_diamond "ge1000"
    run_option1_diamond "len500_999"

    {
      echo -e "sample\t${SAMPLE}"
      echo -e "assembly\t${ASM}"
      echo -e "reads_r1\t${R1}"
      echo -e "reads_r2\t${R2}"
      echo -e "bam\t${BAM}"
      echo -e "cov_tsv\t${COV_TSV}"
    } > "${SAMPLE_OUT}/run_summary.tsv"
  } > "$SAMPLE_LOG" 2>&1

  echo -e "${SAMPLE}\tDONE\tOK" >> "$STATUS_TSV"
}

###############################################################################
# Run batch
###############################################################################
for SAMPLE in "${SAMPLES[@]}"; do
  if ! run_one_sample "$SAMPLE"; then
    echo "[WARN] sample failed: $SAMPLE" | tee -a "$MASTER_LOG"
  fi
done

echo "=== per-sample stage finished: $(date) ===" | tee -a "$MASTER_LOG"

###############################################################################
# Merge stage
###############################################################################
READCOUNT_MERGED="${MERGED_DIR}/all_contig_read_counts.tsv"
echo -e "sample\tlen_bin\tcontig\tcontig_length\tmapped_reads\tunmapped_reads\tbam" > "$READCOUNT_MERGED"
find "$SAMPLES_DIR" -path '*/08_coverage/readcounts.by_lenbin.tsv' -type f | sort | while read -r f; do
  tail -n +2 "$f" >> "$READCOUNT_MERGED"
done

MASTER_COMBINED="${MERGED_DIR}/diamond_best_hits_initial.tsv"
python3 - <<'PY' "$SAMPLES_DIR" "$MASTER_COMBINED"
import os, sys, glob, pandas as pd
samples_dir, out = sys.argv[1], sys.argv[2]
frames = []
for sample_dir in sorted(glob.glob(os.path.join(samples_dir, '*'))):
    sample = os.path.basename(sample_dir)
    for len_bin, p in [
        ('ge1000', os.path.join(sample_dir, '07_diamond_option1_all_depth_ge2_ge1000', 'diamond_out', 'contig_best_hit.filtered.tsv')),
        ('len500_999', os.path.join(sample_dir, '07_diamond_option1_all_depth_ge2_len500_999', 'diamond_out', 'contig_best_hit.filtered.tsv')),
    ]:
        if os.path.exists(p):
            df = pd.read_csv(p, sep='\t', dtype=str)
            if len(df) == 0:
                continue
            df['sample'] = sample
            df['len_bin'] = len_bin
            df['blastx_organism'] = df.get('stitle', '').fillna('').astype(str).str.strip()
            frames.append(df)
if frames:
    outdf = pd.concat(frames, ignore_index=True)
else:
    outdf = pd.DataFrame(columns=['contig','sseqid','pident','aln_len','evalue','bitscore','stitle','n_orfs_strict','sample','len_bin','blastx_organism'])
outdf.to_csv(out, sep='\t', index=False)
print('[OK] wrote', out, 'rows=', len(outdf))
PY

MASTER_WITH_MAPPED="${MERGED_DIR}/diamond_best_hits_with_read_counts.tsv"
python3 - <<'PY' "$MASTER_COMBINED" "$READCOUNT_MERGED" "$MASTER_WITH_MAPPED"
import sys, pandas as pd
MASTER, READS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
m = pd.read_csv(MASTER, sep='\t', dtype=str, low_memory=False)
r = pd.read_csv(READS, sep='\t', dtype=str, low_memory=False)
if 'length' not in m.columns:
    m['length'] = ''
for c in ['sample','contig','len_bin']:
    if c not in m.columns: m[c] = ''
    m[c] = m[c].astype(str).str.strip()
for c in ['sample','contig','len_bin']:
    r[c] = r[c].astype(str).str.strip()
r['mapped_reads'] = pd.to_numeric(r['mapped_reads'], errors='coerce').fillna(0).astype(int)
out = m.merge(r[['sample','len_bin','contig','mapped_reads']], on=['sample','len_bin','contig'], how='left')
out['mapped_reads'] = pd.to_numeric(out['mapped_reads'], errors='coerce').fillna(0).astype(int)
out.to_csv(OUT, sep='\t', index=False)
print('[OK] wrote', OUT, 'rows=', len(out))
PY

TIER1_READS_MATRIX="${MERGED_DIR}/virus_reads_matrix_initial.tsv"
python3 - <<'PY' "$MASTER_WITH_MAPPED" "$TIER1_READS_MATRIX"
import sys, pandas as pd
inp, out = sys.argv[1], sys.argv[2]
df = pd.read_csv(inp, sep='\t', dtype=str, low_memory=False)
if len(df) == 0:
    pd.DataFrame(columns=['blastx_organism']).to_csv(out, sep='\t', index=False)
    print('[OK] wrote empty', out)
    raise SystemExit(0)
if 'blastx_organism' not in df.columns:
    df['blastx_organism'] = df.get('stitle', '').fillna('').astype(str)
if 'mapped_reads' not in df.columns:
    df['mapped_reads'] = 0

df['blastx_organism'] = df['blastx_organism'].fillna('').astype(str).str.strip()
df['sample'] = df['sample'].fillna('').astype(str).str.strip()
df['mapped_reads'] = pd.to_numeric(df['mapped_reads'], errors='coerce').fillna(0)

g = df.groupby(['blastx_organism','sample'], as_index=False)['mapped_reads'].sum()
mat = g.pivot(index='blastx_organism', columns='sample', values='mapped_reads').fillna(0).astype(int).reset_index()
mat.to_csv(out, sep='\t', index=False)
print('[OK] wrote', out, 'rows=', len(mat))
PY

###############################################################################
# Optional taxonomy enrichment stage
###############################################################################
if [[ -n "$MAP_TSV" || ( -n "$WITH_TAXID" && -n "$RANKS_TSV" ) ]]; then
  if [[ -z "$MAP_TSV" ]]; then
    MAP_TSV="${MERGED_DIR}/virus_taxonomy_map.tsv"
    python3 - <<'PY' "$WITH_TAXID" "$RANKS_TSV" "$MAP_TSV"
import sys, pandas as pd
with_taxid, ranks_tsv, out_map = sys.argv[1], sys.argv[2], sys.argv[3]
vt = pd.read_csv(with_taxid, sep='\t', header=None, names=['blastx_organism','taxid'], dtype=str)
rk = pd.read_csv(ranks_tsv, sep='\t', header=None, names=['taxid','kingdom','phylum','class','order','family','genus','species'], dtype=str)
for df in (vt, rk):
    for c in df.columns:
        df[c] = df[c].fillna('').astype(str).str.strip()
m = vt.merge(rk, on='taxid', how='left')
def pick_best_rank(row):
    for rank in ['family','order','class','phylum','kingdom']:
        v = str(row.get(rank, '') or '').strip()
        if v and v.lower() not in {'na','n/a','none'}:
            return pd.Series([v, rank])
    return pd.Series(['Unassigned', 'none'])
res = m.apply(pick_best_rank, axis=1)
m['family_or_best'] = res[0]
m['best_rank'] = res[1]
def mk_lineage(row):
    parts=[]
    for k in ['kingdom','phylum','class','order','family','genus','species']:
        v = str(row.get(k, '') or '').strip()
        if v and v.lower() not in {'na','n/a','none'}:
            parts.append(v)
    return ';'.join(parts)
m['lineage_ranks'] = m.apply(mk_lineage, axis=1)
out = m[['blastx_organism','taxid','family_or_best','best_rank','lineage_ranks']].copy()
out.to_csv(out_map, sep='\t', index=False)
print('[OK] wrote', out_map, 'rows=', len(out))
PY
  fi

  python3 - <<'PY' "$TIER1_READS_MATRIX" "$MAP_TSV" "$MERGED_DIR"
import sys, os, pandas as pd
reads_matrix, map_tsv, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
mat = pd.read_csv(reads_matrix, sep='\t', low_memory=False)
if len(mat) == 0:
    print('[WARN] taxonomy stage skipped because reads matrix is empty')
    raise SystemExit(0)
virus_col = mat.columns[0]
mat = mat.rename(columns={virus_col: 'blastx_organism'})
sample_cols = [c for c in mat.columns if c != 'blastx_organism']
for c in sample_cols:
    mat[c] = pd.to_numeric(mat[c], errors='coerce').fillna(0).astype(int)
mp = pd.read_csv(map_tsv, sep='\t', dtype=str, low_memory=False)
df = mat.merge(mp, on='blastx_organism', how='left')
for c in ['taxid','family_or_best','best_rank','lineage_ranks']:
    if c not in df.columns: df[c] = ''
df['taxid'] = df['taxid'].fillna('').astype(str)
df['family_or_best'] = df['family_or_best'].fillna('Unassigned')
df['best_rank'] = df['best_rank'].fillna('none')
df['lineage_ranks'] = df['lineage_ranks'].fillna('')
meta_cols = ['blastx_organism','taxid','family_or_best','best_rank','lineage_ranks']
df = df[meta_cols + sample_cols]
out_reads = os.path.join(outdir, 'virus_reads_matrix_by_taxonomy.tsv')
out_pa = os.path.join(outdir, 'virus_presence_absence_matrix.tsv')
out_long = os.path.join(outdir, 'virus_reads_by_sample_long.tsv')
df.to_csv(out_reads, sep='\t', index=False)
pa = df.copy(); pa[sample_cols] = (pa[sample_cols] > 0).astype(int); pa.to_csv(out_pa, sep='\t', index=False)
long = df.melt(id_vars=meta_cols, value_vars=sample_cols, var_name='sample', value_name='reads'); long.to_csv(out_long, sep='\t', index=False)
print('[OK] wrote', out_reads)
print('[OK] wrote', out_pa)
print('[OK] wrote', out_long)
PY
fi

echo "=== air full batch finished: $(date) ===" | tee -a "$MASTER_LOG"
echo "Status file: ${STATUS_TSV}" | tee -a "$MASTER_LOG"
