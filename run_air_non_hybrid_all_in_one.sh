#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Air non-hybrid virome pipeline: one-command runner
#
# This wrapper runs the full workflow from per-sample processing to the final
# taxonomy-annotated table. Sankey outputs are intentionally excluded.
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="/home/panda/workspace/projects/nextflow_metagenomic/eDNA/air"
OUT_BASE="${PROJECT_ROOT}/non_hybrid/viral_id_pipeline"
LOG_DIR="${OUT_BASE}/logs"
MERGED_DIR="${OUT_BASE}/merged"
mkdir -p "$LOG_DIR" "$MERGED_DIR"

MASTER_LOG="${LOG_DIR}/air_non_hybrid_all_in_one.log"

RUN_FULL="${SCRIPT_DIR}/01_process_samples.sh"
BUILD_EVIDENCE="${SCRIPT_DIR}/02_build_evidence_table.sh"
TAXONOMY="${SCRIPT_DIR}/03_assign_taxonomy.py"
FINALIZE="${SCRIPT_DIR}/04_finalize_results.sh"
WRITE_README="${SCRIPT_DIR}/05_write_readme.py"

need_file() {
  local f="$1"
  [[ -s "$f" ]] || { echo "[ERROR] missing required script: $f" >&2; exit 1; }
}

need_cmd() {
  local exe="$1"
  command -v "$exe" >/dev/null 2>&1 || { echo "[ERROR] missing command: $exe" >&2; exit 1; }
}

run_step() {
  local name="$1"
  shift

  echo
  echo "============================================================" | tee -a "$MASTER_LOG"
  echo "[$(date)] START: ${name}" | tee -a "$MASTER_LOG"
  echo "============================================================" | tee -a "$MASTER_LOG"

  "$@" 2>&1 | tee -a "$MASTER_LOG"

  echo "[$(date)] DONE: ${name}" | tee -a "$MASTER_LOG"
}

need_file "$RUN_FULL"
need_file "$BUILD_EVIDENCE"
need_file "$TAXONOMY"
need_file "$FINALIZE"
need_file "$WRITE_README"
need_cmd bash
need_cmd python3

echo "=== air non-hybrid all-in-one pipeline started: $(date) ===" | tee "$MASTER_LOG"
echo "[INFO] Recommended command: bash run_air_virome_pipeline.sh" | tee -a "$MASTER_LOG"
echo "[INFO] SCRIPT_DIR=${SCRIPT_DIR}" | tee -a "$MASTER_LOG"
echo "[INFO] OUT_BASE=${OUT_BASE}" | tee -a "$MASTER_LOG"
echo "[INFO] MERGED_DIR=${MERGED_DIR}" | tee -a "$MASTER_LOG"

run_step "1 per-sample viral workflow" bash "$RUN_FULL"

run_step "2 build DIAMOND/evidence/CheckV master table" bash "$BUILD_EVIDENCE"

run_step "3 build taxonomy map and taxonomy matrices" python3 "$TAXONOMY"

run_step "4 merge taxonomy into final annotated tables" bash "$FINALIZE"

run_step "5 write output README" python3 "$WRITE_README"

echo
echo "=== air non-hybrid all-in-one pipeline finished: $(date) ===" | tee -a "$MASTER_LOG"
echo "[DONE] Main outputs:" | tee -a "$MASTER_LOG"
echo "  ${MERGED_DIR}/all_viral_evidence.tsv" | tee -a "$MASTER_LOG"
echo "  ${MERGED_DIR}/final_virus_contigs.tsv" | tee -a "$MASTER_LOG"
echo "  ${MERGED_DIR}/final_virus_contigs_bacteriophage.tsv" | tee -a "$MASTER_LOG"
echo "  ${MERGED_DIR}/final_virus_contigs_eukaryotic.tsv" | tee -a "$MASTER_LOG"
echo "  ${MERGED_DIR}/final_virus_contigs_archaeal.tsv" | tee -a "$MASTER_LOG"
echo "  ${MERGED_DIR}/final_virus_contigs_unidentified.tsv" | tee -a "$MASTER_LOG"
echo "  ${OUT_BASE}/README.md" | tee -a "$MASTER_LOG"
echo "Log file: ${MASTER_LOG}" | tee -a "$MASTER_LOG"
