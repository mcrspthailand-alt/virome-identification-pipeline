#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/config.env}"

[[ -s "$CONFIG" ]] || { echo "[ERROR] missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

mkdir -p "$OUT_BASE/logs" "$OUT_BASE/merged" "$OUT_BASE/samples" "$OUT_BASE/taxonomy"
LOG="${OUT_BASE}/logs/pipeline.log"

run_step() {
  local name="$1"
  shift
  echo
  echo "============================================================" | tee -a "$LOG"
  echo "[$(date)] START: $name" | tee -a "$LOG"
  echo "============================================================" | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
  echo "[$(date)] DONE: $name" | tee -a "$LOG"
}

echo "=== configurable virome pipeline started: $(date) ===" | tee "$LOG"
echo "[INFO] CONFIG=$CONFIG" | tee -a "$LOG"
echo "[INFO] OUT_BASE=$OUT_BASE" | tee -a "$LOG"

run_step "1 process samples and run detection tools" bash "${SCRIPT_DIR}/scripts/01_process_samples.sh" "$CONFIG"
run_step "2 build merged evidence table" python3 "${SCRIPT_DIR}/scripts/02_build_evidence_table.py" --config "$CONFIG"
run_step "3 assign taxonomy" python3 "${SCRIPT_DIR}/scripts/03_assign_taxonomy.py" --config "$CONFIG"
run_step "4 finalize annotated result tables" python3 "${SCRIPT_DIR}/scripts/04_finalize_results.py" --config "$CONFIG"
run_step "5 write output README" python3 "${SCRIPT_DIR}/scripts/05_write_output_readme.py" --config "$CONFIG"

echo
echo "=== configurable virome pipeline finished: $(date) ===" | tee -a "$LOG"
echo "[DONE] Main outputs:" | tee -a "$LOG"
echo "  ${OUT_BASE}/merged/all_viral_evidence.tsv" | tee -a "$LOG"
echo "  ${OUT_BASE}/merged/final_virus_contigs.tsv" | tee -a "$LOG"
echo "  ${OUT_BASE}/README.md" | tee -a "$LOG"
