#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AIR_STAGE=finalize bash "${SCRIPT_DIR}/air_results_engine.sh"
