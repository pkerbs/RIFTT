#!/bin/bash
# RIFTT Part 2 - knowledge-based filtering and evidence scoring over a whole cohort.
# Usage:  ./RIFTT_part2.sh
# Configuration is read from config/params.conf (copy config/params.conf.example).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${RIFTT_CONFIG:-$HERE/config/params.conf}"
SIF="${RIFTT_SIF:-$HERE/RIFTT.sif}"

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; echo "  cp config/params.conf.example config/params.conf  and edit it."; exit 1; }
[[ -f "$SIF" ]]    || { echo "ERROR: container not found: $SIF (build it or set RIFTT_SIF)"; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"

for v in ANNO OUTPUTFOLDER CLINTABLE; do
  [[ -n "${!v:-}" ]] || { echo "ERROR: $v is not set in $CONFIG"; exit 1; }
  [[ -e "${!v}" ]]   || { echo "ERROR: $v path does not exist: ${!v}"; exit 1; }
done

binds="$OUTPUTFOLDER:/outputfolder,$ANNO:/anno.gtf,$CLINTABLE:/clintable.xlsx"
if [[ -n "${USER_BL:-}" && -f "$USER_BL" ]]; then
  binds="$binds,$USER_BL:/blacklist_user.xlsx"
fi

singularity exec \
  --bind "$binds" \
  --env "internal_BL=${INTERNAL_BL:-1},threads=${THREADS:-16}" \
  "$SIF" \
  /scripts/filter.sh
