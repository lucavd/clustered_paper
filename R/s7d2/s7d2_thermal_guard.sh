#!/bin/bash
# s7d2_thermal_guard.sh -- copy of R/s7e/s7e_thermal_guard.sh for the S7d2 run: adjusts results_v4/S7d2/WORKERS
# (read by s7d2_run.R at every scheduling step; running chunks are never killed).
# Every 240 s: Tctl > 90 C -> workers - 4 (min 8); Tctl < 82 C -> workers + 2 (max set in results_v4/S7d2/WORKERS_MAX, default 25).
CTL=~/MMM/results_v4/S7d2/WORKERS; MAXF=~/MMM/results_v4/S7d2/WORKERS_MAX; RUNLOG=~/MMM/results_v4/S7d2/run_progress.log; LOG=~/MMM/reports/morris/S7d2_logs/workers_changes.log
mkdir -p "$(dirname $LOG)"
while ! grep -q " DONE" "$RUNLOG" 2>/dev/null; do
  T=$(sensors 2>/dev/null | awk '/Tctl/ {gsub(/[+°C]/, "", $2); print int($2)}')
  W=$(cat "$CTL"); MX=$(cat "$MAXF" 2>/dev/null || echo 25)
  if [ -n "$T" ]; then
    if [ "$T" -gt 90 ] && [ "$W" -gt 8 ]; then N=$(( W - 4 )); [ $N -lt 8 ] && N=8; echo $N > "$CTL"; echo "$(date -Iseconds) Tctl ${T}C -> WORKERS $W->$N" >> "$LOG"
    elif [ "$T" -lt 82 ] && [ "$W" -lt "$MX" ]; then N=$(( W + 2 )); [ $N -gt $MX ] && N=$MX; echo $N > "$CTL"; echo "$(date -Iseconds) Tctl ${T}C -> WORKERS $W->$N" >> "$LOG"; fi
  fi
  sleep 240
done
echo "$(date -Iseconds) guard stopped: run DONE" >> "$LOG"
