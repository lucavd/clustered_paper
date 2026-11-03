#!/bin/bash
# s7d_thermal_guard.sh -- copy of R/s6b/s6b_thermal_guard.sh for the S7d run by
# adjusting results_v4/S7d/WORKERS (read by s7d_run.R at every scheduling step; running chunks are never killed).
# Every 240 s: Tctl > 90 C -> workers - 4 (min 8); Tctl < 82 C -> workers + 2 (max 28). Stops when the run log says DONE.
CTL=~/MMM/results_v4/S7d/WORKERS; RUNLOG=~/MMM/results_v4/S7d/run_progress.log; LOG=~/MMM/reports/morris/S7d_logs/workers_changes.log
mkdir -p "$(dirname $LOG)"
while ! grep -q " DONE" "$RUNLOG" 2>/dev/null; do
  T=$(sensors 2>/dev/null | awk '/Tctl/ {gsub(/[+°C]/, "", $2); print int($2)}')
  W=$(cat "$CTL")
  if [ -n "$T" ]; then
    if [ "$T" -gt 90 ] && [ "$W" -gt 8 ]; then N=$(( W - 4 )); [ $N -lt 8 ] && N=8; echo $N > "$CTL"; echo "$(date -Iseconds) Tctl ${T}C -> WORKERS $W->$N" >> "$LOG"
    elif [ "$T" -lt 82 ] && [ "$W" -lt 28 ]; then N=$(( W + 2 )); echo $N > "$CTL"; echo "$(date -Iseconds) Tctl ${T}C -> WORKERS $W->$N" >> "$LOG"; fi
  fi
  sleep 240
done
echo "$(date -Iseconds) guard stopped: run DONE" >> "$LOG"
