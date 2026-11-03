#!/bin/bash
# s6b_thermal_guard.sh -- keeps the CPU temperature of the shared host under control during the S6b part B run by
# adjusting results_v4/S6b/ns/WORKERS (read by s6b_ns_run.R at every scheduling step; running chunks are never killed).
# Every 240 s: Tctl > 90 C -> workers - 4 (min 8); Tctl < 82 C -> workers + 2 (max 28). Stops when the run log says DONE.
CTL=~/MMM/results_v4/S6b/ns/WORKERS; RUNLOG=~/MMM/results_v4/S6b/ns/run_progress.log; LOG=~/MMM/reports/morris/S6b_logs/workers_changes.log
while ! grep -q DONE "$RUNLOG" 2>/dev/null; do
  T=$(sensors 2>/dev/null | awk '/Tctl/ {gsub(/[+°C]/, "", $2); print int($2)}')
  W=$(cat "$CTL")
  if [ -n "$T" ]; then
    if [ "$T" -gt 90 ] && [ "$W" -gt 8 ]; then N=$(( W - 4 )); [ $N -lt 8 ] && N=8; echo $N > "$CTL"; echo "$(date -Iseconds) Tctl ${T}C -> WORKERS $W->$N" >> "$LOG"
    elif [ "$T" -lt 82 ] && [ "$W" -lt 28 ]; then N=$(( W + 2 )); echo $N > "$CTL"; echo "$(date -Iseconds) Tctl ${T}C -> WORKERS $W->$N" >> "$LOG"; fi
  fi
  sleep 240
done
echo "$(date -Iseconds) guard stopped: run DONE" >> "$LOG"
