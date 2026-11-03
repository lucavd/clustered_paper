#!/bin/bash
cd ~/MMM
echo "=== 1. decision_guide: ricerca esaustiva ==="
echo "-- file di lavoro (esclusi renv/.git/results_full) --"
grep -rl 'decision_guide' ~/MMM --exclude-dir=renv --exclude-dir=.git --exclude-dir=results_full 2>/dev/null
echo "-- scratch dei job --"
grep -rl 'decision_guide' ~/.claude-science-scratch 2>/dev/null | head -20
echo "-- qualunque .R nel filesystem utente --"
find ~ -name '*.R' -not -path '*/renv/*' 2>/dev/null | xargs grep -l 'decision_guide' 2>/dev/null | head
echo "-- storia git: file mai aggiunti con 'decision' nel nome --"
git log --all --diff-filter=A --name-only --pretty=format: 2>/dev/null | sort -u | grep -i decision
echo "(vuoto = nessun file 'decision' mai committato)"
echo "-- git grep su TUTTI i commit --"
for cm in $(git rev-list --all 2>/dev/null); do
  if git grep -q 'decision_guide' $cm -- '*.R' 2>/dev/null; then echo "TROVATO in $cm"; fi
done | head -5
echo "(nessun output = la stringa non e mai stata in un .R versionato)"
echo "-- file fig_decision_guide: mtime e provenienza --"
ls -l --time-style=+%F_%T ~/MMM/ARTICOLO*/FINAL/fig_decision_guide.* 2>/dev/null
echo "-- il PDF dichiara il produttore? --"
strings ~/MMM/ARTICOLO*/FINAL/fig_decision_guide.pdf 2>/dev/null | grep -iE 'creator|producer|cairo|ggplot|R version' | head -5
echo "-- MMM-package? --"
ls -d ~/MMM-package 2>&1 | head -2

echo ""
echo "=== 2. DFREF_cells.csv: quale script? ==="
echo "-- header --"; head -1 ~/MMM/rev2/DFREF_cells.csv
echo "-- riga N=10,T=4,model=1,mmm,GenBinaryY --"
awk -F, 'NR==1 || ($1=="GenBinaryY" && $2=="mmm" && $3==10 && $4==4 && $5==1)' ~/MMM/rev2/DFREF_cells.csv
echo "-- output effettivi di R3_df.R --"
ls -l ~/MMM/rev2/R3_df_bymethod.csv ~/MMM/rev2/R3_df_worstcase.csv ~/MMM/rev2/R3_df_N25.csv 2>&1
echo "-- header di R3_df_worstcase.csv --"; head -2 ~/MMM/rev2/R3_df_worstcase.csv 2>/dev/null
echo "-- chi nomina DFREF/MMMROB/DGRID in QUALUNQUE .R --"
find ~ -name '*.R' -not -path '*/renv/*' 2>/dev/null | xargs grep -ln 'DFREF_cells\|MMMROB_raw\|DGRID_publishable\|DGRID_raw' 2>/dev/null | head -20

echo ""
echo "=== 3. valori della riesecuzione R3-R* esistono in un file? ==="
for v in 0.5205 0.5636 0.5740 0.5989 0.6351; do
  printf "%s -> " $v
  grep -rl "$v" ~/MMM --exclude-dir=renv --exclude-dir=.git --exclude-dir=results_full 2>/dev/null | tr '\n' ' '
  echo
done
echo "-- file subsample --"
find ~ -iname '*subsample*' 2>/dev/null | head
echo "-- R3_mmm.R: seme e denominatore --"
grep -n 'set.seed\|500\|n_sub\|ratio\|sd_est\|se_' ~/MMM/rev2/R3_mmm.R 2>/dev/null | head -25

echo ""
echo "=== 4. F_paired_mmm_fg.csv esiste? ==="
find ~ -name 'F_paired*' 2>/dev/null | head
echo "=== post_gan.R nel repo? ==="
ls -l ~/MMM/R/postrun/post_gan.R 2>&1 | head -2
ls -1 ~/MMM/R/postrun/ 2>/dev/null | head -30
echo "DONE_JOB3"
