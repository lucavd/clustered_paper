#!/usr/bin/env bash
# S7f pipeline: Table 4 -> extras -> clinical KC tables -> figures -> spec -> registry -> GF3(a) -> render (docx, html, gfm, trace)
# -> GF2, GF4, PF1/PF2, Word format. Run: bash R/s7f/s7f_pipeline.sh   (logs in reports/morris/S7f_logs/)
set -u
cd /home/user/MMM
export OMP_NUM_THREADS=1
L=reports/morris/S7f_logs; F=reports/morris/S7f_figures; mkdir -p $L $F manuscript/v5/figures
run() { if "$@"; then :; else echo "FAILED: $*"; exit 1; fi; }
run Rscript --vanilla R/s7f/s7f_cpflags.R > $L/S7f_cpflags.log 2>&1
run Rscript --vanilla R/s7f/s7f_table4.R > $L/S7f_table4.log 2>&1
run Rscript --vanilla R/s7f/s7f_h2h_boot.R > $L/S7f_h2h_boot.log 2>&1
run Rscript --vanilla R/s7f/s7f_extras.R > $L/S7f_extras.log 2>&1
run Rscript --vanilla R/s7f/s7f_fig1_series.R > $L/S7f_fig1_series.log 2>&1
run Rscript --vanilla R/s7f/s7f_kc_clinical_tables.R > $L/S7f_kc_clinical_tables.log 2>&1
cp py/s7f/fig1_rule_v5.json $F/fig1_rule.json
run python3 -B py/s7f/s7f_figures.py tables/variant_kcfc $F tools/s6_figure.mplstyle tables/V5_S7f_fig1c.csv > $L/S7f_figures.log 2>&1
cp $F/S6_fig1.png manuscript/v5/figures/fig1.png; cp $F/S6_nestedloop.png manuscript/v5/figures/fig2.png
run python3 -B py/s7f_build_spec.py > $L/S7f_spec.log 2>&1
run Rscript --vanilla R/s7f/s7f_numbers.R > $L/S7f_registry.log 2>&1
python3 -B tools/s7f_verify_registry.py > $L/S7f_GF3a.log 2>&1; echo "GF3a exit $?" >> $L/S7f_GF3a.log
cp manuscript/v4/reference.docx manuscript/v5/reference.docx
cd manuscript/v5
for f in MMM_v5 MMM_v5_supplement; do
  quarto render $f.qmd --to docx > /dev/null 2>&1 || echo "render docx $f failed"
  quarto render $f.qmd --to html > /dev/null 2>&1 || echo "render html $f failed"
  quarto render $f.qmd --to gfm --output ${f}_render.md > /dev/null 2>&1 || echo "render gfm $f failed"
  MMM_S7_TRACE=1 quarto render $f.qmd --to gfm --output ${f}_trace.md > /dev/null 2>&1 || echo "trace $f failed"
done
cd /home/user/MMM
python3 -B tools/s7f_check_numbers.py manuscript/v5/MMM_v5_trace.md manuscript/v5/MMM_v5_supplement_trace.md > $L/S7f_GF2.log 2>&1; echo "GF2 exit $?" >> $L/S7f_GF2.log
python3 -B tools/s7f_check_claims.py manuscript/v5/MMM_v5.qmd > $L/S7f_GF4.log 2>&1; echo "GF4 exit $?" >> $L/S7f_GF4.log
python3 -B tools/s7_wordcount.py manuscript/v5/MMM_v5_render.md > $L/S7f_PF1PF2.log 2>&1
python3 -B tools/s7_check_docx_format.py manuscript/v5/MMM_v5.docx manuscript/v5/MMM_v5_supplement.docx > $L/S7f_format.log 2>&1; echo "format exit $?" >> $L/S7f_format.log
tail -3 $L/S7f_GF3a.log; grep -E '^GF2|exit' $L/S7f_GF2.log; head -1 $L/S7f_GF4.log; tail -1 $L/S7f_GF4.log; cat $L/S7f_PF1PF2.log; tail -2 $L/S7f_format.log
