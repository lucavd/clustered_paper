#!/usr/bin/env bash
# S8b pipeline (D-79; copy of R/s7g/s7g_pipeline.sh): Table 4 -> head-to-head (D-74) -> extras -> clinical tables with the
# latent-variable MMM refit (D-70) -> S8b keys -> sub-agent keys -> figures (F01, F02) -> AGQ and table keys -> spec (S8b) ->
# registry -> verification -> render (docx, html, gfm, trace) -> numbers, claims, word count, Word format.
# Prerequisites run once (outputs on disk): R/s8b/s8b_chain_dect1.sh (D-73), R/s8b/s8b_clin_mmm.R, R/s8b/s8b_lus_lateral.R,
# R/s8b/s8b_pb2_collapsed.R. Run: bash R/s8b/s8b_pipeline.sh   (logs in reports/morris/S8b_logs/pipeline/)
set -u
cd /home/user/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=reports/morris/S8b_logs/pipeline; F=reports/morris/S7f_figures; mkdir -p $L $F manuscript/v5/figures
run() { if "$@"; then :; else echo "FAILED: $*"; exit 1; fi; }
run Rscript --vanilla R/s7f/s7f_cpflags.R > $L/cpflags.log 2>&1
run Rscript --vanilla R/s7f/s7f_table4.R > $L/table4.log 2>&1
run Rscript --vanilla R/s7f/s7f_h2h_boot.R > $L/h2h_boot.log 2>&1
run Rscript --vanilla R/s7f/s7f_extras.R > $L/extras.log 2>&1
run Rscript --vanilla R/s7f/s7f_fig1_series.R > $L/fig1_series.log 2>&1
run Rscript --vanilla R/s8b/s8b_clinical_tables.R > $L/clinical_tables.log 2>&1
run Rscript --vanilla R/s7g/s7g_extras.R > $L/s7g_extras.log 2>&1
run Rscript --vanilla R/s8b/s8b_keys.R > $L/s8b_keys.log 2>&1
if [ -f py/s8b/s8b_mcse_captions.py ] && [ -z "${SKIP_SUB:-}" ]; then run python3 -B py/s8b/s8b_mcse_captions.py > $L/mcse_captions.log 2>&1; fi
if [ -f R/s8b/sub/s8b_sub_keys.R ] && [ -z "${SKIP_SUB:-}" ]; then run Rscript --vanilla R/s8b/sub/s8b_sub_keys.R > $L/sub_keys.log 2>&1; fi
cp py/s7f/fig1_rule_v5.json $F/fig1_rule.json
run python3 -B py/s7f/s7f_figures.py tables/variant_kcfc $F tools/s6_figure.mplstyle tables/V5_S7f_fig1c.csv > $L/figures.log 2>&1
cp $F/S6_fig1.png manuscript/v5/figures/fig1.png; cp $F/S6_nestedloop.png manuscript/v5/figures/fig2.png
if [ -f py/s8b/s8b_figS1.py ] && [ -z "${SKIP_SUB:-}" ]; then run python3 -B py/s8b/s8b_figS1.py > $L/figS1.log 2>&1; cp reports/morris/S8b_figures/figS1.png manuscript/v5/figures/figS1.png; fi
run Rscript --vanilla R/s7g/s7g_agq_keys.R > $L/agq_keys.log 2>&1
run Rscript --vanilla R/s7g/s7g_tabkeys.R > $L/tabkeys.log 2>&1
run python3 -B py/s8b_build_spec.py > $L/spec.log 2>&1
run Rscript --vanilla R/s7f/s7f_numbers.R > $L/registry.log 2>&1
python3 -B tools/s7g_verify_registry.py > $L/GB2_verify.log 2>&1; echo "verify exit $?" >> $L/GB2_verify.log
Rscript --vanilla -e '.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); .libPaths(c("/home/user/R/s7d_lib", .libPaths())); suppressMessages({library(fewclust); library(clubSandwich); library(MMLB); library(lme4); library(geepack); library(geesmv)}); writeLines(capture.output(sessionInfo()), "reports/morris/S8b_logs/sessionInfo.txt")' > $L/sessioninfo.log 2>&1
cp manuscript/v4/reference.docx manuscript/v5/reference.docx
cd manuscript/v5
for f in MMM_v5 MMM_v5_supplement; do
  quarto render $f.qmd --to docx > /dev/null 2>&1 || echo "render docx $f failed"
  quarto render $f.qmd --to html > /dev/null 2>&1 || echo "render html $f failed"
  quarto render $f.qmd --to gfm --output ${f}_render.md > /dev/null 2>&1 || echo "render gfm $f failed"
  MMM_S7_TRACE=1 quarto render $f.qmd --to gfm --output ${f}_trace.md > /dev/null 2>&1 || echo "trace $f failed"
done
cd /home/user/MMM
python3 -B tools/s8b_check_numbers.py manuscript/v5/MMM_v5_trace.md manuscript/v5/MMM_v5_supplement_trace.md > $L/GF2.log 2>&1; echo "GF2 exit $?" >> $L/GF2.log
python3 -B tools/s8b_check_claims.py manuscript/v5/MMM_v5_supplement.qmd > $L/GF4_supp.log 2>&1; echo "GF4 supp exit $?" >> $L/GF4_supp.log; cp reports/morris/S8b_logs/S8b_GF4_claims.csv reports/morris/S8b_logs/S8b_GF4_claims_supp.csv
python3 -B tools/s8b_check_claims.py manuscript/v5/MMM_v5.qmd > $L/GF4.log 2>&1; echo "GF4 exit $?" >> $L/GF4.log
python3 -B tools/s7_wordcount.py manuscript/v5/MMM_v5_render.md > $L/wordcount.log 2>&1
python3 -B tools/s7_check_docx_format.py manuscript/v5/MMM_v5.docx manuscript/v5/MMM_v5_supplement.docx > $L/format.log 2>&1; echo "format exit $?" >> $L/format.log
tail -3 $L/GB2_verify.log; grep -E '^GF2|exit' $L/GF2.log; head -1 $L/GF4.log; tail -1 $L/GF4.log; head -1 $L/GF4_supp.log; cat $L/wordcount.log; tail -2 $L/format.log
