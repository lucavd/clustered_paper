#!/usr/bin/env bash
# s9_pipeline.sh -- S9 closing pipeline (D-82, D-83): copy of R/s8b/s8b_pipeline.sh with
#   * the table variant nc (constant-outcome datasets not usable) instead of kcfc;
#   * the scripts that S8b ran once as prerequisites and that read per-replicate GEE rows, rerun in dependency order:
#     sensitivity to the usability limits (R/s7g/s7g_filter.R), average of the MD and KC variances and stability of the
#     rule (R/s8b/provenance/, run in results_v4/S9/prov and copied to tables/ only if the columns are those of the
#     stored file), Tables S36, S37 and the S26 rows (R/s8b/sub/);
#   * a first spec + registry build before the sub-tables, so that their registry checks read the nc paths.
# Prerequisites (outputs on disk): R/s9/s9_nc_pipeline.sh, R/s9/s9_upstream.sh. Run: bash R/s9/s9_pipeline.sh [tag]
# Logs in reports/morris/S9_logs/pipeline_<tag>/ (default tag: run).
set -u
cd /home/user/MMM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TAG=${1:-run}
L=/home/user/MMM/reports/morris/S9_logs/pipeline_$TAG; F=reports/morris/S7f_figures; P=results_v4/S9/prov; mkdir -p $L $F $P manuscript/v5/figures
SKIP=${START:+1}   # START=<step name>: resume from that step (earlier steps are skipped)
run() { b=$1; shift; if [ -n "$SKIP" ]; then if [ "$b" = "$START" ]; then SKIP=; else echo "skip $b"; return 0; fi; fi; t0=$(date +%s); if "$@" > $L/$b.log 2>&1; then echo "OK   $b $(( $(date +%s) - t0 )) s"; else echo "FAIL $b (see $L/$b.log)"; tail -15 $L/$b.log; exit 1; fi; }
run s7_extras    Rscript --vanilla R/s7/s7_extras.R
run cpflags      Rscript --vanilla R/s7f/s7f_cpflags.R
run table4       Rscript --vanilla R/s7f/s7f_table4.R
run s9_keys      Rscript --vanilla R/s9/s9_keys.R
run h2h_boot     Rscript --vanilla R/s7f/s7f_h2h_boot.R
run s7g_filter   Rscript --vanilla R/s7g/s7g_filter.R
( cd $P && run du8_mdkc_avg Rscript --vanilla /home/user/MMM/R/s8b/provenance/du8_mdkc_avg.R ) || exit 1
( cd $P && run lit_pow      Rscript --vanilla /home/user/MMM/R/s8b/provenance/lit_pow_rule_stability_v2.R ) || exit 1
cphdr() { if [ -n "$SKIP" ]; then echo "skip copy $1"; return 0; fi; if [ "$(head -1 $1)" = "$(head -1 $2)" ]; then cp $1 $2; echo "OK   copy $1 -> $2"; else echo "FAIL header of $1 differs from $2"; exit 1; fi; }
cphdr $P/du8_table4_mdkc_avg.csv tables/V5_S8b_mdkc_avg.csv
cphdr $P/lit_pow1_pass_v2.csv tables/V5_S8b_rule_stability.csv
( cd $P && run disc_r2     Rscript --vanilla /home/user/MMM/R/s8b/provenance/disc_r2.R /home/user/MMM/results_v4/S9/prov ) || exit 1
cphdr $P/disc_r2_flip_summary.csv tables/V5_S8b_flip_interval.csv
cphdr $P/disc_r2_flip_lo_points.csv tables/V5_S8b_flip_lo_points.csv
run extras       Rscript --vanilla R/s7f/s7f_extras.R
run fig1_series  Rscript --vanilla R/s7f/s7f_fig1_series.R
run clinical     Rscript --vanilla R/s8b/s8b_clinical_tables.R
run s7g_extras   Rscript --vanilla R/s7g/s7g_extras.R
run agq_keys0    Rscript --vanilla R/s7g/s7g_agq_keys.R
run tabkeys0     Rscript --vanilla R/s7g/s7g_tabkeys.R
run spec0        python3 -B py/s8b_build_spec.py
run registry0    Rscript --vanilla R/s7f/s7f_numbers.R
run s36          Rscript --vanilla R/s8b/sub/s8b_table_s36.R
run s37          Rscript --vanilla R/s8b/sub/s8b_table_s37.R
run mdkc_s26     Rscript --vanilla R/s8b/sub/s8b_mdkc_s26.R
run s8b_keys     Rscript --vanilla R/s8b/s8b_keys.R
run mcse_captions python3 -B py/s8b/s8b_mcse_captions.py
run sub_keys     Rscript --vanilla R/s8b/sub/s8b_sub_keys.R
cp py/s7f/fig1_rule_v5.json $F/fig1_rule.json
run figures      python3 -B py/s7f/s7f_figures.py tables/variant_nc $F tools/s6_figure.mplstyle tables/V5_S7f_fig1c.csv
cp $F/S6_fig1.png manuscript/v5/figures/fig1.png; cp $F/S6_nestedloop.png manuscript/v5/figures/fig2.png
run figS1        python3 -B py/s8b/s8b_figS1.py
cp reports/morris/S8b_figures/figS1.png manuscript/v5/figures/figS1.png
run agq_keys     Rscript --vanilla R/s7g/s7g_agq_keys.R
run tabkeys      Rscript --vanilla R/s7g/s7g_tabkeys.R
run spec         python3 -B py/s8b_build_spec.py
run registry     Rscript --vanilla R/s7f/s7f_numbers.R
python3 -B tools/s7g_verify_registry.py > $L/GB2_verify.log 2>&1; echo "verify exit $?" >> $L/GB2_verify.log
Rscript --vanilla -e '.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); .libPaths(c("/home/user/R/s7d_lib", .libPaths())); suppressMessages({library(fewclust); library(clubSandwich); library(MMLB); library(lme4); library(geepack); library(geesmv)}); writeLines(capture.output(sessionInfo()), "reports/morris/S9_logs/sessionInfo.txt")' > $L/sessioninfo.log 2>&1
cp manuscript/v4/reference.docx manuscript/v5/reference.docx
python3 R/s9/s9_ref_style.py manuscript/v5/reference.docx > /dev/null   # D-88: Title 18 pt left; Heading 1-3 16/14/12 pt
cd manuscript/v5
for f in MMM_v5 MMM_v5_supplement; do
  quarto render $f.qmd --to docx > /dev/null 2>&1 || echo "render docx $f failed"
  quarto render $f.qmd --to html > /dev/null 2>&1 || echo "render html $f failed"
  quarto render $f.qmd --to gfm --output ${f}_render.md > /dev/null 2>&1 || echo "render gfm $f failed"
  MMM_S7_TRACE=1 quarto render $f.qmd --to gfm --output ${f}_trace.md > /dev/null 2>&1 || echo "trace $f failed"
done
cd /home/user/MMM
python3 -B tools/s8b_check_numbers.py manuscript/v5/MMM_v5_trace.md manuscript/v5/MMM_v5_supplement_trace.md > $L/GF2.log 2>&1; echo "GF2 exit $?" >> $L/GF2.log
python3 -B tools/s8b_check_claims.py manuscript/v5/MMM_v5_supplement.qmd > $L/GF4_supp.log 2>&1; echo "GF4 supp exit $?" >> $L/GF4_supp.log
python3 -B tools/s8b_check_claims.py manuscript/v5/MMM_v5.qmd > $L/GF4.log 2>&1; echo "GF4 exit $?" >> $L/GF4.log
python3 -B tools/s7_wordcount.py manuscript/v5/MMM_v5_render.md > $L/wordcount.log 2>&1
python3 -B tools/s7_check_docx_format.py manuscript/v5/MMM_v5.docx manuscript/v5/MMM_v5_supplement.docx > $L/format.log 2>&1; echo "format exit $?" >> $L/format.log
tail -3 $L/GB2_verify.log; grep -E '^GF2|exit' $L/GF2.log; head -1 $L/GF4.log; tail -1 $L/GF4.log; head -1 $L/GF4_supp.log; tail -1 $L/GF4_supp.log; cat $L/wordcount.log; tail -2 $L/format.log
echo "pipeline $TAG done"
