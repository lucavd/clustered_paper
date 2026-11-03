#!/usr/bin/env bash
# S7 pipeline: extras -> registry -> G3 -> render (docx, html) -> trace render -> G2, G4, P1/P2, Word format (Arial, black, 11 pt body).
# Run from anywhere: bash R/s7/s7_pipeline.sh   (logs in reports/morris/S7_logs/)
set -u
cd /home/user/MMM
export OMP_NUM_THREADS=1
L=reports/morris/S7_logs; mkdir -p $L
Rscript --vanilla R/s7/s7_extras.R > $L/S7_extras.log 2>&1 || { echo "extras failed"; tail -20 $L/S7_extras.log; exit 1; }
Rscript --vanilla R/s7/s7_numbers.R > $L/S7_registry.log 2>&1 || { echo "registry failed"; tail -20 $L/S7_registry.log; exit 1; }
python3 -B tools/s7_verify_registry.py > $L/S7_G3.log 2>&1; echo "G3 exit $?" >> $L/S7_G3.log
cd manuscript/v4
for f in MMM_v4 MMM_v4_supplement; do
  quarto render $f.qmd --to docx > /dev/null 2>&1 || echo "render docx $f failed"
  quarto render $f.qmd --to html > /dev/null 2>&1 || echo "render html $f failed"
  quarto render $f.qmd --to gfm --output ${f}_render.md > /dev/null 2>&1 || echo "render gfm $f failed"
  MMM_S7_TRACE=1 quarto render $f.qmd --to gfm --output ${f}_trace.md > /dev/null 2>&1 || echo "trace $f failed"
done
cd /home/user/MMM
python3 -B tools/s7_check_numbers.py manuscript/v4/MMM_v4_trace.md manuscript/v4/MMM_v4_supplement_trace.md > $L/S7_G2.log 2>&1; echo "G2 exit $?" >> $L/S7_G2.log
python3 -B tools/s7_check_claims.py manuscript/v4/MMM_v4.qmd > $L/S7_G4.log 2>&1; echo "G4 exit $?" >> $L/S7_G4.log
python3 -B tools/s7_wordcount.py manuscript/v4/MMM_v4_render.md > $L/S7_P1P2.log 2>&1
python3 -B tools/s7_check_docx_format.py manuscript/v4/MMM_v4.docx manuscript/v4/MMM_v4_supplement.docx > $L/S7_format.log 2>&1; echo "format exit $?" >> $L/S7_format.log
tail -3 $L/S7_G3.log; head -3 $L/S7_G2.log; grep -E '^G2|exit' $L/S7_G2.log; head -1 $L/S7_G4.log; tail -1 $L/S7_G4.log; cat $L/S7_P1P2.log; cat $L/S7_format.log
