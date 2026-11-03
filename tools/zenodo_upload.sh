#!/usr/bin/env bash
# zenodo_upload.sh -- per-replicate results of the paper to the Zenodo deposition (draft); token header in ~/.zenodo_hdr (deleted at the end)
set -uo pipefail
BUCKET="$1"; Z=/mnt/micron/zenodo_staging; mkdir -p $Z; cd /home/user/MMM
trap 'rm -f /home/user/.zenodo_hdr' EXIT
LIST_V4="results_v4/S4/rows_v4.rds results_v4/S4/truth_ar1_rep.rds results_v4/S4/boot_idx.rds results_v4/S6/mmm_boundary_all.rds results_v4/S6cc/cc_all.rds results_v4/S6cc/rows results_v4/S6b/ns/rows_ns.rds results_v4/S6b/ns/boot_idx_67.rds results_v4/S7d/rows results_v4/S7e/rows results_v4/S10/degen/classes.rds results_v4/S11/cv_keys.rds"
t0=$(date +%s)
tar -cf $Z/MMM_results_full.tar results_full || exit 1
tar -cf $Z/MMM_results_v4_per_replicate.tar $LIST_V4 || exit 1
echo "tar done $(( $(date +%s) - t0 )) s"
{ echo "Per-replicate results of the simulation study (deposit 1.0)."
  echo "Code: https://doi.org/10.5281/zenodo.23917381 (lucavd/clustered_paper release 1.0); private source commit $(git rev-parse --short HEAD)."
  echo "Extract both archives in the root of the code repository: paths are relative to it (results_full/, results_v4/)."
  echo "Seeds: 20260805 + 1000*model + i, L'Ecuyer-CMRG; label mmlb = alternative, manual = null."
  echo; echo "Archives (md5):"; (cd $Z && md5sum MMM_results_full.tar MMM_results_v4_per_replicate.tar)
  echo; echo "Contents of MMM_results_v4_per_replicate.tar:"; for p in $LIST_V4; do if [ -d "$p" ]; then echo "  $p/ ($(ls $p | wc -l) files)"; else echo "  $p"; fi; done
  echo; echo "results_full/: $(ls results_full | wc -l) files; md5 of every file in baseline/results_full_md5.txt of the code repository." ; } > $Z/README_deposit.txt
cat $Z/README_deposit.txt | head -12
for f in README_deposit.txt MMM_results_v4_per_replicate.tar MMM_results_full.tar; do
  t1=$(date +%s)
  curl -sS -f -H @/home/user/.zenodo_hdr --upload-file $Z/$f "$BUCKET/$f" -o $Z/$f.resp.json || { echo "upload FAIL $f"; exit 1; }
  echo "uploaded $f $(( $(date +%s) - t1 )) s: $(python3 -c "import json;j=json.load(open('$Z/$f.resp.json'));print(j.get('checksum'), j.get('size'))")"
done
echo "local md5:"; (cd $Z && md5sum README_deposit.txt MMM_results_v4_per_replicate.tar MMM_results_full.tar)
