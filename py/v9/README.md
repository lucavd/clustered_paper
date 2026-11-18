# py/v9: code of the v9 rewrite of the manuscript (2026-11-08/10)

The v9 rewrite (after the SMMR decision of 2026-11-08) added a figure of the four mechanisms (Figure 1), results by
mechanism and region (Figure 2, Figure S1, Figure 3, Tables 3 and 4), a table of method differences (Table 1) and a
clinical table (Table 5). No new simulations: everything starts from the per-design-point tables of the v8 chain.

The scripts are the cells executed in the rewrite session and its two sub-sessions, extracted verbatim and in order
(full record in `provenance/*_cells.json`); only file paths were changed. Run from the repository root:

| Script | Output (tables/v9/) | Check on 2026-11-18 |
|---|---|---|
| 01_build_dp_tidy.py | v8_dp_tidy.csv.gz (45,144 rows = 2,376 design points x 19 strategies) | equal to the file used for v9 |
| 02_results_by_mechanism.py | region_summary_*, mechanism_summary_*, best_*, robust_choice, sensitivity_* | 9 files identical (md5) |
| 03_guidance_by_region.py | guidance_by_region.csv (Table 4) | equal |
| R/v9/mechanisms_generate.R (server, 4 cores, about 2 min) | mech_*.csv | mech_pair_corr.csv identical (md5) |
| 05_mechanisms_figure_data.py | mechanisms_lag_correlation.csv (data of Figure 1) | identical (md5) |
| 06_method_differences_table.py | method_differences.csv, .md (Table 1) | identical (md5) |
| 07_clinical_table.py | clinical_table.csv, .md (Table 5) | identical (md5) |
| 04_figures_heatmap_availability.py | Figures 2, S1, 3 and figC_data_byT.csv | not rerun: needs plotting helpers outside this repository |

Inputs copied from `results_v4/` (not tracked): `S11_degen_table5_points_versions.csv.gz` (md5 e465fd6ec041e697477cce7540e5776f),
`S10_tabelle_s4_all_dp.csv` (md5 e8ea4ce14832ad510ef63aed78ae1afb). Manuscripts: `manuscript/v9/` (SMMR version, sources,
figures, and the Statistics in Medicine package in `manuscript/v9/SiM/`).
