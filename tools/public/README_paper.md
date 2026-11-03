# Choosing an analysis for clustered binary outcomes when clusters are few and series are short

Code and results of the simulation study and of the two clinical applications in the paper
"Choosing an analysis for clustered binary outcomes when clusters are few and series are short:
a simulation study and two clinical applications" (Vedovelli et al., submitted to *Statistical Methods
in Medical Research*).

This repository is a public copy of the authors' working repository, updated in parallel with it.
The R package written for the study, *fewclust*, is at <https://github.com/lucavd/fewclust>.

## Contents

- `R/`, `py/`, `tools/`: simulation, analysis, aggregation and manuscript pipeline. The closing pipeline of the
  paper is `R/s11/s11_pipeline.sh`; every number of the manuscript is read from a file by the registry
  (`tables/V5_S7f_numbers.csv`, built by `R/s7f/s7f_numbers.R` and checked by `tools/s7g_verify_registry.py`).
- `.rundir/`: the code of the original simulation runs.
- `tables/`: per-design-point results and every table behind the paper; `tables/variant_cv/` is the variant
  used in the paper.
- `manuscript/v5/`: manuscript and supplement sources (Quarto) and figures.
- `reports/morris/`, `docs/`: step-by-step reports and project documentation (in Italian), following
  Morris, White and Crowther (2019).
- `renv.lock`: R package versions. *geesmv* 1.3 was archived on CRAN on 12 October 2026; install it with
  `remotes::install_version("geesmv", "1.3")`.

## Not included

- The neonatal cohort data are not publicly available (see the data availability statement of the paper);
  the scripts that read them are included, and the plasmode results are in `tables/` and in the supplement.
- The per-replicate model fits are not in this repository because of their size; they will be archived separately on Zenodo (DOI to be added).

## Licence

MIT (see `LICENSE`).
