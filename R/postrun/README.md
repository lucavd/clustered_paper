# Post-run analysis scripts

Every table in `tables/` is produced by a script in this directory. They were
previously run from per-dispatch scratch directories and existed nowhere else,
so eight tracked tables had no reproducible source and the scratch directories
are swept periodically. That is the same defect this project already hit once,
when the repaired metrics filter lived only in an unversioned script while the
defective version stayed in `R/`.

| script | produces | notes |
|---|---|---|
| `fix_m5_gee.R` | `results_m5gee/` (5,760 chunks) | Refits the uncorrected GEE for scenario 5 only, after the dispatcher was found routing it to scenario 4's fitter. Reuses the production seed `20260805 + 1000*model + i`, so the same datasets are regenerated and every other method's stored rows stay valid. |
| `truth_ar1_conditional.R` | `tables/F_truth_ar1.csv` | Recovers the conditional truth the AR(1) generator calibrates to, per (model, N, T). NOTE: scenarios 3-5 use `~ time * trt` with NO covariate in the mean model; the covariate enters only through the random effects. Two earlier drafts got this wrong and silently produced results for models 1-2 only. |
| `truth_glmmnative_marginal.R` | `tables/F_truth_glmmnative.csv` | Exact Gauss-Hermite marginalisation of the conditional model, giving the marginal effect the Zeger-approximation conditional value actually implies. It is not the flat 0.200 assumed, and it drifts with series length. |
| `reagg_fixed.R` | `tables/G_level_fixed.csv`, `G_coverage_fixed.csv` | Re-aggregates the stored chunks with the scenario-5 refit swapped in, a method-appropriate exclusion filter, both denominators (own and intersection), and the per-cell truths joined. |
| `fix_tables.R` | `tables/F_coverage.csv`, `F_paired_mmm_fg.csv`, `F_worstcase_*.csv` | KNOWN STALE: built before `results_m5gee/` existed, so it still carries the defective scenario-5 GEE rows. Needs the swap from `reagg_fixed.R` before its outputs can be quoted. |
| `verify_plim.R` | (diagnostic, no table) | Tests whether the mixed model coverage trend is an artefact of scoring against a value the estimator does not converge to. It is: at ten occasions the trend reverses from -4.1 to +2.0 points. At two occasions the collapse is real and target-independent. |

All of these read the stored chunks; none re-runs the simulation, except
`fix_m5_gee.R`, which refits one arm of one scenario.

Run them with an explicit library path and `--vanilla`, since `.Rprofile`
activates renv and overwrites `R_LIBS`:

```bash
Rscript --vanilla R/postrun/reagg_fixed.R
```
