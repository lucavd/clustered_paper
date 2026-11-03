# S8b sub-task: snippets for the supplement and the main text

Tokens `r n("f:s9|<key>")` refer to keys of tables/V5_S8b_sub_keys.csv (prefix f:s9| as requested); other tokens are existing registry ids. Part 2 repeats every snippet with values in square brackets, for checking only.

## Part 1. Snippets with tokens

### Table S36 (DEC-T4), new table at the end of the supplement

*Source tables/variant_kcfc/V4_perf_primary.csv via R/s8b/sub/s8b_table_s36.R. The text of Section 2.5 can point here instead of the root mean squared error (P43).*

**Table S36.** Bias in empirical standard errors (Bias/SE), model-based over empirical standard error minus one (SE ratio, %) and coverage of the 95% confidence interval against each method's own estimand (%), under the alternative, on each method's usable replicates, main grid (`r n("x:ndp_grid")` design points) and M1c (`r n("x:ndp_m1c")` design points), marginalized model with the convergence check. Each value is the mean over the design points with at least `r n("x:min_obs_dp")` usable replicates: at least `r n("f:s9|s36|nused_min|main")` of the `r n("f:s9|s36|ndp|main")` design points of a mechanism in the main grid, all `r n("f:s9|s36|ndp|M1c")` in M1c. Not available (n.a.) for the GLMM on GenBinaryY, which has no marginal estimand. Monte Carlo standard errors are at most `r n("f:s9|s36|mcse_max|rb")` for Bias/SE, `r n("f:s9|s36|mcse_max|relerr")` points for the SE ratio and `r n("f:s9|s36|mcse_max|cov")` points for coverage. These are upper bounds: the mean over design points of the per-point Monte Carlo standard errors.

| Method | Grid | Bias/SE GB | Bias/SE NAT | Bias/SE AR1 | SE ratio GB | SE ratio NAT | SE ratio AR1 | Coverage GB | Coverage NAT | Coverage AR1 |
|---|---|---|---|---|---|---|---|---|---|---|
| GEE, Fay-Graubard | Main | `r n("f:s9|s36|gee_fg|GB|main|rb")` | `r n("f:s9|s36|gee_fg|NAT|main|rb")` | `r n("f:s9|s36|gee_fg|AR1|main|rb")` | `r n("f:s9|s36|gee_fg|GB|main|relerr")` | `r n("f:s9|s36|gee_fg|NAT|main|relerr")` | `r n("f:s9|s36|gee_fg|AR1|main|relerr")` | `r n("f:s9|s36|gee_fg|GB|main|cov")` | `r n("f:s9|s36|gee_fg|NAT|main|cov")` | `r n("f:s9|s36|gee_fg|AR1|main|cov")` |
|  | M1c | `r n("f:s9|s36|gee_fg|GB|M1c|rb")` | `r n("f:s9|s36|gee_fg|NAT|M1c|rb")` | `r n("f:s9|s36|gee_fg|AR1|M1c|rb")` | `r n("f:s9|s36|gee_fg|GB|M1c|relerr")` | `r n("f:s9|s36|gee_fg|NAT|M1c|relerr")` | `r n("f:s9|s36|gee_fg|AR1|M1c|relerr")` | `r n("f:s9|s36|gee_fg|GB|M1c|cov")` | `r n("f:s9|s36|gee_fg|NAT|M1c|cov")` | `r n("f:s9|s36|gee_fg|AR1|M1c|cov")` |
| GEE, Mancl-DeRouen | Main | `r n("f:s9|s36|gee_md|GB|main|rb")` | `r n("f:s9|s36|gee_md|NAT|main|rb")` | `r n("f:s9|s36|gee_md|AR1|main|rb")` | `r n("f:s9|s36|gee_md|GB|main|relerr")` | `r n("f:s9|s36|gee_md|NAT|main|relerr")` | `r n("f:s9|s36|gee_md|AR1|main|relerr")` | `r n("f:s9|s36|gee_md|GB|main|cov")` | `r n("f:s9|s36|gee_md|NAT|main|cov")` | `r n("f:s9|s36|gee_md|AR1|main|cov")` |
|  | M1c | `r n("f:s9|s36|gee_md|GB|M1c|rb")` | `r n("f:s9|s36|gee_md|NAT|M1c|rb")` | `r n("f:s9|s36|gee_md|AR1|M1c|rb")` | `r n("f:s9|s36|gee_md|GB|M1c|relerr")` | `r n("f:s9|s36|gee_md|NAT|M1c|relerr")` | `r n("f:s9|s36|gee_md|AR1|M1c|relerr")` | `r n("f:s9|s36|gee_md|GB|M1c|cov")` | `r n("f:s9|s36|gee_md|NAT|M1c|cov")` | `r n("f:s9|s36|gee_md|AR1|M1c|cov")` |
| Marginalized model | Main | `r n("f:s9|s36|mmm|GB|main|rb")` | `r n("f:s9|s36|mmm|NAT|main|rb")` | `r n("f:s9|s36|mmm|AR1|main|rb")` | `r n("f:s9|s36|mmm|GB|main|relerr")` | `r n("f:s9|s36|mmm|NAT|main|relerr")` | `r n("f:s9|s36|mmm|AR1|main|relerr")` | `r n("f:s9|s36|mmm|GB|main|cov")` | `r n("f:s9|s36|mmm|NAT|main|cov")` | `r n("f:s9|s36|mmm|AR1|main|cov")` |
|  | M1c | `r n("f:s9|s36|mmm|GB|M1c|rb")` | `r n("f:s9|s36|mmm|NAT|M1c|rb")` | `r n("f:s9|s36|mmm|AR1|M1c|rb")` | `r n("f:s9|s36|mmm|GB|M1c|relerr")` | `r n("f:s9|s36|mmm|NAT|M1c|relerr")` | `r n("f:s9|s36|mmm|AR1|M1c|relerr")` | `r n("f:s9|s36|mmm|GB|M1c|cov")` | `r n("f:s9|s36|mmm|NAT|M1c|cov")` | `r n("f:s9|s36|mmm|AR1|M1c|cov")` |
| GLMM, Wald | Main | n.a. | `r n("f:s9|s36|glmer|NAT|main|rb")` | `r n("f:s9|s36|glmer|AR1|main|rb")` | `r n("f:s9|s36|glmer|GB|main|relerr")` | `r n("f:s9|s36|glmer|NAT|main|relerr")` | `r n("f:s9|s36|glmer|AR1|main|relerr")` | n.a. | `r n("f:s9|s36|glmer|NAT|main|cov")` | `r n("f:s9|s36|glmer|AR1|main|cov")` |
|  | M1c | n.a. | `r n("f:s9|s36|glmer|NAT|M1c|rb")` | `r n("f:s9|s36|glmer|AR1|M1c|rb")` | `r n("f:s9|s36|glmer|GB|M1c|relerr")` | `r n("f:s9|s36|glmer|NAT|M1c|relerr")` | `r n("f:s9|s36|glmer|AR1|M1c|relerr")` | n.a. | `r n("f:s9|s36|glmer|NAT|M1c|cov")` | `r n("f:s9|s36|glmer|AR1|M1c|cov")` |
| Shrinkage GLM | Main | `r n("f:s9|s36|shrink|GB|main|rb")` | `r n("f:s9|s36|shrink|NAT|main|rb")` | `r n("f:s9|s36|shrink|AR1|main|rb")` | `r n("f:s9|s36|shrink|GB|main|relerr")` | `r n("f:s9|s36|shrink|NAT|main|relerr")` | `r n("f:s9|s36|shrink|AR1|main|relerr")` | `r n("f:s9|s36|shrink|GB|main|cov")` | `r n("f:s9|s36|shrink|NAT|main|cov")` | `r n("f:s9|s36|shrink|AR1|main|cov")` |
|  | M1c | `r n("f:s9|s36|shrink|GB|M1c|rb")` | `r n("f:s9|s36|shrink|NAT|M1c|rb")` | `r n("f:s9|s36|shrink|AR1|M1c|rb")` | `r n("f:s9|s36|shrink|GB|M1c|relerr")` | `r n("f:s9|s36|shrink|NAT|M1c|relerr")` | `r n("f:s9|s36|shrink|AR1|M1c|relerr")` | `r n("f:s9|s36|shrink|GB|M1c|cov")` | `r n("f:s9|s36|shrink|NAT|M1c|cov")` | `r n("f:s9|s36|shrink|AR1|M1c|cov")` |
| GEE, uncorrected | Main | `r n("f:s9|s36|gee|GB|main|rb")` | `r n("f:s9|s36|gee|NAT|main|rb")` | `r n("f:s9|s36|gee|AR1|main|rb")` | `r n("f:s9|s36|gee|GB|main|relerr")` | `r n("f:s9|s36|gee|NAT|main|relerr")` | `r n("f:s9|s36|gee|AR1|main|relerr")` | `r n("f:s9|s36|gee|GB|main|cov")` | `r n("f:s9|s36|gee|NAT|main|cov")` | `r n("f:s9|s36|gee|AR1|main|cov")` |
|  | M1c | `r n("f:s9|s36|gee|GB|M1c|rb")` | `r n("f:s9|s36|gee|NAT|M1c|rb")` | `r n("f:s9|s36|gee|AR1|M1c|rb")` | `r n("f:s9|s36|gee|GB|M1c|relerr")` | `r n("f:s9|s36|gee|NAT|M1c|relerr")` | `r n("f:s9|s36|gee|AR1|M1c|relerr")` | `r n("f:s9|s36|gee|GB|M1c|cov")` | `r n("f:s9|s36|gee|NAT|M1c|cov")` | `r n("f:s9|s36|gee|AR1|M1c|cov")` |
| GEE, Kauermann-Carroll | Main | `r n("f:s9|s36|gee_kc|GB|main|rb")` | `r n("f:s9|s36|gee_kc|NAT|main|rb")` | `r n("f:s9|s36|gee_kc|AR1|main|rb")` | `r n("f:s9|s36|gee_kc|GB|main|relerr")` | `r n("f:s9|s36|gee_kc|NAT|main|relerr")` | `r n("f:s9|s36|gee_kc|AR1|main|relerr")` | `r n("f:s9|s36|gee_kc|GB|main|cov")` | `r n("f:s9|s36|gee_kc|NAT|main|cov")` | `r n("f:s9|s36|gee_kc|AR1|main|cov")` |
|  | M1c | `r n("f:s9|s36|gee_kc|GB|M1c|rb")` | `r n("f:s9|s36|gee_kc|NAT|M1c|rb")` | `r n("f:s9|s36|gee_kc|AR1|M1c|rb")` | `r n("f:s9|s36|gee_kc|GB|M1c|relerr")` | `r n("f:s9|s36|gee_kc|NAT|M1c|relerr")` | `r n("f:s9|s36|gee_kc|AR1|M1c|relerr")` | `r n("f:s9|s36|gee_kc|GB|M1c|cov")` | `r n("f:s9|s36|gee_kc|NAT|M1c|cov")` | `r n("f:s9|s36|gee_kc|AR1|M1c|cov")` |

### MCSE caption sentences (DEC-T5)

*One sentence (or two) to append to each caption. Kind of each value in column how of tables/V5_S8b_sub_keys.csv: registry maximum (S3-S12, S16-S25, S30), upper bound for means over design points (Table 4, S14 usable, S26, S28, S33, S34, S35 means), binomial (Figure 2, S14 largest rate, S15, S27, S32, S35 smallest coverage), block bootstrap (S31), delta-method approximation (S35 bias and SE ratio).*

- Table 4 (main text): Monte Carlo standard errors are at most `r n("f:s9|cap|t4|pu")` points for the usable share, `r n("f:s9|cap|t4|upow")` for power and `r n("f:s9|cap|t4|pu_small")` for the usable share at T ≤ 3 and N ≤ 25. These are upper bounds: the mean over design points of the binomial standard error. Counts of design points rest on the exact interval of each design point.

- Figure 2 (main text): The Monte Carlo standard error of a rate of 5% is `r n("x:mcse_5pct")` points with `r n("x:nsim")` usable replicates and at most `r n("f:s9|cap|fig2|mcse_max")` points at the drawn design point with the fewest.

- Table S3: Largest Monte Carlo standard error: `r n("f:s9|cap|S3|rate")` points for the pooled rate, `r n("f:s9|cap|S3|worst")` for the worst case and `r n("f:s9|cap|S3|mad")` for the MAD.

- Table S4: Largest Monte Carlo standard error: `r n("f:s9|cap|S4|mad")` points.

- Table S5: Largest Monte Carlo standard error: `r n("f:s9|cap|S5|share")` points.

- Table S6: Largest Monte Carlo standard error: `r n("f:s9|cap|S6|rate")` points for the rates and `r n("f:s9|cap|S6|shift")` for the shift.

- Table S7: Largest Monte Carlo standard error: `r n("f:s9|cap|S7|usable")` points.

- Table S8: Largest Monte Carlo standard error: `r n("f:s9|cap|S8|usable")` points for the usable share and `r n("f:s9|cap|S8|uncond")` for unconditional validity.

- Table S9: Largest Monte Carlo standard error: `r n("f:s9|cap|S9|mad")` points.

- Table S10: Largest Monte Carlo standard error: `r n("f:s9|cap|S10|raw")` points for raw power and `r n("f:s9|cap|S10|sa")` for size-adjusted power.

- Table S11: Largest Monte Carlo standard error: `r n("f:s9|cap|S11|rate")` points for the rates and `r n("f:s9|cap|S11|shift")` for the differences.

- Table S12: Largest Monte Carlo standard error: `r n("f:s9|cap|S12|cov")` points.

- Table S14: Monte Carlo standard errors are at most `r n("f:s9|cap|S14|maxrate")` points for the largest rate (binomial, at its design point) and `r n("f:s9|cap|S14|pu")` for the usable share (upper bound: mean over design points of the binomial standard error).

- Table S15: Binomial Monte Carlo standard errors over the `r n("x:lus_nsub")` subsamples are at most `r n("f:s9|cap|S15|share")` points.

- Table S16: Largest Monte Carlo standard error: `r n("f:s9|cap|S16|rate")` points for the pooled rate, `r n("f:s9|cap|S16|worst")` for the worst case and `r n("f:s9|cap|S16|mad")` for the MAD.

- Table S17: Largest Monte Carlo standard error: `r n("f:s9|cap|S17|mad")` points.

- Table S18: Largest Monte Carlo standard error: `r n("f:s9|cap|S18|share")` points.

- Table S19: Largest Monte Carlo standard error: `r n("f:s9|cap|S19|rate")` points for the rates and `r n("f:s9|cap|S19|shift")` for the shift.

- Table S20: Largest Monte Carlo standard error: `r n("f:s9|cap|S20|usable")` points.

- Table S21: Largest Monte Carlo standard error: `r n("f:s9|cap|S21|usable")` points for the usable share and `r n("f:s9|cap|S21|uncond")` for unconditional validity.

- Table S22: Largest Monte Carlo standard error: `r n("f:s9|cap|S22|mad")` points.

- Table S23: Largest Monte Carlo standard error: `r n("f:s9|cap|S23|raw")` points for raw power and `r n("f:s9|cap|S23|sa")` for size-adjusted power.

- Table S24: Largest Monte Carlo standard error: `r n("f:s9|cap|S24|rate")` points for the rates and `r n("f:s9|cap|S24|shift")` for the differences.

- Table S25: Largest Monte Carlo standard error: `r n("f:s9|cap|S25|cov")` points.

- Table S26: The Monte Carlo standard error of the usable share is at most `r n("f:s9|cap|S26|pu")` points (upper bound: mean over design points of the binomial standard error).

- Table S27: Binomial Monte Carlo standard errors are at most `r n("f:s9|cap|S27|finite")` points for finite intervals and `r n("f:s9|cap|S27|gee")` for GEE; they are approximate, because the datasets of a row share replicate seeds across numbers of clusters.

- Table S28: Monte Carlo standard errors are at most `r n("f:s9|cap|S28|type1")` points for the type I error, `r n("f:s9|cap|S28|usable")` for the usable share, `r n("f:s9|cap|S28|power")` for power, `r n("f:s9|cap|S28|biasse")` for bias in empirical standard errors and `r n("f:s9|cap|S28|cov")` points for coverage. These are upper bounds: the mean over design points of the per-point standard error.

- Table S30: Largest Monte Carlo standard error: `r n("f:s9|cap|S30|usable")` points for the usable share, `r n("f:s9|cap|S30|bias")` for the bias and `r n("f:s9|cap|S30|seratio")` for the SE ratio.

- Table S31: Block-bootstrap Monte Carlo standard errors are at most `r n("f:s9|cap|S31|rate")` points for the type I error, `r n("f:s9|cap|S31|mad")` for the departure from 5%, `r n("f:s9|cap|S31|pu")` for the usable share and `r n("f:s9|cap|S31|upow")` for power.

- Table S32: The binomial Monte Carlo standard error of the highest rate is at most `r n("f:s9|cap|S32|maxrate")` points.

- Table S33: The Monte Carlo standard error of the mean usable share is at most `r n("f:s9|cap|S33|pu")` points (upper bound: mean over design points of the binomial standard error).

- Table S34: Monte Carlo standard errors are at most `r n("f:s9|cap|S34|type1")` points for the mean type I error and `r n("f:s9|cap|S34|pu")` for the usable share. These are upper bounds: the mean over design points of the binomial standard error.

- Table S35: Monte Carlo standard errors, with the empirical standard error treated as known, are at most `r n("f:s9|cap|S35|rb")` for the mean bias in standard errors, `r n("f:s9|cap|S35|rbmax")` for the largest, `r n("f:s9|cap|S35|cov")` points for the mean coverage, `r n("f:s9|cap|S35|covmin")` for the smallest coverage and `r n("f:s9|cap|S35|sh")` for the share of fits in which the check intervened. For the SE ratio the typical per-point value is `r n("f:s9|cap|S35|ser")` points. The columns by intervention of the check have no Monte Carlo standard error.

### Table S37 and paragraph for Section S11 (X02)

*Source tables/V5_S8b_rule_stability.csv via R/s8b/sub/s8b_table_s37.R; the calculations are those of R/s8b/provenance/lit_pow_rule_stability_v2.R.*

**Table S37.** Stability of the count of design points clearly above 7.5% over the `r n("f:x|ndp_all5")` design points of the four mechanisms, strategies of Table 4 with an interval. Observed: the count in the study (Table 4). Bootstrap: `r n("f:s9|s37|B")` block-bootstrap resamples of the replicates, the same resample for every design point of a model; share of resamples with no design point clearly above 7.5% (%), mean and 90th percentile of the count. Plug-in: design points treated as independent, with true rate equal to the observed rate (raw) or to the rate shrunk by empirical Bayes within strategy and mechanism (shrunk); probability of no design point clearly above 7.5% (%) and expected count.

| Strategy | Observed | Bootstrap, none (%) | Bootstrap, mean | Bootstrap, 90th percentile | Raw plug-in, none (%) | Raw plug-in, expected | Shrunk plug-in, none (%) | Shrunk plug-in, expected |
|---|---|---|---|---|---|---|---|---|
| GEE AR(1), Fay-Graubard | `r n("f:s9|s37|gee_ar1_fg|hi_obs")` | `r n("f:s9|s37|gee_ar1_fg|p_pass_boot|disp")` | `r n("f:s9|s37|gee_ar1_fg|mean_flags_boot")` | `r n("f:s9|s37|gee_ar1_fg|q90_flags_boot")` | `r n("f:s9|s37|gee_ar1_fg|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_ar1_fg|plugin_exp_flags")` | `r n("f:s9|s37|gee_ar1_fg|eb_p_pass|disp")` | `r n("f:s9|s37|gee_ar1_fg|eb_exp_flags")` |
| Marginalized model with check | `r n("f:s9|s37|mmm_cc|hi_obs")` | `r n("f:s9|s37|mmm_cc|p_pass_boot|disp")` | `r n("f:s9|s37|mmm_cc|mean_flags_boot")` | `r n("f:s9|s37|mmm_cc|q90_flags_boot")` | `r n("f:s9|s37|mmm_cc|plugin_p_pass|disp")` | `r n("f:s9|s37|mmm_cc|plugin_exp_flags")` | `r n("f:s9|s37|mmm_cc|eb_p_pass|disp")` | `r n("f:s9|s37|mmm_cc|eb_exp_flags")` |
| GEE independence, Fay-Graubard | `r n("f:s9|s37|gee_ind_fg|hi_obs")` | `r n("f:s9|s37|gee_ind_fg|p_pass_boot|disp")` | `r n("f:s9|s37|gee_ind_fg|mean_flags_boot")` | `r n("f:s9|s37|gee_ind_fg|q90_flags_boot")` | `r n("f:s9|s37|gee_ind_fg|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_ind_fg|plugin_exp_flags")` | `r n("f:s9|s37|gee_ind_fg|eb_p_pass|disp")` | `r n("f:s9|s37|gee_ind_fg|eb_exp_flags")` |
| GEE exchangeable, Fay-Graubard | `r n("f:s9|s37|gee_exch_fg|hi_obs")` | `r n("f:s9|s37|gee_exch_fg|p_pass_boot|disp")` | `r n("f:s9|s37|gee_exch_fg|mean_flags_boot")` | `r n("f:s9|s37|gee_exch_fg|q90_flags_boot")` | `r n("f:s9|s37|gee_exch_fg|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_exch_fg|plugin_exp_flags")` | `r n("f:s9|s37|gee_exch_fg|eb_p_pass|disp")` | `r n("f:s9|s37|gee_exch_fg|eb_exp_flags")` |
| GEE AR(1), Mancl-DeRouen | `r n("f:s9|s37|gee_ar1_md|hi_obs")` | `r n("f:s9|s37|gee_ar1_md|p_pass_boot|disp")` | `r n("f:s9|s37|gee_ar1_md|mean_flags_boot")` | `r n("f:s9|s37|gee_ar1_md|q90_flags_boot")` | `r n("f:s9|s37|gee_ar1_md|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_ar1_md|plugin_exp_flags")` | `r n("f:s9|s37|gee_ar1_md|eb_p_pass|disp")` | `r n("f:s9|s37|gee_ar1_md|eb_exp_flags")` |
| GEE exchangeable, Mancl-DeRouen | `r n("f:s9|s37|gee_exch_md|hi_obs")` | `r n("f:s9|s37|gee_exch_md|p_pass_boot|disp")` | `r n("f:s9|s37|gee_exch_md|mean_flags_boot")` | `r n("f:s9|s37|gee_exch_md|q90_flags_boot")` | `r n("f:s9|s37|gee_exch_md|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_exch_md|plugin_exp_flags")` | `r n("f:s9|s37|gee_exch_md|eb_p_pass|disp")` | `r n("f:s9|s37|gee_exch_md|eb_exp_flags")` |
| GEE AR(1), Kauermann-Carroll | `r n("f:s9|s37|gee_ar1_kc|hi_obs")` | `r n("f:s9|s37|gee_ar1_kc|p_pass_boot|disp")` | `r n("f:s9|s37|gee_ar1_kc|mean_flags_boot")` | `r n("f:s9|s37|gee_ar1_kc|q90_flags_boot")` | `r n("f:s9|s37|gee_ar1_kc|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_ar1_kc|plugin_exp_flags")` | `r n("f:s9|s37|gee_ar1_kc|eb_p_pass|disp")` | `r n("f:s9|s37|gee_ar1_kc|eb_exp_flags")` |
| GEE exchangeable, Kauermann-Carroll | `r n("f:s9|s37|gee_exch_kc|hi_obs")` | `r n("f:s9|s37|gee_exch_kc|p_pass_boot|disp")` | `r n("f:s9|s37|gee_exch_kc|mean_flags_boot")` | `r n("f:s9|s37|gee_exch_kc|q90_flags_boot")` | `r n("f:s9|s37|gee_exch_kc|plugin_p_pass|disp")` | `r n("f:s9|s37|gee_exch_kc|plugin_exp_flags")` | `r n("f:s9|s37|gee_exch_kc|eb_p_pass|disp")` | `r n("f:s9|s37|gee_exch_kc|eb_exp_flags")` |

Paragraph for S11, after Table S32:

The zero count of the two strategies of the rule is uncertain in a replication of the study. In the block bootstrap, GEE with AR(1) working correlation and Fay-Graubard had at least one design point clearly above 7.5% in `r n("f:s9|s37|fg_pfail_boot")`% of resamples, and the marginalized model with the check in `r n("f:s9|s37|mmm_pfail_boot")`% (Table S37). The ranking by expected count is the same in all three calculations, and in each of them the two strategies of the rule have the fewest expected design points clearly above 7.5%.

### Rows for Table S26 and sentence for S6 (X03)

*Rows go after the row of GEE exchangeable, Kauermann-Carroll. Source R/s8b/sub/s8b_mdkc_s26.R (1,896 points, N - p and N - 2); the text sentence uses tables/V5_S8b_mdkc_avg.csv (2,376 points).*

Version A (recommended; N - 2 computed with the same code, gate passed):

| GEE AR(1), average of Mancl-DeRouen and Kauermann-Carroll | `r n("f:s9|mdkc|gee_ar1_avgvar|N - p|n_hi")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - p|n_lo")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - p|n_lo_N10")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - 2|n_hi")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - 2|n_lo")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - 2|n_lo_N10")` | `r n("f:s9|mdkc|gee_ar1_avgvar|pu")` |
| GEE exchangeable, average of Mancl-DeRouen and Kauermann-Carroll | `r n("f:s9|mdkc|gee_exch_avgvar|N - p|n_hi")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - p|n_lo")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - p|n_lo_N10")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - 2|n_hi")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - 2|n_lo")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - 2|n_lo_N10")` | `r n("f:s9|mdkc|gee_exch_avgvar|pu")` |

Caption sentence for version A: The rows for the average of the Mancl-DeRouen and Kauermann-Carroll variances use the same usability rule, denominators and exact intervals, recomputed from the stored per-replicate results.

Version B (as first instructed; N - 2 columns not filled):

| GEE AR(1), average of Mancl-DeRouen and Kauermann-Carroll | `r n("f:s9|mdkc|gee_ar1_avgvar|N - p|n_hi")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - p|n_lo")` | `r n("f:s9|mdkc|gee_ar1_avgvar|N - p|n_lo_N10")` | n.a. | n.a. | n.a. | `r n("f:s9|mdkc|gee_ar1_avgvar|pu")` |
| GEE exchangeable, average of Mancl-DeRouen and Kauermann-Carroll | `r n("f:s9|mdkc|gee_exch_avgvar|N - p|n_hi")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - p|n_lo")` | `r n("f:s9|mdkc|gee_exch_avgvar|N - p|n_lo_N10")` | n.a. | n.a. | n.a. | `r n("f:s9|mdkc|gee_exch_avgvar|pu")` |

Caption sentence for version B: For the average of the Mancl-DeRouen and Kauermann-Carroll variances the N - 2 columns were not computed (n.a.).

Sentence for the text of S6 (or Section 4.5 of the main text):

Averaging the Mancl-DeRouen and Kauermann-Carroll variances, as studied by Ford and Westgate [14], did not give a strategy that was never clearly above 7.5%. With N - p degrees of freedom, GEE with AR(1) working correlation was clearly above 7.5% at `r n("f:s9|mdkc|gee_ar1_avgvar|hi_all")` of the `r n("f:x|ndp_all5")` design points of the four mechanisms, and GEE with exchangeable working correlation at `r n("f:s9|mdkc|gee_exch_avgvar|hi_all")` (Table S26 gives the three study mechanisms).

### Figure S1 (F02)

*Figure file reports/morris/S8b_figures/figS1.png (copy to manuscript/v5/figures/figS1.png when the supplement is edited). Script py/s8b/s8b_figS1.py.*

![](figures/figS1.png){width=100%}

**Figure S1.** Type I error at every design point of scenarios M6-M7 (`r n("x:n_dp_ns_scen")` design points each) for the strategies of Figure 1, with the corrections from `fewclust` as in Table 4: Wald tests with a t reference (N - p degrees of freedom) and the likelihood-ratio test with a $\chi^2$ reference. Rows: GenBinaryY with latent-SD slope `r n("x:m6_gby")` and `r n("x:m7_gby")`; GLMM-native with random-slope variance `r n("x:m6_nat")` and `r n("x:m7_nat")`. Design points are ordered by number of occasions T (blocks) and number of clusters N. Grey band: 2.5% to 7.5%; line: nominal 5%. Design points with fewer than `r n("x:min_obs_dp")` usable replicates are not drawn (`r n("f:s9|figS1|glmer_lrt|n_hidden")` design points of the likelihood-ratio test, none of the other strategies). Circles: design points clearly above 7.5% (lower limit of the exact 95% interval above 7.5%), as in Table 4: `r n("f:s9|figS1|gee_exch_md|n_hi")` for GEE exchangeable with Mancl-DeRouen, `r n("f:s9|figS1|gee_exch_fg|n_hi")` for GEE exchangeable with Fay-Graubard, all on GenBinaryY; `r n("f:s9|figS1|gee_ar1_fg|n_hi")` for GEE AR(1) with Fay-Graubard, `r n("f:s9|figS1|mmm_cc|n_hi")` for the marginalized model with the check and `r n("f:s9|figS1|glmer_lrt|n_hi")` for the likelihood-ratio test. The binomial Monte Carlo standard error of a drawn rate is at most `r n("f:s9|figS1|mcse_max")` points.

## Part 2. Preview with values (3 significant digits, not the registry rounding)

### Table S36 (DEC-T4), new table at the end of the supplement

**Table S36.** Bias in empirical standard errors (Bias/SE), model-based over empirical standard error minus one (SE ratio, %) and coverage of the 95% confidence interval against each method's own estimand (%), under the alternative, on each method's usable replicates, main grid ([1,440] design points) and M1c ([72] design points), marginalized model with the convergence check. Each value is the mean over the design points with at least [20] usable replicates: at least [416] of the [480] design points of a mechanism in the main grid, all [24] in M1c. Not available (n.a.) for the GLMM on GenBinaryY, which has no marginal estimand. Monte Carlo standard errors are at most [0.0364] for Bias/SE, [5.86] points for the SE ratio and [1.04] points for coverage. These are upper bounds: the mean over design points of the per-point Monte Carlo standard errors.

| Method | Grid | Bias/SE GB | Bias/SE NAT | Bias/SE AR1 | SE ratio GB | SE ratio NAT | SE ratio AR1 | Coverage GB | Coverage NAT | Coverage AR1 |
|---|---|---|---|---|---|---|---|---|---|---|
| GEE, Fay-Graubard | Main | [0.0319] | [0.0233] | [0.0222] | [7.38] | [8.28] | [8.55] | [94.6] | [95.1] | [95.1] |
|  | M1c | [0.0571] | [0.0379] | [0.0515] | [7.01] | [8.52] | [7.74] | [94] | [94.7] | [94.6] |
| GEE, Mancl-DeRouen | Main | [0.032] | [0.0239] | [0.0223] | [6.26] | [7.57] | [8.02] | [94.3] | [94.8] | [94.8] |
|  | M1c | [0.0568] | [0.0379] | [0.0521] | [14.2] | [15.6] | [15.2] | [94.8] | [95.4] | [95.2] |
| Marginalized model | Main | [0.0291] | [0.0219] | [0.0196] | [2.28] | [1.89] | [1.39] | [94.5] | [94.8] | [94.7] |
|  | M1c | [0.0584] | [0.0311] | [0.0426] | [7.39] | [2.18] | [1.01] | [93.1] | [93.2] | [93.2] |
| GLMM, Wald | Main | n.a. | [0.0369] | [-0.0179] | [4.35] | [5.87] | [6.1] | n.a. | [93.6] | [93.3] |
|  | M1c | n.a. | [0.0158] | [0.00346] | [7.81] | [8.47] | [5.32] | n.a. | [95.4] | [95] |
| Shrinkage GLM | Main | [0.124] | [0.0794] | [0.0754] | [2.13] | [4] | [3.87] | [89.9] | [92.7] | [92.5] |
|  | M1c | [0.0383] | [0.0309] | [0.044] | [-4.54] | [2.65] | [0.72] | [88.9] | [92.6] | [92] |
| GEE, uncorrected | Main | [0.0301] | [0.0213] | [0.0201] | [-7.75] | [-6.98] | [-6.75] | [91.5] | [92.3] | [92.2] |
|  | M1c | [0.0559] | [0.0358] | [0.0506] | [-8.67] | [-8.14] | [-8.66] | [90] | [91.1] | [91.1] |
| GEE, Kauermann-Carroll | Main | [0.0304] | [0.0224] | [0.0205] | [-0.77] | [0.155] | [0.427] | [93.3] | [93.9] | [93.8] |
|  | M1c | [0.0595] | [0.0388] | [0.0514] | [0.614] | [2.09] | [1.76] | [92.9] | [93.7] | [93.6] |

### MCSE caption sentences (DEC-T5)

- Table 4 (main text): Monte Carlo standard errors are at most [0.401] points for the usable share, [1.03] for power and [1.42] for the usable share at T ≤ 3 and N ≤ 25. These are upper bounds: the mean over design points of the binomial standard error. Counts of design points rest on the exact interval of each design point.

- Figure 2 (main text): The Monte Carlo standard error of a rate of 5% is [0.69] points with [1,000] usable replicates and at most [6.36] points at the drawn design point with the fewest.

- Table S3: Largest Monte Carlo standard error: [0.0472] points for the pooled rate, [0.303] for the worst case and [0.0471] for the MAD.

- Table S4: Largest Monte Carlo standard error: [0.497] points.

- Table S5: Largest Monte Carlo standard error: [1.11] points.

- Table S6: Largest Monte Carlo standard error: [0.0472] points for the rates and [0.0204] for the shift.

- Table S7: Largest Monte Carlo standard error: [0.126] points.

- Table S8: Largest Monte Carlo standard error: [0.0963] points for the usable share and [0.104] for unconditional validity.

- Table S9: Largest Monte Carlo standard error: [1.35] points.

- Table S10: Largest Monte Carlo standard error: [0.0632] points for raw power and [0.264] for size-adjusted power.

- Table S11: Largest Monte Carlo standard error: [0.0472] points for the rates and [0.035] for the differences.

- Table S12: Largest Monte Carlo standard error: [0.049] points.

- Table S14: Monte Carlo standard errors are at most [1.28] points for the largest rate (binomial, at its design point) and [0.552] for the usable share (upper bound: mean over design points of the binomial standard error).

- Table S15: Binomial Monte Carlo standard errors over the [500] subsamples are at most [2.24] points.

- Table S16: Largest Monte Carlo standard error: [0.216] points for the pooled rate, [0.651] for the worst case and [0.216] for the MAD.

- Table S17: Largest Monte Carlo standard error: [0.284] points.

- Table S18: Largest Monte Carlo standard error: [5.02] points.

- Table S19: Largest Monte Carlo standard error: [0.216] points for the rates and [0.132] for the shift.

- Table S20: Largest Monte Carlo standard error: [0.35] points.

- Table S21: Largest Monte Carlo standard error: [0.26] points for the usable share and [0.277] for unconditional validity.

- Table S22: Largest Monte Carlo standard error: [0.909] points.

- Table S23: Largest Monte Carlo standard error: [0.236] points for raw power and [0.285] for size-adjusted power.

- Table S24: Largest Monte Carlo standard error: [0.216] points for the rates and [0.104] for the differences.

- Table S25: Largest Monte Carlo standard error: [0.203] points.

- Table S26: The Monte Carlo standard error of the usable share is at most [0.397] points (upper bound: mean over design points of the binomial standard error).

- Table S27: Binomial Monte Carlo standard errors are at most [0.323] points for finite intervals and [0.302] for GEE; they are approximate, because the datasets of a row share replicate seeds across numbers of clusters.

- Table S28: Monte Carlo standard errors are at most [1.1] points for the type I error, [0.415] for the usable share, [1.05] for power, [0.0359] for bias in empirical standard errors and [0.896] points for coverage. These are upper bounds: the mean over design points of the per-point standard error.

- Table S30: Largest Monte Carlo standard error: [2.24] points for the usable share, [0.215] for the bias and [0.152] for the SE ratio.

- Table S31: Block-bootstrap Monte Carlo standard errors are at most [0.201] points for the type I error, [0.174] for the departure from 5%, [0.299] for the usable share and [0.0692] for power.

- Table S32: The binomial Monte Carlo standard error of the highest rate is at most [1.66] points.

- Table S33: The Monte Carlo standard error of the mean usable share is at most [0.401] points (upper bound: mean over design points of the binomial standard error).

- Table S34: Monte Carlo standard errors are at most [3.09] points for the mean type I error and [1.01] for the usable share. These are upper bounds: the mean over design points of the binomial standard error.

- Table S35: Monte Carlo standard errors, with the empirical standard error treated as known, are at most [0.0361] for the mean bias in standard errors, [0.192] for the largest, [0.74] points for the mean coverage, [0.892] for the smallest coverage and [1.71] for the share of fits in which the check intervened. For the SE ratio the typical per-point value is [2.19] points. The columns by intervention of the check have no Monte Carlo standard error.

### Table S37 and paragraph for Section S11 (X02)

**Table S37.** Stability of the count of design points clearly above 7.5% over the [2,376] design points of the four mechanisms, strategies of Table 4 with an interval. Observed: the count in the study (Table 4). Bootstrap: [2000] block-bootstrap resamples of the replicates, the same resample for every design point of a model; share of resamples with no design point clearly above 7.5% (%), mean and 90th percentile of the count. Plug-in: design points treated as independent, with true rate equal to the observed rate (raw) or to the rate shrunk by empirical Bayes within strategy and mechanism (shrunk); probability of no design point clearly above 7.5% (%) and expected count.

| Strategy | Observed | Bootstrap, none (%) | Bootstrap, mean | Bootstrap, 90th percentile | Raw plug-in, none (%) | Raw plug-in, expected | Shrunk plug-in, none (%) | Shrunk plug-in, expected |
|---|---|---|---|---|---|---|---|---|
| GEE AR(1), Fay-Graubard | [0] | [38] | [0.949] | [2] | [37] | [0.94] | [93] | [0.0754] |
| Marginalized model with check | [0] | [18] | [1.7] | [3] | [16] | [1.72] | [80] | [0.219] |
| GEE independence, Fay-Graubard | [2] | [1.7] | [3.39] | [5] | [1.4] | [3.35] | [45] | [0.745] |
| GEE exchangeable, Fay-Graubard | [11] | [< 0.1] | [13.3] | [17] | [< 0.1] | [13.3] | [0.9] | [4.16] |
| GEE AR(1), Mancl-DeRouen | [16] | [< 0.1] | [24.7] | [29] | [< 0.1] | [24.8] | [< 0.1] | [9.8] |
| GEE exchangeable, Mancl-DeRouen | [57] | [< 0.1] | [69.6] | [76] | [< 0.1] | [69.6] | [< 0.1] | [46.7] |
| GEE AR(1), Kauermann-Carroll | [132] | [< 0.1] | [136] | [143] | [< 0.1] | [136] | [< 0.1] | [104] |
| GEE exchangeable, Kauermann-Carroll | [199] | [< 0.1] | [202] | [211] | [< 0.1] | [202] | [< 0.1] | [175] |

Paragraph for S11, after Table S32:

The zero count of the two strategies of the rule is uncertain in a replication of the study. In the block bootstrap, GEE with AR(1) working correlation and Fay-Graubard had at least one design point clearly above 7.5% in [62.1]% of resamples, and the marginalized model with the check in [81.8]% (Table S37). The ranking by expected count is the same in all three calculations, and in each of them the two strategies of the rule have the fewest expected design points clearly above 7.5%.

### Rows for Table S26 and sentence for S6 (X03)

Version A (recommended; N - 2 computed with the same code, gate passed):

| GEE AR(1), average of Mancl-DeRouen and Kauermann-Carroll | [37] | [83] | [3] | [62] | [80] | [3] | [92] |
| GEE exchangeable, average of Mancl-DeRouen and Kauermann-Carroll | [74] | [83] | [3] | [104] | [80] | [3] | [92.1] |

Caption sentence for version A: The rows for the average of the Mancl-DeRouen and Kauermann-Carroll variances use the same usability rule, denominators and exact intervals, recomputed from the stored per-replicate results.

Version B (as first instructed; N - 2 columns not filled):

| GEE AR(1), average of Mancl-DeRouen and Kauermann-Carroll | [37] | [83] | [3] | n.a. | n.a. | n.a. | [92] |
| GEE exchangeable, average of Mancl-DeRouen and Kauermann-Carroll | [74] | [83] | [3] | n.a. | n.a. | n.a. | [92.1] |

Caption sentence for version B: For the average of the Mancl-DeRouen and Kauermann-Carroll variances the N - 2 columns were not computed (n.a.).

Sentence for the text of S6 (or Section 4.5 of the main text):

Averaging the Mancl-DeRouen and Kauermann-Carroll variances, as studied by Ford and Westgate [14], did not give a strategy that was never clearly above 7.5%. With N - p degrees of freedom, GEE with AR(1) working correlation was clearly above 7.5% at [51] of the [2,376] design points of the four mechanisms, and GEE with exchangeable working correlation at [113] (Table S26 gives the three study mechanisms).

### Figure S1 (F02)

![](figures/figS1.png){width=100%}

**Figure S1.** Type I error at every design point of scenarios M6-M7 ([96] design points each) for the strategies of Figure 1, with the corrections from `fewclust` as in Table 4: Wald tests with a t reference (N - p degrees of freedom) and the likelihood-ratio test with a $\chi^2$ reference. Rows: GenBinaryY with latent-SD slope [0.1] and [0.25]; GLMM-native with random-slope variance [0.01] and [0.04]. Design points are ordered by number of occasions T (blocks) and number of clusters N. Grey band: 2.5% to 7.5%; line: nominal 5%. Design points with fewer than [20] usable replicates are not drawn ([60] design points of the likelihood-ratio test, none of the other strategies). Circles: design points clearly above 7.5% (lower limit of the exact 95% interval above 7.5%), as in Table 4: [17] for GEE exchangeable with Mancl-DeRouen, [4] for GEE exchangeable with Fay-Graubard, all on GenBinaryY; [0] for GEE AR(1) with Fay-Graubard, [0] for the marginalized model with the check and [0] for the likelihood-ratio test. The binomial Monte Carlo standard error of a drawn rate is at most [3.39] points.
