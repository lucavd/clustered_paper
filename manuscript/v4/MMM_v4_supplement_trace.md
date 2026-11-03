# Supplementary material: choosing an analysis for clustered binary
outcomes when clusters are few and series are short


This supplement accompanies the manuscript. Every number is read from
the same registry of values as the main text. GB = GenBinaryY, NAT =
GLMM-native, AR1 = GLMM-AR1; MMM = marginalized multilevel model.

## S1. Parameters of the data-generating mechanisms

**Table S1.** Parameters under the alternative, read from the table
extracted from the executed code. Under the null only the interaction
coefficient is set to zero. Coefficients are in the order (Intercept),
time, trt, covariate if any, time:trt; they are marginal for GenBinaryY,
conditional for GLMM-native and the marginal target for GLMM-AR1, whose
conditional coefficients are calibrated on each replicate.

| Mechanism | Scenario | Mean model | Coefficients | Dependence |
|----|----|----|----|----|
| GenBinaryY | M1 | § | § | lv §, σ = (§); transition §, γ = (§) |
| GenBinaryY | M2 | § | § | lv §, σ = (§); transition §, γ = (§) |
| GenBinaryY | M3 | § | § | lv §, σ = (§); transition §, γ = (§) |
| GenBinaryY | M4 | § | § | lv §, σ = (§); transition §, γ = (§) |
| GenBinaryY | M5 | § | § | lv §, σ = (§); transition §, γ = (§) |
| GLMM-native | M1 | § | § | random effects §, D = \[§\] |
| GLMM-native | M2 | § | § | random effects §, D = \[§\] |
| GLMM-native | M3 | § | § | random effects §, D = \[§\] |
| GLMM-native | M4 | § | § | random effects §, D = \[§\] |
| GLMM-native | M5 | § | § | random effects §, D = \[§\] |
| GLMM-AR1 | M1 | § | § | random effects §, D = \[§\]; AR(1) ρ = §, variance § |
| GLMM-AR1 | M2 | § | § | random effects §, D = \[§\]; AR(1) ρ = §, variance § |
| GLMM-AR1 | M3 | § | § | random effects §, D = \[§\]; AR(1) ρ = §, variance § |
| GLMM-AR1 | M4 | § | § | random effects §, D = \[§\]; AR(1) ρ = §, variance § |
| GLMM-AR1 | M5 | § | § | random effects §, D = \[§\]; AR(1) ρ = §, variance § |

## S2. Software

**Table S2.** Versions in the project library used for every rerun and
for the tables of this paper. The original run was executed under an
earlier patch release of R; its datasets were reproduced bit for bit
under the current release. The full `sessionInfo()` is provided with the
code (`results_v4/S7/sessionInfo.txt`).

| Software   | Version |
|------------|---------|
| R          | §       |
| MMLB       | §       |
| lme4       | §       |
| geepack    | §       |
| geesmv     | §       |
| data.table | §       |
| crew       | §       |
| mirai      | §       |
| Matrix     | §       |

## S3. Tables of the main grid

**Table S3.** Type I error at the nominal 5% (%): pooled rate, worst
case and mean absolute departure (MAD) over cluster counts (percentage
points), main grid, with the marginalized model analysed as in the
original study, without the convergence check (sensitivity analysis).
Only the row of the marginalized model differs from Table 2. GB =
GenBinaryY, NAT = GLMM-native, AR1 = GLMM-AR1.

|  | rate GB | rate NAT | rate AR1 | worst GB | worst NAT | worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | § | § | § | § | § | § | § | § | § |
| GEE, Mancl-DeRouen | § | § | § | § | § | § | § | § | § |
| Marginalized model | § | § | § | § | § | § | § | § | § |
| GLMM, Wald | § | § | § | § | § | § | § | § | § |
| Shrinkage GLM | § | § | § | § | § | § | § | § | § |
| GEE, uncorrected | § | § | § | § | § | § | § | § | § |
| GLMM, likelihood ratio | § | § | § | § | § | § | § | § | § |
| GEE, Kauermann-Carroll | § | § | § | § | § | § | § | § | § |

**Table S4.** Mean absolute departure from 5% by number of occasions
(percentage points), cell-weighted, main grid, marginalized model with
the convergence check. Monte Carlo standard errors by block bootstrap;
values given to the digit they support.

|                        | T = 2 | T = 3 | T \>= 4 |
|------------------------|-------|-------|---------|
| GEE, Fay-Graubard      | §     | §     | §       |
| GEE, Mancl-DeRouen     | §     | §     | §       |
| Marginalized model     | §     | §     | §       |
| GLMM, Wald             | §     | §     | §       |
| Shrinkage GLM          | §     | §     | §       |
| GEE, uncorrected       | §     | §     | §       |
| GLMM, likelihood ratio | §     | §     | §       |
| GEE, Kauermann-Carroll | §     | §     | §       |

**Table S5.** Method-cells detectably off nominal (%), exact binomial
test with Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY)
corrections within method, main grid, marginalized model with the
convergence check. Monte Carlo standard errors by block bootstrap;
values given to the digit they support.

|                        | BH  | BY  |
|------------------------|-----|-----|
| GEE, Fay-Graubard      | §   | §   |
| GEE, Mancl-DeRouen     | §   | §   |
| Marginalized model     | §   | §   |
| GLMM, Wald             | §   | §   |
| Shrinkage GLM          | §   | §   |
| GEE, uncorrected       | §   | §   |
| GLMM, likelihood ratio | §   | §   |
| GEE, Kauermann-Carroll | §   | §   |

**Table S6.** Type I error on each method’s own usable replicates and on
the replicates usable for all methods (common set) (%), main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.

|  | own GB | own NAT | own AR1 | common GB | common NAT | common AR1 | shift GB | shift NAT | shift AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | § | § | § | § | § | § | § | § | § |
| GEE, Mancl-DeRouen | § | § | § | § | § | § | § | § | § |
| Marginalized model | § | § | § | § | § | § | § | § | § |
| GLMM, Wald | § | § | § | § | § | § | § | § | § |
| Shrinkage GLM | § | § | § | § | § | § | § | § | § |
| GEE, uncorrected | § | § | § | § | § | § | § | § | § |
| GLMM, likelihood ratio | § | § | § | § | § | § | § | § | § |
| GEE, Kauermann-Carroll | § | § | § | § | § | § | § | § | § |

**Table S7.** Replicates returning a usable result by number of
occasions (%), every attempted replicate in the denominator, main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.

|  | GEE, Fay-Graubard | GEE, Mancl-DeRouen | Marginalized model | GLMM, Wald | Shrinkage GLM | GEE, uncorrected | GLMM, likelihood ratio | GEE, Kauermann-Carroll |
|----|----|----|----|----|----|----|----|----|
| random intercept, T = 2 | § | § | § | § | § | § | § | § |
| random intercept, T = 3 | § | § | § | § | § | § | § | § |
| random intercept, T = 4 | § | § | § | § | § | § | § | § |
| random intercept, T = 6 | § | § | § | § | § | § | § | § |
| random intercept, T = 8 | § | § | § | § | § | § | § | § |
| random intercept, T = 10 | § | § | § | § | § | § | § | § |
| random slopes, T = 2 | § | § | § | § | § | § | § | § |
| random slopes, T = 3 | § | § | § | § | § | § | § | § |
| random slopes, T = 4 | § | § | § | § | § | § | § | § |
| random slopes, T = 6 | § | § | § | § | § | § | § | § |
| random slopes, T = 8 | § | § | § | § | § | § | § | § |
| random slopes, T = 10 | § | § | § | § | § | § | § | § |

**Table S8.** Usable result and unconditional validity at 20 clusters or
fewer (%), main grid, marginalized model with the convergence check.
Monte Carlo standard errors by block bootstrap; values given to the
digit they support.

|  | usable | unconditional validity | design points with a usable replicate |
|----|----|----|----|
| GEE, Fay-Graubard | § | § | § |
| GEE, Mancl-DeRouen | § | § | § |
| Marginalized model | § | § | § |
| GLMM, Wald | § | § | § |
| Shrinkage GLM | § | § | § |
| GEE, uncorrected | § | § | § |
| GLMM, likelihood ratio | § | § | § |
| GEE, Kauermann-Carroll | § | § | § |

**Table S9.** Mean absolute departure from 5% by estimability stratum
(percentage points); strata by the mean usable proportion of the five
marginal methods, main grid, marginalized model with the convergence
check. Monte Carlo standard errors by block bootstrap; values given to
the digit they support.

|                        | §   | §   | §   | §   | §   |
|------------------------|-----|-----|-----|-----|-----|
| GEE, Fay-Graubard      | §   | §   | §   | §   | §   |
| GEE, Mancl-DeRouen     | §   | §   | §   | §   | §   |
| Marginalized model     | §   | §   | §   | §   | §   |
| GLMM, Wald             | §   | §   | §   | §   | §   |
| Shrinkage GLM          | §   | §   | §   | §   | §   |
| GEE, uncorrected       | §   | §   | §   | §   | §   |
| GLMM, likelihood ratio | §   | §   | §   | §   | §   |
| GEE, Kauermann-Carroll | §   | §   | §   | §   | §   |
| design points          | §   | §   | §   | §   | §   |

**Table S10.** Raw and size-adjusted power at the interaction contrast
(%), main grid, marginalized model with the convergence check. Monte
Carlo standard errors by block bootstrap; values given to the digit they
support.

|                                     | raw | size-adjusted |
|-------------------------------------|-----|---------------|
| GenBinaryY, GEE, Fay-Graubard       | §   | §             |
| GenBinaryY, GEE, Mancl-DeRouen      | §   | §             |
| GenBinaryY, Marginalized model      | §   | §             |
| GenBinaryY, GLMM, Wald              | §   | §             |
| GenBinaryY, Shrinkage GLM           | §   | §             |
| GenBinaryY, GEE, uncorrected        | §   | §             |
| GenBinaryY, GLMM, likelihood ratio  | §   | §             |
| GenBinaryY, GEE, Kauermann-Carroll  | §   | §             |
| GLMM-native, GEE, Fay-Graubard      | §   | §             |
| GLMM-native, GEE, Mancl-DeRouen     | §   | §             |
| GLMM-native, Marginalized model     | §   | §             |
| GLMM-native, GLMM, Wald             | §   | §             |
| GLMM-native, Shrinkage GLM          | §   | §             |
| GLMM-native, GEE, uncorrected       | §   | §             |
| GLMM-native, GLMM, likelihood ratio | §   | §             |
| GLMM-native, GEE, Kauermann-Carroll | §   | §             |
| GLMM-AR1, GEE, Fay-Graubard         | §   | §             |
| GLMM-AR1, GEE, Mancl-DeRouen        | §   | §             |
| GLMM-AR1, Marginalized model        | §   | §             |
| GLMM-AR1, GLMM, Wald                | §   | §             |
| GLMM-AR1, Shrinkage GLM             | §   | §             |
| GLMM-AR1, GEE, uncorrected          | §   | §             |
| GLMM-AR1, GLMM, likelihood ratio    | §   | §             |
| GLMM-AR1, GEE, Kauermann-Carroll    | §   | §             |

**Table S11.** Type I error under a normal and a t reference (%), shift
in percentage points, main grid, marginalized model with the convergence
check. Monte Carlo standard errors by block bootstrap; values given to
the digit they support.

|  | normal GB | normal NAT | normal AR1 | t GB | t NAT | t AR1 | shift GB | shift NAT | shift AR1 | shift at N \<= 25 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | § | § | § | § | § | § | § | § | § | § |
| GEE, Mancl-DeRouen | § | § | § | § | § | § | § | § | § | § |
| Marginalized model | § | § | § | § | § | § | § | § | § | § |
| GLMM, Wald | § | § | § | § | § | § | § | § | § | § |
| Shrinkage GLM | § | § | § | § | § | § | § | § | § | § |
| GEE, uncorrected | § | § | § | § | § | § | § | § | § | § |
| GEE, Kauermann-Carroll | § | § | § | § | § | § | § | § | § | § |

**Table S12.** Coverage of the 95% confidence interval against each
method’s own estimand (%); not available (n.a.) for the GLMM on
GenBinaryY, main grid, marginalized model with the convergence check.
Monte Carlo standard errors by block bootstrap; values given to the
digit they support.

|                        | GenBinaryY | GLMM-native | GLMM-AR1 |
|------------------------|------------|-------------|----------|
| GEE, Fay-Graubard      | §          | §           | §        |
| GEE, Mancl-DeRouen     | §          | §           | §        |
| Marginalized model     | §          | §           | §        |
| GLMM, Wald             | n.a.       | §           | §        |
| Shrinkage GLM          | §          | §           | §        |
| GEE, uncorrected       | §          | §           | §        |
| GEE, Kauermann-Carroll | §          | §           | §        |

**Table S13.** Pooled type I error (%) on the main grid with fits that
returned warnings kept (U0, primary) and excluded (U1),
replicate-weighted over the null cells. Binomial Monte Carlo standard
error at most § points.

| Method                         | GB, U0 | GB, U1 | NAT, U0 | NAT, U1 | AR1, U0 | AR1, U1 |
|--------------------------------|--------|--------|---------|---------|---------|---------|
| GEE, Fay-Graubard              | §      | §      | §       | §       | §       | §       |
| GEE, Mancl-DeRouen             | §      | §      | §       | §       | §       | §       |
| Marginalized model, with check | §      | §      | §       | §       | §       | §       |
| GLMM, Wald                     | §      | §      | §       | §       | §       | §       |
| Shrinkage GLM                  | §      | §      | §       | §       | §       | §       |
| GEE, uncorrected               | §      | §      | §       | §       | §       | §       |
| GLMM, likelihood ratio         | §      | §      | §       | §       | §       | §       |
| GEE, Kauermann-Carroll         | §      | §      | §       | §       | §       | §       |

**Table S14.** Scenarios M6-M7 (§ design points each): design points
clearly above 7.5% and clearly below 2.5% under the null, largest
rejection rate (%) and replicates returning a usable result (%).

| Scenario | Strategy | Above 7.5% | Below 2.5% | Largest rate (%) | Usable (%) |
|----|----|----|----|----|----|
| GenBinaryY, slope § | Marginalized model with check, t | § | § | § | § |
|  | Marginalized model with check, normal | § | § | § | § |
|  | Marginalized model without check, t | § | § | § | § |
|  | GEE, Fay-Graubard, t | § | § | § | § |
|  | GEE, Mancl-DeRouen, t | § | § | § | § |
|  | GLMM, likelihood ratio | § | § | § | § |
| GenBinaryY, slope § | Marginalized model with check, t | § | § | § | § |
|  | Marginalized model with check, normal | § | § | § | § |
|  | Marginalized model without check, t | § | § | § | § |
|  | GEE, Fay-Graubard, t | § | § | § | § |
|  | GEE, Mancl-DeRouen, t | § | § | § | § |
|  | GLMM, likelihood ratio | § | § | § | § |
| GLMM-native, slope variance § | Marginalized model with check, t | § | § | § | § |
|  | Marginalized model with check, normal | § | § | § | § |
|  | Marginalized model without check, t | § | § | § | § |
|  | GEE, Fay-Graubard, t | § | § | § | § |
|  | GEE, Mancl-DeRouen, t | § | § | § | § |
|  | GLMM, likelihood ratio | § | § | § | § |
| GLMM-native, slope variance § | Marginalized model with check, t | § | § | § | § |
|  | Marginalized model with check, normal | § | § | § | § |
|  | Marginalized model without check, t | § | § | § | § |
|  | GEE, Fay-Graubard, t | § | § | § | § |
|  | GEE, Mancl-DeRouen, t | § | § | § | § |
|  | GLMM, likelihood ratio | § | § | § | § |

**Table S15.** Clinical application: outcome of the fits in § subsamples
of whole patients at each subsample size N (% of subsamples). Not
returned: no point estimate; no SE: estimate without a standard error;
separated: \|estimate\| or SE above §; usable: the remainder. For each
method and N the four columns add up to all subsamples.

| Method                 | N   | Not returned | No SE | Separated | Usable |
|------------------------|-----|--------------|-------|-----------|--------|
| GLMM                   | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
| Marginalized model     | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
| GEE, uncorrected       | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
| GEE, Fay-Graubard      | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
| GEE, Kauermann-Carroll | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
| GEE, Mancl-DeRouen     | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |
|                        | §   | §            | §     | §         | §      |

## S4. Numerical behaviour of the small-sample corrections

The corrections of Mancl and DeRouen, Kauermann and Carroll and Fay and
Graubard were computed with `geesmv`. On sparse binary data with few
clusters the three functions can loop indefinitely in compiled code, so
each call ran in a separate process killed after § seconds; the largest
share of time-outs in a cell was §%. In the earlier version of this work
the Kauermann-Carroll correction appeared to shrink rather than inflate
the variance; recomputing the three corrections from the
estimating-equation building blocks showed that this was a property of
the implementation, which refits the model internally and applies a
symmetric matrix square root to a non-symmetric matrix, and not of the
estimator. The hand computation of Fay-Graubard agreed with the package
at every grid point checked. These checks are documented with the code.

## S5. Tables of the scenario anchored to the cohort (M1c)

**Table S16.** Type I error at the nominal 5% (%): pooled rate, worst
case and mean absolute departure (MAD) over cluster counts (percentage
points), scenario M1c (§ design points).

|  | rate GB | rate NAT | rate AR1 | worst GB | worst NAT | worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | § | § | § | § | § | § | § | § | § |
| GEE, Mancl-DeRouen | § | § | § | § | § | § | § | § | § |
| Marginalized model | § | § | § | § | § | § | § | § | § |
| GLMM, Wald | § | § | § | § | § | § | § | § | § |
| Shrinkage GLM | § | § | § | § | § | § | § | § | § |
| GEE, uncorrected | § | § | § | § | § | § | § | § | § |
| GLMM, likelihood ratio | § | § | § | § | § | § | § | § | § |
| GEE, Kauermann-Carroll | § | § | § | § | § | § | § | § | § |

**Table S17.** Mean absolute departure from 5% by number of occasions
(percentage points), cell-weighted, scenario M1c (§ design points).

|                        | T = 2 | T = 3 | T \>= 4 |
|------------------------|-------|-------|---------|
| GEE, Fay-Graubard      | §     | §     | §       |
| GEE, Mancl-DeRouen     | §     | §     | §       |
| Marginalized model     | §     | §     | §       |
| GLMM, Wald             | §     | §     | §       |
| Shrinkage GLM          | §     | §     | §       |
| GEE, uncorrected       | §     | §     | §       |
| GLMM, likelihood ratio | §     | §     | §       |
| GEE, Kauermann-Carroll | §     | §     | §       |

**Table S18.** Method-cells detectably off nominal (%), exact binomial
test with Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY)
corrections within method, scenario M1c (§ design points).

|                        | BH  | BY  |
|------------------------|-----|-----|
| GEE, Fay-Graubard      | §   | §   |
| GEE, Mancl-DeRouen     | §   | §   |
| Marginalized model     | §   | §   |
| GLMM, Wald             | §   | §   |
| Shrinkage GLM          | §   | §   |
| GEE, uncorrected       | §   | §   |
| GLMM, likelihood ratio | §   | §   |
| GEE, Kauermann-Carroll | §   | §   |

**Table S19.** Type I error on each method’s own usable replicates and
on the replicates usable for all methods (common set) (%), scenario M1c
(§ design points).

|  | own GB | own NAT | own AR1 | common GB | common NAT | common AR1 | shift GB | shift NAT | shift AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | § | § | § | § | § | § | § | § | § |
| GEE, Mancl-DeRouen | § | § | § | § | § | § | § | § | § |
| Marginalized model | § | § | § | § | § | § | § | § | § |
| GLMM, Wald | § | § | § | § | § | § | § | § | § |
| Shrinkage GLM | § | § | § | § | § | § | § | § | § |
| GEE, uncorrected | § | § | § | § | § | § | § | § | § |
| GLMM, likelihood ratio | § | § | § | § | § | § | § | § | § |
| GEE, Kauermann-Carroll | § | § | § | § | § | § | § | § | § |

**Table S20.** Replicates returning a usable result by number of
occasions (%), every attempted replicate in the denominator, scenario
M1c (§ design points).

|  | GEE, Fay-Graubard | GEE, Mancl-DeRouen | Marginalized model | GLMM, Wald | Shrinkage GLM | GEE, uncorrected | GLMM, likelihood ratio | GEE, Kauermann-Carroll |
|----|----|----|----|----|----|----|----|----|
| random intercept, T = 2 | § | § | § | § | § | § | § | § |
| random intercept, T = 3 | § | § | § | § | § | § | § | § |
| random intercept, T = 4 | § | § | § | § | § | § | § | § |
| random intercept, T = 6 | § | § | § | § | § | § | § | § |

**Table S21.** Usable result and unconditional validity at 20 clusters
or fewer (%), scenario M1c (§ design points).

|  | usable | unconditional validity | design points with a usable replicate |
|----|----|----|----|
| GEE, Fay-Graubard | § | § | § |
| GEE, Mancl-DeRouen | § | § | § |
| Marginalized model | § | § | § |
| GLMM, Wald | § | § | § |
| Shrinkage GLM | § | § | § |
| GEE, uncorrected | § | § | § |
| GLMM, likelihood ratio | § | § | § |
| GEE, Kauermann-Carroll | § | § | § |

**Table S22.** Mean absolute departure from 5% by estimability stratum
(percentage points); strata by the mean usable proportion of the five
marginal methods, scenario M1c (§ design points).

|                        | §   | §   | §   | §   | §   |
|------------------------|-----|-----|-----|-----|-----|
| GEE, Fay-Graubard      | §   | §   | §   | §   | §   |
| GEE, Mancl-DeRouen     | §   | §   | §   | §   | §   |
| Marginalized model     | §   | §   | §   | §   | §   |
| GLMM, Wald             | §   | §   | §   | §   | §   |
| Shrinkage GLM          | §   | §   | §   | §   | §   |
| GEE, uncorrected       | §   | §   | §   | §   | §   |
| GLMM, likelihood ratio | §   | §   | §   | §   | §   |
| GEE, Kauermann-Carroll | §   | §   | §   | §   | §   |
| design points          | §   | §   | §   | §   | §   |

**Table S23.** Raw and size-adjusted power at the interaction contrast
(%), scenario M1c (§ design points).

|                                     | raw | size-adjusted |
|-------------------------------------|-----|---------------|
| GenBinaryY, GEE, Fay-Graubard       | §   | §             |
| GenBinaryY, GEE, Mancl-DeRouen      | §   | §             |
| GenBinaryY, Marginalized model      | §   | §             |
| GenBinaryY, GLMM, Wald              | §   | §             |
| GenBinaryY, Shrinkage GLM           | §   | §             |
| GenBinaryY, GEE, uncorrected        | §   | §             |
| GenBinaryY, GLMM, likelihood ratio  | §   | §             |
| GenBinaryY, GEE, Kauermann-Carroll  | §   | §             |
| GLMM-native, GEE, Fay-Graubard      | §   | §             |
| GLMM-native, GEE, Mancl-DeRouen     | §   | §             |
| GLMM-native, Marginalized model     | §   | §             |
| GLMM-native, GLMM, Wald             | §   | §             |
| GLMM-native, Shrinkage GLM          | §   | §             |
| GLMM-native, GEE, uncorrected       | §   | §             |
| GLMM-native, GLMM, likelihood ratio | §   | §             |
| GLMM-native, GEE, Kauermann-Carroll | §   | §             |
| GLMM-AR1, GEE, Fay-Graubard         | §   | §             |
| GLMM-AR1, GEE, Mancl-DeRouen        | §   | §             |
| GLMM-AR1, Marginalized model        | §   | §             |
| GLMM-AR1, GLMM, Wald                | §   | §             |
| GLMM-AR1, Shrinkage GLM             | §   | §             |
| GLMM-AR1, GEE, uncorrected          | §   | §             |
| GLMM-AR1, GLMM, likelihood ratio    | §   | §             |
| GLMM-AR1, GEE, Kauermann-Carroll    | §   | §             |

**Table S24.** Type I error under a normal and a t reference (%), shift
in percentage points, scenario M1c (§ design points).

|  | normal GB | normal NAT | normal AR1 | t GB | t NAT | t AR1 | shift GB | shift NAT | shift AR1 | shift at N \<= 25 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | § | § | § | § | § | § | § | § | § | § |
| GEE, Mancl-DeRouen | § | § | § | § | § | § | § | § | § | § |
| Marginalized model | § | § | § | § | § | § | § | § | § | § |
| GLMM, Wald | § | § | § | § | § | § | § | § | § | § |
| Shrinkage GLM | § | § | § | § | § | § | § | § | § | § |
| GEE, uncorrected | § | § | § | § | § | § | § | § | § | § |
| GEE, Kauermann-Carroll | § | § | § | § | § | § | § | § | § | § |

**Table S25.** Coverage of the 95% confidence interval against each
method’s own estimand (%); not available (n.a.) for the GLMM on
GenBinaryY, scenario M1c (§ design points).

|                        | GenBinaryY | GLMM-native | GLMM-AR1 |
|------------------------|------------|-------------|----------|
| GEE, Fay-Graubard      | §          | §           | §        |
| GEE, Mancl-DeRouen     | §          | §           | §        |
| Marginalized model     | §          | §           | §        |
| GLMM, Wald             | n.a.       | §           | §        |
| Shrinkage GLM          | §          | §           | §        |
| GEE, uncorrected       | §          | §           | §        |
| GEE, Kauermann-Carroll | §          | §           | §        |
