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
| GenBinaryY | M1 | ~time \* trt | -2, 0.1, 0.2, 0.2 | lv ~1, σ = (1); transition ~1, γ = (1) |
| GenBinaryY | M2 | ~time \* trt + age | -2, 0.1, 0.2, 0.01, 0.2 | lv ~1, σ = (1); transition ~1, γ = (1) |
| GenBinaryY | M3 | ~time \* trt | -2, 0.1, 0.2, 0.2 | lv ~time, σ = (1, 0.01); transition ~1, γ = (1) |
| GenBinaryY | M4 | ~time \* trt | -2, 0.2, 0.2, 0.05 | lv ~time, σ = (1, 0.01); transition ~trt, γ = (1, 0.01) |
| GenBinaryY | M5 | ~time \* trt + dose | -2, 0.1, 0.2, 0.05, 0.2 | lv ~time + dose, σ = (1, 0.01, 0.01); transition ~trt + steroids, γ = (1, 0.1, 0.05) |
| GLMM-native | M1 | ~time \* trt | -2.32, 0.116, 0.232, 0.232 | random effects ~1, D = \[1\] |
| GLMM-native | M2 | ~time \* trt + age | -2.32, 0.116, 0.232, 0.012, 0.232 | random effects ~1, D = \[1\] |
| GLMM-native | M3 | ~time \* trt | -2.32, 0.116, 0.232, 0.232 | random effects ~1 + time, D = \[1, 0.01; 0.01, 0.001\] |
| GLMM-native | M4 | ~time \* trt | -2.32, 0.232, 0.232, 0.058 | random effects ~1 + time, D = \[1, 0.01; 0.01, 0.001\] |
| GLMM-native | M5 | ~time \* trt + dose | -2.32, 0.116, 0.232, 0.058, 0.232 | random effects ~1 + time + dose, D = \[1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.001\] |
| GLMM-AR1 | M1 | ~time \* trt | -2, 0.1, 0.2, 0.2 | random effects ~1, D = \[1\]; AR(1) ρ = 0.6, variance 0.5 |
| GLMM-AR1 | M2 | ~time \* trt + age | -2, 0.1, 0.2, 0.01, 0.2 | random effects ~1, D = \[1\]; AR(1) ρ = 0.6, variance 0.5 |
| GLMM-AR1 | M3 | ~time \* trt | -2, 0.1, 0.2, 0.2 | random effects ~1 + time, D = \[1, 0.01; 0.01, 0.001\]; AR(1) ρ = 0.6, variance 0.5 |
| GLMM-AR1 | M4 | ~time \* trt | -2, 0.2, 0.2, 0.05 | random effects ~1 + time, D = \[1, 0.01; 0.01, 0.001\]; AR(1) ρ = 0.6, variance 0.5 |
| GLMM-AR1 | M5 | ~time \* trt + dose | -2, 0.1, 0.2, 0.05, 0.2 | random effects ~1 + time + dose, D = \[1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.01\]; AR(1) ρ = 0.6, variance 0.5 |

## S2. Software

**Table S2.** Versions in the project library used for every rerun and
for the tables of this paper. The original run was executed under an
earlier patch release of R; its datasets were reproduced bit for bit
under the current release. The full `sessionInfo()` is provided with the
code (`results_v4/S7/sessionInfo.txt`).

| Software   | Version |
|------------|---------|
| R          | 4.6.1   |
| MMLB       | 1.0     |
| lme4       | 2.0.6   |
| geepack    | 1.3.13  |
| geesmv     | 1.3     |
| data.table | 1.18.4  |
| crew       | 1.3.2   |
| mirai      | 2.7.2   |
| Matrix     | 1.7.6   |

## S3. Tables of the main grid

**Table S3.** Type I error at the nominal 5% (%): pooled rate, worst
case and mean absolute departure (MAD) over cluster counts (percentage
points), main grid, with the marginalized model analysed as in the
original study, without the convergence check (sensitivity analysis).
Only the row of the marginalized model differs from Table 2. GB =
GenBinaryY, NAT = GLMM-native, AR1 = GLMM-AR1.

|  | rate GB | rate NAT | rate AR1 | worst GB | worst NAT | worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 4.7 | 2.8 | 2.9 | 1.0 | 0.5 | 0.6 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 8.6 | 5.9 | 6.0 | 1.5 | 0.9 | 1.0 |
| Marginalized model | 6.6 | 6.2 | 6.4 | 2.2 | 1.5 | 1.6 | 1.6 | 1.2 | 1.4 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 4.5 | 1.8 | 2.2 | 3.5 | 1.4 | 1.6 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 4.9 | 2.3 | 2.4 | 4.2 | 1.6 | 1.9 |
| GEE, uncorrected | 8.3 | 7.6 | 7.7 | 21.0 | 17.4 | 18.3 | 3.8 | 3.0 | 3.1 |
| GLMM, likelihood ratio | 7.8 | 5.8 | 6.1 | 3.6 | 1.8 | 1.8 | 2.8 | 0.8 | 1.1 |
| GEE, Kauermann-Carroll | 13.1 | 10.4 | 10.7 | 9.2 | 6.5 | 6.7 | 7.9 | 5.3 | 5.6 |

**Table S4.** Mean absolute departure from 5% by number of occasions
(percentage points), cell-weighted, main grid, marginalized model with
the convergence check. Monte Carlo standard errors by block bootstrap;
values given to the digit they support.

|                        | T = 2 | T = 3 | T \>= 4 |
|------------------------|-------|-------|---------|
| GEE, Fay-Graubard      | 2.4   | 1.0   | 1.2     |
| GEE, Mancl-DeRouen     | 2.5   | 1.2   | 1.6     |
| Marginalized model     | 1.5   | 1.0   | 0.7     |
| GLMM, Wald             | 5.6   | 4.0   | 2.1     |
| Shrinkage GLM          | 10.8  | 4.8   | 2.2     |
| GEE, uncorrected       | 2.2   | 3.6   | 3.8     |
| GLMM, likelihood ratio | 10.4  | 2.4   | 1.4     |
| GEE, Kauermann-Carroll | 2.5   | 3.5   | 7.8     |

**Table S5.** Method-cells detectably off nominal (%), exact binomial
test with Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY)
corrections within method, main grid, marginalized model with the
convergence check. Monte Carlo standard errors by block bootstrap;
values given to the digit they support.

|                        | BH  | BY  |
|------------------------|-----|-----|
| GEE, Fay-Graubard      | 23  | 14  |
| GEE, Mancl-DeRouen     | 28  | 20  |
| Marginalized model     | 3   | 1   |
| GLMM, Wald             | 51  | 42  |
| Shrinkage GLM          | 53  | 45  |
| GEE, uncorrected       | 46  | 35  |
| GLMM, likelihood ratio | 35  | 28  |
| GEE, Kauermann-Carroll | 87  | 80  |

**Table S6.** Type I error on each method’s own usable replicates and on
the replicates usable for all methods (common set) (%), main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.

|  | own GB | own NAT | own AR1 | common GB | common NAT | common AR1 | shift GB | shift NAT | shift AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 5.8 | 5.4 | 5.4 | -0.1 | 0.1 | 0.1 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 6.1 | 5.7 | 5.7 | -0.0 | 0.2 | 0.1 |
| Marginalized model | 5.5 | 5.2 | 5.4 | 5.3 | 5.1 | 5.3 | -0.3 | -0.1 | -0.1 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 8.6 | 6.3 | 6.6 | -0.0 | -0.0 | 0.0 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 8.6 | 6.3 | 6.6 | -0.7 | -0.4 | -0.4 |
| GEE, uncorrected | 8.3 | 7.6 | 7.7 | 7.8 | 7.4 | 7.5 | -0.5 | -0.2 | -0.2 |
| GLMM, likelihood ratio | 7.8 | 5.8 | 6.1 | 6.8 | 5.0 | 5.3 | -1.0 | -0.8 | -0.8 |
| GEE, Kauermann-Carroll | 13.1 | 10.4 | 10.7 | 13.5 | 10.9 | 11.2 | 0.4 | 0.5 | 0.5 |

**Table S7.** Replicates returning a usable result by number of
occasions (%), every attempted replicate in the denominator, main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.

|  | GEE, Fay-Graubard | GEE, Mancl-DeRouen | Marginalized model | GLMM, Wald | Shrinkage GLM | GEE, uncorrected | GLMM, likelihood ratio | GEE, Kauermann-Carroll |
|----|----|----|----|----|----|----|----|----|
| random intercept, T = 2 | 72.7 | 72.7 | 68.2 | 68.6 | 74.0 | 78.7 | 98.9 | 72.6 |
| random intercept, T = 3 | 88.3 | 88.3 | 87.1 | 87.8 | 88.9 | 91.7 | 99.6 | 88.3 |
| random intercept, T = 4 | 93.8 | 93.7 | 92.9 | 93.7 | 94.1 | 95.9 | 99.9 | 93.7 |
| random intercept, T = 6 | 97.7 | 97.6 | 97.4 | 97.8 | 97.9 | 98.6 | 100.0 | 97.7 |
| random intercept, T = 8 | 98.9 | 98.8 | 98.9 | 99.0 | 99.0 | 99.4 | 100.0 | 98.8 |
| random intercept, T = 10 | 99.5 | 99.4 | 99.5 | 99.6 | 99.6 | 99.7 | 100.0 | 99.4 |
| random slopes, T = 2 | 73.1 | 73.1 | 71.0 | 0.0 | 0.0 | 78.8 | 0.1 | 73.1 |
| random slopes, T = 3 | 88.4 | 88.4 | 87.9 | 56.7 | 59.3 | 91.8 | 66.5 | 88.4 |
| random slopes, T = 4 | 94.1 | 94.0 | 93.8 | 90.8 | 93.7 | 96.1 | 99.4 | 94.0 |
| random slopes, T = 6 | 98.0 | 97.9 | 97.8 | 97.7 | 98.1 | 98.8 | 100.0 | 97.9 |
| random slopes, T = 8 | 99.1 | 99.0 | 99.0 | 99.2 | 99.2 | 99.5 | 100.0 | 99.1 |
| random slopes, T = 10 | 99.5 | 99.4 | 99.4 | 99.6 | 99.6 | 99.7 | 100.0 | 99.5 |

**Table S8.** Usable result and unconditional validity at 20 clusters or
fewer (%), main grid, marginalized model with the convergence check.
Monte Carlo standard errors by block bootstrap; values given to the
digit they support.

|  | usable | unconditional validity | design points with a usable replicate |
|----|----|----|----|
| GEE, Fay-Graubard | 68.9 | 63.8 | 270 |
| GEE, Mancl-DeRouen | 68.7 | 62.3 | 270 |
| Marginalized model | 67.9 | 64.4 | 270 |
| GLMM, Wald | 64.5 | 61.0 | 244 |
| Shrinkage GLM | 66.9 | 62.6 | 251 |
| GEE, uncorrected | 78.7 | 65.2 | 270 |
| GLMM, likelihood ratio | 85.8 | 79.9 | 255 |
| GEE, Kauermann-Carroll | 68.7 | 62.2 | 270 |

**Table S9.** Mean absolute departure from 5% by estimability stratum
(percentage points); strata by the mean usable proportion of the five
marginal methods, main grid, marginalized model with the convergence
check. Monte Carlo standard errors by block bootstrap; values given to
the digit they support.

|                        | \<50% | 50-80% | 80-95% | 95-99.9% | 100% |
|------------------------|-------|--------|--------|----------|------|
| GEE, Fay-Graubard      | 3.9   | 2.8    | 2.8    | 1.4      | 0.7  |
| GEE, Mancl-DeRouen     | 4.3   | 4.2    | 4.2    | 1.7      | 0.6  |
| Marginalized model     | 2.6   | 1.2    | 1.0    | 0.9      | 0.7  |
| GLMM, Wald             | 5     | 3      | 2.6    | 3.1      | 2.3  |
| Shrinkage GLM          | 10    | 7      | 4.0    | 3.5      | 2.3  |
| GEE, uncorrected       | 8.4   | 9.9    | 8.8    | 3.6      | 1.1  |
| GLMM, likelihood ratio | 4.6   | 5.2    | 5.0    | 3.2      | 1.2  |
| GEE, Kauermann-Carroll | 3.0   | 1.9    | 3.0    | 4.7      | 8.4  |
| design points          | 83    | 100    | 140    | 346      | 771  |

**Table S10.** Raw and size-adjusted power at the interaction contrast
(%), main grid, marginalized model with the convergence check. Monte
Carlo standard errors by block bootstrap; values given to the digit they
support.

|                                     | raw  | size-adjusted |
|-------------------------------------|------|---------------|
| GenBinaryY, GEE, Fay-Graubard       | 23.3 | 22.0          |
| GenBinaryY, GEE, Mancl-DeRouen      | 23.5 | 21.9          |
| GenBinaryY, Marginalized model      | 24.1 | 23.4          |
| GenBinaryY, GLMM, Wald              | 31.9 | 24.2          |
| GenBinaryY, Shrinkage GLM           | 32.2 | 23.6          |
| GenBinaryY, GEE, uncorrected        | 26.2 | 21.5          |
| GenBinaryY, GLMM, likelihood ratio  | 30.1 | 26.4          |
| GenBinaryY, GEE, Kauermann-Carroll  | 34.6 | 23.0          |
| GLMM-native, GEE, Fay-Graubard      | 23.7 | 23.1          |
| GLMM-native, GEE, Mancl-DeRouen     | 23.9 | 23.1          |
| GLMM-native, Marginalized model     | 24.4 | 24.2          |
| GLMM-native, GLMM, Wald             | 28.0 | 23.5          |
| GLMM-native, Shrinkage GLM          | 28.1 | 23.2          |
| GLMM-native, GEE, uncorrected       | 26.5 | 22.8          |
| GLMM-native, GLMM, likelihood ratio | 26.8 | 26.3          |
| GLMM-native, GEE, Kauermann-Carroll | 31.9 | 24.1          |
| GLMM-AR1, GEE, Fay-Graubard         | 24.1 | 23.5          |
| GLMM-AR1, GEE, Mancl-DeRouen        | 24.4 | 23.4          |
| GLMM-AR1, Marginalized model        | 25.2 | 24.6          |
| GLMM-AR1, GLMM, Wald                | 29.2 | 25.4          |
| GLMM-AR1, Shrinkage GLM             | 29.3 | 25.0          |
| GLMM-AR1, GEE, uncorrected          | 27.0 | 23.0          |
| GLMM-AR1, GLMM, likelihood ratio    | 27.9 | 26.7          |
| GLMM-AR1, GEE, Kauermann-Carroll    | 32.7 | 24.5          |

**Table S11.** Type I error under a normal and a t reference (%), shift
in percentage points, main grid, marginalized model with the convergence
check. Monte Carlo standard errors by block bootstrap; values given to
the digit they support.

|  | normal GB | normal NAT | normal AR1 | t GB | t NAT | t AR1 | shift GB | shift NAT | shift AR1 | shift at N \<= 25 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 5.1 | 4.6 | 4.6 | -0.8 | -0.7 | -0.7 | 2.2 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 5.3 | 4.7 | 4.8 | -0.8 | -0.8 | -0.8 | 2.4 |
| Marginalized model | 5.5 | 5.2 | 5.4 | 4.8 | 4.5 | 4.6 | -0.7 | -0.7 | -0.8 | 2.1 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 8.0 | 5.8 | 6.0 | -0.7 | -0.5 | -0.6 | 1.4 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 8.6 | 6.1 | 6.3 | -0.7 | -0.5 | -0.6 | 1.5 |
| GEE, uncorrected | 8.3 | 7.6 | 7.7 | 7.2 | 6.5 | 6.6 | -1.1 | -1.1 | -1.1 | 3.7 |
| GEE, Kauermann-Carroll | 13.1 | 10.4 | 10.7 | 12.0 | 9.5 | 9.7 | -1.1 | -1.0 | -1.0 | 2.7 |

**Table S12.** Coverage of the 95% confidence interval against each
method’s own estimand (%); not available (n.a.) for the GLMM on
GenBinaryY, main grid, marginalized model with the convergence check.
Monte Carlo standard errors by block bootstrap; values given to the
digit they support.

|                        | GenBinaryY | GLMM-native | GLMM-AR1 |
|------------------------|------------|-------------|----------|
| GEE, Fay-Graubard      | 94.6       | 94.9        | 94.9     |
| GEE, Mancl-DeRouen     | 94.2       | 94.7        | 94.6     |
| Marginalized model     | 94.4       | 94.7        | 94.5     |
| GLMM, Wald             | n.a.       | 93.5        | 93.2     |
| Shrinkage GLM          | 90.1       | 92.9        | 92.6     |
| GEE, uncorrected       | 91.9       | 92.5        | 92.5     |
| GEE, Kauermann-Carroll | 86.4       | 89.1        | 88.8     |

**Table S13.** Pooled type I error (%) on the main grid with fits that
returned warnings kept (U0, primary) and excluded (U1),
replicate-weighted over the null cells. Binomial Monte Carlo standard
error at most 0.051 points.

| Method                         | GB, U0 | GB, U1 | NAT, U0 | NAT, U1 | AR1, U0 | AR1, U1 |
|--------------------------------|--------|--------|---------|---------|---------|---------|
| GEE, Fay-Graubard              | 5.8    | 5.8    | 5.3     | 5.3     | 5.3     | 5.3     |
| GEE, Mancl-DeRouen             | 6.1    | 6.1    | 5.5     | 5.5     | 5.6     | 5.6     |
| Marginalized model, with check | 5.5    | 5.4    | 5.2     | 5.2     | 5.4     | 5.3     |
| GLMM, Wald                     | 8.6    | 5.8    | 6.3     | 4.1     | 6.6     | 4.4     |
| Shrinkage GLM                  | 9.3    | 6.3    | 6.6     | 4.4     | 6.9     | 4.6     |
| GEE, uncorrected               | 8.3    | 8.2    | 7.6     | 7.5     | 7.7     | 7.6     |
| GLMM, likelihood ratio         | 7.8    | 7.2    | 5.8     | 5.3     | 6.1     | 5.6     |
| GEE, Kauermann-Carroll         | 13.1   | 13.1   | 10.4    | 10.4    | 10.7    | 10.7    |

**Table S14.** Scenarios M6-M7 (96 design points each): design points
clearly above 7.5% and clearly below 2.5% under the null, largest
rejection rate (%) and replicates returning a usable result (%).

| Scenario | Strategy | Above 7.5% | Below 2.5% | Largest rate (%) | Usable (%) |
|----|----|----|----|----|----|
| GenBinaryY, slope 0.1 | Marginalized model with check, t | 0 | 3 | 7.0 | 90.5 |
|  | Marginalized model with check, normal | 0 | 2 | 7.1 | 90.5 |
|  | Marginalized model without check, t | 0 | 3 | 7.5 | 88.8 |
|  | GEE, Fay-Graubard, t | 0 | 6 | 9.1 | 91.2 |
|  | GEE, Mancl-DeRouen, t | 8 | 7 | 13.3 | 91.2 |
|  | GLMM, likelihood ratio | 0 | 3 | 7.3 | 83.3 |
| GenBinaryY, slope 0.25 | Marginalized model with check, t | 0 | 2 | 7.8 | 90.1 |
|  | Marginalized model with check, normal | 2 | 1 | 9.7 | 90.1 |
|  | Marginalized model without check, t | 1 | 2 | 9.4 | 88.3 |
|  | GEE, Fay-Graubard, t | 2 | 6 | 10.4 | 90.8 |
|  | GEE, Mancl-DeRouen, t | 10 | 6 | 14.1 | 90.8 |
|  | GLMM, likelihood ratio | 0 | 4 | 7.5 | 83.2 |
| GLMM-native, slope variance 0.01 | Marginalized model with check, t | 0 | 1 | 6.9 | 91.5 |
|  | Marginalized model with check, normal | 0 | 0 | 7.2 | 91.5 |
|  | Marginalized model without check, t | 1 | 2 | 10.0 | 90.2 |
|  | GEE, Fay-Graubard, t | 0 | 6 | 7.0 | 92.3 |
|  | GEE, Mancl-DeRouen, t | 0 | 6 | 9.3 | 92.2 |
|  | GLMM, likelihood ratio | 0 | 0 | 6.1 | 83.3 |
| GLMM-native, slope variance 0.04 | Marginalized model with check, t | 0 | 3 | 7.8 | 91.7 |
|  | Marginalized model with check, normal | 1 | 0 | 10.6 | 91.7 |
|  | Marginalized model without check, t | 6 | 2 | 11.1 | 90.4 |
|  | GEE, Fay-Graubard, t | 0 | 4 | 6.7 | 92.3 |
|  | GEE, Mancl-DeRouen, t | 0 | 6 | 9.0 | 92.3 |
|  | GLMM, likelihood ratio | 0 | 1 | 6.6 | 83.3 |

**Table S15.** Clinical application: outcome of the fits in 500
subsamples of whole patients at each subsample size N (% of subsamples).
Not returned: no point estimate; no SE: estimate without a standard
error; separated: \|estimate\| or SE above 10; usable: the remainder.
For each method and N the four columns add up to all subsamples.

| Method                 | N   | Not returned | No SE | Separated | Usable |
|------------------------|-----|--------------|-------|-----------|--------|
| GLMM                   | 10  | 21.8         | 0.0   | 0.0       | 78.2   |
|                        | 15  | 14.6         | 0.0   | 0.0       | 85.4   |
|                        | 20  | 5.4          | 0.0   | 0.0       | 94.6   |
|                        | 25  | 3.4          | 0.0   | 0.0       | 96.6   |
|                        | 50  | 0.0          | 0.0   | 0.0       | 100.0  |
| Marginalized model     | 10  | 23.0         | 0.0   | 0.0       | 77.0   |
|                        | 15  | 14.8         | 0.0   | 0.0       | 85.2   |
|                        | 20  | 5.6          | 0.0   | 0.0       | 94.4   |
|                        | 25  | 3.4          | 0.0   | 0.0       | 96.6   |
|                        | 50  | 0.0          | 0.0   | 0.0       | 100.0  |
| GEE, uncorrected       | 10  | 24.0         | 0.0   | 6.0       | 70.0   |
|                        | 15  | 14.8         | 0.0   | 4.4       | 80.8   |
|                        | 20  | 5.6          | 0.0   | 3.4       | 91.0   |
|                        | 25  | 3.4          | 0.0   | 1.2       | 95.4   |
|                        | 50  | 0.0          | 0.0   | 0.0       | 100.0  |
| GEE, Fay-Graubard      | 10  | 24.0         | 41.4  | 0.4       | 34.2   |
|                        | 15  | 14.8         | 23.8  | 0.2       | 61.2   |
|                        | 20  | 5.6          | 14.2  | 0.2       | 80.0   |
|                        | 25  | 3.4          | 8.0   | 0.2       | 88.4   |
|                        | 50  | 0.0          | 0.2   | 0.0       | 99.8   |
| GEE, Kauermann-Carroll | 10  | 24.0         | 41.6  | 0.4       | 34.0   |
|                        | 15  | 14.8         | 24.0  | 0.2       | 61.0   |
|                        | 20  | 5.6          | 14.2  | 0.2       | 80.0   |
|                        | 25  | 3.4          | 8.0   | 0.2       | 88.4   |
|                        | 50  | 0.0          | 0.2   | 0.0       | 99.8   |
| GEE, Mancl-DeRouen     | 10  | 24.0         | 60.4  | 0.6       | 15.0   |
|                        | 15  | 14.8         | 47.4  | 0.6       | 37.2   |
|                        | 20  | 5.6          | 31.8  | 0.4       | 62.2   |
|                        | 25  | 3.4          | 16.6  | 0.0       | 80.0   |
|                        | 50  | 0.0          | 1.2   | 0.0       | 98.8   |

## S4. Numerical behaviour of the small-sample corrections

The corrections of Mancl and DeRouen, Kauermann and Carroll and Fay and
Graubard were computed with `geesmv`. On sparse binary data with few
clusters the three functions can loop indefinitely in compiled code, so
each call ran in a separate process killed after 5 seconds; the largest
share of time-outs in a cell was 17%. In the earlier version of this
work the Kauermann-Carroll correction appeared to shrink rather than
inflate the variance; recomputing the three corrections from the
estimating-equation building blocks showed that this was a property of
the implementation, which refits the model internally and applies a
symmetric matrix square root to a non-symmetric matrix, and not of the
estimator. The hand computation of Fay-Graubard agreed with the package
at every grid point checked. These checks are documented with the code.

## S5. Tables of the scenario anchored to the cohort (M1c)

**Table S16.** Type I error at the nominal 5% (%): pooled rate, worst
case and mean absolute departure (MAD) over cluster counts (percentage
points), scenario M1c (72 design points).

|  | rate GB | rate NAT | rate AR1 | worst GB | worst NAT | worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.1 | 5.0 | 5.0 | 1.8 | 1.2 | 1.8 | 0.6 | 0.5 | 0.5 |
| GEE, Mancl-DeRouen | 4.5 | 4.5 | 4.4 | 1.2 | 1.1 | 1.1 | 0.8 | 0.7 | 0.8 |
| Marginalized model | 7.0 | 6.8 | 6.9 | 2.4 | 2.4 | 2.9 | 2.0 | 1.8 | 2.0 |
| GLMM, Wald | 7.7 | 6.0 | 6.2 | 4.5 | 2.8 | 2.5 | 2.5 | 1.3 | 1.5 |
| Shrinkage GLM | 11.2 | 7.9 | 8.3 | 6.6 | 3.3 | 4.5 | 6.2 | 2.9 | 3.3 |
| GEE, uncorrected | 8.3 | 8.1 | 8.2 | 10 | 9 | 10 | 3.6 | 3.3 | 3.5 |
| GLMM, likelihood ratio | 10.5 | 7.5 | 8.0 | 7.0 | 4.8 | 5.6 | 5.5 | 2.5 | 3.0 |
| GEE, Kauermann-Carroll | 6.5 | 8.1 | 8.1 | 2.3 | 4.3 | 4.1 | 1.5 | 3.0 | 3.0 |

**Table S17.** Mean absolute departure from 5% by number of occasions
(percentage points), cell-weighted, scenario M1c (72 design points).

|                        | T = 2 | T = 3 | T \>= 4 |
|------------------------|-------|-------|---------|
| GEE, Fay-Graubard      | 2.5   | 0.7   | 1.2     |
| GEE, Mancl-DeRouen     | 3.0   | 1.0   | 0.9     |
| Marginalized model     | 1.9   | 1.8   | 1.9     |
| GLMM, Wald             | 1.5   | 1.3   | 3.6     |
| Shrinkage GLM          | 7.1   | 2.6   | 3.9     |
| GEE, uncorrected       | 0.9   | 3.1   | 4.7     |
| GLMM, likelihood ratio | 6.7   | 3.2   | 2.5     |
| GEE, Kauermann-Carroll | 2.6   | 3.0   | 3.0     |

**Table S18.** Method-cells detectably off nominal (%), exact binomial
test with Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY)
corrections within method, scenario M1c (72 design points).

|                        | BH  | BY  |
|------------------------|-----|-----|
| GEE, Fay-Graubard      | 28  | 25  |
| GEE, Mancl-DeRouen     | 26  | 19  |
| Marginalized model     | 57  | 42  |
| GLMM, Wald             | 54  | 38  |
| Shrinkage GLM          | 78  | 72  |
| GEE, uncorrected       | 61  | 53  |
| GLMM, likelihood ratio | 67  | 64  |
| GEE, Kauermann-Carroll | 82  | 68  |

**Table S19.** Type I error on each method’s own usable replicates and
on the replicates usable for all methods (common set) (%), scenario M1c
(72 design points).

|  | own GB | own NAT | own AR1 | common GB | common NAT | common AR1 | shift GB | shift NAT | shift AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.1 | 5.0 | 5.0 | 4.7 | 4.7 | 4.7 | -0.3 | -0.3 | -0.3 |
| GEE, Mancl-DeRouen | 4.5 | 4.5 | 4.4 | 4.2 | 4.2 | 4.2 | -0.2 | -0.2 | -0.2 |
| Marginalized model | 7.0 | 6.8 | 6.9 | 6.4 | 6.2 | 6.4 | -0.6 | -0.6 | -0.5 |
| GLMM, Wald | 7.7 | 6.0 | 6.2 | 7.7 | 6.0 | 6.3 | 0.0 | 0.0 | 0.1 |
| Shrinkage GLM | 11.2 | 7.9 | 8.3 | 7.7 | 6.0 | 6.3 | -3.5 | -1.9 | -2.0 |
| GEE, uncorrected | 8.3 | 8.1 | 8.2 | 7.5 | 7.4 | 7.6 | -0.8 | -0.7 | -0.7 |
| GLMM, likelihood ratio | 10.5 | 7.5 | 8.0 | 7.3 | 4.9 | 5.3 | -3.3 | -2.6 | -2.6 |
| GEE, Kauermann-Carroll | 6.5 | 8.1 | 8.1 | 6.5 | 7.9 | 8.0 | -0.1 | -0.2 | -0.1 |

**Table S20.** Replicates returning a usable result by number of
occasions (%), every attempted replicate in the denominator, scenario
M1c (72 design points).

|  | GEE, Fay-Graubard | GEE, Mancl-DeRouen | Marginalized model | GLMM, Wald | Shrinkage GLM | GEE, uncorrected | GLMM, likelihood ratio | GEE, Kauermann-Carroll |
|----|----|----|----|----|----|----|----|----|
| random intercept, T = 2 | 81.3 | 81.3 | 69.6 | 68.8 | 82.9 | 84.1 | 99.9 | 81.3 |
| random intercept, T = 3 | 93.9 | 93.8 | 89.1 | 90.8 | 94.6 | 95.1 | 100.0 | 93.9 |
| random intercept, T = 4 | 96.6 | 96.4 | 94.3 | 96.0 | 97.0 | 97.5 | 100.0 | 96.5 |
| random intercept, T = 6 | 98.4 | 98.1 | 97.5 | 98.6 | 98.7 | 98.9 | 100.0 | 98.2 |

**Table S21.** Usable result and unconditional validity at 20 clusters
or fewer (%), scenario M1c (72 design points).

|  | usable | unconditional validity | design points with a usable replicate |
|----|----|----|----|
| GEE, Fay-Graubard | 85.9 | 81.5 | 36 |
| GEE, Mancl-DeRouen | 85.5 | 81.7 | 36 |
| Marginalized model | 78.8 | 73.3 | 36 |
| GLMM, Wald | 79.8 | 75.6 | 36 |
| Shrinkage GLM | 87.2 | 79.3 | 36 |
| GEE, uncorrected | 88.6 | 79.7 | 36 |
| GLMM, likelihood ratio | 99.9 | 90.0 | 36 |
| GEE, Kauermann-Carroll | 85.7 | 79.8 | 36 |

**Table S22.** Mean absolute departure from 5% by estimability stratum
(percentage points); strata by the mean usable proportion of the five
marginal methods, scenario M1c (72 design points).

|                        | \<50% | 50-80% | 80-95% | 95-99.9% | 100% |
|------------------------|-------|--------|--------|----------|------|
| GEE, Fay-Graubard      | 4.7   | 2.3    | 2.4    | 0.8      | 0.6  |
| GEE, Mancl-DeRouen     | 5.0   | 2.3    | 2.4    | 0.8      | 0.6  |
| Marginalized model     | 1     | 1.6    | 2.5    | 2.1      | 1.1  |
| GLMM, Wald             | 1     | 1.9    | 1.5    | 2.2      | 5.0  |
| Shrinkage GLM          | 8     | 9      | 4.2    | 2.9      | 5.0  |
| GEE, uncorrected       | 1     | 6      | 5.4    | 2.8      | 1.8  |
| GLMM, likelihood ratio | 7     | 8.1    | 5.5    | 2.7      | 1.1  |
| GEE, Kauermann-Carroll | 4.3   | 2.1    | 2.0    | 3.1      | 3.5  |
| design points          | 3     | 7      | 15     | 33       | 14   |

**Table S23.** Raw and size-adjusted power at the interaction contrast
(%), scenario M1c (72 design points).

|                                     | raw  | size-adjusted |
|-------------------------------------|------|---------------|
| GenBinaryY, GEE, Fay-Graubard       | 9.7  | 9.4           |
| GenBinaryY, GEE, Mancl-DeRouen      | 8.6  | 9.2           |
| GenBinaryY, Marginalized model      | 11.2 | 8.7           |
| GenBinaryY, GLMM, Wald              | 10.6 | 4.7           |
| GenBinaryY, Shrinkage GLM           | 14.2 | 4.1           |
| GenBinaryY, GEE, uncorrected        | 13.4 | 9.0           |
| GenBinaryY, GLMM, likelihood ratio  | 15.7 | 8.9           |
| GenBinaryY, GEE, Kauermann-Carroll  | 10.0 | 8.0           |
| GLMM-native, GEE, Fay-Graubard      | 9.5  | 9.1           |
| GLMM-native, GEE, Mancl-DeRouen     | 8.6  | 9.0           |
| GLMM-native, Marginalized model     | 11.1 | 8.8           |
| GLMM-native, GLMM, Wald             | 8.6  | 4.8           |
| GLMM-native, Shrinkage GLM          | 10.4 | 4.8           |
| GLMM-native, GEE, uncorrected       | 13.1 | 8.9           |
| GLMM-native, GLMM, likelihood ratio | 11.8 | 9.1           |
| GLMM-native, GEE, Kauermann-Carroll | 11.9 | 8.1           |
| GLMM-AR1, GEE, Fay-Graubard         | 9.8  | 9.7           |
| GLMM-AR1, GEE, Mancl-DeRouen        | 8.7  | 9.5           |
| GLMM-AR1, Marginalized model        | 11.2 | 8.6           |
| GLMM-AR1, GLMM, Wald                | 8.9  | 4.8           |
| GLMM-AR1, Shrinkage GLM             | 11.0 | 4.8           |
| GLMM-AR1, GEE, uncorrected          | 13.5 | 9.3           |
| GLMM-AR1, GLMM, likelihood ratio    | 12.5 | 9.1           |
| GLMM-AR1, GEE, Kauermann-Carroll    | 12.2 | 8.0           |

**Table S24.** Type I error under a normal and a t reference (%), shift
in percentage points, scenario M1c (72 design points).

|  | normal GB | normal NAT | normal AR1 | t GB | t NAT | t AR1 | shift GB | shift NAT | shift AR1 | shift at N \<= 25 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.1 | 5.0 | 5.0 | 3.6 | 3.6 | 3.5 | -1.5 | -1.4 | -1.5 | 1.9 |
| GEE, Mancl-DeRouen | 4.5 | 4.5 | 4.4 | 3.2 | 3.2 | 3.2 | -1.3 | -1.3 | -1.2 | 1.6 |
| Marginalized model | 7.0 | 6.8 | 6.9 | 5.1 | 5.0 | 5.1 | -1.9 | -1.8 | -1.9 | 2.4 |
| GLMM, Wald | 7.7 | 6.0 | 6.2 | 6.6 | 5.1 | 5.2 | -1.1 | -0.9 | -1.0 | 1.1 |
| Shrinkage GLM | 11.2 | 7.9 | 8.3 | 9.7 | 6.7 | 7.0 | -1.5 | -1.2 | -1.3 | 1.6 |
| GEE, uncorrected | 8.3 | 8.1 | 8.2 | 6.0 | 5.9 | 5.8 | -2.3 | -2.2 | -2.4 | 3.1 |
| GEE, Kauermann-Carroll | 6.5 | 8.1 | 8.1 | 4.9 | 6.1 | 6.2 | -1.7 | -1.9 | -1.9 | 2.2 |

**Table S25.** Coverage of the 95% confidence interval against each
method’s own estimand (%); not available (n.a.) for the GLMM on
GenBinaryY, scenario M1c (72 design points).

|                        | GenBinaryY | GLMM-native | GLMM-AR1 |
|------------------------|------------|-------------|----------|
| GEE, Fay-Graubard      | 93.9       | 94.5        | 94.5     |
| GEE, Mancl-DeRouen     | 94.7       | 95.2        | 95.1     |
| Marginalized model     | 93.0       | 93.2        | 93.1     |
| GLMM, Wald             | n.a.       | 95.1        | 94.9     |
| Shrinkage GLM          | 89.4       | 93.0        | 92.4     |
| GEE, uncorrected       | 90.2       | 91.2        | 91.2     |
| GEE, Kauermann-Carroll | 93.7       | 92.6        | 92.3     |
