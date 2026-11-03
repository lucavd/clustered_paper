# Supplementary material for: Choosing an analysis for clustered binary
outcomes when clusters are few and series are short: a simulation study
and two clinical applications


This supplement accompanies the manuscript. Every number is read from
the same registry of values as the main text. GB = GenBinaryY, NAT =
GLMM-native, AR1 = GLMM-AR1, COP = Copula-AR1; MMM = marginalized
multilevel model. In every table the Kauermann-Carroll correction is
computed with *fewclust* (Section S4). Numbers in square brackets refer
to the reference list of the main text.

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
| GLMM-AR1 | M5 | ~time \* trt + dose | -2, 0.1, 0.2, 0.05, 0.2 | random effects ~1 + time + dose, D = \[1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.001\]; AR(1) ρ = 0.6, variance 0.5 |

Scenarios and mechanisms not in the table use the same code. M1c has the
mean model of M1 with the prevalence (0.59) and between-subject variance
(3.5) of the clinical cohort. M6 and M7 let the latent standard
deviation of GenBinaryY grow by 0.1 or 0.25 per occasion, and give
GLMM-native a random-slope variance of 0.01 or 0.04. Copula-AR1 has the
marginal coefficients of each scenario and latent correlation ρ = 0.6
(Section S8). Covariates are drawn as follows: age, per subject, rnorm(n
= N, mean = 50, sd = 10), then standardized; dose, per observation,
rnorm(n = length(id)), then standardized; steroids, per subject,
sample(x = 1:3, size = N, replace = TRUE), in the transition model of
GenBinaryY M5 only. Treatment is rbinom(n = N, size = 1, prob = 0.5) per
subject.

The conditional coefficients of GLMM-AR1 are calibrated on each
replicate, so that the implied marginal coefficients equal the target.
In M6 and M7, with 10 occasions the latent standard deviation of
GenBinaryY grows from 1 to 1.9 or 3.25, against 1.09 in M3. The
autoregressive parameters of GLMM-AR1 add a serial dependence to the
cluster effect.<!-- noclaim: justification of parameter values --> In 44
design points of the main grid and 6 of M1c, all in GLMM-AR1, only 971
to 999 datasets could be generated, because the per-replicate
calibration encountered a singular system; every replicate that could
not be generated had a degenerate design, defined in Section 2.2 of the
main text.

## S2. Software

**Table S2.** Versions in the project library and, for *fewclust*, in a
separate library built from the commit shown, used for every rerun and
for the tables of this paper. The simulations were first run under an
earlier patch release of R; their datasets were reproduced bit for bit
under the current release. The full *sessionInfo()* is provided with the
code (*reports/morris/S8b_logs/sessionInfo.txt*, with *fewclust* and
*clubSandwich* loaded). The marginalized model used 20 Gauss-Hermite
points; Bie et al. used more with large random-effect variances
(reference \[19\] of the main text), and more points were not examined.

| Software         | Version |
|------------------|---------|
| R                | 4.6.1   |
| MMLB             | 1.0     |
| lme4             | 2.0.6   |
| geepack          | 1.3.13  |
| geesmv           | 1.3     |
| fewclust         | 0.2.0   |
| fewclust, commit | 8538cad |
| clubSandwich     | 0.7.0   |
| data.table       | 1.18.4  |
| crew             | 1.3.2   |
| mirai            | 2.7.2   |
| Matrix           | 1.7.6   |

Replicates were run in parallel, each from its own seed, so that any
dataset can be regenerated alone; with the L’Ecuyer-CMRG generator the
sequences started from different seeds do not overlap in practice, and
separate streams \[21\] were not used. The 2,376 design points comprise
4,752,000 datasets, analysed with about 52.1 million model fits over all
runs (a likelihood-ratio test counted once; the logistic fits of the
sign-flip test and of its interval search not counted), and every
analysis added after the main run regenerated them from their seeds.

## S3. Tables of the main grid

The marginalized model without the check is reported as a sensitivity
analysis (Table S3). Fits that returned with warnings were kept;
excluding them is a sensitivity analysis (Table S13). Counts of
miscalibrated cells used exact binomial tests with the
Benjamini-Hochberg and Benjamini-Yekutieli corrections.

With the Benjamini-Hochberg correction within method, the marginalized
model was detectably off nominal in 1% of its cells, Fay-Graubard in
23%, Mancl-DeRouen in 28% and Kauermann-Carroll in 38% (Table
S5).<!-- claim: cl:C10b|BH as=larger --> The uncorrected GEE was above
7.5% by 20 points at N = 10 on GenBinaryY. Its smallest excess over 5%
across cluster counts was under half a point on every
mechanism.<!-- claim: cl:C12a|GB as=descr; cl:C12c|GB as=similar; cl:C12c|NAT as=similar; cl:C12c|AR1 as=similar -->
Above 45 clusters, Fay-Graubard and the marginalized model were both
within about half a point of 5% on the random-effects
mechanisms.<!-- claim: cl:C25a|NAT as=descr; cl:C25b|NAT as=similar; cl:C25b|AR1 as=similar -->
The ordering of the methods was stable across mechanisms, with Spearman
correlations of 0.83 to 0.98 between pairs of mechanisms. GLMM-AR1
differed from GLMM-native by 0.16 percentage points on average across
methods.<!-- claim: cl:C08|rate as=descr; cl:C09|all as=descr -->

Section 3.5 of the main text compares Fay-Graubard with exchangeable,
independence and AR(1) working correlations over the three study
mechanisms. Per design point, 90% of the differences lay between -1.4
and 0.40 points, and 18 of the 1,896 differences with independence and
11 with AR(1) exceeded 2.5 points in absolute value (at most 5.1 and
3.7).<!-- noclaim: descriptive values -->

**Table S3.** Type I error at the nominal 5% (%): pooled rate, worst
case and mean absolute departure (MAD) over cluster counts (percentage
points), main grid, with the marginalized model without the convergence
check (sensitivity analysis). Only the row of the marginalized model
differs from Table 3. GB = GenBinaryY, NAT = GLMM-native, AR1 =
GLMM-AR1. Largest Monte Carlo standard error: 0.047 points for the
pooled rate, 0.32 for the worst case and 0.047 for the MAD.

|  | Rate GB | Rate NAT | Rate AR1 | Worst GB | Worst NAT | Worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 4.7 | 2.8 | 2.9 | 1.0 | 0.5 | 0.6 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 8.6 | 5.9 | 6.0 | 1.5 | 0.9 | 1.0 |
| Marginalized model | 6.6 | 6.2 | 6.4 | 2.2 | 1.5 | 1.6 | 1.6 | 1.2 | 1.4 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 4.5 | 1.8 | 2.2 | 3.5 | 1.4 | 1.6 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 4.9 | 2.3 | 2.4 | 4.2 | 1.6 | 1.9 |
| GEE, uncorrected | 8.4 | 7.6 | 7.7 | 22.1 | 18.0 | 19.1 | 3.9 | 3.0 | 3.2 |
| GLMM, likelihood ratio | 7.8 | 5.8 | 6.1 | 3.6 | 1.8 | 1.8 | 2.8 | 0.8 | 1.1 |
| GEE, Kauermann-Carroll | 6.9 | 6.3 | 6.4 | 14.1 | 10.7 | 11.2 | 2.5 | 1.7 | 1.8 |

**Table S4.** Mean absolute departure from 5% by number of occasions
(percentage points), cell-weighted, main grid, marginalized model with
the convergence check. Monte Carlo standard errors by block bootstrap;
values given to the digit they support. Largest Monte Carlo standard
error: 0.50 points.

|                        | T = 2 | T = 3 | T ≥ 4 |
|------------------------|-------|-------|-------|
| GEE, Fay-Graubard      | 2.4   | 1.0   | 1.2   |
| GEE, Mancl-DeRouen     | 2.5   | 1.2   | 1.6   |
| Marginalized model     | 1.6   | 1.1   | 0.7   |
| GLMM, Wald             | 5.6   | 4.0   | 2.1   |
| Shrinkage GLM          | 10.8  | 4.8   | 2.2   |
| GEE, uncorrected       | 2.6   | 3.7   | 3.8   |
| GLMM, likelihood ratio | 10.4  | 2.4   | 1.4   |
| GEE, Kauermann-Carroll | 2.4   | 2.0   | 2.5   |

**Table S5.** Method-cells detectably off nominal (%), exact binomial
test with Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY)
corrections within method, main grid, marginalized model with the
convergence check. Monte Carlo standard errors by block bootstrap;
values given to the digit they support. Largest Monte Carlo standard
error: 1.1 points.

|                        | BH  | BY  |
|------------------------|-----|-----|
| GEE, Fay-Graubard      | 23  | 14  |
| GEE, Mancl-DeRouen     | 28  | 20  |
| Marginalized model     | 1   | 0.0 |
| GLMM, Wald             | 51  | 42  |
| Shrinkage GLM          | 53  | 45  |
| GEE, uncorrected       | 46  | 35  |
| GLMM, likelihood ratio | 35  | 28  |
| GEE, Kauermann-Carroll | 38  | 29  |

**Table S6.** Type I error on each method’s own usable replicates and on
the replicates usable for all methods (common set) (%), main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.
Largest Monte Carlo standard error: 0.047 points for the rates and 0.023
for the shift.

|  | Own GB | Own NAT | Own AR1 | Common GB | Common NAT | Common AR1 | Shift GB | Shift NAT | Shift AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 5.6 | 5.4 | 5.3 | -0.2 | 0.0 | 0.0 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 5.8 | 5.5 | 5.5 | -0.2 | 0.0 | -0.1 |
| Marginalized model | 5.3 | 5.1 | 5.3 | 5.2 | 5.1 | 5.2 | -0.1 | 0.0 | 0.0 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 8.6 | 6.3 | 6.5 | 0.0 | 0.0 | 0.0 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 8.6 | 6.3 | 6.5 | -0.7 | -0.3 | -0.4 |
| GEE, uncorrected | 8.4 | 7.6 | 7.7 | 7.5 | 7.1 | 7.1 | -0.9 | -0.5 | -0.6 |
| GLMM, likelihood ratio | 7.8 | 5.8 | 6.1 | 6.9 | 5.0 | 5.4 | -0.9 | -0.8 | -0.8 |
| GEE, Kauermann-Carroll | 6.9 | 6.3 | 6.4 | 6.6 | 6.2 | 6.3 | -0.3 | -0.1 | -0.1 |

**Table S7.** Replicates returning a usable result by number of
occasions (%), under the null, averaged over design points, every
replicate with an estimable design in the denominator, main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.
Largest Monte Carlo standard error: 0.13 points.

|  | GEE, Fay-Graubard | GEE, Mancl-DeRouen | Marginalized model | GLMM, Wald | Shrinkage GLM | GEE, uncorrected | GLMM, likelihood ratio | GEE, Kauermann-Carroll |
|----|----|----|----|----|----|----|----|----|
| Random intercept, T = 2 | 72.7 | 72.7 | 68.2 | 68.6 | 74.0 | 77.8 | 98.9 | 73.0 |
| Random intercept, T = 3 | 88.3 | 88.3 | 87.1 | 87.8 | 88.9 | 91.4 | 99.6 | 88.7 |
| Random intercept, T = 4 | 93.8 | 93.7 | 92.9 | 93.7 | 94.1 | 95.8 | 99.9 | 93.9 |
| Random intercept, T = 6 | 97.7 | 97.6 | 97.4 | 97.8 | 97.9 | 98.6 | 100.0 | 97.7 |
| Random intercept, T = 8 | 98.9 | 98.8 | 98.9 | 99.0 | 99.0 | 99.4 | 100.0 | 98.9 |
| Random intercept, T = 10 | 99.5 | 99.4 | 99.5 | 99.6 | 99.6 | 99.7 | 100.0 | 99.4 |
| Random slopes, T = 2 | 73.1 | 73.1 | 42.2 | 0.0 | 0.0 | 78.0 | 0.1 | 73.2 |
| Random slopes, T = 3 | 88.4 | 88.4 | 69.1 | 56.7 | 59.3 | 91.5 | 66.5 | 88.8 |
| Random slopes, T = 4 | 94.1 | 94.0 | 80.6 | 90.8 | 93.7 | 96.0 | 99.4 | 94.1 |
| Random slopes, T = 6 | 98.0 | 97.9 | 90.4 | 97.7 | 98.1 | 98.8 | 100.0 | 97.9 |
| Random slopes, T = 8 | 99.1 | 99.0 | 94.6 | 99.2 | 99.2 | 99.5 | 100.0 | 99.1 |
| Random slopes, T = 10 | 99.5 | 99.4 | 96.8 | 99.6 | 99.6 | 99.7 | 100.0 | 99.4 |

At 20 clusters or fewer, the likelihood-ratio test returned a usable
result in 86% of attempts and the uncorrected GEE in 78%. The three
corrected sandwiches returned a usable result in 69% to 69% of attempts,
and the marginalized model, which gives no result when the refit of its
convergence check does not converge (Section 2.4 of the main text), in
fewer, 54%.<!-- claim: cl:C14a|all as=larger; cl:C14d|all as=larger -->

**Table S8.** Usable result and unconditional validity at 20 clusters or
fewer (%), under the null, averaged over design points, main grid,
marginalized model with the convergence check. Monte Carlo standard
errors by block bootstrap; values given to the digit they support.
Largest Monte Carlo standard error: 0.096 points for the usable share
and 0.10 for unconditional validity.

|  | Usable | Unconditional validity | Design points with a usable replicate |
|----|----|----|----|
| GEE, Fay-Graubard | 68.9 | 63.8 | 270 |
| GEE, Mancl-DeRouen | 68.7 | 62.3 | 270 |
| Marginalized model | 54.1 | 51.4 | 270 |
| GLMM, Wald | 64.5 | 61.0 | 244 |
| Shrinkage GLM | 66.9 | 62.6 | 251 |
| GEE, uncorrected | 77.6 | 64.1 | 270 |
| GLMM, likelihood ratio | 85.8 | 79.9 | 255 |
| GEE, Kauermann-Carroll | 69.3 | 60.6 | 270 |

**Table S9.** Mean absolute departure from 5% by estimability stratum
(percentage points); strata by the mean usable proportion of the five
marginal methods, main grid, marginalized model with the convergence
check. Monte Carlo standard errors by block bootstrap; values given to
the digit they support. Largest Monte Carlo standard error: 1.2 points.

|                        | \<50% | 50-80% | 80-95% | 95-99.9% | 100% |
|------------------------|-------|--------|--------|----------|------|
| GEE, Fay-Graubard      | 3.8   | 2.8    | 2.2    | 1.1      | 0.6  |
| GEE, Mancl-DeRouen     | 4.3   | 4.2    | 3.1    | 1.3      | 0.6  |
| Marginalized model     | 3.0   | 1.6    | 0.9    | 0.7      | 0.6  |
| GLMM, Wald             | 6     | 2.3    | 2.5    | 3.1      | 2.2  |
| Shrinkage GLM          | 11    | 5      | 3.7    | 3.4      | 2.2  |
| GEE, uncorrected       | 10.4  | 9.3    | 6.4    | 2.8      | 1.1  |
| GLMM, likelihood ratio | 4.3   | 4.6    | 4.5    | 2.7      | 1.2  |
| GEE, Kauermann-Carroll | 6.1   | 6.2    | 4.5    | 2.0      | 0.8  |
| design points          | 93    | 123    | 196    | 398      | 630  |

The size-adjusted critical value was the 95th percentile of
\|estimate/SE\| over the usable null replicates of the same cell; for
the likelihood-ratio test it was the 5th percentile of the null
p-values. It uses the null distribution of each design point, unknown in
practice, and serves only to compare methods at equal level. Table S10
gives raw and size-adjusted power.<!-- noclaim: description --> After
size adjustment, power ranged from 22% to 27% across methods and
mechanisms.<!-- noclaim: range --> On GenBinaryY the Wald test of the
mixed model had more size-adjusted power than Fay-Graubard, by 3.3
points with the default optimizer of lme4 and by 3.7 points with the
optimizer bobyqa (Section S15). On GLMM-AR1 the difference was about two
points with the default optimizer, and the two were similar with bobyqa.
On GLMM-native the two were similar with the default optimizer, and with
bobyqa Fay-Graubard had more power, by about two
points.<!-- claim: cl:C16a|GB as=larger; cl:C16a|GB:bobyqa as=larger; cl:C16a|AR1 as=about1; cl:C16a|AR1:bobyqa as=similar; cl:C16a|NAT as=similar; cl:C16a|NAT:bobyqa as=about1 -->
Part of the size-adjusted power of the Wald test of the mixed model came
from fits with a collapsed standard error
($|\hat\beta/\mathrm{SE}| > 37$, a p-value numerically zero), which the
usability limits keep: averaged over the 1,254 cells of the main grid
with a size-adjusted value, they made up 2.3 of its 24 percentage
points, and in 193 cells the size-adjusted critical value itself
exceeded 37.<!-- claim: reg:pb2_mean as=descr -->

**Table S10.** Raw and size-adjusted power at the interaction contrast
(%), main grid, pooled over design points, by mechanism, marginalized
model with the convergence check; Kauermann-Carroll computed with
*fewclust*. Size adjustment uses the 95th percentile of \|estimate/SE\|
under the null of the same cell (for the likelihood-ratio test, the 5th
percentile of the null p-values), and needs at least two usable null
replicates. Monte Carlo standard errors by block bootstrap, at most 0.26
points; values given to the digit they support.

|                                     | Raw  | Size-adjusted |
|-------------------------------------|------|---------------|
| GenBinaryY, GEE, Fay-Graubard       | 23.3 | 22.0          |
| GenBinaryY, GEE, Mancl-DeRouen      | 23.5 | 21.9          |
| GenBinaryY, Marginalized model      | 25.3 | 24.8          |
| GenBinaryY, GLMM, Wald              | 31.9 | 25.3          |
| GenBinaryY, Shrinkage GLM           | 32.2 | 24.9          |
| GenBinaryY, GEE, uncorrected        | 26.3 | 21.5          |
| GenBinaryY, GLMM, likelihood ratio  | 30.1 | 26.4          |
| GenBinaryY, GEE, Kauermann-Carroll  | 24.8 | 21.9          |
| GLMM-native, GEE, Fay-Graubard      | 23.7 | 23.1          |
| GLMM-native, GEE, Mancl-DeRouen     | 23.9 | 23.1          |
| GLMM-native, Marginalized model     | 25.8 | 25.7          |
| GLMM-native, GLMM, Wald             | 28.0 | 23.9          |
| GLMM-native, Shrinkage GLM          | 28.1 | 23.7          |
| GLMM-native, GEE, uncorrected       | 26.6 | 22.8          |
| GLMM-native, GLMM, likelihood ratio | 26.8 | 26.3          |
| GLMM-native, GEE, Kauermann-Carroll | 25.1 | 23.1          |
| GLMM-AR1, GEE, Fay-Graubard         | 24.1 | 23.5          |
| GLMM-AR1, GEE, Mancl-DeRouen        | 24.4 | 23.4          |
| GLMM-AR1, Marginalized model        | 26.7 | 26.2          |
| GLMM-AR1, GLMM, Wald                | 29.2 | 25.8          |
| GLMM-AR1, Shrinkage GLM             | 29.3 | 25.6          |
| GLMM-AR1, GEE, uncorrected          | 27.0 | 23.0          |
| GLMM-AR1, GLMM, likelihood ratio    | 27.9 | 26.7          |
| GLMM-AR1, GEE, Kauermann-Carroll    | 25.6 | 23.3          |

Under the t reference, Fay-Graubard with exchangeable working
correlation had a smaller worst-case departure than the marginalized
model with the check: 0.9 to 1.1 points against 3.3 to
3.6.<!-- claim: cl:C17d|GB as=larger; cl:C17d|NAT as=larger; cl:C17d|AR1 as=larger -->

**Table S11.** Type I error under a normal and a t reference (%), shift
in percentage points, main grid, marginalized model with the convergence
check. Monte Carlo standard errors by block bootstrap; values given to
the digit they support. Largest Monte Carlo standard error: 0.047 points
for the rates and 0.035 for the differences.

|  | Normal GB | Normal NAT | Normal AR1 | t GB | t NAT | t AR1 | t minus normal GB | t minus normal NAT | t minus normal AR1 | Normal minus t at N ≤ 25 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 5.1 | 4.6 | 4.6 | -0.8 | -0.7 | -0.7 | 2.2 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 5.3 | 4.7 | 4.8 | -0.8 | -0.8 | -0.8 | 2.4 |
| Marginalized model | 5.3 | 5.1 | 5.3 | 4.7 | 4.4 | 4.6 | -0.7 | -0.7 | -0.7 | 2.0 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 8.0 | 5.8 | 6.0 | -0.7 | -0.5 | -0.6 | 1.4 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 8.6 | 6.1 | 6.3 | -0.7 | -0.5 | -0.6 | 1.5 |
| GEE, uncorrected | 8.4 | 7.6 | 7.7 | 7.2 | 6.5 | 6.6 | -1.1 | -1.1 | -1.1 | 3.8 |
| GEE, Kauermann-Carroll | 6.9 | 6.3 | 6.4 | 6.0 | 5.4 | 5.5 | -0.9 | -0.9 | -0.9 | 3.1 |

**Table S12.** Coverage of the 95% confidence interval against each
method’s own estimand (%); not available (n.a.) for the GLMM on
GenBinaryY, main grid, marginalized model with the convergence check.
Monte Carlo standard errors by block bootstrap; values given to the
digit they support. Largest Monte Carlo standard error: 0.049 points.

|                        | GenBinaryY | GLMM-native | GLMM-AR1 |
|------------------------|------------|-------------|----------|
| GEE, Fay-Graubard      | 94.6       | 94.9        | 94.9     |
| GEE, Mancl-DeRouen     | 94.2       | 94.7        | 94.6     |
| Marginalized model     | 94.6       | 94.8        | 94.7     |
| GLMM, Wald             | n.a.       | 93.5        | 93.2     |
| Shrinkage GLM          | 90.1       | 92.9        | 92.6     |
| GEE, uncorrected       | 91.9       | 92.5        | 92.5     |
| GEE, Kauermann-Carroll | 93.4       | 93.8        | 93.8     |

**Table S13.** Pooled type I error (%) on the main grid with fits that
returned warnings kept (U0, primary) and excluded (U1),
replicate-weighted over the null cells. Binomial Monte Carlo standard
error at most 0.044 points.

| Method                         | GB, U0 | GB, U1 | NAT, U0 | NAT, U1 | AR1, U0 | AR1, U1 |
|--------------------------------|--------|--------|---------|---------|---------|---------|
| GEE, Fay-Graubard              | 5.8    | 5.8    | 5.3     | 5.3     | 5.3     | 5.3     |
| GEE, Mancl-DeRouen             | 6.1    | 6.1    | 5.5     | 5.5     | 5.6     | 5.6     |
| Marginalized model, with check | 5.3    | 5.2    | 5.1     | 5.0     | 5.3     | 5.2     |
| GLMM, Wald                     | 8.6    | 5.8    | 6.3     | 4.1     | 6.6     | 4.4     |
| Shrinkage GLM                  | 9.3    | 6.3    | 6.6     | 4.4     | 6.9     | 4.6     |
| GEE, uncorrected               | 8.4    | 8.2    | 7.6     | 7.5     | 7.7     | 7.6     |
| GLMM, likelihood ratio         | 7.8    | 7.2    | 5.8     | 5.3     | 6.1     | 5.6     |
| GEE, Kauermann-Carroll         | 6.9    | 6.9    | 6.3     | 6.3     | 6.4     | 6.4     |

**Table S14.** Scenarios M6-M7 (96 design points each): design points
clearly above 7.5% and clearly below 2.5% under the null, largest
rejection rate (%) over the design points with at least 20 usable null
replicates, and replicates returning a usable result under the null (%),
averaged over design points. The GEE rows use the exchangeable working
correlation with the corrections computed by *geesmv*, as in Table 3, so
their counts can differ slightly from those of Table 5 and Figure S1,
which use *fewclust*. Monte Carlo standard errors are at most 1.3 points
for the largest rate (binomial, at its design point) and 0.57 for the
usable share (upper bound: mean over design points of the binomial
standard error).

| Scenario | Strategy | Above 7.5% | Below 2.5% | Largest rate (%) | Usable (%) |
|----|----|----|----|----|----|
| GenBinaryY, slope 0.1 | Marginalized model with check, t | 0 | 2 | 7.0 | 85.3 |
|  | Marginalized model with check, normal | 0 | 0 | 7.1 | 85.3 |
|  | Marginalized model without check, t | 0 | 1 | 7.5 | 88.8 |
|  | GEE, Fay-Graubard, t | 0 | 4 | 9.1 | 91.2 |
|  | GEE, Mancl-DeRouen, t | 8 | 4 | 13.3 | 91.2 |
|  | GLMM, likelihood ratio | 0 | 1 | 7.3 | 83.3 |
| GenBinaryY, slope 0.25 | Marginalized model with check, t | 0 | 1 | 7.6 | 84.4 |
|  | Marginalized model with check, normal | 1 | 0 | 9.3 | 84.4 |
|  | Marginalized model without check, t | 1 | 0 | 9.4 | 88.3 |
|  | GEE, Fay-Graubard, t | 3 | 4 | 10.4 | 90.8 |
|  | GEE, Mancl-DeRouen, t | 10 | 4 | 14.1 | 90.8 |
|  | GLMM, likelihood ratio | 0 | 1 | 7.5 | 83.2 |
| GLMM-native, slope variance 0.01 | Marginalized model with check, t | 0 | 1 | 6.9 | 86.6 |
|  | Marginalized model with check, normal | 0 | 0 | 7.2 | 86.6 |
|  | Marginalized model without check, t | 2 | 1 | 10.0 | 90.2 |
|  | GEE, Fay-Graubard, t | 0 | 4 | 7.0 | 92.3 |
|  | GEE, Mancl-DeRouen, t | 0 | 4 | 9.3 | 92.2 |
|  | GLMM, likelihood ratio | 0 | 0 | 6.1 | 83.3 |
| GLMM-native, slope variance 0.04 | Marginalized model with check, t | 0 | 4 | 7.8 | 87.0 |
|  | Marginalized model with check, normal | 1 | 0 | 10.3 | 87.0 |
|  | Marginalized model without check, t | 6 | 0 | 11.1 | 90.4 |
|  | GEE, Fay-Graubard, t | 0 | 3 | 6.7 | 92.3 |
|  | GEE, Mancl-DeRouen, t | 0 | 4 | 9.0 | 92.3 |
|  | GLMM, likelihood ratio | 0 | 0 | 6.6 | 83.3 |

**Table S15.** Clinical application: outcome of the fits in 500
subsamples of whole patients at each subsample size N (% of subsamples).
Not returned: no point estimate; no SE: estimate without a standard
error; extreme: \|estimate\| or SE above 10; usable: the remainder. For
each method and N the four columns add up to all subsamples. Binomial
Monte Carlo standard errors over the 500 subsamples are at most 2.2
points.

| Method                 | N   | Not returned | No SE | Extreme | Usable |
|------------------------|-----|--------------|-------|---------|--------|
| GLMM                   | 10  | 21.8         | 0.0   | 0.0     | 78.2   |
|                        | 15  | 14.6         | 0.0   | 0.0     | 85.4   |
|                        | 20  | 5.4          | 0.0   | 0.0     | 94.6   |
|                        | 25  | 3.4          | 0.0   | 0.0     | 96.6   |
|                        | 50  | 0.0          | 0.0   | 0.0     | 100.0  |
| Marginalized model     | 10  | 22.8         | 0.0   | 0.0     | 77.2   |
|                        | 15  | 14.6         | 0.0   | 0.0     | 85.4   |
|                        | 20  | 5.4          | 0.0   | 0.0     | 94.6   |
|                        | 25  | 3.4          | 0.0   | 0.0     | 96.6   |
|                        | 50  | 0.0          | 0.0   | 0.0     | 100.0  |
| GEE, uncorrected       | 10  | 24.0         | 0.0   | 6.0     | 70.0   |
|                        | 15  | 14.8         | 0.0   | 4.4     | 80.8   |
|                        | 20  | 5.6          | 0.0   | 3.4     | 91.0   |
|                        | 25  | 3.4          | 0.0   | 1.2     | 95.4   |
|                        | 50  | 0.0          | 0.0   | 0.0     | 100.0  |
| GEE, Fay-Graubard      | 10  | 24.0         | 41.4  | 0.4     | 34.2   |
|                        | 15  | 14.8         | 23.8  | 0.2     | 61.2   |
|                        | 20  | 5.6          | 14.2  | 0.2     | 80.0   |
|                        | 25  | 3.4          | 8.0   | 0.2     | 88.4   |
|                        | 50  | 0.0          | 0.2   | 0.0     | 99.8   |
| GEE, Kauermann-Carroll | 10  | 24.0         | 63.6  | 1.2     | 11.2   |
|                        | 15  | 14.8         | 49.0  | 1.4     | 34.8   |
|                        | 20  | 5.6          | 31.4  | 1.8     | 61.2   |
|                        | 25  | 3.4          | 16.0  | 1.0     | 79.6   |
|                        | 50  | 0.0          | 1.2   | 0.0     | 98.8   |
| GEE, Mancl-DeRouen     | 10  | 24.0         | 60.4  | 0.6     | 15.0   |
|                        | 15  | 14.8         | 47.4  | 0.6     | 37.2   |
|                        | 20  | 5.6          | 31.8  | 0.4     | 62.2   |
|                        | 25  | 3.4          | 16.6  | 0.0     | 80.0   |
|                        | 50  | 0.0          | 1.2   | 0.0     | 98.8   |

## S4. Numerical behaviour of the small-sample corrections

The corrections of Mancl and DeRouen, Kauermann and Carroll and Fay and
Graubard were computed with two implementations. *geesmv* was used for
the exchangeable working correlation. On sparse binary data with few
clusters its three functions can loop indefinitely in compiled code, so
each call ran in a separate process killed after 5 seconds; the largest
share of time-outs at a design point under one hypothesis was
20%.<!-- noclaim: descriptive values --> *geesmv* also takes the
clusters from any column named *id*, whatever cluster variable is
passed: in the simulations and in the cohort the cluster column has that
name, and in the respiratory trial it was renamed before every call.
With AR(1) working correlation *geesmv* fails as soon as one subject has
a single observation, so it was not used with that correlation.

*fewclust* (Table S2) computes the three corrections from the
estimating-equation building blocks. On 666 or more datasets for each
working correlation and correction on which both computations returned a
standard error, its standard errors agreed with a computation by hand;
the largest difference over the three corrections was below
$10^{-12}$.<!-- noclaim: descriptive values --> Against *geesmv* on the
same exchangeable fits of the first 50 replicates of each design point
and hypothesis, the subsample on which *geesmv* was run in this
comparison (174,804 datasets), Mancl-DeRouen and Fay-Graubard agreed to
within one per cent in 98% and 97% of datasets, and Kauermann-Carroll in
only 2.0%. *geesmv::GEE.var.kc* rebuilds $(I - H_{ii})^{-1/2}$ with the
transpose of the eigenvectors of a non-symmetric matrix, an identity
that holds only for symmetric matrices. With the published formula,
Kauermann-Carroll with exchangeable correlation and a t reference had a
mean type I error of 5.6% over the 1,896 design points of the three
study mechanisms, against 9.7% with *geesmv*, and it was clearly above
7.5% at 146 instead of 963 design points. With *geesmv* the correction
appears to shrink rather than inflate the variance; this is a property
of the implementation, not of the
estimator.<!-- noclaim: interpretation --> The CR2 estimator of
*clubSandwich* is the Bell-McCaffrey version of Kauermann-Carroll; it
was used here with a t reference on N - p degrees of freedom, not with
the Satterthwaite-type degrees of freedom for which it was proposed
(reference \[27\] of the main text), and with those its results could
differ.

With exchangeable correlation, Mancl-DeRouen and Fay-Graubard were also
computed with *geesmv*, under a hard time limit of 5 seconds; a time-out
counted as an unusable fit, and these values enter Table 3, Table 4 and
Table S10. Every Kauermann-Carroll result is from *fewclust*, because
*geesmv* applies a symmetric matrix square root to a non-symmetric
matrix (above). Computed with *fewclust* instead of *geesmv*,
Fay-Graubard with exchangeable correlation rejected a true null at
similar rates (difference 0.061 points) and had a similar usable share
(difference 0.21
points).<!-- claim: cv:impl_fg as=similar; cv:impl_fg_pu as=similar -->
The Kauermann-Carroll correction computed from its published formula
rejected a true null 4.1 points less often than the *geesmv* version on
average under a t reference.<!-- claim: cv:impl_kc as=larger --> Its
rejection rates in Table 3 and Table S10 are therefore lower than those
obtained with *geesmv*.<!-- noclaim: consequence of the correction -->

## S5. Tables of the scenario anchored to the cohort (M1c)

In scenario M1c, the size-adjusted power of the Wald test of the mixed
model fell below that of Fay-Graubard, by 4.1 to 4.8 points with the
default optimizer. With bobyqa the deficit was about two points on
GenBinaryY and GLMM-AR1, and the two were similar on GLMM-native. The
default optimizer returned degenerate fits with a near-zero standard
error.<!-- claim: cm:C16a|GB as=larger; cm:C16a|NAT as=larger; cm:C16a|AR1 as=larger; cm:C16a|GB:bobyqa as=about1; cm:C16a|NAT:bobyqa as=similar; cm:C16a|AR1:bobyqa as=about1 -->
In M1c, the method with the fewest miscalibrated cells could not be
identified.<!-- claim: cm:C10b|BH as=notdist -->

**Table S16.** Type I error at the nominal 5% (%): pooled rate, worst
case and mean absolute departure (MAD) over cluster counts (percentage
points), scenario M1c (72 design points). Largest Monte Carlo standard
error: 0.22 points for the pooled rate, 0.65 for the worst case and 0.22
for the MAD.

|  | Rate GB | Rate NAT | Rate AR1 | Worst GB | Worst NAT | Worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.1 | 5.0 | 5.0 | 1.8 | 1.2 | 1.8 | 0.6 | 0.5 | 0.5 |
| GEE, Mancl-DeRouen | 4.5 | 4.5 | 4.4 | 1.2 | 1.1 | 1.1 | 0.8 | 0.7 | 0.8 |
| Marginalized model | 7.0 | 6.8 | 6.9 | 2.4 | 2.4 | 2.9 | 2.0 | 1.8 | 2.0 |
| GLMM, Wald | 7.7 | 6.0 | 6.2 | 4.5 | 2.8 | 2.5 | 2.5 | 1.3 | 1.5 |
| Shrinkage GLM | 11.2 | 7.9 | 8.3 | 6.6 | 3.3 | 4.5 | 6.2 | 2.9 | 3.3 |
| GEE, uncorrected | 8.3 | 8.1 | 8.2 | 10 | 9 | 10 | 3.6 | 3.3 | 3.5 |
| GLMM, likelihood ratio | 10.5 | 7.5 | 8.0 | 7.0 | 4.8 | 5.6 | 5.5 | 2.5 | 3.0 |
| GEE, Kauermann-Carroll | 6.1 | 6.0 | 6.0 | 5 | 4 | 4 | 1.3 | 1.1 | 1.2 |

**Table S17.** Mean absolute departure from 5% by number of occasions
(percentage points), cell-weighted, scenario M1c (72 design points).
Largest Monte Carlo standard error: 0.28 points.

|                        | T = 2 | T = 3 | T ≥ 4 |
|------------------------|-------|-------|-------|
| GEE, Fay-Graubard      | 2.5   | 0.7   | 1.2   |
| GEE, Mancl-DeRouen     | 3.0   | 1.0   | 0.9   |
| Marginalized model     | 1.9   | 1.8   | 1.9   |
| GLMM, Wald             | 1.5   | 1.3   | 3.6   |
| Shrinkage GLM          | 7.1   | 2.6   | 3.9   |
| GEE, uncorrected       | 0.9   | 3.1   | 4.7   |
| GLMM, likelihood ratio | 6.7   | 3.2   | 2.5   |
| GEE, Kauermann-Carroll | 2.1   | 1.1   | 2.3   |

**Table S18.** Method-cells detectably off nominal (%), exact binomial
test with Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY)
corrections within method, scenario M1c (72 design points). Largest
Monte Carlo standard error: 5.0 points.

|                        | BH  | BY  |
|------------------------|-----|-----|
| GEE, Fay-Graubard      | 28  | 25  |
| GEE, Mancl-DeRouen     | 26  | 19  |
| Marginalized model     | 57  | 42  |
| GLMM, Wald             | 54  | 38  |
| Shrinkage GLM          | 78  | 72  |
| GEE, uncorrected       | 61  | 53  |
| GLMM, likelihood ratio | 67  | 64  |
| GEE, Kauermann-Carroll | 44  | 36  |

**Table S19.** Type I error on each method’s own usable replicates and
on the replicates usable for all methods (common set) (%), scenario M1c
(72 design points). Largest Monte Carlo standard error: 0.22 points for
the rates and 0.13 for the shift.

|  | Own GB | Own NAT | Own AR1 | Common GB | Common NAT | Common AR1 | Shift GB | Shift NAT | Shift AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.1 | 5.0 | 5.0 | 4.7 | 4.7 | 4.7 | -0.3 | -0.3 | -0.3 |
| GEE, Mancl-DeRouen | 4.5 | 4.5 | 4.4 | 4.2 | 4.2 | 4.2 | -0.2 | -0.2 | -0.2 |
| Marginalized model | 7.0 | 6.8 | 6.9 | 6.4 | 6.2 | 6.4 | -0.6 | -0.6 | -0.5 |
| GLMM, Wald | 7.7 | 6.0 | 6.2 | 7.7 | 6.0 | 6.3 | 0.0 | 0.0 | 0.1 |
| Shrinkage GLM | 11.2 | 7.9 | 8.3 | 7.7 | 6.0 | 6.3 | -3.5 | -1.9 | -2.0 |
| GEE, uncorrected | 8.3 | 8.1 | 8.2 | 7.5 | 7.4 | 7.5 | -0.8 | -0.7 | -0.7 |
| GLMM, likelihood ratio | 10.5 | 7.5 | 8.0 | 7.3 | 4.9 | 5.3 | -3.3 | -2.6 | -2.6 |
| GEE, Kauermann-Carroll | 6.1 | 6.0 | 6.0 | 5.7 | 5.6 | 5.7 | -0.4 | -0.3 | -0.3 |

**Table S20.** Replicates returning a usable result by number of
occasions (%), under the null, averaged over design points, every
replicate with an estimable design in the denominator, scenario M1c (72
design points). Largest Monte Carlo standard error: 0.35 points.

|  | GEE, Fay-Graubard | GEE, Mancl-DeRouen | Marginalized model | GLMM, Wald | Shrinkage GLM | GEE, uncorrected | GLMM, likelihood ratio | GEE, Kauermann-Carroll |
|----|----|----|----|----|----|----|----|----|
| Random intercept, T = 2 | 81.3 | 81.3 | 69.6 | 68.8 | 82.9 | 84.1 | 99.9 | 79.9 |
| Random intercept, T = 3 | 93.9 | 93.8 | 89.1 | 90.8 | 94.6 | 95.1 | 100.0 | 93.1 |
| Random intercept, T = 4 | 96.6 | 96.4 | 94.3 | 96.0 | 97.0 | 97.5 | 100.0 | 96.1 |
| Random intercept, T = 6 | 98.4 | 98.1 | 97.5 | 98.6 | 98.7 | 98.9 | 100.0 | 98.1 |

**Table S21.** Usable result and unconditional validity at 20 clusters
or fewer (%), under the null, averaged over design points, scenario M1c
(72 design points). Largest Monte Carlo standard error: 0.26 points for
the usable share and 0.28 for unconditional validity.

|  | Usable | Unconditional validity | Design points with a usable replicate |
|----|----|----|----|
| GEE, Fay-Graubard | 85.9 | 81.5 | 36 |
| GEE, Mancl-DeRouen | 85.5 | 81.7 | 36 |
| Marginalized model | 78.8 | 73.3 | 36 |
| GLMM, Wald | 79.8 | 75.6 | 36 |
| Shrinkage GLM | 87.2 | 79.3 | 36 |
| GEE, uncorrected | 88.6 | 79.7 | 36 |
| GLMM, likelihood ratio | 99.9 | 90.0 | 36 |
| GEE, Kauermann-Carroll | 84.5 | 78.7 | 36 |

**Table S22.** Mean absolute departure from 5% by estimability stratum
(percentage points); strata by the mean usable proportion of the five
marginal methods, scenario M1c (72 design points). Largest Monte Carlo
standard error: 0.91 points.

|                        | \<50% | 50-80% | 80-95% | 95-99.9% | 100% |
|------------------------|-------|--------|--------|----------|------|
| GEE, Fay-Graubard      | 4.7   | 2.3    | 2.4    | 0.8      | 0.6  |
| GEE, Mancl-DeRouen     | 5.0   | 2.3    | 2.4    | 0.8      | 0.6  |
| Marginalized model     | 1     | 1.6    | 2.5    | 2.1      | 1.1  |
| GLMM, Wald             | 1     | 1.9    | 1.5    | 2.2      | 5.0  |
| Shrinkage GLM          | 8     | 9      | 4.2    | 2.9      | 5.0  |
| GEE, uncorrected       | 1     | 6      | 5.4    | 2.8      | 1.8  |
| GLMM, likelihood ratio | 7     | 8.1    | 5.5    | 2.7      | 1.1  |
| GEE, Kauermann-Carroll | 4.3   | 3.8    | 3.3    | 1.2      | 0.9  |
| design points          | 3     | 7      | 15     | 33       | 14   |

**Table S23.** Raw and size-adjusted power at the interaction contrast
(%), scenario M1c (72 design points). Largest Monte Carlo standard
error: 0.24 points for raw power and 0.29 for size-adjusted power.

|                                     | Raw  | Size-adjusted |
|-------------------------------------|------|---------------|
| GenBinaryY, GEE, Fay-Graubard       | 9.7  | 9.4           |
| GenBinaryY, GEE, Mancl-DeRouen      | 8.6  | 9.2           |
| GenBinaryY, Marginalized model      | 11.2 | 8.7           |
| GenBinaryY, GLMM, Wald              | 10.6 | 4.8           |
| GenBinaryY, Shrinkage GLM           | 14.2 | 4.8           |
| GenBinaryY, GEE, uncorrected        | 13.4 | 9.0           |
| GenBinaryY, GLMM, likelihood ratio  | 15.7 | 8.9           |
| GenBinaryY, GEE, Kauermann-Carroll  | 10.8 | 9.2           |
| GLMM-native, GEE, Fay-Graubard      | 9.5  | 9.1           |
| GLMM-native, GEE, Mancl-DeRouen     | 8.6  | 9.0           |
| GLMM-native, Marginalized model     | 11.1 | 8.8           |
| GLMM-native, GLMM, Wald             | 8.6  | 5.1           |
| GLMM-native, Shrinkage GLM          | 10.4 | 5.0           |
| GLMM-native, GEE, uncorrected       | 13.1 | 8.9           |
| GLMM-native, GLMM, likelihood ratio | 11.8 | 9.1           |
| GLMM-native, GEE, Kauermann-Carroll | 10.7 | 9.0           |
| GLMM-AR1, GEE, Fay-Graubard         | 9.8  | 9.7           |
| GLMM-AR1, GEE, Mancl-DeRouen        | 8.7  | 9.5           |
| GLMM-AR1, Marginalized model        | 11.2 | 8.6           |
| GLMM-AR1, GLMM, Wald                | 8.9  | 4.9           |
| GLMM-AR1, Shrinkage GLM             | 11.0 | 5.0           |
| GLMM-AR1, GEE, uncorrected          | 13.5 | 9.3           |
| GLMM-AR1, GLMM, likelihood ratio    | 12.5 | 9.1           |
| GLMM-AR1, GEE, Kauermann-Carroll    | 10.9 | 9.4           |

**Table S24.** Type I error under a normal and a t reference (%), shift
in percentage points, scenario M1c (72 design points). Largest Monte
Carlo standard error: 0.22 points for the rates and 0.10 for the
differences.

|  | Normal GB | Normal NAT | Normal AR1 | t GB | t NAT | t AR1 | t minus normal GB | t minus normal NAT | t minus normal AR1 | Normal minus t at N ≤ 25 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.1 | 5.0 | 5.0 | 3.6 | 3.6 | 3.5 | -1.5 | -1.4 | -1.5 | 1.9 |
| GEE, Mancl-DeRouen | 4.5 | 4.5 | 4.4 | 3.2 | 3.2 | 3.2 | -1.3 | -1.3 | -1.2 | 1.6 |
| Marginalized model | 7.0 | 6.8 | 6.9 | 5.1 | 5.0 | 5.1 | -1.9 | -1.8 | -1.9 | 2.4 |
| GLMM, Wald | 7.7 | 6.0 | 6.2 | 6.6 | 5.1 | 5.2 | -1.1 | -0.9 | -1.0 | 1.1 |
| Shrinkage GLM | 11.2 | 7.9 | 8.3 | 9.7 | 6.7 | 7.0 | -1.5 | -1.2 | -1.3 | 1.6 |
| GEE, uncorrected | 8.3 | 8.1 | 8.2 | 6.0 | 5.9 | 5.8 | -2.3 | -2.2 | -2.4 | 3.1 |
| GEE, Kauermann-Carroll | 6.1 | 6.0 | 6.0 | 4.3 | 4.2 | 4.2 | -1.9 | -1.7 | -1.8 | 2.4 |

**Table S25.** Coverage of the 95% confidence interval against each
method’s own estimand (%); not available (n.a.) for the GLMM on
GenBinaryY, scenario M1c (72 design points). Largest Monte Carlo
standard error: 0.20 points.

|                        | GenBinaryY | GLMM-native | GLMM-AR1 |
|------------------------|------------|-------------|----------|
| GEE, Fay-Graubard      | 93.9       | 94.5        | 94.5     |
| GEE, Mancl-DeRouen     | 94.7       | 95.2        | 95.1     |
| Marginalized model     | 93.0       | 93.2        | 93.1     |
| GLMM, Wald             | n.a.       | 95.1        | 94.9     |
| Shrinkage GLM          | 89.4       | 93.0        | 92.4     |
| GEE, uncorrected       | 90.2       | 91.2        | 91.2     |
| GEE, Kauermann-Carroll | 92.9       | 93.6        | 93.5     |

## S6. Degrees of freedom of the t reference

Degrees of freedom N - 2, as used by Li and Redden \[15\], were examined
as a sensitivity analysis. Averaging the Mancl-DeRouen and
Kauermann-Carroll variances, as studied by Ford and Westgate \[13\], did
not give a strategy that was never clearly above 7.5%. With N - p
degrees of freedom, GEE with AR(1) working correlation and the averaged
variance was clearly above 7.5% at 56 of the 2,376 design points of the
four mechanisms, and GEE with exchangeable working correlation and the
averaged variance at 121 (Table S26 gives the three study
mechanisms).<!-- claim: reg:mdkc as=descr -->

**Table S26.** The decision rule with N - p and with N - 2 degrees of
freedom, as used by Li and Redden, on the 1,896 design points of the
three study mechanisms (the fourth mechanism was analysed with N - p
only). Above: design points clearly above 7.5%; below: clearly below
2.5%; usable: replicates returning a usable result under the null (%),
averaged over design points. With both choices the rule selects the
marginalized model with the convergence check; GEE with AR(1) working
correlation and Fay-Graubard is clearly above 7.5% at 2 design points
with N - p and at 2 with N - 2. The rows for the average of the
Mancl-DeRouen and Kauermann-Carroll variances use the same usability
rule, denominators and exact intervals, recomputed from the stored
per-replicate results. The Monte Carlo standard error of the usable
share is at most 0.55 points (upper bound: mean over design points of
the binomial standard error).

| Strategy | Above, N - p | Below, N - p | Below at ten clusters, N - p | Above, N - 2 | Below, N - 2 | Below at ten clusters, N - 2 | Usable (%) |
|----|----|----|----|----|----|----|----|
| GEE AR(1), Fay-Graubard | 2 | 89 | 5 | 2 | 86 | 4 | 92.0 |
| Marginalized model with check | 0 | 39 | 30 | 0 | 21 | 13 | 84.2 |
| GEE independence, Fay-Graubard | 4 | 84 | 7 | 4 | 80 | 5 | 92.2 |
| GEE exchangeable, Fay-Graubard | 7 | 87 | 3 | 12 | 85 | 3 | 92.1 |
| GEE AR(1), Mancl-DeRouen | 13 | 85 | 4 | 31 | 82 | 3 | 92.0 |
| GEE exchangeable, Mancl-DeRouen | 35 | 83 | 3 | 66 | 82 | 3 | 92.0 |
| GEE exchangeable, Kauermann-Carroll | 146 | 59 | 3 | 180 | 55 | 3 | 91.9 |
| GEE AR(1), average of Mancl-DeRouen and Kauermann-Carroll | 41 | 74 | 3 | 69 | 71 | 3 | 91.8 |
| GEE exchangeable, average of Mancl-DeRouen and Kauermann-Carroll | 80 | 74 | 3 | 108 | 71 | 3 | 91.9 |
| GLMM, Wald | 292 | 15 | 12 | 292 | 10 | 8 | 80.7 |

## S7. Sign-flip score test and the interval obtained by inverting it

The model without the interaction was fitted as a logistic regression
with the mean model of the scenario, and the effective score
contribution of each subject to the interaction was obtained by
projecting out the other coefficients. The p-value compares the observed
sum of contributions with its distribution under random changes of sign:
all sign assignments up to 12 subjects, and otherwise 1,999 random
assignments plus the identity. The test gives a p-value only; the
interval obtained by inverting it is described below.

The sign-flip test of the main text tests the null of no interaction.
Testing instead the null that the interaction equals a value $b_0$, by
entering $b_0 \cdot t \cdot \mathrm{trt}$ as an offset in the model
without the interaction, and collecting the values that are not rejected
at 5%, gives a confidence interval. Its ends were searched from the
independence GEE estimate by a walk in fixed steps followed by
bisection, with the same sign vectors at every value and every fit
started from the default initial value; a value at which the test was
not computable counted as not rejected, so the interval can be
unbounded. Under the alternative the test evaluated at the true value
kept its level: it covered the true value in 95% of datasets over the
2,376 design points, and its non-coverage was clearly outside Bradley’s
band at no design point. The interval, which counts as covering when it
is unbounded, covered the true value in 95% of datasets, against 96% for
GEE with AR(1) working correlation, Fay-Graubard and a t reference; its
non-coverage was clearly below 2.5% at 97 design points under the
alternative and 112 under the null, 90 of the former with two occasions,
where it is almost always unbounded. Both ends were finite in 79% of
datasets, and in 14% with two or three occasions and at most 25 clusters
(Table S27). Where both intervals were finite it was wider than the GEE
interval by 17% on average over design points (geometric mean of the
ratio of widths at each design point, then over design points), and by
6.4% with 75 clusters or more; over datasets, which gives more weight to
the larger designs, the geometric mean was 8.7% (5.6% with 75 clusters
or more).<!-- claim: reg:ci as=descr -->

With the rules of Section 2.4 of the main text the test was not
computable in some of the small designs, often because one arm had only
zero outcomes. With two or three occasions and at most 25 clusters it
answered in 72% of datasets, against 99.9% when the rules were not
applied, with a mean type I error of 5.2%. A version that also excluded
fits with a coefficient of at least 10 is a sensitivity analysis. With
this coefficient filter, on the 1,896 design points of the three study
mechanisms, it rejected a true null in 5.5% of datasets on average, was
clearly above 7.5% at 0 design points and answered in 97% of
datasets.<!-- claim: reg:flip as=descr -->

**Table S27.** Interval obtained by inverting the sign-flip test, under
the alternative, by number of occasions T and number of clusters N:
datasets with an estimable design, datasets with both ends of the
interval finite (%), and datasets in which GEE with AR(1) working
correlation, Fay-Graubard and a t reference returned a usable interval
(%). Binomial Monte Carlo standard errors are at most 0.32 points for
finite intervals and 0.30 for GEE; they are approximate, because the
datasets of a row share replicate seeds across numbers of clusters.

| Occasions | Clusters | Datasets | Interval finite (%) | GEE AR(1), Fay-Graubard usable (%) |
|----|----|----|----|----|
| T = 2 | N = 10 | 26,343 | 0.9 | 10.2 |
|  | N = 15-25 | 80,898 | 6.2 | 33.3 |
|  | N = 35-65 | 102,000 | 22.7 | 79.6 |
|  | N = 75-145 | 192,000 | 72.0 | 98.6 |
| T = 3 | N = 10 | 26,794 | 6.7 | 36.4 |
|  | N = 15-25 | 80,977 | 27.6 | 71.4 |
|  | N = 35-65 | 102,000 | 72.7 | 97.3 |
|  | N = 75-145 | 192,000 | 97.3 | 100.0 |
| T = 4 | N = 10 | 26,883 | 15.6 | 56.5 |
|  | N = 15-25 | 80,996 | 51.6 | 87.1 |
|  | N = 35-65 | 102,000 | 90.9 | 99.4 |
|  | N = 75-145 | 192,000 | 99.7 | 100.0 |
| T = 6 | N = 10 | 26,922 | 34.6 | 80.6 |
|  | N = 15-25 | 80,996 | 78.0 | 96.6 |
|  | N = 35-65 | 102,000 | 98.8 | 99.9 |
|  | N = 75-145 | 192,000 | 100.0 | 100.0 |
| T = 8 | N = 10 | 23,930 | 47.0 | 89.5 |
|  | N = 15-25 | 72,000 | 87.5 | 98.7 |
|  | N = 35-65 | 96,000 | 99.8 | 100.0 |
|  | N = 75-145 | 192,000 | 100.0 | 100.0 |
| T = 10 | N = 10 | 23,944 | 55.9 | 93.9 |
|  | N = 15-25 | 71,996 | 92.0 | 99.4 |
|  | N = 35-65 | 96,000 | 99.9 | 100.0 |
|  | N = 75-145 | 192,000 | 100.0 | 100.0 |

## S8. The fourth mechanism

Copula-AR1 draws, for each subject, a vector of normal scores with
correlation $\rho^{|t_j - t_k|}$ (ρ = 0.6) and sets each outcome to one
when its score is below the normal quantile of the marginal mean of the
scenario. The marginal means are therefore exactly those of the
scenario, and the correlation between occasions of a subject decreases
with their distance in time, with no component that is constant within
subject. Without the convergence check the marginalized model was
clearly above 7.5% at 62 of the 480 design points. To check whether the
excess of the mixed model comes from a non-zero null value of its
conditional coefficient, a single dataset of 20,000 subjects was
generated under the null for T = 2, 4 and 10 and fitted with a random
intercept: the estimates ranged from 0.0046 to 0.060, with standard
errors up to 0.12, none distinguishable from zero. On the latent scale
the correlation of Copula-AR1 is that of the AR(1) working correlation
of Section 2.4 of the main text; on the binary outcomes it was 0.36
between adjacent occasions and 0.20 two occasions apart, so it decays
more slowly than an AR(1) correlation with the same first-lag
value.<!-- noclaim: checked in R/s7f/s7f_extras.R (lag-2 and lag-3 correlations above the AR(1) values) -->
The fourth mechanism was simulated on the main grid only, with the seeds
of the other three.

With a random intercept only, and N ≥ 85, under the null and averaged
over the design points of each scenario and number of occasions, the
model-based standard error of the mixed model (root mean square over
replicates) was 0.64 to 0.81 times the standard deviation of its
estimates, while the estimates were centred.<!-- noclaim: mechanism -->
The marginalized model estimated the interaction with a mean bias of
0.026 empirical standard errors and covered the true value in 94% of
replicates under a normal reference, although its model of the
dependence is misspecified on this mechanism (Table
S28).<!-- noclaim: descriptive -->

**Table S28.** Copula-AR1 (480 design points), Wald tests with a t
reference: mean type I error (%), design points clearly above 7.5% and
clearly below 2.5%, replicates returning a usable result under the null
(%, averaged over design points), power (%), and, for the marginal
estimators, mean bias in empirical standard errors and coverage of the
95% interval with a normal and a t reference (%); n.a. = not applicable.
Monte Carlo standard errors are at most 1.2 points for the type I error,
0.58 for the usable share, 1.1 for power, 0.041 for bias in empirical
standard errors and 0.91 points for coverage. These are upper bounds:
the mean over design points of the per-point standard error.

| Strategy | Type I error | Above 7.5% | Below 2.5% | Usable | Power | Bias / empirical SE | Coverage, normal | Coverage, t |
|----|----|----|----|----|----|----|----|----|
| GEE AR(1), Fay-Graubard | 4.6 | 0 | 17 | 91.4 | 17.4 | 0.0 | 94.7 | 95.5 |
| Marginalized model with check | 4.9 | 0 | 3 | 83.2 | 18.6 | 0.0 | 94.0 | 94.9 |
| Marginalized model without check | 6.3 | 62 | 2 | 89.6 | 19.4 | 0.0 | 92.3 | 93.3 |
| GEE independence, Fay-Graubard | 4.7 | 1 | 20 | 91.7 | 17.0 | 0.0 | 94.8 | 95.6 |
| GEE exchangeable, Fay-Graubard | 5.0 | 8 | 17 | 91.5 | 17.3 | 0.0 | 94.3 | 95.1 |
| GEE exchangeable, Mancl-DeRouen | 5.3 | 28 | 18 | 91.4 | 17.5 | 0.0 | 93.9 | 94.8 |
| GEE exchangeable, Kauermann-Carroll | 6.1 | 61 | 11 | 91.2 | 18.5 | 0.0 | 93.0 | 94.1 |
| GEE exchangeable, uncorrected | 7.6 | 100 | 0 | 93.5 | 20.1 | 0.0 | 91.2 | 92.5 |
| GLMM, Wald | 10.2 | 223 | 0 | 79.2 | 23.2 | n.a. | n.a. | n.a. |
| Shrinkage GLM | 11.4 | 244 | 0 | 81.3 | 24.1 | n.a. | n.a. | n.a. |
| GLMM, likelihood ratio | 9.4 | 189 | 10 | 86.4 | 24.7 | n.a. | n.a. | n.a. |
| Sign-flip score test | 5.4 | 0 | 0 | 96.1 | 18.2 | n.a. | n.a. | n.a. |

## S9. Second clinical application: the respiratory trial

The respiratory trial distributed with *geepack* has 111 patients
observed at four visits. The full mean model had treatment, baseline
status, visit, age, sex and centre; the reduced model had treatment and
baseline status only. Kauermann-Carroll was computed with *fewclust* on
the same GEE fits and the same draws, in both applications;
Mancl-DeRouen and Fay-Graubard were computed with *geesmv*. Recomputed
with *fewclust*, these two agreed with *geesmv* to within one per cent
in 99% and 99% of the draws of the cohort and in 99.9% and 99.6% of
those of the trial. In both applications the marginalized model had a
latent variable with an intercept only (20 quadrature points) and the
main text reports its model-based standard error. With its robust
standard error the corrected ratio of Section 2.7 of the main text was
0.85 to 0.99 in the cohort and 0.55 to 0.87 in the trial; on the full
data the robust standard error was 0.17 in the cohort and 0.35 in the
trial. In the respiratory trial the contrast of interest was the
treatment effect, adjusted for baseline status, visit, age, sex and
centre; the latent variable of the marginalized model had an intercept
only, as in the cohort, since with four visits per patient a time trend
of the latent standard deviation would rest on four observations per
patient.

In the trial, the bias of the marginalized model against its full-data
estimate was within two Monte Carlo standard errors of zero at 10 and at
50 patients and reached 0.40 at 15 (Table S30). Its corrected ratio,
with the model-based standard error, was 0.70 at 10 patients and 0.85 at
50.<!-- noclaim: descriptive -->

**Table S29.** Full-data estimate (SE) of the treatment effect, full and
reduced mean model; GEE with exchangeable working correlation;
model-based SE for the marginalized model and the mixed model. The mixed
model estimates a conditional effect.

| Method                 | Full model | Reduced model |
|------------------------|------------|---------------|
| Marginalized model     | 1.3 (0.34) | 1.2 (0.33)    |
| GEE, uncorrected       | 1.3 (0.35) | 1.2 (0.33)    |
| GEE, Mancl-DeRouen     | 1.3 (0.37) | 1.2 (0.34)    |
| GEE, Kauermann-Carroll | 1.3 (0.36) | 1.2 (0.33)    |
| GEE, Fay-Graubard      | 1.3 (0.35) | 1.2 (0.33)    |
| GLMM                   | 2.1 (0.54) | 2.1 (0.53)    |

**Table S30.** Subsamples of whole patients, full mean model: draws not
returned, with an estimate but no SE, extreme (\|estimate\| or SE above
10) and returning a usable result (%), bias of the estimate against the
full-data estimate, and ratio of the mean standard error (model-based
for the marginalized model and the mixed model, sandwich for the GEEs)
to the standard deviation of the estimates multiplied by the
finite-population factor of Section 2.7 of the main text. Largest Monte
Carlo standard error: 2.2 points for the usable share, 0.21 for the bias
and 0.15 for the SE ratio.

| Method | N | Not returned | No SE | Extreme | Usable (%) | Bias | SE ratio, corrected |
|----|----|----|----|----|----|----|----|
| Marginalized model | 10 | 30.6 | 1.4 | 18.2 | 50 | 0.2 | 0.7 |
|  | 15 | 8.2 | 2.6 | 5.6 | 84 | 0.4 | 0.7 |
|  | 20 | 1.8 | 2.0 | 1.8 | 94 | 0.2 | 0.7 |
|  | 25 | 0.0 | 0.8 | 0.0 | 99.2 | 0.2 | 0.7 |
|  | 50 | 0.0 | 0.0 | 0.0 | 100.0 | 0.0 | 0.9 |
| GEE, uncorrected | 10 | 10.2 | 14.6 | 28.4 | 47 | 0.2 | 0.5 |
|  | 15 | 1.8 | 3.6 | 11.6 | 83 | 0.5 | 0.6 |
|  | 20 | 0.8 | 1.4 | 3.8 | 94 | 0.2 | 0.7 |
|  | 25 | 0.0 | 0.0 | 0.8 | 99.2 | 0.2 | 0.7 |
|  | 50 | 0.0 | 0.0 | 0.0 | 100.0 | 0.0 | 0.9 |
| GEE, Mancl-DeRouen | 10 | 10.2 | 63.4 | 3.4 | 23 | 0.1 | 2.2 |
|  | 15 | 1.8 | 30.0 | 0.6 | 68 | 0.4 | 1.4 |
|  | 20 | 0.8 | 11.4 | 0.2 | 88 | 0.2 | 1.2 |
|  | 25 | 0.0 | 2.6 | 0.0 | 97 | 0.2 | 1.0 |
|  | 50 | 0.0 | 0.0 | 0.0 | 100.0 | 0.0 | 1.0 |
| GEE, Kauermann-Carroll | 10 | 10.2 | 65.4 | 2.8 | 22 | 0.4 | 1.0 |
|  | 15 | 1.8 | 27.6 | 1.4 | 69 | 0.6 | 0.9 |
|  | 20 | 0.8 | 9.2 | 0.6 | 89 | 0.2 | 0.9 |
|  | 25 | 0.0 | 2.2 | 0.6 | 97 | 0.2 | 0.9 |
|  | 50 | 0.0 | 0.0 | 0.0 | 100.0 | 0.0 | 0.9 |
| GEE, Fay-Graubard | 10 | 10.2 | 61.8 | 4.6 | 23 | 0.1 | 2.2 |
|  | 15 | 1.8 | 28.4 | 2.0 | 68 | 0.3 | 1.3 |
|  | 20 | 0.8 | 11.0 | 0.0 | 88 | 0.2 | 1.0 |
|  | 25 | 0.0 | 2.4 | 0.0 | 98 | 0.2 | 0.9 |
|  | 50 | 0.0 | 0.0 | 0.0 | 100.0 | 0.0 | 0.9 |
| GLMM | 10 | 9.6 | 0.0 | 37.0 | 53 | -0.6 | 0.7 |
|  | 15 | 1.8 | 0.0 | 10.0 | 88 | -0.1 | 0.7 |
|  | 20 | 0.8 | 0.0 | 2.4 | 97 | -0.1 | 0.6 |
|  | 25 | 0.0 | 0.0 | 0.0 | 100.0 | 0.0 | 0.7 |
|  | 50 | 0.0 | 0.0 | 0.0 | 100.0 | -0.1 | 0.7 |

## S10. The marginalized model and the GEE strategy compared

The strategy selected by the decision rule, the marginalized model with
the convergence check, was compared with the GEE strategy (AR(1) working
correlation, Fay-Graubard correction) by number of clusters and of
occasions after the results had been seen, both with a t reference, with
the block bootstrap and the one-point margin for mean differences of
Section 2.5 of the main text. The marginalized model was never clearly
above 7.5%; the GEE strategy was clearly above 7.5% at 2 design points,
both with ten clusters and two occasions (Section 3.4 of the main text).
Before the datasets with a constant outcome were excluded (Section 2.5
of the main text), the GEE strategy was usable in 172 and 136 replicates
at these two design points, with the same 12 and 11 rejections: every
excluded dataset had been a usable non-rejection, so the exclusion
raised the rate. The comparison does not enter the rule (Table 5). With
two or three occasions and at most 25 clusters the GEE strategy returned
a usable result more often, by 16 points.<!-- claim: reg:hb as=descr -->

**Table S31.** GEE with AR(1) working correlation and Fay-Graubard (GEE)
and the marginalized multilevel model with the convergence check (MMM),
both with a t reference (N - p degrees of freedom), over the 2,376
design points of the four mechanisms, by number of clusters N and of
occasions T: mean type I error (%), mean absolute departure from 5%
(percentage points), their difference GEE minus MMM with its Monte Carlo
standard error (block bootstrap), its reading against the one-point
margin (the strategy closer to 5% when the interval of ± 2 standard
errors lies beyond the margin, similar when it lies within it, about the
margin when it straddles it), datasets returning a usable result under
the null (%) and datasets in which the simulated effect was detected
among those attempted (%), both averaged over design points.
Block-bootstrap Monte Carlo standard errors are at most 0.47 points for
the type I error, 0.42 for the departure from 5%, 0.30 for the usable
share and 0.069 for power.

| Clusters | Occasions | Type I error, GEE / MMM | Departure from 5%, GEE / MMM | Difference (MCSE) | Reading | Usable, GEE / MMM | Power, GEE / MMM |
|----|----|----|----|----|----|----|----|
| N = 10 | T = 2 | 3.7 / 0.5 | 5.4 / 5.1 | 0.3 (0.54) | about the margin | 10.2 / 4.3 | 0.2 / 0.0 |
| N = 10 | T = 3 | 4.6 / 1.2 | 1.9 / 3.9 | -1.9 (0.25) | GEE | 35.0 / 21.0 | 1.3 / 0.4 |
| N = 10 | T ≥ 4 | 4.1 / 1.8 | 1.6 / 3.2 | -1.6 (0.088) | GEE | 76.7 / 59.4 | 3.6 / 2.5 |
| N = 15-25 | T = 2 | 0.8 / 1.9 | 4.3 / 3.3 | 1.1 (0.13) | about the margin | 32.3 / 19.3 | 0.2 / 0.6 |
| N = 15-25 | T = 3 | 4.5 / 2.7 | 1.4 / 2.5 | -1.1 (0.12) | about the margin | 69.2 / 45.9 | 3.1 / 1.7 |
| N = 15-25 | T ≥ 4 | 4.7 / 3.9 | 1.0 / 1.4 | -0.43 (0.047) | similar | 94.3 / 83.7 | 9.3 / 9.9 |
| N \> 25 | T = 2 | 2.9 / 4.1 | 2.1 / 1.3 | 0.78 (0.045) | similar | 91.3 / 67.7 | 3.2 / 3.3 |
| N \> 25 | T = 3 | 4.7 / 5.0 | 0.6 / 0.7 | -0.12 (0.036) | similar | 98.9 / 91.2 | 6.6 / 6.7 |
| N \> 25 | T ≥ 4 | 4.8 / 5.2 | 0.6 / 0.6 | -0.08 (0.017) | similar | 99.9 / 99.0 | 34.0 / 35.9 |
| all | all | 4.4 / 4.4 | 1.1 / 1.2 | -0.10 (0.016) | similar | 91.9 / 84.0 | 19.4 / 20.3 |

## S11. What the decision rule guarantees

The rule replaces a criterion that selected at each cluster count the
method with the highest unconditional validity, and was fixed after the
results of the version that used that criterion had been seen (Section
4.1 of the main
text).<!-- noclaim: description of the earlier criterion -->

A design point is clearly above 7.5% when the lower limit of the exact
95% interval of its rejection rate exceeds 7.5%. With 1,000 usable
replicates this requires an observed rate of at least 9.3%, with half as
many 10%, and with a tenth as many 14%. A strategy that is never clearly
above 7.5% is therefore one in which no excess was detected. Table S32
gives, for every strategy of Table 5, the design points at which an
excess could not be excluded (upper limit of the exact interval above
7.5%, but not clearly above) and the highest rate observed. An exact
interval was used because with no rejections the binomial Monte Carlo
standard error is zero. The highest rate at a design point was 20% for
the GEE strategy and 8.5% for the marginalized model, and an excess
could not be excluded, the upper limit of the exact interval lying above
7.5%, at 227 and 402 design points (Table
S32).<!-- noclaim: descriptive values -->

**Table S32.** Strategies of Table 5 over the 2,376 design points of the
four mechanisms, under the null: design points clearly above 7.5%,
design points at which an excess could not be excluded, design points
not assessed (fewer than 20 usable replicates), and the highest
rejection rate at an assessed design point (%) with its number of usable
replicates. Wald tests with a t reference (N - p degrees of freedom),
likelihood-ratio test with a $\chi^2$ reference. The binomial Monte
Carlo standard error of the highest rate is at most 6.0 points.

| Strategy | Clearly above 7.5% | Excess not excluded | Not assessed | Highest rate (%) | Usable replicates |
|----|----|----|----|----|----|
| GEE AR(1), Fay-Graubard | 2 | 227 | 0 | 19.7 | 61 |
| Marginalized model with check | 0 | 402 | 15 | 8.5 | 991 |
| GEE independence, Fay-Graubard | 5 | 271 | 0 | 19.3 | 57 |
| GEE exchangeable, Fay-Graubard | 15 | 404 | 0 | 19.7 | 61 |
| GEE AR(1), Mancl-DeRouen | 20 | 369 | 0 | 19.0 | 58 |
| GEE independence, Mancl-DeRouen | 58 | 402 | 0 | 19.6 | 56 |
| GEE exchangeable, Mancl-DeRouen | 63 | 453 | 0 | 19.0 | 58 |
| GEE AR(1), Kauermann-Carroll | 138 | 589 | 0 | 22.4 | 49 |
| GEE independence, Kauermann-Carroll | 189 | 595 | 0 | 22.6 | 53 |
| GEE exchangeable, Kauermann-Carroll | 207 | 651 | 0 | 22.4 | 49 |
| GEE independence, CR2 | 251 | 561 | 0 | 22.2 | 677 |
| GEE exchangeable, CR2 | 254 | 660 | 0 | 22.3 | 677 |
| GEE AR(1), uncorrected | 363 | 737 | 0 | 24.7 | 648 |
| GEE independence, uncorrected | 383 | 770 | 0 | 25.7 | 835 |
| GEE exchangeable, uncorrected | 411 | 799 | 0 | 26.7 | 681 |
| Shrinkage GLM | 608 | 687 | 320 | 25.3 | 949 |
| GLMM, Wald | 515 | 702 | 322 | 19.8 | 979 |
| GLMM, likelihood ratio | 412 | 450 | 304 | 28.4 | 1,000 |
| Sign-flip score test | 0 | 760 | 0 | 8.8 | 296 |

The zero count of the marginalized model, and the two design points of
the GEE strategy, are uncertain in a replication of the study. In the
block bootstrap, GEE with AR(1) working correlation and Fay-Graubard had
at least one design point clearly above 7.5% in 99% of resamples, and
the marginalized model with the check in 58% (Table S37 in Section S18).
Among the eight strategies of Table S37 the marginalized model has the
fewest expected design points clearly above 7.5% in all three
calculations: in the bootstrap (0.86, against 3.7 for the GEE strategy),
in the plug-in calculation with the observed rates (0.84 against 3.7)
and with the rates shrunk by empirical Bayes (normal-normal method of
moments on the rate scale, within strategy and mechanism; 0.12 against
0.15). Shrinkage reduces the GEE strategy’s count most, because its two
design points rest on 61 and 63 usable
replicates.<!-- claim: reg:s37 as=descr -->

The usability limits of Section 2.5 of the main text (\|estimate\| and
standard error below 10) were varied to 5 and 20 for the strategies of
Table 5 with an estimate and an interval, on the stored per-replicate
results. With the published limit the recomputation reproduced Table 5
exactly. With usability limits of 5 or 20 in place of 10, the
marginalized model was never clearly above 7.5%. The GEE strategy was
never clearly above 7.5% with the limit of 5, so that the rule of
Section 2.6 of the main text would have selected it (usable share 92%
against 84%), and was clearly above 7.5% at 6 design points with the
limit of 20 (Table S33). Its two design points with the published limit
therefore depend on datasets with an estimate or a standard error
between 5 and 10 in absolute
value.<!-- claim: reg:filter_t4 as=descr -->

**Table S33.** Design points clearly above 7.5% under the null and
usable share under the null averaged over design points (%) with
usability limits of 5, 10 (published) and 20, 2,376 design points, Wald
tests with a t reference. The Monte Carlo standard error of the mean
usable share is at most 0.56 points (upper bound: mean over design
points of the binomial standard error).

| Strategy | Above 7.5%, limit 5 | Above 7.5%, limit 10 | Above 7.5%, limit 20 | Usable, limit 5 | Usable, limit 10 | Usable, limit 20 |
|----|----|----|----|----|----|----|
| GEE AR(1), Fay-Graubard | 0 | 2 | 6 | 91.9 | 91.9 | 91.9 |
| Marginalized model with check | 0 | 0 | 0 | 84.0 | 84.0 | 84.0 |
| GEE independence, Fay-Graubard | 2 | 5 | 8 | 92.1 | 92.1 | 92.1 |
| GEE exchangeable, Fay-Graubard | 13 | 15 | 23 | 92.0 | 92.0 | 92.0 |
| GEE AR(1), Mancl-DeRouen | 19 | 20 | 29 | 91.8 | 91.8 | 91.9 |
| GEE independence, Mancl-DeRouen | 55 | 58 | 60 | 92.0 | 92.0 | 92.0 |
| GEE exchangeable, Mancl-DeRouen | 61 | 63 | 78 | 91.9 | 91.9 | 92.0 |
| GEE AR(1), Kauermann-Carroll | 133 | 138 | 161 | 91.7 | 91.7 | 91.8 |
| GEE independence, Kauermann-Carroll | 186 | 189 | 193 | 92.0 | 92.0 | 92.0 |
| GEE exchangeable, Kauermann-Carroll | 200 | 207 | 255 | 91.8 | 91.8 | 91.9 |
| GEE independence, CR2 | 241 | 251 | 252 | 93.8 | 93.9 | 93.9 |
| GEE exchangeable, CR2 | 241 | 254 | 284 | 93.4 | 93.5 | 93.6 |
| GEE AR(1), uncorrected | 331 | 363 | 400 | 93.6 | 93.7 | 93.8 |
| GEE independence, uncorrected | 376 | 383 | 394 | 93.9 | 94.0 | 94.0 |
| GEE exchangeable, uncorrected | 378 | 411 | 451 | 93.7 | 93.8 | 94.0 |
| Shrinkage GLM | 602 | 608 | 622 | 81.5 | 81.7 | 81.9 |
| GLMM, Wald | 434 | 515 | 553 | 79.0 | 80.4 | 81.1 |

## S12. Adaptive quadrature for the mixed model

This analysis was added after the main run. The mixed model with a
random intercept (scenarios M1, M2 and M1c) was refitted with
*lme4::glmer* and *nAGQ* = 10 on the same null datasets, at 123 design
points: N ∈ {10, 25, 55, 145} and T ∈ {2, 4, 10} in M1 and M2 on the
four mechanisms, and N ∈ {10, 25, 45} and T ∈ {2, 4, 6} in M1c on the
three study mechanisms. Only the number of quadrature nodes changed;
formula, optimizer and usability limits were those of the study.
Refitted by the Laplace approximation, the 122,904 replicate fits
reproduced the stored estimates and standard errors exactly. On the
first 100 replicates of each design point the fit was repeated with 25
nodes, which gave a different decision in 10 of the 9,844 fits usable
for both.
<!-- noclaim: description of a sensitivity analysis added after the main run -->
Quadrature lowered the rejection rate most at two occasions with 145
clusters, from 10% to 20% under the Laplace approximation to 2.6% to
4.5%.<!-- noclaim: descriptive values of the sensitivity analysis -->

**Table S34.** Wald test of the mixed model with a random intercept and
a t reference, under the null, Laplace approximation (L) against
adaptive quadrature (Q), by mechanism and number of occasions T: design
points, mean type I error (%), design points clearly above 7.5% and
clearly below 2.5% (exact 95% interval), and mean usable share (%). Last
column: mean difference Q minus L over the design points of the
mechanism (percentage points), with its 95% block-bootstrap interval.
Monte Carlo standard errors are at most 3.1 points for the mean type I
error and 1.1 for the usable share. These are upper bounds: the mean
over design points of the binomial standard error.

| Mechanism | Occasions | Design points | Type I error, L | Type I error, Q | Above 7.5%, L | Above 7.5%, Q | Below 2.5%, L | Below 2.5%, Q | Usable, L | Usable, Q | Q minus L (mechanism) |
|----|----|----|----|----|----|----|----|----|----|----|----|
| GenBinaryY | T = 2 | 11 | 8.0 | 1.9 | 3 | 0 | 0 | 4 | 53.4 | 56.1 | -3.4 (-3.7 to -3.2) |
|  | T = 4 | 11 | 6.0 | 3.4 | 1 | 0 | 0 | 4 | 86.0 | 86.6 |  |
|  | T = 6 | 3 | 9.8 | 4.8 | 2 | 0 | 0 | 1 | 97.1 | 97.4 |  |
|  | T = 10 | 8 | 7.6 | 6.4 | 3 | 2 | 0 | 2 | 98.1 | 98.1 |  |
| GLMM-native | T = 2 | 11 | 6.3 | 1.8 | 2 | 0 | 0 | 3 | 57.7 | 58.6 | -3.3 (-3.6 to -3.0) |
|  | T = 4 | 11 | 3.8 | 2.3 | 0 | 0 | 1 | 4 | 89.0 | 89.3 |  |
|  | T = 6 | 3 | 8.0 | 2.8 | 2 | 0 | 0 | 1 | 97.8 | 97.8 |  |
|  | T = 10 | 8 | 7.1 | 3.0 | 2 | 0 | 0 | 2 | 98.8 | 98.8 |  |
| GLMM-AR1 | T = 2 | 11 | 5.7 | 1.9 | 2 | 0 | 0 | 6 | 56.2 | 57.6 | -3.1 (-3.5 to -2.7) |
|  | T = 4 | 11 | 3.8 | 2.5 | 0 | 0 | 2 | 4 | 88.6 | 88.9 |  |
|  | T = 6 | 3 | 8.4 | 3.3 | 2 | 0 | 0 | 1 | 97.8 | 97.9 |  |
|  | T = 10 | 8 | 7.1 | 4.1 | 2 | 0 | 0 | 2 | 98.7 | 98.7 |  |
| Copula-AR1 | T = 2 | 8 | 12.1 | 3.1 | 3 | 0 | 0 | 1 | 48.0 | 51.2 | -4.0 (-4.6 to -3.5) |
|  | T = 4 | 8 | 7.0 | 4.6 | 3 | 0 | 0 | 2 | 82.1 | 82.4 |  |
|  | T = 10 | 8 | 12.8 | 10.8 | 6 | 6 | 0 | 0 | 98.3 | 98.3 |  |

## S13. Estimation by the marginalized model with the check in scenarios M6-M7

In M6-M7 the fit that the convergence check substitutes, without a time
trend of the latent standard deviation, is misspecified. Table S35 gives
bias, standard error and coverage of the marginalized model with the
check against the marginal true value, and the same coverage and
standard-error ratio separately for usable fits that are the fit without
the time trend and for those that are the original fit; a fit at the
boundary whose refit did not converge gives no usable result (Section
2.4 of the main text). Intervention depends on the data, so the split is
a description, not a criterion. The estimates of the marginalized model
with the check had a bias of at most 0.55 empirical standard errors at
any design point, its model-based standard error, averaged over
replicates, was a median 2.0% to 4.2% below the empirical one across
scenarios and hypotheses, and the t interval covered the true value in
95% to 96% of replicates on average, with non-coverage clearly above
7.5% at no design point (Table S35). Usable fits that are the fit
without the time trend covered the true value as often as the others
(Table S35). Under a normal reference its mean type I error was 5.2% to
5.8% across the four scenarios. GEE with exchangeable working
correlation and Fay-Graubard with a t reference was clearly above 7.5%
at 4 design points, 3 of them on GenBinaryY with the steeper slope (M7)
(Table 5 and Figure S1).<!-- claim: reg:t4 as=descr -->

<img src="figures/figS1.png" style="width:100.0%" />

**Figure S1.** Type I error at every design point of scenarios M6-M7 (96
design points each) for the strategies of Figure 1, with the corrections
from *fewclust* as in Table 5: Wald tests with a t reference (N - p
degrees of freedom) and the likelihood-ratio test with a $\chi^2$
reference. Rows: GenBinaryY with latent-SD slope 0.1 and 0.25;
GLMM-native with random-slope variance 0.01 and 0.04. Design points are
ordered by number of occasions T (blocks) and number of clusters N. Grey
band: 2.5% to 7.5%; line: nominal 5%. Design points with fewer than 20
usable replicates are not drawn (60 design points of the
likelihood-ratio test, none of the other strategies). Circles: design
points clearly above 7.5% (lower limit of the exact 95% interval above
7.5%), as in Table 5: 18 for GEE exchangeable with Mancl-DeRouen, 4 for
GEE exchangeable with Fay-Graubard, all on GenBinaryY; 0 for GEE AR(1)
with Fay-Graubard, 0 for the marginalized model with the check and 0 for
the likelihood-ratio test. The binomial Monte Carlo standard error of a
drawn rate is at most 3.4 points.

**Table S35.** Marginalized model with the convergence check in
scenarios M6-M7, by hypothesis, over design points with at least 20
usable replicates: mean and largest absolute bias in empirical standard
errors, median of model-based over empirical standard error minus one
(%; model-based standard error averaged over replicates), mean and
smallest coverage of the 95% t interval (%), design points with
non-coverage clearly above 7.5%, mean share of usable fits that are the
fit without the time trend (%), and mean coverage and median
standard-error ratio minus one for usable fits that are the fit without
the time trend (with) and that are the original fit (without). Monte
Carlo standard errors, with the empirical standard error treated as
known, are at most 0.038 for the mean bias in standard errors, 0.22 for
the largest, 0.73 points for the mean coverage, 0.86 for the smallest
coverage and 1.8 for the share of usable fits that are the fit without
the time trend. For the SE ratio the typical per-point value is 2.2
points. The columns by intervention of the check have no Monte Carlo
standard error.

| Scenario | Hypothesis | Design points | Bias / SE, mean | Bias / SE, largest | SE ratio minus one (%) | Coverage (%) | Smallest coverage (%) | Non-coverage above 7.5% | Fit without trend (%) | Coverage, with / without | SE ratio minus one, with / without |
|----|----|----|----|----|----|----|----|----|----|----|----|
| GenBinaryY, slope 0.1 | H0 | 95 | 0.000053 | 0.14 | -2.5 | 95.9 | 93.0 | 0 | 44.9 | 96.1 / 95.5 | -3.1 / -1.8 |
|  | H1 | 96 | 0.036 | 0.55 | -2.0 | 95.7 | 93.1 | 0 | 43.1 | 96.0 / 95.3 | -3.2 / -1.6 |
| GenBinaryY, slope 0.25 | H0 | 95 | 0.0018 | 0.13 | -3.8 | 95.5 | 92.4 | 0 | 35.6 | 95.9 / 95.0 | -4.4 / -3.6 |
|  | H1 | 95 | 0.043 | 0.17 | -4.2 | 95.2 | 92.2 | 0 | 33.7 | 95.6 / 94.8 | -4.2 / -3.1 |
| GLMM-native, slope variance 0.01 | H0 | 95 | 0.0017 | 0.14 | -3.6 | 95.8 | 93.0 | 0 | 53.1 | 95.9 / 95.5 | -4.7 / -2.8 |
|  | H1 | 95 | 0.026 | 0.18 | -3.0 | 95.8 | 93.4 | 0 | 51.5 | 95.9 / 95.6 | -4.3 / -2.7 |
| GLMM-native, slope variance 0.04 | H0 | 95 | -0.00091 | 0.081 | -3.9 | 95.5 | 92.8 | 0 | 45.1 | 95.7 / 95.2 | -4.7 / -3.5 |
|  | H1 | 96 | 0.012 | 0.45 | -3.5 | 95.5 | 92.6 | 0 | 43.5 | 95.8 / 95.0 | -3.4 / -3.6 |

## S14. Previous evaluations of the small-sample corrections

Cluster-randomized trials with ten to thirty clusters are common. Among
100 reviewed trials, the 99 that reported the number of clusters had a
median of 25 (interquartile range 15 to 44) \[5\], and Barker et al.,
citing an earlier review, reported that 45% of stepped-wedge trials had
fewer than ten clusters \[6\].<!-- noclaim: literature --> In the
evaluations cited here, the unit within a cluster is a site within a
mouth \[7\], a patient within a general practice \[9\], a participant
within a cluster \[15,19\] or a youth within a community \[20\]. None is
a subject measured on several occasions. Cluster sizes range from 4
observations per cluster \[20\] to 150 on average \[15\].

Three results are established elsewhere and are cited rather than
claimed here: the uncorrected sandwich is anticonservative with few
clusters \[7,8\], the corrections trade level for power \[15,17\], and
marginal and conditional estimands differ
\[4\].<!-- noclaim: cited results -->

## S15. Further details of the methods and performance measures

For the likelihood-ratio test, both nested fits had to return and the
p-value had to be finite. Its full fit used the optimizer bobyqa in both
stages with a limit of $10^5$ function evaluations, whereas the Wald
test used the default of *lme4* 2.0.6 (bobyqa, then Nelder-Mead); the
Wald test of that full fit is reported as the mixed model with bobyqa
(Section S3 and Section S5).

The Wald statistic of the shrinkage GLM equals that of the mixed model
on the same fit; its rejection rate can still differ, because the
usability limits of Section 2.5 of the main text apply to the rescaled
estimate and standard error.

In the block bootstrap of the direct comparison of the two strategies of
the rule (Table S31), scenario M1c, whose seeds are those of M1, was
resampled with its own indices rather than jointly with M1; this affects
only the Monte Carlo standard errors of that
comparison.<!-- noclaim: definition -->

In the text, values are given to two significant figures, or fewer where
the MCSE does not support them.<!-- noclaim: definition --> Tables give
the digit that the MCSE supports, with the MCSE in the caption.

## S16. The convergence check of the marginalized model: diagnosis and sensitivity analyses

*MMLB::mm* writes the latent standard deviation with an exponentiated
slope, so a slope near zero lies on the boundary of the parameter space,
and in scenarios M3-M5 the generating slope was close to zero. At the
boundary *nlm* either stopped without converging or converged with a
collapsed standard error. Both failures produced excess rejections that
grew with the number of occasions. This is why, without the check, the
marginalized model appears miscalibrated with long series.

The check was specified after the diagnosis and before the reanalysis of
the 1,685,204 fits of M3-M5 that returned an estimate. The check changed
raw power at eight or ten occasions by 0.26 points on average.

Thresholds of 0.0001 and 0.01 were examined as sensitivity analyses, and
scenarios M6-M7 were simulated afterwards to test the check. With the
thresholds 0.0001 and 0.01 in place of 0.001, and the same rules for the
refit and for constant outcomes, the tables were recomputed from the
stored fits. The reading (difference, similar or about the margin;
Section 2.5 of the main text) changed in 0 and 0 of the 156 comparisons
of the claim table of the main grid; of the 37 cited in this paper, 0
and 0 changed. With either threshold, the marginalized model with a t
reference was clearly above 7.5% at 0 and 0 design points of the main
grid and
M1c.<!-- noclaim: descriptive values of the sensitivity analysis -->

## S17. The neonatal cohort: subsampling details

Subsampling whole patients preserves the dependence structure exactly,
and the covariate distribution only in expectation. Of the 25 patients
ventilated invasively, a draw of 10 contains none with probability 0.27,
and the mean model is then not identifiable.

In subsamples of 10 patients, 22% of draws contained no invasively
ventilated patient (expected proportion 0.27 under the simple random
draws of Section 2.7 of the main text, a deviation of 2.5 binomial
standard errors at this size and of at most 1.1 at the others), so that
the mean model could not be fitted (Table S15). Among the fits that were
returned, a usable result was obtained by the mixed model in 100%, by
the marginalized model in 100% and by the uncorrected GEE in 92%. It was
obtained by Fay-Graubard in 45% and by Mancl-DeRouen in 20%; most
failures of the corrections returned a point estimate without a standard
error.<!-- noclaim: descriptive values per method -->

Section 3.11 of the main text reports, for the mixed model in the
cohort, the ratio of its model-based standard error to the empirical
standard deviation of its estimates. A random intercept for the
patient-timepoint did not change this ratio (0.50 to 0.57), because its
variance was estimated at the boundary in the subsamples. The
posterior-versus-anterior contrast varied strongly between patients,
which a random intercept ignores: with a random slope for posterior
zones (variance 10 in the full cohort) the ratio rose to 0.66 to 0.73 at
25 patients or more, so this heterogeneity explains part of the
shortfall.<!-- noclaim: explanation examined after the results -->

The bias of the marginalized model against its full-cohort estimate fell
from 0.17 at 10 patients to a value within two Monte Carlo standard
errors of zero at 50.

## S18. Standard errors and coverage, and stability of the decision rule

Table S36 gives, under the alternative, the bias in standard errors, the
ratio of the model-based to the empirical standard error and the
coverage of each method, cited in Section 2.5 of the main text. Table
S37 gives the stability of the count of design points clearly above 7.5%
in a replication of the study, discussed in Section
S11.<!-- noclaim: pointer to the tables -->

**Table S36.** Bias in empirical standard errors (Bias/SE), model-based
over empirical standard error minus one (SE ratio, %; model-based
standard error as the root mean square over replicates \[21\]) and
coverage of the 95% confidence interval against each method’s own
estimand (%), under the alternative, on each method’s usable replicates,
main grid (1,440 design points) and M1c (72 design points), marginalized
model with the convergence check. Each value is the mean over the design
points with at least 20 usable replicates: at least 416 of the 480
design points of a mechanism in the main grid, all 24 in M1c. Not
available (n.a.) for the GLMM on GenBinaryY, which has no marginal
estimand. Monte Carlo standard errors are at most 0.039 for Bias/SE, 5.9
points for the SE ratio and 1.0 points for coverage. These are upper
bounds: the mean over design points of the per-point Monte Carlo
standard errors.

| Method | Grid | Bias/SE GB | Bias/SE NAT | Bias/SE AR1 | SE ratio GB | SE ratio NAT | SE ratio AR1 | Coverage GB | Coverage NAT | Coverage AR1 |
|----|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | Main | 0.03 | 0.02 | 0.02 | 7 | 8 | 9 | 95 | 95 | 95 |
|  | M1c | 0.06 | 0.04 | 0.05 | 7 | 9 | 8 | 94 | 95 | 95 |
| GEE, Mancl-DeRouen | Main | 0.03 | 0.02 | 0.02 | 6 | 8 | 8 | 94 | 95 | 95 |
|  | M1c | 0.06 | 0.04 | 0.05 | 14 | 16 | 15 | 95 | 95 | 95 |
| Marginalized model | Main | 0.03 | 0.02 | 0.02 | 2 | 2 | 2 | 95 | 95 | 95 |
|  | M1c | 0.06 | 0.03 | 0.04 | 7 | 2 | 1 | 93 | 93 | 93 |
| GLMM, Wald | Main | n.a. | 0.04 | -0.02 | 4 | 6 | 6 | n.a. | 94 | 93 |
|  | M1c | n.a. | 0.02 | 0.00 | 8 | 8 | 5 | n.a. | 95 | 95 |
| Shrinkage GLM | Main | 0.12 | 0.08 | 0.08 | 2 | 4 | 4 | 90 | 93 | 92 |
|  | M1c | 0.04 | 0.03 | 0.04 | -5 | 3 | 1 | 89 | 93 | 92 |
| GEE, uncorrected | Main | 0.03 | 0.02 | 0.02 | -8 | -7 | -7 | 91 | 92 | 92 |
|  | M1c | 0.06 | 0.04 | 0.05 | -9 | -8 | -9 | 90 | 91 | 91 |
| GEE, Kauermann-Carroll | Main | 0.03 | 0.02 | 0.02 | -1 | 0 | 0 | 93 | 94 | 94 |
|  | M1c | 0.06 | 0.04 | 0.05 | 1 | 2 | 2 | 93 | 94 | 94 |

**Table S37.** Stability of the count of design points clearly above
7.5% over the 2,376 design points of the four mechanisms, for the eight
strategies of Table 5 with a small-sample correction or the convergence
check. Observed: the count in the study (Table 5). Bootstrap: 2,000
block-bootstrap resamples of the replicates, the same resample for every
design point of a model; share of resamples with no design point clearly
above 7.5% (%), mean and 90th percentile of the count. Plug-in: design
points treated as independent, with true rate equal to the observed rate
(raw) or to the rate shrunk by empirical Bayes within strategy and
mechanism (shrunk); probability of no design point clearly above 7.5%
(%) and expected count.

| Strategy | Observed | Bootstrap, none (%) | Bootstrap, mean | Bootstrap, 90th percentile | Raw plug-in, none (%) | Raw plug-in, expected | Shrunk plug-in, none (%) | Shrunk plug-in, expected |
|----|----|----|----|----|----|----|----|----|
| GEE AR(1), Fay-Graubard | 2 | 0.5 | 3.7 | 6 | 0.6 | 3.7 | 86 | 0.15 |
| Marginalized model with check | 0 | 42 | 0.86 | 2 | 42 | 0.84 | 88 | 0.12 |
| GEE independence, Fay-Graubard | 5 | \< 0.1 | 5.7 | 8 | \< 0.1 | 5.6 | 41 | 0.82 |
| GEE exchangeable, Fay-Graubard | 15 | \< 0.1 | 17 | 21 | \< 0.1 | 17 | 0.3 | 5.0 |
| GEE AR(1), Mancl-DeRouen | 20 | \< 0.1 | 28 | 33 | \< 0.1 | 28 | \< 0.1 | 11 |
| GEE exchangeable, Mancl-DeRouen | 63 | \< 0.1 | 73 | 80 | \< 0.1 | 74 | \< 0.1 | 49 |
| GEE AR(1), Kauermann-Carroll | 138 | \< 0.1 | 142 | 149 | \< 0.1 | 142 | \< 0.1 | 110 |
| GEE exchangeable, Kauermann-Carroll | 207 | \< 0.1 | 208 | 217 | \< 0.1 | 209 | \< 0.1 | 180 |

## S19. Separated datasets

In a separated dataset a linear combination of the covariates of the
mean model predicts the outcome perfectly, so the maximum likelihood
estimate of the logistic mean model does not exist \[31\]. Every dataset
with an estimable design and a non-constant outcome was regenerated from
its seed and checked by linear programming \[32\] with the mean model of
the GEEs (time, treatment and their interaction, plus age in M2 and dose
in M5); the R package *detectseparation* \[33\] identified which
coefficients had infinite estimates. On a stratified sample of 9,504
datasets the linear-programming check written for the study and
*detectseparation* agreed on every dataset.<!-- noclaim: methods -->

Of the 4,739,647 datasets with an estimable design and a non-constant
outcome, 357,701 were separated (7.5%): 55% of those with two or three
occasions and at most 25 clusters, and 0.98% of those with 45 clusters
or more. In 357,545 the estimate of the interaction itself was infinite.
The usability limit of Section 2.5 of the main text does not identify
these datasets: in the sample, an ordinary logistic regression gave an
interaction above 10 in absolute value in 503 of the 729 separated
datasets, so the limit removed only part of
them.<!-- noclaim: descriptive values -->

The primary analysis keeps separated datasets and judges each method’s
result on them by the usability rule (Section 2.5 of the main text).
Over all design points, 0.19% of the usable results of the GEE strategy
came from datasets with an infinite interaction, against 4 of the
4,014,389 usable results of the marginalized model with the check. At
the two design points where the GEE strategy was clearly above 7.5% (M2,
N = 10, T = 2; GLMM-AR1 and GLMM-native), it rejected a true null in 12
of 61 and 11 of 63 usable replicates; 12 and 11 of these rejections were
on datasets with an infinite interaction, which were 27 and 26 of its
usable replicates. Refitting these separated datasets with 200 instead
of the 25 iterations of the study moved the GEE estimate further from
zero, by more than 0.5 on the logit scale, in 35 of
53.<!-- noclaim: descriptive values -->

Table S38 repeats the rule of Table 5 with separated datasets removed
from the usable results, in four sensitivity versions with unchanged
denominators: (a) datasets with an infinite interaction removed for the
strategies with a Wald interval (GEEs, marginalized model, Wald test of
the mixed model, shrinkage GLM); (b) the same datasets removed for every
strategy, including the likelihood-ratio and sign-flip tests; (c) every
separated dataset removed for the strategies with a Wald interval; (d)
as (a), plus the datasets in which the outcome varies within at most one
subject, all events or all non-events fall in one subject, one arm has
no events or only events, or one arm has a single subject observed after
the first occasion. In every version the GEE strategy is no longer
clearly above 7.5% at any design point and, answering more often, is the
strategy the rule selects; the marginalized model with the check is
selected only in the primary
analysis.<!-- noclaim: descriptive values of the sensitivity analysis -->

In the clinical applications, an infinite treatment contrast occurred in
365 of the 5,000 subsamples of the respiratory trial and in 0 subsamples
of the cohort.<!-- noclaim: descriptive values -->

**Table S38.** The decision rule of Table 5 with separated datasets
removed from the usable results, 2,376 design points of the four
mechanisms, Wald tests with a t reference and the likelihood-ratio test
with a $\chi^2$ reference; strategies of the rule and the two tests
without an interval. Above and below: design points clearly above 7.5%
and clearly below 2.5% (exact 95% interval, at least 20 usable null
replicates); usable: replicates returning a usable result under the null
(%), averaged over design points, every replicate with an estimable
design in the denominator; power: rejections of a false null per
attempted dataset (%). Versions as in Section S19.

| Version | Strategy | Above 7.5% | Below 2.5% | Usable (%) | Usable, T ≤ 3 and N ≤ 25 (%) | Power (%) |
|----|----|----|----|----|----|----|
| Primary (no removal) | Marginalized model with check | 0 | 42 | 84 | 28 | 20 |
|  | GEE AR(1), Fay-Graubard | 2 | 106 | 92 | 44 | 19 |
|  | GLMM, likelihood ratio (no interval) | 412 | 13 | 86 | 62 | 24 |
|  | Sign-flip score test (no interval) | 0 | 0 | 97 | 72 | 21 |
| \(a\) Infinite interaction removed, Wald strategies | Marginalized model with check | 0 | 42 | 84 | 28 | 20 |
|  | GEE AR(1), Fay-Graubard | 0 | 105 | 92 | 42 | 19 |
|  | GLMM, likelihood ratio (no interval) | 412 | 13 | 86 | 62 | 24 |
|  | Sign-flip score test (no interval) | 0 | 0 | 97 | 72 | 21 |
| \(b\) Infinite interaction removed, all strategies | Marginalized model with check | 0 | 42 | 84 | 28 | 20 |
|  | GEE AR(1), Fay-Graubard | 0 | 105 | 92 | 42 | 19 |
|  | GLMM, likelihood ratio (no interval) | 337 | 14 | 82 | 33 | 24 |
|  | Sign-flip score test (no interval) | 0 | 32 | 92 | 43 | 20 |
| \(c\) Any separation removed, Wald strategies | Marginalized model with check | 0 | 42 | 84 | 28 | 20 |
|  | GEE AR(1), Fay-Graubard | 0 | 105 | 92 | 42 | 19 |
|  | GLMM, likelihood ratio (no interval) | 412 | 13 | 86 | 62 | 24 |
|  | Sign-flip score test (no interval) | 0 | 0 | 97 | 72 | 21 |
| \(d\) As (a), plus the other degenerate classes | Marginalized model with check | 0 | 38 | 84 | 27 | 20 |
|  | GEE AR(1), Fay-Graubard | 0 | 100 | 91 | 40 | 19 |
|  | GLMM, likelihood ratio (no interval) | 412 | 13 | 86 | 62 | 24 |
|  | Sign-flip score test (no interval) | 0 | 0 | 97 | 72 | 21 |

## S20. Versions of the convergence check

The check declares a fit of the marginalized model at the boundary when
*nlm* does not report convergence (return code other than 1 or 2) or an
estimated slope of the latent standard deviation is below 0.001 (Section
2.4 of the main text). The primary analysis replaces such a fit by the
fit without the time trend only when *nlm* reports convergence for that
refit and its estimate and standard error are finite; otherwise the
dataset gives no usable result. Under the null, the check declared 71%
of the fits of M3-M5 at the boundary on the main grid, 45% in M6-M7 and
76% on the fourth mechanism, and left 16%, 8.2% and 16% without a
result.<!-- noclaim: descriptive values -->

Table S39 compares the primary rule with four other versions, computed
on the same stored fits by a script written independently of the
pipeline of the tables; for version (iii) it reproduces Table 5. Keeping
the original fit when the refit does not converge, version (ii), leaves
the marginalized model clearly above 7.5% at 4 design points, so that no
strategy with an interval would meet the rule. Accepting a refit that
did not converge, version (i), which an earlier version of this analysis
used, raises the usable share to
91%.<!-- noclaim: descriptive values of the versions of the check -->

**Table S39.** The marginalized model with the convergence check and a t
reference under five versions of the rule for a fit at the boundary,
2,376 design points of the four mechanisms: design points clearly above
7.5% (main grid and M1c, M6-M7, Copula-AR1) and clearly below 2.5%
(all), usable share under the null averaged over design points (%), the
same with two or three occasions and at most 25 clusters, and power per
attempted dataset (%). *nlm* code 3: the last global step failed to
locate a point lower than the estimate, which may be an approximate
minimum.

| Version | Above 7.5%, main and M1c | Above 7.5%, M6-M7 | Above 7.5%, Copula-AR1 | Below 2.5% | Usable (%) | Usable, T ≤ 3 and N ≤ 25 (%) | Power (%) |
|----|----|----|----|----|----|----|----|
| \(iii\) Refit used only when converged; otherwise no result (primary) | 0 | 0 | 0 | 42 | 84 | 28 | 20 |
| \(i\) Refit used when finite, any *nlm* code; otherwise original fit | 0 | 0 | 0 | 24 | 91 | 39 | 21 |
| \(ii\) Refit used only when converged; otherwise original fit | 3 | 0 | 1 | 28 | 91 | 39 | 21 |
| (ii’) As (ii), *nlm* code 3 also accepted | 2 | 0 | 1 | 29 | 91 | 39 | 21 |
| (iii’) As (iii), *nlm* code 3 also accepted | 0 | 0 | 0 | 30 | 85 | 29 | 20 |
