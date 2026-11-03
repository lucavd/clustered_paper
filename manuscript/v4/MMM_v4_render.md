# Choosing an analysis for clustered binary outcomes when clusters are
few and series are short: a simulation study and a clinical application


## Abstract

Binary outcomes measured repeatedly on few subjects are analysed with
generalized estimating equations (GEE), generalized linear mixed models
or marginalized multilevel models. Small-sample corrections to the
sandwich variance have been evaluated mainly by the number of clusters,
with little attention to the number of measurement occasions. We
compared eight analysis strategies over 1,896 simulated design points
(10 to 145 clusters, 2 to 10 occasions, three data-generating
mechanisms), including a scenario anchored to a neonatal lung-ultrasound
cohort and four scenarios in which the latent variability grows strongly
over time, and in subsamples of that cohort. The marginalized multilevel
model was fitted with a convergence check that drops the time trend of
the latent standard deviation when the fit reaches the boundary. With
this check and a t reference distribution, it was the only strategy
whose type I error was never clearly above 7.5% at any design
point.<!-- claim: reg:ra as=descr --> With two or three occasions and up
to 25 clusters no strategy was both available and
calibrated.<!-- claim: reg:ra as=descr --> We recommend the marginalized
model with the convergence check and a t reference, and GEE with the
Fay-Graubard correction with a t reference as a simpler alternative.

**Keywords.** clustered binary data; few clusters; sandwich variance;
small-sample correction; marginalized multilevel model; estimand;
simulation study.

## 1. Introduction

Binary outcomes measured repeatedly within clusters, such as
observations within patients or zones within a lung, are analysed in one
of three ways. A generalized estimating equation (GEE) targets a
population-averaged contrast and treats the dependence as a nuisance. A
generalized linear mixed model (GLMM) targets a cluster-specific
contrast and models the dependence explicitly. A marginalized multilevel
model (MMM) specifies both: a marginal mean model of scientific interest
and a conditional model that generates the dependence \[1,2\]. In
applied work the choice usually follows the conventions of the field
rather than evidence about the design at
hand.<!-- noclaim: background -->

The three approaches target different quantities. Under a logit link the
marginal and conditional contrasts differ. With a normal random
intercept of variance $\sigma^2$, the marginal contrast is approximately
the conditional one multiplied by

$$\left(1 + k^2\sigma^2\right)^{-1/2}, \qquad k = \frac{16\sqrt{3}}{15\pi},$$

so the gap widens with the variance of the random effects \[3,4\]. A
further difference arises when the two are adjusted for different
covariate sets, because the odds ratio is not collapsible \[5\]. Method
comparisons nevertheless often score a mixed model against a marginal
target it does not estimate.

This paper concerns designs with few clusters and short series, where
the asymptotic arguments behind all three approaches are weakest.
Cluster-randomized trials with ten to thirty clusters are common. Among
100 reviewed trials, the 99 that reported the number of clusters had a
median of 25 (interquartile range 15 to 44) \[6\], and 45% of
stepped-wedge trials in another review had fewer than ten clusters
\[7\]. Single-centre cohorts with two or three measurement occasions per
subject are equally common in perinatal and paediatric
research.<!-- noclaim: literature --> The sandwich variance estimator is
anticonservative with few clusters, and a family of small-sample
corrections exists to repair it \[8-10\]. These corrections have been
evaluated on subjects followed over time, with as few as ten subjects
and four occasions \[11-15\]. They have also been evaluated on parallel
cluster-randomized trials, where the simulation evidence has been
reviewed systematically \[16,17\]. In the clustered designs of these
evaluations, the unit within a cluster is a site within a mouth \[8\], a
patient within a general practice \[10\], a participant within a cluster
\[16,18\] or a youth within a community \[19\]. It is never a subject
measured on several occasions. Cluster sizes range from 4 observations
per cluster \[19\] to 150 on average \[16\]. None of these studies fits
a random slope, and only \[10\] generates random regression
coefficients, as a misspecification check. Whether their recommendations
hold with two or three occasions per subject has not been established.

We report a simulation study planned and reported following Morris,
White and Crowther \[20\], and a clinical application. The study has
three components. The first is a simulation over three data-generating
mechanisms and a grid that extends to ten clusters and two occasions.
One scenario takes its prevalence and between-subject variance from a
neonatal lung-ultrasound cohort. The second is a recommendation stated
as a decision rule that must hold at every simulated design point,
because an analyst knows the number of clusters and occasions but not
the mechanism that generated the data. The third is a plasmode
application, in which the cohort is repeatedly subsampled into the
few-cluster regime with the full-cohort estimate as reference. Three
results are established elsewhere and are cited rather than claimed
here: the uncorrected sandwich is anticonservative with few clusters
\[8,9\], the corrections trade level for power \[16,21\], and marginal
and conditional estimands differ
\[5\].<!-- noclaim: aims and cited results -->

## 2. Methods

The study is reported in the order aims, data-generating mechanisms,
estimands, methods and performance measures \[20\]. The revision
reported here followed a protocol in this structure, versioned before
any rerun and released with the code. Table 1 summarizes it. The full
parameter table, extracted from the executed code, is in Table S1.

### 2.1 Aims

The aims concern the treatment-by-time interaction on the logit scale,
in longitudinal binary data with few clusters and short series. We
evaluated (i) the type I error of the test of no interaction; (ii) its
power, raw and size-adjusted; (iii) bias, precision and
confidence-interval coverage of the estimate against the estimand of
each method; and (iv) the probability that a method returns a usable
result. The clinical application and the decision rule of Section 2.6
build on these aims but are not simulation aims.

### 2.2 Data-generating mechanisms

The main grid crossed 16 cluster counts, N ∈ {10, 15, 20, 25, 35, 45,
55, 65, 75, 85, 95, 105, 115, 125, 135, 145}, with 6 maximum numbers of
occasions, T ∈ {2, 3, 4, 6, 8, 10}. It also crossed five mean-model
scenarios (M1-M5) and three mechanisms, giving 1,440 design points. Each
design point was simulated under the alternative (H1) and the null (H0),
and the two arms differed only in the tested coefficient. T is a
maximum. Each subject was observed at all T occasions with probability
0.7, and otherwise at a number of occasions drawn uniformly from $1$ to
$T-1$. At T = 2, between 29% and 31% of subjects therefore contributed a
single observation. Treatment was assigned to subjects with probability
0.5.

The mean model of scenario M1 is
$\mathrm{logit}\,\mu = \beta_0 + \beta_1 t + \beta_2 \mathrm{trt} + \beta_3\, t \times \mathrm{trt}$.
Its marginal coefficients under H1 are (-2, 0.1, 0.2, 0.2), which give a
marginal prevalence of about 0.14 at T = 2. M2 adds a subject-level
covariate (age). M3 and M4 add a random slope on time; in GenBinaryY its
counterpart is a latent standard deviation that changes over time. M4
uses a smaller interaction (0.05). M5 adds a second covariate that
varies over time (dose), with its own random slope.

The three mechanisms differ in how the dependence is generated.

1.  GenBinaryY is a marginalized model with a latent variable and a
    first-order transition term (R package MMLB) \[1,2,22\]. The latent
    standard deviation is a linear predictor, which in M3-M5 grows by
    0.01 per occasion. This is the mechanism of the marginalized model.
2.  GLMM-native is a logistic mixed model with normal random effects.
    The random intercept has variance 1 and, in M3-M5, the random slope
    on time has variance 0.001. The conditional coefficients are the
    marginal ones multiplied by 1.16. This is the mechanism of the mixed
    model.
3.  GLMM-AR1 adds a latent first-order autoregressive process to the
    random effects of GLMM-native (ρ = 0.6, variance 0.5). Its
    conditional coefficients are calibrated on each replicate, so that
    the implied marginal coefficients equal the target. It is the
    mechanism of no method.

A sixth scenario, M1c, uses the mean model of M1 with the prevalence
(0.59) and the between-subject variance (3.5) of the clinical cohort of
Section 2.7. It covers 72 design points, with N up to 45 and T up to 6.
Four further scenarios, M6 and M7 on GenBinaryY and on GLMM-native, test
the convergence check of Section 2.4 where the latent variability grows
strongly over time. In GenBinaryY the latent standard deviation grew by
0.1 or 0.25 per occasion; in GLMM-native the random-slope variance was
0.01 or 0.04. They cover the full grid (384 design points). With 10
occasions the latent standard deviation of GenBinaryY then grows from 1
to 1.9 or 3.25, against 1.09 in M3.

The values of the main scenarios were compared with the clinical cohort
of Section 2.7. Their marginal prevalence (0.14 at T = 2 and 0.28 at T =
10 in M1) is much lower than in the cohort (0.59), and their
random-intercept variance (1) is lower than the patient variance of the
cohort (3.5); scenario M1c was added for this reason. The cohort has no
treatment and few patients with repeated timepoints, so the interaction
per occasion has no direct anchor in it. Its cumulative value at ten
occasions (1.8 on the logit scale) is of the order of the regional
contrasts estimated in the cohort, for example 1.1 between lateral and
anterior zones. The autoregressive parameters of GLMM-AR1 were chosen to
give a serial dependence that is detectable without dominating the
cluster effect.<!-- noclaim: justification of parameter values -->

Each design point received 1,000 replicates per arm. In 44 design
points, all in GLMM-AR1, only 971 to 999 datasets could be generated,
because the per-replicate calibration met a singular system. Replicates
in which one treatment arm was never observed after the first occasion
were excluded from every denominator (247 per arm, the same in the three
mechanisms). Replicate i of scenario m used the seed 20260805 + 1000m +
i under the L’Ecuyer-CMRG generator. The seed was the same for both arms
and for every N, T and mechanism, so the datasets are common random
numbers across cells \[20\].

### 2.3 Estimands and true values

The estimand is the coefficient of the treatment-by-time interaction on
the logit scale. It is marginal for the marginalized model, the GEEs and
the shrinkage estimator, and conditional on the random effects for the
mixed model. In GenBinaryY the marginal true value is the generating
coefficient. There is no conditional true value, so bias and coverage of
the mixed model are not defined on that mechanism. In GLMM-native the
conditional true value is the generating coefficient. The marginal value
was obtained by Gauss-Hermite quadrature on the realized design; under
H0 it is not exactly zero, and it ranges from -0.0030 to -0.0011 across
scenarios and T. In GLMM-AR1 the marginal true value is the target, and
the conditional value was computed for each replicate.

### 2.4 Methods compared

Eight strategies were applied to every dataset (Table 1). The
marginalized multilevel model was fitted with `MMLB::mm` (20 quadrature
points), with the latent-variable and transition formulas of the
GenBinaryY generator. GEEs with exchangeable working correlation were
fitted with `geepack` \[23\]. They were combined with the uncorrected
sandwich variance and with the corrections of Mancl and DeRouen \[8\],
Kauermann and Carroll \[9\] and Fay and Graubard \[10\], computed by
`geesmv`. Each correction ran under a hard time limit of 5 seconds, and
a time-out counted as an unusable fit. The mixed model was fitted with
`lme4::glmer` \[24\] (Laplace approximation). It had a random intercept
in M1-M2 and the random slopes of the generating model in M3-M5. It was
tested by Wald and by a likelihood-ratio test against the model without
the interaction, with the same random effects. The shrinkage estimator
rescales the mixed-model estimate and standard error by
$(1 + k^2 D_{11})^{-1/2}$, with $D_{11}$ the random-intercept variance.
Its Wald statistic therefore equals that of the mixed model, and only
its interval differs. Software versions are in Table S2.

`MMLB::mm` writes the latent standard deviation as
$\exp(a_0) + \exp(a_1)\,t$, so a slope of zero lies on the boundary of
the parameter space. In the main scenarios with a latent standard
deviation that changes over time (M3-M5), the generating slope is close
to that boundary. In an earlier analysis of these data, two failure
modes of the optimizer (`nlm`) produced excess rejections with long
series: non-convergence with a shifted estimate, and convergence with a
collapsed standard error. We therefore added a convergence check. A fit
is declared at the boundary if `nlm` does not report convergence (return
code other than 1 or 2) or if an estimated slope $\exp(a_k)$ is below
0.001. Such a fit is replaced by the fit without the time trend of the
latent standard deviation, and kept if that refit
fails.<!-- noclaim: definition of the check --> The check was specified
after the diagnosis and before the reanalysis of the 1,685,204 fits of
M3-M5 that returned an estimate. Thresholds of 0.0001 and 0.01 were
examined as sensitivity analyses, and scenarios M6-M7 were simulated
afterwards to test the check. The marginalized model without the check
is reported as a sensitivity analysis (Table S3).

Wald tests were referred to the standard normal in the primary tables
and, as an alternative, to a t distribution with N - p degrees of
freedom, p being the number of mean-model coefficients. The
likelihood-ratio test was referred to a $\chi^2$ distribution with one
degree of freedom.

### 2.5 Performance measures

A replicate was usable for a Wald-based method when the fit returned
without error and the estimate, its standard error and its p-value were
finite, with \|estimate\| \< 10 and $0 <$ SE \< 10. Fits that returned
with warnings were kept; excluding them is a sensitivity analysis (Table
S13). For the likelihood-ratio test, both nested fits had to return and
the p-value had to be finite. The probability of a usable result was
reported for every method and cell, with every attempted replicate in
the denominator.

Type I error was the rejection rate at the nominal 5% under H0. Over the
grid we report the pooled rate and two summaries computed within each
mechanism. For each of the 16 cluster counts we took the departure from
5% of the replicate-weighted rate. The worst case is the largest of
these departures, and the mean absolute departure is their mean. At each
design point, a strategy was clearly above 7.5% when its rejection rate
minus two Monte Carlo standard errors exceeded 7.5%. It was clearly
below 2.5% when its rate plus two standard errors fell below 2.5%. These
limits are Bradley’s liberal criterion
\[25\].<!-- noclaim: definitions --> Unconditional validity is the
probability of a usable result multiplied by the probability of not
rejecting a true null given one. Power was reported raw and
size-adjusted. The size-adjusted critical value was the 95th percentile
of \|estimate/SE\| over the usable null replicates of the same cell; for
the likelihood-ratio test it was the 5th percentile of the null
p-values. Bias, empirical and model-based standard errors, root mean
squared error and coverage of the 95% interval were computed against the
estimand of each method.

Monte Carlo standard errors (MCSE) were computed for every measure
\[20\]. With 1,000 usable replicates, the MCSE of a rejection rate at 5%
is 0.69 percentage points. This number of replicates was chosen so that
a single cell resolves a departure from 5% of about 1.4 points. The
datasets are common random numbers across cells. The MCSE of every
summary over cells, and of every difference between methods, was
therefore obtained by a bootstrap that resamples replicate indices
jointly across all cells (2,000 resamples). A comparison between methods
is written as a difference only when its bootstrap interval, the
estimate ± 2 MCSE, lies entirely beyond a relevance margin fixed in
advance. The margin is one percentage point for differences in type I
error or coverage averaged over cells, 2.5 points for a single cell, and
two points for differences in power. Intervals entirely within the
margin are reported as similar. Intervals that straddle it are reported
as differences of about the size of the
margin.<!-- noclaim: definitions --> Counts of miscalibrated cells used
exact binomial tests with the Benjamini-Hochberg and Benjamini-Yekutieli
corrections. In the text, values are given to two significant figures.
Tables give the digit that the MCSE supports, with the MCSE in the
caption.

### 2.6 Decision rule

The recommendation of Figure 1 comes from a rule applied to every
strategy, a strategy being a method with its reference distribution. A
strategy was excluded if its type I error was clearly above 7.5% at any
design point. Among the remaining strategies, the one that returned a
usable result most often was
recommended.<!-- noclaim: definition of the rule --> The rule was
applied to all 1,896 design points. It replaces the criterion of an
earlier version of this work, which selected at each cluster count the
method with the highest unconditional validity (Section
4.1).<!-- noclaim: description of the earlier criterion -->

### 2.7 Clinical application

The cohort comprises 207 neonates with ten lung-ultrasound zones scored
at one to three timepoints. It gives 230 patient-timepoints and 2,296
zone-level outcomes; 22 patients contribute more than one
timepoint.<!-- noclaim: cohort description --> The outcome is severe
loss of aeration (prevalence 0.59). In a logistic mixed model of the
outcome on zone region with a patient random intercept, the patient
variance is 3.5. Adding a patient-timepoint level drives the patient
variance to the boundary and places a variance of 3.8 at the new level
(likelihood-ratio $\chi^2$ = 47 on one degree of freedom). The contrast
of interest is posterior against anterior zones, adjusted for disease,
invasive ventilation and gestational age. It was estimated by the
marginalized model, the GEEs and the mixed model. With 207 clusters the
cohort lies above the simulated range, so it was used as a plasmode. We
drew 500 subsamples of whole patients at each of 10, 15, 20, 25 and 50
patients, with the full-cohort estimate of each method as its reference.
Subsampling whole patients preserves the dependence structure exactly,
and the covariate distribution only in expectation. Of the 25 patients
ventilated invasively, a draw of 10 contains none with probability 0.27,
and the mean model is then not identifiable.

## 3. Results

### 3.1 Type I error on the main grid

Table 2 gives the type I error of the eight methods on the 1,440 design
points of the main grid, under a normal reference. The pooled rejection
rate of the marginalized model with the convergence check ranged from
5.2% to 5.5% across mechanisms. It ranged from 5.3% to 5.8% for
Fay-Graubard and from 5.5% to 6.1% for Mancl-DeRouen. On GenBinaryY it
reached 13% for Kauermann-Carroll and 8.3% for the uncorrected
GEE.<!-- noclaim: descriptive values per method -->

The marginalized model had the smallest worst-case departure from 5% on
all three mechanisms: 0.83, 0.44 and 0.56 percentage points on
GenBinaryY, GLMM-native and GLMM-AR1. Its lead over the next method was
2.8 points on GenBinaryY and about one point on the two random-effects
mechanisms.<!-- claim: cl:C01|GB as=larger; cl:C01|NAT as=about1; cl:C01|AR1 as=about1 -->
Its mean absolute departure (0.24 to 0.54 points) and that of
Fay-Graubard (0.52 to 1.0) were
similar.<!-- claim: cl:C02|GB as=similar; cl:C02|NAT as=similar; cl:C02|AR1 as=similar -->
On matched replicates usable for both methods, the two rejected a true
null at similar
rates.<!-- claim: cl:C03|GB as=similar; cl:C03|NAT as=similar; cl:C03|AR1 as=similar -->

Each method-cell was tested against 5% with the Benjamini-Hochberg
correction within method. The marginalized model was detectably off
nominal in 3% of its cells, Fay-Graubard in 23%, Mancl-DeRouen in 28%
and Kauermann-Carroll in 87% (Table
S5).<!-- claim: cl:C10b|BH as=larger; cl:C10c|BH as=larger --> The
uncorrected GEE was above 5% by 19 points at N = 10 on GenBinaryY. Its
smallest excess over the cluster counts was under half a point on every
mechanism.<!-- claim: cl:C12a|GB as=descr; cl:C12c|GB as=similar; cl:C12c|NAT as=similar; cl:C12c|AR1 as=similar -->
Above 45 clusters, Fay-Graubard and the marginalized model were both
within about half a point of 5% on the random-effects
mechanisms.<!-- claim: cl:C25a|NAT as=descr; cl:C25b|NAT as=similar; cl:C25b|AR1 as=similar -->

The ordering of the methods was stable across mechanisms, with Spearman
correlations of 0.90 to 0.98 between pairs of mechanisms. GLMM-AR1
differed from GLMM-native by 0.19 percentage points on average across
methods.<!-- claim: cl:C08|rate as=descr; cl:C09|all as=descr --> The
mean absolute departure of the marginalized model fell from 1.5 points
at two occasions to 0.71 at four or more (Table S4). At two occasions
the five marginal methods departed by 1.5 to 2.5 points on average. The
conditional methods departed by 5.6 (mixed model, Wald) to 11 points
(shrinkage GLM).<!-- claim: cl:C13d|all as=larger -->

### 3.2 The convergence check of the marginalized model

In scenarios M3-M5 the check replaced 68% of the fits of the
marginalized model with the fit without the time trend of the latent
standard deviation. Without the check, the marginalized model was
clearly above 7.5% at 138 of the 864 null design points of these
scenarios, 84 of them with eight or ten occasions. With the check it was
clearly above 7.5% at
0.<!-- claim: reg:cc_hi_study as=descr; reg:cc_hi_cc as=descr --> Its
mean absolute departure at ten occasions fell from 2.8 to 0.67 points.
Size-adjusted power at eight or ten occasions changed by 0.31 points on
average. The thresholds 0.0001 and 0.01 in place of 0.001 changed the
class of none of the 156 comparisons classified for this paper. With
either threshold, the marginalized model with a t reference was not
clearly above 7.5% at any design
point.<!-- claim: reg:thr_pa6_t4 as=descr; reg:thr_pa6_t2 as=descr -->

In scenarios M6-M7 the fit without the time trend is misspecified. There
the check intervened in 36% to 51% of fits under H0, and in 58% to 60%
at two or three occasions, where `nlm` most often failed to converge.
With a t reference, the marginalized model with the check was never
clearly above 7.5% at the 384 design points. Under a normal reference
its mean type I error was 5.2% to 5.8% across the four scenarios.
Without the check, and with a t reference, it was clearly above 7.5% at
8 design points. Fay-Graubard with a t reference was clearly above 7.5%
at 2 design points, all on GenBinaryY with the steeper slope (M7); its
largest rate was 10% (Table
S14).<!-- claim: reg:ns_mmmt_hi as=descr; reg:ns_nocc_t_hi as=descr; reg:ns_fgt_hi as=descr; reg:ns_fgt_g7_hi equals reg:ns_fgt_hi -->

### 3.3 Availability of a result

The probability of a usable result depended mostly on the number of
occasions (Table 3). With random slopes and two occasions, the mixed
model returned a usable Wald test in 0.002% of attempts. The
likelihood-ratio test with the same random effects was usable in 0.14%,
against 71% for the marginalized model and 73% for
Fay-Graubard.<!-- claim: cl:C13a|all as=descr --> With a random
intercept only, the likelihood-ratio test was usable in 99% of attempts
at two occasions, against 69% for the Wald
test.<!-- claim: cl:C13c|all as=larger --> From four occasions every
method returned a usable result in more than 90% of
attempts.<!-- claim: reg:pu_T4_min as=descr -->

At 20 clusters or fewer, the likelihood-ratio test returned a usable
result in 86% of attempts and the uncorrected GEE in 79%. The
marginalized model (68%) and the three corrected sandwiches (each 69%)
were similar.<!-- claim: cl:C14a|all as=larger; cl:C14d|all as=descr -->
Unconditional validity ranged from 80% for the likelihood-ratio test to
61% for the Wald test of the mixed
model.<!-- claim: cl:C14b|all as=larger -->

### 3.4 Reference distribution and the decision rule

Referring the Wald statistic to a t distribution with N - p degrees of
freedom lowered every pooled rejection rate. The shift was 0.53 to 1.1
percentage points over the grid and 1.4 to 3.7 points at 25 clusters or
fewer (Table S11).<!-- claim: cl:C17b|all as=larger --> Under the t
reference the marginalized model became conservative at ten clusters.
There it was clearly below 2.5% at 47 of 102 design points, against 18
for
Fay-Graubard.<!-- claim: reg:rs_mmmt_N10_lo as=descr; reg:rs_fgt_N10_lo as=descr -->
Fay-Graubard also had the smaller worst-case departure under the t
reference, 0.9 to 1.1 points against 3.0 to
3.4.<!-- claim: cl:C17d|GB as=larger; cl:C17d|NAT as=larger; cl:C17d|AR1 as=larger -->

Table 4 applies the rule of Section 2.6 to the 1,896 design points of
the main grid, M1c and M6-M7. The marginalized model with the check and
a t reference was never clearly above 7.5%. Fay-Graubard with a t
reference was clearly above 7.5% at 2 design points and Mancl-DeRouen
with a t reference at 32. The likelihood-ratio test was clearly above
7.5% at 216 design points, and the Wald test of the mixed model with a t
reference at 265.<!-- claim: reg:ra as=descr --> The recommended
strategy returned a usable result in 91% of datasets, against 92% for
Fay-Graubard with a t reference. It was clearly below 2.5% at 77 design
points, against 116. With two or three occasions and at most 25
clusters, both returned a usable result in fewer than half of the
datasets (40% and 44%). Both rejected a false null in fewer than 2% of
attempted
datasets.<!-- claim: reg:ra as=descr; reg:upow_T23_ceiling as=descr -->

### 3.5 Power and coverage

Raw power was highest for Kauermann-Carroll (35% on GenBinaryY), which
was also the most anticonservative method. The ordering by raw power did
not survive size adjustment (Table
5).<!-- claim: cl:C16c|GB as=larger; cl:C16d|GB as=larger --> After size
adjustment, power ranged from 21% to 27% across methods and
mechanisms.<!-- noclaim: range --> On GenBinaryY the Wald test of the
mixed model had more size-adjusted power than Fay-Graubard, by 2.2
points with the default optimizer of lme4 and by 3.1 points with bobyqa.
On GLMM-AR1 the difference was about two points with the default
optimizer, and the two were similar with bobyqa. On GLMM-native the two
were similar with the default optimizer, and with bobyqa Fay-Graubard
had more power, by about two
points.<!-- claim: cl:C16a|GB as=larger; cl:C16a|GB:bobyqa as=larger; cl:C16a|AR1 as=about1; cl:C16a|AR1:bobyqa as=similar; cl:C16a|NAT as=similar; cl:C16a|NAT:bobyqa as=about1 -->

Coverage of the 95% interval was 94% to 95% for Fay-Graubard,
Mancl-DeRouen and the marginalized model on every mechanism (Table S12).
It was furthest from nominal for Kauermann-Carroll, at 86% to
89%.<!-- claim: cl:C20a|GB as=similar; cl:C20a|NAT as=similar; cl:C20a|AR1 as=similar; cl:C20b|GB as=larger -->

### 3.6 Scenario anchored to the cohort

Scenario M1c has the prevalence and between-subject variance of the
cohort (72 design points). Two results of the main grid did not hold
there. Under a normal reference the marginalized model rejected a true
null in 6.8% to 7.0% of usable replicates, against 5.0% to 5.1% for
Fay-Graubard. Its worst-case departure was larger than that of the best
method by about one to two
points.<!-- claim: cm:C03|GB as=larger; cm:C03|NAT as=larger; cm:C03|AR1 as=larger; cm:C01|GB as=about1; cm:C01|NAT as=about1; cm:C01|AR1 as=about1 -->
Under a t reference its rate was 5.0% to 5.1%, and that of Fay-Graubard
3.5% to 3.6%.<!-- noclaim: descriptive values --> The size-adjusted
power of the Wald test of the mixed model fell below that of
Fay-Graubard, by 4.3 to 4.9 points with the default optimizer. With
bobyqa the deficit was about two points on GenBinaryY and GLMM-AR1, and
the two were similar on GLMM-native. The default optimizer returned
degenerate fits with a near-zero standard
error.<!-- claim: cm:C16a|GB as=larger; cm:C16a|NAT as=larger; cm:C16a|AR1 as=larger; cm:C16a|GB:bobyqa as=about1; cm:C16a|NAT:bobyqa as=similar; cm:C16a|AR1:bobyqa as=about1 -->
Which method had the fewest or the most miscalibrated cells could not be
distinguished in
M1c.<!-- claim: cm:C10b|BH as=notdist; cm:C10c|BH as=notdist -->

### 3.7 Clinical application

In the full cohort the posterior-versus-anterior contrast was 0.38 (SE
0.17) by the marginalized model and 0.38 by GEE. The mixed model gave
0.48 (SE 0.12) on its conditional scale. In subsamples of 10 patients,
22% of draws contained no invasively ventilated patient, so that the
mean model could not be fitted (Table S15). Among the fits that were
returned, a usable result was obtained by the mixed model in 100%, by
the marginalized model in 100% and by the uncorrected GEE in 92%. It was
obtained by Fay-Graubard in 45% and by Mancl-DeRouen in 20%; most
failures of the corrections returned a point estimate without a standard
error.<!-- noclaim: descriptive values per method --> Across subsample
sizes, the standard error of the mixed model was 0.53 to 0.65 times the
empirical standard deviation of its estimates. The corresponding ratio
was 0.87 to 1.1 for the marginalized model. The bias of the marginalized
model against its full-cohort estimate fell from 0.16 at 10 patients to
0.015 at 50.

## 4. Discussion

### 4.1 The recommendation

The rule of Section 2.6 selects the marginalized multilevel model with
the convergence check and a t reference with N - p degrees of freedom,
at every number of clusters and occasions (Figure 1). Across the 1,896
simulated design points it was the only strategy whose type I error was
never clearly above 7.5%. It returned a usable result in 91% of
datasets.<!-- claim: reg:ra as=descr --> Fay-Graubard with a t reference
returned a usable result about as often and needs no convergence check,
so it is a reasonable alternative. It was clearly above 7.5% only with
at most 25 clusters and a latent variability that grows steeply over
time.<!-- claim: reg:fgt_hi_Nmax as=descr --> The two strategies are
conservative in different places: the marginalized model at ten
clusters, Fay-Graubard more often elsewhere (Table
4).<!-- claim: reg:rs_mmmt_N10_lo as=descr -->

The rule asks that a recommendation hold at every design point, because
an analyst knows the number of clusters and occasions but not the
mechanism that generated the data. It replaces a criterion that
maximized unconditional validity at each cluster count. That criterion
rewards a method for answering, whatever its level. On the present data
it would have recommended the likelihood-ratio test of the mixed model
up to 25 clusters, a test that was clearly above 7.5% at 216 design
points.<!-- claim: reg:switch_GB as=descr; reg:ra as=descr --> The
change of criterion also reverses the general recommendation of the
earlier version, which was Fay-Graubard. The reason is that scenarios
M6-M7 exposed design points at which Fay-Graubard with a t reference was
clearly above
7.5%.<!-- noclaim: statement of the change of criterion -->

With two or three occasions and up to 25 clusters, no strategy was both
available and calibrated. The recommended strategy and Fay-Graubard
returned a usable result in fewer than half of the datasets and almost
never detected the simulated effect.<!-- claim: reg:ra as=descr --> With
a random slope and two occasions, the random-effects covariance of the
mixed model is not identifiable, since two subject-level parameters
saturate two observations, and its Wald test was essentially never
usable. The case for the maximal random-effects structure \[26\] is thus
bounded by power \[27\] and, with two occasions, by identifiability. In
this corner the design determines what can be learned more than the
analysis does. Four or more occasions per subject, or more clusters,
should be planned where
possible.<!-- noclaim: interpretation and design advice -->

### 4.2 The convergence check

The check addresses a property of the parametrization rather than of the
marginalized model. `MMLB::mm` writes the latent standard deviation with
an exponentiated slope, so a slope near zero lies on the boundary of the
parameter space, and in all the main scenarios the generating slope was
close to zero. At the boundary `nlm` either stopped without converging,
with an estimate shifted by an amount that did not shrink with the
standard error, or converged with a collapsed standard error. Both
failures produced excess rejections that grew with the number of
clusters and occasions. This is why the marginalized model appeared
miscalibrated with long series in an earlier analysis. Replacing such
fits with the fit without the time trend removed the
excess.<!-- noclaim: interpretation of the mechanism --> The check
intervened far more often than we had anticipated, mainly because `nlm`
failed to converge with two or three occasions or few
clusters.<!-- claim: reg:pb1_n6 as=descr --> Where it intervened, the
type I error stayed close to nominal. It also remained calibrated in
M6-M7, where the fit without the time trend is misspecified. Analysts
who use the marginalized model with a latent standard deviation that
changes over time should apply the check and report how often it
intervened.<!-- noclaim: recommendation -->

### 4.3 Relation to previous work

We reproduce the anticonservatism of the uncorrected sandwich with few
clusters \[8,9\]. Li and Redden found that Kauermann-Carroll attains the
nominal level at ten clusters under a t reference \[16\]. In our data
its rate at ten clusters was 3.9% to 5.5%. Over the whole grid, however,
it was clearly above 7.5% at 951 of 1,896 design points, so the finding
does not extend beyond the cluster range in which it was
established.<!-- claim: cl:C24a|AR1 as=descr; cl:C24a|GB as=descr; reg:ra as=descr -->
Mancl and DeRouen reported sizes close to nominal with ten subjects
\[8\], and Bie et al. favoured GEE with their correction \[18\]. Under a
t reference its pooled rate was 4.7% to 5.3%, inside their band. It was
nevertheless clearly above 7.5% at 32 design points, 18 of them in
scenarios
M6-M7.<!-- claim: cl:C24c|NAT as=descr; cl:C24c|GB as=descr; reg:ra as=descr -->
The marginalized model has been studied for moderate to long series
\[22\], and the short series studied here complement that work.
Estimation methods for mixed logistic models with few clusters are
reviewed by McNeish \[28\].

### 4.4 Limitations

In the main scenarios the random-slope variance and the time trend of
the latent standard deviation were small. Conclusions about random
slopes therefore rest mainly on scenarios M6-M7, which cover one
structure (a slope on time) and two mechanisms. The working correlation
of the GEEs was exchangeable throughout, and misspecified working
correlations were not simulated. The degrees of freedom of the t
reference were N - p and were not varied, whereas Li and Redden used K -
2 \[16\]. GLMM-AR1 behaved almost exactly like GLMM-native, so the study
rests on two distinct mechanisms. Each of these is the home model of one
method: GenBinaryY of the marginalized model and GLMM-native of the
mixed model. In GLMM-native the marginal value under the null is not
exactly zero. It ranges from -0.0030 to -0.0011 in the main scenarios
and reaches -0.0079 in M6-M7. The convergence check was introduced after
part of the results had been seen, and was then tested on scenarios
simulated afterwards. Small-sample corrections were computed with
`geesmv`, whose numerical failures at small samples are documented in
Section S4. An independent implementation, such as the CR2 estimator of
`clubSandwich` \[29,30\], would avoid them. The clinical application
uses one cohort and one outcome. The subsampling resamples whole
patients, and the patient-timepoint level that carries most of the
dependence is identified by only 22 patients measured more than
once.<!-- noclaim: limitations -->

## 5. Conclusions

With few clusters and short series, the choice of analysis for a
clustered binary outcome decides whether the test keeps its nominal
level and whether it returns an answer at all. The marginalized
multilevel model was the only strategy that kept its level at every
simulated design point, provided that it was fitted with a check that
drops the time trend of the latent standard deviation at the boundary
and referred to a t distribution. With two or three occasions and up to
25 clusters, no strategy was both available and calibrated, so in that
regime the design limits what any analysis can
show.<!-- claim: reg:ra as=descr -->

We recommend the marginalized model with the convergence check and a t
reference. GEE with the Fay-Graubard correction with a t reference is a
simpler alternative, to be used with caution when clusters are few and
the between-subject variability grows over
time.<!-- noclaim: recommendation -->

## Code and data availability

The simulation code, the per-replicate results and the scripts that
produce every number, table and figure of this paper are available at
https://github.com/lucavd/MMM (commit f9ab317). They will be archived on
Zenodo with a DOI at submission. The manuscript source reads every
number in the text from a registry of values, each with its source file
and checksum. The clinical cohort cannot be shared publicly; the
plasmode results are included as summary tables.

## Declarations

\[Authors: funding, conflicts of interest, ethics approval of the cohort
and an ICMJE-conformant statement on the use of any software tool in
drafting, language editing or consistency checking are to be completed
by the authors.\]

## References

1.  Heagerty PJ. Marginally specified logistic-normal models for
    longitudinal binary data. *Biometrics* 1999;55(3):688-698.
    doi:10.1111/j.0006-341X.1999.00688.x
2.  Heagerty PJ, Zeger SL. Marginalized multilevel models and likelihood
    inference. *Statistical Science* 2000;15(1):1-26.
    doi:10.1214/ss/1009212671
3.  Zeger SL, Liang KY, Albert PS. Models for longitudinal data: a
    generalized estimating equation approach. *Biometrics*
    1988;44(4):1049-1060. doi:10.2307/2531734
4.  Heagerty PJ, Kurland BF. Misspecified maximum likelihood estimates
    and generalised linear mixed models. *Biometrika*
    2001;88(4):973-985. doi:10.1093/biomet/88.4.973
5.  Daniel R, Zhang J, Farewell D. Making apples from oranges: comparing
    noncollapsible effect estimators and their standard errors after
    adjustment for different covariate sets. *Biometrical Journal*
    2021;63(3):528-557. doi:10.1002/bimj.201900297
6.  Kahan BC, Forbes G, Ali Y, et al. Increased risk of type I errors in
    cluster randomized trials with small or medium numbers of clusters:
    a review, reanalysis, and simulation study. *Trials* 2016;17:438.
    doi:10.1186/s13063-016-1571-2
7.  Barker D, D’Este C, Campbell MJ, McElduff P. Minimum number of
    clusters and comparison of analysis methods for cross sectional
    stepped wedge cluster randomized trials with binary outcomes: a
    simulation study. *Trials* 2017;18:119.
    doi:10.1186/s13063-017-1862-2
8.  Mancl LA, DeRouen TA. A covariance estimator for GEE with improved
    small-sample properties. *Biometrics* 2001;57(1):126-134.
    doi:10.1111/j.0006-341X.2001.00126.x
9.  Kauermann G, Carroll RJ. A note on the efficiency of sandwich
    covariance matrix estimation. *Journal of the American Statistical
    Association* 2001;96(456):1387-1396. doi:10.1198/016214501753382309
10. Fay MP, Graubard BI. Small-sample adjustments for Wald-type tests
    using sandwich estimators. *Biometrics* 2001;57(4):1198-1206.
    doi:10.1111/j.0006-341X.2001.01198.x
11. Pan W, Wall MM. Small-sample adjustments in using the sandwich
    variance estimator in generalized estimating equations. *Statistics
    in Medicine* 2002;21(10):1429-1441. doi:10.1002/sim.1142
12. Morel JG, Bokossa MC, Neerchal NK. Small sample correction for the
    variance of GEE estimators. *Biometrical Journal*
    2003;45(4):395-409. doi:10.1002/bimj.200390021
13. Westgate PM, Burchett WW. Improving power in small-sample
    longitudinal studies when using generalized estimating equations.
    *Statistics in Medicine* 2016;35(21):3733-3744. doi:10.1002/sim.6967
14. Ford WP, Westgate PM. A comparison of bias-corrected empirical
    covariance estimators with generalized estimating equations in
    small-sample longitudinal study settings. *Statistics in Medicine*
    2018;37(28):4318-4329. doi:10.1002/sim.7917
15. Gosho M, Ishii R, Noma H, Maruo K. A comparison of bias-adjusted
    generalized estimating equations for sparse binary data in
    small-sample longitudinal studies. *Statistics in Medicine*
    2023;42(15):2711-2727. doi:10.1002/sim.9744
16. Li P, Redden DT. Small sample performance of bias-corrected sandwich
    estimators for cluster-randomized trials with binary outcomes.
    *Statistics in Medicine* 2015;34(2):281-296. doi:10.1002/sim.6344
17. Hemming K, Thompson J, Kristunas C, Forbes A, Watson S, Bowden R,
    McKenzie JE, Kasza J. The performance of small sample correction
    methods for controlling type I error when analyzing parallel cluster
    randomized trials: a systematic review of simulation studies.
    *Journal of Clinical Epidemiology* 2025;185:111838.
    doi:10.1016/j.jclinepi.2025.111838
18. Bie R, Haneuse S, Huey N, Schildcrout J, Shepherd BE. Fitting
    marginal models in small samples: a simulation study of marginalized
    multilevel models and generalized estimating equations. *Statistics
    in Medicine* 2021;40(24):5298-5312. doi:10.1002/sim.9126
19. Lu B, Preisser JS, Qaqish BF, Suchindran C, Bangdiwala SI,
    Wolfson M. A comparison of two bias-corrected covariance estimators
    for generalized estimating equations. *Biometrics*
    2007;63(3):935-941. doi:10.1111/j.1541-0420.2007.00764.x
20. Morris TP, White IR, Crowther MJ. Using simulation studies to
    evaluate statistical methods. *Statistics in Medicine*
    2019;38(11):2074-2102. doi:10.1002/sim.8086
21. Leyrat C, Morgan KE, Leurent B, Kahan BC. Cluster randomized trials
    with a small number of clusters: which analyses should be used?
    *International Journal of Epidemiology* 2018;47(1):321-331.
    doi:10.1093/ije/dyx169
22. Schildcrout JS, Heagerty PJ. Marginalized models for moderate to
    long series of longitudinal binary response data. *Biometrics*
    2007;63(2):322-331. doi:10.1111/j.1541-0420.2006.00680.x
23. Halekoh U, Højsgaard S, Yan J. The R package geepack for generalized
    estimating equations. *Journal of Statistical Software*
    2006;15(2):1-11. doi:10.18637/jss.v015.i02
24. Bates D, Mächler M, Bolker B, Walker S. Fitting linear mixed-effects
    models using lme4. *Journal of Statistical Software*
    2015;67(1):1-48. doi:10.18637/jss.v067.i01
25. Bradley JV. Robustness? *British Journal of Mathematical and
    Statistical Psychology* 1978;31(2):144-152.
    doi:10.1111/j.2044-8317.1978.tb00581.x
26. Barr DJ, Levy R, Scheepers C, Tily HJ. Random effects structure for
    confirmatory hypothesis testing: keep it maximal. *Journal of Memory
    and Language* 2013;68(3):255-278. doi:10.1016/j.jml.2012.11.001
27. Matuschek H, Kliegl R, Vasishth S, Baayen H, Bates D. Balancing type
    I error and power in linear mixed models. *Journal of Memory and
    Language* 2017;94:305-315. doi:10.1016/j.jml.2017.01.001
28. McNeish D. Estimation methods for mixed logistic models with few
    clusters. *Multivariate Behavioral Research* 2016;51(6):790-804.
    doi:10.1080/00273171.2016.1236237
29. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
    variance estimation and hypothesis testing in fixed effects models.
    *Journal of Business & Economic Statistics* 2018;36(4):672-683.
    doi:10.1080/07350015.2016.1247004
30. Pustejovsky JE, Tipton E. Corrigendum: small sample methods for
    cluster-robust variance estimation and hypothesis testing in fixed
    effects models. *Journal of Business & Economic Statistics*
    2023;41(2):650-652. doi:10.1080/07350015.2023.2174123

## Tables and figures

**Table 1.** Summary of the simulation study in the ADEMP structure
\[20\].

|  |  |
|----|----|
| Aims | Type I error and power of the test of no treatment-by-time interaction, bias and coverage of its estimate, and probability of a usable result, with few clusters and short series |
| Data-generating mechanisms | GenBinaryY (marginalized latent-variable transition model), GLMM-native (logistic mixed model), GLMM-AR1 (mixed model with latent first-order autoregression); scenarios M1-M5 on the full grid of 16 cluster counts (10 to 145) and 6 series lengths (2 to 10), 1,440 design points; M1c anchored to the clinical cohort (72 design points); M6-M7 with a steep time trend of the latent variability (384 design points); null and alternative arm for each, 1,000 replicates per arm |
| Estimands | Interaction coefficient on the logit scale: marginal (marginalized model, GEEs, shrinkage GLM) or conditional on the random effects (GLMM) |
| Methods | Marginalized multilevel model (`MMLB::mm`) with convergence check; GEE with uncorrected sandwich and with the Mancl-DeRouen, Kauermann-Carroll and Fay-Graubard corrections (`geepack`, `geesmv`); GLMM with Wald and likelihood-ratio tests (`lme4`); shrinkage-rescaled GLMM. Wald tests with normal and t (N - p degrees of freedom) reference |
| Performance measures | Probability of a usable result; type I error (pooled, worst case and mean absolute departure over cluster counts, clearly above 7.5% or below 2.5% per design point); raw and size-adjusted power; bias, empirical and model standard error, coverage; Monte Carlo standard errors by block bootstrap over replicates |

**Table 2.** Type I error at the nominal 5% on the main grid (1,440
design points), normal reference. Rate: pooled rejection rate (%).
Worst: largest absolute departure from 5% over the 16 cluster counts
(percentage points). MAD: mean absolute departure over the cluster
counts. GB = GenBinaryY, NAT = GLMM-native, AR1 = GLMM-AR1. Monte Carlo
standard errors (block bootstrap) at most 0.052 for rates and MAD and at
most 0.30 for the worst case; values of the marginalized model without
the check in Table S3.

| Method | Rate GB | Rate NAT | Rate AR1 | Worst GB | Worst NAT | Worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 4.7 | 2.8 | 2.9 | 1.0 | 0.5 | 0.6 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 8.6 | 5.9 | 6.0 | 1.5 | 0.9 | 1.0 |
| Marginalized model, with check | 5.5 | 5.2 | 5.4 | 0.8 | 0.4 | 0.6 | 0.5 | 0.2 | 0.4 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 4.5 | 1.8 | 2.2 | 3.5 | 1.4 | 1.6 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 4.9 | 2.3 | 2.4 | 4.2 | 1.6 | 1.9 |
| GEE, uncorrected | 8.3 | 7.6 | 7.7 | 21.0 | 17.4 | 18.3 | 3.8 | 3.0 | 3.1 |
| GLMM, likelihood ratio | 7.8 | 5.8 | 6.1 | 3.6 | 1.8 | 1.8 | 2.8 | 0.8 | 1.1 |
| GEE, Kauermann-Carroll | 13.1 | 10.4 | 10.7 | 9.2 | 6.5 | 6.7 | 7.9 | 5.3 | 5.6 |

**Table 3.** Replicates returning a usable result (%), by structure of
the random effects and number of occasions T, pooled over mechanisms and
cluster counts, every attempted replicate in the denominator. LR =
likelihood-ratio test with the random-effects structure of the Wald
test. Monte Carlo standard errors at most 0.13 points. All methods in
Table S7.

| Random effects | Occasions | GLMM, Wald | GLMM, LR | Marginalized model | GEE, Fay-Graubard | GEE, uncorrected |
|----|----|----|----|----|----|----|
| Intercept (M1-M2) | T = 2 | 68.6 | 98.9 | 68.2 | 72.7 | 78.7 |
|  | T = 3 | 87.8 | 99.6 | 87.1 | 88.3 | 91.7 |
|  | T = 4 | 93.7 | 99.9 | 92.9 | 93.8 | 95.9 |
|  | T = 6 | 97.8 | 100.0 | 97.4 | 97.7 | 98.6 |
|  | T = 10 | 99.6 | 100.0 | 99.5 | 99.5 | 99.7 |
| Intercept and slopes (M3-M5) | T = 2 | 0.0 | 0.1 | 71.0 | 73.1 | 78.8 |
|  | T = 3 | 56.7 | 66.5 | 87.9 | 88.4 | 91.8 |
|  | T = 4 | 90.8 | 99.4 | 93.8 | 94.1 | 96.1 |
|  | T = 6 | 97.7 | 100.0 | 97.8 | 98.0 | 98.8 |
|  | T = 10 | 99.6 | 100.0 | 99.4 | 99.5 | 99.7 |

**Table 4.** The decision rule applied to 1,896 design points: the main
grid and M1c (1,512) and scenarios M6-M7 (384). Clearly above 7.5%:
rejection rate minus two Monte Carlo standard errors above 7.5% under
the null; clearly below 2.5%: rate plus two standard errors below 2.5%;
design points with fewer than 20 usable replicates are not assessed.
Usable: replicates returning a usable result (%), averaged over design
points; power: rejections of a false null per attempted dataset (%). The
column at ten clusters covers the main grid and M1c (102 design points);
the last column covers two or three occasions and at most 25 clusters.

| Strategy | Above 7.5%, main and M1c | Above 7.5%, M6-M7 | Below 2.5%, all | Below 2.5%, ten clusters | Usable (%) | Power (%) | Usable, T ≤ 3 and N ≤ 25 (%) |
|----|----|----|----|----|----|----|----|
| Marginalized model with check, t | 0 | 0 | 77 | 47 | 91.0 | 21.2 | 40.4 |
| Marginalized model with check, normal | 1 | 3 | 18 | 7 | 91.0 | 22.5 | 40.4 |
| GEE, Fay-Graubard, t | 0 | 2 | 116 | 18 | 91.9 | 20.5 | 43.6 |
| GEE, Fay-Graubard, normal | 36 | 28 | 88 | 12 | 91.9 | 21.7 | 43.6 |
| GEE, Mancl-DeRouen, t | 14 | 18 | 124 | 17 | 91.8 | 20.6 | 43.5 |
| GEE, Kauermann-Carroll, t | 727 | 224 | 48 | 14 | 91.9 | 29.1 | 43.5 |
| GLMM, Wald, t | 187 | 78 | 18 | 14 | 80.7 | 23.2 | 30.2 |
| GLMM, likelihood ratio | 216 | 0 | 13 | 4 | 86.4 | 24.4 | 62.7 |

**Table 5.** Raw and size-adjusted power (%) at the interaction contrast
on the main grid, pooled over design points, by mechanism. Size
adjustment uses the 95th percentile of \|estimate/SE\| under the null of
the same cell (for the likelihood-ratio test, the 5th percentile of the
null p-values). Monte Carlo standard errors at most 0.27 points.

| Method | Raw GB | Raw NAT | Raw AR1 | Adjusted GB | Adjusted NAT | Adjusted AR1 |
|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 23.3 | 23.7 | 24.1 | 22.0 | 23.1 | 23.5 |
| GEE, Mancl-DeRouen | 23.5 | 23.9 | 24.4 | 21.9 | 23.1 | 23.4 |
| Marginalized model, with check | 24.1 | 24.4 | 25.2 | 23.4 | 24.2 | 24.6 |
| GLMM, Wald | 31.9 | 28.0 | 29.2 | 24.2 | 23.5 | 25.4 |
| Shrinkage GLM | 32.2 | 28.1 | 29.3 | 23.6 | 23.2 | 25.0 |
| GEE, uncorrected | 26.2 | 26.5 | 27.0 | 21.5 | 22.8 | 23.0 |
| GLMM, likelihood ratio | 30.1 | 26.8 | 27.9 | 26.4 | 26.3 | 26.7 |
| GEE, Kauermann-Carroll | 34.6 | 31.9 | 32.7 | 23.0 | 24.1 | 24.5 |

<img src="figures/fig1.png" style="width:100.0%" />

**Figure 1.** Recommended analysis for the test of the treatment-by-time
interaction with few clusters. (a) The rule: the marginalized multilevel
model with the convergence check and a t reference (N - p degrees of
freedom) at any number of clusters and occasions, with Fay-Graubard and
a t reference as alternative; N = number of clusters, p = number of
mean-model coefficients. (b) Mean absolute departure from the nominal 5%
(percentage points) by number of occasions, averaged over scenarios
M1-M5, cluster counts and mechanisms, logarithmic axis; bars are
estimate ± 2 Monte Carlo standard errors. (c) Percentage of design
points at which each strategy is clearly above 7.5%, by number of
clusters, over the 1,512 design points of the main grid and M1c;
scenarios M6-M7 are not drawn (Table 4). (d) Worst-case against mean
departure from 5% on GenBinaryY, over the 96 combinations of N and T,
normal reference. FG = Fay-Graubard, MD = Mancl-DeRouen, LRT =
likelihood-ratio test of the GLMM, KC = Kauermann-Carroll, MMM =
marginalized multilevel model.

<img src="figures/fig2.png" style="width:100.0%" />

**Figure 2.** Type I error at every design point of the main grid and
M1c (nested-loop plot), normal reference, marginalized model with the
convergence check. Design points are ordered by scenario (outer loop),
number of occasions T and number of clusters N (inner loop, bottom
panel); rows are the mechanisms. Grey band: 2.5% to 7.5%; line: nominal
5%. Points with fewer than 20 usable replicates are not drawn and the
axis is truncated at 30%.
