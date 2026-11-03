# Choosing an analysis for clustered binary outcomes when clusters are
few and series are short: a simulation study and two clinical
applications


## Abstract

Binary outcomes measured repeatedly on few subjects are analysed with
generalized estimating equations (GEE), mixed models or marginalized
multilevel models. Small-sample corrections to the sandwich variance
have been evaluated mainly by number of clusters, rarely by number of
occasions. We compared analysis strategies over 2,376 simulated design
points (10-145 clusters, 2-10 occasions, four mechanisms, one purely
serial) and in subsamples of two clinical datasets. GEE with a
first-order autoregressive working correlation, the Fay-Graubard
correction and a t reference, and the marginalized multilevel model with
a convergence check and a t reference, were never clearly above 7.5%
type I error at any design point and returned a usable result in 92% and
91% of datasets.<!-- claim: reg:t4 as=descr --> Their calibration
differed materially only at ten clusters, where GEE was closer to
nominal.<!-- claim: reg:hb as=descr --> The mixed-model Wald test was
anticonservative under the Laplace approximation; with adaptive
quadrature its excess remained only with ten
occasions.<!-- noclaim: summary of Sections 3.6 and 3.10 --> With two or
three occasions and up to 25 clusters both strategies answered in fewer
than half of datasets; a sign-flip score test answered in 72% and kept
its level, but gives no estimate.<!-- claim: reg:t4 as=descr --> We
recommend the GEE strategy, which needs no convergence check; the
marginalized model is a comparable
alternative.<!-- noclaim: recommendation -->

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
covariate sets, because the odds ratio is not collapsible \[5\]. A
method comparison can nevertheless score a mixed model against a
marginal target it does not estimate.

This paper concerns designs with few clusters and short series, where the asymptotic arguments behind all three approaches are weakest. Cluster-randomized trials with ten to thirty clusters are common \[6\], and many stepped-wedge trials have fewer than ten \[7\].<!-- noclaim: literature --> The sandwich variance estimator is anticonservative with few clusters, and a family of small-sample corrections exists to repair it \[8-10\]. These corrections have been evaluated on subjects followed over time, with as few as ten subjects and four occasions \[11-15\]. They have also been evaluated on parallel cluster-randomized trials, where the simulation evidence has been reviewed systematically \[16,17\]. In the original evaluations and in the cluster-randomized comparisons, the unit within a cluster is never a subject measured on several occasions \[8,10,16,18,19\]. None of these studies fits a random slope, and only Fay and Graubard \[10\] generate random regression coefficients, as a misspecification check. Whether their recommendations hold with two or three occasions per subject has not been established.

We report a simulation study planned and reported following Morris, White and Crowther \[20\], and two clinical applications. The study has three components. The first is a simulation over four data-generating mechanisms and a grid that extends to ten clusters and two occasions. Three mechanisms contain a component that is constant within subject, as the methods compared assume, and the fourth has only serial dependence. One scenario takes its prevalence and between-subject variance from a neonatal lung-ultrasound cohort. The second is a recommendation stated as a decision rule that must hold at every simulated design point. The third is a plasmode application to the cohort and to a public respiratory trial, each repeatedly subsampled into the few-cluster regime with the full-data estimate as reference.<!-- noclaim: aims -->

## 2. Methods

The study is reported in the order aims, data-generating mechanisms,
estimands, methods and performance measures (ADEMP) \[20\]. The analysis
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
result. We also evaluated a sign-flip score test, which returns a
p-value only, and the interval obtained by inverting it. The clinical
applications and the decision rule of Section 2.6 build on these aims
but are not simulation aims.

### 2.2 Data-generating mechanisms

The main grid crossed 16 cluster counts, N ∈ {10, 15, 20, 25, 35, 45,
55, 65, 75, 85, 95, 105, 115, 125, 135, 145}, with 6 maximum numbers of
occasions, T ∈ {2, 3, 4, 6, 8, 10}. It also crossed five mean-model
scenarios (M1-M5) and the first three mechanisms below, giving 1,440
design points, and the same five scenarios with the fourth mechanism
(480 design points). Each design point was simulated under the
alternative (H1) and the null (H0), and the two hypotheses differed only
in the tested coefficient. T is a maximum. Each subject was observed at
all T occasions with probability 0.7, and otherwise at a number of
occasions drawn uniformly from $1$ to $T-1$. At T = 2, between 29% and
31% of subjects therefore contributed a single observation. Treatment
was assigned to subjects with probability 0.5.

The mean model of scenario M1 is
$\mathrm{logit}\,\mu = \beta_0 + \beta_1 t + \beta_2 \mathrm{trt} + \beta_3\, t \times \mathrm{trt}$.
Its marginal coefficients under H1 are (-2, 0.1, 0.2, 0.2), which give a
marginal prevalence of about 0.14 at T = 2. M2 adds a subject-level
covariate (age). M3 and M4 add a random slope on time; in GenBinaryY its
counterpart is a latent standard deviation that changes over time. M4
uses a smaller interaction (0.05). M5 adds a second covariate that
varies over time (dose), with its own random slope.

The four mechanisms differ in how the dependence is generated.

1. GenBinaryY is a marginalized model with a latent variable and a first-order transition term (R package MMLB) \[1,2,22\]. The latent standard deviation is a linear predictor, which in M3-M5 grows by 0.01 per occasion. This is the mechanism of the marginalized model.
2. GLMM-native is a logistic mixed model with normal random effects. The random intercept has variance 1 and, in M3-M5, the random slope on time has variance 0.001. The conditional coefficients are the marginal ones multiplied by 1.16. This is the mechanism of the mixed model.
3. GLMM-AR1 adds a latent first-order autoregressive process to the random effects of GLMM-native (ρ = 0.6, variance 0.5). Its random intercept is still constant within subject.
4. Copula-AR1 is a Gaussian copula with logistic margins. The latent normal scores of a subject have correlation $\rho^{|t_j - t_k|}$, with ρ = 0.6, and no component is constant within subject. None of the methods compared assumes this mechanism.

A sixth scenario, M1c, simulated with the first three mechanisms, uses the mean model of M1 with the prevalence (0.59) and the between-subject variance (3.5) of the clinical cohort of Section 2.7. It covers 72 design points, with N up to 45 and T up to 6. Four further scenarios, M6 and M7 on GenBinaryY and on GLMM-native, test the convergence check of Section 2.4 where the latent variability grows strongly over time. In GenBinaryY the latent standard deviation grows by 0.1 or 0.25 per occasion; in GLMM-native the random-slope variance is 0.01 or 0.04. They cover the full grid (384 design points).

The values of the main scenarios were compared with the clinical cohort of Section 2.7. Their marginal prevalence (0.14 at T = 2 and 0.28 at T = 10 in M1) is much lower than in the cohort (0.59), and their random-intercept variance (1) is lower than the patient variance of the cohort (3.5); scenario M1c was added for this reason. The cohort has no treatment and few patients with repeated timepoints, so the interaction per occasion has no direct anchor in it. The cumulative value of the simulated interaction at ten occasions (1.8 on the logit scale) is of the order of the regional contrasts estimated in the cohort, for example 1.1 between lateral and anterior zones.<!-- noclaim: justification of parameter values -->

Each design point received 1,000 replicates per hypothesis. Replicates in which one treatment arm was never observed after the first occasion were excluded from every denominator (247 per hypothesis in each of the three study mechanisms). Replicate i of scenario m used the seed 20260805 + 1000m + i under the L’Ecuyer-CMRG generator. The seed was the same for both hypotheses and for every N, T and mechanism, so the datasets are common random numbers across cells \[20\].

### 2.3 Estimands and true values

The estimand is the coefficient of the treatment-by-time interaction on
the logit scale. It is marginal for the marginalized model, the GEEs and
the shrinkage GLM, and conditional on the random effects for the mixed
model. In GenBinaryY the marginal true value is the generating
coefficient. There is no conditional true value, so bias and coverage of
the mixed model are not defined on that mechanism. In GLMM-native the
conditional true value is the generating coefficient. The marginal value
was obtained by Gauss-Hermite quadrature on the realized design; under
H0 it is not exactly zero, and it ranges from -0.0030 to -0.0011 across
scenarios and T. In GLMM-AR1 the marginal true value is the target, and
the conditional value was computed for each replicate. In Copula-AR1 the
marginal true value is the generating coefficient, and there is no
conditional true value.

### 2.4 Methods compared

The eight methods of Table 2 and the further strategies of Table 1 were applied to every dataset. The marginalized multilevel model was fitted with `MMLB::mm` (20 quadrature points), with the latent-variable and transition formulas of the GenBinaryY generator. GEEs were fitted with `geepack` \[23\] with exchangeable, independence and first-order autoregressive (AR(1)) working correlations. They were combined with the uncorrected sandwich variance and with the corrections of Mancl and DeRouen \[8\], Kauermann and Carroll \[9\] and Fay and Graubard \[10\] (bound b = 0.75, as in \[10\]), computed with the R package `fewclust` (Table S2), written for this study and released with the code, which implements the published formulas and agreed with a computation by hand (Section S4). The CR2 estimator of `clubSandwich` \[24,25\] was added with exchangeable and independence correlations. The mixed model was fitted with `lme4::glmer` \[26\] (Laplace approximation). After the main run, the mixed model with a random intercept was refitted by adaptive Gauss-Hermite quadrature with 10 nodes on a subset of null design points (Section 3.10). It had a random intercept in M1-M2 and the random slopes of the generating model in M3-M5. It was tested by a Wald test and by a likelihood-ratio test against the model without the interaction, with the same random effects. The shrinkage GLM rescales the mixed-model estimate and standard error by $(1 + k^2 D_{11})^{-1/2}$, with $D_{11}$ the random-intercept variance. Its Wald statistic therefore equals that of the mixed model on the same fit; its rejection rate can still differ, because the usability limits of Section 2.5 apply to the rescaled estimate and standard error. Software versions are in Table S2.

`MMLB::mm` writes the latent standard deviation as $\exp(a_0) + \exp(a_1)\,t$, so a slope of zero lies on the boundary of the parameter space. In the main scenarios with a latent standard deviation that changes over time (M3-M5), the generating slope is close to that boundary. Without a check, two failure modes of the optimizer (`nlm`) produced excess rejections with long series: non-convergence with a shifted estimate, and convergence with a collapsed standard error. We therefore added a convergence check. A fit was declared at the boundary if `nlm` did not report convergence (return code other than 1 or 2) or if an estimated slope $\exp(a_k)$ was below 0.001. Such a fit was replaced by the fit without the time trend of the latent standard deviation, and kept if that refit failed.<!-- noclaim: definition of the check -->

Wald tests were referred to the standard normal in the primary tables and, as an alternative, to a t distribution with N - p degrees of freedom, p being the number of mean-model coefficients. The likelihood-ratio test was referred to a $\chi^2$ distribution with one degree of freedom.

A sign-flip score test \[27\] was applied to every dataset. It assumes independent clusters; it is exact when the score contributions are symmetric about zero under the null, asymptotically valid otherwise, and allows clusters to differ in variance \[27\]. The test was declared not computable, and returned no answer, when the model without the interaction diverged (a coefficient of at least 1,000 in absolute value) or every score contribution was below $10^{-6}$ in absolute value, because the interaction is then not identified. The test gives a p-value only; its computation and the interval obtained by inverting it are described in Section S7.

### 2.5 Performance measures

A replicate was usable for a Wald-based method when the fit returned without error and the estimate, its standard error and its p-value were finite, with \|estimate\| \< 10 and $0 <$ SE \< 10. The probability of a usable result was reported for every method and cell. Its denominator, here and in Table 3, Table 4 and Table S7, is every replicate with an estimable design (Section 2.2), including those whose dataset could not be generated or whose fit failed. Rates under the null and the alternative were computed on each method’s usable replicates, so a method that fails on difficult datasets is assessed on easier ones; Table S6 repeats the type I error on the replicates usable for all methods, and Table S9 by stratum of estimability.<!-- noclaim: definitions -->

Type I error was the rejection rate at the nominal 5% under H0. Over the grid we report the pooled rate and two summaries computed within each mechanism. A design point is a combination of scenario, mechanism, number of clusters and number of occasions; a cell is a design point under one hypothesis. For each of the 16 cluster counts we took the absolute departure from 5% of the replicate-weighted rate. The worst case is the largest of these absolute departures, and the mean absolute departure is their mean. In the direct comparison of the two strategies of the rule (Section 3.4, Figure 1, Table S31) the absolute departure is instead averaged over design points. At each design point, a strategy was clearly above 7.5% when the lower limit of the exact (Clopper-Pearson) 95% interval of its rejection rate exceeded 7.5%, and clearly below 2.5% when the upper limit fell below 2.5%. These limits are Bradley’s liberal criterion, 0.5 to 1.5 times the nominal level \[28\]; for the coverage of 95% intervals they were applied to the rate of non-coverage.<!-- noclaim: definitions --> Unconditional validity is the probability of a usable result multiplied by the probability of not rejecting a true null given one. Power was reported raw and size-adjusted. Bias, empirical and model-based standard errors, root mean squared error and coverage of the 95% interval were computed against the estimand of each method.

Monte Carlo standard errors (MCSE) were computed for every measure \[20\]. With 1,000 usable replicates, the MCSE of a rejection rate at 5% is 0.69 percentage points. This number of replicates was chosen so that a single cell resolves a departure from 5% of about 1.4 points. The datasets are common random numbers across cells. The MCSE of every summary over cells, and of every difference between methods, was therefore obtained by a bootstrap that resamples replicate indices jointly across all cells (2,000 resamples). A comparison between methods is written as a difference only when its bootstrap interval, the estimate ± 2 MCSE, lies entirely beyond a relevance margin fixed in advance. The margin is one percentage point for differences in type I error or coverage averaged over cells, 2.5 points for a single cell, two points for differences in power and 2.5 points for differences in usable shares. Intervals entirely within the margin are reported as similar. Intervals that straddle it are reported as differences of about the size of the margin.<!-- noclaim: definitions -->

### 2.6 Decision rule

The recommendation of Figure 1 comes from a rule applied to every strategy that returns an estimate and an interval, a strategy being a method with its reference distribution. A strategy was excluded if its type I error was clearly above 7.5% at any design point. With 1,000 usable replicates a design point is flagged only from an observed rate of 9.3%, and with half as many from 10%: a strategy never clearly above 7.5% is one in which no excess was detected, not one shown to be calibrated. Among the remaining strategies, the one that returned a usable result most often was recommended.<!-- noclaim: definition of the rule --> The rule was applied to all 2,376 design points of the four mechanisms. The likelihood-ratio test and the sign-flip test give no interval and are reported for comparison.<!-- noclaim: definition of the rule -->

### 2.7 Clinical applications

The cohort comprises 207 neonates with ten lung-ultrasound zones scored at one to three timepoints. It gives 230 patient-timepoints and 2,296 zone-level outcomes; 22 patients contribute more than one timepoint.<!-- noclaim: cohort description --> The outcome is severe loss of aeration (prevalence 0.59). In a logistic mixed model of the outcome on zone region with a patient random intercept, the patient variance is 3.5. Adding a patient-timepoint level drives the patient variance to the boundary and places a variance of 3.8 at the new level (likelihood-ratio $\chi^2$ = 47 on one degree of freedom). The contrast of interest is posterior against anterior zones, adjusted for disease, invasive ventilation and gestational age. It was estimated by the marginalized model, the GEEs and the mixed model. With 207 clusters the cohort lies above the simulated range, so it was used as a plasmode. We drew 500 subsamples of whole patients at each of 10, 15, 20, 25 and 50 patients, with the full-cohort estimate of each method as its reference. In both applications the ratio of the mean model-based standard error to the standard deviation of the estimates over subsamples was multiplied by $\sqrt{1 - N/P}$, with P the number of patients in the full dataset, because subsamples drawn without replacement from a finite cohort vary less than independent samples.<!-- noclaim: definition -->

The same analysis was repeated on the respiratory trial distributed with `geepack` \[23\] (111 patients, four visits, active treatment against placebo). The contrast of interest was the treatment effect, adjusted for baseline status, visit, age, sex and centre.

## 3. Results

### 3.1 Type I error on the main grid

Table 2 gives the type I error of the eight methods on the 1,440 design
points of the main grid with the first three mechanisms, under a normal
reference, and Figure 2 shows the type I error of the marginalized model
at every design point. The pooled rejection rate of the marginalized
model with the convergence check ranged from 5.2% to 5.5% across
mechanisms. It ranged from 5.3% to 5.8% for Fay-Graubard and from 5.5%
to 6.1% for Mancl-DeRouen. On GenBinaryY it reached 6.9% for
Kauermann-Carroll and 8.3% for the uncorrected
GEE.<!-- noclaim: descriptive values per method -->

The marginalized model had the smallest worst-case departure from 5% on
all three mechanisms: 0.83, 0.44 and 0.56 percentage points on
GenBinaryY, GLMM-native and GLMM-AR1, respectively. Its worst-case
departure was smaller than that of the next method by 2.8 points on
GenBinaryY and by about one point on the two random-effects
mechanisms.<!-- claim: cl:C01|GB as=larger; cl:C01|NAT as=about1; cl:C01|AR1 as=about1 -->
Its mean absolute departure (0.24 to 0.54 points) and that of
Fay-Graubard (0.52 to 1.0) were
similar.<!-- claim: cl:C02|GB as=similar; cl:C02|NAT as=similar; cl:C02|AR1 as=similar -->
On matched replicates usable for both methods, the two rejected a true
null at similar
rates.<!-- claim: cl:C03|GB as=similar; cl:C03|NAT as=similar; cl:C03|AR1 as=similar -->

The uncorrected GEE was above 5% by 19 points at N = 10 on GenBinaryY.<!-- claim: cl:C12a|GB as=descr --> Above 45 clusters, Fay-Graubard and the marginalized model were both within about half a point of 5% on the random-effects mechanisms.<!-- claim: cl:C25a|NAT as=descr; cl:C25b|NAT as=similar; cl:C25b|AR1 as=similar -->

The mean absolute departure of the marginalized model fell from 1.5 points at two occasions to 0.71 at four or more (Table S4). At two occasions the five marginal methods departed by 1.5 to 2.5 points on average. The conditional methods departed by 5.6 (mixed model, Wald) to 11 points (shrinkage GLM).<!-- claim: cl:C13d|all as=larger -->

### 3.2 The convergence check of the marginalized model

In scenarios M3-M5 the check replaced 68% of the fits of the marginalized model with the fit without the time trend of the latent standard deviation. Without the check, the marginalized model was clearly above 7.5% at 149 of the 864 null design points of these scenarios, 90 of them with eight or ten occasions. With the check it was clearly above 7.5% at 0.<!-- claim: reg:cc_hi_study as=descr; reg:cc_hi_cc as=descr --> Its mean absolute departure at ten occasions fell from 2.8 to 0.67 points.<!-- noclaim: descriptive values -->

In scenarios M6-M7 the fit without the time trend is misspecified. There the check intervened in 36% to 51% of fits under H0, and in 58% to 60% at two or three occasions, where `nlm` most often failed to converge. With a t reference, the marginalized model with the check was never clearly above 7.5% at the 384 design points. Without the check, and with a t reference, it was clearly above 7.5% at 9 design points.<!-- claim: reg:ns_mmmt_hi as=descr; reg:ns_nocc_t_hi as=descr --> Bias, standard errors and coverage of the marginalized model with the check in these scenarios are in Section S13.

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
marginalized model (68%) and the three corrected sandwiches (69% to 70%)
were similar.<!-- claim: cl:C14a|all as=larger; cl:C14d|all as=descr -->
Unconditional validity ranged from 61% for the Wald test of the mixed
model to 80% for the likelihood-ratio
test.<!-- claim: cl:C14b|all as=larger -->

### 3.4 Reference distribution and the decision rule

Referring the Wald statistic to a t distribution with N - p degrees of freedom lowered every pooled rejection rate. The shift was 0.53 to 1.1 percentage points over the grid and 1.4 to 3.7 points at 25 clusters or fewer (Table S11).<!-- claim: cl:C17b|all as=larger --> Under the t reference the marginalized model became conservative at ten clusters. There it was clearly below 2.5% at 16 of 102 design points, against 4 for Fay-Graubard.<!-- claim: reg:rs_mmmt_N10_lo as=descr; reg:rs_fgt_N10_lo as=descr -->

Table 4 applies the rule of Section 2.6 to the 2,376 design points of
the four mechanisms. Among the strategies with an interval, only two
were never clearly above 7.5%: GEE with AR(1) working correlation,
Fay-Graubard and a t reference; and the marginalized model with the
check and a t reference.<!-- claim: reg:t4 as=descr --> The rule
selected the GEE strategy, which returned a usable result in 92% of
datasets, against 91% for the marginalized model, a difference within
the margin for usable shares (Section
2.5).<!-- claim: reg:t4 as=descr; reg:hb as=descr --> Never clearly
above 7.5% means that no excess was detected. The highest rate at a
design point was 9% for the GEE strategy and 8.4% for the marginalized
model, and an excess could not be excluded, the upper limit of the exact
interval lying above 7.5%, at 218 and 438 design points (Table
S32).<!-- noclaim: descriptive values --> The GEE strategy was clearly
below 2.5% at 116 design points and the marginalized model at 24, 18 of
them at ten clusters.<!-- claim: reg:t4 as=descr -->

With exchangeable or independence working correlation, Fay-Graubard with
a t reference was clearly above 7.5% at 11 and 2 design points, none of
them in the main grid or M1c. Mancl-DeRouen with a t reference was
clearly above 7.5% at 16 to 57 design points across working
correlations, Kauermann-Carroll at 132 to 199, and the Wald test of the
mixed model at 515, 55 of them on data generated by its own model
(GLMM-native). The likelihood-ratio test was clearly above 7.5% at 412
design points.<!-- claim: reg:t4 as=descr -->

With two or three occasions and at most 25 clusters, the GEE strategy
and the marginalized model returned a usable result in fewer than half
of the datasets (46% and 39%).<!-- claim: reg:t4 as=descr --> With N - 2
degrees of freedom \[16\], on the 1,896 design points of the three study
mechanisms, the rule selected the same strategy; the marginalized model
was again never clearly above 7.5%, and clearly below 2.5% at 8 instead
of 23 design points (Table S26).<!-- claim: reg:t2 as=descr -->

Compared directly, with the block bootstrap and the margins of Section 2.5 (Figure 1, Table S31), the two strategies had similar mean absolute departures from 5% over all design points (both 1.1 percentage points), similar usable shares (difference 1.2 points) and similar power (difference 1.4 points).<!-- claim: reg:hb as=descr --> The GEE strategy was closer to 5% at ten clusters with three or more occasions (1.6 against 2.9 points), where the marginalized model was conservative (mean type I error 2.1% against 4.0%).<!-- claim: reg:hb as=descr --> With two occasions and 15 or more clusters the GEE strategy was conservative (2.5%) and the marginalized model was closer to 5% by about the margin (difference 1.0 points); with 15 or more clusters and three or more occasions the two were similar. With two occasions and 15 or more clusters the GEE strategy returned a usable result more often, by 4.2 points.<!-- claim: reg:hb as=descr -->

### 3.5 Working correlation and implementation of the corrections

Changing the working correlation of Fay-Graubard from exchangeable to independence or AR(1) changed its type I error by -0.29 and -0.31 percentage points on average over the 1,896 design points of the three study mechanisms, within the one-point margin.<!-- claim: cv:wc_ind as=similar; cv:wc_ar1 as=similar --> The two implementations of the corrections are compared in Section S4.

### 3.6 The fourth mechanism

On Copula-AR1 (480 design points) neither the GEE strategy nor the marginalized model with the check was clearly above 7.5% at any design point (Table 4). Their mean type I errors with a t reference were 4.6% and 5.1%, respectively.<!-- claim: reg:t4 as=descr --> For both, the difference from the three study mechanisms was within the one-point margin.<!-- claim: cv:cop_gee as=similar; cv:cop_mmm as=similar --> The marginalized model estimated the interaction with a mean bias of 0.027 empirical standard errors and covered the true value in 94% of replicates under a normal reference, although its model of the dependence is misspecified on this mechanism (Table S28).<!-- noclaim: descriptive --> The Wald test of the mixed model (Laplace approximation) with a t reference rejected a true null in 10% of usable replicates on average, 3.2 points more than on the three study mechanisms, and was clearly above 7.5% at 223 of the 480 design points.<!-- claim: cv:cop_glmm as=larger --> With a random intercept only, and N ≥ 85, its model-based standard error was 0.64 to 0.81 times the standard deviation of its estimates, while the estimates were centred.<!-- noclaim: mechanism --> GEE with exchangeable working correlation and Fay-Graubard was clearly above 7.5% at 7 design points, all with T = 4 and N = 10 to 25.

### 3.7 Sign-flip score test

Over the 2,376 design points the sign-flip test rejected a true null in 5.4% of datasets on average, and it was clearly above 7.5% or clearly below 2.5% at no design point.<!-- claim: reg:flip as=descr --> It answered in 97% of datasets with an estimable design and, with two or three occasions and at most 25 clusters, in 72%, against 46% for the GEE strategy (Table 4).<!-- claim: reg:t4 as=descr --> Its type I error in those small designs was 5.2%, and its power per attempted dataset over all design points was 21%.<!-- noclaim: descriptive -->

### 3.8 Power and coverage

Raw and size-adjusted power of the eight methods are in Table S10.

Coverage of the 95% interval was 94% to 95% for Fay-Graubard,
Mancl-DeRouen and the marginalized model on every mechanism, and 93% to
94% for Kauermann-Carroll (Table
S12).<!-- claim: cl:C20a|GB as=similar; cl:C20a|NAT as=similar; cl:C20a|AR1 as=similar -->

### 3.9 Scenario anchored to the cohort

Scenario M1c has the prevalence and between-subject variance of the cohort (72 design points). Under a normal reference the marginalized model rejected a true null in 6.8% to 7.0% of usable replicates, against 5.0% to 5.1% for Fay-Graubard. Its worst-case departure was larger than that of the best method by about the one-point margin.<!-- claim: cm:C03|GB as=larger; cm:C03|NAT as=larger; cm:C03|AR1 as=larger; cm:C01|GB as=about1; cm:C01|NAT as=about1; cm:C01|AR1 as=about1 --> Under a t reference its rate was 5.0% to 5.1%, and that of Fay-Graubard 3.5% to 3.6%.<!-- noclaim: descriptive values -->

### 3.10 Adaptive quadrature for the mixed model

The refit covered 123 null design points of scenarios M1, M2 and M1c (Section S12, Table S34). With a t reference, quadrature lowered the rejection rate under the null by 3.1 to 4.0 percentage points on average, depending on the mechanism. The change was largest at two occasions with 145 clusters, from 10% to 20% under the Laplace approximation to 2.6% to 4.5%. The Wald test was clearly above 7.5% at 33 of these design points under the Laplace approximation and at 8 with quadrature, all with T = 10: 6 on Copula-AR1 (12% to 16%) and 2 on GenBinaryY (9.3% and 9.6%). On data generated by its own model it was no longer clearly above 7.5% at any of them. Quadrature also made the test conservative, mostly with few clusters: it was clearly below 2.5% at 37 design points, against 3 under the Laplace approximation (Table S34).<!-- noclaim: descriptive values of the sensitivity analysis -->

### 3.11 Clinical applications

In the full cohort the posterior-versus-anterior contrast was 0.38 (SE 0.17) by the marginalized model and 0.38 by GEE. The mixed model gave 0.48 (SE 0.12) on its conditional scale. In subsamples of 10 patients, 22% of draws contained no invasively ventilated patient, so that the mean model could not be fitted (Table S15). Over all subsamples of this size, including those, a usable result was obtained by the mixed model in 78%, by the marginalized model in 77%, by the uncorrected GEE in 70%, by Fay-Graubard in 34% and by Mancl-DeRouen in 15%.<!-- noclaim: descriptive values per method --> Most failures of the corrections returned a point estimate without a standard error (Table S15). Across subsample sizes, after the finite-population correction, the standard error of the mixed model was 0.51 to 0.58 times the empirical standard deviation of its estimates.<!-- noclaim: descriptive --> The corresponding ratio was 0.84 to 0.99 for the marginalized model. The bias of the marginalized model against its full-cohort estimate fell from 0.16 at 10 patients to a value within two Monte Carlo standard errors of zero at 50.

In the respiratory trial the treatment effect was 1.3 (SE 0.35) by the marginalized model and 1.3 by GEE, and 2.1 on the conditional scale of the mixed model. In subsamples of 10 patients, with the full mean model, a usable result was returned by the mixed model in 53% of draws, by the marginalized model in 51%, by the uncorrected GEE in 47%, by Fay-Graubard in 23% and by Mancl-DeRouen in 23%; as in the cohort, the corrections mostly failed with an estimate but no standard error (Table S30).<!-- noclaim: descriptive values per method -->

## 4. Discussion

### 4.1 The recommendation

Two strategies met the rule of Section 2.6: GEE with an AR(1) working
correlation, the Fay-Graubard correction and a t reference with N - p
degrees of freedom; and the marginalized model with the convergence
check and a t reference (Figure 1). Across the 2,376 design points of
the four mechanisms neither was ever clearly above 7.5%, and they
returned a usable result in 92% and 91% of
datasets.<!-- claim: reg:t4 as=descr --> The rule selects the GEE
strategy by usable share, but the difference is within the margin. In
the direct comparison, which was not pre-specified, the two had similar
level, usable share and power over the design as a
whole.<!-- claim: reg:hb as=descr --> We place the GEE strategy first
for two reasons: where the two differed beyond the margin, at ten
clusters, it was the one closer to nominal, and it needs no convergence
check. With two occasions and 15 or more clusters the two are
comparable; there the GEE strategy is conservative and returns a usable
result more often.<!-- claim: reg:hb as=descr --> Because the fourth
mechanism has AR(1) dependence on the latent scale, it could favour the
AR(1) working correlation. The working correlation changed the level of
Fay-Graubard by less than half a point on average, and the rule selected
the same strategy on the three study mechanisms alone (Table
S26).<!-- claim: cv:wc_ar1 as=similar; reg:t2 as=descr -->

The rule asks that a recommendation hold at every design point, because
an analyst knows the number of clusters and occasions but not the
mechanism that generated the data. It replaces a criterion that
maximized unconditional validity at each cluster count. That criterion
rewards a method for answering, whatever its level. On the present data
it would have recommended the likelihood-ratio test of the mixed model
up to 25 clusters, a test that was clearly above 7.5% at 412 design
points.<!-- claim: reg:switch_GB as=descr; reg:t4 as=descr -->
Fay-Graubard with an exchangeable working correlation would not have met
the rule: with a t reference it was clearly above 7.5% at 11 design
points, 4 in scenarios M6-M7 and 7 on the fourth mechanism (Table
4).<!-- claim: reg:t4 as=descr -->

With two or three occasions and up to 25 clusters, no strategy with an
interval was both available and calibrated. The recommended strategy and
the marginalized model returned a usable result in fewer than half of
the datasets and detected the simulated effect in fewer than 2% of
attempted
datasets.<!-- claim: reg:t4 as=descr; reg:upow_small_ceiling as=descr -->
The sign-flip test answered in 72% of those datasets and kept its level,
but it detected the simulated effect in only 4.2% of attempted datasets:
it gives a valid test where the estimators return nothing, not a
powerful
one.<!-- claim: reg:flip as=descr; reg:flip_upow_small as=descr -->
Inverting it gives an interval with nominal coverage (95% under the
alternative, against 96% for the GEE strategy, and its non-coverage
never clearly outside Bradley’s limits at any design point), but in
those designs both of its ends were finite in only 14% of datasets, and
over all design points, where both intervals existed, it was 17% wider
than the GEE interval on average (Section S7). We therefore do not
propose it as an interval to
use.<!-- claim: cv:ci_cov as=larger; reg:ci as=descr --> With a random
slope and two occasions, the random-effects covariance of the mixed
model is not identifiable, since two subject-level parameters saturate
two observations, and its Wald test was essentially never usable. The
case for the maximal random-effects structure \[29\] is thus bounded by
power \[30\] and, with two occasions, by identifiability. In these
designs, what can be learned depends more on the number of occasions and
clusters than on the analysis. Four or more occasions per subject, or
more clusters, should be planned where
possible.<!-- noclaim: interpretation and design advice -->

### 4.2 The convergence check

The check addresses a property of the parametrization rather than of the marginalized model. Because the slope is exponentiated, a slope near zero lies on the boundary, where `nlm` fails in the two ways described in Section 2.4. Replacing such fits with the fit without the time trend removed the excess.<!-- noclaim: interpretation of the mechanism --> The check intervened in a large share of fits (Section 3.2), mainly because `nlm` failed to converge with two or three occasions.<!-- claim: reg:pb1_n6 as=descr --> Where it intervened, the type I error stayed close to nominal. It also remained calibrated in M6-M7, where the fit without the time trend is misspecified. Analysts who use the marginalized model with a latent standard deviation that changes over time should apply the check and report how often it intervened.<!-- noclaim: recommendation -->

### 4.3 Relation to previous work

We reproduce the anticonservatism of the uncorrected sandwich with few
clusters \[8,9\]. Li and Redden found that Kauermann-Carroll attains the
nominal level at ten clusters under a t reference \[16\]. In our data,
with the correction computed from its published formula, its rate at ten
clusters under a t reference was 9.1% to 12%, above Bradley’s limit, and
over the four mechanisms it was clearly above 7.5% at 132 to 199 of
2,376 design points across working
correlations.<!-- claim: cv:kc10_NAT as=descr; cv:kc10_GB as=descr; cv:kc10_AR1 as=descr; reg:t4 as=descr -->
Their clusters were groups of participants \[16\], whereas ours are
subjects measured on several occasions, so the difference may reflect
the structure of the clusters rather than the estimator. Our
Kauermann-Carroll values come from `fewclust`, which agreed with a
computation by hand (Section S4); `geesmv` computes a different
quantity, so published values obtained with it are not directly
comparable with ours.<!-- noclaim: interpretation --> Mancl and DeRouen
reported sizes close to nominal with ten subjects \[8\], and Bie et
al. favoured GEE with the Mancl-DeRouen correction \[18\]. Under a t
reference its pooled rate here was 4.7% to 5.3%, inside Bradley’s
limits. With exchangeable working correlation it was nevertheless
clearly above 7.5% at 57 design points, 17 of them in scenarios M6-M7
and 26 on the fourth
mechanism.<!-- claim: cl:C24c|NAT as=descr; cl:C24c|GB as=descr; reg:t4 as=descr -->
The marginalized model has been studied for moderate to long series
\[22\], and the short series studied here complement that work.
Estimation methods for mixed logistic models with few clusters are
reviewed by McNeish \[31\].

### 4.4 The clinical applications and the fourth mechanism

In the neonatal cohort the standard error of the mixed model fell short
of the spread of its estimates at every subsample size, and that of the
marginalized model was closer to it (Section 3.11). In the respiratory
trial both fell short with few patients, the marginalized model more
than the mixed model up to 15 patients (Table S30). With seven
mean-model coefficients and few patients its estimates are unstable, and
their spread is inflated by large estimates that the usability filter
keeps. On the fourth mechanism with ten occasions, the excess rejection
of the mixed model persisted under adaptive quadrature and comes from
its standard error, which ignores the stronger correlation between
neighbouring occasions. Elsewhere much of its excess came from the
Laplace approximation, which is inaccurate with few binary observations
per cluster
\[31\].<!-- noclaim: interpretation of Sections 3.6 and 3.10 -->

### 4.5 Limitations

In the main scenarios the random-slope variance and the time trend of
the latent standard deviation were small. Conclusions about random
slopes therefore rest mainly on scenarios M6-M7, which cover one
structure (a slope on time) and two mechanisms. Trends of the latent
standard deviation between those of the main scenarios (0.01 per
occasion) and of M6 (0.1) were not simulated; reporting how often the
check intervenes (Section 4.2) is the safeguard. The number of occasions
observed did not depend on the outcomes, so informative dropout was not
studied, and the results concern a linear treatment-by-time interaction.
Three working correlations were compared; a misspecified variance
function was not simulated. Three of the four mechanisms contain a
component that is constant within subject, and two of them are the home
models of the marginalized model and of the mixed model. The fourth has
only serial dependence, AR(1) on the latent scale, like the working
correlation of the GEE strategy; the rule selected the same strategy
without it. In GLMM-native the marginal value under the null is not
exactly zero. The excess rejection it can cause is at most 0.029
percentage points in the main scenarios and 0.26 in M6-M7 for the
marginalized model and the GEEs, and subtracting it changes none of the
flags clearly above 7.5% or clearly below 2.5% of these
methods.<!-- claim: reg:null_flips as=descr --> The convergence check
was introduced after part of the results had been seen, and was then
tested on scenarios simulated afterwards. The decision rule itself
replaced the criterion of an earlier version after that version’s
results had been seen, and the direct comparison of the two strategies
that met the rule was not pre-specified. The tables report the mixed
model fitted by the Laplace approximation, the default of `lme4`;
adaptive quadrature, which `lme4` offers only with a single random
intercept, was examined after the main run on a subset of design points
(Section 3.10) and could not be examined with random slopes \[31\].
Penalized GEE \[15\] and the degrees of freedom proposed with the
Fay-Graubard correction \[10\] were not evaluated. Treatment was
assigned by simple randomization with probability one half, so with ten
subjects some datasets had unbalanced arms, and datasets with a
degenerate design were excluded. Mancl-DeRouen and Fay-Graubard were
computed with two implementations that agreed; the Kauermann-Carroll
results rest on one. The sign-flip test gives no estimate, and its
interval is rarely bounded where it would be most needed. The clinical
applications use two datasets with one outcome each. The subsampling
resamples whole patients, and in the neonatal cohort the
patient-timepoint level that carries most of the dependence is
identified by only 22 patients measured more than
once.<!-- noclaim: limitations -->

## 5. Conclusions

With few clusters and short series, the choice of analysis for a
clustered binary outcome decides whether the test keeps its nominal
level and whether it returns an answer at all. GEE with an AR(1) working
correlation, the Fay-Graubard correction and a t reference, and the
marginalized multilevel model fitted with a convergence check and a t
reference, were never clearly above 7.5% at any simulated design point
of four mechanisms, and their calibration differed beyond the relevance
margin only at ten clusters. The Wald test of the mixed model fitted by
the Laplace approximation was anticonservative; with adaptive
quadrature, examined for a random intercept, the excess remained only
with ten occasions, mainly under serial dependence, and with few
clusters the test became conservative. With two or three occasions and
up to 25 clusters, neither strategy returned a usable result in half of
the datasets. A sign-flip score test can still test for an effect there,
with little power, and no analysis estimates the effect reliably; in
that regime the design limits what any analysis can
show.<!-- claim: reg:t4 as=descr; reg:hb as=descr; cv:cop_glmm as=larger -->

In the designs studied we recommend GEE with an AR(1) working
correlation, the Fay-Graubard correction and a t reference with N - p
degrees of freedom, which needs no convergence check and was closer to
nominal at ten clusters. The marginalized model with the convergence
check and a t reference is a competitive alternative, slightly closer to
nominal with two occasions and 15 or more
clusters.<!-- noclaim: recommendation -->

## Code and data availability

The simulation code, the per-replicate results and the scripts that
produce every number, table and figure of this paper are available at
https://github.com/lucavd/MMM (commit ab90156). They will be archived on
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
24. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
    variance estimation and hypothesis testing in fixed effects models.
    *Journal of Business & Economic Statistics* 2018;36(4):672-683.
    doi:10.1080/07350015.2016.1247004
25. Pustejovsky JE, Tipton E. Corrigendum: small sample methods for
    cluster-robust variance estimation and hypothesis testing in fixed
    effects models. *Journal of Business & Economic Statistics*
    2023;41(2):650-652. doi:10.1080/07350015.2023.2174123
26. Bates D, Mächler M, Bolker B, Walker S. Fitting linear mixed-effects
    models using lme4. *Journal of Statistical Software*
    2015;67(1):1-48. doi:10.18637/jss.v067.i01
27. Hemerik J, Goeman JJ, Finos L. Robust testing in generalized linear
    models by sign flipping score contributions. *Journal of the Royal
    Statistical Society Series B* 2020;82(3):841-864.
    doi:10.1111/rssb.12369
28. Bradley JV. Robustness? *British Journal of Mathematical and
    Statistical Psychology* 1978;31(2):144-152.
    doi:10.1111/j.2044-8317.1978.tb00581.x
29. Barr DJ, Levy R, Scheepers C, Tily HJ. Random effects structure for
    confirmatory hypothesis testing: keep it maximal. *Journal of Memory
    and Language* 2013;68(3):255-278. doi:10.1016/j.jml.2012.11.001
30. Matuschek H, Kliegl R, Vasishth S, Baayen H, Bates D. Balancing type
    I error and power in linear mixed models. *Journal of Memory and
    Language* 2017;94:305-315. doi:10.1016/j.jml.2017.01.001
31. McNeish D. Estimation methods for mixed logistic models with few
    clusters. *Multivariate Behavioral Research* 2016;51(6):790-804.
    doi:10.1080/00273171.2016.1236237

## Tables and figures

**Table 1.** Summary of the simulation study in the ADEMP structure
\[20\].

|  |  |
|----|----|
| Aims | Type I error and power of the test of no treatment-by-time interaction, bias and coverage of its estimate, and probability of a usable result, with few clusters and short series |
| Data-generating mechanisms | GenBinaryY (marginalized latent-variable transition model), GLMM-native (logistic mixed model), GLMM-AR1 (mixed model with latent first-order autoregression), Copula-AR1 (Gaussian copula with logistic margins and serial correlation only); scenarios M1-M5 on the full grid of 16 cluster counts (10 to 145) and 6 series lengths (2 to 10), 1,440 design points with the first three mechanisms and 480 with the fourth; M1c anchored to the clinical cohort (72 design points); M6-M7 with a steep time trend of the latent variability (384 design points); null and alternative arm for each, 1,000 replicates per arm |
| Estimands | Interaction coefficient on the logit scale: marginal (marginalized model, GEEs, shrinkage GLM) or conditional on the random effects (GLMM) |
| Methods | Marginalized multilevel model (`MMLB::mm`) with convergence check; GEE with exchangeable, independence and AR(1) working correlation, with uncorrected sandwich and the Mancl-DeRouen, Kauermann-Carroll and Fay-Graubard corrections (`geepack`, `fewclust`, `geesmv`) and CR2 (`clubSandwich`); GLMM with Wald and likelihood-ratio tests (`lme4`); shrinkage GLM (rescaled GLMM); sign-flip score test. Wald tests with normal and t (N - p degrees of freedom) reference |
| Performance measures | Probability of a usable result; type I error (pooled, worst case and mean absolute departure over cluster counts, clearly above 7.5% or below 2.5% per design point); raw and size-adjusted power; bias, empirical and model standard error, coverage; Monte Carlo standard errors by block bootstrap over replicates |

**Table 2.** Type I error at the nominal 5% on the main grid (1,440
design points), normal reference, first three mechanisms;
Kauermann-Carroll computed with `fewclust`. Rate: pooled rejection rate
(%). Worst: largest absolute departure from 5% over the 16 cluster
counts (percentage points). MAD: mean absolute departure over the
cluster counts (percentage points). GB = GenBinaryY, NAT = GLMM-native,
AR1 = GLMM-AR1. Monte Carlo standard errors (block bootstrap) at most
0.052 percentage points for rates and MAD and at most 0.30 for the worst
case; values of the marginalized model without the check in Table S3.

| Method | Rate GB | Rate NAT | Rate AR1 | Worst GB | Worst NAT | Worst AR1 | MAD GB | MAD NAT | MAD AR1 |
|----|----|----|----|----|----|----|----|----|----|
| GEE, Fay-Graubard | 5.8 | 5.3 | 5.3 | 4.7 | 2.8 | 2.9 | 1.0 | 0.5 | 0.6 |
| GEE, Mancl-DeRouen | 6.1 | 5.5 | 5.6 | 8.6 | 5.9 | 6.0 | 1.5 | 0.9 | 1.0 |
| Marginalized model, with check | 5.5 | 5.2 | 5.4 | 0.8 | 0.4 | 0.6 | 0.5 | 0.2 | 0.4 |
| GLMM, Wald | 8.6 | 6.3 | 6.6 | 4.5 | 1.8 | 2.2 | 3.5 | 1.4 | 1.6 |
| Shrinkage GLM | 9.3 | 6.6 | 6.9 | 4.9 | 2.3 | 2.4 | 4.2 | 1.6 | 1.9 |
| GEE, uncorrected | 8.3 | 7.6 | 7.7 | 21.0 | 17.4 | 18.3 | 3.8 | 3.0 | 3.1 |
| GLMM, likelihood ratio | 7.8 | 5.8 | 6.1 | 3.6 | 1.8 | 1.8 | 2.8 | 0.8 | 1.1 |
| GEE, Kauermann-Carroll | 6.9 | 6.3 | 6.4 | 13.3 | 10.3 | 10.7 | 2.4 | 1.7 | 1.8 |

**Table 3.** Replicates returning a usable result (%), by structure of
the random effects and number of occasions T, pooled over mechanisms and
cluster counts, every replicate with an estimable design in the
denominator; Kauermann-Carroll computed with `fewclust`. LR =
likelihood-ratio test with the random-effects structure of the Wald
test. Monte Carlo standard errors at most 0.13 points. T = 8 is omitted
for space; all methods and all T are in Table S7.

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

**Table 4.** The decision rule applied to the 2,376 design points of the
four mechanisms: main grid and M1c (1,512), scenarios M6-M7 (384) and
Copula-AR1 (480). Wald tests with a t reference (N - p degrees of
freedom); likelihood-ratio test with a $\chi^2$ reference. Clearly above
7.5% and clearly below 2.5%: lower limit of the exact 95% interval of
the rejection rate under the null above 7.5%, upper limit below 2.5%;
design points with fewer than 20 usable replicates are not assessed.
Usable: replicates returning a usable result (%), averaged over design
points, denominator the replicates with an estimable design; power:
rejections of a false null per attempted dataset (%). The column at ten
clusters covers 156 design points; the last column covers two or three
occasions and at most 25 clusters (216 design points). The rule
considers only strategies with an interval.

| Strategy | Above 7.5%, main and M1c | Above 7.5%, M6-M7 | Above 7.5%, Copula-AR1 | Below 2.5%, all | Below 2.5%, ten clusters | Usable (%) | Power (%) | Usable, T ≤ 3 and N ≤ 25 (%) |
|----|----|----|----|----|----|----|----|----|
| GEE AR(1), Fay-Graubard | 0 | 0 | 0 | 116 | 12 | 92.1 | 19.4 | 45.9 |
| Marginalized model with check | 0 | 0 | 0 | 24 | 18 | 90.9 | 20.8 | 39.2 |
| GEE independence, Fay-Graubard | 0 | 1 | 1 | 120 | 16 | 92.3 | 19.3 | 46.4 |
| GEE exchangeable, Fay-Graubard | 0 | 4 | 7 | 114 | 10 | 92.2 | 19.9 | 46.1 |
| GEE AR(1), Mancl-DeRouen | 4 | 6 | 6 | 116 | 11 | 92.0 | 19.6 | 45.3 |
| GEE exchangeable, Mancl-DeRouen | 14 | 17 | 26 | 114 | 10 | 92.1 | 20.0 | 45.5 |
| GEE AR(1), Kauermann-Carroll | 58 | 37 | 37 | 82 | 4 | 91.9 | 20.7 | 44.4 |
| GEE exchangeable, Kauermann-Carroll | 97 | 44 | 58 | 82 | 4 | 92.0 | 21.1 | 44.6 |
| GEE exchangeable, CR2 | 133 | 52 | 69 | 10 | 0 | 93.5 | 21.4 | 54.3 |
| GEE exchangeable, uncorrected | 232 | 68 | 97 | 1 | 0 | 94.0 | 22.7 | 58.3 |
| Shrinkage GLM | 271 | 93 | 244 | 6 | 5 | 81.7 | 23.7 | 33.7 |
| GLMM, Wald | 206 | 86 | 223 | 15 | 12 | 80.4 | 23.2 | 28.2 |
| GLMM, likelihood ratio (no interval) | 223 | 0 | 189 | 13 | 0 | 86.4 | 24.5 | 61.9 |
| Sign-flip score test (no interval) | 0 | 0 | 0 | 0 | 0 | 96.5 | 20.9 | 72.3 |

<img src="{{artifact:art_3c56d66b-d70e-4dd7-b4c2-b4c56cb89565}}" style="width:100.0%" />

**Figure 1.** Recommended analysis for the test of the treatment-by-time
interaction with few clusters. (a) Comparison of the two strategies of
the rule in each region of clusters and occasions, by the difference in
mean absolute departure from the nominal 5% (block bootstrap; Table
S31): the strategy is named where it is closer by more than the
one-point margin of Section 2.5 over the whole interval, “GEE = MMM”
elsewhere, with a dotted frame where the interval crosses the margin;
“usable”: datasets returning a usable result (%), GEE/MMM. GEE = AR(1)
working correlation, Fay-Graubard correction and a t reference with N -
p degrees of freedom, the recommended strategy; MMM = marginalized
multilevel model with the convergence check and a t reference, the
alternative; N = number of clusters, p = number of mean-model
coefficients. Dashed frame: two or three occasions and at most 25
clusters, where on average both answer in fewer than half the datasets
(Table 4) and the sign-flip test is a complement. (b) Mean absolute
departure from the nominal 5% (percentage points) by number of
occasions, over the 2,376 design points of the four mechanisms
(likelihood-ratio test: the 2,094 design points at which it could be
assessed), logarithmic axis. (c) Percentage of design points at which
each strategy is clearly above 7.5% (exact 95% interval), by number of
clusters (Table 4). (d) Mean absolute departure from 5% of the two
strategies of the rule by number of clusters, with two occasions
(dashed) and three or more (solid). (e) Datasets returning a usable
result (%) and (f) datasets in which the simulated effect was detected,
among those attempted (%), by number of occasions, for the two
strategies and the sign-flip test. Wald tests with a t reference,
likelihood-ratio test with a $\chi^2$ reference. Bars in (b) and (d):
estimate ± 2 Monte Carlo standard errors by block bootstrap, for the two
strategies of the rule. AR(1) = first-order autoregressive working
correlation, exch. = exchangeable, FG = Fay-Graubard, MD =
Mancl-DeRouen, LRT = likelihood-ratio test of the GLMM, MMM =
marginalized multilevel model.

<img src="figures/fig2.png" style="width:100.0%" />

**Figure 2.** Type I error at every design point of the main grid and
M1c (nested-loop plot), normal reference, marginalized model with the
convergence check. Design points are ordered by scenario (outer loop),
number of occasions T and number of clusters N (inner loop, bottom
panel); rows are the mechanisms. Grey band: 2.5% to 7.5%; line: nominal
5%. Points with fewer than 20 usable replicates are not drawn and the
axis is truncated at 30%.
