---
title: "Analyzing clustered binary outcomes with few clusters and short series: which method works under which clustering mechanism? A simulation study and two clinical applications"
---

## Abstract

Analyses of clustered binary outcomes with few clusters have been evaluated mainly by number of clusters, seldom by number of occasions or by dependence mechanism. We compared six strategies for testing a treatment-by-time interaction: uncorrected generalized estimating equations (GEE), GEE with the Mancl-DeRouen or Fay-Graubard correction, the Wald test of a generalized linear mixed model (GLMM), a marginalized multilevel model (MMM) and a sign-flip score test. We simulated 10 to 145 clusters and 2 to 10 occasions under four mechanisms: random intercept, random intercept plus serial correlation, latent-variable transition model and purely serial correlation. Validity depended on the mechanism. Uncorrected GEE was anticonservative with 25 or fewer clusters under every mechanism. The GLMM was close to nominal on its own model but anticonservative under serial dependence, increasingly with more clusters (mean type I error up to 17%). Exchangeable GEE with a correction failed mainly with serial dependence and 25 or fewer clusters; an AR(1) working correlation removed this excess. The MMM was the only strategy with an interval never clearly above a 7.5% type I error, under any mechanism or design. AR(1) GEE with Fay-Graubard kept its level from three occasions and answered more often. Availability depended on the design: with two occasions and 25 or fewer clusters no valid strategy with an interval answered in more than 26% of datasets. Two subsampled clinical datasets agreed. We recommend the MMM with a convergence check; AR(1) GEE with Fay-Graubard, after a separation check, is a simpler alternative.

**Keywords.** clustered binary data; few clusters; generalized estimating equations; small-sample correction; marginalized multilevel model; generalized linear mixed model; simulation study.

## 1. Introduction

Binary outcomes measured repeatedly within a small number of subjects or clusters are common in clinical research: lung zones scored repeatedly within a patient, visits within a trial participant, periods within a cluster in a stepped-wedge trial. Three families of analysis are used. Generalized estimating equations (GEE) model the marginal (population-averaged) mean and treat the dependence as a nuisance, with inference from a sandwich variance. Generalized linear mixed models (GLMM) model the outcome conditionally on random effects and represent the dependence explicitly. Marginalized multilevel models (MMM) combine the two: a marginal mean model of scientific interest and a conditional model that generates the dependence [1,2].

Designs with few clusters are common: cluster-randomized trials with ten to thirty clusters are frequent [3], and many stepped-wedge trials have fewer than ten [4]. With few clusters, the asymptotic arguments behind all three methods are weak. The sandwich variance is biased downwards, and a family of small-sample corrections exists to repair it [5-7]. These corrections have been compared in longitudinal data [8-12], in parallel cluster-randomized trials [13-15] and in stepped-wedge trials [16], and the MMM has been compared with GEE in small samples [17]. Two features of these evaluations limit their use for short longitudinal series. First, they vary the number of clusters, rarely the number of occasions per cluster. In the original evaluations of the corrections and in the cluster-randomized comparisons [5,7,13,17,18], the unit within a cluster is never a subject measured on several occasions, and whether their conclusions hold with two or three occasions per subject is not known. Second, a method whose model matches the data-generating mechanism has an advantage in any comparison, so a conclusion obtained under one mechanism need not transfer to another; to separate the effect of the mechanism from that of sample size, the mechanism itself has to be varied.

This paper addresses three questions for a within-subject treatment-by-time interaction in longitudinal binary data.

1. With few clusters and short series, which analyses keep the nominal type I error, and does the answer depend on how the within-subject dependence is generated?
2. How often does each analysis return a usable result, and with what power?
3. When the mechanism is unknown, as it is in practice, which analysis can be used for a given number of clusters and occasions?

To answer them we simulated four dependence mechanisms over a grid that extends to ten clusters and two occasions, and repeated the comparison in two clinical datasets subsampled into the few-cluster regime. The study was planned and reported following Morris, White and Crowther [19].

## 2. Methods

### 2.1 The methods compared and how they differ

Table 1 summarizes the six strategies of the main text. They differ in four respects, each of which turns out to explain part of the results.

GEE and the MMM target a marginal log odds ratio; the GLMM targets a conditional one. Under a logit link the conditional coefficient is larger, approximately by a factor $(1 + k^{2}\sigma^{2})^{1/2}$ with $k = 16\sqrt{3}/(15\pi)$ and $\sigma^{2}$ the random-intercept variance [20]. Under the null of no interaction both are zero, so type I error is comparable across methods; bias and coverage are assessed against each method's own estimand. GEE does not need the dependence to be modeled correctly: the working correlation (exchangeable, independence or first-order autoregressive, AR(1)) affects efficiency, and the sandwich variance protects validity, but only asymptotically. The GLMM with a random intercept assumes that, given the subject effect, occasions are independent, which implies an exchangeable correlation; serial correlation that decays with lag is outside the model, and its model-based standard error has no protection against this. The MMM of Schildcrout and Heagerty [21] has both a latent subject effect and a first-order transition term, so it can represent a constant component and a short-range serial component. The Mancl-DeRouen correction (MD) inflates each cluster's residuals by the inverse of one minus its leverage; Fay-Graubard (FG) applies a bounded inflation (bound 0.75) [5,7]. The GLMM and the MMM rely on likelihood theory with no correction. All Wald tests in the main text use a t reference with $N - p$ degrees of freedom, $N$ being the number of clusters and $p$ the number of mean-model coefficients. The sign-flip score test [22] needs no variance estimate: it compares the observed effective score of the interaction with its distribution under random sign flips of the subject contributions.

Each method has its own failure mode in small designs. The corrected sandwich can be undefined when a cluster has leverage close to one. The GLMM with a random slope is not identified with two occasions per subject. The MMM's latent standard deviation, when it trends over time, is parametrized so that a zero trend lies on the boundary of the parameter space, where the optimizer fails; we therefore used a convergence check (Section 2.4). And when the data are separated, so that the maximum likelihood estimate of the interaction is infinite [23], GEE stops at a finite estimate with no warning, whereas the likelihood-based methods usually fail to converge.

Further strategies (Kauermann-Carroll and CR2 corrections, independence working correlation, the GLMM likelihood-ratio test, a shrinkage GLM, GLMM fitted by adaptive quadrature) were evaluated and are reported in the Supplementary Material, Section S3.

### 2.2 Data-generating mechanisms

Four mechanisms generated the within-subject dependence (Table 2, Figure 1). They share the marginal mean model and differ only in the dependence.

1. *GLMM-native*: a logistic mixed model with a normal random intercept (variance 1); the mechanism of the GLMM.
2. *GLMM-AR1*: GLMM-native plus a latent AR(1) process ($\rho = 0.6$, variance 0.5).
3. *GenBinaryY*: a marginalized model with a latent variable and a first-order transition term, simulated with the *MMLB* package [21,24]; the mechanism of the MMM.
4. *Copula-AR1*: a Gaussian copula with logistic margins and latent correlation $0.6^{|t_j - t_k|}$, with no component constant within subject; the mechanism of none of the likelihood-based methods.

Figure 1 shows what distinguishes them. Under the first three, the correlation between two occasions of the same subject levels off at 0.11 to 0.16 however far apart they are, because part of the dependence is constant within subject. Under Copula-AR1 it falls from 0.37 at lag 1 to 0.02 or less from lag 6. GenBinaryY combines both features: a short-range excess from the transition term and a plateau from the latent variable. The marginal prevalence was nearly identical across mechanisms (from 0.13 at the first occasion to 0.27 at the tenth).

The marginal mean model was $\mathrm{logit}\,\mu = \beta_{0} + \beta_{1}t + \beta_{2}\,\mathrm{trt} + \beta_{3}\,t \times \mathrm{trt}$, with coefficients $(-2, 0.1, 0.2, 0.2)$ under the alternative and $\beta_{3} = 0$ under the null (scenario M1). Scenarios M2 to M5 added a subject-level covariate, a small random slope on time (with a latent standard deviation growing by 0.01 per occasion in GenBinaryY), a smaller interaction (0.05) and a time-varying covariate with its own random slope. The main grid crossed the five scenarios and the four mechanisms with 16 numbers of clusters, $N \in \{10, 15, 20, 25, 35, \ldots, 145\}$, and six maximum numbers of occasions, $T \in \{2, 3, 4, 6, 8, 10\}$: 1,920 design points. Each subject was observed at all $T$ occasions with probability 0.7, and otherwise at a number drawn uniformly from 1 to $T - 1$. Treatment was assigned to subjects with probability 0.5.

Two further sets tested the conclusions outside the main grid. Scenario M1c (72 design points, first three mechanisms, $N \le 45$, $T \le 6$) took its prevalence (0.59) and between-subject variance (3.5) from the neonatal cohort of Section 2.5. Scenarios M6 and M7 (384 design points, GenBinaryY and GLMM-native) had a random-slope variance 10 to 40 times larger and a latent variability growing 10 to 25 times faster than in M3 to M5. Each design point received 1,000 replicates under the null and 1,000 under the alternative, with common random numbers across design points. Full parameters are in Table S1.

### 2.3 Performance measures

A replicate was *usable* for a Wald-based strategy when the fit returned a finite estimate, standard error and p-value with $|\hat\beta| < 10$ and $0 < \mathrm{SE} < 10$; for the sign-flip test when it returned a finite p-value. Replicates in which one treatment arm was never observed after the first occasion had no estimable design and were excluded from every denominator. A dataset in which the outcome was constant gave no usable result for any method. The usable share is reported with every replicate with an estimable design in the denominator, so that failures count against a method.

Type I error is the rejection rate at the nominal 5% among usable replicates under the null. At each design point a strategy was *clearly above 7.5%* when the lower limit of the exact 95% interval of its rejection rate exceeded 7.5%, and *clearly below 2.5%* when the upper limit fell below 2.5%; these are Bradley's liberal limits [25]. With 1,000 usable replicates a design point is flagged above only from an observed rate of 9.3%, so a strategy never clearly above 7.5% is one in which no excess was detected; design points with fewer than 20 usable replicates were not assessed. Power is the proportion of attempted datasets, usable or not, in which the null was rejected under the alternative, so that it reflects both availability and sensitivity.

To display the results by design, the grid was divided into twelve regions: three bands of occasions ($T = 2$; 3 to 4; 6 to 10) crossed with four bands of clusters ($N = 10$; 15 to 25; 35 to 65; 75 to 145). Within each region and mechanism we report the mean type I error over assessed design points, the number of design points clearly above 7.5%, the mean usable share and the mean power. With 1,000 replicates the Monte Carlo standard error of a rejection rate of 5% is 0.69 percentage points at a single design point. The datasets are common random numbers across design points, so the standard error of a region mean is bounded by the mean of the single-point standard errors and is not reduced as for independent cells (Section S2).

### 2.4 Implementation

GEEs were fitted with *geepack* [26] and the corrections computed with the R package *fewclust*, written for this study and checked against a computation by hand (Section S2). The GLMM was fitted with *lme4* [27] (Laplace approximation), with the random effects of the generating scenario. The MMM was fitted with *MMLB::mm* [24]. Its convergence check declared a fit at the boundary when the optimizer did not report convergence or an estimated slope of the latent standard deviation was below 0.001. Such a fit was replaced by the fit without the time trend only if that refit converged; otherwise the dataset gave no MMM result (Section S2). The sign-flip test used the effective score of the logistic model without the interaction, with subjects as units [22].

### 2.5 Clinical applications

Two datasets were repeatedly subsampled into the few-cluster regime, with the full-data estimate of each method as reference. The first is a prospective cohort of 207 neonates with ten lung-ultrasound zones scored at one to three timepoints [28]; the outcome is severe loss of aeration (prevalence 0.59), and the contrast is posterior against anterior zones, adjusted for disease, invasive ventilation and gestational age. The second is the respiratory trial distributed with *geepack* (111 patients, four visits); the contrast is active treatment against placebo, adjusted for baseline status, visit, age, sex and center. From each, 500 subsamples of whole patients were drawn at 10, 15, 20, 25 and 50 patients. We report the usable share and the ratio of the mean standard error to the standard deviation of the estimates over subsamples, multiplied by $\sqrt{1 - n/P}$, with $n$ the subsample size and $P$ the number of patients, to account for sampling without replacement.

## 3. Results

### 3.1 Which methods keep their level, and under which mechanism

Table 3 counts the design points at which each strategy was clearly above 7.5%, by mechanism, and Figure 2 shows the mean type I error in each region. Three patterns emerged: failures associated with few clusters, with the mechanism, and with the two together.

*Failures associated with few clusters*

Uncorrected GEE was clearly above 7.5% at 70 to 100 of 480 design points per mechanism, almost all with 25 or fewer clusters. With ten clusters and three or more occasions its mean type I error was 14% to 21% under every mechanism (Figure 2). Its excess was similar under the four mechanisms.

*Failures associated with the mechanism*

The GLMM Wald test behaved differently under each mechanism. On its own model (GLMM-native) and on GLMM-AR1 it was clearly above 7.5% at 42 and 43 of 416 assessed design points, with mean type I errors of 5.7% and 5.9%. Under GenBinaryY it was clearly above at 110, and under Copula-AR1 at 223 (54%), with mean type I error 9.7%. Under serial dependence the excess grew with the number of clusters. On Copula-AR1 with two occasions the mean rejection rate rose from 6.8% with 15 to 25 clusters to 17% with 75 to 145; with three or four occasions it rose from 7.0% with ten clusters to 12% with 75 to 145 (Figure 2d). On Copula-AR1 with 85 or more clusters and a random intercept, its standard error was on average 0.64 to 0.81 times the spread of its estimates under the null, while the estimates were centerd on zero (Section S4). Refitting the GLMM by adaptive quadrature on 123 null design points lowered its rejection rate by 3.1 to 4.0 percentage points and removed every flag on its own model. It left 8 design points clearly above 7.5%, all with ten occasions: 6 on Copula-AR1 and 2 on GenBinaryY (Section S4).

*Failures associated with the mechanism and few clusters together*

GEE with exchangeable working correlation and the Mancl-DeRouen or Fay-Graubard correction kept its level under the two random-intercept mechanisms, with at most 2 of 480 design points clearly above 7.5%. Under the two mechanisms with serial dependence it did not: Mancl-DeRouen was clearly above at 14 (GenBinaryY) and 28 (Copula-AR1) design points, and Fay-Graubard at 1 and 8, all but one with 25 or fewer clusters and three or more occasions. On Copula-AR1 with 10 to 25 clusters and three or four occasions their mean type I errors were 7.7% to 9.0% (Figure 2d). With an AR(1) working correlation, the Fay-Graubard strategy was clearly above 7.5% at no design point under GenBinaryY and Copula-AR1; its mean type I error in the same Copula-AR1 regions was 6.2% to 6.4%. The same pattern appeared in scenarios M6-M7 on GenBinaryY, where exchangeable Mancl-DeRouen and Fay-Graubard were clearly above 7.5% at 18 and 4 of 192 design points and the AR(1) version at none (Section S4). Averaged over each mechanism, the AR(1) working correlation lowered the type I error of Fay-Graubard by 0.19 to 0.23 percentage points under the random-intercept mechanisms and by 0.43 under the two serial ones; the difference was concentrated where clusters were few (Figure 2).

*Strategies that kept their level under every mechanism*

The MMM with the convergence check was never clearly above 7.5% under any mechanism (0 of 1,909 assessed design points in the main grid). It was conservative with few clusters: with ten clusters its mean type I error was 1.1% to 2.0% from three occasions upwards (Figure 2). GEE with AR(1) working correlation and Fay-Graubard was clearly above 7.5% at 2 of 1,920 design points, both with ten clusters, two occasions and a subject-level covariate (scenario M2) under the random-intercept mechanisms. At those two points every one of its 23 false positives occurred in a separated dataset, in which the estimate of the interaction does not exist and GEE reports the value at which its iterations stopped; on its non-separated datasets it rejected none (Section S4). The sign-flip test was never clearly above 7.5% and its mean type I error was 5.4% to 5.5% under each mechanism.

### 3.2 Availability and power depend on the design

Whether a method returned a result depended on the number of occasions and clusters, and hardly on the mechanism (Figure 3; Figure S1). With ten occasions every strategy answered in 97% or more of datasets. With two occasions the GEE strategies answered in 73% to 78% of datasets, the MMM in 52% and the GLMM Wald test in 27% (Figure 3a). The GLMM loss came from the random-slope scenarios, in which its Wald test was essentially never usable with two occasions (Table S6). The MMM loss came from the convergence check, which left 16% of null fits without a result in the random-slope scenarios, 32% with two or three occasions (Section S2).

With two or three occasions and 25 or fewer clusters, the GEE strategies with a correction returned a usable result in 39% of datasets, the MMM in 22% and the GLMM in 23%, and they detected the simulated interaction in fewer than 2% of attempted datasets. The sign-flip test answered in 69% of these datasets, kept its level (mean 5.1%), but had power of 3.8%. Over the whole grid, the strategies that kept their level had similar power, 20% to 21% per attempted dataset (Table 3); the higher power of uncorrected GEE and of the GLMM (23% for both) reflects their excess rejections.

Separated datasets were common in small designs. Of the datasets with an estimable design and a non-constant outcome, 7.5% were separated, and 55% of those with two or three occasions and at most 25 clusters (Section S4). The MMM rarely returned a result on separated data; GEE returned a plausible estimate and standard error, with no sign of the problem.

### 3.3 Choosing an analysis when the mechanism is unknown

An analyst knows the number of clusters and occasions but not the mechanism. Table 4 therefore lists, for each region, the strategies that were never clearly above 7.5% under any of the four mechanisms, and the availability and power of those that remained valid.

The MMM was the only strategy with an interval that remained valid in all twelve regions. With three or more occasions GEE with AR(1) working correlation and Fay-Graubard also remained valid in every region; with ten clusters and two occasions it did not. The GEE answered more often with 25 or fewer clusters (40% against 22% of datasets with ten clusters and three or four occasions; 75% against 50% with 15 to 25 clusters) and was closer to nominal, while the MMM was conservative. Their power was similar. Exchangeable Fay-Graubard was also valid from 35 clusters with three or four occasions and from 15 clusters with six or more, and from 75 clusters every GEE strategy, including the uncorrected one, was valid. The GLMM Wald test with the Laplace approximation was valid in only three regions: the two with two occasions and 25 or fewer clusters, where it answered in at most 8% of datasets, and ten clusters with six or more occasions, where it was conservative.

With two occasions and 25 or fewer clusters, no strategy with an interval was both valid and informative: among the valid ones the best usable share was 26% and the power 0.4% or less. The sign-flip test answered in 32% to 68% of these datasets, with power of 2.0% to 3.5%. With two or three occasions and at most 25 clusters, inverting it gave an interval with both ends finite in only 14% of datasets (Section S5).

### 3.4 Sensitivity analyses

The conclusions of Sections 3.1 to 3.3 held in the two additional sets of scenarios. In the cohort-anchored scenario M1c, with prevalence 0.59 and between-subject variance 3.5, GEE with AR(1) working correlation and Fay-Graubard, and the MMM, were never clearly above 7.5% (72 design points). The corrected GEEs were conservative (mean type I error 2.8% to 3.4%), and the GLMM Wald test was clearly above 7.5% at 11 design points. In M6-M7, with steep random slopes and growing latent variability, the same two strategies were again never clearly above 7.5% (384 design points), while the GLMM was clearly above at 86 (Section S4). Using a normal instead of a t reference raised every rejection rate, by 1.4 to 3.8 percentage points with 25 or fewer clusters; using $N - 2$ instead of $N - p$ degrees of freedom did not change which strategies remained valid (Section S3). The Kauermann-Carroll and CR2 corrections were clearly above 7.5% at 138 to 254 of the 2,376 design points across working correlations (Section S3).

### 3.5 Clinical applications

Table 5 summarizes the two applications. In the full datasets the two marginal methods agreed (posterior against anterior zones: 0.35 by the MMM and 0.38 by GEE; treatment effect: 1.3 by both), and the GLMM gave larger values on its conditional scale (0.48 and 2.1).

The subsamples reproduced three findings of the simulation. Availability fell with fewer patients for every method. The uncorrected sandwich understated the spread of the GEE estimates with 10 to 15 patients in the trial (ratio 0.5 to 0.6). The corrections returned a result less often than the other methods: with 10 patients Fay-Graubard was usable in 34% of cohort subsamples and 23% of trial subsamples, against 77% and 50% for the MMM. Most of their failures returned an estimate without a standard error. The applications also showed what the simulation did not: in both datasets the model-based standard errors of the MMM and the GLMM fell short of the spread of their estimates with up to 25 patients (ratios 0.5 to 0.7). In the cohort the robust standard error of the MMM was close to the spread (0.85 to 0.99).

## 4. Discussion

### 4.1 Principal findings

Different methods performed better under different mechanisms, and the reason can be traced to how each method treats the dependence. Uncorrected GEE failed because of few clusters, whatever the mechanism, which is the known downward bias of the sandwich [5,6]. The GLMM Wald test failed because of the mechanism. On its own random-intercept model it was nearly calibrated from three occasions. A random intercept, however, implies a correlation that does not decay with lag. When the true correlation does decay, the model-based standard error of a contrast over time is too small, and adding clusters does not correct it: the standard error shrinks while the error in it persists, so the excess grew with the number of clusters. A second, smaller source of excess was the Laplace approximation, which is inaccurate with few binary observations per cluster [29]. Exchangeable GEE with a small-sample correction failed mainly when the two problems met: serial dependence and 25 or fewer clusters. Replacing the exchangeable by an AR(1) working correlation removed that failure, and under random-intercept mechanisms cost nothing. The MMM kept its level under every mechanism, possibly because its dependence model has both a constant and a first-order serial part; it answered less often and was conservative with few clusters. In both clinical datasets the model-based standard errors of the MMM and the GLMM fell short of the spread of their estimates, consistent with dependence that their models did not capture.

Availability, in contrast, depended on the design and hardly on the mechanism. With two or three occasions and 25 or fewer clusters, more than half of the datasets were separated, and none of the interval methods answered reliably.

### 4.2 Practical recommendations

For a within-subject contrast over time with few clusters, we suggest the following, in the order of the questions of Section 1.

1. Do not use the uncorrected sandwich with fewer than about 75 clusters, nor the GLMM Wald test with the Laplace approximation when serial dependence is plausible, which in longitudinal data it usually is. Adaptive quadrature reduces but does not remove the excess under serial dependence.
2. Use the MMM with a convergence check and a t reference with $N - p$ degrees of freedom. It was the only strategy with an interval never clearly above 7.5% at any of the 2,376 design points, under any mechanism. It is conservative with ten clusters and returns no result more often than GEE; report how often the check intervened and how often it left no result.
3. GEE with an AR(1) working correlation, the Fay-Graubard correction and a t reference is a simpler alternative from three occasions, and the analysis to use when the MMM returns no result. It kept its level under all four mechanisms from three occasions and answered more often. Check for separation first, for example with the *detectseparation* package [30], and do not report a Wald test of a separated coefficient: at both design points where it exceeded 7.5%, all its false positives came from separated datasets. The R package *fewclust* applies this sequence, fitting the MMM first, GEE when the MMM gives no usable result, and stopping with an explanation when the mean model is separated.
4. With two occasions and 25 or fewer clusters, no analysis evaluated here estimates the effect reliably. The sign-flip test still gives a valid p-value, with little power. Where possible, plan four or more occasions per subject.

### 4.3 Relation to previous work

Our results agree with the anticonservatism of the uncorrected sandwich with few clusters [5,6] and with the good performance of the Fay-Graubard correction in stepped-wedge trials with binary outcomes [16]. A systematic review of parallel cluster-randomized trials recommended Fay-Graubard for continuous outcomes and, for binary outcomes, a mixed model with a between-within correction or Mancl-DeRouen [15]; in our longitudinal setting, Mancl-DeRouen with an exchangeable working correlation exceeded the level under serial dependence with few clusters. Li and Redden found Kauermann-Carroll at nominal level with ten clusters [13]; here it exceeded Bradley's limit. Their trials had large clusters, low intracluster correlation and a cluster-level treatment, so the difference may reflect the design. Bie et al. found the MMM's standard error nearly correct for within-cluster contrasts with few clusters [17], as here. McNeish preferred restricted penalized quasi-likelihood to quadrature for mixed logistic models with few clusters [29]; it was not evaluated here. The dependence of the GLMM's validity on the mechanism has a parallel in the debate on random-effects structures in mixed models [31,32]: a model that omits a source of variation in the contrast tested understates its uncertainty. The present results add that this source can be serial correlation, and that it is not visible from the number of clusters.

### 4.4 Limitations

The treatment-by-time interaction was linear and the mean depended only on current covariate values; with lagged effects of time-varying exposures, non-independence working correlations can bias GEE [33,34]. Occasions were missing completely at random. Random-slope variances in the main scenarios were small; M6-M7 cover steeper slopes for one structure and two mechanisms. No mechanism had a random-effect variance that differed between arms. Each Wald strategy was judged on its own usable replicates, so a method that fails on difficult datasets is assessed on easier ones; Table S7 repeats the comparison on replicates usable for all methods. The criterion "clearly above 7.5%" detects excess only from an observed rate of 9.3% at a single design point. Three elements were added after part of the results had been seen: the convergence check of the MMM, adaptive quadrature for the GLMM, and the analysis of separated datasets; the check was then tested on scenarios M6-M7, simulated afterwards. The AR(1) and independence versions of the corrections and Kauermann-Carroll rest on one implementation, checked against a computation by hand. The clinical applications use two datasets with one outcome each. They show that the availability problems arise in real data of this size, but they cannot measure type I error, because the truth is unknown.

### 4.5 Conclusion

With few clusters and short series, whether a test of a clustered binary outcome keeps its level depends on the mechanism that generates the dependence, and the methods differ in which mechanisms they tolerate. The marginalized multilevel model with a convergence check kept its level under all four mechanisms and at every design studied, and is the analysis we recommend. GEE with an AR(1) working correlation, the Fay-Graubard correction and a t reference kept its level under all four mechanisms from three occasions, answered more often, and is a simpler alternative when preceded by a check for separation. Whether any method returns a usable answer depends on the design. With two occasions and 25 or fewer clusters, no analysis does so reliably.

## Code and data availability

The simulation code, the per-replicate results and the scripts that produce every number, table and figure are available in a public repository and archived on Zenodo; the links are given on the title page. The companion R package *fewclust* implements the corrections and the convergence check. The neonatal cohort cannot be shared; the subsampling results are included as summary tables. The respiratory trial is distributed with the R package *geepack*.


## References

1. Heagerty PJ. Marginally specified logistic-normal models for longitudinal binary data. Biometrics 1999; 55: 688-698. DOI: 10.1111/j.0006-341X.1999.00688.x.
2. Heagerty PJ and Zeger SL. Marginalized multilevel models and likelihood inference. Stat Sci 2000; 15: 1-26. DOI: 10.1214/ss/1009212671.
3. Kahan BC, Forbes G, Ali Y, et al. Increased risk of type I errors in cluster randomised trials with small or medium numbers of clusters: a review, reanalysis, and simulation study. Trials 2016; 17: 438. DOI: 10.1186/s13063-016-1571-2.
4. Barker D, D'Este C, Campbell MJ, et al. Minimum number of clusters and comparison of analysis methods for cross sectional stepped wedge cluster randomised trials with binary outcomes: a simulation study. Trials 2017; 18: 119. DOI: 10.1186/s13063-017-1862-2.
5. Mancl LA and DeRouen TA. A covariance estimator for GEE with improved small-sample properties. Biometrics 2001; 57: 126-134. DOI: 10.1111/j.0006-341X.2001.00126.x.
6. Kauermann G and Carroll RJ. A note on the efficiency of sandwich covariance matrix estimation. J Am Stat Assoc 2001; 96: 1387-1396. DOI: 10.1198/016214501753382309.
7. Fay MP and Graubard BI. Small-sample adjustments for Wald-type tests using sandwich estimators. Biometrics 2001; 57: 1198-1206. DOI: 10.1111/j.0006-341X.2001.01198.x.
8. Pan W and Wall MM. Small-sample adjustments in using the sandwich variance estimator in generalized estimating equations. Stat Med 2002; 21: 1429-1441. DOI: 10.1002/sim.1142.
9. Morel JG, Bokossa MC and Neerchal NK. Small sample correction for the variance of GEE estimators. Biom J 2003; 45: 395-409. DOI: 10.1002/bimj.200390021.
10. Westgate PM and Burchett WW. Improving power in small-sample longitudinal studies when using generalized estimating equations. Stat Med 2016; 35: 3733-3744. DOI: 10.1002/sim.6967.
11. Ford WP and Westgate PM. A comparison of bias-corrected empirical covariance estimators with generalized estimating equations in small-sample longitudinal study settings. Stat Med 2018; 37: 4318-4329. DOI: 10.1002/sim.7917.
12. Gosho M, Ishii R, Noma H, et al. A comparison of bias-adjusted generalized estimating equations for sparse binary data in small-sample longitudinal studies. Stat Med 2023; 42: 2711-2727. DOI: 10.1002/sim.9744.
13. Li P and Redden DT. Small sample performance of bias-corrected sandwich estimators for cluster-randomized trials with binary outcomes. Stat Med 2015; 34: 281-296. DOI: 10.1002/sim.6344.
14. Leyrat C, Morgan KE, Leurent B, et al. Cluster randomized trials with a small number of clusters: which analyses should be used? Int J Epidemiol 2018; 47: 321-331. DOI: 10.1093/ije/dyx169.
15. Hemming K, Thompson J, Kristunas C, et al. The performance of small sample correction methods for controlling type I error when analyzing parallel cluster randomized trials: a systematic review of simulation studies. J Clin Epidemiol 2025; 185: 111838. DOI: 10.1016/j.jclinepi.2025.111838.
16. Thompson JA, Hemming K, Forbes A, et al. Comparison of small-sample standard-error corrections for generalised estimating equations in stepped wedge cluster randomised trials with a binary outcome: a simulation study. Stat Methods Med Res 2021; 30: 425-439. DOI: 10.1177/0962280220958735.
17. Bie R, Haneuse S, Huey N, et al. Fitting marginal models in small samples: a simulation study of marginalized multilevel models and generalized estimating equations. Stat Med 2021; 40: 5298-5312. DOI: 10.1002/sim.9126.
18. Lu B, Preisser JS, Qaqish BF, et al. A comparison of two bias-corrected covariance estimators for generalized estimating equations. Biometrics 2007; 63: 935-941. DOI: 10.1111/j.1541-0420.2007.00764.x.
19. Morris TP, White IR and Crowther MJ. Using simulation studies to evaluate statistical methods. Stat Med 2019; 38: 2074-2102. DOI: 10.1002/sim.8086.
20. Zeger SL, Liang KY and Albert PS. Models for longitudinal data: a generalized estimating equation approach. Biometrics 1988; 44: 1049-1060. DOI: 10.2307/2531734.
21. Schildcrout JS and Heagerty PJ. Marginalized models for moderate to long series of longitudinal binary response data. Biometrics 2007; 63: 322-331. DOI: 10.1111/j.1541-0420.2006.00680.x.
22. Hemerik J, Goeman JJ and Finos L. Robust testing in generalized linear models by sign flipping score contributions. J R Stat Soc Series B Stat Methodol 2020; 82: 841-864. DOI: 10.1111/rssb.12369.
23. Heinze G and Schemper M. A solution to the problem of separation in logistic regression. Stat Med 2002; 21: 2409-2419. DOI: 10.1002/sim.1047.
24. Mercaldo N and Schildcrout J. MMLB: marginalized models for longitudinal binary data. R package version 1.0, https://github.com/mercaldo/MMLB (2017, accessed 2 November 2026).
25. Bradley JV. Robustness? Br J Math Stat Psychol 1978; 31: 144-152. DOI: 10.1111/j.2044-8317.1978.tb00581.x.
26. Halekoh U, Højsgaard S and Yan J. The R package geepack for generalized estimating equations. J Stat Softw 2006; 15: 1-11. DOI: 10.18637/jss.v015.i02.
27. Bates D, Mächler M, Bolker B, et al. Fitting linear mixed-effects models using lme4. J Stat Softw 2015; 67: 1-48. DOI: 10.18637/jss.v067.i01.
28. Loi B, Sartorius V, Vivalda L, et al. Global and regional heterogeneity of lung aeration in neonates with different respiratory disorders: a physiologic observational study. Anesthesiology 2024; 141: 719-731. DOI: 10.1097/ALN.0000000000005026.
29. McNeish D. Estimation methods for mixed logistic models with few clusters. Multivariate Behav Res 2016; 51: 790-804. DOI: 10.1080/00273171.2016.1236237.
30. Kosmidis I, Schumacher D, Schwendinger F, et al. detectseparation: detect and check for separation and infinite maximum likelihood estimates. R package version 0.4.0, https://CRAN.R-project.org/package=detectseparation (accessed 3 November 2026).
31. Barr DJ, Levy R, Scheepers C, et al. Random effects structure for confirmatory hypothesis testing: keep it maximal. J Mem Lang 2013; 68: 255-278. DOI: 10.1016/j.jml.2012.11.001.
32. Matuschek H, Kliegl R, Vasishth S, et al. Balancing type I error and power in linear mixed models. J Mem Lang 2017; 94: 305-315. DOI: 10.1016/j.jml.2017.01.001.
33. Pepe MS and Anderson GL. A cautionary note on inference for marginal regression models with longitudinal data and general correlated response data. Commun Stat Simul Comput 1994; 23: 939-951. DOI: 10.1080/03610919408813210.
34. Schildcrout JS and Heagerty PJ. Regression analysis of longitudinal binary data with time-dependent environmental covariates: bias and efficiency. Biostatistics 2005; 6: 633-652. DOI: 10.1093/biostatistics/kxi033.


## Tables and figures

**Table 1.** Key differences between the strategies of the main text. N = number of clusters; p = number of mean-model coefficients.

| | GEE, uncorrected | GEE with Mancl-DeRouen or Fay-Graubard | GLMM, Wald test | MMM with convergence check | Sign-flip score test |
|---|---|---|---|---|---|
| Estimand | Marginal log odds ratio | Marginal log odds ratio | Conditional on the random effects; larger than marginal | Marginal log odds ratio | Tests zero marginal interaction; no estimate |
| Within-subject dependence | Working correlation as a nuisance; sandwich variance | As uncorrected GEE | Modelled: random intercept (and slopes); occasions independent given them | Modelled: latent subject effect plus first-order transition | Subjects as independent units; no model for dependence |
| Small-sample device | None | Leverage-based inflation of the sandwich | None (Laplace approximation) | None; check for boundary fits | Sign-flip distribution of the score |
| Reference | t, N - p df | t, N - p df | t, N - p df | t, N - p df | Randomization |
| Returns | Estimate, interval | Estimate, interval | Estimate, interval (conditional scale) | Estimate, interval | p-value only |
| Typical failure in small designs | Anticonservative with few clusters | Variance not computable; finite estimate on separated data | Random slope not identified at T = 2; misses serial dependence | No result after boundary fit; rarely answers on separated data | Not computable when the null fit diverges; low power |

**Table 2.** Data-generating mechanisms. All four share the marginal mean model of each scenario. Main grid: scenarios M1-M5, 16 numbers of clusters (10 to 145) and 6 numbers of occasions (2 to 10), 480 design points per mechanism; 1,000 replicates under the null and under the alternative at each. Full parameters in Table S1.

| Mechanism | How the dependence is generated | Constant within-subject component | Serial component | Model it favours | Additional scenarios |
|---|---|---|---|---|---|
| GLMM-native | Logistic mixed model, normal random intercept (variance 1), small random slopes in M3-M5 | Yes | No | GLMM | M1c, M6-M7 |
| GLMM-AR1 | GLMM-native plus latent AR(1) process (ρ = 0.6, variance 0.5) | Yes | Yes (latent) | None | M1c |
| GenBinaryY | Marginalized model: latent variable plus first-order transition term (*MMLB*) | Yes | Yes (transition) | MMM | M1c, M6-M7 |
| Copula-AR1 | Gaussian copula with logistic margins, latent correlation 0.6^\|lag\| | No | Yes | None | None |

**Table 3.** Type I error by mechanism, main grid (scenarios M1-M5, 480 design points per mechanism), Wald tests with a t reference (N - p degrees of freedom). Each mechanism cell: design points clearly above 7.5% / design points assessed (at least 20 usable null replicates); mean type I error (%) over assessed design points. Usable: replicates returning a usable result under the null (%); power: rejections under the alternative per attempted dataset (%); both averaged over the 1,920 design points. Results for all 19 strategies in Table S4.

| Strategy | GLMM-native | GLMM-AR1 | GenBinaryY | Copula-AR1 | Usable (%) | Power (%) |
|---|---|---|---|---|---|---|
| GEE, uncorrected | 70/480; 6.7 | 76/480; 6.8 | 91/480; 7.5 | 100/480; 7.6 | 93.8 | 22.7 |
| GEE, Mancl-DeRouen | 2/480; 4.6 | 1/480; 4.7 | 14/480; 5.2 | 28/480; 5.3 | 91.9 | 20.1 |
| GEE, Fay-Graubard | 1/480; 4.5 | 1/480; 4.5 | 1/480; 5.0 | 8/480; 5.0 | 92.0 | 19.9 |
| GEE AR(1), Fay-Graubard | 1/480; 4.3 | 1/480; 4.3 | 0/480; 4.5 | 0/480; 4.6 | 91.9 | 19.5 |
| GLMM, Wald (Laplace) | 42/416; 5.7 | 43/416; 5.9 | 110/416; 7.9 | 223/414; 9.7 | 80.4 | 23.1 |
| MMM with check | 0/477; 4.1 | 0/477; 4.3 | 0/478; 4.4 | 0/477; 4.9 | 83.5 | 20.3 |
| Sign-flip score test | 0/480; 5.5 | 0/480; 5.4 | 0/480; 5.5 | 0/480; 5.4 | 96.5 | 20.9 |

**Table 4.** Strategies whose type I error was never clearly above 7.5% under any of the four mechanisms, by region of the main grid, and the usable share and power (%) of the three strategies that were valid in most regions. Values are means over the design points and mechanisms of each region. GEE unc = uncorrected; MD = Mancl-DeRouen; FG = Fay-Graubard (exchangeable working correlation unless stated); GEE AR(1)-FG = AR(1) working correlation with Fay-Graubard.

| Occasions | Clusters | Strategies with an interval never clearly above 7.5% | Usable (%): GEE AR(1)-FG / MMM / sign-flip | Power (%): GEE AR(1)-FG / MMM / sign-flip |
|---|---|---|---|---|
| 2 | 10 | GLMM, MMM | 6 / 2 / 32 | 0.3 / 0.0 / 2.0 |
| 2 | 15-25 | GEE MD, GEE FG, GEE AR(1)-FG, GLMM, MMM | 26 / 13 / 68 | 0.2 / 0.3 / 3.5 |
| 2 | 35-65 | GEE unc, GEE MD, GEE FG, GEE AR(1)-FG, MMM | 77 / 48 / 97 | 1.5 / 1.8 / 6.0 |
| 2 | 75-145 | GEE unc, GEE MD, GEE FG, GEE AR(1)-FG, MMM | 98 / 76 / 100 | 4.1 / 4.0 / 6.2 |
| 3-4 | 10 | GEE AR(1)-FG, MMM | 40 / 22 / 64 | 1.6 / 0.4 / 3.6 |
| 3-4 | 15-25 | GEE AR(1)-FG, MMM | 75 / 50 / 91 | 3.8 / 2.0 / 5.5 |
| 3-4 | 35-65 | GEE FG, GEE AR(1)-FG, MMM | 98 / 85 / 100 | 6.4 / 5.9 / 7.4 |
| 3-4 | 75-145 | GEE unc, GEE MD, GEE FG, GEE AR(1)-FG, MMM | 100 / 98 / 100 | 8.9 / 9.6 / 9.6 |
| 6-10 | 10 | GEE AR(1)-FG, GLMM, MMM | 85 / 67 / 92 | 4.1 / 2.8 / 9.1 |
| 6-10 | 15-25 | GEE FG, GEE AR(1)-FG, MMM | 98 / 90 / 99 | 11.0 / 12.1 / 14.3 |
| 6-10 | 35-65 | GEE MD, GEE FG, GEE AR(1)-FG, MMM | 100 / 99 / 100 | 27.4 / 29.7 / 28.3 |
| 6-10 | 75-145 | GEE unc, GEE MD, GEE FG, GEE AR(1)-FG, MMM | 100 / 100 / 100 | 49.0 / 50.9 / 48.9 |

**Table 5.** Clinical applications: full-data estimate (standard error) and usable results in 500 subsamples of whole patients at each size (%). GLMM estimates are on the conditional scale. FG and MD: Fay-Graubard and Mancl-DeRouen corrections of exchangeable GEE. SE/SD: mean standard error divided by the standard deviation of the estimates over subsamples, finite-population corrected, range over subsample sizes from 10 to 25 patients. Cohort, 10 patients: 22% of subsamples contained no ventilated patient and could not be fitted by any method; they are in the denominator. Details in Section S6.

| Dataset | Method | Full-data estimate (SE) | Usable, 10 | 15 | 20 | 25 | 50 | SE/SD |
|---|---|---|---|---|---|---|---|---|
| Neonatal cohort (207 patients) | MMM | 0.35 (0.10) | 77 | 85 | 95 | 97 | 100 | 0.56-0.62 (robust SE 0.85-0.99)^a^ |
| | GEE, uncorrected | 0.38 (0.17) | 70 | 81 | 91 | 95 | 100 | n.a. |
| | GEE, FG | 0.38 | 34 | 61 | 80 | 88 | 100 | n.a. |
| | GEE, MD | 0.38 | 15 | 37 | 62 | 80 | 99 | n.a. |
| | GLMM | 0.48 (0.12) | 78 | 85 | 95 | 97 | 100 | 0.51-0.58^a^ |
| Respiratory trial (111 patients) | MMM | 1.3 (0.34) | 50 | 84 | 94 | 99 | 100 | 0.7 |
| | GEE, uncorrected | 1.3 (0.35) | 47 | 83 | 94 | 99 | 100 | 0.5-0.7 |
| | GEE, FG | 1.3 (0.35) | 23 | 68 | 88 | 98 | 100 | 0.9-2.2 |
| | GEE, MD | 1.3 (0.37) | 23 | 68 | 88 | 97 | 100 | 1.0-2.2 |
| | GLMM | 2.1 (0.54) | 53 | 88 | 97 | 100 | 100 | 0.6-0.7 |

^a^ Range over all subsample sizes, 10 to 50 patients. n.a. = not computed.





![]({{artifact:art_caf50544-6c23-4520-9e9d-4a2582cdc075}}){width=6.5in}

**Figure 1.** The four data-generating mechanisms, scenario M1 under the null, ten datasets of 2,000 subjects with up to ten occasions per mechanism. (a) Correlation between the binary outcomes of two occasions of the same subject, by lag, computed within treatment arm and occasion pair and averaged; bands: ± 2 Monte Carlo standard errors. (b) Marginal prevalence by occasion.

![]({{artifact:art_8f565bb6-3625-4d6b-a262-e04329534a6a}}){width=6.5in}

**Figure 2.** Mean type I error (%) of the seven strategies of Table 3 in each region of the main grid, by mechanism: (a) GenBinaryY, (b) GLMM-native, (c) GLMM-AR1, (d) Copula-AR1. Columns: number of occasions T (top) and of clusters N (bottom). Colour: white at the nominal 5%, red above, blue below, clipped at 0% and 15%. Boxed cells contain at least one design point clearly above 7.5% (lower limit of the exact 95% interval above 7.5%). Wald tests with a t reference; n.a.: no design point with 20 or more usable replicates.

![]({{artifact:art_87b2b623-700c-45cb-8bf3-fa408f26a0c6}}){width=6.5in}

**Figure 3.** (a) Replicates returning a usable result under the null and (b) power per attempted dataset, by number of occasions, averaged over the 16 numbers of clusters, five scenarios and four mechanisms of the main grid. The power of uncorrected GEE and of the GLMM Wald test includes their excess rejections.

