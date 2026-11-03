"""s7f_build_supplement.py -- S7f: supplement v5 from the v4 source. Every number is a registry token.
Usage: python s7f_build_supplement.py <v4 supplement qmd> <v5 supplement qmd>"""
import sys
src, dst = sys.argv[1], sys.argv[2]
T = lambda i: '`r n("%s")`' % i
sp = open(src, encoding="utf-8").read().split("\n\n")
Q = list(sp)
Q[1] = Q[1].replace("~/MMM/tables/V4_S7_numbers.csv", "~/MMM/tables/V5_S7f_numbers.csv")
Q[2] = ("This supplement accompanies the manuscript. Every number is read from the same registry of values as the main text. "
        "GB = GenBinaryY, NAT = GLMM-native, AR1 = GLMM-AR1, COP = Copula-AR1; MMM = marginalized multilevel model. "
        "In every table the Kauermann-Carroll correction is computed with `fewclust` (Section S4).")
Q[7] = Q[7].replace("Versions in the project library", "Versions in the project library, with the commit of `fewclust`")
SW = [("R", "R"), ("MMLB", "MMLB"), ("lme4", "lme4"), ("geepack", "geepack"), ("geesmv", "geesmv"), ("fewclust", "fewclust"),
      ("fewclust_commit", "fewclust, commit"), ("clubSandwich", "clubSandwich"), ("data.table", "data.table"), ("crew", "crew"),
      ("mirai", "mirai"), ("Matrix", "Matrix")]
# Tables S11 and S24: the last column is normal minus t, the other shift columns t minus normal (second review): say so in the header
for k_ in range(len(Q)):
    Q[k_] = Q[k_].replace("| shift GB | shift NAT | shift AR1 | shift at N <= 25 |", "| t minus normal GB | t minus normal NAT | t minus normal AR1 | normal minus t at N <= 25 |")
Q[8] = "| Software | Version |\n|---|---|\n" + "\n".join(f"| {lab} | {T('f:sw|' + k)} |" for k, lab in SW)
assert Q[36].startswith("## S4."), Q[36][:30]
Q[37] = (
    "The corrections of Mancl and DeRouen, Kauermann and Carroll and Fay and Graubard were computed with two implementations. "
    "`geesmv` was used in the original study with exchangeable working correlation. On sparse binary data with few clusters its three functions "
    "can loop indefinitely in compiled code, so each call ran in a separate process killed after `r n(\"x:geesmv_timeout\")` seconds; the largest "
    "share of time-outs in a cell was `r n(\"x:to_max\")`%. `geesmv` also takes the clusters from any column named `id`, whatever cluster variable "
    "is passed: in the simulations and in the cohort the cluster column has that name, and in the respiratory trial it was renamed before every call. "
    "With AR(1) working correlation `geesmv` fails as soon as one subject has a single observation, so it was not used with that correlation.\n\n"
    f"`fewclust` (Table S2) computes the three corrections from the estimating-equation building blocks. On {T('f:ge4|nmin')} or more datasets for each "
    f"working correlation its standard errors agreed with a computation by hand; the largest difference over the three corrections was below $10^{-12}$. "
    f"Against `geesmv` on the same exchangeable fits ({T('f:ec|md|n')} datasets), Mancl-DeRouen and Fay-Graubard agreed to within one per cent in "
    f"{T('f:ec|md|share')}% and {T('f:ec|fg|share')}% of datasets, and Kauermann-Carroll in only {T('f:ec|kc|share')}%. "
    "`geesmv::GEE.var.kc` rebuilds $(I - H_{ii})^{-1/2}$ with the transpose of the eigenvectors of a non-symmetric matrix, an identity that holds only "
    "for symmetric matrices. "
    f"With the published formula, Kauermann-Carroll with exchangeable correlation and a t reference had a mean type I error of "
    f"{T('f:r7d|gee_exch_kc|mean_rate')}% over the {T('f:x|ndp_three')} design points of the three study mechanisms, against "
    f"{T('f:r7d|study_gee_kc|mean_rate')}% with `geesmv`, and it was clearly above 7.5% at {T('f:r7d|gee_exch_kc|n_hi')} instead of "
    f"{T('f:r7d|study_gee_kc|n_hi')} design points. In the earlier version of this work the `geesmv` value appeared to shrink rather than inflate "
    "the variance; this was a property of the implementation, not of the estimator.")
# ---- S6. Degrees of freedom ------------------------------------------------------------------------------------------------
S26 = [("gee_ar1_fg", "GEE AR(1), Fay-Graubard"), ("mmm_cc", "Marginalized model with check"), ("gee_ind_fg", "GEE independence, Fay-Graubard"),
       ("gee_exch_fg", "GEE exchangeable, Fay-Graubard"), ("gee_ar1_md", "GEE AR(1), Mancl-DeRouen"), ("gee_exch_md", "GEE exchangeable, Mancl-DeRouen"),
       ("gee_exch_kc", "GEE exchangeable, Kauermann-Carroll"), ("glmer", "GLMM, Wald")]
sec6 = ("## S6. Degrees of freedom of the t reference\n\n"
        f"**Table S26.** The decision rule with N - p and with N - 2 degrees of freedom, as used by Li and Redden, on the {T('f:x|ndp_three')} design points "
        "of the three study mechanisms (the fourth mechanism was analysed with N - p only). Above: design points clearly above 7.5%; below: clearly below 2.5%; "
        "usable: replicates returning a usable result (%), averaged over design points. With both choices the rule selects GEE with AR(1) working correlation "
        "and Fay-Graubard.\n\n"
        "| Strategy | Above, N - p | Below, N - p | Below at ten clusters, N - p | Above, N - 2 | Below, N - 2 | Below at ten clusters, N - 2 | Usable (%) |\n"
        "|---|---|---|---|---|---|---|---|\n"
        + "\n".join(f"| {lab} | {T(f'f:t2|{s}|N - p|n_hi')} | {T(f'f:t2|{s}|N - p|n_lo')} | {T(f'f:t2|{s}|N - p|n_lo_N10')} | {T(f'f:t2|{s}|N - 2|n_hi')} | "
                    f"{T(f'f:t2|{s}|N - 2|n_lo')} | {T(f'f:t2|{s}|N - 2|n_lo_N10')} | {T(f'f:t2|{s}|N - p|pu@tab')} |" for s, lab in S26))
# ---- S7. Sign-flip test and its inversion ----------------------------------------------------------------------------------
rows27 = []
for t in [2, 3, 4, 6, 8, 10]:
    for k, g in enumerate(["10", "15-25", "35-65", "75-145"]):
        rows27.append(f"| {('T = %d' % t) if k == 0 else ''} | N = {g} | {T(f'f:ef|{t}|{g}|datasets')} | {T(f'f:ef|{t}|{g}|finite@tab')} | "
                      f"{T(f'f:ef|{t}|{g}|gee_ar1_fg_usable@tab')} |")
sec7 = ("## S7. Sign-flip score test and the interval obtained by inverting it\n\n"
        "The sign-flip test of the main text tests the null of no interaction. Testing instead the null that the interaction equals a value $b_0$, by entering "
        "$b_0 \\cdot t \\cdot \\mathrm{trt}$ as an offset in the model without the interaction, and collecting the values that are not rejected at 5%, gives a "
        "confidence interval. Its ends were searched from the independence GEE estimate by a walk in fixed steps followed by bisection, with the same sign "
        "vectors at every value and every fit started from the default initial value; a value at which the test was not computable counted as not rejected, "
        "so the interval can be unbounded. "
        f"Under the alternative the interval covered the true value in {T('f:ci|PG1')}% of datasets over the {T('f:x|ndp_all5')} design points, against "
        f"{T('f:ci|PG8')}% for GEE with AR(1) working correlation, Fay-Graubard and a t reference, and its coverage was clearly outside Bradley's band at no "
        f"design point. Both ends were finite in {T('f:ci|PG5')}% of datasets, and in {T('f:ci|PG6')}% with two or three occasions and at most 25 clusters "
        f"(Table S27). Where both intervals were finite it was wider than the GEE interval by {T('f:x|ci_width_pct')}% on average (geometric mean of the "
        f"ratio of widths), and by {T('f:x|ci_width_pct_N75')}% with 75 clusters or more.\n\n"
        "With the rules of Section 2.4 of the main text the test was not computable in part of the small designs, mostly because one arm had only zero "
        f"outcomes. With two or three occasions and at most 25 clusters it answered in {T('f:fl|small|H0|mean_pu')}% of datasets, against "
        f"{T('f:fl7|small|H0|mean_pu')}% when the rules were not applied, with a mean type I error of {T('f:fl|small|H0|mean_rate')}%.\n\n"
        "**Table S27.** Interval obtained by inverting the sign-flip test, under the alternative, by number of occasions T and number of clusters N: datasets "
        "with an estimable design, datasets with both ends of the interval finite (%), and datasets in which GEE with AR(1) working correlation, "
        "Fay-Graubard and a t reference returned a usable interval (%).\n\n"
        "| Occasions | Clusters | Datasets | Interval finite (%) | GEE AR(1), Fay-Graubard usable (%) |\n|---|---|---|---|---|\n" + "\n".join(rows27))
# ---- S8. Fourth mechanism --------------------------------------------------------------------------------------------------
S28 = [("gee_ar1_fg", "GEE AR(1), Fay-Graubard"), ("mmm_cc", "Marginalized model with check"), ("mmm_study", "Marginalized model without check"),
       ("gee_ind_fg", "GEE independence, Fay-Graubard"), ("gee_exch_fg", "GEE exchangeable, Fay-Graubard"), ("gee_exch_md", "GEE exchangeable, Mancl-DeRouen"),
       ("gee_exch_kc", "GEE exchangeable, Kauermann-Carroll"), ("gee_exch_unc", "GEE exchangeable, uncorrected"), ("glmer", "GLMM, Wald"),
       ("shrink", "Shrinkage GLM"), ("glmer_lrt", "GLMM, likelihood ratio"), ("flip", "Sign-flip score test")]
BIAS = {"gee_ar1_fg", "mmm_cc", "mmm_study", "gee_ind_fg", "gee_exch_fg", "gee_exch_md", "gee_exch_kc", "gee_exch_unc"}
bb = lambda s, f: T(f"f:bs|{s}|{f}@tab") if s in BIAS else "n.a."
sec8 = ("## S8. The fourth mechanism\n\n"
        f"Copula-AR1 draws, for each subject, a vector of normal scores with correlation $\\rho^{{|t_j - t_k|}}$ (ρ = {T('f:x|cop_rho')}) and sets each "
        "outcome to one when its score is below the normal quantile of the marginal mean of the scenario. The marginal means are therefore exactly those of "
        "the scenario, and the correlation between occasions of a subject fades with their distance, with no component that is constant within subject. "
        f"Without the convergence check the marginalized model was clearly above 7.5% at {T('f:s7e|mmm_study|clearly_above')} of the {T('f:x|ndp_cop')} "
        "design points.\n\n"
        f"**Table S28.** Copula-AR1 ({T('f:x|ndp_cop')} design points), Wald tests with a t reference: mean type I error (%), design points clearly above "
        "7.5% and clearly below 2.5%, replicates returning a usable result (%), power (%), and, for the marginal estimators, mean bias in empirical standard "
        "errors and coverage of the 95% interval with a normal and a t reference (%); n.a. = not applicable.\n\n"
        "| Strategy | Type I error | Above 7.5% | Below 2.5% | Usable | Power | Bias / empirical SE | Coverage, normal | Coverage, t |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
        + "\n".join(f"| {lab} | {T(f'f:s7e|{s}|mean_type1@tab')} | {T(f'f:s7e|{s}|clearly_above')} | {T(f'f:s7e|{s}|clearly_below')} | "
                    f"{T(f'f:s7e|{s}|mean_usable@tab')} | {T(f'f:s7e|{s}|mean_power@tab')} | {bb(s, 'mean_bias_rel')} | {bb(s, 'mean_cov_n')} | "
                    f"{bb(s, 'mean_cov_t')} |" for s, lab in S28))
# ---- S9. Respiratory trial --------------------------------------------------------------------------------------------------
RM = [("MMM", "Marginalized model"), ("GEE", "GEE, uncorrected"), ("GEE-MD", "GEE, Mancl-DeRouen"), ("GEE-KC", "GEE, Kauermann-Carroll"),
      ("GEE-FG", "GEE, Fay-Graubard"), ("GLMM", "GLMM")]
rows30 = [f"| {lab if k == 0 else ''} | {T(f'f:rN|{n}')} | {T(f'f:ru|{m}|{n}@tab')} | {T(f'f:rb|{m}|{n}@tab')} | {T(f'f:sr|resp|{m}|{n}@tab')} |"
          for m, lab in RM for k, n in enumerate([10, 15, 20, 25, 50])]
sec9 = ("## S9. Second clinical application: the respiratory trial\n\n"
        f"The respiratory trial distributed with `geepack` has {T('f:rs|P')} patients observed at four visits. The full mean model had treatment, baseline "
        "status, visit, age, sex and centre; the reduced model had treatment and baseline status only. Kauermann-Carroll was computed with `fewclust` on "
        "the same GEE fits and the same draws, in both applications; Mancl-DeRouen and Fay-Graubard were computed with `geesmv`. Recomputed with "
        f"`fewclust`, these two agreed with `geesmv` to within one per cent in {T('f:ag|neonatal cohort|md')}% and {T('f:ag|neonatal cohort|fg')}% of "
        f"the draws of the cohort and in {T('f:ag|respiratory|md')}% and {T('f:ag|respiratory|fg')}% of those of the trial.\n\n"
        "**Table S29.** Full-data estimate (SE) of the treatment effect, full and reduced mean model. The mixed model estimates a conditional effect.\n\n"
        "| Method | Full model | Reduced model |\n|---|---|---|\n"
        + "\n".join(f"| {lab} | {T(f'f:rr|{m}|estimate')} ({T(f'f:rr|{m}|se')}) | {T(f'f:rrr|{m}|estimate')} ({T(f'f:rrr|{m}|se')}) |" for m, lab in RM)
        + "\n\n**Table S30.** Subsamples of whole patients, full mean model: draws returning a usable result (%), bias of the estimate against the full-data "
        "estimate, and ratio of the mean model standard error to the standard deviation of the estimates multiplied by the finite-population factor of "
        "Section 2.7 of the main text.\n\n"
        "| Method | N | Usable (%) | Bias | SE ratio, corrected |\n|---|---|---|---|---|\n" + "\n".join(rows30))
# ---- S10. The two strategies of the rule compared (D-61; block bootstrap and margins after the second review) -----------
REG = [("map|10|2", "N = 10", "T = 2"), ("map|10|3", "N = 10", "T = 3"), ("map|10|4+", "N = 10", "T ≥ 4"),
       ("map|15-25|2", "N = 15-25", "T = 2"), ("map|15-25|3", "N = 15-25", "T = 3"), ("map|15-25|4+", "N = 15-25", "T ≥ 4"),
       ("map|>25|2", "N > 25", "T = 2"), ("map|>25|3", "N > 25", "T = 3"), ("map|>25|4+", "N > 25", "T ≥ 4"), ("all", "all", "all")]
hb = lambda r, q: T(f"f:hb|{r}|{q}@tab"); hs = lambda r, q: T(f"f:hbs|{r}|{q}@tab")
rows31 = [f"| {ln} | {lt} | {hb(r, 'rate_gee')} / {hb(r, 'rate_mmm')} | {hb(r, 'mad_gee')} / {hb(r, 'mad_mmm')} | {T(f'f:hb|{r}|mad_gee_minus_mmm')} ({T(f'f:hbs|{r}|mad_gee_minus_mmm')}) | "
          f"{T('f:hbc|' + r)} | {hb(r, 'pu_gee')} / {hb(r, 'pu_mmm')} | {hb(r, 'upow_gee')} / {hb(r, 'upow_mmm')} |" for r, ln, lt in REG]
sec10 = ("## S10. The two strategies of the rule compared\n\n"
         "The two strategies that met the decision rule were compared by number of clusters and of occasions after the results had been seen, with the block bootstrap and the "
         "one-point margin for mean differences of Section 2.5 of the main text. Neither was clearly above 7.5% at any design point. The comparison does not enter the rule, which "
         "selects the GEE strategy (Table 4).\n\n"
         f"**Table S31.** GEE with AR(1) working correlation and Fay-Graubard (GEE) and the marginalized multilevel model with the convergence check (MMM), both with a t "
         f"reference (N - p degrees of freedom), over the {T('f:x|ndp_all5')} design points of the four mechanisms, by number of clusters N and of occasions T: mean type I error (%), "
         "mean absolute departure from 5% (percentage points), their difference GEE minus MMM with its Monte Carlo standard error (block bootstrap), its reading against the "
         "one-point margin (the strategy closer to 5% when the interval of ± 2 standard errors lies beyond the margin, similar when it lies within it, about the margin when it "
         "straddles it), datasets returning a usable result (%) and datasets in which the simulated effect was detected among those attempted (%).\n\n"
         "| Clusters | Occasions | Type I error, GEE / MMM | Departure from 5%, GEE / MMM | Difference (MCSE) | Reading | Usable, GEE / MMM | Power, GEE / MMM |\n"
         "|---|---|---|---|---|---|---|---|\n" + "\n".join(rows31))
out = "\n\n".join(Q[:59] + [sec6, sec7, sec8, sec9, sec10]) + "\n"
open(dst, "w", encoding="utf-8").write(out)
print(len(out))
