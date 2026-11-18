Padova, 10 November 2026

Dr Chris Metcalfe\
*Statistical Methods in Medical Research*

Dear Dr Metcalfe,

**Resubmission of manuscript SMM-26-1128, now entitled "Analyzing clustered binary outcomes with few clusters and short series: which method works under which clustering mechanism? A simulation study and two clinical applications"**

Thank you for your decision of 8 November 2026 and for inviting a resubmission. We agree with your assessment. The paper had grown through successive additions to the point where its questions and findings were hard to follow. We have rewritten it from the beginning, using the same simulation results and clinical data, and organized it around the points you raised.

*A concise statement of the questions.* The Introduction now ends with three numbered questions: which analyses keep the nominal type I error with few clusters and short series, and whether this depends on the dependence mechanism; how often each analysis returns a usable result, and with what power; and which analysis can be used when the mechanism is unknown. Sections 3.1, 3.2 and 3.3 answer them in that order, and the recommendations of Section 4.2 follow the same order.

*The key differences between the methods.* A new Section 2.1 and Table 1 set out, side by side, what each method estimates, how it treats the within-subject dependence, how it handles few clusters, what it returns and how it typically fails in small designs. The Results and Discussion refer back to these differences to explain the findings. The main text now has six strategies (uncorrected GEE, GEE with the Mancl-DeRouen or Fay-Graubard correction, the GLMM Wald test, the marginalized multilevel model and a sign-flip score test), plus the AR(1) working-correlation version of Fay-Graubard. The remaining strategies are reported in the Supplementary Material.

*The clustering mechanisms, and which method works under which.* The comparison of mechanisms is now the core of the paper. A new Figure 1 shows how the four mechanisms differ: the within-subject correlation levels off under the three mechanisms with a constant component and decays to zero under purely serial dependence, while the marginal prevalence is the same. Results are reported by mechanism (Table 3) and by region of clusters and occasions within each mechanism (Figure 2). This makes the mechanism-specific pattern you anticipated explicit:

1.  uncorrected GEE fails with few clusters under every mechanism;
2.  the GLMM Wald test keeps its level on its own model but not under serial dependence, and the excess grows with the number of clusters;
3.  exchangeable GEE with a small-sample correction fails only when serial dependence and 25 or fewer clusters occur together, and an AR(1) working correlation removes this failure;
4.  the marginalized multilevel model with a convergence check is the only strategy with an interval that never exceeds the level under any mechanism or design, and GEE with AR(1) working correlation and Fay-Graubard keeps its level under all four mechanisms from three occasions upwards and answers more often.

Our recommendation is therefore the marginalized multilevel model, with GEE with AR(1) working correlation and the Fay-Graubard correction, preceded by a check for separation, as a simpler alternative from three occasions.

*A summary of the findings that helps understanding.* We replaced the single decision rule of the previous version, which compressed the comparison into counts over 2,376 design points, with a table that a reader can consult by design (Table 4). For each combination of clusters and occasions it lists the strategies that kept their level under every mechanism, with their availability and power. A separate result shows that availability depends on the design and hardly on the mechanism (Figure 3). The two clinical applications are summarized in one table (Table 5) and linked to the simulation findings they confirm.

The main text has been shortened from about 8,000 to about 4,700 words, with five tables and three figures. The supplement has been reorganized into seven sections. Sensitivity analyses (cohort-anchored and steep-variability scenarios, reference distribution, degrees of freedom, adaptive quadrature, separated datasets, versions of the convergence check) are collected in one section and summarized in one paragraph of the main text. No new simulations were run; every number comes from the per-replicate results of the previous submission, which are publicly archived together with the code.

The manuscript has not been published and is not under consideration elsewhere. All authors have read and approved the resubmitted version. The use of generative artificial intelligence tools is declared on the title page. The manuscript and the supplementary material are anonymized for double-anonymized review; the title page contains the authors, the declarations and the link to the repository.

Yours sincerely,

Dario Gregori, on behalf of the authors\
Unit of Biostatistics, Epidemiology and Public Health, University of Padova\
dario.gregori@unipd.it
