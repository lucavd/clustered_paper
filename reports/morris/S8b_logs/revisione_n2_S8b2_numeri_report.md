# Reviewer 2, second pass (S8b changes): numbers, tables and figures

**Scope.** Numbers, tables and figures changed in S8b, checked against the repository on lesexp-server (branch morris-v4, commit 459457a). The working tree has one uncommitted change, `manuscript/v5/MMM_v5_supplement.qmd`, and that file is byte-identical to the supplied supplement artifact apart from trailing whitespace. Every value was recomputed with my own code from the source files. Server scratch is under `/home/user/MMM_rev_S8b2/numbers/`; two read-only single-core R jobs were run. Log: 413 checks, 400 PASS, 12 MISMATCH, 1 UNVERIFIABLE. **No blocking and no major findings.**

## Verdict by task item

1. **Table S36**: all 126 cells (120 numeric, 6 n.a.) were recomputed from `tables/variant_kcfc/V4_perf_primary.csv` (H1, U0, own set, design points with at least 20 usable replicates) and all agree after rounding. The caption values also agree: 1,440/72 design points, at least 416 of 480 per mechanism, all 24 in M1c, and the MCSE maxima 0.0364, 5.86 and 1.04. The "upper bound" wording matches the computation (mean of per-point MCSE). One caveat: the Bias/SE MCSE treats the empirical SE as known.
2. **Table S37 and S11**: the 32 plug-in cells (raw and empirical-Bayes) were recomputed independently from the per-point n and k. The observed counts equal Table 4. The 24 bootstrap cells match `V5_S8b_rule_stability.csv`; the resampling was not re-run. The 62% and 82% (62.1 and 81.8) and the identical ranking in the three calculations are confirmed. **Table S26**: the seven existing rows were recomputed from the per-point exact intervals. The two new average rows (37/83/3/62/80/3 and 74/83/3/104/80/3) and the **S6 sentence (51 and 113 of 2,376)** were reproduced by an independent rerun on `results_v4/S7d+S7e/rows`.
3. **Figure S1**: the counts 17 (7+10), 4 (1+3, all GenBinaryY), 0, 0 and 0 equal Table 4, and 60 likelihood-ratio points are hidden. However, 4 further likelihood-ratio points (T = 2, N = 10) are isolated and plotted as lines without markers, so they are invisible. The caption MCSE maximum of 3.4 lies at one of these invisible points; over the visible points the maximum is 1.6.
4. **MCSE captions**: MCSE sentences of 18 captions were recomputed (Tables 2, 3, 4, S10, S12, S14, S15, S26, S28, S30, S31, S32, S33, S34, S35 and S36, Figures 2 and S1). The wording matches the computation everywhere except two cases: Figure 2 (wrong attribution) and Table S28 (one row omitted).
5. **Chain corrections**:
   - Table S10 (D-73): all 48 values were reproduced exactly by an independent rerun on `rows_v4.rds`.
   - Table S31 (D-74): all point estimates were recomputed, and the reading classes are consistent with the margin.
   - Table S1 (D-75): every entry matches `V4_dgm_parameters.csv`, including GLMM-AR1 M5 D[3,3] = 0.001.
   - Collapsed standard errors (the PB2 sentence, now in Section S3): 1,254 cells, 2.3 of 24 points.
   - GLMM-native marginal truth under H1: 0.0436 to 0.2018, which matches the printed 0.044 to 0.20.
6. **Clinical re-fits (D-70)**: all marginalized-model numbers in main Sections 3.11 and 4.4, supplement S9 and S17, and Tables S15, S29 and S30 were recomputed from `V5_S8b_clin_mmm_fits.csv` and `_ref.csv`. These cover the usable shares, outcome categories, bias, the model-based and robust SE ratios with the finite-population factor, and the full-data estimates and SEs. All agree.
7. **Cross-references**: all references resolve; Tables S1 to S37 are numbered in order; Table 5 no longer exists. The numbers shared by the abstract, Section 5 and the tables are consistent.

## Mismatches (all minor)

| Location | Finding | Proposed fix |
|---|---|---|
| Table S28 caption | The usable-share MCSE "at most 0.41" omits the printed row *Marginalized model without check* (mmm_study), whose bound is 0.49. The caption script builds the Copula frame without mmm_study. | 0.41 -> 0.49 |
| Figure 2 caption | "at most 6.4 points at the drawn design point with the fewest": the 6.36 maximum is at shrinkage GLM, GLMM-AR1 M2, N=10, T=2 (n=47, rate 25.5%). The point with the fewest replicates (n=20) has 4.87. | "at most 6.4 points at any drawn design point" |
| S3, collapsed-SE sentence | "in 193 cells they exceeded 5% of the usable null replicates": 193 counts cells whose type-7 critical value exceeds 37. Cells with more than 5% of null replicates at \|z\| > 37 number 199. | Reword to "the size-adjusted critical value exceeded 37", or report 199 |
| S11 paragraph | Highest rate "9%" for GEE AR(1) FG, while Table S32 prints 8.6 for the same quantity (8.605). The MMM value is 8.4 in both. | Print 8.6 in the text |
| Figure S1 | 4 likelihood-ratio points (T=2, N=10, one per panel) are drawn but invisible, and the MCSE maximum of 3.4 refers to one of them (the visible maximum is 1.6). | Add markers, or count them as not drawn (64) and give MCSE 1.6 |
| Table S37 caption | "strategies of Table 4 with an interval": 8 of 12 are shown; GEE exchangeable CR2, GEE exchangeable uncorrected, shrinkage GLM and GLMM Wald are missing. The S11 ranking claim holds over these 8. | State the 8 strategies, or add the missing rows |
| Table S26 | The two new average rows print usable "92" with no decimal; the other rows use one decimal (values 92.0 and 92.1). | 92.0 and 92.1 |
| Section S6 | "GEE with AR(1) working correlation was clearly above 7.5% at 51" means the MD/KC-average variance. Read alone, it contradicts the zero for GEE AR(1) FG. | Name the averaged-variance strategy |
| Table S14 vs Table 4, Figure S1 and S13 | M6-M7 counts differ: GEE FG t 3 vs 4, GEE MD t 18 vs 17, because S14 uses a different implementation (recomputed from `V4_S6b_ns_rule_safety_dp.csv`). The S14 caption does not say which implementation. | State the implementation in the S14 caption, or align it with fewclust |
| Section 2.3 | H0 marginal truth "-0.0030 to -0.0011 across scenarios": holds for M1-M5. M1c is +0.00057. | Add "(M1-M5; +0.0006 in M1c)" |
| Section S9 | The implementation and agreement sentence appears twice; the second copy cites "(Section S9)" from inside S9. | Remove the duplicate |

**UNVERIFIABLE (1):** the S1 count of 6 M1c design points with fewer than 1,000 datasets. Its source rds was not read.

**Notes on checks that passed:**
- The trial robust SE ratio is computed on the usable set defined by the robust SE (n = 257 vs 249 at N = 10); the rounded range 0.55 to 0.87 is unchanged.
- In main Section 4.4, the range 0.64 to 0.81 spans all T (minimum at T = 2), although the sentence concerns ten occasions (0.71 to 0.73 at T = 10).
- The S1 covariate note omits that dose is standardized.
- The Table S26 MCSE bound does not include the two average rows.

## What was not re-run
- The block bootstraps behind S37, S31 and the registry MCSE column: values were recomputed from their output files.
- The GLMM and GEE clinical fits: only the marginalized-model re-fits (D-70) were recomputed from fit-level output.
