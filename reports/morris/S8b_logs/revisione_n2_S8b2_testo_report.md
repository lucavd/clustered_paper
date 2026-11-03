# Reviewer n.2, second pass on S8b: text, support and citations

Scope: the parts of the main text and supplement that changed in S8b (the 64-item revision, then the shortening). Inputs: main qmd v9d1ebaff and render v0ccf23e1, supplement qmd ve562e8dc and render v2a604661, diffs v76bf3b35 / v5e4714b0, registry vaddb01b1, server repo /home/user/MMM at commit 459457a, read only.

Result: 115 checks, of which 70 PASS, 35 MISMATCH (4 major, 31 minor), 5 UNSUPPORTED (all blocking) and 5 UNVERIFIABLE (minor). The full log is in `revisione_n2_S8b2_text_log.csv`.

## Blocking: claims in Section 4.1 and Section 5 without main-text support

After the shortening, these values appear only in the supplement or in server tables. No sentence in Section 3, and no main table or figure, carries them.

1. **4.1, earlier criterion** (P41-11). "Would have recommended the likelihood-ratio test up to 25 clusters." The only source is `tables/V4_S5_switch_points.csv` (switch at N = 25 on all three mechanisms). The 3.3 sentence that partly supported it has also lost its scope (see MT-01). *Fix:* restore "At 20 clusters or fewer, unconditional validity ranged from 61% ... to 80% (Table S8)" in 3.3 and add the up-to-25 statement there.
2. **4.1, power of the two strategies in small designs** (P41-14). The text says "fewer than 2%". Recomputed: 1.45% (GEE) and 1.42% (MMM), from `V5_S7f_rule_dp.csv` and `V5_S7f_h2h_boot.csv`. Table 4 gives power over all design points only, and Figure 1f gives it by T only. *Fix:* add one sentence to 3.4.
3. **4.1 and Section 5, sign-flip power** (P41-15, C5-05). The text says "only 4.2%" and "little power". Recomputed: 4.15% from `V5_S7f_rule_dp.csv`. Section 3.7 gives only the 21% over all design points. *Fix:* add the small-design power to 3.7.
4. **4.1, inverted sign-flip interval** (P41-16). The values 95% vs 96%, 14% finite, 101 / 91 design points and 17% wider are in Section S7 and Table S27 only. *Fix:* add a short paragraph to 3.7, or keep only "(Section S7)" in 4.1.

Every other sentence of the abstract, Section 4.1 and Section 5 has main-text support. In particular, "differed beyond the margin only at ten clusters" holds: in `V5_S7f_h2h_boot.csv`, only the cells N = 10 with T = 3 and N = 10 with T ≥ 4 read "GEE"; the T = 2 cells read "about the margin". The quadrature claims also hold: all 8 design points still above 7.5% are at T = 10, and 31 of the 37 design points clearly below 2.5% are at N ≤ 25.

## Major

- **I-09, Introduction.** "In the evaluations cited here, the unit within a cluster is never a subject measured on several occasions" is false as worded. Reference [12], cited one sentence earlier, simulates 10 to 100 subjects with four occasions each. The original wording ("the original evaluations of the corrections and the cluster-randomized comparisons") was accurate.
- **D45-04, Section 4.5.** "The zero count is uncertain; the ranking is not" misstates Table S37. In 62% (GEE) and 82% (MMM) of bootstrap resamples, the strategy is clearly above 7.5% at one or more design points. The ranking is identical across the three calculations, but there is no interval for it. Quote the numbers; the PI may also wish to qualify "never clearly above 7.5%" in the abstract.
- **MT-01, Section 3.3.** "Unconditional validity ranged from 61% to 80%" lost its scope ("at 20 clusters or fewer", v:val20).
- **MT-02, Sections 2.4 and 3.5.** The main text no longer says that the exchangeable Mancl-DeRouen and Fay-Graubard values in Tables 2 and 3 come from geesmv (Section S4). Section 2.4 therefore misdescribes those tables, and "the two implementations" in 3.5 has no antecedent.

## Minor (summary)

**Citations.**
- [15] is grouped with the systematic review although it is a primary simulation study.
- Thompson et al. [18] did not single out Fay-Graubard; they found approximate Kauermann-Carroll comparable. They also reported *unstable* geesmv Kauermann-Carroll standard errors, not a *different quantity*.
- McNeish [31] recommends restricted penalized quasi-likelihood mainly for many observations per cluster.
- Bie et al. [19] found the marginalized model robust to non-normal random effects; it was sensitive only to a random-effect variance that differs between arms.
- Bates et al. [26] (the lme4 paper) does not cover adaptive quadrature.
- The first author of Thompson et al. should be "Thompson JA".

**Section 4.4.**
- "By similar amounts" hides a consistent gap: the marginalized model's ratio is higher by 0.04 in the cohort and 0.06 in the trial, and in the trial at N = 50 it is 0.85 against 0.70.
- The trial ratios are no longer in the main-text Results.
- The 0.64-0.81 range needs a pointer to Section S8.
- "Seven coefficients" describes the trial model only.

**Section 4.2.** "Most often because nlm did not converge" is not true in GLMM-native M6 (31.9% small slope vs 31.0% non-convergence). There is no breakdown by reason for M3-M5.

**Main-text coherence.**
- 3.10 lost the list of T values refitted.
- 2.6 no longer introduces the post hoc comparison.
- 2.5 no longer states that fits with warnings were kept.

**Supplement, moved sentences.**
- Dangling "also" and "its" (S3, S8, S13), and "The change" without an antecedent (S12).
- Self-references to the same section (S4, S7, S8, S9).
- Duplicated sentences (S8, S9, S10, S11, S12, S13).
- Scope missing (S5 "In M1c"; S9 "In the trial").
- "With bobyqa" used before its definition in S15.
- "Often" vs "usually" for the same fact in S7.
- "Its smallest excess" now reads as relative to 7.5% (S3).
- S16 does not say what changed raw power.

S8b's two descriptor changes are correct: "above 7.5% by 19 points" (C12a = rate − 7.5) and "raw power" (P11 is a raw-power difference).

Each log row carries a proposed wording.

## Verified without findings

- Numbering by first citation runs 1 to 33 in order, with no orphan or missing entries.
- All entries follow the Vancouver rule (at most three authors, otherwise three plus et al.).
- Citation numbers in Section S14 were renumbered correctly.
- Sources confirm what they are cited for: Kahan (median 25, IQR 15-44), Barker (45%), Westgate & Burchett (10 subjects, 4 occasions), Leyrat (12 analysis approaches compared), Hemming (recommendations), Li & Redden (ICC 0.001-0.05, mean cluster size 10-150, K − 2 df), Zeger (constant 16√3/(15π)), Bie (120 quadrature points; within- and between-cluster SEs).

## Deviations and limits

- Pepe & Anderson 1994 and Schildcrout & Heagerty 2005 were not accessible; only abstracts and secondary descriptions were read, so both are UNVERIFIABLE.
- Morel, Ford & Westgate and Gosho were read as abstracts only.
- Bradley was not searched (known to be unreadable).
- Numerical values owned by other tracks were taken from the registry unless they decided a support call; those were recomputed from server tables.
- References whose citing sentences did not change ([1], [2], [22]-[25], [27], [29], [30]) were not re-read.
