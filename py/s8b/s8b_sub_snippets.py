"""s8b_sub_snippets.py -- S8b sub-task: markdown snippets for the supplement and the main text (Table S36, MCSE caption
sentences, Table S37 and the S11 paragraph, rows for Table S26 and the S6 sentence, Figure S1 caption). Every number is a
registry token `r n("f:s9|<key>")` (keys of tables/V5_S8b_sub_keys.csv) or an existing registry id; nothing is typed.
Sentences whose wording depends on a result (ranking in S11, hidden points in Figure S1) are chosen from the keys.
A second part repeats each snippet with the tokens resolved (preview only: 3 significant digits; the registry applies
its own rounding). Output: reports/morris/S8b_logs/sub/snippets.md. Run from ~/MMM after R/s8b/sub/s8b_sub_keys.R.
"""
import os, re, pandas as pd
os.chdir(os.path.expanduser("~/MMM"))
K = pd.read_csv("tables/V5_S8b_sub_keys.csv", dtype=str, keep_default_na=False).set_index("key")
R = pd.read_csv("tables/V5_S7f_numbers.csv", dtype=str, keep_default_na=False).set_index("id")
PFX = "f:s9|"
def k(key):
    assert key in K.index, f"key not found: {key}"
    return f'`r n("{PFX}{key}")`'
def r(i):
    assert i in R.index, f"registry id not found: {i}"
    return f'`r n("{i}")`'
def kv(key): return K.loc[key, "value"]
def ki(key): return int(float(kv(key)))
LAB = {"gee_fg": "GEE, Fay-Graubard", "gee_md": "GEE, Mancl-DeRouen", "mmm": "Marginalized model", "glmer": "GLMM, Wald",
       "shrink": "Shrinkage GLM", "gee": "GEE, uncorrected", "gee_kc": "GEE, Kauermann-Carroll"}
T4L = {"gee_ar1_fg": "GEE AR(1), Fay-Graubard", "mmm_cc": "Marginalized model with check", "gee_ind_fg": "GEE independence, Fay-Graubard",
       "gee_exch_fg": "GEE exchangeable, Fay-Graubard", "gee_ar1_md": "GEE AR(1), Mancl-DeRouen", "gee_exch_md": "GEE exchangeable, Mancl-DeRouen",
       "gee_ar1_kc": "GEE AR(1), Kauermann-Carroll", "gee_exch_kc": "GEE exchangeable, Kauermann-Carroll"}
MECH = ["GB", "NAT", "AR1"]
S = []   # (title, note, markdown)
# ---- 1. Table S36 -----------------------------------------------------------------------------------------------------
cap = (f"**Table S36.** Bias in empirical standard errors (Bias/SE), model-based over empirical standard error minus one (SE ratio, %) "
       f"and coverage of the 95% confidence interval against each method's own estimand (%), under the alternative, on each method's usable replicates, "
       f"main grid ({r('x:ndp_grid')} design points) and M1c ({r('x:ndp_m1c')} design points), marginalized model with the convergence check. "
       f"Each value is the mean over the design points with at least {r('x:min_obs_dp')} usable replicates: at least {k('s36|nused_min|main')} of the "
       f"{k('s36|ndp|main')} design points of a mechanism in the main grid, all {k('s36|ndp|M1c')} in M1c. "
       f"Not available (n.a.) for the GLMM on GenBinaryY, which has no marginal estimand. "
       f"Monte Carlo standard errors are at most {k('s36|mcse_max|rb')} for Bias/SE, {k('s36|mcse_max|relerr')} points for the SE ratio "
       f"and {k('s36|mcse_max|cov')} points for coverage. These are upper bounds: the mean over design points of the per-point Monte Carlo standard errors.")
hdr = ["Method", "Grid"] + [f"Bias/SE {m}" for m in MECH] + [f"SE ratio {m}" for m in MECH] + [f"Coverage {m}" for m in MECH]
rows = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
for me in ["gee_fg", "gee_md", "mmm", "glmer", "shrink", "gee", "gee_kc"]:
    for j, (g, gl) in enumerate([("main", "Main"), ("M1c", "M1c")]):
        c = []
        for ms in ["rb", "relerr", "cov"]:
            for mc in MECH:
                key = f"s36|{me}|{mc}|{g}|{ms}"
                c.append("n.a." if kv(key) == "n.a." else k(key))
        rows.append("| " + " | ".join([LAB[me] if j == 0 else "", gl] + c) + " |")
S.append(("Table S36 (DEC-T4), new table at the end of the supplement",
          "Source tables/variant_cv/V4_perf_primary.csv via R/s8b/sub/s8b_table_s36.R. The text of Section 2.5 can point here instead of the root mean squared error (P43).",
          cap + "\n\n" + "\n".join(rows)))
# ---- 2. MCSE sentences (DEC-T5) -----------------------------------------------------------------------------------------
B = "These are upper bounds: the mean over design points of the binomial standard error."
CAP = [
 ("Table 4 (main text)", f"Monte Carlo standard errors are at most {k('cap|t4|pu')} points for the usable share, {k('cap|t4|upow')} for power and {k('cap|t4|pu_small')} for the usable share at T ≤ 3 and N ≤ 25. {B} Counts of design points rest on the exact interval of each design point."),
 ("Figure 2 (main text)", f"The Monte Carlo standard error of a rate of 5% is {r('x:mcse_5pct')} points with {r('x:nsim')} usable replicates and at most {k('cap|fig2|mcse_max')} points at the drawn design point with the fewest."),
 ("Table S3", f"Largest Monte Carlo standard error: {k('cap|S3|rate')} points for the pooled rate, {k('cap|S3|worst')} for the worst case and {k('cap|S3|mad')} for the MAD."),
 ("Table S4", f"Largest Monte Carlo standard error: {k('cap|S4|mad')} points."),
 ("Table S5", f"Largest Monte Carlo standard error: {k('cap|S5|share')} points."),
 ("Table S6", f"Largest Monte Carlo standard error: {k('cap|S6|rate')} points for the rates and {k('cap|S6|shift')} for the shift."),
 ("Table S7", f"Largest Monte Carlo standard error: {k('cap|S7|usable')} points."),
 ("Table S8", f"Largest Monte Carlo standard error: {k('cap|S8|usable')} points for the usable share and {k('cap|S8|uncond')} for unconditional validity."),
 ("Table S9", f"Largest Monte Carlo standard error: {k('cap|S9|mad')} points."),
 ("Table S10", f"Largest Monte Carlo standard error: {k('cap|S10|raw')} points for raw power and {k('cap|S10|sa')} for size-adjusted power."),
 ("Table S11", f"Largest Monte Carlo standard error: {k('cap|S11|rate')} points for the rates and {k('cap|S11|shift')} for the differences."),
 ("Table S12", f"Largest Monte Carlo standard error: {k('cap|S12|cov')} points."),
 ("Table S14", f"Monte Carlo standard errors are at most {k('cap|S14|maxrate')} points for the largest rate (binomial, at its design point) and {k('cap|S14|pu')} for the usable share (upper bound: mean over design points of the binomial standard error)."),
 ("Table S15", f"Binomial Monte Carlo standard errors over the {r('x:lus_nsub')} subsamples are at most {k('cap|S15|share')} points."),
 ("Table S16", f"Largest Monte Carlo standard error: {k('cap|S16|rate')} points for the pooled rate, {k('cap|S16|worst')} for the worst case and {k('cap|S16|mad')} for the MAD."),
 ("Table S17", f"Largest Monte Carlo standard error: {k('cap|S17|mad')} points."),
 ("Table S18", f"Largest Monte Carlo standard error: {k('cap|S18|share')} points."),
 ("Table S19", f"Largest Monte Carlo standard error: {k('cap|S19|rate')} points for the rates and {k('cap|S19|shift')} for the shift."),
 ("Table S20", f"Largest Monte Carlo standard error: {k('cap|S20|usable')} points."),
 ("Table S21", f"Largest Monte Carlo standard error: {k('cap|S21|usable')} points for the usable share and {k('cap|S21|uncond')} for unconditional validity."),
 ("Table S22", f"Largest Monte Carlo standard error: {k('cap|S22|mad')} points."),
 ("Table S23", f"Largest Monte Carlo standard error: {k('cap|S23|raw')} points for raw power and {k('cap|S23|sa')} for size-adjusted power."),
 ("Table S24", f"Largest Monte Carlo standard error: {k('cap|S24|rate')} points for the rates and {k('cap|S24|shift')} for the differences."),
 ("Table S25", f"Largest Monte Carlo standard error: {k('cap|S25|cov')} points."),
 ("Table S26", f"The Monte Carlo standard error of the usable share is at most {k('cap|S26|pu')} points (upper bound: mean over design points of the binomial standard error)."),
 ("Table S27", f"Binomial Monte Carlo standard errors are at most {k('cap|S27|finite')} points for finite intervals and {k('cap|S27|gee')} for GEE; they are approximate, because the datasets of a row share replicate seeds across numbers of clusters."),
 ("Table S28", f"Monte Carlo standard errors are at most {k('cap|S28|type1')} points for the type I error, {k('cap|S28|usable')} for the usable share, {k('cap|S28|power')} for power, {k('cap|S28|biasse')} for bias in empirical standard errors and {k('cap|S28|cov')} points for coverage. These are upper bounds: the mean over design points of the per-point standard error."),
 ("Table S30", f"Largest Monte Carlo standard error: {k('cap|S30|usable')} points for the usable share, {k('cap|S30|bias')} for the bias and {k('cap|S30|seratio')} for the SE ratio."),
 ("Table S31", f"Block-bootstrap Monte Carlo standard errors are at most {k('cap|S31|rate')} points for the type I error, {k('cap|S31|mad')} for the departure from 5%, {k('cap|S31|pu')} for the usable share and {k('cap|S31|upow')} for power."),
 ("Table S32", f"The binomial Monte Carlo standard error of the highest rate is at most {k('cap|S32|maxrate')} points."),
 ("Table S33", f"The Monte Carlo standard error of the mean usable share is at most {k('cap|S33|pu')} points (upper bound: mean over design points of the binomial standard error)."),
 ("Table S34", f"Monte Carlo standard errors are at most {k('cap|S34|type1')} points for the mean type I error and {k('cap|S34|pu')} for the usable share. {B}"),
 ("Table S35", f"Monte Carlo standard errors, with the empirical standard error treated as known, are at most {k('cap|S35|rb')} for the mean bias in standard errors, {k('cap|S35|rbmax')} for the largest, {k('cap|S35|cov')} points for the mean coverage, {k('cap|S35|covmin')} for the smallest coverage and {k('cap|S35|sh')} for the share of fits in which the check intervened. For the SE ratio the typical per-point value is {k('cap|S35|ser')} points. The columns by intervention of the check have no Monte Carlo standard error."),
]
S.append(("MCSE caption sentences (DEC-T5)", "One sentence (or two) to append to each caption. Kind of each value in column how of tables/V5_S8b_sub_keys.csv: registry maximum (S3-S12, S16-S25, S30), upper bound for means over design points (Table 4, S14 usable, S26, S28, S33, S34, S35 means), binomial (Figure 2, S14 largest rate, S15, S27, S32, S35 smallest coverage), block bootstrap (S31), delta-method approximation (S35 bias and SE ratio).",
          "\n\n".join(f"- {t}: {s}" for t, s in CAP)))
# ---- 3. Table S37 and the S11 paragraph ------------------------------------------------------------------------------
st = [s for s in T4L if f"s37|{s}|hi_obs" in K.index]
cap37 = (f"**Table S37.** Stability of the count of design points clearly above 7.5% over the {r('f:x|ndp_all5')} design points of the four mechanisms, "
         f"strategies of Table 4 with an interval. Observed: the count in the study (Table 4). Bootstrap: {k('s37|B')} block-bootstrap resamples of the replicates, "
         f"the same resample for every design point of a model; share of resamples with no design point clearly above 7.5% (%), mean and 90th percentile of the count. "
         f"Plug-in: design points treated as independent, with true rate equal to the observed rate (raw) or to the rate shrunk by empirical Bayes within strategy and mechanism (shrunk); "
         f"probability of no design point clearly above 7.5% (%) and expected count.")
h = ["Strategy", "Observed", "Bootstrap, none (%)", "Bootstrap, mean", "Bootstrap, 90th percentile", "Raw plug-in, none (%)", "Raw plug-in, expected", "Shrunk plug-in, none (%)", "Shrunk plug-in, expected"]
rw = ["| " + " | ".join(h) + " |", "|" + "---|" * len(h)]
for s in st:
    rw.append("| " + " | ".join([T4L[s], k(f"s37|{s}|hi_obs"), k(f"s37|{s}|p_pass_boot|disp"), k(f"s37|{s}|mean_flags_boot"), k(f"s37|{s}|q90_flags_boot"),
                                 k(f"s37|{s}|plugin_p_pass|disp"), k(f"s37|{s}|plugin_exp_flags"), k(f"s37|{s}|eb_p_pass|disp"), k(f"s37|{s}|eb_exp_flags")]) + " |")
p1 = (f"The zero count of the two strategies of the rule is uncertain in a replication of the study. In the block bootstrap, GEE with AR(1) working correlation and "
      f"Fay-Graubard had at least one design point clearly above 7.5% in {k('s37|fg_pfail_boot')}% of resamples, and the marginalized model with the check in {k('s37|mmm_pfail_boot')}% (Table S37).")
if ki("s37|rank_same|table4") == 1 and ki("s37|rule_top2|table4") == 1 and ki("s37|rule_top2|file") == 1:
    p2 = "The ranking by expected count is the same in all three calculations, and in each of them the two strategies of the rule have the fewest expected design points clearly above 7.5%."
elif ki("s37|rule_top2|table4") == 1:
    p2 = "In each of the three calculations the two strategies of the rule have the fewest expected design points clearly above 7.5%, although the order of the other strategies changes."
else:
    p2 = "The two strategies of the rule do not have the fewest expected design points clearly above 7.5% in every calculation (Table S37)."
S.append(("Table S37 and paragraph for Section S11 (X02)", "Source tables/V5_S8b_rule_stability.csv via R/s8b/sub/s8b_table_s37.R; the calculations are those of R/s8b/provenance/lit_pow_rule_stability_v2.R.",
          cap37 + "\n\n" + "\n".join(rw) + "\n\nParagraph for S11, after Table S32:\n\n" + p1 + " " + p2))
# ---- 4. Table S26 rows and the S6 sentence ------------------------------------------------------------------------------
def rowA(s, lab):
    c = [k(f"mdkc|{s}|{d}|{f}") for d in ["N - p"] for f in ["n_hi", "n_lo", "n_lo_N10"]] + [k(f"mdkc|{s}|{d}|{f}") for d in ["N - 2"] for f in ["n_hi", "n_lo", "n_lo_N10"]]
    return "| " + " | ".join([lab] + c + [k(f"mdkc|{s}|pu")]) + " |"
def rowB(s, lab):
    c = [k(f"mdkc|{s}|N - p|{f}") for f in ["n_hi", "n_lo", "n_lo_N10"]] + ["n.a."] * 3
    return "| " + " | ".join([lab] + c + [k(f"mdkc|{s}|pu")]) + " |"
labs = [("gee_ar1_avgvar", "GEE AR(1), average of Mancl-DeRouen and Kauermann-Carroll"), ("gee_exch_avgvar", "GEE exchangeable, average of Mancl-DeRouen and Kauermann-Carroll")]
s26 = ("Version A (recommended; N - 2 computed with the same code, gate passed):\n\n" + "\n".join(rowA(s, l) for s, l in labs) +
       "\n\nCaption sentence for version A: The rows for the average of the Mancl-DeRouen and Kauermann-Carroll variances use the same usability rule, denominators and exact intervals, "
       "recomputed from the stored per-replicate results.\n\n"
       "Version B (as first instructed; N - 2 columns not filled):\n\n" + "\n".join(rowB(s, l) for s, l in labs) +
       "\n\nCaption sentence for version B: For the average of the Mancl-DeRouen and Kauermann-Carroll variances the N - 2 columns were not computed (n.a.).\n\n"
       "Sentence for the text of S6 (or Section 4.5 of the main text):\n\n"
       f"Averaging the Mancl-DeRouen and Kauermann-Carroll variances, as studied by Ford and Westgate [14], did not give a strategy that was never clearly above 7.5%. "
       f"With N - p degrees of freedom, GEE with AR(1) working correlation was clearly above 7.5% at {k('mdkc|gee_ar1_avgvar|hi_all')} of the {r('f:x|ndp_all5')} design points "
       f"of the four mechanisms, and GEE with exchangeable working correlation at {k('mdkc|gee_exch_avgvar|hi_all')} (Table S26 gives the three study mechanisms).")
S.append(("Rows for Table S26 and sentence for S6 (X03)", "Rows go after the row of GEE exchangeable, Kauermann-Carroll. Source R/s8b/sub/s8b_mdkc_s26.R (1,896 points, N - p and N - 2); the text sentence uses tables/V5_S8b_mdkc_avg.csv (2,376 points).", s26))
# ---- 5. Figure S1 caption ----------------------------------------------------------------------------------------------
others0 = all(ki(f"figS1|{s}|n_hidden") == 0 for s in ["gee_ar1_fg", "mmm_cc", "gee_exch_fg", "gee_exch_md"])
hid = f"({k('figS1|glmer_lrt|n_hidden')} design points of the likelihood-ratio test{', none of the other strategies' if others0 else ''})"
nat0 = all(ki(f"figS1|{s}|NAT|{sc}|n_hi") == 0 for s in ["gee_ar1_fg", "mmm_cc", "gee_exch_fg", "gee_exch_md", "glmer_lrt"] for sc in ["M6", "M7"])
clip = "" if ki("figS1|n_clipped") == 0 else f" {k('figS1|n_clipped')} values above the axis are drawn at its limit."
figcap = ("![](figures/figS1.png){width=100%}\n\n"
          f"**Figure S1.** Type I error at every design point of scenarios M6-M7 ({r('x:n_dp_ns_scen')} design points each) for the strategies of Figure 1, with the corrections from `fewclust` as in Table 4: "
          f"Wald tests with a t reference (N - p degrees of freedom) and the likelihood-ratio test with a $\\chi^2$ reference. Rows: GenBinaryY with latent-SD slope {r('x:m6_gby')} and {r('x:m7_gby')}; "
          f"GLMM-native with random-slope variance {r('x:m6_nat')} and {r('x:m7_nat')}. Design points are ordered by number of occasions T (blocks) and number of clusters N. "
          f"Grey band: 2.5% to 7.5%; line: nominal 5%. Design points with fewer than {r('x:min_obs_dp')} usable replicates are not drawn {hid}. "
          f"Circles: design points clearly above 7.5% (lower limit of the exact 95% interval above 7.5%), as in Table 4: {k('figS1|gee_exch_md|n_hi')} for GEE exchangeable with Mancl-DeRouen, "
          f"{k('figS1|gee_exch_fg|n_hi')} for GEE exchangeable with Fay-Graubard{', all on GenBinaryY' if nat0 else ''}; {k('figS1|gee_ar1_fg|n_hi')} for GEE AR(1) with Fay-Graubard, "
          f"{k('figS1|mmm_cc|n_hi')} for the marginalized model with the check and {k('figS1|glmer_lrt|n_hi')} for the likelihood-ratio test.{clip} "
          f"The binomial Monte Carlo standard error of a drawn rate is at most {k('figS1|mcse_max')} points.")
S.append(("Figure S1 (F02)", "Figure file reports/morris/S8b_figures/figS1.png (copy to manuscript/v5/figures/figS1.png when the supplement is edited). Script py/s8b/s8b_figS1.py.", figcap))
# ---- write -----------------------------------------------------------------------------------------------------------
def fmt(v, kind):
    if kind == "int": return str(int(float(v)))
    if kind == "text": return v
    x = float(v); return f"{x:.3g}"
def resolve(t):
    def one(m):
        i = m.group(1)
        if i.startswith(PFX): kk = i[len(PFX):]; return "[" + fmt(K.loc[kk, "value"], K.loc[kk, "kind"]) + "]"
        return "[" + R.loc[i, "text"] + "]"
    return re.sub(r'`r n\("([^"]+)"\)`', one, t)
L = ["# S8b sub-task: snippets for the supplement and the main text", "",
     f"Tokens `r n(\"{PFX}<key>\")` refer to keys of tables/V5_S8b_sub_keys.csv (prefix {PFX} as requested); other tokens are existing registry ids. "
     "Part 2 repeats every snippet with values in square brackets, for checking only.", "", "## Part 1. Snippets with tokens", ""]
for t, n_, m in S: L += [f"### {t}", "", f"*{n_}*", "", m, ""]
L += ["## Part 2. Preview with values (3 significant digits, not the registry rounding)", ""]
for t, n_, m in S: L += [f"### {t}", "", resolve(m), ""]
txt = "\n".join(L)
assert "\u2014" not in txt and "\u2013" not in txt, "em or en dash found"
os.makedirs("reports/morris/S8b_logs/sub", exist_ok=True)
open("reports/morris/S8b_logs/sub/snippets.md", "w", encoding="utf-8").write(txt)
print("snippets written:", len(S), "blocks,", len(re.findall(r'`r n\(', txt)), "tokens")
