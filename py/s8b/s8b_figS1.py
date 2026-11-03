"""s8b_figS1.py -- S8b sub-task, plan item F02: Figure S1 of the supplement, type I error at every design point of
scenarios M6-M7 for the strategies of Figure 1 (Table 4 versions: fewclust corrections, t reference with N - p df for the
Wald tests, chi-squared reference for the likelihood-ratio test), with the exact-interval flags of Table 4 (D-63).
Adapted from the S8 reviewer's prototype ~/MMM_rev_S8/Figure/fig_typeI_m67.py (copied, not moved): same layout, no
panel titles (D-31), style tools/s6_figure.mplstyle, colours and line weights of Figure 1 (py/s7f/s7f_figures.py).
Inputs (read-only): tables/V5_S7f_rule_dp.csv (flags, as Table 4), tables/V5_S7f_cp_S7d_dp_rates.csv (rates and usable n
of the GEE strategies and of the marginalized model with the check, ref t), tables/V5_S7f_cp_S7b_rule_dp.csv (likelihood-
ratio test, ref n), tables/V5_S7f_numbers.csv (scenario parameters x:m6_gby ... x:m7_nat, x:min_obs_dp, Table 4 counts).
Outputs: reports/morris/S8b_figures/figS1.png (300 dpi) and .pdf, reports/morris/S8b_figures/figS1_data.csv (every
plotted value), tables/V5_S8b_sub_figS1_keys.csv (keys figS1|...), tables/V5_S8b_sub_figS1_checks.csv.
Run from ~/MMM: python3 py/s8b/s8b_figS1.py
"""
import os, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

os.chdir(os.path.expanduser("~/MMM"))
OUT = "reports/morris/S8b_figures"; os.makedirs(OUT, exist_ok=True)
plt.style.use("tools/s6_figure.mplstyle")
REG = pd.read_csv("tables/V5_S7f_numbers.csv", dtype=str, keep_default_na=False).set_index("id")
def rtext(i): return REG.loc[i, "text"]
def rval(i): return float(REG.loc[i, "value"])
MIN_OBS = int(rval("x:min_obs_dp"))
FOCAL = {"mmm": "#1b4f9c", "gee_fg": "#c2531b", "glmer_lrt": "#6a3d9a"}     # py/s7f/s7f_figures.py
GREYS = {"gee_md": "#8c8c8c"}
# (strategy, legend label, colour, line width, line style, zorder): the five series of Figure 1, panel c
STRAT = [("gee_ar1_fg", "GEE AR(1), Fay-Graubard, t", FOCAL["gee_fg"], 1.5, "-", 5),
         ("mmm_cc", "Marginalized model with check, t", FOCAL["mmm"], 1.5, "-", 4),
         ("gee_exch_fg", "GEE exchangeable, Fay-Graubard, t", FOCAL["gee_fg"], 0.8, "--", 3),
         ("gee_exch_md", "GEE exchangeable, Mancl-DeRouen, t", GREYS["gee_md"], 0.8, "-", 2),
         ("glmer_lrt", "GLMM, likelihood ratio", FOCAL["glmer_lrt"], 0.8, "-", 2)]
ROWS = [("GenBinaryY", "M6", f"GenBinaryY\nlatent-SD slope {rtext('x:m6_gby')}"), ("GenBinaryY", "M7", f"GenBinaryY\nlatent-SD slope {rtext('x:m7_gby')}"),
        ("GLMM-native", "M6", f"GLMM-native\nslope variance {rtext('x:m6_nat')}"), ("GLMM-native", "M7", f"GLMM-native\nslope variance {rtext('x:m7_nat')}")]
KP = ["dgp", "scen", "N", "ntime"]
RD = pd.read_csv("tables/V5_S7f_rule_dp.csv"); RD = RD[RD.scen.isin(["M6", "M7"])]
A = pd.read_csv("tables/V5_S7f_cp_S7d_dp_rates.csv"); A = A[(A.ref == "t") & A.scen.isin(["M6", "M7"])].copy()
A.loc[A.strategy == "study_mmm", "strategy"] = "mmm_cc"
C = pd.read_csv("tables/V5_S7f_cp_S7b_rule_dp.csv"); C = C[(C.ref == "n") & (C.method == "glmer_lrt") & C.scen.isin(["M6", "M7"])].copy()
C["strategy"] = "glmer_lrt"; C["kk"] = np.where(C.n > 0, np.round(C.n * C.rate.fillna(0) / 100), 0)
cols = KP + ["strategy", "n", "nd", "kk", "rate", "assess", "hi"]
D = pd.concat([A[A.strategy.isin([s[0] for s in STRAT])][cols], C[cols]], ignore_index=True)
assert D.loc[D.kk.isna(), "n"].eq(0).all(); D["kk"] = D.kk.fillna(0)
D = D.merge(RD[KP + ["strategy", "hi"]].rename(columns={"hi": "hi_t4"}), on=KP + ["strategy"], how="left", validate="one_to_one")
CHK = []
def chk(name, ok, detail=""):
    CHK.append((name, bool(ok), detail)); print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
chk("every plotted point has a Table 4 flag", D.hi_t4.notna().all(), f"{len(D)} rows")
chk("per-point flags equal V5_S7f_rule_dp.csv", (D.hi.astype(bool) == D.hi_t4.astype(bool)).all())
chk("assessed = at least x:min_obs_dp usable replicates", ((D.n >= MIN_OBS) == D.assess.astype(bool)).all())
chk("384 points per strategy", (D.groupby("strategy").size() == 384).all())
D["hi"] = D.hi_t4.astype(bool); D["assess"] = D.assess.astype(bool)
TS = sorted(D.ntime.unique()); Ns = sorted(D.N.unique()); gap = 3
xpos = {(t, N): i * (len(Ns) + gap) + j for i, t in enumerate(TS) for j, N in enumerate(Ns)}
D["x"] = [xpos[(t, N)] for t, N in zip(D.ntime, D.N)]
YMAX = 15.0
fig, axes = plt.subplots(4, 1, figsize=(7.2, 8.2), sharex=True, gridspec_kw=dict(hspace=0.12))
fig.subplots_adjust(left=0.13, right=0.98, top=0.98, bottom=0.10)
REC = []
for ax, (g, sc, lab) in zip(axes, ROWS):
    ax.axhspan(2.5, 7.5, color="#e9e9e9", lw=0, zorder=0); ax.axhline(5, color="#555555", lw=0.6, zorder=1)
    for s, slab, col, lw, ls, z in STRAT[::-1]:
        d = D[(D.dgp == g) & (D.scen == sc) & (D.strategy == s)].sort_values("x")
        y = np.where(d.assess, d.rate, np.nan)
        for t in TS:
            m = (d.ntime == t).values
            ax.plot(d.x.values[m], np.minimum(y[m], YMAX), color=col, lw=lw, ls=ls, zorder=z)
            yb = y[m]; fin = np.isfinite(yb); prv = np.r_[False, fin[:-1]]; nxt = np.r_[fin[1:], False]
            iso = fin & ~prv & ~nxt  # S8b review N283: a drawn point between undrawn neighbours gets a marker
            if iso.any(): ax.scatter(d.x.values[m][iso], np.minimum(yb[iso], YMAX), s=7, color=col, zorder=z)
        h = d[d.hi]
        ax.scatter(h.x, np.minimum(h.rate, YMAX), s=16, facecolors="none", edgecolors="#b2182b", linewidths=0.9, zorder=7)
        for r_, yy in zip(d.itertuples(), y):
            REC.append((g, sc, s, int(r_.ntime), int(r_.N), int(r_.n), int(r_.nd), int(r_.kk), yy, bool(r_.hi)))
    for i in range(1, len(TS)):
        ax.axvline(i * (len(Ns) + gap) - (gap + 1) / 2, color="#c8c8c8", lw=0.5, zorder=0)
    ax.set_ylim(0, YMAX); ax.set_yticks([0, 2.5, 5, 7.5, 10, 15]); ax.set_yticklabels(["0", "2.5", "5", "7.5", "10", "15"])
    ax.set_ylabel(f"{lab}\nType I error (%)")
centres = [i * (len(Ns) + gap) + (len(Ns) - 1) / 2 for i in range(len(TS))]
axes[-1].set_xticks(centres); axes[-1].set_xticklabels([f"T = {t}" for t in TS]); axes[-1].tick_params(axis="x", length=0)
axes[-1].set_xlabel(f"Occasions T (blocks); within each block the number of clusters N rises from {Ns[0]} to {Ns[-1]} ({len(Ns)} values)")
handles = [mpl.lines.Line2D([], [], color=c, lw=lw, ls=ls, label=l) for _, l, c, lw, ls, _ in STRAT]
handles.append(mpl.lines.Line2D([], [], ls="none", marker="o", mfc="none", mec="#b2182b", label="clearly above 7.5% (exact 95% interval)"))
fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.55, 0.0), fontsize=6.5)
assert not any(ax.get_title() for ax in axes) and fig._suptitle is None   # D-31: no titles
fig.savefig(f"{OUT}/figS1.png", dpi=300); fig.savefig(f"{OUT}/figS1.pdf")
S = pd.DataFrame(REC, columns=["dgp", "scen", "strategy", "ntime", "N", "n", "nd", "kk", "rate_plotted", "clearly_above"])
S.to_csv(f"{OUT}/figS1_data.csv", index=False)
## ---- keys and checks ------------------------------------------------------------------------------------------------
K = []
def put(key, value, kind, how, mcse="", source="tables/V5_S7f_rule_dp.csv"): K.append((key, str(int(value)) if kind == "int" else repr(float(value)), mcse, kind, source, how))
SH = {"GenBinaryY": "GB", "GLMM-native": "NAT"}
for s, *_ in STRAT:
    tot = int(S[S.strategy == s].clearly_above.sum())
    t4 = int(float(REG.loc[f"f:t4|{s}|hi_m67", "value"]))
    chk(f"circled points of {s} equal Table 4 f:t4|{s}|hi_m67", tot == t4, f"{tot} vs {t4}")
    put(f"figS1|{s}|n_hi", tot, "int", "circled points over the four panels (= f:t4|<strategy>|hi_m67)")
    for g, sc, _ in ROWS:
        put(f"figS1|{s}|{SH[g]}|{sc}|n_hi", int(S[(S.strategy == s) & (S.dgp == g) & (S.scen == sc)].clearly_above.sum()), "int", "circled points in the panel")
drawn = S.rate_plotted.notna()
put("figS1|n_points_panel", int(len(S) / len(STRAT) / len(ROWS)), "int", "design points per panel")
put("figS1|n_hidden", int((~drawn).sum()), "int", f"strategy x point pairs not drawn (fewer than {MIN_OBS} usable replicates)", source="tables/V5_S7f_cp_S7d_dp_rates.csv")
for s, *_ in STRAT:
    put(f"figS1|{s}|n_hidden", int((~drawn & (S.strategy == s)).sum()), "int", "points of the strategy not drawn", source="tables/V5_S7f_cp_S7d_dp_rates.csv")
put("figS1|n_clipped", int((S.rate_plotted > YMAX).sum()), "int", "drawn values above the axis limit (drawn at the limit)", source="tables/V5_S7f_cp_S7d_dp_rates.csv")
put("figS1|ymax", YMAX, "plain", "upper limit of the vertical axis (%)", source="py/s8b/s8b_figS1.py")
put("figS1|max_rate", float(S.rate_plotted.max()), "sig2", "largest drawn rate (%)", source="tables/V5_S7f_cp_S7d_dp_rates.csv")
mc_max = float(np.nanmax(np.where(drawn, 100 * np.sqrt(S.rate_plotted / 100 * (1 - S.rate_plotted / 100) / S.n), np.nan)))
put("figS1|mcse_max", mc_max, "sig2", "largest binomial MCSE of a drawn rate, 100 sqrt(p(1-p)/n) (pp)", source="tables/V5_S7f_cp_S7d_dp_rates.csv")
pd.DataFrame(K, columns=["key", "value", "mcse", "kind", "source", "how"]).to_csv("tables/V5_S8b_sub_figS1_keys.csv", index=False)
pd.DataFrame(CHK, columns=["check", "pass", "detail"]).to_csv("tables/V5_S8b_sub_figS1_checks.csv", index=False)
print(S.groupby(["dgp", "scen", "strategy"]).agg(drawn=("rate_plotted", lambda v: int(v.notna().sum())), max_rate=("rate_plotted", "max"), n_hi=("clearly_above", "sum")).to_string())
assert all(c[1] for c in CHK), "a check failed"
print("figS1 written; keys", len(K))
