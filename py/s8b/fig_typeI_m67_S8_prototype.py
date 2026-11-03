"""fig_typeI_m67.py -- proposal for the next step (doubt 2 received by S8): redraw of
reports/morris/S6b_figures/S6b_fig_typeI_newscen.png (false positives in scenarios M6-M7) WITHOUT panel titles
(D-31: what each panel shows goes in the caption) and with the exact Clopper-Pearson flags of D-63, from
tables/V4_S6b_ns_dp.csv only. Style: tools/s6_figure.mplstyle; colours and line weights as in Figures 1-2 of v5.
Rows (top to bottom): GenBinaryY M6 (latent-SD slope 0.1), GenBinaryY M7 (0.25), GLMM-native M6 (random-slope
variance 0.01), GLMM-native M7 (0.04); scenario parameters from tables/V4_S7_extra.csv (x:m6_gby ... x:m7_nat).
Rate under H0 = 100 * kt / n (Wald test with t at N - p, column kt), n = usable replicates; points with n < 20
are not drawn and not flagged (as in Figure 2 and Table S14). Clearly above 7.5%: lower limit of the exact
two-sided 95% Clopper-Pearson interval above 7.5%.
Usage: python3 fig_typeI_m67.py <V4_S6b_ns_dp.csv> <s6_figure.mplstyle> <out_dir>
Outputs: <out_dir>/fig_typeI_m67.png/.pdf and fig_typeI_m67_data.csv (every plotted value, for a read-back check).
"""
import sys, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import beta

SRC, STYLE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
plt.style.use(STYLE)
MIN_OBS = 20
ROWS = [("GenBinaryY", "M6", "GenBinaryY\nlatent-SD slope 0.1"), ("GenBinaryY", "M7", "GenBinaryY\nlatent-SD slope 0.25"),
        ("GLMM-native", "M6", "GLMM-native\nslope variance 0.01"), ("GLMM-native", "M7", "GLMM-native\nslope variance 0.04")]
STRAT = [("mmm", "MMM with convergence check, t", "#1b4f9c", 1.5, "-", 4),
         ("gee_fg", "GEE, Fay-Graubard, t", "#c2531b", 1.0, "-", 3),
         ("mmm_nocc", "MMM without check, t", "#8c8c8c", 0.8, "-", 2)]
TS = [2, 3, 4, 6, 8, 10]
D = pd.read_csv(SRC)
D = D[D.arm == "H0"].copy()
D["rate"] = np.where(D.n > 0, 100 * D.kt / D.n.where(D.n > 0, 1), np.nan)
k, n = D.kt.values, D.n.values
lo = np.where(k > 0, beta.ppf(0.025, np.maximum(k, 1e-12), n - k + 1), 0.0)
D["cp_lo"] = 100 * np.where(n > 0, lo, np.nan)
D["assess"] = D.n >= MIN_OBS
D["hi"] = D.assess & (D.cp_lo > 7.5)
Ns = sorted(D.N.unique()); gap = 3
xpos = {(t, N): i * (len(Ns) + gap) + j for i, t in enumerate(TS) for j, N in enumerate(Ns)}
D["x"] = [xpos[(t, N)] for t, N in zip(D.ntime, D.N)]
fig, axes = plt.subplots(4, 1, figsize=(7.2, 8.2), sharex=True, gridspec_kw=dict(hspace=0.12))
fig.subplots_adjust(left=0.13, right=0.98, top=0.98, bottom=0.12)
REC = []
for ax, (g, sc, lab) in zip(axes, ROWS):
    ax.axhspan(2.5, 7.5, color="#e9e9e9", lw=0, zorder=0); ax.axhline(5, color="#555555", lw=0.6, zorder=1)
    for m, mlab, col, lw, ls, z in STRAT[::-1]:
        d = D[(D.dgp == g) & (D.scen == sc) & (D.method == m)].sort_values("x")
        y = np.where(d.assess, d.rate, np.nan)
        for t in TS:
            s = (d.ntime == t).values
            ax.plot(d.x.values[s], y[s], color=col, lw=lw, ls=ls, zorder=z)
        h = d[d.hi]
        ax.scatter(h.x, h.rate, s=16, facecolors="none", edgecolors="#b2182b", linewidths=0.9, zorder=6)
        for r_, yy in zip(d.itertuples(), y):
            REC.append((g, sc, m, int(r_.ntime), int(r_.N), int(r_.n), int(r_.kt), yy, bool(r_.hi)))
    for i in range(1, len(TS)):
        ax.axvline(i * (len(Ns) + gap) - (gap + 1) / 2, color="#c8c8c8", lw=0.5, zorder=0)
    ax.set_ylim(0, 15); ax.set_yticks([0, 2.5, 5, 7.5, 10, 15]); ax.set_yticklabels(["0", "2.5", "5", "7.5", "10", "15"])
    ax.set_ylabel(f"{lab}\nType I error (%)")
centres = [i * (len(Ns) + gap) + (len(Ns) - 1) / 2 for i in range(len(TS))]
axes[-1].set_xticks(centres); axes[-1].set_xticklabels([f"T = {t}" for t in TS])
axes[-1].tick_params(axis="x", length=0)
axes[-1].set_xlabel(f"Occasions T (blocks); within each block the number of clusters N rises from {Ns[0]} to {Ns[-1]} ({len(Ns)} values)")
handles = [mpl.lines.Line2D([], [], color=c, lw=lw, label=l) for _, l, c, lw, _, _ in STRAT]
handles.append(mpl.lines.Line2D([], [], ls="none", marker="o", mfc="none", mec="#b2182b", label="clearly above 7.5% (exact 95% interval)"))
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.55, 0.0), fontsize=6.5)
assert not any(ax.get_title() for ax in axes) and fig._suptitle is None   # D-31: no titles
fig.savefig(f"{OUT}/fig_typeI_m67.png"); fig.savefig(f"{OUT}/fig_typeI_m67.pdf")
pd.DataFrame(REC, columns=["dgp", "scen", "method", "ntime", "N", "n", "kt", "rate_plotted", "clearly_above"]).to_csv(f"{OUT}/fig_typeI_m67_data.csv", index=False)
S = pd.DataFrame(REC, columns=["dgp", "scen", "method", "ntime", "N", "n", "kt", "rate", "hi"])
print(S.groupby(["dgp", "scen", "method"]).agg(n_pts=("rate", "size"), n_drawn=("rate", lambda v: int(np.isfinite(v).sum())), max_rate=("rate", "max"), n_hi=("hi", "sum")).to_string())
