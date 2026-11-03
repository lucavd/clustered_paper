"""s7f_figures.py -- S7f: Figure 1 (v5) and nested-loop plot; copy of py/s6/s6_figures.py with panel c read from
tables/V5_S7f_fig1c.csv (2,376 design points of the four mechanisms). Usage: python s7f_figures.py <tables_dir> <out_dir> <style> <fig1c csv>
(original docstring follows)

Every plotted value is read from tables/V4_S6_values.csv, V4_S6_nestedloop.csv or
V4_S6_zip.csv and written back, as plotted, to reports/morris/S6_figures/*_data.csv
(re-read by R/s6/s6_check_figures.R, gate G3). Style: tools/s6_figure.mplstyle
(rcParams of figure-style, sizes 8/7/6). Monte Carlo intervals: estimate +/- 2 MCSE
(block-bootstrap SE, B = 2000). Usage: python s6_figures.py <tables_dir> <out_dir> <style>
"""
import sys, json, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

TAB, OUT, STYLE, FIG1C = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
plt.style.use(STYLE)
LAB = {"gee_fg": "FG", "gee_md": "MD", "mmm": "MMM", "glmer": "GLMM", "shrink": "Shrink",
       "gee": "GEE", "glmer_lrt": "LRT", "gee_kc": "KC"}
LONG = {"gee_fg": "GEE, Fay-Graubard", "gee_md": "GEE, Mancl-DeRouen", "mmm": "Marginalized model",
        "glmer": "GLMM, Wald", "shrink": "Shrinkage GLM", "gee": "GEE, uncorrected",
        "glmer_lrt": "GLMM, likelihood ratio", "gee_kc": "GEE, Kauermann-Carroll"}
ORD = ["gee_fg", "gee_md", "mmm", "glmer", "shrink", "gee", "glmer_lrt", "gee_kc"]
FOCAL = {"mmm": "#1b4f9c", "gee_fg": "#c2531b", "glmer_lrt": "#6a3d9a"}
GREYS = {"gee_md": "#8c8c8c", "glmer": "#6e6e6e", "shrink": "#a8a8a8", "gee": "#b5b5b5", "gee_kc": "#7d7d7d"}
COL = {**FOCAL, **GREYS}
MIN_OBS = 20   # nested-loop plot: design points with fewer usable replicates are not drawn
MECH = {"GB": "GenBinaryY", "NAT": "GLMM-native", "AR1": "GLMM-AR1"}
V = pd.read_csv(f"{TAB}/V4_S6_values.csv")
VM = V.set_index(["variant", "name"])
RECORD = []   # (figure, panel, series, x, est, lo, hi, source names)

def val(name, variant="main"):
    r = VM.loc[(variant, name)]
    return float(r.est), float(r.se)

def letter(ax, s):
    ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")

def bbox_check(fig):
    fig.canvas.draw(); r = fig.canvas.get_renderer()   # renderer at the figure dpi, after any savefig at 300 dpi
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes if ax.axison for s in ax.spines.values() if s.get_visible()]
    tick = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both")) for ax in fig.axes}
    ov = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    ov += [(t.get_text(), "spine") for t, bt in texts for s, bs in spines if bt.overlaps(bs) and t not in tick[s.axes]]
    fb = fig.bbox
    out = [t.get_text() for t, bt in texts if bt.x0 < fb.x0 - 1 or bt.y0 < fb.y0 - 1 or bt.x1 > fb.x1 + 1 or bt.y1 > fb.y1 + 1]
    return ov, out

def place_labels(ax, x, ys, labels, colors, mingap):
    """direct labels at the right end, staggered so that no two are closer than mingap (data units)."""
    o = np.argsort(ys); pos = np.array(ys, float)[o].copy()
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + mingap)
    shift = (np.array(ys, float)[o].mean() - pos.mean()); pos = pos + min(0, shift) * 0
    for k, idx in enumerate(o):
        ax.annotate(labels[idx], xy=(x, ys[idx]), xytext=(x * 1.0, pos[k]), textcoords="data", color=colors[idx],
                    fontsize=7, va="center", ha="left",
                    arrowprops=dict(arrowstyle="-", color=colors[idx], lw=0.5, shrinkA=0, shrinkB=0) if abs(pos[k] - ys[idx]) > 1e-9 else None)

# ------------------------------------------------------------------------------------------
def figure1(rule_lines, fname, switch_N=None):
    # v5 (D-61): every panel on the strategies of the decision rule, t reference, 2,376 design points of the four mechanisms
    fig = plt.figure(figsize=(7.2, 9.6))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.62, left=0.14, right=0.87, top=0.96, bottom=0.07)
    S5 = [("gee_ar1_fg", "FG, AR(1), t", FOCAL["gee_fg"], 1.8, "-"), ("mmm_cc", "MMM, t", FOCAL["mmm"], 1.8, "-"),
          ("gee_exch_fg", "FG, exch., t", FOCAL["gee_fg"], 0.9, "--"), ("gee_exch_md", "MD, exch., t", GREYS["gee_md"], 0.9, "-"),
          ("glmer_lrt", "LRT", FOCAL["glmer_lrt"], 1.0, "-")]
    FLIPC = "#2a7f62"
    def endlabels(ax, items, x0, x1, gap, log=False):
        items = sorted(items); v = [np.log10(e[0]) if log else e[0] for e in items]; pos = list(v)
        for i in range(1, len(pos)): pos[i] = max(pos[i], pos[i - 1] + gap)
        for (y, lab, col), p_ in zip(items, pos):
            yy = 10 ** p_ if log else p_
            ax.plot([x0, x1], [y, yy], color=col, lw=0.5, clip_on=False)
            ax.text(x1 + (x1 - x0) * 0.15, yy, lab, color=col, fontsize=7, va="center", ha="left", clip_on=False)
    # a (D-61, second review): design map. Reading of the GEE - MMM difference in mean absolute departure against the
    # one-point margin (block bootstrap, tables/V5_S7f_h2h_boot.csv); small text = usable results (%) of GEE / MMM
    HB = pd.read_csv("tables/V5_S7f_h2h_boot.csv")
    hb = lambda reg, q: HB[(HB.region == reg) & (HB.q == q)].iloc[0]
    ax = fig.add_subplot(gs[0, 0]); letter(ax, "a")
    NBS = ["10", "15-25", ">25"]; TMS = ["4+", "3", "2"]
    for i_, nb in enumerate(NBS):
        for j_, tm in enumerate(TMS):
            reg = f"map|{nb}|{tm}"; d_ = hb(reg, "mad_gee_minus_mmm"); cl = d_["class"]
            lean = "GEE" if d_["est"] < 0 else "MMM"
            # PI 2026-10-27: three labels only. A strategy is named only when the whole interval lies beyond the margin;
            # otherwise "GEE = MMM" (dotted frame where the interval crosses the margin)
            lab, key = {"GEE": ("GEE", "GEE"), "MMM": ("MMM", "MMM"), "similar": ("GEE = MMM", "sim")}.get(cl, ("GEE = MMM", "ab"))
            fc = {"GEE": "#fbeee6", "MMM": "#e6edf8", "sim": "#f2f2f2", "ab": "#f2f2f2"}[key]
            ec = {"GEE": FOCAL["gee_fg"], "MMM": FOCAL["mmm"], "sim": "#9a9a9a", "ab": "#9a9a9a"}[key]
            tc = {"GEE": FOCAL["gee_fg"], "MMM": FOCAL["mmm"], "sim": "#555555", "ab": "#555555"}[key]
            ax.add_patch(mpl.patches.Rectangle((i_ + 0.04, j_ + 0.04), 0.92, 0.92, fc=fc, ec=ec, lw=0.8, ls="-" if key in ("GEE", "MMM", "sim") else (0, (2, 1.5))))
            ax.text(i_ + 0.5, j_ + 0.6, lab, ha="center", va="center", fontsize=8.5 if key in ("GEE", "MMM") else 6.6, linespacing=1.1, fontweight="bold" if key in ("GEE", "MMM") else "normal", color=tc)
            pg, pm = hb(reg, "pu_gee")["est"], hb(reg, "pu_mmm")["est"]
            pct = lambda v: ">99" if 99.5 <= v < 100 else f"{v:.0f}"
            ax.text(i_ + 0.5, j_ + 0.22, f"usable {pct(pg)}/{pct(pm)}%", ha="center", va="center", fontsize=5.1, color="#444444")
            RECORD.append(("fig1", "a", reg, "mad_gee_minus_mmm", d_["est"], d_["est"] - 2 * d_["mcse"], d_["est"] + 2 * d_["mcse"], f"V5_S7f_h2h_boot:{reg}|mad_gee_minus_mmm|{cl}"))
            RECORD.append(("fig1", "a", reg, "pu_gee", pg, np.nan, np.nan, f"V5_S7f_h2h_boot:{reg}|pu_gee"))
            RECORD.append(("fig1", "a", reg, "pu_mmm", pm, np.nan, np.nan, f"V5_S7f_h2h_boot:{reg}|pu_mmm"))
    ax.add_patch(mpl.patches.Rectangle((0.0, 1.0), 2.0, 2.0, fc="none", ec="#2a7f62", lw=1.3, ls="--", zorder=5))
    ax.text(1.0, 3.06, "+ sign-flip test (p-value only)", ha="center", va="bottom", fontsize=6.5, color="#2a7f62")
    ax.set_xlim(0, 3); ax.set_ylim(0, 3.32); ax.set_aspect("auto")
    ax.set_xticks([0.5, 1.5, 2.5]); ax.set_xticklabels(["10", "15-25", "> 25"])
    ax.set_yticks([0.5, 1.5, 2.5]); ax.set_yticklabels(["4 or more", "3", "2"])
    ax.tick_params(length=0); [s_.set_visible(False) for s_ in ax.spines.values()]
    ax.set_xlabel("Clusters (e.g. patients)"); ax.set_ylabel("Occasions per cluster")
    # S8b (F01): no title above the panel (D-31); the question is in the caption
    # b: mean absolute departure by exact T
    BT = pd.read_csv("tables/V5_S7f_fig1_byT.csv")
    SB = pd.read_csv("tables/V5_S7f_h2h_series_boot.csv")
    ax = fig.add_subplot(gs[0, 1]); letter(ax, "b"); ends = []
    for key, lab, col, lw, ls in S5:
        d = BT[BT.strategy == key].sort_values("ntime")
        ax.plot(d.ntime, d.mad, color=col, lw=lw, ls=ls, marker="o", ms=3 if lw > 1.5 else 2, zorder=4 if lw > 1.5 else 2)
        if key in ("gee_ar1_fg", "mmm_cc"):   # block-bootstrap MCSE (s7f_h2h_boot.R), the two strategies of the rule only
            se_ = np.array([SB[(SB.series == f"T|{int(t_)}") & (SB.strategy == key)].mad_se.iloc[0] for t_ in d.ntime])
            ax.vlines(d.ntime, np.maximum(d.mad - 2 * se_, 0.3), d.mad + 2 * se_, color=col, lw=0.8, zorder=2)
        else:
            se_ = np.full(len(d), np.nan)
        for r, e_ in zip(d.itertuples(), se_):
            RECORD.append(("fig1", "b", key, int(r.ntime), r.mad, r.mad - 2 * e_, r.mad + 2 * e_, f"V5_S7f_fig1_byT:{key}|{int(r.ntime)}|mad"))
        ends.append((float(d.mad.iloc[-1]), lab, col))
    ax.set_yscale("log"); ax.set_yticks([0.5, 1, 2, 5, 10]); ax.set_yticklabels(["0.5", "1", "2", "5", "10"]); ax.set_ylim(0.4, 14)
    ax.set_xticks([2, 3, 4, 6, 8, 10]); ax.set_xlim(1.6, 10.4)
    endlabels(ax, ends, 10.15, 10.9, 0.085, log=True)
    ax.set_xlabel("Occasions per cluster"); ax.set_ylabel("Mean absolute departure\nfrom 5% (pp), lower = better")
    # c: share of design points clearly above 7.5% by N
    ax = fig.add_subplot(gs[1, 0]); letter(ax, "c")
    SF = pd.read_csv(FIG1C); ends = []; ymax = 0
    for key, lab, col, lw, ls in S5:
        d = SF[SF.strategy == key].sort_values("N"); yv = 100 * d.n_hi / d.n_dp; ymax = max(ymax, yv.max())
        ax.plot(d.N, yv, color=col, lw=lw, ls=ls, marker="o", ms=3 if key == "gee_ar1_fg" else 2, zorder=4 if key == "gee_ar1_fg" else 2)
        for n_, yy in zip(d.N, yv):
            RECORD.append(("fig1", "c", key, int(n_), yy, np.nan, np.nan, f"V5_S7f_fig1c:{key}|{int(n_)}"))
        ends.append((float(yv.iloc[-1]), lab, col))
    ax.set_xlim(5, 150); ax.set_xticks([20, 40, 60, 80, 100, 120, 140]); ax.set_ylim(-1, max(10, 5 * np.ceil((ymax + 2) / 5)))
    endlabels(ax, ends, 146, 151, 2.6)
    ax.set_xlabel("Clusters (e.g. patients)"); ax.set_ylabel("Design points with false positives\nclearly above 7.5% (%), lower = better")
    # d: mean absolute departure by N, the two strategies, T = 2 (dashed) and T >= 3 (solid)
    BN = pd.read_csv("tables/V5_S7f_fig1_byN.csv")
    ax = fig.add_subplot(gs[1, 1]); letter(ax, "d"); ends = []
    for key, lab, col in [("gee_ar1_fg", "FG, AR(1), t", FOCAL["gee_fg"]), ("mmm_cc", "MMM, t", FOCAL["mmm"])]:
        for tb, ls in [("2", "--"), ("3+", "-")]:
            d = BN[(BN.strategy == key) & (BN.Tb == tb)].sort_values("N")
            ax.plot(d.N, d.mad, color=col, lw=1.4, ls=ls, marker="o", ms=2.2, zorder=3)
            se_ = np.array([SB[(SB.series == f"NT|{int(n_)}|{tb}") & (SB.strategy == key)].mad_se.iloc[0] for n_ in d.N])
            ax.vlines(d.N, d.mad - 2 * se_, d.mad + 2 * se_, color=col, lw=0.7, zorder=2)
            for r, e_ in zip(d.itertuples(), se_):
                RECORD.append(("fig1", "d", f"{key}|{tb}", int(r.N), r.mad, r.mad - 2 * e_, r.mad + 2 * e_, f"V5_S7f_fig1_byN:{key}|{int(r.N)}|{tb}|mad"))
            ends.append((float(d.mad.iloc[-1]), ("FG" if key == "gee_ar1_fg" else "MMM") + (", T = 2" if tb == "2" else ", T ≥ 3"), col))
    ax.set_xlim(5, 150); ax.set_xticks([20, 40, 60, 80, 100, 120, 140]); ax.set_ylim(0, 5.3)
    endlabels(ax, ends, 146, 151, 0.36)
    ax.set_xlabel("Clusters (e.g. patients)"); ax.set_ylabel("Mean absolute departure\nfrom 5% (pp), lower = better")
    # e-f: usable results and power by T
    EF = pd.read_csv("tables/V5_S7f_fig1_ef.csv")
    S3 = [("gee_ar1_fg", "FG, AR(1), t", FOCAL["gee_fg"], "-"), ("mmm_cc", "MMM, t", FOCAL["mmm"], "-"), ("flip", "Sign-flip", FLIPC, "--")]
    for pos, col_, lab_y, let, ylim, gap in [(gs[2, 0], "pu", "Datasets with a usable\nresult (%), higher = better", "e", (65, 101), 2.2),
                                            (gs[2, 1], "upow", "Effect detected among\nattempted datasets (%)", "f", (0, 58), 3.2)]:
        ax = fig.add_subplot(pos); letter(ax, let); ends = []
        for key, lab, col, ls in S3:
            d = EF[EF.strategy == key].sort_values("ntime")
            ax.plot(d.ntime, d[col_], color=col, lw=1.6 if key != "flip" else 1.1, ls=ls, marker="o", ms=2.6, zorder=3)
            for r in d.itertuples():
                v_ = getattr(r, col_); s_ = getattr(r, col_ + "_se")
                RECORD.append(("fig1", let, key, int(r.ntime), v_, v_ - 2 * s_, v_ + 2 * s_, f"V5_S7f_fig1_ef:{key}|{int(r.ntime)}|{col_}"))
            ends.append((float(d[col_].iloc[-1]), lab, col))
        ax.set_xticks([2, 3, 4, 6, 8, 10]); ax.set_xlim(1.6, 10.4); ax.set_ylim(*ylim)
        if let == "e": ax.set_yticks([70, 80, 90, 100])
        endlabels(ax, ends, 10.15, 10.9, gap)
        ax.set_xlabel("Occasions per cluster"); ax.set_ylabel(lab_y)
    fig.savefig(f"{OUT}/{fname}.png"); fig.savefig(f"{OUT}/{fname}.pdf")
    return fig

# ------------------------------------------------------------------------------------------
def nested_loop(fname):
    NL = pd.read_csv(f"{TAB}/V4_S6_nestedloop.csv")
    SC = ["M1", "M2", "M3", "M4", "M5", "M1c"]
    fig, axes = plt.subplots(4, 1, figsize=(7.2, 8.4), sharex=True, gridspec_kw=dict(height_ratios=[1, 1, 1, 0.45], hspace=0.14))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.95, bottom=0.06)
    keys = NL[["dgp", "scen", "ntime", "N"]].drop_duplicates()
    order = []
    xpos = 0
    for sc in SC:
        k = keys[(keys.scen == sc) & (keys.dgp == "GenBinaryY")].sort_values(["ntime", "N"])
        for _, r in k.iterrows():
            order.append((sc, r.ntime, r.N, xpos)); xpos += 1
        xpos += 12
    X = pd.DataFrame(order, columns=["scen", "ntime", "N", "x"])
    ymax = 0
    for ax, (g, gl) in zip(axes[:3], [("GenBinaryY", "GenBinaryY"), ("GLMM-native", "GLMM-native"), ("GLMM-AR1", "GLMM-AR1")]):
        ax.axhspan(2.5, 7.5, color="#e9e9e9", lw=0, zorder=0); ax.axhline(5, color="#555555", lw=0.6, zorder=1)
        for m in ORD[::-1]:
            d = NL[(NL.dgp == g) & (NL.method == m)].merge(X, on=["scen", "ntime", "N"]).sort_values("x")
            y = np.where(d.n_obs >= MIN_OBS, d.rej, np.nan)
            foc = m in FOCAL
            for sc in SC:
                s = d.scen.values == sc
                ax.plot(d.x.values[s], y[s], color=COL[m], lw=0.9 if foc else 0.45, drawstyle="steps-mid", zorder=3 if foc else 2)
            for r, yy in zip(d.itertuples(), y):
                RECORD.append(("nestedloop", g, m, f"{r.scen}|{r.ntime}|{r.N}", yy, yy - 2 * r.mcse_rej if r.n_obs >= MIN_OBS else np.nan,
                               yy + 2 * r.mcse_rej if r.n_obs >= MIN_OBS else np.nan, f"V4_S6_nestedloop:{g}|{r.scen}|{r.N}|{r.ntime}|{m}"))
            ymax = max(ymax, np.nanmax(y))
        ax.set_ylabel(f"{gl}\nType I error (%)")
    for ax in axes[:3]:
        ax.set_ylim(0, 30); ax.set_yticks([0, 5, 10, 20, 30])
    clipped = int(((NL.rej > 30) & (NL.n_obs >= MIN_OBS)).sum()); hidden = int((NL.n_obs < MIN_OBS).sum())
    for sc in SC:
        xs = X[X.scen == sc].x
        axes[0].text(xs.mean(), 31.0, sc, ha="center", va="bottom", fontsize=7)
    # legend by direct labels in the first panel's whitespace
    lx = X[X.scen == "M1c"].x.max() + 8
    for i, m in enumerate(ORD):
        axes[0].text(lx, 28 - 2.6 * i, LAB[m], color=COL[m], fontsize=6.5, va="center", ha="left")
    # loop indicators
    ax = axes[3]
    lN = np.log(X.N) ; nN = (lN - np.log(10)) / (np.log(145) - np.log(10))
    for sc in SC:
        s = X.scen == sc
        ax.plot(X.x[s], 2.2 + nN[s], color="#333333", lw=0.6, drawstyle="steps-mid")
        ax.plot(X.x[s], 0.2 + (X.ntime[s] - 2) / 8, color="#333333", lw=0.6, drawstyle="steps-mid")
    ax.set_yticks([0.2, 1.2, 2.2, 3.2]); ax.set_yticklabels(["T = 2", "T = 10", "N = 10", "N = 145"])
    ax.set_ylim(-0.1, 3.5); ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
    ax.set_xlabel("Design points: scenario (outer loop), occasions T, clusters N (inner loop)")
    ax.set_xlim(-5, lx + 40)
    fig.savefig(f"{OUT}/{fname}.png"); fig.savefig(f"{OUT}/{fname}.pdf")
    return fig, (clipped, hidden)

# ------------------------------------------------------------------------------------------
def zip_plot(fname):
    Z = pd.read_csv(f"{TAB}/V4_S6_zip.csv"); Z = Z[Z.usable]
    Ms = ["gee_fg", "gee_md", "gee_kc", "gee", "mmm", "glmer", "shrink"]; Ns = [10, 45, 145]
    fig, axes = plt.subplots(3, 7, figsize=(7.2, 5.6), sharey=True, gridspec_kw=dict(hspace=0.55, wspace=0.12))
    fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.06)
    for i, n in enumerate(Ns):
        for j, m in enumerate(Ms):
            ax = axes[i, j]; d = Z[(Z.N == n) & (Z.method == m)].copy()
            z = (d.estimate - d.truth) / d.se; d["rk"] = 100 * (np.argsort(np.argsort(np.abs(z.values))) + 0.5) / len(d)
            cov = np.abs(z.values) <= 1.959964
            lo = d.estimate - 1.959964 * d.se; hi = d.estimate + 1.959964 * d.se
            ax.hlines(d.rk[cov], lo[cov], hi[cov], color="#9db4d6", lw=0.35)
            ax.hlines(d.rk[~cov], lo[~cov], hi[~cov], color="#c2531b", lw=0.35)
            t = d.truth.iloc[0]; ax.axvline(t, color="#222222", lw=0.6)
            p = cov.mean(); se = np.sqrt(p * (1 - p) / len(d))
            ax.axhspan(100 * (p - 1.96 * se), 100 * (p + 1.96 * se), color="#222222", alpha=0.25, lw=0)
            ax.axhline(95, color="#222222", lw=0.5, ls=":")
            w = np.quantile(np.abs(np.r_[lo - t, hi - t]), 0.9)
            ax.set_xlim(t - w, t + w); ax.set_xticks([round(t, 2)]); ax.tick_params(axis="x", labelsize=6)
            RECORD.append(("zip", f"N={n}", m, "coverage", 100 * p, 100 * (p - 1.96 * se), 100 * (p + 1.96 * se), f"V4_S6_zip:{n}|{m}"))
            if i == 0: ax.set_title(LAB[m], loc="center", fontsize=7)
            if j == 0: ax.set_ylabel(f"N = {n}\ncentile of |z|")
    fig.savefig(f"{OUT}/{fname}.png"); fig.savefig(f"{OUT}/{fname}.pdf")
    return fig

if __name__ == "__main__":
    rule = json.load(open(f"{OUT}/fig1_rule.json"))
    rep = {}
    f1 = figure1([tuple(x) for x in rule["lines"]], "S6_fig1"); rep["S6_fig1"] = bbox_check(f1)
    f2, clipped = nested_loop("S6_nestedloop"); rep["S6_nestedloop"] = bbox_check(f2)
    f3 = zip_plot("S6_zip"); rep["S6_zip"] = bbox_check(f3)
    pd.DataFrame(RECORD, columns=["figure", "panel", "series", "x", "est", "lo", "hi", "source"]).to_csv(f"{OUT}/S6_figures_data.csv", index=False)
    json.dump({k: {"overlaps": v[0], "outside": v[1]} for k, v in rep.items()} | {"nestedloop_clipped_above_30_and_hidden_below_min_obs": clipped}, open(f"{OUT}/S6_figures_bbox.json", "w"), indent=1)
    for k, v in rep.items(): print(k, "overlaps", len(v[0]), v[0][:6], "outside", v[1][:6])
    print("nested loop: values above 30% (clipped), points with n_obs < MIN_OBS (hidden):", clipped)
