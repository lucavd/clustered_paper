"""s6_figures.py -- S6: Figure 1 (v4), nested-loop plot and zip plot.

Every plotted value is read from tables/V4_S6_values.csv, V4_S6_nestedloop.csv or
V4_S6_zip.csv and written back, as plotted, to reports/morris/S6_figures/*_data.csv
(re-read by R/s6/s6_check_figures.R, gate G3). Style: tools/s6_figure.mplstyle
(rcParams of figure-style, sizes 8/7/6). Monte Carlo intervals: estimate +/- 2 MCSE
(block-bootstrap SE, B = 2000). Usage: python s6_figures.py <tables_dir> <out_dir> <style>
"""
import sys, json, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

TAB, OUT, STYLE = sys.argv[1], sys.argv[2], sys.argv[3]
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
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.62, left=0.14, right=0.89, top=0.95, bottom=0.11)
    # a: rule (text boxes)
    ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off(); letter(ax, "a")
    y = 0.95
    for kind, txt in rule_lines:
        fc = {"q": "#e8edf5", "b": "#f0f0f0", "lrt": "#efe9f5", "fg": "#fbeee6", "mmm": "#e6edf8", "end": "#f6f6f6"}[kind]
        ec = {"q": "#555555", "b": "#555555", "lrt": FOCAL["glmer_lrt"], "fg": FOCAL["gee_fg"], "mmm": FOCAL["mmm"], "end": "#555555"}[kind]
        ax.text(0.45, y, txt, transform=ax.transAxes, ha="center", va="top", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.35", fc=fc, ec=ec, lw=0.7), fontweight="bold" if kind == "q" else "normal")
        y -= 0.11 + 0.068 * txt.count("\n")   # D-35: four boxes
    # b: MAD by exact T, all eight, log scale
    ax = fig.add_subplot(gs[0, 1]); letter(ax, "b")
    Ts = [2, 3, 4, 6, 8, 10]; ends = []
    for m in ORD:
        e = np.array([val(f"madT1|{m}|{t}")[0] for t in Ts]); s = np.array([val(f"madT1|{m}|{t}")[1] for t in Ts])
        foc = m in FOCAL
        ax.plot(Ts, e, color=COL[m], lw=1.6 if foc else 0.8, marker="o", ms=3 if foc else 2, zorder=3 if foc else 2)
        ax.vlines(Ts, np.maximum(e - 2 * s, 0.3), e + 2 * s, color=COL[m], lw=0.8, zorder=2)
        for t, a, b in zip(Ts, e, s):
            RECORD.append(("fig1", "b", m, t, a, a - 2 * b, a + 2 * b, f"madT1|{m}|{t}"))
        ends.append(e[-1])
    ax.set_yscale("log"); ax.set_yticks([1, 2, 5, 10]); ax.set_yticklabels(["1", "2", "5", "10"])
    ax.set_xticks(Ts); ax.set_xlim(1.6, 10.6)
    ax.set_xlabel("Occasions per cluster"); ax.set_ylabel("Mean absolute departure\nfrom 5% (pp), lower = better")
    le = np.log10(ends); o = np.argsort(le); pos = le[o].copy()
    for i in range(1, len(pos)): pos[i] = max(pos[i], pos[i - 1] + 0.115)
    for k, i in enumerate(o):
        m = ORD[i]
        ax.plot([10.1, 10.9], [ends[i], 10 ** pos[k]], color=COL[m], lw=0.5, clip_on=False)
        ax.text(11.0, 10 ** pos[k], LAB[m], color=COL[m], fontsize=7, va="center", ha="left", clip_on=False)
    # c: share of simulated scenarios clearly above 7.5% false positives, by N (all series lengths)
    ax = fig.add_subplot(gs[1, 0]); letter(ax, "c")
    SF = pd.read_csv(f"{TAB}/V4_S6_rule_safety.csv"); SF = SF[(SF.band == "all") & SF.N.notna()]
    SC = [("mmm|t", "MMM, t", FOCAL["mmm"], 1.8, "-"), ("gee_fg|t", "FG, t", FOCAL["gee_fg"], 1.2, "-"), ("gee_fg|n", "FG, normal", FOCAL["gee_fg"], 0.9, "--"),
          ("gee_md|t", "MD, t", GREYS["gee_md"], 0.9, "-"), ("glmer_lrt|n", "LRT", FOCAL["glmer_lrt"], 1.2, "-")]
    ends = []
    for key, lab, col, lw, ls in SC:
        d = SF[SF.strategy == key].sort_values("N"); y = 100 * d.n_hi / d.n_dp
        ax.plot(d.N, y, color=col, lw=lw, ls=ls, marker="o", ms=3 if key == "mmm|t" else 2, zorder=4 if key == "mmm|t" else 2)
        for n_, yy, k_ in zip(d.N, y, d.n_hi):
            RECORD.append(("fig1", "c", key, int(n_), yy, np.nan, np.nan, f"V4_S6_rule_safety:{key}|{int(n_)}"))
        ends.append((y.iloc[-1], lab, col))
    ends.sort(); pos = [e[0] for e in ends]
    for i in range(1, len(pos)): pos[i] = max(pos[i], pos[i - 1] + 2.6)
    for (v_, lab, col), p_ in zip(ends, pos):
        ax.plot([146, 151], [v_, p_], color=col, lw=0.5, clip_on=False)
        ax.text(152, p_, lab, color=col, fontsize=7, va="center", ha="left", clip_on=False)
    ax.set_xlim(5, 150); ax.set_xticks([20, 40, 60, 80, 100, 120, 140]); ax.set_ylim(-1, 45)
    ax.set_xlabel("Clusters (e.g. patients)"); ax.set_ylabel("Scenarios clearly above 7.5%\nfalse positives (%), lower = better")
    # d: worst case vs mean departure on GenBinaryY (96 N x T combinations)
    ax = fig.add_subplot(gs[1, 1]); letter(ax, "d")
    for m in ORD:
        we, ws = val(f"wc96|GB|{m}"); me, ms_ = val(f"mad96|GB|{m}")
        foc = m in FOCAL
        ax.errorbar(me, we, xerr=2 * ms_, yerr=2 * ws, fmt="o", color=COL[m], ms=5 if foc else 3.5, lw=0.8, capsize=0, zorder=3)
        RECORD.append(("fig1", "d", m, "wc", we, we - 2 * ws, we + 2 * ws, f"wc96|GB|{m}"))
        RECORD.append(("fig1", "d", m, "mad", me, me - 2 * ms_, me + 2 * ms_, f"mad96|GB|{m}"))
        dx, dy, ha = 0.18, 0.9, "left"
        if m == "mmm": dx, dy = 0.2, -1.6
        if m == "gee_fg": dx, dy = -0.2, 1.2; ha = "right"
        if m in ("glmer", "gee"): dx, dy, ha = -0.2, 0.9, "right"
        ax.text(me + dx, we + dy, LAB[m], color=COL[m], fontsize=7, ha=ha, va="center")
    ax.set_xlabel("Mean departure from 5% (pp)"); ax.set_ylabel("Worst-case departure (pp)\nlower = better on both")
    ax.set_xlim(0.5, 8.6); ax.set_ylim(0, 28)
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
