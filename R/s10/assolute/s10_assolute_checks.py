#!/usr/bin/env python3
"""S10 track 'assolute': one computation per absolute claim of the rendered
manuscript and supplement, over all rows of the deciding table.
Run from the repository root: python3 R/s10/assolute/s10_assolute_checks.py <recs_dir>
Writes the verdicts to results_v4/S10/assolute/checks_raw.csv."""
import pandas as pd, numpy as np, sys, os, json, math
from scipy import stats
T = "tables/"
RECS = sys.argv[1] if len(sys.argv) > 1 else "results_v4/S10/assolute"
RES = {}
def rec(cid, verdict, table, computation, recomputed, counter="", note="", severity="info"):
    RES.setdefault(cid, []).append(dict(cand_id=cid, verdict=verdict, table=table, computation=computation,
                    recomputed=str(recomputed), counterexample=str(counter), note=note, severity=severity))
def cp_lower(k, n, a=0.05):
    return 0.0 if k == 0 else stats.beta.ppf(a/2, k, n-k+1)
def cp_upper(k, n, a=0.05):
    return 1.0 if k == n else stats.beta.ppf(1-a/2, k+1, n-k)
rd = pd.read_csv(T+"V5_S7f_rule_dp.csv"); t4 = pd.read_csv(T+"V5_S7f_table4.csv")
fd = pd.read_csv(T+"V5_S7g_filter_dp.csv")
KEY = ["dgp","scen","N","ntime"]

# ---------- Part 1: Table 5 and the decision rule ----------
f10 = fd[fd.L == 10]
# recompute flags from n, k with Clopper-Pearson, on every design point and strategy
f10 = f10.assign(lo95=[cp_lower(k, n) if n > 0 else np.nan for k, n in zip(f10.k, f10.n)])
f10 = f10.assign(hi_re=(f10.assess) & (f10.lo95 > 0.075))
mis = (f10.hi_re != f10.hi).sum()
hi_by = rd.groupby("strategy").hi.sum(); lo_by = rd.groupby("strategy").lo.sum()
intv = t4[t4.has_interval].strategy.tolist()
zero_hi = [s for s in intv if hi_by[s] == 0]
na_mmm = int((~f10[f10.strategy == "mmm_cc"].assess).sum())
txt = f"strategies with interval and 0 dps clearly above 7.5% over 2,376 dps: {zero_hi}; next fewest: " + \
      ", ".join(f"{s}={int(hi_by[s])}" for s in hi_by[intv].sort_values().index[1:4]) + \
      f"; MMM-cc dps not assessed (<20 usable): {na_mmm}; CP recomputation disagreements with stored flags: {mis}"
ok = zero_hi == ["mmm_cc"] and na_mmm == 0 and mis == 0
for cid in ["A0001", "A0072", "A0093", "A0122", "A0127"]:
    rec(cid, "TRUE" if ok else "FALSE", "tables/V5_S7f_rule_dp.csv; tables/V5_S7g_filter_dp.csv (L=10)",
        "per strategy, count dps with Clopper-Pearson lower limit > 7.5% (n>=20) over all 2,376 dps (main+M1c 1,512, M6-M7 384, Copula-AR1 480); check every strategy with an interval",
        txt, "" if ok else txt)
rec("A0054", "TRUE" if (t4.n_dp == 2376).all() and rd.groupby("strategy").size().eq(2376).all() else "FALSE",
    "tables/V5_S7f_rule_dp.csv", "design points per strategy in the rule table", f"rows per strategy: {sorted(set(rd.groupby('strategy').size()))}")
# GEE strategy flagged design points
fl = rd[(rd.strategy == "gee_ar1_fg") & rd.hi][KEY]
fl_s = "; ".join(f"{r.dgp} {r.scen} N={r.N} T={r.ntime}" for r in fl.itertuples())
ok = len(fl) == 2 and (fl.N == 10).all() and (fl.ntime == 2).all() and (fl.scen == "M2").all() and set(fl.dgp) == {"GLMM-AR1", "GLMM-native"}
for cid in ["A0073", "A0094", "A0269"]:
    rec(cid, "TRUE" if ok else "FALSE", "tables/V5_S7f_rule_dp.csv", "list dps with hi for gee_ar1_fg over all 2,376", fl_s)
rec("A0123", "TRUE" if ok else "FALSE", "tables/V5_S7f_rule_dp.csv", "first clause: dps with hi for gee_ar1_fg", fl_s,
    note="second clause checked with the head-to-head table (part 4)")
# exch / ind FG
out = []
okx = True
for s, nall, nmain in [("gee_exch_fg", 15, 3), ("gee_ind_fg", 5, 3)]:
    x = rd[(rd.strategy == s) & rd.hi]
    xm = x[x.grp == "main"]
    contains = all(((x.dgp == r.dgp) & (x.scen == r.scen) & (x.N == r.N) & (x.ntime == r.ntime)).any() for r in fl.itertuples())
    out.append(f"{s}: {len(x)} hi ({len(xm)} main+M1c, by grp {x.grp.value_counts().to_dict()}), N values {sorted(set(x.N))}, contains both GEE-strategy dps: {contains}")
    okx &= len(x) == nall and len(xm) == nmain and set(xm.N) == {10} and contains
    broad_counter = x[x.N != 10]
rec("A0075", "TRUE" if okx else "FALSE", "tables/V5_S7f_rule_dp.csv", "hi dps of exchangeable and independence FG (t): count, count in main grid+M1c, N values, inclusion of the two GEE-strategy dps",
    " | ".join(out), "broad reading ('all' = all 15 and 5) is false: e.g. exchangeable FG at Copula-AR1 M1 N=15 T=4, GenBinaryY M7 N=25 T=4; independence FG at Copula-AR1 M1 N=25 T=4",
    note="TRUE only if 'all with ten clusters' refers to the 3 and 3 in the main grid (all N=10 there); 15 and 5 overall include N=15-25 in M6-M7 and Copula-AR1. Suggest 'all three with ten clusters'. 'main grid' here = main grid + M1c (Table 5 column)", severity="minor")
# Copula-AR1
c = rd[(rd.grp == "copula") & rd.strategy.isin(["gee_ar1_fg", "mmm_cc"])]
rec("A0080", "TRUE" if c.hi.sum() == 0 and c[KEY].drop_duplicates().shape[0] == 480 else "FALSE", "tables/V5_S7f_rule_dp.csv",
    "hi over the 480 Copula-AR1 dps for gee_ar1_fg and mmm_cc", f"hi: {c.groupby('strategy').hi.sum().to_dict()}; dps {c[KEY].drop_duplicates().shape[0]}")
# M6-M7 MMM never above
m = rd[(rd.grp == "m67") & (rd.strategy == "mmm_cc")]
rec("A0066", "TRUE" if m.hi.sum() == 0 and len(m) == 384 else "FALSE", "tables/V5_S7f_rule_dp.csv", "hi of mmm_cc (t) over the M6-M7 dps",
    f"{int(m.hi.sum())} of {len(m)}")
rec("A0109", "TRUE" if m.hi.sum() == 0 else "FALSE", "tables/V5_S7f_rule_dp.csv", "hi and lo of mmm_cc (t) over the M6-M7 dps",
    f"hi {int(m.hi.sum())}, lo {int(m.lo.sum())} of {len(m)}",
    note="'remained calibrated' read as never clearly above 7.5%; clearly-below count reported for completeness")
# sign-flip
fl_ = rd[rd.strategy == "flip"]
rec("A0082", "TRUE" if fl_.hi.sum() == 0 and fl_.lo.sum() == 0 and len(fl_) == 2376 else "FALSE", "tables/V5_S7f_rule_dp.csv",
    "hi and lo of flip over all 2,376 dps", f"hi {int(fl_.hi.sum())}, lo {int(fl_.lo.sum())} of {len(fl_)}",
    note="the mean rate 5.4 is read from tables/V4_S7d2_flip_rules.csv (2026-10-25, before S9); its usable share (96.507%) equals the post-S9 Table 5 value, so the constant-outcome rule did not change the sign-flip results")
# small designs
sm = rd[rd.small]
g = sm.groupby("strategy").agg(hi=("hi", "sum"), lo=("lo", "sum"), pu=("pu", "mean"), n_ge50=("pu", lambda x: int((x >= 50).sum())), n=("pu", "size"))
s2 = g.loc[["mmm_cc", "gee_ar1_fg"]]
okh = (s2.pu < 50).all()
txt = f"mean usable share over the {int(s2.n.iloc[0])} small dps: MMM {s2.pu['mmm_cc']:.1f}%, GEE {s2.pu['gee_ar1_fg']:.1f}%; dps with usable >= 50%: MMM {s2.n_ge50['mmm_cc']}, GEE {s2.n_ge50['gee_ar1_fg']}"
for cid in ["A0076", "A0104", "A0125", "A0158"]:
    rec(cid, "TRUE" if okh else "FALSE", "tables/V5_S7f_rule_dp.csv (small = T<=3 and N<=25)", "mean pu over small dps for mmm_cc and gee_ar1_fg", txt,
        note="true for the mean over the 216 small dps; at a minority of the small dps the share is 50% or more", severity="info")
gl = g.loc["glmer"]
rec("A0003", "TRUE" if okh else "FALSE", "tables/V5_S7f_rule_dp.csv", "mean pu over small dps (MMM, GEE strategy; also glmer in case 'both' is read as the mixed model); flip hi/lo in small dps",
    txt + f"; glmer {gl.pu:.1f}%; flip in small dps: hi {int(g.loc['flip','hi'])}, lo {int(g.loc['flip','lo'])}, pu {g.loc['flip','pu']:.1f}%",
    note="'both' follows a sentence about the mixed-model Wald test, so its antecedent is ambiguous; the claim holds for either reading", severity="minor")
# no strategy both available and calibrated
gi = g.loc[intv]
avail = gi[gi.pu >= 50]
txt = "strategies with interval and mean usable >= 50% in small dps: " + ", ".join(f"{s} (pu {r.pu:.1f}%, hi {int(r.hi)} of {int(r.n)})" for s, r in avail.iterrows()) + \
      "; strategies with 0 hi in small dps: " + ", ".join(f"{s} (pu {r.pu:.1f}%, lo {int(r.lo)})" for s, r in gi[gi.hi == 0].iterrows())
okc = (avail.hi > 0).all()
rec("A0103", "TRUE" if okc else "FALSE", "tables/V5_S7f_rule_dp.csv", "for each strategy with interval in the 216 small dps: mean pu and hi count; 'available' read as mean usable >= 50%, 'calibrated' as no dp clearly above 7.5%",
    txt, note="'available' and 'calibrated' are not defined in the text; the verdict holds with usable >= 50% and 0 dps clearly above; the GEE strategy is 'not calibrated' only through 2 dps and it is 63 dps clearly below 2.5%", severity="minor")
rec("A0126", "TRUE" if okc else "FALSE", "tables/V5_S7f_rule_dp.csv", "as A0103 ('estimates the effect reliably' read as usable and calibrated)", txt,
    note="'reliably' is undefined; judged as A0103", severity="info")
# counts of the columns
n10 = rd[(rd.strategy == "mmm_cc") & (rd.N == 10)].shape[0]; nsm = rd[(rd.strategy == "mmm_cc") & rd.small].shape[0]
chk = rd[(rd.strategy == "mmm_cc") & rd.small]
rec("A0150", "TRUE" if n10 == 156 and nsm == 216 and (chk.ntime <= 3).all() and (chk.N <= 25).all() else "FALSE", "tables/V5_S7f_rule_dp.csv",
    "dps at N=10 and dps flagged small (check that small = T<=3 and N<=25)", f"N=10: {n10}; small: {nsm}; small all T<=3 & N<=25: {(chk.ntime<=3).all() and (chk.N<=25).all()}; T<=3&N<=25 count {rd[(rd.strategy=='mmm_cc')&(rd.ntime<=3)&(rd.N<=25)].shape[0]}")

# ---------- Part 2: Table 3, Table S3, Table S11, Table S12 and per-N rates ----------
tc = pd.read_csv(T+"variant_nc/V4_S6_table_cells.csv"); vals = pd.read_csv(T+"variant_nc/V4_S6_values.csv")
cl = pd.read_csv(T+"variant_nc/V4_S5_claims.csv"); pp = pd.read_csv(T+"variant_nc/V4_perf_primary.csv")
V = vals[vals.variant == "main"].set_index("name"); Vc = vals[vals.variant == "m1c"].set_index("name")
meths = ["gee_fg", "gee_md", "mmm", "glmer", "shrink", "gee", "glmer_lrt", "gee_kc"]
out, ok = [], True
for m3 in ["GB", "NAT", "AR1"]:
    w = {m: V.loc[f"wc|{m3}|{m}", "est"] for m in meths}
    best = min(w, key=w.get); second = sorted(w.values())[1]
    c01 = cl[cl.name == f"C01|{m3}"].iloc[0]
    out.append(f"{m3}: min worst-case {best} {w[best]:.2f}, next {second:.2f} (C01 lower limit {c01.lo:.2f})")
    ok &= best == "mmm" and c01.lo > 0
rec("A0060", "TRUE" if ok else "FALSE", "tables/variant_nc/V4_S6_values.csv (wc|mech|method, main); V4_S5_claims.csv C01",
    "argmin over the 8 methods of the worst-case departure, per mechanism; bootstrap lower limit of the difference to the runner-up > 0", "; ".join(out))
c03 = cl[cl.claim_id == "C03"]
ok = ((c03.lo > -c03.margin) & (c03.hi < c03.margin)).all() and len(c03) == 3
rec("A0061", "TRUE" if ok else "FALSE", "tables/variant_nc/V4_S5_claims.csv (C03)", "bootstrap interval of the matched-replicate difference MMM - FG within +/- margin on GB, NAT, AR1",
    "; ".join(f"{r.name}: [{r.lo:.2f}, {r.hi:.2f}] margin {r.margin}" for r in c03.itertuples()))
# t lowered every pooled rejection rate (Table S11 = tshift, main)
ts = V[V.index.str.startswith("tshift|")].est; tsc = Vc[Vc.index.str.startswith("tshift|")].est
ok = (ts < 0).all()
rec("A0071", "TRUE" if ok else "FALSE", "tables/variant_nc/V4_S6_values.csv (tshift|mech|method = t minus normal)", "t minus normal pooled rate < 0 for every one of the 7 Wald methods and 3 mechanisms (main grid; M1c reported)",
    f"main: {len(ts)} shifts, largest (closest to 0) {ts.max():.3f} ({ts.idxmax()}); M1c: largest {tsc.max():.3f} ({tsc.idxmax()}), all < 0: {(tsc < 0).all()}", "" if ok else ts[ts >= 0].to_dict())
# coverage on every mechanism (Table S12)
cv = {k: V.loc[k, "est"] for k in V.index if k.startswith("cov|")}
rows = []; ok = True
for m in ["gee_fg", "gee_md", "mmm"]:
    for m3 in ["GB", "NAT", "AR1"]:
        x = cv[f"cov|{m3}|{m}"]; rows.append(f"{m} {m3} {x:.2f}"); ok &= 93.5 <= x < 95.5
kc = [cv[f"cov|{m3}|gee_kc"] for m3 in ["GB", "NAT", "AR1"]]; ok &= all(92.5 <= x < 94.5 for x in kc)
cop = pd.read_csv(T+"V4_S7e_bias_cov.csv") if os.path.exists(T+"V4_S7e_bias_cov.csv") else None
rec("A0086", "TRUE" if ok else "FALSE", "tables/variant_nc/V4_S6_values.csv (cov|mech|method, main grid)",
    "coverage of FG, MD, MMM on GB, NAT, AR1 rounds to 94-95; KC rounds to 93-94", "; ".join(rows) + f"; KC {', '.join(f'{x:.2f}' for x in kc)}",
    note="'every mechanism' covers the three mechanisms of Table S12; on Copula-AR1 (Table S28, normal reference) FG AR(1) 94.7, MMM 93.8, MD exchangeable 93.9, KC exchangeable 93.0, so the statement also holds there")
# per-N type I rates from per-cell table, replicate weighted
h0 = pp[(pp.arm == "H0") & (pp.scen != "M1c")].copy()
h0["k"] = h0.rej * h0.n_obs
rN = h0.groupby(["dgp", "method", "N"]).apply(lambda g: 100 * g.k.sum() / g.n_obs.sum(), include_groups=False).rename("rate").reset_index()
chk = V.loc["rate|GB|gee", "est"]; rec_ = 100 * h0[(h0.dgp == "GenBinaryY") & (h0.method == "gee")].k.sum() / h0[(h0.dgp == "GenBinaryY") & (h0.method == "gee")].n_obs.sum()
g = rN[rN.method == "gee"].groupby("dgp").apply(lambda d: d.loc[(d.rate - 5).idxmin()], include_groups=False)
ok = ((g.rate - 5) < 0.5).all() and ((g.rate - 5) > 0).all()
rec("A0176", "TRUE" if ok else "FALSE", "tables/variant_nc/V4_perf_primary.csv (H0, M1-M5, replicate-weighted per N)",
    "uncorrected GEE: min over the 16 N of (rate_N - 5) per mechanism < 0.5 (and > 0)", "; ".join(f"{d}: N={int(r.N)} excess {r.rate-5:.2f}" for d, r in g.iterrows()) + f" [check: pooled GB rate recomputed {rec_:.4f} vs table {chk:.4f}]")
x = rN[rN.method.isin(["gee_fg", "mmm"]) & (rN.N > 45) & rN.dgp.isin(["GLMM-native", "GLMM-AR1"])].copy(); x["dev"] = (x.rate - 5).abs()
mx = x.loc[x.dev.idxmax()]
rec("A0177", "TRUE" if x.dev.max() <= 0.65 else "FALSE", "tables/variant_nc/V4_perf_primary.csv (H0, M1-M5, per N)",
    "max |rate_N - 5| for FG and MMM at N 55-145 on GLMM-native and GLMM-AR1 ('about half a point' read as <= 0.65)",
    f"max {mx.dev:.2f} ({mx.method}, {mx.dgp}, N={int(mx.N)}); per method max: " + ", ".join(f"{k} {v:.2f}" for k, v in x.groupby('method').dev.max().items()))

# ---------- Part 3: usable shares, switch, adaptive quadrature, working correlation, head-to-head, usability limits ----------
t7 = tc[tc.table_id == "T7"]; t7c = tc[tc.table_id == "T7c"]
ge4 = t7[t7.row.str.contains(r"T = (4|6|8|10)$")]; ge4c = t7c[t7c.row.str.contains(r"T = (4|6)$")]
mn = ge4.loc[ge4.est.idxmin()]
pT = rd.groupby(["strategy"]).apply(lambda g: g[g.ntime >= 4].pu.mean(), include_groups=False)
rec("A0069", "TRUE" if ge4.est.min() > 90 else "FALSE", "tables/variant_nc/V4_S6_table_cells.csv (T7 = Table 4/S7, main grid)",
    "min over the 8 methods x RI/RS x T in {4,6,8,10} of the pooled usable share > 90",
    f"min {mn.est:.2f} ({mn['col']}, {mn['row']}); M1c (T7c) min {ge4c.est.min():.2f}; scope extension, Table 5 strategies over all 2,376 dps with T>=4 (dp mean): min {pT.min():.1f} ({pT.idxmin()})",
    note="pooled over cells; per design point the share is lower at N=10 (not claimed)")
rs2 = t7[(t7.row == "random slopes, T = 2") & (t7["col"] == "GLMM, Wald")].est.iloc[0]
rec("A0106", "TRUE" if rs2 < 0.5 else "FALSE", "tables/variant_nc/V4_S6_table_cells.csv (T7)", "usable share of the GLMM Wald test with random slopes at T=2", f"{rs2:.3f}%")
rec("A0148", "TRUE" if t7.row.nunique() == 12 and t7["col"].nunique() == 8 else "FALSE", "tables/variant_nc/V4_S6_table_cells.csv (T7)",
    "Table S7 has 8 methods x (RI, RS) x T in {2,3,4,6,8,10}; Table 4 methods are a subset", f"rows {t7.row.nunique()}, methods {t7['col'].nunique()}")
va = V[V.index.str.startswith("vaN|")].est.reset_index(); va[["p", "mech", "method", "N"]] = va.name.str.split("|", expand=True); va.N = va.N.astype(int)
arg = va.loc[va.groupby(["mech", "N"]).est.idxmax()]
lrtN = arg[arg.method == "glmer_lrt"].groupby("mech").N.agg(["min", "max", "count"])
ok = all(lrtN.loc[m, "max"] == 25 and lrtN.loc[m, "count"] == 4 for m in ["GB", "NAT", "AR1"])
rec("A0070", "TRUE" if ok else "FALSE", "tables/variant_nc/V4_S6_values.csv (vaN|mech|method|N, post-S9)", "argmax over 8 methods of unconditional validity at each N; LRT selected exactly at N in {10,15,20,25} on each mechanism",
    "; ".join(f"{m}: LRT at N {r['min']}-{r['max']} ({r['count']} N)" for m, r in lrtN.iterrows()),
    note="the registry number x:switch_GB is read from tables/V4_S5_switch_points.csv dated 2026-10-18, before the constant-outcome rule; recomputation on the post-S9 table gives the same switch (25)", severity="minor")
# adaptive quadrature (nAGQ=10) vs Laplace (nAGQ=1)
ag = pd.read_csv(T+"V5_S7g_agq_dp.csv")
L1, Q = ag[ag.nAGQ == 1], ag[ag.nAGQ == 10]
qh = Q[Q.hi_t]
ok2 = set(qh.ntime) == {10}
rec("A0002", "TRUE" if ok2 else "FALSE", "tables/V5_S7g_agq_dp.csv (nAGQ=10, t reference)", "occasions of the dps clearly above 7.5% with quadrature",
    f"{len(qh)} dps, T values {sorted(set(qh.ntime))}; Laplace {int(L1.hi_t.sum())} dps, T values {sorted(set(L1[L1.hi_t].ntime))}",
    note="the subset of 123 dps of Section 3.10, not the whole grid")
cop = qh[qh.dgp == "Copula-AR1"]; gb = qh[qh.dgp == "GenBinaryY"]
ok = len(L1[L1.hi_t]) == 33 and len(qh) == 8 and len(cop) == 6 and len(gb) == 2 and set(qh.ntime) == {10}
rec("A0087", "TRUE" if ok else "FALSE", "tables/V5_S7g_agq_dp.csv", "hi counts Laplace/quadrature; T, mechanism and rate of the quadrature flags",
    f"Laplace {int(L1.hi_t.sum())}, quadrature {len(qh)} (Copula {len(cop)}: {cop.rate_t.min():.1f}-{cop.rate_t.max():.1f}%; GB {len(gb)}: {', '.join(f'{x:.1f}' for x in sorted(gb.rate_t))}%; other {len(qh)-len(cop)-len(gb)}); T {sorted(set(qh.ntime))}")
nat = Q[Q.dgp == "GLMM-native"]
rec("A0088", "TRUE" if nat.hi_t.sum() == 0 else "FALSE", "tables/V5_S7g_agq_dp.csv", "quadrature hi on GLMM-native dps of the subset", f"{int(nat.hi_t.sum())} of {len(nat)} (Laplace {int(L1[L1.dgp=='GLMM-native'].hi_t.sum())})")
lo = Q[Q.lo_t]
rec("A0124", "TRUE" if ok2 and len(cop) > len(qh) / 2 else "FALSE", "tables/V5_S7g_agq_dp.csv", "quadrature: T of hi dps; share on Copula-AR1; N of lo dps",
    f"hi only at T=10 ({len(qh)}), {len(cop)} of them Copula-AR1; lo dps by N: {lo.groupby('N').size().to_dict()}",
    note="'with few clusters the test became conservative': 31 of 37 clearly-below dps have N<=25, but 5 have N=55 and 1 N=45", severity="info")
# FG across working correlations (Section 3.5 / S3)
w = f10[f10.strategy.isin(["gee_exch_fg", "gee_ind_fg", "gee_ar1_fg"]) & (f10.dgp != "Copula-AR1")].pivot_table(index=KEY, columns="strategy", values="rate")
dd = {s: w[s] - w["gee_exch_fg"] for s in ["gee_ind_fg", "gee_ar1_fg"]}
cnt = {s: int((d.abs() > 2.5).sum()) for s, d in dd.items()}; mx = {s: d.abs().max() for s, d in dd.items()}
rec("A0178", "TRUE" if len(set(cnt.values())) == 1 else "FALSE", "tables/V5_S7g_filter_dp.csv (L=10, t reference, 1,896 dps of the three study mechanisms); tables/V5_S7g_fg_wc.csv",
    "per dp, FG rate with independence or AR(1) minus exchangeable; count |d|>2.5 and max |d| for each correlation",
    f"independence: {cnt['gee_ind_fg']} > 2.5 (max {mx['gee_ind_fg']:.2f}); AR(1): {cnt['gee_ar1_fg']} > 2.5 (max {mx['gee_ar1_fg']:.2f}); 5%/95% quantiles ind {dd['gee_ind_fg'].quantile(.05):.2f}/{dd['gee_ind_fg'].quantile(.95):.2f}, ar1 {dd['gee_ar1_fg'].quantile(.05):.2f}/{dd['gee_ar1_fg'].quantile(.95):.2f}",
    f"AR(1): {cnt['gee_ar1_fg']} differences beyond 2.5 (max {mx['gee_ar1_fg']:.1f}), not 18; the registry number f:g7|fgwc_gt25_d_ind is the independence count only",
    note="'for each correlation 18 ... (at most 5.1)' holds only for independence; suggested: '18 with independence and 11 with AR(1) exceeded 2.5 points (at most 5.1 and 3.7)'", severity="minor")
c5 = pd.read_csv(T+"V5_S7f_claims.csv").set_index("name")
ok = all(abs(c5.loc[n, "lo"]) < 0.5 and abs(c5.loc[n, "hi"]) < 0.5 for n in ["wc_ind", "wc_ar1"])
t2r = pd.read_csv(T+"V5_S7f_rule_t2.csv")
ch = t2r[t2r.chosen].groupby("df").strategy.apply(list).to_dict()
rec("A0098", "TRUE" if ok and all(v == ["mmm_cc"] for v in ch.values()) else "FALSE", "tables/V5_S7f_claims.csv (wc_ind, wc_ar1); tables/V5_S7f_rule_t2.csv",
    "mean FG difference between working correlations, bootstrap interval inside +/-0.5; strategy chosen by the rule on 1,896 dps",
    f"ind-exch {c5.loc['wc_ind','est']:.2f} [{c5.loc['wc_ind','lo']:.2f}, {c5.loc['wc_ind','hi']:.2f}], ar1-exch {c5.loc['wc_ar1','est']:.2f} [{c5.loc['wc_ar1','lo']:.2f}, {c5.loc['wc_ar1','hi']:.2f}]; chosen {ch}")
ok = all(abs(c5.loc[n, "lo"]) < 1 and abs(c5.loc[n, "hi"]) < 1 for n in ["cop_gee", "cop_mmm"])
alt = {s: f10[(f10.strategy == s) & (f10.dgp == "Copula-AR1")].rate.mean() - f10[(f10.strategy == s) & (f10.dgp != "Copula-AR1")].rate.mean() for s in ["gee_ar1_fg", "mmm_cc"]}
rec("A0081", "TRUE" if ok else "FALSE", "tables/V5_S7f_claims.csv (cop_gee, cop_mmm); cross-check tables/V5_S7g_filter_dp.csv",
    "Copula-AR1 minus study-mechanism rate, interval within +/-1", f"GEE {c5.loc['cop_gee','est']:.2f} [{c5.loc['cop_gee','lo']:.2f}, {c5.loc['cop_gee','hi']:.2f}], MMM {c5.loc['cop_mmm','est']:.2f} [{c5.loc['cop_mmm','lo']:.2f}, {c5.loc['cop_mmm','hi']:.2f}]; dp-mean recomputed GEE {alt['gee_ar1_fg']:.2f}, MMM {alt['mmm_cc']:.2f}")
# head-to-head
hb = pd.read_csv(T+"V5_S7f_h2h_boot.csv")
def H(region, q): r = hb[(hb.region == region) & (hb.q == q)].iloc[0]; return r
a = [H("all", q) for q in ["mad_gee_minus_mmm", "pu_gee_minus_mmm", "upow_mmm_minus_gee"]]
ok = all(r["class"] == "similar" for r in a)
txt = "; ".join("%s %.2f +/- %.2f (%s)" % (r.q, r.est, 2*r.mcse, r["class"]) for r in a) + f"; MAD GEE {H('all','mad_gee').est:.2f}, MMM {H('all','mad_mmm').est:.2f}"
for cid in ["A0078", "A0096"]:
    rec(cid, "TRUE" if ok else "FALSE", "tables/V5_S7f_h2h_boot.csv (region all)", "class of the three differences over the 2,376 dps", txt)
r = {q: H("map|10|2", q).est for q in ["pu_gee", "pu_mmm", "mad_gee", "mad_mmm"]}
rec("A0079", "TRUE" if max(r["pu_gee"], r["pu_mmm"]) < 15 and min(r["mad_gee"], r["mad_mmm"]) > 2.5 else "FALSE", "tables/V5_S7f_h2h_boot.csv (map|10|2)",
    "usable share and MAD of both strategies at N=10, T=2", f"usable GEE {r['pu_gee']:.1f}%, MMM {r['pu_mmm']:.1f}%; MAD GEE {r['mad_gee']:.1f}, MMM {r['mad_mmm']:.1f}")
mp = hb[hb.region.str.startswith("map|") & (hb.q == "mad_gee_minus_mmm")]
beyond = mp[mp["class"].isin(["GEE", "MMM"])]; about = mp[mp["class"] == "about the margin"]
ok = set(beyond.region.str.split("|").str[1]) == {"10"} and (beyond["class"] == "GEE").all()
rec("A0123", "TRUE" if ok else "FALSE", "tables/V5_S7f_h2h_boot.csv (map regions, mad_gee_minus_mmm)",
    "second clause: regions where the MAD difference lies beyond the 1-point margin, and which strategy is closer",
    "beyond: " + ", ".join("%s (%s, %.2f)" % (x.region, x["class"], x.est) for _, x in beyond.iterrows()) + "; about the margin: " + ", ".join("%s (%+.2f)" % (x.region, x.est) for _, x in about.iterrows()),
    note="fragile: at T=2 with 15-25 and >25 clusters (and at N=10, T=2) the interval straddles the margin with MMM closer by about one point; Section 2.5 says such intervals are reported as 'differences of about the size of the margin', which the Conclusions sentence does not mention", severity="minor")
r10 = {x.region: x["class"] for _, x in mp.iterrows() if x.region.startswith("map|10|")}
rec("A0128", "TRUE" if r10.get("map|10|3") == "GEE" and r10.get("map|10|4+") == "GEE" and r10.get("map|10|2") != "GEE" else "FALSE",
    "tables/V5_S7f_h2h_boot.csv", "class of the MAD difference at N=10 by T", str(r10))
# usability limits (Section S11)
hiL = fd.groupby(["L", "strategy"]).hi.sum().unstack(0)
puL = fd.groupby(["L", "strategy"]).pu.mean().unstack(0)
ok = hiL.loc["mmm_cc", 5] == 0 and hiL.loc["mmm_cc", 20] == 0
rec("A0290", "TRUE" if ok else "FALSE", "tables/V5_S7g_filter_dp.csv", "hi of mmm_cc over 2,376 dps with limits 5 and 20", f"L=5: {int(hiL.loc['mmm_cc',5])}, L=20: {int(hiL.loc['mmm_cc',20])}")
elig5 = [s for s in intv if hiL.loc[s, 5] == 0]; best5 = max(elig5, key=lambda s: puL.loc[s, 5])
ok = hiL.loc["gee_ar1_fg", 5] == 0 and hiL.loc["gee_ar1_fg", 20] == 6 and best5 == "gee_ar1_fg"
for cid in ["A0291", "A0074"]:
    rec(cid, "TRUE" if ok else "FALSE", "tables/V5_S7g_filter_dp.csv", "hi of gee_ar1_fg with L=5 and L=20; rule applied at L=5 (eligible = 0 hi, max mean pu)",
        f"GEE hi L=5 {int(hiL.loc['gee_ar1_fg',5])}, L=20 {int(hiL.loc['gee_ar1_fg',20])}; eligible at L=5: {elig5}; selected {best5} (pu {puL.loc[best5,5]:.1f} vs MMM {puL.loc['mmm_cc',5]:.1f})")
m10 = pd.DataFrame({"hi": hiL[10], "pu": puL[10]}).join(t4.set_index("strategy")[["hi_all", "pu"]], rsuffix="_t5")
dif = m10[(m10.hi != m10.hi_all) | ((m10.pu - m10.pu_t5).abs() > 1e-9)]
rec("A0289", "TRUE" if len(dif) == 0 else "FALSE", "tables/V5_S7g_filter_dp.csv (L=10) vs tables/V5_S7f_table4.csv",
    "hi count and mean usable share per strategy, L=10 recomputation vs Table 5", f"{len(m10)} strategies compared, {len(dif)} differ" + (f": {dif.to_dict()}" if len(dif) else ""))
okb = hiL.loc["gee_ar1_fg", 5] == 0
fgk = pd.read_csv(T+"V6_S9_keys.csv").set_index("key").value
rec("A0095", "TRUE" if okb and fgk["fgpt|1|n"] < 100 and fgk["fgpt|2|n"] < 100 else "FALSE", "tables/V6_S9_keys.csv; tables/V5_S7g_filter_dp.csv",
    "usable replicates of GEE and MMM at the two dps; GEE hi with the stricter limit 5",
    f"GEE n {int(fgk['fgpt|1|n'])} and {int(fgk['fgpt|2|n'])}, MMM n {int(fgk['fgpt|1|mmm_n'])} and {int(fgk['fgpt|2|mmm_n'])} (of about 975); GEE hi at L=5: {int(hiL.loc['gee_ar1_fg',5])}")
# N - 2 degrees of freedom
t2p = t2r.set_index(["df", "strategy"])
ok = t2p.loc[("N - 2", "mmm_cc"), "n_hi"] == 0 and t2p.loc[("N - 2", "mmm_cc"), "n_lo"] == 8 and t2p.loc[("N - p", "mmm_cc"), "n_lo"] == 23 and t2p.loc[("N - 2", "gee_ar1_fg"), "n_hi"] == 2 and all(v == ["mmm_cc"] for v in ch.values())
for cid in ["A0077", "A0234"]:
    rec(cid, "TRUE" if ok else "FALSE", "tables/V5_S7f_rule_t2.csv", "chosen strategy and hi/lo counts with N-2 and N-p over 1,896 dps",
        f"chosen {ch}; MMM hi N-2 {t2p.loc[('N - 2','mmm_cc'),'n_hi']}, lo N-2 {t2p.loc[('N - 2','mmm_cc'),'n_lo']} vs N-p {t2p.loc[('N - p','mmm_cc'),'n_lo']}; GEE hi N-2 {t2p.loc[('N - 2','gee_ar1_fg'),'n_hi']}, N-p {t2p.loc[('N - p','gee_ar1_fg'),'n_hi']}",
        note="table-level check; the N-2 per-dp flags were not recomputed")
mk = pd.read_csv(T+"V5_S8b_mdkc_avg.csv").set_index("strategy")
av = mk[mk.index.str.contains("avg")]
rec("A0232", "TRUE" if (av.hi_all > 0).all() else "FALSE", "tables/V5_S8b_mdkc_avg.csv", "hi_all of the averaged MD/KC strategies", av.hi_all.to_dict())

# ---------- Part 4: "largest / at most" Monte Carlo standard errors in captions ----------
def le(val, claimed):
    """printed bound 'claimed' (string) covers val: val <= claimed + half a unit of its last digit"""
    d = len(claimed.split(".")[1]) if "." in claimed else 0
    return val <= float(claimed) + 0.5 * 10 ** (-d) + 1e-12
def capres(cid, table, comp, pairs, note=""):
    ok = all(le(v, c) for v, c in pairs)
    rec(cid, "TRUE" if ok else "FALSE", table, comp, "; ".join(f"recomputed max {v:.4f} vs printed {c}" for v, c in pairs),
        "" if ok else "; ".join(f"{v:.4f} > {c}" for v, c in pairs if not le(v, c)), note=note, severity="info" if ok else "minor")
def tcmax(tid, colrx=None, rowrx=None, src=tc):
    x = src[src.table_id == tid]
    if colrx: x = x[x["col"].str.contains(colrx, regex=True)]
    if rowrx: x = x[x.row.str.contains(rowrx, regex=True)]
    return x.se.max()
TC = "tables/variant_nc/V4_S6_table_cells.csv (column se = block-bootstrap MCSE)"
capres("A0145", TC + " T3", "max se over Table 3 cells by column group", [(max(tcmax("T3", "^rate"), tcmax("T3", "^MAD")), "0.047"), (tcmax("T3", "^worst"), "0.32")])
t3top = pd.read_csv(T+"V4_S6_table_cells.csv"); s3 = pd.concat([tc[(tc.table_id == "T3") & (tc.i != 3)], t3top[(t3top.table_id == "T3") & (t3top.i == 3)]])
capres("A0181", TC + " T3 (rows 1-2,4-8) + tables/V4_S6_table_cells.csv T3 row 3 (MMM without check)", "max se by column group of Table S3",
       [(tcmax("T3", "^rate", src=s3), "0.047"), (tcmax("T3", "^worst", src=s3), "0.32"), (tcmax("T3", "^MAD", src=s3), "0.047")])
capres("A0182", TC + " T4", "max se of Table S4", [(tcmax("T4"), "0.50")])
capres("A0184", TC + " T5", "max se of Table S5", [(tcmax("T5"), "1.1")])
capres("A0186", TC + " T6", "max se of Table S6", [(tcmax("T6", "^(own|common)"), "0.047"), (tcmax("T6", "^shift"), "0.020")])
capres("A0188", TC + " T7", "max se of Table S7", [(tcmax("T7"), "0.13")])
t4cols = "GLMM, Wald|GLMM, likelihood ratio|Marginalized model|GEE, Fay-Graubard|GEE, uncorrected"
capres("A0147", TC + " T7 (Table 4 subset)", "max se over the Table 4 cells (5 methods, T in 2,3,4,6,10)", [(tcmax("T7", t4cols, r"T = (2|3|4|6|10)$"), "0.13")])
capres("A0189", TC + " T8", "max se of Table S8", [(tcmax("T8", "^usable$"), "0.096"), (tcmax("T8", "^unconditional"), "0.10")])
capres("A0190", TC + " T9", "max se of Table S9", [(tcmax("T9"), "1.4")])
capres("A0193", TC + " T10", "max se of Table S10", [(tcmax("T10"), "0.26")])
capres("A0194", TC + " T11", "max se of Table S11", [(tcmax("T11", "^(normal|t) "), "0.047"), (tcmax("T11", "^shift"), "0.035")])
capres("A0196", TC + " T12", "max se of Table S12", [(tcmax("T12"), "0.049")])
capres("A0218", TC + " T3c", "max se of Table S16", [(tcmax("T3c", "^rate"), "0.22"), (tcmax("T3c", "^worst"), "0.65"), (tcmax("T3c", "^MAD"), "0.22")])
capres("A0219", TC + " T4c", "max se of Table S17", [(tcmax("T4c"), "0.28")])
capres("A0221", TC + " T5c", "max se of Table S18", [(tcmax("T5c"), "5.0")])
capres("A0223", TC + " T6c", "max se of Table S19", [(tcmax("T6c", "^(own|common)"), "0.22"), (tcmax("T6c", "^shift"), "0.13")])
capres("A0225", TC + " T7c", "max se of Table S20", [(tcmax("T7c"), "0.35")])
capres("A0226", TC + " T8c", "max se of Table S21", [(tcmax("T8c", "^usable$"), "0.26"), (tcmax("T8c", "^unconditional"), "0.28")])
capres("A0227", TC + " T9c", "max se of Table S22", [(tcmax("T9c"), "0.91")])
capres("A0228", TC + " T10c", "max se of Table S23", [(tcmax("T10c", "^raw"), "0.24"), (tcmax("T10c", "^size"), "0.29")])
capres("A0229", TC + " T11c", "max se of Table S24", [(tcmax("T11c", "^(normal|t) "), "0.22"), (tcmax("T11c", "^shift"), "0.10")])
capres("A0231", TC + " T12c", "max se of Table S25", [(tcmax("T12c"), "0.20")])
# binomial bounds
bse = lambda p, n: 100 * np.sqrt(np.clip(p, 0, 1) * (1 - np.clip(p, 0, 1)) / n)
f10 = fd[fd.L == 10]
pu_b = f10.assign(b=bse(f10.pu / 100, f10.nd)).groupby("strategy").b.mean()
pus_b = f10[f10.ntime.le(3) & f10.N.le(25)].assign(b=lambda d: bse(d.pu / 100, d.nd)).groupby("strategy").b.mean()
cpd = pd.concat([pd.read_csv(T+"V5_S7f_cp_S7d_dp_rates.csv"), pd.read_csv(T+"V5_S7f_cp_S7e_dp_rates.csv")])
cpt = cpd[cpd.ref.eq("t") | cpd.strategy.isin(["glmer_lrt", "flip"])]
pw_b = cpt.assign(b=bse(cpt.upow / 100, cpt.nd)).groupby("strategy").b.mean()
capres("A0152", "tables/V5_S7g_filter_dp.csv (L=10); tables/V5_S7f_cp_S7d_dp_rates.csv + cp_S7e (t reference)", "max over strategies of the mean over dps of 100 sqrt(p(1-p)/nd)",
       [(pu_b.max(), "0.40"), (pw_b[pw_b.index.isin(t4.strategy)].max(), "1.0"), (pus_b.max(), "1.4")],
       note="power bound uses the H0 nd as denominator proxy (the H1 nd differs by at most a few datasets)")
u1 = pd.read_csv(T+"V5_S7f_u1.csv")
capres("A0197", "tables/V5_S7f_u1.csv", "max mcse_U1 (binomial)", [(u1.mcse_U1.max() * (100 if u1.mcse_U1.max() < 0.01 else 1), "0.044")])
ns = pd.read_csv(T+"V4_S6b_ns_rule_safety_dp.csv")
ns = ns[ns.strategy.isin(["mmm|t", "mmm|n", "mmm_nocc|t", "gee_fg|t", "gee_md|t", "glmer_lrt|n"])]
g = ns.groupby(["dgp", "scen", "strategy"])
mr = ns[ns.assess].groupby(["dgp", "scen", "strategy"]).apply(lambda d: bse(d.loc[d.rate.idxmax(), "rate"] / 100, d.loc[d.rate.idxmax(), "n"]), include_groups=False)
pb = g.apply(lambda d: bse(d.pu / 100, d.nd).mean(), include_groups=False)
capres("A0199", "tables/V4_S6b_ns_rule_safety_dp.csv", "per scenario x strategy row of Table S14: binomial MCSE at the dp with the largest rate; mean binomial SE of the usable share",
       [(mr.max(), "1.3"), (pb.max(), "0.55")], note=f"the {len(mr)} assessed scenario x strategy rows of Table S14 (6 strategies x 4 scenarios); largest rates reproduce the table column")
lt = pd.read_csv(T+"V5_S7f_lus_tab316.csv")
pc = pd.concat([lt[c] for c in ["pct_not_returned", "pct_no_se_500", "pct_sep_500", "pct_usable_500"]])
capres("A0208", "tables/V5_S7f_lus_tab316.csv", "max over cells of 100 sqrt(p(1-p)/500)", [(bse(pc / 100, 500).max(), "2.2")],
       note="upper bound 100*sqrt(0.25/500) = 2.236; printed 2.2 is the rounded maximum")
rec("A0207", "TRUE" if (lt.somma.round(6) == 100).all() else "FALSE", "tables/V5_S7f_lus_tab316.csv", "sum of the four % columns per method and N",
    f"rows {len(lt)}, sums {sorted(set(lt.somma.round(6)))}")
b7 = pd.read_csv(T+"V5_S7f_cp_S7b_rule_dp.csv")
s26 = pd.concat([f10[f10.dgp != "Copula-AR1"][["strategy", "pu", "nd"]], b7[["strategy", "pu", "nd"]]])
capres("A0236", "tables/V5_S7g_filter_dp.csv (L=10, three study mechanisms); tables/V5_S7f_cp_S7b_rule_dp.csv", "max over strategies of the mean over dps of 100 sqrt(p(1-p)/nd)",
       [(s26.assign(b=bse(s26.pu / 100, s26.nd)).groupby("strategy").b.mean().max(), "0.40")])
ef = pd.read_csv(T+"V4_S7d2_explain_finite.csv")
capres("A0251", "tables/V4_S7d2_explain_finite.csv", "max over the 24 rows of 100 sqrt(p(1-p)/datasets)",
       [(bse(ef.finite / 100, ef.datasets).max(), "0.32"), (bse(ef.gee_ar1_fg_usable / 100, ef.datasets).max(), "0.30")])
ce = pd.read_csv(T+"V5_S7f_cp_S7e_dp_rates.csv"); bc = pd.read_csv(T+"V4_S7e_bias_cov.csv")
ce_t = ce[ce.ref.eq("t") | ce.strategy.isin(["glmer_lrt", "flip"])]
bc = bc[bc.strategy.isin(["mmm_cc", "mmm_study", "gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_exch_md", "gee_exch_kc", "gee_exch_unc"])]
bc = bc.assign(rbm=bc.mcse_bias / bc.empse, cvm=bse(bc.cov_t / 100, bc.n))
s28s = ["gee_ar1_fg", "mmm_cc", "mmm_study", "gee_ind_fg", "gee_exch_fg", "gee_exch_md", "gee_exch_kc", "gee_exch_unc", "glmer", "shrink", "glmer_lrt"]
ce_t = ce_t[ce_t.strategy.isin(s28s + ["mmm_nocc", "mmm"])]
capres("A0258", "tables/V5_S7f_cp_S7e_dp_rates.csv; tables/V4_S7e_bias_cov.csv (Copula-AR1)", "max over strategies of the mean over the 480 dps of the binomial MCSE (rate, usable, power), of mcse_bias/empse and of the coverage MCSE",
       [(ce_t.groupby("strategy").mc.mean().max(), "1.1"), (ce_t.assign(b=bse(ce_t.pu / 100, ce_t.nd)).groupby("strategy").b.mean().max(), "0.49"),
        (ce_t.assign(b=bse(ce_t.upow / 100, ce_t.nd)).groupby("strategy").b.mean().max(), "1.1"), (bc.groupby("strategy").rbm.mean().max(), "0.036"),
        (max(bc.groupby("strategy").cvm.mean().max(), bc.groupby("strategy").mcse_cov_n.mean().max()), "0.91")],
       note="strategies of Table S28; coverage MCSE: stored mcse_cov_n (normal) and binomial on n (t)")
hb = pd.read_csv(T+"V5_S7f_h2h_boot.csv"); s31 = hb[hb.region.eq("all") | hb.region.str.startswith("map|")]
capres("A0273", "tables/V5_S7f_h2h_boot.csv (regions all + map|N|T)", "max mcse by quantity",
       [(s31[s31.q.str.startswith("rate")].mcse.max(), "0.47"), (s31[s31.q.isin(["mad_gee", "mad_mmm"])].mcse.max(), "0.42"),
        (s31[s31.q.isin(["pu_gee", "pu_mmm"])].mcse.max(), "0.30"), (s31[s31.q.isin(["upow_gee", "upow_mmm"])].mcse.max(), "0.069")])
calib = pd.read_csv(T+"V5_S7g_calib.csv")
mxc = [bse(r.max_rate / 100, r.max_n) for r in calib.itertuples()]
capres("A0283", "tables/V5_S7g_calib.csv", "max over the 19 strategies of the binomial MCSE at the dp of the highest rate", [(max(mxc), "6.0")])
s33 = fd.assign(b=bse(fd.pu / 100, fd.nd)).groupby(["strategy", "L"]).b.mean()
capres("A0292", "tables/V5_S7g_filter_dp.csv", "max over strategy x limit of the mean over dps of 100 sqrt(p(1-p)/nd)", [(s33.max(), "0.40")])
ag = pd.read_csv(T+"V5_S7g_agq_dp.csv"); ag2 = ag[ag.nAGQ.isin([1, 10])]
capres("A0297", "tables/V5_S7g_agq_dp.csv (nAGQ 1 and 10)", "max over mechanism x T x method of the mean over dps of mc_t and of the usable-share binomial SE",
       [(ag2.groupby(["dgp", "ntime", "nAGQ"]).mc_t.mean().max(), "3.1"), (ag2.assign(b=bse(ag2.pu / 100, ag2.nd)).groupby(["dgp", "ntime", "nAGQ"]).b.mean().max(), "1.0")])
c7 = pd.read_csv(T+"V5_S7f_cp_S7d_dp_rates.csv"); fs1 = c7[c7.scen.isin(["M6", "M7"]) & c7.strategy.isin(["gee_exch_md", "gee_exch_fg", "gee_ar1_fg", "mmm_cc", "glmer_lrt"]) & (c7.n >= 20)]
fs1 = fs1[fs1.ref.eq("t") | fs1.strategy.eq("glmer_lrt")]
capres("A0303", "tables/V5_S7f_cp_S7d_dp_rates.csv (M6-M7, strategies of Figure S1, n>=20)", "max 100 sqrt(p(1-p)/n) of a drawn rate", [(bse(fs1.rate / 100, fs1.n).max(), "3.4")],
       note="the bound holds; this recomputation gives a smaller maximum than the registry key figS1|mcse_max (3.39), so the drawn set used by the figure script is not exactly the one reconstructed here")
mc = pd.read_csv(T+"V5_S7g_m67_cc_dp.csv"); mc20 = mc[mc.n >= 20]
gm = mc20.groupby(["dgp", "scen", "arm"])
rbmax = gm.apply(lambda d: 1 / np.sqrt(d.loc[d.rb.abs().idxmax(), "n"]), include_groups=False)
covmin = gm.apply(lambda d: bse(d.loc[d["cov"].idxmin(), "cov"] / 100, d.loc[d["cov"].idxmin(), "n"]), include_groups=False)
capres("A0305", "tables/V5_S7g_m67_cc_dp.csv (n>=20)", "per row (mechanism x scenario x hypothesis): mean 1/sqrt(n); 1/sqrt(n) at the largest |bias/SE|; mean binomial SE of coverage; SE at the smallest coverage; mean binomial SE of the intervention share",
       [(gm.apply(lambda d: (1 / np.sqrt(d.n)).mean(), include_groups=False).max(), "0.036"), (rbmax.max(), "0.19"),
        (gm.apply(lambda d: bse(d["cov"] / 100, d.n).mean(), include_groups=False).max(), "0.74"), (covmin.max(), "0.89"),
        (gm.apply(lambda d: bse(d.share_B, d.n).mean(), include_groups=False).max(), "1.7")])
pp1 = pp[(pp.arm == "H1") & (pp.n_obs >= 20) & ~pp.method.isin(["glmer_bobyqa"])].copy()
pp1["grid"] = np.where(pp1.scen == "M1c", "M1c", "main")
g36 = pp1.assign(rbm=pp1.mcse_bias / pp1.empse).groupby(["grid", "dgp", "method"])
capres("A0326", "tables/variant_nc/V4_perf_primary.csv (H1, n_obs>=20)", "max over method x mechanism x grid of the mean per-point MCSE of bias/SE, SE ratio, coverage",
       [(g36.rbm.mean().max(), "0.036"), (g36.mcse_relerr.mean().max(), "5.9"), (g36.mcse_cover.mean().max() * (100 if pp1.mcse_cover.max() < 1 else 1), "1.0")],
       note="computed over all methods of the file except glmer_bobyqa (Table S36 rows)")
nl = pd.read_csv(T+"variant_nc/V4_S6_nestedloop.csv")
m5 = bse(0.05, 1000)
capres("A0162", "tables/variant_nc/V4_S6_nestedloop.csv (n_obs>=20)", "100 sqrt(.05*.95/1000); max mcse_rej among drawn points",
       [(m5, "0.69"), (nl[nl.n_obs >= 20].mcse_rej.max(), "6.4")],
       note="the 6.4 is the binomial MCSE at the observed rate of the point, not at 5%: the sentence reads as if it were the MCSE of a 5% rate (that would be at most 4.9 with n=20)", )

# ---------- Part 5: per-replicate facts, design, supplementary claims ----------
rb = pd.read_csv(os.path.join(RECS, "recs_by_dp.csv")); cc_ = pd.read_csv(os.path.join(RECS, "const_compare.csv")).iloc[0]
fg_ = pd.read_csv(os.path.join(RECS, "fgpt_check.csv"))
rb["grp"] = np.where(rb.scen.isin(["M6", "M7"]), "m67", np.where(rb.scen == "M1c", "m1c", "main"))
RS = "results_v4/S10/assolute/recs_by_dp.csv (from all 47,520 chunks of results_v4/S7d/rows and S7e/rows)"
ok = (rb.n_recs == 1000).all() and len(rb) == 4752
for cid in ["A0020", "A0133"]:
    rec(cid, "TRUE" if ok else "FALSE", RS, "records per design point x hypothesis", f"{len(rb)} dp x hypothesis rows, records per row {sorted(set(rb.n_recs))}")
dg_ = rb.groupby(["grp", "scen", "dgp", "label"]).n_deg.sum().reset_index()
main_ = dg_[dg_.grp == "main"].groupby(["dgp", "label"]).n_deg.sum()
m1c_ = dg_[dg_.grp == "m1c"].groupby(["dgp", "label"]).n_deg.sum()
m6 = dg_[dg_.scen == "M6"].groupby(["dgp", "label"]).n_deg.sum(); m7 = dg_[dg_.scen == "M7"].groupby(["dgp", "label"]).n_deg.sum()
pp_ = pp.merge(rb[["dgp", "scen", "N", "ntime", "label", "n_deg"]], on=["dgp", "scen", "N", "ntime", "label"], how="left")
den_ok = (pp_.n_den == 1000 - pp_.n_deg).all()
fdm = fd.merge(rb[rb.label == "manual"][["dgp", "scen", "N", "ntime", "n_deg"]], on=KEY, how="left")
nd_ok = (fdm.nd == 1000 - fdm.n_deg).all()
ok = set(main_) == {247} and set(m1c_) == {43} and set(m6) == {47} and set(m7) == {55} and den_ok and nd_ok
rec("A0039", "TRUE" if den_ok and nd_ok else "FALSE", "tables/variant_nc/V4_perf_primary.csv (n_den); tables/V5_S7g_filter_dp.csv (nd); " + RS,
    "denominator = every replicate with an estimable design (1000 - degenerate), whatever the fit returned", f"perf_primary cells {len(pp_)}: {den_ok}; rule rows {len(fdm)}: {nd_ok}")
rec("A0021", "TRUE" if ok else "FALSE", RS + "; tables/variant_nc/V4_perf_primary.csv (n_den); tables/V5_S7g_filter_dp.csv (nd)",
    "count design_degenerate per mechanism x hypothesis and grid; check every denominator = 1000 - degenerate",
    f"main grid per mechanism x hypothesis: {sorted(set(main_))} (incl. Copula-AR1); M1c {sorted(set(m1c_))}; M6 {sorted(set(m6))}; M7 {sorted(set(m7))}; perf n_den = 1000 - deg in all {len(pp_)} cells: {den_ok}; rule nd = 1000 - deg in all {len(fdm)} rows: {nd_ok}",
    note="'never observed after the first occasion' is the definition of design_degenerate; not re-derived from the data")
ar = rb[(rb.dgp == "GLMM-AR1")]
sh = ar[ar.n_gen < 1000]
n_sh = {g: sh[sh.grp == g][["scen", "N", "ntime"]].drop_duplicates().shape[0] for g in ["main", "m1c"]}
ok = n_sh == {"main": 44, "m1c": 6} and sh.n_gen.min() == 971 and sh.n_gen.max() == 999 and rb[rb.dgp != "GLMM-AR1"].n_gen.min() == 1000 and rb.n_notgen_nondeg.sum() == 0
rec("A0170", "TRUE" if ok else "FALSE", RS, "dps with fewer than 1,000 generated datasets, by mechanism and grid; range; non-generated datasets that are not flagged degenerate",
    f"short dps: {n_sh} (GLMM-AR1 only; other mechanisms min generated {rb[rb.dgp!='GLMM-AR1'].n_gen.min()}); generated range {sh.n_gen.min()}-{sh.n_gen.max()}; non-generated and not degenerate: {rb.n_notgen_nondeg.sum()}")
e = rb[rb.run == "S7e"]
ok_gee = rb.gee_missing_on_ok.sum() == 0 and rb.flip_seen_on_ok.sum() == rb.n_ok_design.sum()
miss = {c.replace("_missing_on_ok", ""): int(e[c].sum()) for c in e.columns if c.endswith("_missing_on_ok") and c.startswith("n_")}
rec("A0027", "TRUE" if ok_gee and all(v == 0 for k, v in miss.items() if k != "n_shrink") else "FALSE", RS,
    "every dataset with an estimable design has GEE rows for the 3 working correlations and a sign-flip record (S7d+S7e); in the Copula-AR1 run, a stored row of each study method",
    f"GEE rows missing on estimable datasets: {rb.gee_missing_on_ok.sum()}; flip record missing: {rb.n_ok_design.sum()-rb.flip_seen_on_ok.sum()}; Copula-AR1 study rows missing: {miss}",
    note="for the main grid, M1c and M6-M7 the study methods (MMM, GLMM, shrinkage, LRT) live in results_v4/variant_cc/S4/rows_v4.rds (5.2 GB) and results_v4/S6b/ns, not read here; in the Copula-AR1 run the shrinkage GLM has no stored row for 134,122 estimable datasets (28%), which the stored rows cannot distinguish from 'not applied' (they count as not usable)", severity="minor")
rec("A0032", "TRUE" if rb.flip_seen_on_ok.sum() == rb.n_ok_design.sum() else "FALSE", RS, "sign-flip record (p-value or reason) on every dataset with an estimable design",
    f"{rb.flip_seen_on_ok.sum()} of {rb.n_ok_design.sum()}; also on {rb.flip_on_notok.sum()} generated datasets with a degenerate design")
ok = cc_.recomputed == 9711 and cc_.both == 9711 and cc_.only_keys == 0 and cc_.all_datasets == 4752000 and round(100 * cc_.small_recomputed / cc_.recomputed) == 92
rec("A0037", "TRUE" if ok else "FALSE", "results_v4/S10/assolute/const_compare.csv vs results_v4/S9/const_keys.rds",
    "constant-outcome datasets among estimable designs recomputed from sum_y and n_obs; share with T<=3 and N<=25",
    f"recomputed {int(cc_.recomputed)}, in keys {int(cc_['keys'])}, both {int(cc_.both)}; datasets {int(cc_.all_datasets)}; small share {100*cc_.small_recomputed/cc_.recomputed:.1f}%; constant including degenerate designs {int(cc_.const_any_design)}",
    note="'no usable result for any method' is the rule (D-82) applied in the tables, not a property of the fits: before S9 uncorrected GEE, geesmv KC and the S7d GEE strategies returned 'usable' results on these datasets (e.g. 111 at GLMM-AR1 M2 N=10 T=2)")
f1 = fg_.set_index("dgp")
ok = (f1.loc["GLMM-AR1", ["n_pre", "k_pre", "n_post", "k_post"]].tolist() == [172, 12, 61, 12] and f1.loc["GLMM-native", ["n_pre", "k_pre", "n_post", "k_post"]].tolist() == [136, 11, 63, 11]
      and (f1.const_usable_rej == 0).all())
ncst = rb[(rb.scen == "M2") & (rb.N == 10) & (rb.ntime == 2) & (rb.label == "manual")].set_index("dgp").n_const
ok &= all(f1.loc[d, "const_usable"] == ncst[d] for d in ["GLMM-AR1", "GLMM-native"])
rec("A0270", "TRUE" if ok else "FALSE", "results_v4/S10/assolute/fgpt_check.csv (S7d rows, AR(1), FG, t with N-5 df, usability |est|<10, 0<SE<10)",
    "at the two dps: usable and rejections before/after removing constant datasets; constant datasets that were usable and their rejections",
    "; ".join(f"{d}: pre {f1.loc[d,'n_pre']}/{f1.loc[d,'k_pre']}, post {f1.loc[d,'n_post']}/{f1.loc[d,'k_post']}, constant {ncst[d]}, of which usable {f1.loc[d,'const_usable']} with {f1.loc[d,'const_usable_rej']} rejections" for d in ["GLMM-AR1", "GLMM-native"]))
# parameters: hypotheses differ only in the tested coefficient
dgm = pd.read_csv(T+"V4_dgm_parameters.csv")
a = dgm[dgm.arm == "H0"].merge(dgm[dgm.arm == "H1"], on=["dgp", "model", "block", "parameter"], suffixes=("_0", "_1"), how="outer")
dif = a[a.value_0.astype(str) != a.value_1.astype(str)]
def lastonly(r):
    v0 = [float(x) for x in str(r.value_0).split(",")]; v1 = [float(x) for x in str(r.value_1).split(",")]
    return len(v0) == len(v1) and all(abs(x - y) < 1e-12 for x, y in zip(v0[:-1], v1[:-1]))
dif = dif.assign(last_only=[lastonly(r) for r in dif.itertuples()])
bad = dif[~dif.last_only]
rec("A0015", "TRUE" if len(bad) == 0 else "DEFINITION_DIVERGENCE", "tables/V4_dgm_parameters.csv", "H0 vs H1 parameter rows: parameters that differ, and whether only the last (interaction) element differs",
    f"{len(dif)} differing parameter rows; differing beyond the interaction: " + "; ".join(f"{r.dgp} M{r.model} {r.parameter}" for r in bad.itertuples()),
    "GLMM-AR1 beta_c_calibrated M1: H0 (-2.499, 0.119, 0.236, 0.0018) vs H1 (-2.499, 0.119, 0.227, 0.255)",
    note="true for the coefficients stated in the paper (marginal targets of GLMM-AR1, conditional of GLMM-native, marginal of GenBinaryY); the conditional coefficients that actually generate GLMM-AR1 are calibrated per hypothesis and also differ in the treatment main effect (and slightly in intercept/time in M2, M5)", severity="minor")
rec("A0165", RES["A0015"][-1]["verdict"], "tables/V4_dgm_parameters.csv", "as A0015", RES["A0015"][-1]["recomputed"], RES["A0015"][-1]["counterexample"], note=RES["A0015"][-1]["note"], severity="minor")
pf = dgm[dgm.parameter.isin(["p_full_followup", "followup_rule"])].value.unique()
rec("A0016", "TRUE" if set(pf) == {"0.7", "ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))"} else "FALSE",
    "tables/V4_dgm_parameters.csv", "follow-up rule recorded for every mechanism, scenario and hypothesis", "; ".join(pf), note="recorded generator rule; generator code not re-read")
sr_ = dgm[dgm.parameter == "seed_rule"].value.unique()
rec("A0022", "TRUE" if list(sr_) == ["20260805L + 1000L * model + i"] else "FALSE", "tables/V4_dgm_parameters.csv (seed_rule)", "seed rule across mechanisms, scenarios and hypotheses",
    f"distinct rules: {list(sr_)}; it depends on scenario (model) and replicate only, not on N, T, mechanism or hypothesis",
    note="the per-dataset rows of S7d/S7e share sim_id across mechanisms; seeds themselves are not stored in the rows")
ra = dgm[dgm.block.isin(["random_effects", "association", "serial"])]
trt_re = ra[ra.parameter.isin(["re.formula", "lv.formula"]) & ra.value.str.contains("trt")]
trt_t = ra[(ra.parameter == "t.formula") & ra.value.str.contains("trt")][["dgp", "model"]].drop_duplicates()
rec("A0112", "TRUE" if len(trt_re) == 0 else "FALSE", "tables/V4_dgm_parameters.csv", "random-effect / latent-variable formulas containing trt",
    f"random-effect or latent-SD formulas with trt: {len(trt_re)}; transition formulas with trt: " + ", ".join(f"{r.dgp} M{r.model}" for r in trt_t.itertuples()),
    note="the GenBinaryY transition (serial association) term depends on treatment in M4 and M5 (t.formula ~trt, ~trt + steroids), so the within-subject association does differ by arm there; Table 2 describes M4 as 'As M3', which has no such term", severity="info")
tg = pd.read_csv(T+"V4_S2_truth_glmmnative.csv")
h0t = tg[tg.label == "manual"].truth_marginal_weighted; h1t = tg[tg.label == "mmlb"].truth_marginal_weighted
ok = (h0t != 0).all() and round(h0t.min(), 4) == -0.0030 and round(h0t.max(), 4) == -0.0011 and round(h1t.min(), 3) == 0.044 and round(h1t.max(), 2) == 0.20
rec("A0024", "TRUE" if ok else "FALSE", "tables/V4_S2_truth_glmmnative.csv (GLMM-native, M1-M5 x T)", "range of the marginal true value under H0 and H1; any exactly zero under H0",
    f"H0 {h0t.min():.4f} to {h0t.max():.4f} ({(h0t == 0).sum()} exactly zero of {len(h0t)}); H1 {h1t.min():.4f} to {h1t.max():.4f}", note="M1c range (0.19-0.20) read from tables/V5_S8b_keys.csv, not recomputed")
nanp = pp.p_usable.isna().sum(); nanr = rd.pu.isna().sum()
rec("A0038", "TRUE" if nanp == 0 and nanr == 0 else "FALSE", "tables/variant_nc/V4_perf_primary.csv; tables/V5_S7f_rule_dp.csv", "missing usable share over all method x cell rows",
    f"perf_primary: {nanp} missing of {len(pp)}; rule_dp: {nanr} missing of {len(rd)}")
thr = {}
for n in [1000, 500, 100]:
    k = next(k for k in range(n + 1) if cp_lower(k, n) > 0.075); thr[n] = 100 * k / n
ok = round(thr[1000], 1) == 9.3 and round(thr[500]) == 10 and round(thr[100]) == 14
for cid in ["A0052", "A0277"]:
    rec(cid, "TRUE" if ok else "FALSE", "Clopper-Pearson (scipy beta quantiles)", "smallest k/n with exact 95% lower limit > 7.5% for n = 1000, 500, 100", str({n: round(v, 2) for n, v in thr.items()}))
nl = pd.read_csv(T+"variant_nc/V4_S6_nestedloop.csv"); hid = nl[nl.n_obs < 20].groupby("method").size().to_dict()
for cid in ["A0059", "A0161"]:
    rec(cid, "TRUE", "tables/variant_nc/V4_S6_nestedloop.csv", "design points per method in the figure data; points not drawn (n_obs < 20)",
        f"{nl[KEY].drop_duplicates().shape[0]} dps x {nl.method.nunique()} methods; not drawn: {hid}",
        note="every dp is in the figure data, but 567 method x dp points (GLMM Wald, LRT, shrinkage) are not drawn; the Figure 2 caption says so, the sentence in Section 3.1 does not", severity="info")
fb = pd.read_csv(T+"V4_S6b_ns_fallback.csv"); f0 = fb[(fb.arm == "H0") & (fb.band != "all")]
top = f0.loc[f0.groupby(["dgp", "scen"]).nonconv.idxmax()]
rec("A0065", "TRUE" if (top.band == "T2-3").all() else "FALSE", "tables/V4_S6b_ns_fallback.csv (H0)", "band of occasions with the highest nlm non-convergence share, per mechanism x scenario",
    "; ".join(f"{r.dgp} {r.scen}: {r.band} {100*r.nonconv:.0f}%" for r in top.itertuples()))
s6 = pd.read_csv(T+"V5_S7f_cp_S6cc_cells.csv")
rec("A0108", "TRUE" if s6[s6.arm == "H0"].hi_c.sum() == 0 else "FALSE", "tables/V5_S7f_cp_S6cc_cells.csv (M3-M5, H0, normal reference)", "null dps clearly above 7.5% with and without the check",
    f"without {int(s6[s6.arm=='H0'].hi_s.sum())}, with {int(s6[s6.arm=='H0'].hi_c.sum())} of {int((s6.arm=='H0').sum())}")
sr = pd.read_csv(T+"V5_S7f_seratio_fpc.csv"); full = sr[(sr.spec == "pieno") & sr.method.isin(["GLMM", "MMM"])]; red = sr[(sr.spec == "ridotto") & sr.method.isin(["GLMM", "MMM"])]
rec("A0110", "TRUE" if (full.ratio_fpc < 1).all() else "FALSE", "tables/V5_S7f_seratio_fpc.csv (full mean model, as Table S30)", "corrected SE ratio < 1 for GLMM and MMM at every N in both applications",
    "full model: " + "; ".join(f"{a} {m} {g.ratio_fpc.min():.2f}-{g.ratio_fpc.max():.2f}" for (a, m), g in full.groupby(["application", "method"])) + f"; reduced model max {red.ratio_fpc.max():.2f}",
    note="with the reduced mean model of the respiratory trial the GLMM ratio is 1.07 at N=10 and 1.02 at N=15 (not short); the sentence cites Table S30 (full model)", severity="info")
s37 = pd.read_csv(T+"V5_S8b_sub_s37.csv").set_index("strategy")
f_m, f_g = 100 - s37.loc["mmm_cc", "p_pass_boot"], 100 - s37.loc["gee_ar1_fg", "p_pass_boot"]
for cid in ["A0114", "A0286"]:
    rec(cid, "TRUE" if round(f_m) == 82 and round(f_g) == 99 else "FALSE", "tables/V5_S8b_sub_s37.csv", "share of bootstrap resamples with at least one dp clearly above 7.5%", f"MMM {f_m:.2f}%, GEE {f_g:.2f}%")
mins = {c: s37[c].idxmin() for c in ["mean_flags_boot", "plugin_exp_flags", "eb_exp_flags"]}
ok = mins["mean_flags_boot"] == "mmm_cc" and mins["plugin_exp_flags"] == "mmm_cc" and mins["eb_exp_flags"] == "gee_ar1_fg"
for cid in ["A0115", "A0287", "A0288"]:
    rec(cid, "TRUE" if ok else "FALSE", "tables/V5_S8b_sub_s37.csv (8 strategies)", "argmin over the 8 strategies of expected flags (bootstrap, plug-in, shrunk plug-in)",
        "; ".join(f"{c}: {s37[c].idxmin()} {s37[c].min():.2f} (next {s37[c].nsmallest(2).iloc[1]:.2f})" for c in mins))
nf = pd.read_csv(T+"V5_S7f_cp_null_flips.csv").value.iloc[0]
rec("A0117", "TRUE" if nf == 0 else "FALSE", "tables/V5_S7f_cp_null_flips.csv; tables/V5_S7f_extra.csv", "flags changed after subtracting the null excess",
    f"flips {nf}", note="the bounds 0.029 and 0.26 are read from tables/V5_S7f_extra.csv, not recomputed here")
gt = pd.read_csv(T+"geesmv_timeouts.csv")
rec("A0210", "TRUE" if round(gt.pct_timeout.max()) == 17 else "FALSE", "tables/geesmv_timeouts.csv", "max share of time-outs over cells", f"{gt.pct_timeout.max()}% ({gt.loc[gt.pct_timeout.idxmax(), ['dgp','N','ntime','method']].to_dict()})")
ge = pd.read_csv(T+"V4_S7d_GE4.csv")
mx = max(ge[["max_md", "max_kc", "max_fg"]].max()); both_min = ge[["both_md", "both_kc", "both_fg"]].min().min(); one = ge[["one_md", "one_kc", "one_fg"]].sum().sum()
rec("A0212", "FALSE" if both_min < 670 else "TRUE", "tables/V4_S7d_GE4.csv", "datasets compared per working correlation, datasets where both computations returned, max |difference|",
    f"n {ge.n.min()}-{ge.n.max()}; both returned min {both_min}; one side only {one}; max |d| {mx:.2e}",
    "exchangeable MD: SEs agree on 666 datasets (4 where only one computation returned); AR(1) KC 667 (4); independence KC 669 (2)",
    note="max difference < 1e-12 holds; '670 or more datasets' counts datasets compared, not datasets on which both computations returned an SE", severity="minor")
m1 = pd.read_csv(T+"variant_nc/V4_S6_m1c_claims.csv")
c10 = m1[m1.name.str.startswith("C10b")]
rec("A0217", "TRUE" if (c10["class"] == "non risolta").all() else "FALSE", "tables/variant_nc/V4_S6_m1c_claims.csv (C10b)", "class of the difference in miscalibrated-cell share between the two best methods in M1c",
    "; ".join("%s: %.2f [%.2f, %.2f] %s" % (r["name"], r.est, r.lo, r.hi, r["class"]) for _, r in c10.iterrows()))
fi = pd.read_csv(T+"V5_S8b_flip_interval.csv")
rec("A0243", "TRUE" if (fi.hi_t == 0).all() and (fi.lo_t == 0).all() and (fi.points == 2376).all() else "FALSE", "tables/V5_S8b_flip_interval.csv", "test at the true value: dps with non-coverage clearly outside Bradley's band; mean coverage",
    f"hi {fi.hi_t.tolist()}, lo {fi.lo_t.tolist()} over {fi.points.tolist()} dps; mean coverage {fi.mean_cov_t.round(2).tolist()}")
lpz = pd.read_csv(T+"V5_S8b_flip_lo_points.csv")
rec("A0244", "TRUE", "tables/V5_S8b_flip_lo_points.csv; tables/V4_S7d2_explain_finite.csv", "share of unbounded intervals at the clearly-below points, and at T=2 overall",
    f"clearly-below points: median share unbounded {lpz.median_share_unb.round(1).tolist()}% (min {lpz.min_share_unb.round(1).tolist()}); all T=2 datasets: finite {100*0+ (ef[ef.ntime==2].finite*ef[ef.ntime==2].datasets).sum()/ef[ef.ntime==2].datasets.sum():.1f}%",
    note="true if 'where' refers to the 90 design points (unbounded in 86-100% of datasets); at two occasions in general the interval is finite in 42% of datasets (72% with 75-145 clusters)", severity="minor")
pt = pd.read_csv(T+"V4_S7e_glmm_pseudotrue.csv"); z = (pt.est_timetrt / pt.se_timetrt).abs()
p1_ = pt[pt.model == 1]
rec("A0255", "TRUE" if (z < 1.96).all() else "FALSE", "tables/V4_S7e_glmm_pseudotrue.csv", "|estimate/SE| < 1.96 for every fitted large dataset; range of estimates and SEs (random intercept, model 1)",
    f"max |z| {z.max():.2f}; model 1 estimates {p1_.est_timetrt.min():.4f}-{p1_.est_timetrt.max():.3f}, SE max {p1_.se_timetrt.max():.3f}; file also has 2 rows of scenario M3 (SE up to {pt.se_timetrt.max():.3f}) not described in the text", severity="info")
cdp = rd[rd.grp == "copula"][KEY].drop_duplicates()
rec("A0256", "TRUE" if len(cdp) == 480 and set(cdp.scen) == {"M1", "M2", "M3", "M4", "M5"} else "FALSE", "tables/V5_S7f_rule_dp.csv", "Copula-AR1 dps and scenarios", f"{len(cdp)} dps, scenarios {sorted(set(cdp.scen))}")
fa = f10[f10.assess].assign(up=[cp_upper(k, n) for k, n in zip(f10[f10.assess].k, f10[f10.assess].n)])
mxr = fa.groupby("strategy").rate.max(); ne = fa[(fa.up > 0.075) & ~fa.hi].groupby("strategy").size()
ok = round(mxr["gee_ar1_fg"]) == 20 and round(mxr["mmm_cc"], 1) == 8.4 and ne["gee_ar1_fg"] == 227 and ne["mmm_cc"] == 438
rec("A0281", "TRUE" if ok else "FALSE", "tables/V5_S7g_filter_dp.csv (L=10, assessed dps)", "highest rate and count of dps with CP upper > 7.5% but not clearly above",
    f"max rate GEE {mxr['gee_ar1_fg']:.2f}, MMM {mxr['mmm_cc']:.2f}; not excluded GEE {ne['gee_ar1_fg']}, MMM {ne['mmm_cc']}")
agg_ = pd.read_csv(T+"V5_S7g_agq_gates.csv").set_index("gate")
rec("A0294", "TRUE", "tables/V5_S7g_agq_gates.csv (GG2)", "gate detail", agg_.loc["GG2", "detail"],
    note="the 122,904 refits were compared with 122,688 stored rows: 216 refits have no stored counterpart, so 'reproduced ... exactly' holds for the 122,688 comparable rows", severity="minor")
ag = pd.read_csv(T+"V5_S7g_agq_dp.csv"); a25 = ag[ag.nAGQ == 25]
rec("A0295", "TRUE" if len(a25) == 123 and a25.nd.max() <= 100 else "FALSE", "tables/V5_S7g_agq_dp.csv (nAGQ=25)", "dps with the 25-node check and datasets per dp",
    f"{len(a25)} dps, estimable datasets per dp {a25.nd.min()}-{a25.nd.max()}", note="the 10 of 9,844 decision changes are read from tables/V5_S7g_agq_keys.csv")
mc = pd.read_csv(T+"V5_S7g_m67_cc_dp.csv"); mc20 = mc[mc.n >= 20].copy()
mc20["ncov_lo"] = [cp_lower(n - k, n) for k, n in zip(mc20.kcov, mc20.n)]
rbmax = mc20.rb.abs().max(); nhi = int((mc20.ncov_lo > 0.075).sum())
rec("A0299", "TRUE" if round(rbmax, 2) <= 0.40 and nhi == 0 else "FALSE", "tables/V5_S7g_m67_cc_dp.csv (n>=20)", "max |bias/SE| over dps x hypotheses; dps with non-coverage CP lower > 7.5%",
    f"max |bias/SE| {rbmax:.4f}; non-coverage clearly above 7.5%: {nhi} of {len(mc20)}")
ns_ = pd.read_csv(T+"V4_S6b_ns_rule_safety_dp.csv")
hid = {"gee_" + s: int(((c7.scen.isin(['M6', 'M7'])) & (c7.strategy == "gee_" + s) & (c7.ref == "t") & (c7.n < 20)).sum()) for s in ["exch_md", "exch_fg", "ar1_fg"]}
hid["mmm_cc"] = int(((ns_.strategy == "mmm|t") & (ns_.n < 20)).sum()); hid["glmer_lrt"] = int(((ns_.strategy == "glmer_lrt|n") & (ns_.n < 20)).sum())
rec("A0301", "TRUE" if hid["glmer_lrt"] == 60 and sum(v for k, v in hid.items() if k != "glmer_lrt") == 0 else "FALSE", "tables/V5_S7f_cp_S7d_dp_rates.csv (M6-M7, t); tables/V4_S6b_ns_rule_safety_dp.csv (MMM, LRT)",
    "dps with fewer than 20 usable replicates per strategy of Figure S1", str(hid))
m67 = rd[rd.grp == "m67"]; hm = m67[m67.hi].groupby(["strategy", "dgp"]).size()
want = {"gee_exch_md": 18, "gee_exch_fg": 4}
ok = all(m67[(m67.strategy == s) & m67.hi].shape[0] == v for s, v in want.items()) and set(m67[m67.strategy.isin(list(want)) & m67.hi].dgp) == {"GenBinaryY"} and \
     all(m67[(m67.strategy == s) & m67.hi].shape[0] == 0 for s in ["gee_ar1_fg", "mmm_cc", "glmer_lrt"])
rec("A0302", "TRUE" if ok else "FALSE", "tables/V5_S7f_rule_dp.csv (M6-M7)", "hi per strategy and mechanism in M6-M7", hm.to_dict())
thr_ = pd.read_csv(T+"V4_S6b_thr_rule.csv"); td = pd.read_csv(T+"V4_S6b_thr_claims_diff.csv")
rec("A0315", "UNCHECKABLE", "tables/V4_S6b_thr_claims_diff.csv (2026-10-19)", "rows of the claims-difference table between thresholds",
    f"{len(td)} changed readings in the file", note="the deciding table predates the constant-outcome rule (S9, 2026-11-02) and the fewclust KC; the 156 comparisons were not re-evaluated with the thresholds after S9, so the claim is verified only on the pre-S9 tables", severity="minor")
x = thr_[thr_.strategy == "mmm|t"]
rec("A0316", "TRUE" if (x.n_hi == 0).all() else "FALSE", "tables/V4_S6b_thr_rule.csv (2026-10-19)", "hi of mmm|t with thresholds 1e-4 (t4) and 1e-2 (t2)",
    "; ".join(f"{r.variant}: {r.n_hi} of {r.n_dp}" for r in x.itertuples()),
    note="'any design point' covers the 1,512 dps of the main grid and M1c only (not M6-M7 or Copula-AR1); table predates S9, but S9 changed no MMM result (tables/V6_S9_nc_usable_change.csv, dU0 = 0 for mmm)", severity="minor")
P, K_, n_ = 207, 25, 10
ph = math.comb(P - K_, n_) / math.comb(P, n_)
rec("A0318", "TRUE" if round(ph, 2) == 0.27 else "FALSE", "hypergeometric with P=207 (tables/V5_S7f_seratio_fpc.csv), 25 ventilated", "P(no ventilated patient in a draw of 10)", f"{ph:.4f}")
ld = pd.read_csv(T+"V5_S7g_lus_draws.csv")
zo = ld[ld.N != 10].z.abs().max()
rec("A0319", "TRUE" if round(abs(ld[ld.N == 10].z.iloc[0]), 1) == 2.5 and round(zo, 1) <= 1.1 else "FALSE", "tables/V5_S7g_lus_draws.csv", "z at N=10 and max |z| at the other sizes",
    f"N=10 z {ld[ld.N==10].z.iloc[0]:.2f}; others max |z| {zo:.2f}")
c36 = pp1.groupby(["grid", "dgp", "method"]).apply(lambda d: d[["scen", "N", "ntime"]].drop_duplicates().shape[0], include_groups=False)
mm = c36.loc["main"]; m1 = c36.loc["M1c"]
rec("A0324", "TRUE" if mm.min() == 416 and (m1 == 24).all() else "FALSE", "tables/variant_nc/V4_perf_primary.csv (H1, n_obs>=20)", "dps with at least 20 usable replicates per method x mechanism (main grid of 480, M1c of 24)",
    f"main min {mm.min()} ({mm.idxmin()}); M1c min {m1.min()} ({m1.idxmin()})")

# ---------- Part 6: remaining items (literature, definitions, process) ----------
LIT = "literature claim about the cited studies; not decidable from the study tables (would need the cited papers)"
rec("A0007", "UNCHECKABLE", "Section S14 text", "internal consistency with Section S14 only", "S14 lists the units of [7],[9],[15],[19],[20]: site within mouth, patient within practice, participant within cluster, youth within community; none is a subject measured on several occasions",
    note=LIT + "; consistent with Section S14", severity="info")
rec("A0308", "UNCHECKABLE", "Section S14 text", "internal consistency with the list of units in the same paragraph", "consistent", note=LIT, severity="info")
rec("A0008", "UNCHECKABLE", "-", "-", "-", note=LIT + "; 'these studies' has no explicit list in the text (the preceding sentence cites [7,9,15,19,20]); Section S14 does not say which studies fit random slopes", severity="minor")
rec("A0018", "DEFINITION_DIVERGENCE", "Table 1/2 of the paper; method definitions", "reading of 'assumes'",
    "GEE with AR(1) or independence working correlation assumes only the marginal mean model, which Copula-AR1 satisfies exactly (logistic margins, Section S8); the sign-flip test assumes no dependence model either",
    note="true only if 'assumes this mechanism' means 'is a likelihood model of it' (MMM, GLMM); the GEEs and the sign-flip test are valid under Copula-AR1 by construction, which is the point of Section 3.6. Suggest 'None of the likelihood-based methods compared models this mechanism'", severity="minor")
rec("A0055", "DEFINITION_DIVERGENCE", "Section 2.4, Section S7", "LRT and sign-flip: interval available?",
    "the sign-flip test has an interval obtained by inverting it (Section S7, Table S27), excluded from the rule; the LRT has none",
    note="'give no interval' is true of the rule as defined (Wald-type intervals), not literally of the sign-flip test", severity="info")
sm = rd[rd.small]; g = sm.groupby("strategy").pu.mean()
rec("A0105", "TRUE", "tables/V5_S7f_rule_dp.csv (216 small dps)", "flip usable share, hi/lo in small dps; lower bound of the share of datasets where flip answers and the GEE strategy does not",
    f"flip {g['flip']:.1f}%, GEE {g['gee_ar1_fg']:.1f}%, MMM {g['mmm_cc']:.1f}%: flip answers without a GEE result in at least {g['flip']-g['gee_ar1_fg']:.1f}% of datasets; flip hi {int(sm[sm.strategy=='flip'].hi.sum())}, lo {int(sm[sm.strategy=='flip'].lo.sum())}",
    note="'where the estimators return nothing' is supported as a bound from the dp shares; the joint per-dataset overlap with the MMM was not computed", severity="info")
rec("A0130", "TRUE", "manuscript/v5/MMM_v5.qmd", "numerals in prose outside `r n()` calls, after removing references, citations, cross-references, scenario and mechanism labels",
    "only design constants and labels remain (7.5%, 2.5%, 5%, 95%, cluster counts and occasions such as N = 10, T = 2)",
    note="'every number' holds for computed results; design constants are literal", severity="info")
rec("A0163", "TRUE", "manuscript/v5/MMM_v5_supplement.qmd", "as A0130", "only design constants remain, plus the literal threshold 37 (|estimate/SE| > 37, Section S3) and the figure width 100%",
    note="the threshold 37 is a fixed constant of the analysis written as a literal", severity="info")
t3top = pd.read_csv(T+"V4_S6_table_cells.csv")
kc_old = t3top[(t3top.table_id == "T3") & (t3top.i == 8)].est.values; kc_new = tc[(tc.table_id == "T3") & (tc.i == 8)].est.values
for cid in ["A0164", "A0215"]:
    rec(cid, "UNCHECKABLE", "tables/variant_nc/V4_S6_table_cells.csv vs tables/V4_S6_table_cells.csv (pre-S9, geesmv KC)", "partial: KC row of Table 3 differs from the geesmv version",
        f"Table 3 KC row: max |fewclust - geesmv| {np.abs(kc_new - kc_old).max():.2f} points (fewclust in use); per-table provenance of KC for all other tables not traced",
        note="verified for Table 3 and Table 5 (strategies gee_*_kc come from fewclust per tables/V5_S7f_claims.csv impl_kc); not traced table by table for the supplement", severity="info")
rec("A0173", "TRUE", "results_v4/S10/assolute/recs_by_dp.csv", "2,376 x 2 x 1,000 and records read",
    f"{2376*2*1000} = records read {int(rb.n_recs.sum())}", note="'52.1 million model fits' and 'every analysis added after the main run regenerated them from their seeds' were not recounted (process claims)", severity="info")
rec("A0180", "TRUE", "rendered Table 3 and Table S3; tables/variant_nc/V4_S6_table_cells.csv; tables/V4_S6_table_cells.csv", "row-by-row comparison of the rendered tables; source of each row",
    "rows 1-2 and 4-8 identical in the rendering and read from the same variant_nc cells; the MMM row of Table S3 comes from tables/V4_S6_table_cells.csv (2026-10-18)",
    note="the MMM-without-check row is read from a table that predates the constant-outcome rule; S9 changed no MMM result (V6_S9_nc_usable_change.csv, dU0 = 0 for mmm), so the row is still current", severity="info")
n_ok = int((rd[rd.strategy == "glmer_lrt"].pu > 0).sum())
rec("A0159", "TRUE" if n_ok == 2094 else "FALSE", "tables/V5_S7f_rule_dp.csv", "dps where the LRT has a positive usable share", f"{n_ok} of 2376")
rbias = pd.read_csv(T+"V5_S7f_resp_bias.csv"); ru = pd.read_csv(T+"V5_S7f_resp_usable.csv"); srr = pd.read_csv(T+"V5_S7f_seratio_fpc.csv")
rf, uf, sf = rbias[rbias.spec == "pieno"], ru[ru.spec == "pieno"], srr[(srr.application == "respiratory") & (srr.spec == "pieno")]
capres("A0266", "tables/V5_S7f_resp_usable.csv, V5_S7f_resp_bias.csv, V5_S7f_seratio_fpc.csv (respiratory, full model)", "max MCSE of usable share, bias and SE ratio",
       [(uf.mcse_usable_500.max(), "2.2"), (rf.mcse_bias.max(), "0.21"), (sf.mcse_ratio.max(), "0.15")])

# ---------- write the verdict table ----------
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s10_assolute_classify as C
cand = pd.read_csv(os.path.join(RECS, "candidates_mapped.csv"))
rows = []
for r in cand.itertuples():
    base = dict(cand_id=r.cand_id, sentence=r.sentence, file="manuscript/v5/" + r.file, render_line=r.line_start, render_line_end=r.line_end,
                section=r.section, words=r.words, claim_marker=bool(r.has_claim_marker), marker=r.marker if isinstance(r.marker, str) else "",
                qmd_line=r.qmd_line, registry_ids=r.reg_ids if isinstance(r.reg_ids, str) else "")
    if r.cand_id in C.DISCARD:
        rows.append({**base, "table": "", "computation": "", "verdict": "NOT_A_RESULT_CLAIM", "definition_divergence": False, "recomputed": "",
                     "counterexample": "", "discard_reason": C.DISCARD[r.cand_id], "note": "", "severity": ""})
        continue
    rs = RES.get(r.cand_id)
    if not rs:
        rows.append({**base, "table": "", "computation": "", "verdict": "UNCHECKABLE", "definition_divergence": False, "recomputed": "", "counterexample": "",
                     "discard_reason": "", "note": "no computation written", "severity": "minor"}); continue
    v = rs[-1]
    verdict = v["verdict"]; dd = verdict == "DEFINITION_DIVERGENCE"
    if dd: verdict = "TRUE"
    rows.append({**base, "table": v["table"], "computation": v["computation"], "verdict": verdict, "definition_divergence": dd, "recomputed": v["recomputed"],
                 "counterexample": v["counterexample"], "discard_reason": "", "note": v["note"], "severity": v["severity"]})
out = pd.DataFrame(rows)
os.makedirs("results_v4/S10/assolute", exist_ok=True)
out.to_csv(os.environ.get("S10A_OUT", "tables/V7_S10_assolute_check.csv"), index=False)
print(out.verdict.value_counts().to_string()); print("definition divergences:", int(out.definition_divergence.sum()))
