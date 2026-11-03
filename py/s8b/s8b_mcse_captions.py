"""s8b_mcse_captions.py -- S8b sub-task, plan item DEC-T5 (technical decision): one Monte Carlo standard error (MCSE)
value per quantity for every table or figure whose caption lacks one (list: reports/morris/S8_logs/H2_mcse_by_table.csv).
Four kinds of value, written in the column `how` of each key:
  (a) registry: maximum of the registry column `mcse` (block bootstrap, B = 2000) over the registry ids used in the
      columns of that quantity, ids parsed from the table block of the qmd (read-only);
  (b) bound: for an average over design points of binomial proportions, mean over the design points of the per-point
      binomial MCSE 100 sqrt(p (1 - p) / n); since sd(mean) <= mean(sd) whatever the correlation between design points,
      this is an upper bound for the MCSE of the average; the caption gives its maximum over the rows of the table;
  (c) binomial: a single binomial share (largest rate at a point, shares over subsamples or datasets): 100 sqrt(p (1-p)/n),
      maximum over the table;
  (d) bootstrap: the block-bootstrap mcse column of tables/V5_S7f_h2h_boot.csv (Table S31), maximum per quantity.
  (e) approx: Table S35 measures that are not proportions: per-point delta-method MCSE with the empirical SE treated as
      known (bias / SE: 1 / sqrt(n); SE ratio: (1 + ratio) / sqrt(2 (n - 1)), empirical-SE component only).
Every recomputed point estimate is checked against the value the table prints (gate), so the MCSE refers to the same
design points and denominators. No bootstrap is run here.
Inputs read-only. Outputs: tables/V5_S8b_sub_cap_keys.csv (keys cap|<table>|<quantity>), tables/V5_S8b_sub_cap_checks.csv.
Run from ~/MMM: python3 py/s8b/s8b_mcse_captions.py
"""
import os, re, json, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/MMM"))
REG = pd.read_csv("tables/V5_S7f_numbers.csv", dtype=str, keep_default_na=False)
RI = REG.set_index("id")
TOK = re.compile(r'`r n\("([^"]+)"\)`')
KEYS, CHK = [], []
def chk(name, ok, detail=""):
    CHK.append((name, bool(ok), str(detail))); print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
def ceil_sig2(x):   # S11 (finding F15): an upper bound is rounded up to the two significant digits it is printed with
    import math
    x = float(x)
    if not math.isfinite(x) or x <= 0: return x
    e = math.floor(math.log10(x)) - 1
    return math.ceil(round(x / 10 ** e, 9)) * 10 ** e
def put(key, value, kind_mc, source, how, kind="sig2"):
    if kind_mc == "bound" and kind == "sig2": value = ceil_sig2(value)   # S11 (F15)
    assert np.isfinite(value), key
    KEYS.append((key, str(int(value)) if kind == "int" else repr(float(value)), "", kind, source, f"{kind_mc}: {how}"))
def b(p, n):   # binomial MCSE in percentage points, p in %
    p = np.asarray(p, float) / 100; n = np.asarray(n, float)
    return 100 * np.sqrt(p * (1 - p) / n)
def num(x): return pd.to_numeric(pd.Series(x), errors="coerce").to_numpy()
# ---- table blocks ----------------------------------------------------------------------------------------------------
def blocks(path):
    L = open(path, encoding="utf-8").read().split("\n"); out = {}; cur = None
    for l in L:
        m = re.match(r"\*\*(Table S?\d+|Figure S?\d+)\.\*\*", l)
        if m: cur = m.group(1); out[cur] = []; continue
        if cur and l.startswith("|"): out[cur].append(l)
        elif cur and out[cur] and not l.startswith("|"): cur = None
    return out
def split_row(r):
    ids = []
    s = TOK.sub(lambda m: (ids.append(m.group(1)), f"@@{len(ids) - 1}@@")[1], r)
    return [(re.sub(r"@@(\d+)@@", lambda m: ids[int(m.group(1))], c).strip(), [ids[int(x)] for x in re.findall(r"@@(\d+)@@", c)])
            for c in s.strip().strip("|").split("|")]
def cells(rows):
    rows = [r for r in rows if not re.match(r"^\|[\s\-:|]+$", r)]
    hdr = [c for c, _ in split_row(rows[0])]; out = []
    for r in rows[1:]:
        cs = split_row(r)
        for j, (c, ii) in enumerate(cs):
            for i in ii: out.append((hdr[j] if j < len(hdr) else f"col{j}", i))
    D = pd.DataFrame(out, columns=["col", "id"])
    D = D[~D.id.str.startswith("f:s9|")]  # S8b: rows built from this script's own keys (MD/KC average) are not inputs
    D = D.merge(REG[["id", "file", "filter", "field", "kind", "mcse", "value"]], on="id", how="left", validate="many_to_one")
    assert D.file.notna().all(), "registry id missing"
    D["m"] = num(D.mcse.replace("", np.nan)); D["v"] = num(D.value)
    D["flt"] = [json.loads(f) if f else {} for f in D["filter"]]
    return D
SB = blocks("manuscript/v5/MMM_v5_supplement.qmd"); MB = blocks("manuscript/v5/MMM_v5.qmd")
# ---- (a) registry maxima ---------------------------------------------------------------------------------------------
GROUPS = {  # table: {quantity: regex on the column header}
    "S3": {"rate": r"^Rate", "worst": r"^Worst", "mad": r"^MAD"}, "S16": {"rate": r"^Rate", "worst": r"^Worst", "mad": r"^MAD"},
    "S4": {"mad": r"."}, "S17": {"mad": r"."}, "S5": {"share": r"."}, "S18": {"share": r"."},
    "S6": {"rate": r"^(Own|Common)", "shift": r"^Shift"}, "S19": {"rate": r"^(Own|Common)", "shift": r"^Shift"},
    "S7": {"usable": r"."}, "S20": {"usable": r"."}, "S8": {"usable": r"^Usable", "uncond": r"^Unconditional"},
    "S21": {"usable": r"^Usable", "uncond": r"^Unconditional"}, "S9": {"mad": r"."}, "S22": {"mad": r"."},
    "S10": {"raw": r"^Raw", "sa": r"^Size"}, "S23": {"raw": r"^Raw", "sa": r"^Size"},
    "S11": {"rate": r"^(Normal|t) (GB|NAT|AR1)$", "shift": r"minus"}, "S24": {"rate": r"^(Normal|t) (GB|NAT|AR1)$", "shift": r"minus"},
    "S12": {"cov": r"."}, "S25": {"cov": r"."}, "S30": {"usable": r"^Usable", "bias": r"^Bias"}}
for t, G in GROUPS.items():
    D = cells(SB[f"Table {t}"]); D = D[D.kind != "int"]
    for q, rx in G.items():
        d = D[D.col.str.contains(rx, regex=True)]
        chk(f"{t} {q}: every registry id carries an mcse", len(d) > 0 and d.m.notna().all(), f"{d.m.notna().sum()}/{len(d)} ids")
        put(f"cap|{t}|{q}", d.m.max(), "registry", "tables/V5_S7f_numbers.csv",
            f"max of column mcse over the {len(d)} registry ids of columns {sorted(set(d.col))} of Table {t} (block bootstrap)")
# Table S30, SE ratio: mcse_ratio of tables/V5_S7f_seratio_fpc.csv (rows of the registry ids)
D = cells(SB["Table S30"]); d = D[D.col.str.startswith("SE ratio")]
SR = pd.read_csv("tables/V5_S7f_seratio_fpc.csv", dtype=str)
mm = []
for f in d.flt:
    r = SR
    for k, v in f.items(): r = r[r[k] == str(v)]
    assert len(r) == 1; mm.append(float(r.mcse_ratio.iloc[0]))
chk("S30 SE ratio: one row of V5_S7f_seratio_fpc.csv per id", len(mm) == len(d), len(mm))
put("cap|S30|seratio", max(mm), "registry", "tables/V5_S7f_seratio_fpc.csv", f"max of mcse_ratio over the {len(mm)} rows of the respiratory trial in Table S30 (ratio before the finite-population factor, which is below 1, so conservative)")
# ---- (d) Table S31: block bootstrap ------------------------------------------------------------------------------------
HB = pd.read_csv("tables/V5_S7f_h2h_boot.csv", dtype={"region": str, "q": str})
D = cells(SB["Table S31"]); D = D[D.field == "est"]
D["region"] = [f.get("region") for f in D.flt]; D["q"] = [f.get("q") for f in D.flt]
D = D.merge(HB[["region", "q", "mcse", "est"]].rename(columns={"mcse": "mc_hb"}), on=["region", "q"], how="left")
chk("S31: every id found in V5_S7f_h2h_boot.csv", D.mc_hb.notna().all(), len(D))
for q, rx in {"rate": r"^rate_", "mad": r"^mad_(gee|mmm)$", "diff": r"minus", "pu": r"^pu_", "upow": r"^upow_"}.items():
    d = D[D.q.str.contains(rx, regex=True)]
    put(f"cap|S31|{q}", d.mc_hb.max(), "bootstrap", "tables/V5_S7f_h2h_boot.csv", f"max of the block-bootstrap mcse over the {len(d)} cells ({q})")
# ---- per-point files for Table 4, S26, S28 ------------------------------------------------------------------------------
KP = ["dgp", "scen", "N", "ntime"]
A = pd.read_csv("tables/V5_S7f_cp_S7d_dp_rates.csv"); Bc = pd.read_csv("tables/V5_S7f_cp_S7e_dp_rates.csv")
C = pd.read_csv("tables/V5_S7f_cp_S7b_rule_dp.csv"); FL = pd.read_csv("tables/V5_S7f_cp_S7d2_flip_rules_dp.csv")
EST = ["gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_ar1_md", "gee_ind_md", "gee_exch_md", "gee_ar1_kc", "gee_ind_kc", "gee_exch_kc",
       "gee_ind_cr2", "gee_exch_cr2", "gee_ar1_unc", "gee_ind_unc", "gee_exch_unc"]
FLv = FL[FL.version == "with D-57 rules"]
ND1 = FLv[FLv.label == "mmlb"][KP + ["nd"]].rename(columns={"nd": "nd1"})      # estimable designs under H1 (same for every strategy)
keep = KP + ["strategy", "n", "nd", "pu", "upow", "rate", "mc"]
s1 = A[(A.ref == "t") & A.strategy.isin(EST + ["study_mmm"])][keep].replace({"strategy": {"study_mmm": "mmm_cc"}})
s2 = C[((C.ref == "t") & C.method.isin(["glmer", "shrink"])) | ((C.ref == "n") & (C.method == "glmer_lrt"))].assign(strategy=lambda x: x.method)[keep]
s3 = Bc[((Bc.ref == "t") & Bc.strategy.isin(EST + ["mmm_cc", "mmm_study", "glmer", "shrink"])) | ((Bc.ref == "n") & (Bc.strategy == "glmer_lrt"))][keep]
f0 = FLv[FLv.label == "manual"][KP + ["nd", "n", "pu", "rate", "mcse"]].rename(columns={"mcse": "mc"})
f1 = FLv[FLv.label == "mmlb"][KP + ["pu", "rate"]].assign(upow=lambda x: x.pu * x.rate / 100)[KP + ["upow"]]
fl = f0.merge(f1, on=KP).assign(strategy="flip")[keep]
DP = pd.concat([s1, s2, s3, fl], ignore_index=True).merge(ND1, on=KP, how="left", validate="many_to_one")
chk("H1 estimable-design count found for every point", DP.nd1.notna().all())
DP["small"] = (DP.ntime <= 3) & (DP.N <= 25)
DP["mc_pu"] = b(DP.pu, DP.nd); DP["mc_upow"] = b(DP.upow, DP.nd1)
chk("pu = 100 n / nd on every row", np.allclose(DP.pu, 100 * DP.n / DP.nd, atol=1e-9))
chk("2,376 points per strategy (480 Copula-AR1 points for the model without the check)", (DP[DP.strategy != "mmm_study"].groupby("strategy").size() == 2376).all() and (DP.strategy == "mmm_study").sum() == 480)
# ---- Table 4 (main text) ---------------------------------------------------------------------------------------------
T4 = pd.read_csv("tables/V5_S7f_table4.csv").set_index("strategy")
D4 = cells(MB["Table 5"]); S4 = sorted(set(i.split("|")[1] for i in D4.id if i.startswith("f:t4|")))   # S9: v6 numbering (the rule table is Table 5)
agg = DP[DP.strategy.isin(S4)].groupby("strategy").apply(lambda g: pd.Series({
    "pu": g.pu.mean(), "upow": g.upow.mean(), "pu_small": g.pu[g.small].mean(),
    "b_pu": g.mc_pu.mean(), "b_upow": g.mc_upow.mean(), "b_pu_small": g.mc_pu[g.small].mean()}), include_groups=False)
for q in ["pu", "upow", "pu_small"]:
    chk(f"Table 4 {q}: recomputed means equal V5_S7f_table4.csv", np.allclose(agg[q], T4.loc[agg.index, q], atol=1e-9), f"{len(agg)} strategies")
    w = agg[f"b_{q}"].idxmax()
    put(f"cap|t4|{q}", agg[f"b_{q}"].max(), "bound", "tables/V5_S7f_cp_S7d_dp_rates.csv; V5_S7f_cp_S7e_dp_rates.csv; V5_S7f_cp_S7b_rule_dp.csv; V5_S7f_cp_S7d2_flip_rules_dp.csv",
        f"max over the {len(agg)} strategies of Table 4 of the mean over design points of 100 sqrt(p(1-p)/nd) (at {w}; power with the H1 estimable designs)")
# ---- Table S26: usable share on the 1,896 points ------------------------------------------------------------------------
D = cells(SB["Table S26"]); s26 = sorted(set(f["strategy"] for f in D.flt if "strategy" in f))  # S8b: skip the f:s9 rows of the MD/KC average
T2 = pd.read_csv("tables/V5_S7f_rule_t2.csv"); T2 = T2[T2.df == "N - p"].set_index("strategy")
d = DP[DP.strategy.isin(s26) & (DP.dgp != "Copula-AR1")]
g = d.groupby("strategy").agg(pu=("pu", "mean"), bnd=("mc_pu", "mean"), npts=("pu", "size"))
chk("S26 usable: recomputed means equal V5_S7f_rule_t2.csv on 1,896 points", (g.npts == 1896).all() and np.allclose(g.pu, T2.loc[g.index, "pu"], atol=1e-9), list(g.index))
put("cap|S26|pu", g.bnd.max(), "bound", "tables/V5_S7f_cp_S7d_dp_rates.csv; V5_S7f_cp_S7b_rule_dp.csv", f"max over the {len(g)} strategies of Table S26 of the mean over the 1,896 points of 100 sqrt(p(1-p)/nd)")
# ---- Table S28: Copula-AR1 --------------------------------------------------------------------------------------------
D = cells(SB["Table S28"]); s28 = sorted(set(f["strategy"] for f in D.flt if "dgp" in f) | {"mmm_study"})  # S8b review N303: the printed row of the model without the check
S7e = pd.read_csv("tables/V5_S7f_cp_S7e_summary.csv"); S7e = S7e[S7e.dgp == "Copula-AR1"].set_index("strategy")
d = DP[(DP.dgp == "Copula-AR1") & DP.strategy.isin(s28)]
g = d.groupby("strategy").apply(lambda x: pd.Series({"t1": x.rate.mean(), "b_t1": x.mc[x.rate.notna()].mean(), "pu": x.pu.mean(), "b_pu": x.mc_pu.mean(),
                                                     "pw": x.upow.mean(), "b_pw": x.mc_upow.mean(), "npts": len(x)}), include_groups=False)
for q, c in [("t1", "mean_type1"), ("pu", "mean_usable"), ("pw", "mean_power")]:
    chk(f"S28 {c}: recomputed means equal V5_S7f_cp_S7e_summary.csv", np.allclose(g[q], S7e.loc[g.index, c], atol=1e-9), list(g.index))
put("cap|S28|type1", g.b_t1.max(), "bound", "tables/V5_S7f_cp_S7e_dp_rates.csv; V5_S7f_cp_S7d2_flip_rules_dp.csv", "max over the strategies of the mean over the 480 points of the per-point binomial MCSE of the rate")
put("cap|S28|usable", g.b_pu.max(), "bound", "tables/V5_S7f_cp_S7e_dp_rates.csv", "max over the strategies of the mean over the 480 points of 100 sqrt(p(1-p)/nd)")
put("cap|S28|power", g.b_pw.max(), "bound", "tables/V5_S7f_cp_S7e_dp_rates.csv", "max over the strategies of the mean over the 480 points of 100 sqrt(q(1-q)/nd_H1)")
BC = pd.read_csv("tables/V4_S7e_bias_cov.csv"); BS = pd.read_csv("tables/V4_S7e_bias_summary.csv").set_index("strategy")
sb = sorted(set(f["strategy"] for f, c in zip(D.flt, D.col) if c.startswith(("Bias", "Coverage"))))
bc = BC[(BC.label == "mmlb") & BC.strategy.isin(sb)].copy()
bc["mc_rb"] = bc.mcse_bias / bc.empse; bc["mc_ct"] = b(bc.cov_t, bc.n)
g = bc.groupby("strategy").apply(lambda x: pd.Series({"rb": x.bias_rel.mean(), "cn": x.cov_n.mean(), "ct": x.cov_t.mean(),
    "b_rb": x.mc_rb[x.bias_rel.notna()].mean(), "b_cn": x.mcse_cov_n[x.cov_n.notna()].mean(), "b_ct": x.mc_ct[x.cov_t.notna()].mean()}), include_groups=False)
for q, c in [("rb", "mean_bias_rel"), ("cn", "mean_cov_n"), ("ct", "mean_cov_t")]:
    chk(f"S28 {c}: recomputed means equal V4_S7e_bias_summary.csv", np.allclose(g[q], BS.loc[g.index, c], atol=1e-9), list(g.index))
put("cap|S28|biasse", g.b_rb.max(), "bound", "tables/V4_S7e_bias_cov.csv", "max over the marginal estimators of the mean over points of mcse_bias / empse")
put("cap|S28|cov", max(g.b_cn.max(), g.b_ct.max()), "bound", "tables/V4_S7e_bias_cov.csv", "max over estimators and references of the mean over points of the binomial MCSE of coverage")
# ---- Table S14: M6-M7 nested-loop strategies ----------------------------------------------------------------------------
NS = pd.read_csv("tables/V4_S6b_ns_rule_safety_dp.csv", keep_default_na=False, na_values=[""])
D = cells(SB["Table S14"]); d = D[D.col.isin(["Largest rate (%)", "Usable (%)"])]
mx, bp, rows = [], [], set()
for f, col, v in zip(d.flt, d.col, d.v):
    x = NS[(NS.strategy == f["strategy"]) & (NS.dgp == f["dgp"]) & (NS.scen == f["scen"])]
    assert len(x) == 96, (f, len(x)); rows.add((f["strategy"], f["dgp"], f["scen"]))
    if col.startswith("Largest"):
        a = x[x.assess.astype(str).str.upper() == "TRUE"]; top = a[np.isclose(a.rate, a.rate.max())]
        chk(f"S14 largest rate {f['strategy']} {f['dgp']} {f['scen']} equals the table", np.isclose(a.rate.max(), v, atol=1e-9))
        mx.append(float(b(top.rate, top.n).max()))
    else:
        chk(f"S14 usable {f['strategy']} {f['dgp']} {f['scen']} equals the table", np.isclose(x.pu.mean(), v, atol=1e-9))
        bp.append(float(b(x.pu, x.nd).mean()))
put("cap|S14|maxrate", max(mx), "binomial", "tables/V4_S6b_ns_rule_safety_dp.csv", f"max over the {len(mx)} rows of 100 sqrt(p(1-p)/n) at the design point with the largest rate")
put("cap|S14|pu", max(bp), "bound", "tables/V4_S6b_ns_rule_safety_dp.csv", f"max over the {len(bp)} rows of the mean over the 96 points of 100 sqrt(p(1-p)/nd)")
# ---- Table S15: clinical shares over the subsamples ---------------------------------------------------------------------
LT = pd.read_csv("tables/V5_S7f_lus_tab316.csv", dtype={"N": str}); NSUB = int(float(RI.loc["x:lus_nsub", "value"]))
D = cells(SB["Table S15"]); d = D[D.col != "N"]; mc = []
for f, fld, v in zip(d.flt, d.field, d.v):
    x = LT
    for k, val in f.items(): x = x[x[k].astype(str) == str(val)]
    assert len(x) == 1; mc.append(float(b(x[fld].iloc[0], NSUB)))
chk("S15: every cell found once", len(mc) == len(d), len(mc))
put("cap|S15|share", max(mc), "binomial", "tables/V5_S7f_lus_tab316.csv", f"max over the {len(mc)} cells of 100 sqrt(p(1-p)/{NSUB}), {NSUB} subsamples (x:lus_nsub)")
# ---- Table S27: shares of datasets -------------------------------------------------------------------------------------
EF = pd.read_csv("tables/V4_S7d2_explain_finite.csv", dtype=str)
D = cells(SB["Table S27"]); mf, mg = [], []
for f, fld in zip(D.flt, D.field):
    x = EF
    for k, val in f.items(): x = x[x[k] == str(val)]
    assert len(x) == 1
    if fld == "finite": mf.append(float(b(float(x.finite.iloc[0]), float(x.datasets.iloc[0]))))
    if fld == "gee_ar1_fg_usable": mg.append(float(b(float(x.gee_ar1_fg_usable.iloc[0]), float(x.datasets.iloc[0]))))
put("cap|S27|finite", max(mf), "binomial", "tables/V4_S7d2_explain_finite.csv", f"max over the {len(mf)} rows of 100 sqrt(p(1-p)/datasets); approximate, the datasets of a row share replicate seeds across N")
put("cap|S27|gee", max(mg), "binomial", "tables/V4_S7d2_explain_finite.csv", f"max over the {len(mg)} rows of 100 sqrt(p(1-p)/datasets); approximate, as above")
# ---- Table S32: highest rate per strategy -------------------------------------------------------------------------------
CA = pd.read_csv("tables/V5_S7g_calib.csv").set_index("strategy")
D = cells(SB["Table S32"]); s32 = sorted(set(re.sub(r"^f:gt\|cal\|([a-z0-9_]+)\|.*$", r"\1", i) for i in D.id))
ca = CA.loc[s32]
chk("S32: max_mc = 100 sqrt(p(1-p)/n) at the highest rate", np.allclose(ca.max_mc, b(ca.max_rate, ca.max_n), atol=1e-9), len(ca))
put("cap|S32|maxrate", ca.max_mc.max(), "binomial", "tables/V5_S7g_calib.csv", f"max of max_mc over the {len(ca)} strategies (binomial MCSE at the point of the highest rate, at {ca.max_mc.idxmax()})")
# ---- Table S33: usable share by usability limit -------------------------------------------------------------------------
FD = pd.read_csv("tables/V5_S7g_filter_dp.csv"); FT = pd.read_csv("tables/V5_S7g_filter_t4.csv").set_index(["strategy", "L"])
D = cells(SB["Table S33"]); d = D[D.col.str.startswith("Usable")]
sl = sorted(set((m.group(1), int(m.group(2))) for m in (re.match(r"^f:gt\|flt\|([a-z0-9_]+)\|(\d+)\|pu", i) for i in d.id)))
g = FD.assign(mc=b(FD.pu, FD.nd)).groupby(["strategy", "L"]).agg(pu=("pu", "mean"), bnd=("mc", "mean"), npts=("pu", "size")).loc[sl]
chk("S33: recomputed mean usable equals V5_S7g_filter_t4.csv", np.allclose(g.pu, FT.loc[g.index, "pu"], atol=1e-9), f"{len(g)} strategy x limit")
put("cap|S33|pu", g.bnd.max(), "bound", "tables/V5_S7g_filter_dp.csv", f"max over the {len(g)} strategy x limit cells of the mean over points of 100 sqrt(p(1-p)/nd)")
# ---- Table S34: Laplace against adaptive quadrature ---------------------------------------------------------------------
AQ = pd.read_csv("tables/V5_S7g_agq_dp.csv"); AT = pd.read_csv("tables/V5_S7g_agq_tab.csv")
MS = {"GB": "GenBinaryY", "NAT": "GLMM-native", "AR1": "GLMM-AR1", "COP": "Copula-AR1"}
D = cells(SB["Table S34"]); used = sorted(set((m.group(1), int(m.group(2))) for m in (re.match(r"^f:gq\|t_([A-Z0-9]+)_(\d+)_", i) for i in D.id) if m))
bt, bu = [], []
for mech, T in used:
    for nag, sfx in [(1, "lap"), (10, "agq")]:
        x = AQ[(AQ.dgp == MS[mech]) & (AQ.ntime == T) & (AQ.nAGQ == nag)]
        ref = AT[(AT.dgp == MS[mech]) & (AT.ntime == T)].iloc[0]
        ok_r = np.isclose(x.rate_t.mean(), ref[f"r{sfx}"], atol=1e-9); ok_p = np.isclose(x.pu.mean(), ref[f"p{sfx}"], atol=1e-9)
        chk(f"S34 {mech} T={T} {sfx}: mean rate and usable equal V5_S7g_agq_tab.csv", ok_r and ok_p and len(x) == ref.n, len(x))
        bt.append(float(x.mc_t[x.rate_t.notna()].mean())); bu.append(float(b(x.pu, x.nd).mean()))
put("cap|S34|type1", max(bt), "bound", "tables/V5_S7g_agq_dp.csv", f"max over the {len(bt)} cells (mechanism x T x method) of the mean over points of the binomial MCSE of the rate (mc_t)")
put("cap|S34|pu", max(bu), "bound", "tables/V5_S7g_agq_dp.csv", f"max over the {len(bu)} cells of the mean over points of 100 sqrt(p(1-p)/nd)")
# ---- Table S35: marginalized model with the check in M6-M7 -------------------------------------------------------------
MP = pd.read_csv("tables/V5_S7g_m67_cc_dp.csv"); M67 = pd.read_csv("tables/V5_S7g_m67_cc.csv")
MIN_OBS = int(float(RI.loc["x:min_obs_dp", "value"])); MP = MP[MP.n >= MIN_OBS].copy()
MP["mc_rb"] = 1 / np.sqrt(MP.n); MP["mc_cov"] = b(MP["cov"], MP.n); MP["mc_sh"] = b(100 * MP.share_B, MP.n)
MP["mc_ser"] = 100 * (1 + MP.ser) / np.sqrt(2 * (MP.n - 1))
r = {"rb": [], "rbmax": [], "cov": [], "covmin": [], "sh": [], "ser": []}
for _, m in M67.iterrows():
    x = MP[(MP.dgp == m.dgp) & (MP.scen == m.scen) & (MP.arm == m.arm)]
    ok = len(x) == m.n_points and np.isclose(x.rb.mean(), m.mean_rb) and np.isclose(x["cov"].mean(), m.mean_cov) and np.isclose(100 * x.share_B.mean(), m.mean_share_B)
    chk(f"S35 {m.dgp} {m.scen} {m.arm}: recomputed means equal V5_S7g_m67_cc.csv", ok, len(x))
    r["rb"].append(x.mc_rb.mean()); r["cov"].append(x.mc_cov.mean()); r["sh"].append(x.mc_sh.mean())
    r["rbmax"].append(float(x.mc_rb[np.isclose(x.rb.abs(), x.rb.abs().max())].max())); r["covmin"].append(float(x.mc_cov[np.isclose(x["cov"], x["cov"].min())].max()))
    r["ser"].append(x.mc_ser.median())
put("cap|S35|rb", max(r["rb"]), "approx", "tables/V5_S7g_m67_cc_dp.csv", "max over the 8 rows of the mean over points of 1 / sqrt(n) (bias / SE, empirical SE treated as known; bound for the mean)")
put("cap|S35|rbmax", max(r["rbmax"]), "approx", "tables/V5_S7g_m67_cc_dp.csv", "max over rows of 1 / sqrt(n) at the point with the largest |bias / SE|")
put("cap|S35|ser", max(r["ser"]), "approx", "tables/V5_S7g_m67_cc_dp.csv", "max over the 8 rows of the median over points of 100 (1 + ratio) / sqrt(2 (n - 1)) (typical per-point MCSE of the ratio, empirical-SE component only; not a bound for the median over points)")
put("cap|S35|cov", max(r["cov"]), "bound", "tables/V5_S7g_m67_cc_dp.csv", "max over rows of the mean over points of the binomial MCSE of coverage")
put("cap|S35|covmin", max(r["covmin"]), "binomial", "tables/V5_S7g_m67_cc_dp.csv", "max over rows of the binomial MCSE at the point with the smallest coverage")
put("cap|S35|sh", max(r["sh"]), "bound", "tables/V5_S7g_m67_cc_dp.csv", "max over rows of the mean over points of the binomial MCSE of the share of fits in which the check intervened")
# ---- Figure 2: reuse x:mcse_5pct; largest per-point MCSE among drawn points, for information ----------------------------
chk("Figure 2: registry id x:mcse_5pct exists", "x:mcse_5pct" in RI.index, RI.loc["x:mcse_5pct", "text"])
NLf = "tables/variant_cv/V4_S6_nestedloop.csv"
if os.path.exists(NLf):
    NL = pd.read_csv(NLf); NL = NL[NL.n_obs >= MIN_OBS]
    put("cap|fig2|mcse_max", NL.mcse_rej.max(), "binomial", NLf, f"largest per-point MCSE of a drawn rate (n_obs >= {MIN_OBS}), all methods (percentage points)")
K = pd.DataFrame(KEYS, columns=["key", "value", "mcse", "kind", "source", "how"])
assert K.key.is_unique
K.to_csv("tables/V5_S8b_sub_cap_keys.csv", index=False)
pd.DataFrame(CHK, columns=["check", "pass", "detail"]).to_csv("tables/V5_S8b_sub_cap_checks.csv", index=False)
with pd.option_context("display.max_colwidth", 60, "display.width", 200): print(K[["key", "value", "how"]].to_string())
nf = sum(not c[1] for c in CHK); print(f"keys {len(K)} | checks {len(CHK)} | failed {nf}")
assert nf == 0, "a check failed"
