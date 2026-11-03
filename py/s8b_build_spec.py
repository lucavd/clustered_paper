"""s8b_build_spec.py -- S8b (copy of py/s7g_build_spec.py with the families f:s8, f:s9, f:pd, f:rux and the D-70 remaps of the cohort values of the marginalized model). s7g_build_spec.py -- S7g (copy of py/s7f_build_spec.py with the rules for f:g7, f:gq, f:gt). S7f: manuscript/v5/numbers_spec.csv from the tokens of the v5 manuscript and supplement.
Tokens already in the v4 specification keep their selector, with the source moved to the tables of the variant with
Kauermann-Carroll from fewclust (D-60) where the v4 source was the variant with the convergence check; new tokens get
their selector from the rules below. Every row: id, file, filter (JSON), field, op, scale, abs, fmt, kind, mcse_field.
Usage (from ~/MMM): python3 -B py/s7f_build_spec.py"""
import csv, json, re, sys, os
V4 = {r["id"]: r for r in csv.DictReader(open("manuscript/v4/numbers_spec.csv", encoding="utf-8"))}
toks = []
for f in ("manuscript/v5/MMM_v5.qmd", "manuscript/v5/MMM_v5_supplement.qmd"):
    toks += re.findall(r'`r n\("([^"]+)"\)`', open(f, encoding="utf-8").read())
for f in ("manuscript/v5/MMM_v5.qmd",):   # registry values named by reg: tags of the claim check (as in the v4 registry)
    for tag in re.findall(r"<!-- claim: (.*?) -->", open(f, encoding="utf-8").read()):
        toks += ["x:" + k for k in re.findall(r"reg:(\w+)", tag) if "x:" + k in V4]
ids = sorted(set(toks))
REMAP = {"tables/variant_cc/V4_S6_values.csv": "tables/variant_cv/V4_S6_values.csv",
         "tables/variant_cc/V4_S5_claims.csv": "tables/variant_cv/V4_S5_claims.csv",
         "tables/variant_cc/V4_S6_table_cells.csv": "tables/variant_cv/V4_S6_table_cells.csv",
         "tables/variant_cc/V4_S6_rule_safety.csv": "tables/variant_cv/V4_S6_rule_safety.csv",
         "tables/V4_S6_m1c_claims.csv": "tables/variant_cv/V4_S6_m1c_claims.csv",
         "rev2/T5_tab316.csv": "tables/V5_S7f_lus_tab316.csv",
         "tables/V4_S7_u1.csv": "tables/V5_S7f_u1.csv",
         # exact-interval flags (PI, 2026-10-27; R/s7f/s7f_cpflags.R)
         "tables/V4_S6b_ns_rule_safety.csv": "tables/V5_S7f_cp_S6b_ns_rule_safety.csv",
         "tables/V4_S6cc_cells.csv": "tables/V5_S7f_cp_S6cc_cells.csv",
         "tables/V4_S7d_rule.csv": "tables/V5_S7f_cp_S7d_rule.csv",
         "tables/V4_S7e_summary.csv": "tables/V5_S7f_cp_S7e_summary.csv"}
REMAP2 = {"tables/variant_cv/V4_S6_rule_safety.csv": "tables/V5_S7f_cp_kcfc_rule_safety.csv"}
S7G_FILES = {"f:g7": "tables/V5_S7g_extra.csv", "f:gq": "tables/V5_S7g_agq_keys.csv", "f:gt": "tables/V5_S7g_tab_keys.csv"}
S7G_KINDS = {f: {r["key"]: r["kind"] for r in csv.DictReader(open(f, encoding="utf-8"))} for f in S7G_FILES.values()}
S8B_FILES = {"f:s8": "tables/V5_S8b_keys.csv", "f:s9": "tables/V5_S8b_sub_keys.csv", "f:n9": "tables/V6_S9_keys.csv", "f:n11": "tables/V8_S11_keys.csv"}   # S9: f:n9 (R/s9/s9_keys.R); S11: f:n11 (R/s11/s11_keys.R)
S8B_KINDS = {f: ({r["key"]: r["kind"] for r in csv.DictReader(open(f, encoding="utf-8"))} if os.path.exists(f) else {}) for f in S8B_FILES.values()}
# D-70 (S8b): cohort values of the marginalized model from the latent-variable refit with the model-based SE
ID_REMAP = {"x:lus_est_mmm": ("tables/V5_S7f_lus_ref.csv", {"spec": "pieno", "method": "MMM"}, "estimate"),
            "x:lus_se_mmm": ("tables/V5_S7f_lus_ref.csv", {"spec": "pieno", "method": "MMM"}, "se"),
            "x:lus_bias_mmm_10": ("tables/V5_S7f_lus_bias.csv", {"spec": "pieno", "method": "MMM", "N": "10"}, "bias"),
            "x:lus_lateral": ("tables/V5_S8b_keys.csv", {"key": "lus_lateral_lv"}, "value")}
def row(i, file, filt, field, op="value", kind="sig2", fmt="text", scale="1", ab="FALSE", mcse=""):
    return {"id": i, "file": file, "filter": json.dumps(filt), "field": field, "op": op, "scale": scale, "abs": ab, "fmt": fmt, "kind": kind, "mcse_field": mcse}
INTX = {"ndp_all5", "ndp_cop", "ndp_three", "flip_B", "flip_enum_maxn", "guard_div", "cop_exfg_T", "cop_exfg_Nmin", "cop_exfg_Nmax", "upow_small_ceiling", "cop_exfg_nT", "glmer_hi_native", "glmer_ndp_native"}
PLAINX = {"cop_rho", "flip_coefmax", "guard_u0", "fg_bound"}
FLSET = {"all": "all 2,376 points", "small": "small designs (T 2-3, N <= 25)"}
CIQ = {"PG1": "PG1 coverage of the test at the truth, H1", "PG5": "PG5 finite interval share, H1",
       "PG6": "PG6 finite interval share, H1, small points", "PG8": "PG8 coverage GEE AR(1) FG t, H1"}
def new(i):
    base, tab = (i[:-4], True) if i.endswith("@tab") else (i, False)
    ab = "FALSE"
    if base.endswith("@abs"): base, ab = base[:-4], "TRUE"
    fmt = "tab" if tab else "text"
    p = base.split("|"); h = p[0]
    if h == "f:t4":
        k = "int" if re.match(r"^(hi_|lo_|n_)", p[2]) else "sig2"
        return row(i, "tables/V5_S7f_table4.csv", {"strategy": p[1]}, p[2], kind=k, fmt=fmt)
    if h == "f:x":
        k = "int" if p[1] in INTX else "plain" if p[1] in PLAINX else "text" if p[1] == "git_commit" else "sig2"
        return row(i, "tables/V5_S7f_extra.csv", {"key": p[1]}, "value", kind=k, fmt="text" if k == "text" else fmt)
    if h == "f:ge4":
        if p[1] == "nmin": return row(i, "tables/V4_S7d_GE4.csv", {}, "n", op="min", kind="int")
        return row(i, "tables/V5_S7f_extra.csv", {"key": "ge4_max"}, "value", kind="plain")
    if h == "f:cv": return row(i, "tables/V5_S7f_claims.csv", {"name": p[1]}, p[2], ab=ab)
    if h == "f:t2":
        k = "int" if p[3].startswith("n_") else "sig2"
        return row(i, "tables/V5_S7f_rule_t2.csv", {"df": p[2], "strategy": p[1]}, p[3], kind=k, fmt=fmt)
    if h in ("f:fl", "f:fl7"):
        return row(i, "tables/V4_S7d2_flip_rules.csv", {"label": p[2], "version": "with D-57 rules" if h == "f:fl" else "as in S7d", "set": FLSET[p[1]]}, p[3])
    if h == "f:ci": return row(i, "tables/V4_S7d2_boot_means.csv", {"q": CIQ[p[1]]}, "est")
    if h == "f:bs": return row(i, "tables/V4_S7e_bias_summary.csv", {"strategy": p[1]}, p[2], fmt=fmt)
    if h == "f:s7e":
        k = "int" if p[2].startswith("clearly_") else "sig2"
        return row(i, "tables/V5_S7f_cp_S7e_summary.csv", {"strategy": p[1], "dgp": "Copula-AR1"}, p[2], kind=k, fmt=fmt)
    if h == "f:sr":
        if p[1] == "neo": return row(i, "tables/V5_S7f_seratio_fpc.csv", {"application": "neonatal cohort", "spec": "pieno", "method": p[2]}, "ratio_fpc", op=p[3])
        return row(i, "tables/V5_S7f_seratio_fpc.csv", {"application": "respiratory", "spec": "pieno", "method": p[2], "N": p[3]}, "ratio_fpc", fmt=fmt)
    if h in ("f:rr", "f:rrr"):
        return row(i, "tables/V5_S7f_resp_ref.csv", {"method": p[1], "spec": "pieno" if h == "f:rr" else "ridotto"}, p[2])
    if h == "f:rs": return row(i, "tables/V5_S7f_seratio_fpc.csv", {"application": "respiratory", "spec": "pieno", "method": "MMM", "N": "10"}, "P", kind="int")
    if h == "f:ru":
        if p[1] in ("Nmin", "Nmax"): return row(i, "tables/V5_S7f_resp_usable.csv", {"spec": "pieno", "method": "MMM"}, "N", op=p[1][1:], kind="int")
        return row(i, "tables/V5_S7f_resp_usable.csv", {"spec": "pieno", "method": p[1], "N": p[2]}, "pct_usable_500", fmt=fmt, mcse="mcse_usable_500")
    if h == "f:rN": return row(i, "tables/V5_S7f_resp_usable.csv", {"spec": "pieno", "method": "MMM", "N": p[1]}, "N", kind="int")
    if h == "f:rb": return row(i, "tables/V5_S7f_resp_bias.csv", {"spec": "pieno", "method": p[1], "N": p[2]}, "bias", fmt=fmt, mcse="mcse_bias")
    if h == "f:ef":
        k = "int" if p[3] == "datasets" else "sig2"
        return row(i, "tables/V4_S7d2_explain_finite.csv", {"ntime": p[1], "Ng": p[2]}, p[3], kind=k, fmt=fmt)
    if h == "f:ec":
        if p[2] == "n": return row(i, "tables/V4_S7d_exch_coincidence.csv", {"corr": p[1]}, "n", kind="int")
        return row(i, "tables/V4_S7d_exch_coincidence.csv", {"corr": p[1]}, "share_within_1pct", scale="100")
    if h == "f:r7d":
        k = "int" if p[2] == "n_hi" else "sig2"
        return row(i, "tables/V5_S7f_cp_S7d_rule.csv", {"strategy": p[1]}, p[2], kind=k)
    if h == "f:sw": return row(i, "tables/V5_S7f_software.csv", {"package": p[1]}, "version", kind="text")
    if h == "f:ag":
        if p[1] == "min": return row(i, "tables/V5_S7f_kc_clinical_agree.csv", {}, "share_within_1pct", op="min")
        return row(i, "tables/V5_S7f_kc_clinical_agree.csv", {"app": p[1], "corr": p[2]}, "share_within_1pct")
    if h == "v:use20":
        return row(i, "tables/variant_cv/V4_S6_values.csv", {"variant": "main", "name": base[2:]}, "est", fmt=fmt, mcse="se")
    if h == "f:hh":
        k = "int" if p[4] in ("n_lo", "n_ok", "n_pts") else "sig2"
        ms = {"mad": "mad_se", "pu": "pu_se", "upow": "upow_se"}.get(p[4], "")
        return row(i, "tables/V5_S7f_h2h_grid.csv", {"strategy": p[1], "Nb": p[2], "Tb": p[3]}, p[4], kind=k, fmt=fmt, mcse=ms)
    if h == "f:hd": return row(i, "tables/V5_S7f_h2h_diff.csv", {"measure": p[1]}, p[2], mcse="se" if (p[2] == "est" and p[1] != "max_mcse_ef") else "")
    if h in ("f:hb", "f:hbs"):   # block-bootstrap direct comparison (s7f_h2h_boot.R): f:hb|<region>|<quantity>; region may contain "|"
        q, region = p[-1], "|".join(p[1:-1])
        if h == "f:hbs": return row(i, "tables/V5_S7f_h2h_boot.csv", {"region": region, "q": q}, "mcse", fmt=fmt)
        return row(i, "tables/V5_S7f_h2h_boot.csv", {"region": region, "q": q}, "est", fmt=fmt, mcse="mcse")
    if h == "f:hbc": return row(i, "tables/V5_S7f_h2h_boot.csv", {"region": "|".join(p[1:]), "q": "mad_gee_minus_mmm"}, "class", kind="text")
    if h == "f:dp": return row(i, "tables/V4_S7e_dependence_profiles.csv", {"mech": "copula floor=0.00 rho=0.6", "lag": p[1]}, "cor")
    if h in ("f:g7", "f:gq", "f:gt"):   # S7g registry sources (R/s7g/s7g_extras.R, s7g_agq_keys.R, s7g_tabkeys.R)
        FILE = S7G_FILES[h]; key = "|".join(p[1:]); kd = S7G_KINDS[FILE][key]
        k = "sig2" if kd == "tab1" else kd
        return row(i, FILE, {"key": key}, "value", kind=k, fmt="text" if k == "text" else fmt, ab=ab, mcse="mcse")
    if h in ("f:s8", "f:s9", "f:n9", "f:n11"):   # S8b registry sources (R/s8b/s8b_keys.R, R/s8b/sub/s8b_sub_keys.R)
        FILE = S8B_FILES[h]; key = "|".join(p[1:]); kd = S8B_KINDS[FILE][key]
        k = "sig2" if kd == "tab1" else kd
        return row(i, FILE, {"key": key}, "value", kind=k, fmt="text" if k == "text" else fmt, ab=ab, mcse="mcse")
    if h == "f:pd": return row(i, "tables/V4_dgm_parameters.csv", {"dgp": p[1], "model": p[2], "arm": "H1", "parameter": p[3]}, "value", kind="text", fmt="text")
    if h == "f:rux": return row(i, "tables/V5_S7f_resp_usable.csv", {"spec": "pieno", "method": p[1], "N": p[2]}, p[3], fmt=fmt)
    if i == "x:lus_se_gee": return row(i, "rev2/T2_fullcohort.csv", {"tag": "A:pieno", "method": "GEE"}, "se")
    raise SystemExit("no rule for token " + i)
OUT, n_old, n_new = [], 0, 0
for i in ids:
    if i in V4:
        r = dict(V4[i]); r["file"] = REMAP.get(r["file"], r["file"]); r["file"] = REMAP2.get(r["file"], r["file"])
        if i.startswith("x:ts_T3|") and json.loads(r["filter"]).get("i") != "3":
            r["file"] = "tables/variant_cv/V4_S6_table_cells.csv"   # Table S3: only the row of the marginalized model differs from Table 2
        if i in ID_REMAP: f_, flt_, fld_ = ID_REMAP[i]; r["file"], r["filter"], r["field"], r["mcse_field"] = f_, json.dumps(flt_), fld_, ""
        OUT.append(r); n_old += 1
    else:
        OUT.append(new(i)); n_new += 1
cols = ["id", "file", "filter", "field", "op", "scale", "abs", "fmt", "kind", "mcse_field"]
w = csv.DictWriter(open("manuscript/v5/numbers_spec.csv", "w", newline="", encoding="utf-8"), fieldnames=cols); w.writeheader()
for r in OUT: w.writerow({c: r[c] for c in cols})
print(f"spec: {len(OUT)} tokens ({n_old} from v4, {n_new} new) -> manuscript/v5/numbers_spec.csv")
