#!/usr/bin/env python3
"""S7 gate G4: claim tags in the manuscript source. Each <!-- claim: REF as=W; ... --> refers to a row of the claims
files (cl: main grid, variant with convergence check; cm: scenario M1c; suffix :bobyqa uses the bobyqa columns) or
to a registry value (reg:KEY). The wording W is checked against the bootstrap interval (estimate +/- 2 MCSE) and the
relevance margin: larger = interval beyond the margin, similar = interval within the margin, about1 = interval
excludes 0 and straddles the margin, notdist = interval includes 0, descr = descriptive (class confermata,
descrittiva or risolta). Sentences with comparative words that no tag covers are listed. Usage: s7_check_claims.py file.qmd"""
import re, sys, csv
CL = {r["name"]: r for r in csv.DictReader(open("tables/variant_kcfc/V4_S5_claims.csv"))}
CM = {r["name"]: r for r in csv.DictReader(open("tables/variant_kcfc/V4_S6_m1c_claims.csv"))}
CV = {r["name"]: r for r in csv.DictReader(open("tables/V5_S7f_claims.csv"))}
EXK = {r["key"]: r for r in csv.DictReader(open("tables/V5_S7f_extra.csv"))}
REG = {r["id"]: r for r in csv.DictReader(open("tables/V5_S7f_numbers.csv"))}
def f(x):
    try: return float(x)
    except Exception: return float("nan")
def judge(row, w, bob):
    lo, hi = (f(row["lo_bobyqa"]), f(row["hi_bobyqa"])) if bob else (f(row["lo"]), f(row["hi"]))
    mg = f(row["margin"]); cls = row["class_bobyqa"] if bob else row["class"]
    excl0 = lo * hi > 0; a, b = sorted([abs(lo), abs(hi)])
    if w == "larger": ok = excl0 and (mg != mg or a > mg)
    elif w == "similar": ok = mg == mg and -mg <= lo and hi <= mg
    elif w == "about1": ok = excl0 and mg == mg and a <= mg < b
    elif w == "notdist": ok = lo <= 0 <= hi
    elif w == "descr": ok = cls in ("confermata", "descrittiva", "risolta")
    else: ok = False
    return ok, f"[{lo:.3g}, {hi:.3g}] margin {mg} class {cls}"
src = open(sys.argv[1], encoding="utf-8").read()
src = src[:src.find("## References")] + (src[src.find("## Tables and figures"):] if "## Tables and figures" in src else "") if "## References" in src else src
rows = []; bad = 0
for m in re.finditer(r"<!-- claim: (.*?) -->", src):
    for part in m.group(1).split(";"):
        part = part.strip()
        mm = re.match(r"(cl|cm|cv):([^ ]+?)(:bobyqa)? as=(\w+)$", part)
        if mm:
            tab = {"cl": CL, "cm": CM, "cv": CV}[mm.group(1)]; name = mm.group(2)
            if name not in tab: rows.append((part, "MISSING", "")); bad += 1; continue
            ok, why = judge(tab[name], mm.group(4), bool(mm.group(3))); rows.append((part, "ok" if ok else "VIOLATION", why)); bad += (not ok); continue
        mr = re.match(r"reg:(\w+) as=descr$", part)
        if mr:
            ok = ("x:" + mr.group(1)) in REG or ("f:x|" + mr.group(1)) in REG or mr.group(1) in EXK or mr.group(1) in ("t4", "t2", "flip", "ci", "hh", "hb"); rows.append((part, "ok" if ok else "MISSING", REG.get("x:" + mr.group(1), REG.get("f:x|" + mr.group(1), {})).get("text", ""))); bad += (not ok); continue
        me = re.match(r"reg:(\w+) equals reg:(\w+)$", part)
        if me:
            a, b = REG.get("x:" + me.group(1)), REG.get("x:" + me.group(2)); ok = a is not None and b is not None and a["value"] == b["value"]
            rows.append((part, "ok" if ok else "VIOLATION", f"{a and a['value']} vs {b and b['value']}")); bad += (not ok); continue
        rows.append((part, "UNPARSED", "")); bad += 1
# uncovered comparative sentences (abstract, results, discussion, conclusions)
CUE = re.compile(r"\b(than|smallest|largest|highest|lowest|fewer|more often|less often|best|worst-case departure was larger|only strategy|lead|similar|furthest|closest|outperform|exceed)\b", re.I)
unc = []
txt = re.sub(r"`r n\(\"[^\"]*\"\)`", "§", src)
for para in re.split(r"\n\s*\n", txt):
    if para.lstrip().startswith(("|", "#", "**Table", "**Figure", "```", "!", "---")): continue
    pieces = re.split(r"(<!-- (?:claim|noclaim):.*?-->)", para)
    for k in range(0, len(pieces), 2):
        seg = pieces[k]; covered = k + 1 < len(pieces)
        if covered: continue
        for s in re.split(r"(?<=[.;])\s+(?=[A-Z])", seg):
            if CUE.search(s): unc.append(s.strip()[:220])
w = csv.writer(open("reports/morris/S7g_logs/S7g_GF4_claims.csv", "w", newline="")); w.writerow(["tag", "outcome", "interval"]); w.writerows(rows)
print(f"GF4: {len(rows)} tags, {bad} not ok; {len(unc)} comparative sentences without a tag")
for r in rows:
    if r[1] != "ok": print("  ", r)
for s in unc: print("  UNTAGGED:", s)
sys.exit(1 if (bad or unc) else 0)
