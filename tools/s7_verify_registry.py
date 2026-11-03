#!/usr/bin/env python3
"""S7 gate G3: re-read every number of tables/V4_S7_numbers.csv from its source file, independently of the
R builder (stdlib only). Checks value (|diff| < 1e-9, or equal text), printed text (rounding rule of D-39)
and md5 of the source file. Exit 1 on any discordance. Run from ~/MMM: python3 -B tools/s7_verify_registry.py"""
import csv, json, hashlib, math, sys
REG = list(csv.DictReader(open("tables/V4_S7_numbers.csv", newline="", encoding="utf-8")))
cache = {}
def rows(f):
    if f not in cache: cache[f] = list(csv.DictReader(open(f, newline="", encoding="utf-8")))
    return cache[f]
def num(x):
    try: return float(x)
    except Exception: return None
def eq(col, v):
    if v == "": return col in ("", "NA", None)
    a, b = num(col), num(v)
    if b is not None: return a is not None and a == b
    return col == v
def sel(f, filt):
    out = []
    for r in rows(f):
        ok = True
        for k, v in filt.items():
            v = str(v)
            if k.endswith("__lt"):
                a = num(r[k[:-4]]); ok = a is not None and a < float(v)
            elif k.endswith("__gt"):
                a = num(r[k[:-4]]); ok = a is not None and a > float(v)
            elif k.endswith("__in"):
                a = num(r[k[:-4]]); ok = a is not None and a in [float(z) for z in v.split(",")]
            else:
                ok = eq(r[k], v)
            if not ok: break
        if ok: out.append(r)
    return out
def sig2(v, mcse):
    if v == 0: return "0"
    e = math.floor(math.log10(abs(v))); d = max(0, 1 - e)
    if mcse is not None and math.isfinite(mcse) and mcse > 0: d = min(d, max(0, -math.floor(math.log10(mcse))))
    s = "%.*f" % (d, v)
    if abs(float(s)) >= 10 ** (e + 1) and d > 0: s = "%.*f" % (d - 1, v)
    return s
def fmt(v, mcse, kind, f):
    if kind == "int": return "{:,}".format(int(round(v)))
    if kind == "plain": return "%.15g" % v
    if f == "tab": return "%.*f" % (0 if (mcse is not None and math.isfinite(mcse) and mcse >= 0.5) else 1, v)
    return sig2(v, mcse)
md5 = {}
bad = []
for r in REG:
    f = r["file"]
    if f not in md5: md5[f] = hashlib.md5(open(f, "rb").read()).hexdigest()
    if md5[f] != r["md5"]: bad.append((r["id"], "md5")); continue
    R = sel(f, json.loads(r["filter"])); op = r["op"]; mc = None
    if op == "count": val = len(R)
    elif op == "value":
        if len(R) != 1: bad.append((r["id"], f"{len(R)} rows")); continue
        val = R[0][r["field"]]
        if r["mcse_field"]: mc = num(R[0][r["mcse_field"]])
    else:
        xs = [float(x[r["field"]]) for x in R]; val = {"min": min, "max": max, "sum": sum}[op](xs)
    if r["kind"] == "text":
        if str(val) != r["value"] or str(val) != r["text"]: bad.append((r["id"], "text")); continue
        continue
    v = float(val) * float(r["scale"])
    if r["abs"] == "TRUE": v = abs(v)
    if mc is not None: mc = mc * float(r["scale"])
    if abs(v - float(r["value"])) > 1e-9 * max(1.0, abs(v)): bad.append((r["id"], f"value {v} vs {r['value']}")); continue
    t = fmt(v, mc, r["kind"], r["fmt"])
    if t != r["text"]: bad.append((r["id"], f"text {t} vs {r['text']}"))
print(f"G3: {len(REG)} numbers, {len(md5)} source files, {len(bad)} discordances")
for b in bad[:50]: print("  ", b)
sys.exit(1 if bad else 0)
