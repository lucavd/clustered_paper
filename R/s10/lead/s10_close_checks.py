
# s10_close_checks.py -- S10 closing checks: GD4 coverage of the 524 main-text ids and the 367 sampled supplement ids,
# supplement YAML title vs main title, source of the fewclust coherence ids. Read-only.
import csv, re, collections
inv = list(csv.DictReader(open("results_v4/S10/ids_inventory.csv")))
main_ids = {r["id"] for r in inv if r["part"] == "main"}
sup_ids = {r["id"] for r in inv if r["part"] == "supplement" and r["sampled"] == "True"}
def ids_of(path):
    rows = list(csv.DictReader(open(path))); cols = rows[0].keys() if rows else []
    idc = [c for c in cols if c.lower() in ("id", "reg_id", "registry_id", "ids")]
    out = set()
    for r in rows:
        for c in idc:
            for x in re.split(r"[;, ]+", r[c] or ""):
                x = x.strip().replace("@tab", "") if False else x.strip()
                if x: out.add(x)
    return out, idc, len(rows)
seen = set()
for f in ["tables/V7_S10_tabelle_check.csv", "tables/V7_S10_testo_check.csv"]:
    s, idc, n = ids_of(f); print(f, "rows", n, "id columns", idc, "distinct ids", len(s)); seen |= s
miss = sorted(main_ids - seen)
print("main ids", len(main_ids), "covered", len(main_ids & seen), "missing", len(miss), miss[:15])
s, idc, n = ids_of("tables/V7_S10_supplemento_check.csv"); print("supp check rows", n, "id cols", idc)
ms = sorted(sup_ids - s); print("supp sampled", len(sup_ids), "covered", len(sup_ids & s), "missing", len(ms), ms[:10])
a = list(csv.DictReader(open("tables/V7_S10_assolute_check.csv")))
print("assolute rows", len(a), collections.Counter((r.get("status") or r.get("esito") or "?") for r in a).most_common(8))
for f in ["manuscript/v5/MMM_v5.qmd", "manuscript/v5/MMM_v5_supplement.qmd"]:
    t = open(f).read(); m = re.search(r"^title:\s*(.+)$", t, re.M); print(f, "YAML title:", m.group(1)[:220] if m else None)
spec = {r["id"]: r for r in csv.DictReader(open("manuscript/v5/numbers_spec.csv"))}
for i in ["f:n9|fc|n_cmp", "f:n9|fc|n_coh"]:
    print(i, {k: spec.get(i, {}).get(k) for k in ("file", "filter", "field")})
