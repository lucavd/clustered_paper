
# s10_inbreve_sources.py -- S10: reads back from the files every number of the In breve summary and writes
# tables/V7_S10_inbreve_sources.csv (sentence, file, row, value). Read-only on the inputs.
import csv
T = "tables/"
def rows(f): return list(csv.DictReader(open(f)))
out = []
def add(sentence, f, row, value): out.append(dict(sentence=sentence, file=f, row=row, value=value))
g = rows(T + "V7_S10_degen_gates.csv"); add("4.752.000 dataset rigenerati", T + "V7_S10_degen_gates.csv", "GD2 datasets", g[0]["detail"])
inv = rows("results_v4/S10/ids_inventory.csv")
add("524 numeri del testo principale", "results_v4/S10/ids_inventory.csv", "part == main, distinct id", len({r["id"] for r in inv if r["part"] == "main"}))
add("367 numeri del supplemento", "results_v4/S10/ids_inventory.csv", "supplement, sampled", len({r["id"] for r in inv if r["part"] == "supplement" and r["sampled"] == "True"}))
a = rows(T + "V7_S10_assolute_check.csv"); from collections import Counter; ca = Counter(r["verdict"] for r in a)
add("151 frasi assolute", T + "V7_S10_assolute_check.csv", "verdict != NOT_A_RESULT_CLAIM", sum(v for k, v in ca.items() if k != "NOT_A_RESULT_CLAIM"))
add("3 frasi false (rilettura)", T + "V7_S10_assolute_check.csv", "verdict == FALSE", ca["FALSE"])
tb = rows(T + "V7_S10_tabelle_check.csv"); ct = Counter(r["status"] for r in tb if r["id"])
add("234 celle, 64 con scelte non scritte", T + "V7_S10_tabelle_check.csv", "rows with id, by status", dict(ct))
tx = rows(T + "V7_S10_testo_check.csv"); cx = Counter(r["status"] for r in tx)
add("tranne 4 (2 diversi, 2 non ricalcolabili); 1 frase falsa", T + "V7_S10_testo_check.csv", "status counts", dict(cx))
cl = [r for r in rows(T + "V7_S10_degen_classes.csv")]
al = [r for r in cl if r["scope"] == "all"][0]; add("quasi 8 su 100 separati", T + "V7_S10_degen_classes.csv", "scope all, pct_S", al["pct_S"])
for r in cl:
    if r["scope"] == "region": add("separati per regione", T + "V7_S10_degen_classes.csv", "region " + r["region"] + ", pct_S", r["pct_S"])
u = rows(T + "V7_S10_degen_usable.csv")
for st in ("mmm_cc", "gee_ar1_fg"):
    rr = [r for r in u if r["family"] == "T5" and r["strategy"] == st and r["class"] == "Si" and "all" not in (r["mechanism"], r["region"], r["arm"])]
    add("risultati usabili su dataset separati, " + st, T + "V7_S10_degen_usable.csv", "T5 %s Si, sum over %d mechanism x region x arm rows" % (st, len(rr)),
        "%d of %d" % (sum(int(r["usable_class"]) for r in rr), sum(int(r["usable_total"]) for r in rr)))
v = rows(T + "V7_S10_degen_table5_versions.csv")
for r in v:
    if r["strategy"] in ("gee_ar1_fg", "mmm_cc", "flip") and r["version"] in ("base", "a", "b"):
        add("regola per versione", T + "V7_S10_degen_table5_versions.csv", r["version"] + " " + r["strategy"], "hi_all %s pu %.2f pu_small %.2f chosen %s" % (r["hi_all"], float(r["pu"]), float(r["pu_small"]), r["chosen"]))
s = rows(T + "V7_S10_degen_sensitivity.csv")
t4 = [r for r in s if r["version"] == "a" and r["table"] == "Table 4"]
m = min(t4, key=lambda r: float(r["diff"])); add("fino a 5 risultati usabili su 100 in meno (Tabella 4)", T + "V7_S10_degen_sensitivity.csv", "version a Table 4 min diff: " + m["row"] + " " + m["column"], m["diff"])
cv = rows(T + "V7_S10_testo_cc_variants.csv")
for r in cv:
    if r["variant"] in ("ii", "iii") and r["quantity"] in ("hi_all", "pu", "pu_small") and "Table 5" in r["table"]:
        add("varianti del controllo", T + "V7_S10_testo_cc_variants.csv", r["variant"] + " " + r["quantity"], r["value"])
to = rows(T + "V7_S10_supplemento_geesmv_timeouts.csv"); mx = max(to, key=lambda r: float(r["pct_to"]))
add("time-out di geesmv fino al 20%", T + "V7_S10_supplemento_geesmv_timeouts.csv", "max pct_to: %s %s N%s T%s %s" % (mx["dgp"], mx["scen"], mx["N"], mx["ntime"], mx["label"]), mx["pct_to"])
c = rows(T + "V7_S10_degen_clinical.csv")
add("365 sottocampioni su 5.000 separati (studio respiratorio)", T + "V7_S10_degen_clinical.csv", "app respiratory, method (draws), sum Si / sum draws",
    "%d / %d" % (sum(int(r["Si"]) for r in c if r["method"] == "(draws)" and r["app"].startswith("resp")), sum(int(r["draws"]) for r in c if r["method"] == "(draws)" and r["app"].startswith("resp"))))
with open(T + "V7_S10_inbreve_sources.csv", "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=["sentence", "file", "row", "value"]); w.writeheader(); w.writerows(out)
for r in out: print(r["sentence"], "|", r["row"], "|", r["value"])
