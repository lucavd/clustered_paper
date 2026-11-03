
# s10_ids.py -- S10: inventory of the registry ids used in the manuscript source (main and supplement), with line,
# section, table context and assignment to the S10 counter-check tracks. Read-only on the source.
import re, csv, random
out = []
for f, part in [("manuscript/v5/MMM_v5.qmd", "main"), ("manuscript/v5/MMM_v5_supplement.qmd", "supplement")]:
    sec = ""; tab = ""; lines = open(f).read().splitlines()
    for ln, l in enumerate(lines, 1):
        if l.startswith("#"): sec = l.lstrip("#").strip(); tab = ""
        m = re.match(r"^\s*:?\s*\**(Table S?\d+)[\.:]", l)
        if m: tab = m.group(1)
        in_table = l.lstrip().startswith("|")
        for k, i in enumerate(re.findall(r'n\("([^"]+)"\)', l)):
            out.append(dict(part=part, line=ln, section=sec, table_context=tab if in_table else "", in_table=in_table, id=i))
# assignment: main Tables 3-5 cells -> Tabelle; other main -> Testo; supplement -> Supplemento (sample)
seen = {}
for r in out:
    if r["part"] == "main":
        r["track"] = "Tabelle" if r["in_table"] and r["table_context"] in ("Table 3", "Table 4", "Table 5") else "Testo"
    else:
        r["track"] = "Supplemento"
# supplement sample: 10% of distinct ids per table (text counts as one stratum per section), min 5, seed 20261102
random.seed(20261102)
strata = {}
for r in out:
    if r["part"] == "supplement":
        s = r["table_context"] or ("text: " + r["section"])
        strata.setdefault(s, set()).add(r["id"])
pick = set()
for s in sorted(strata):
    ids = sorted(strata[s]); k = min(len(ids), max(5, -(-len(ids) // 10)))
    pick |= {(s, i) for i in random.sample(ids, k)}
for r in out:
    if r["part"] == "supplement":
        s = r["table_context"] or ("text: " + r["section"])
        r["stratum"] = s; r["sampled"] = (s, r["id"]) in pick
    else:
        r["stratum"] = r["table_context"] or ("text: " + r["section"]); r["sampled"] = True
with open("results_v4/S10/ids_inventory.csv", "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=["part", "line", "section", "table_context", "in_table", "id", "track", "stratum", "sampled"]); w.writeheader(); w.writerows(out)
from collections import Counter
print("rows", len(out))
print(Counter((r["part"], r["track"]) for r in out))
print("main distinct", len({r["id"] for r in out if r["part"] == "main"}),
      "Tabelle distinct", len({r["id"] for r in out if r["track"] == "Tabelle"}),
      "Testo distinct", len({r["id"] for r in out if r["track"] == "Testo"}))
print("supp strata", len(strata), "supp distinct", len({r["id"] for r in out if r["part"] == "supplement"}), "sampled distinct", len({i for _, i in pick}))
print(sorted(Counter(r["table_context"] for r in out if r["part"] == "main" and r["in_table"]).items()))
