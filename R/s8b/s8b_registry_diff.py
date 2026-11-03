
import csv, json, collections
old = {r["id"]: r for r in csv.DictReader(open("/mnt/micron/MMM_backup_S8b/V5_S7f_numbers.csv", newline="", encoding="utf-8"))}
new = {r["id"]: r for r in csv.DictReader(open("tables/V5_S7f_numbers.csv", newline="", encoding="utf-8"))}
ch = [(i, old[i]["text"], new[i]["text"], new[i]["file"]) for i in old if i in new and old[i]["text"] != new[i]["text"]]
fam = collections.Counter(i.split("|")[0] for i,_,_,_ in ch)
print("ids old", len(old), "new", len(new), "common", len(set(old)&set(new)), "removed", len(set(old)-set(new)), "added", len(set(new)-set(old)))
print("changed printed text:", len(ch), dict(fam))
with open("reports/morris/S8b_logs/registry_diff_S8b.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["id", "old_text", "new_text", "file"]); w.writerows(ch)
# GB3a: adjusted power ids vs reviewer recomputation
A = list(csv.DictReader(open("reports/morris/S8_logs/adjpower_fix_S8.csv")))
mx = 0; nmatch = 0; dig = 0; miss = []
for a in A:
    r = new.get(a["id"])
    if r is None: miss.append(a["id"]); continue
    d = abs(float(r["value"]) - float(a["corrected_value"])); mx = max(mx, d); nmatch += 1
    if old.get(a["id"]) and old[a["id"]]["text"] != r["text"]: dig += 1
print("GB3a: adjpower ids", len(A), "found", nmatch, "missing", miss[:5], "max |diff|", mx, "printed changes", dig)
# list of changed ids outside expected families
for x in ch: print(x)
