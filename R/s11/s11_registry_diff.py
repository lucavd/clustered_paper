import csv, subprocess, io, re, json
old = {r["id"]: r["text"] for r in csv.DictReader(io.StringIO(subprocess.run(["git", "show", "7a68326:tables/V5_S7f_numbers.csv"], capture_output=True, text=True, cwd="/home/user/MMM").stdout))}
new = {r["id"]: r["text"] for r in csv.DictReader(open("/home/user/MMM/tables/V5_S7f_numbers.csv"))}
ch = {i: (old[i], new[i]) for i in old if i in new and old[i] != new[i]}
print("ids before", len(old), "after", len(new), "changed", len(ch), "removed", len(set(old) - set(new)), "added", len(set(new) - set(old)))
out = []
for f in ("MMM_v5.qmd", "MMM_v5_supplement.qmd"):
    src = open("/home/user/MMM/manuscript/v5/" + f).read()
    for para in re.split(r"\n\s*\n", src):
        ids = re.findall(r'n\("([^"]+)"\)', para)
        hit = [i for i in ids if i in ch]
        if hit and not para.lstrip().startswith("|"):
            out.append({"file": f, "ids": {i: ch[i] for i in hit}, "text": re.sub(r"<!--.*?-->", "", para)[:3000]})
json.dump({"changed": ch, "paras": out}, open("/home/user/MMM/results_v4/S11/registry_diff.json", "w"), indent=1)
print("paragraphs with changed values (tables excluded):", len(out))
