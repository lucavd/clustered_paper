#!/usr/bin/env python3
"""Map each candidate sentence (rendered text) to the registry ids that
produce its numbers in the qmd source. Usage: map_ids.py cand.csv main.qmd supp.qmd numbers.csv out.csv"""
import re, sys, csv, pandas as pd
cand, mq, sq, reg, out = sys.argv[1:6]
R = pd.read_csv(reg, dtype=str).set_index("id")["text"].to_dict()
def norm(s):
    s = re.sub(r"<!--.*?-->", "", s)
    s = s.replace("\\", "").replace("*", "").replace("’", "'").replace("“", '"').replace("”", '"')
    return re.sub(r"\s+", " ", s).strip()
TOK = re.compile(r"`r n\(\"([^\"]+)\"\)`")
def render_lines(path):
    res = []
    for ln, l in enumerate(open(path, encoding="utf-8"), 1):
        txt, ids, pos = "", [], 0
        for m in TOK.finditer(l):
            txt += l[pos:m.start()]
            v = R.get(m.group(1), "??")
            ids.append((len(norm(txt)), m.group(1)))
            txt += v; pos = m.end()
        txt += l[pos:]
        res.append((ln, norm(txt), ids, l))
    return res
SRC = {"main": render_lines(mq), "supplement": render_lines(sq)}
d = pd.read_csv(cand)
src_line, src_ids = [], []
for r in d.itertuples():
    s = norm(r.sentence)
    hit = None
    for ln, t, ids, raw in SRC[r.part]:
        # recompute ids positions on norm(text before token) is approximate; search sentence
        k = t.find(s[:60]) if len(s) >= 20 else t.find(s)
        if k >= 0:
            hit = (ln, [i for p, i in ids if k - 2 <= p <= k + len(s) + 2]); break
    src_line.append(hit[0] if hit else None); src_ids.append(";".join(hit[1]) if hit else "")
d["qmd_line"] = src_line; d["reg_ids"] = src_ids
d.to_csv(out, index=False)
print(d.qmd_line.isna().sum(), "unmatched;", (d.reg_ids != "").sum(), "with ids")
print(d[d.qmd_line.isna()][["cand_id","part","line_start"]].to_string())
