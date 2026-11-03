#!/usr/bin/env python3
"""S10 track 'assolute': extract every sentence of the rendered manuscript and
supplement that contains an absolute / universal quantifier.
Usage: python3 s10_assolute_extract.py <render_main.md> <render_supp.md> <out.csv>
Sentences are built per paragraph (blank-line separated); table rows (lines
starting with '|') are treated one row = one unit. HTML comments are removed
from the sentence text, and the <!-- claim: ... --> / <!-- noclaim: ... -->
marker that closes a sentence is recorded with it."""
import re, sys, csv

PATTERNS = [
 ("none", r"\bnone\b"), ("never", r"\bnever\b"), ("only", r"\bonly\b"),
 ("every", r"\bevery(?:one|thing|where)?\b"), ("all", r"\ball\b"),
 ("no", r"\bno\b"), ("always", r"\balways\b"), ("any", r"\bany(?:one|thing|where)?\b"),
 ("neither", r"\bneither\b"), ("nor", r"\bnor\b"), ("each", r"\beach\b"),
 ("exactly", r"\bexact(?:ly)?\b"), ("solely", r"\bsolely\b"),
 ("without exception", r"without exception"), ("nowhere", r"\bnowhere\b"),
 ("nothing", r"\bnothing\b"), ("not a single", r"not a single"),
 ("half", r"\b(?:fewer|less|under|below|more|over|above) than half\b|\bunder half\b|\bhalf of\b"),
 ("at most", r"\bat most\b"), ("at least", r"\bat least\b"),
 ("both", r"\bboth\b"), ("zero", r"\bzero\b"),
 ("whole/entire", r"\b(?:whole|entire|entirely)\b"),
 ("invariably/uniformly", r"\b(?:invariably|uniformly|universally|consistently)\b"),
 ("superlative", r"\b(?:smallest|largest|lowest|highest|fewest|most often|least often)\b"),
 ("without", r"\bwithout\b"),
]
RX = [(k, re.compile(p, re.I)) for k, p in PATTERNS]
ABBR = r"(?<!\be\.g)(?<!\bi\.e)(?<!\bet al)(?<!\bvs)(?<!\bFig)(?<!\bEq)(?<!\bcf)(?<!\bSect)(?<!\bapprox)"
SPLIT = re.compile(ABBR + r"(?<=[.!?])(\s*<!--.*?-->)?\s+(?=[A-Z*(\[`$])")
COMMENT = re.compile(r"<!--(.*?)-->", re.S)

def blocks(lines):
    cur, start = [], None
    for i, l in enumerate(lines, 1):
        if l.strip() == "":
            if cur: yield start, cur
            cur, start = [], None
        elif l.lstrip().startswith("|"):
            if cur: yield start, cur
            cur, start = [], None
            yield i, [(i, l.rstrip("\n"))]
        else:
            if start is None: start = i
            cur.append((i, l.rstrip("\n")))
    if cur: yield start, cur

def section_of(lines):
    sec, out = "", {}
    for i, l in enumerate(lines, 1):
        if l.startswith("#"): sec = l.strip("# \n")
        out[i] = sec
    return out

def units(path, part):
    lines = open(path, encoding="utf-8").read().split("\n")
    secs = section_of(lines)
    for start, blk in blocks(lines):
        text, pos2line = "", []
        for ln, l in blk:
            piece = l.strip() + " "
            text += piece; pos2line += [ln] * len(piece)
        is_table = blk[0][1].lstrip().startswith("|")
        if is_table:
            spans = [(0, len(text))]
        else:
            spans, s = [], 0
            for m in SPLIT.finditer(text):
                e = m.end(1) if m.group(1) else m.start()
                spans.append((s, e)); s = m.end()
            spans.append((s, len(text)))
        for s, e in spans:
            raw = text[s:e].strip()
            if not raw: continue
            markers = " | ".join(c.strip() for c in COMMENT.findall(raw))
            sent = re.sub(r"\s+", " ", COMMENT.sub("", raw)).strip()
            words = [k for k, rx in RX if rx.search(sent)]
            if not words: continue
            if is_table and re.fullmatch(r"[|\s:\-+=]*", sent): continue
            ln0 = pos2line[min(s, len(pos2line)-1)]; ln1 = pos2line[min(e-1, len(pos2line)-1)]
            yield dict(part=part, file=path.split("/")[-1], line_start=ln0, line_end=ln1,
                       section=secs.get(ln0, ""), is_table_row=is_table, words=";".join(words),
                       marker=markers, has_claim_marker=("claim:" in markers and "noclaim" not in markers.replace("claim:", "", 0) or bool(re.search(r"(?<!no)claim:", markers))),
                       sentence=sent)

if __name__ == "__main__":
    a, b, out = sys.argv[1:4]
    rows = list(units(a, "main")) + list(units(b, "supplement"))
    for k, r in enumerate(rows, 1):
        r["has_claim_marker"] = bool(re.search(r"(?<!no)claim:", r["marker"]))
        r["has_noclaim_marker"] = "noclaim:" in r["marker"]
        r["cand_id"] = f"A{k:04d}"
    cols = ["cand_id","part","file","line_start","line_end","section","is_table_row","words","has_claim_marker","has_noclaim_marker","marker","sentence"]
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({c: r[c] for c in cols})
    print(len(rows))
