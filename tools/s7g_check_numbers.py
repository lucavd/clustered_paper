#!/usr/bin/env python3
"""S7 gate G2: digits left in the manuscript and supplement rendered in trace mode (MMM_S7_TRACE=1: every number
read from the registry becomes the sign "§"). Every remaining digit run must belong to a class fixed in the
pre-registration (docs/steps/S7.md): (a) references and citation markers, (b) labels, (c) grid levels written
as levels, (d) nominal levels and design thresholds, (e) mathematical notation; digits inside identifiers
(lme4, AR1, H0) are reported as identifiers. Anything else is a violation. Usage: s7_check_numbers.py trace.md ..."""
import re, sys, csv, collections
NSET = {10, 15, 20, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145}
TSET = {2, 3, 4, 6, 8, 10}
NOMINAL = {"5", "95", "0.05", "2.5", "7.5"}
def classify(text, s, e):
    tok = text[s:e]; before = text[max(0, s - 22):s]; after = text[e:e + 22]
    if s > 0 and re.match(r"[A-Za-z_]", text[s - 1]): return "identifier"
    if before.endswith("AR("): return "identifier"
    lb = text.rfind("{", 0, s); rb = text.find("}", lb) if lb >= 0 else -1
    if lb >= 0 and rb >= e and re.search(r"[NT] ∈ $", text[max(0, lb - 4):lb]) and re.fullmatch(r"[\d,\s…]*", text[lb + 1:rb]):
        return "c-grid-N" if text[lb - 4] == "N" else "c-grid-T"
    if re.match(r"[A-Za-z]", after[:1] or " "): return "identifier"
    if re.search(r"(Section|Sections|Table|Tables|Figure|Figures|Panel|panel)\s+S?[\d.\-, S]*$", before): return "b-label"
    if re.search(r"(^|\n)#{1,4} S?[\d.]*$", before) or re.search(r"(^|\n)\s*$", before) and re.match(r"\.\s", after): return "b-label"
    if re.search(r"\bS$", before) and re.match(r"\.", after): return "b-label"
    if re.search(r"\b[MK]$|\bM\d-M$|\bM$", before): return "c-scenario"
    if re.search(r"\b[A-Z]\s-\s$", before): return "e-math"
    try: v = float(tok)
    except ValueError: v = None
    if tok in NOMINAL and (re.match(r"%|th percentile", after) or (tok == "2.5" and re.match(r" points", after))): return "d-nominal"
    ctxN = re.search(r"(N\s*(=|≤|≥|<|>|∈)|N up to|\bN\b|cluster|subject)", before + " " + after, re.I) is not None
    ctxT = re.search(r"(T\s*(=|≤|≥|<|>|∈)|T up to|\bT\b|occasion|series length)", before + " " + after) is not None
    if v is not None and v == int(v):
        if int(v) in NSET and ctxN: return "c-grid-N"
        if int(v) in TSET and ctxT: return "c-grid-T"
        if re.search(r"N ∈ \{[\d, …]*$", before) and int(v) in NSET: return "c-grid-N"
        if re.search(r"T ∈ \{[\d, ]*$", before) and int(v) in TSET: return "c-grid-T"
    if re.search(r"±\s*$", before) and tok == "2": return "d-nominal"
    return "VIOLATION"
def strip(text):
    i = text.find("## References")
    if i < 0: i = text.find("# References")
    j = text.find("## Tables and figures")
    body = text if i < 0 else (text[:i] + (text[j:] if j > i else ""))
    body = re.sub(r"\$\$.*?\$\$", " MATH ", body, flags=re.S); body = re.sub(r"\$[^$]+?\$", " MATH ", body)
    body = re.sub(r"\\?\[(\d+(?:[,\-–]\d+)*)\\?\]", " CITE ", body)
    body = re.sub(r"<!--.*?-->", " ", body, flags=re.S)
    body = re.sub(r"https?://\S+", " URL ", body)
    body = re.sub(r"!\[[^\]]*\]\([^)]*\)(\{[^}]*\})?", " IMG ", body)
    body = re.sub(r"<img[^>]*>", " IMG ", body)
    body = re.sub(r"^---\n.*?\n---\n", "", body, flags=re.S)
    return body
out = []; counts = collections.Counter()
for f in sys.argv[1:]:
    t = strip(open(f, encoding="utf-8").read())
    for m in re.finditer(r"\d+(?:[.,]\d+)*", t):
        cl = classify(t, m.start(), m.end()); counts[(f, cl)] += 1
        out.append((f, cl, m.group(0), t[max(0, m.start() - 40):m.end() + 30].replace("\n", " ")))
w = csv.writer(open("reports/morris/S7g_logs/S7g_GF2_digits.csv", "w", newline="")); w.writerow(["file", "class", "digits", "context"]); w.writerows(out)
viol = [o for o in out if o[1] == "VIOLATION"]
for k, v in sorted(counts.items()): print(k, v)
print(f"GF2: {len(viol)} violations")
for v in viol[:80]: print("  ", v)
sys.exit(1 if viol else 0)
