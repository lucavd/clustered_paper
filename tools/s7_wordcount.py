#!/usr/bin/env python3
"""S7 prediction P1/P2: words of the main text (Introduction to Conclusions; abstract, tables, captions, references,
availability and declarations excluded) plus 200 per table or figure; sentence statistics of the same text.
Usage: s7_wordcount.py MMM_v4.md (rendered, numbers substituted)"""
import re, sys, statistics
t = open(sys.argv[1], encoding="utf-8").read()
i0 = t.find("## 1. Introduction"); i1 = t.find("## Code and data availability")
body = t[i0:i1]
body = re.sub(r"<!--.*?-->", " ", body, flags=re.S)
lines = [l for l in body.splitlines() if not l.startswith(("|", "#", "**Table", "**Figure", "<img", "![")) and l.strip()]
txt = " ".join(lines)
txt = re.sub(r"\$\$.*?\$\$", " ", txt); txt = re.sub(r"\$[^$]+\$", " x ", txt)
words = re.findall(r"[A-Za-z0-9][A-Za-z0-9\-'.,%:/()]*", txt)
n_disp = len(re.findall(r"^\*\*(Table|Figure) \d+\.\*\*", t, re.M))
sents = [s for s in re.split(r"(?<=[.;?])\s+(?=[A-Z(`])", txt) if len(s.split()) > 2]
L = [len(s.split()) for s in sents]
abst = t[t.find("## Abstract"):t.find("**Keywords")]
abst = re.sub(r"<!--.*?-->", " ", abst).replace("## Abstract", "")
print(f"P1: words {len(words)} + 200 x {n_disp} display items = {len(words) + 200 * n_disp}")
print(f"P2: sentences {len(L)}, median {statistics.median(L)}, mean {statistics.mean(L):.1f}, p90 {sorted(L)[int(0.9 * len(L))]}, max {max(L)}")
print(f"abstract words {len(re.findall(r'[A-Za-z0-9][A-Za-z0-9\-.,%]*', abst))}")
print(f"em dashes {t.count(chr(8212))}, en dashes {t.count(chr(8211))}")
