#!/usr/bin/env python3
"""tools/okf_registers.py -- S8b: prototype of the S8 review applied to docs/. Original: du1_okf_registers.py -- S8 review, track Dubbi, doubt 1 (codes of other registers in docs/).
Works ONLY on a copy of ~/MMM/docs and ~/MMM/tools/okf.py in ./okf_copy (okf.py check rewrites the indexes).
1. runs `python3 -B tools/okf.py check --bundle docs` on the copy (output: du1_okf_check_before.txt);
2. counts, per file, the codes of other registers that okf.py flags (its own OTHER_CODES_RE, outside code spans
   and links): V-nn (verification log), Qn (expectations), Morris audit requirements ([B-G]n, H2-H6);
   each code is checked against its register (VERIFICATION_LOG_v4.md rows, audit/MORRIS2019_conformita.csv ids);
3. prototype generator: one concept file per verification entry (docs/verification/V-nn.md) and per Morris requirement
   (docs/audit/requirements/<id>.md), then a linker for V-nn (everywhere) and Morris codes (only in lines that speak of
   Morris / requisit* / audit, to avoid S3 block names F1/F2, hypotheses H2... and table cells); Qn are left as they are;
4. reruns the check on the copy and recounts (du1_after_*). Nothing is written in ~/MMM."""
import os, re, csv, sys, shutil, subprocess, importlib.util, collections, pathlib, datetime as dt
SRC = "/home/user/MMM"; W = SRC   # S8b: applied to the real bundle (doubt 1 received from S8)
os.chdir(os.path.join(SRC, "reports/morris/S8b_logs/du1"))
def okf_check(tag):
    r = subprocess.run([sys.executable, "-B", "tools/okf.py", "check", "--bundle", "docs"], cwd=W, capture_output=True, text=True)
    open(f"du1_okf_check_{tag}.txt", "w").write(r.stdout + r.stderr); return r
spec = importlib.util.spec_from_file_location("okf", os.path.join(W, "tools/okf.py")); okf = importlib.util.module_from_spec(spec); spec.loader.exec_module(okf)
B = pathlib.Path(W, "docs")
# ---- registers ----------------------------------------------------------------------------------------------------
VL = {}
for line in (B / "VERIFICATION_LOG_v4.md").read_text(encoding="utf-8").split("\n"):
    m = re.match(r"^\|\s*(V-\d+)\s*\|(.*)$", line)
    if m:
        cells = [c.strip() for c in m.group(2).split("|")]
        VL[m.group(1)] = dict(step=cells[0], number=cells[1], source=cells[2] if len(cells) > 2 else "", command=cells[3] if len(cells) > 3 else "")
MR = {r["id"]: r for r in csv.DictReader(open(B / "audit/MORRIS2019_conformita.csv", encoding="utf-8"))}
def register(code):
    if code.startswith("V-"): return "verification_log"
    if code.startswith("Q"): return "expectation"
    return "morris_requirement"
def count(tag):
    rows, ctx = [], []
    for p in sorted(B.rglob("*.md")):
        fm, body = okf.read_frontmatter(p)
        if body is None: body = p.read_text(encoding="utf-8")
        s = okf._strip_protected(body)
        C = collections.Counter(okf.OTHER_CODES_RE.findall(s))
        for code, n in sorted(C.items()):
            reg = register(code)
            valid = (code in VL) if reg == "verification_log" else (code in MR) if reg == "morris_requirement" else None
            rows.append(dict(file=str(p.relative_to(B)), register=reg, code=code, n=n, in_register=valid))
            for m in list(re.finditer(r"\b" + re.escape(code) + r"\b", s))[:2]:
                ctx.append(dict(file=str(p.relative_to(B)), code=code, context=s[max(0, m.start() - 70):m.end() + 50].replace("\n", " ")))
    with open(f"du1_{tag}_by_file_code.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["file", "register", "code", "n", "in_register"]); w.writeheader(); w.writerows(rows)
    agg = collections.defaultdict(lambda: collections.Counter())
    for r in rows: agg[r["file"]][r["register"]] += r["n"]
    with open(f"du1_{tag}_by_file.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["file", "verification_log", "expectation", "morris_requirement", "total"])
        for k in sorted(agg, key=lambda k: -sum(agg[k].values())):
            a = agg[k]; w.writerow([k, a["verification_log"], a["expectation"], a["morris_requirement"], sum(a.values())])
    if tag == "before":
        with open("du1_contexts.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["file", "code", "context"]); w.writeheader(); w.writerows(ctx)
    tot = collections.Counter(); 
    for r in rows: tot[r["register"]] += r["n"]
    print(tag, "files with codes:", len(agg), "occurrences:", dict(tot), "distinct codes:", len({r['code'] for r in rows}))
    return rows
r0 = okf_check("before"); print(r0.stdout.strip().split("\n")[-2:])
before = count("before")
# ---- prototype generator -----------------------------------------------------------------------------------------
NOW = dt.datetime.now().astimezone().isoformat(timespec="seconds"); GEN = "{ by: mmm-agent/skill-v4, at: %s }" % NOW
q = lambda s: '"' + re.sub(r"\s+", " ", s.replace('"', "'")).strip() + '"'
short = lambda s, n=90: (s if len(s) <= n else s[:n].rsplit(" ", 1)[0] + "...")
os.makedirs(B / "verification", exist_ok=True)
for code, v in VL.items():
    step = re.sub(r"\]\((\.\./)?steps/", "](../steps/", v["step"])
    txt = f"""---
type: Verification Entry
title: {q(f"{code}: " + short(re.sub(r'[`*]', '', v['number'])))}
description: {q("Voce del log di verifica: un numero dei report o del manoscritto, il file che lo ha prodotto e il comando che lo rilegge.")}
status: stable
entry_id: "{code}"
generated: {GEN}
sources:
  - id: verification_log
    resource: "../VERIFICATION_LOG_v4.md"
---
# Numero
{v['number']}

# Step
{step}

# Fonte
{v['source']}

# Comando
{v['command']}
"""
    (B / "verification" / f"{code}.md").write_text(txt, encoding="utf-8")
os.makedirs(B / "audit/requirements", exist_ok=True)
for code, r in MR.items():
    txt = f"""---
type: Requirement
title: {q(f"{code}: " + r['requisito'])}
description: {q(f"Requisito {code} di Morris, White e Crowther (2019), area {r['area_ADEMP']} ({r['rif_Morris2019']}): esito dell'audit iniziale e stato nella versione 4.")}
status: stable
requirement_id: "{code}"
generated: {GEN}
sources:
  - id: morris2019
    resource: https://doi.org/10.1002/sim.8086
  - id: audit_csv
    resource: "../MORRIS2019_conformita.csv"
---
# Requisito
{r['requisito']} (Morris 2019, {r['rif_Morris2019']})[^morris2019]

# Esito dell'audit iniziale
{r['esito']}. {r['evidenza']}

# Azione
{r['azione']}

# Stato nella versione 4
{r['esito_v4']} ({r['step_v4']}). {r['nota_v4']}

[^morris2019]: Morris, White, Crowther (2019), Statistics in Medicine 38:2074-2102.
"""
    (B / "audit/requirements" / f"{code}.md").write_text(txt, encoding="utf-8")
# ---- prototype linker -------------------------------------------------------------------------------------------
VRE = re.compile(r"(?<![\w\[/-])(V-\d{2})(?![\w])")
MRE = re.compile(r"(?<![\w\[/-])([A-H]\d)(?![\w.])")
MCTX = re.compile(r"Morris|requisit|audit|conformit", re.I)
def link_body(text, fdir):
    rel = lambda *a: os.path.relpath(os.path.join(B, *a), fdir)
    out, fence = [], False
    for line in text.split("\n"):
        if line.lstrip().startswith("```"): fence = not fence
        if fence or line.lstrip().startswith("```") or line.lstrip().startswith("#"): out.append(line); continue
        # bracketed citations "[V-53]" or "[V-54, V-55]" (not links) become plain codes first, then are linked
        line = re.sub(r"\[((?:V-\d{2})(?:,\s*V-\d{2})*)\](?!\()", lambda m: m.group(1), line)
        parts = okf.PROTECT_RE.split(line)
        for k in range(0, len(parts), 2):
            p = VRE.sub(lambda m: f"[{m.group(1)}]({rel('verification', m.group(1) + '.md')})" if m.group(1) in VL else m.group(1), parts[k])
            if MCTX.search(line):
                p = MRE.sub(lambda m: f"[{m.group(1)}: {MR[m.group(1)]['requisito']}]({rel('audit/requirements', m.group(1) + '.md')})" if m.group(1) in MR else m.group(1), p)
            parts[k] = p
        out.append("".join(parts))
    return "\n".join(out)
skip = {"index.md"}
TITLES = okf.decision_titles(str(B))
for p in sorted(list(B.glob("verification/V-*.md")) + list(B.glob("audit/requirements/*.md"))):
    # generated concepts: decision codes and steps linked with okf's own linker; codes removed from titles
    txt = p.read_text(encoding="utf-8"); fm, body = txt[4:].split("\n---\n", 1)
    fm = "\n".join(re.sub(r"\bD-\d{2}[bc]?\b", "", l).replace("()", "") if l.startswith(("title:", "description:")) else l for l in fm.split("\n"))
    p.write_text("---\n" + fm + "\n---\n" + okf.link_text(body, TITLES, str(B), str(p.parent)), encoding="utf-8")
for p in sorted(B.rglob("*.md")):
    if p.name in skip: continue
    txt = p.read_text(encoding="utf-8")
    if txt.startswith("---\n") and "\n---\n" in txt[4:]:
        fm, body = txt[4:].split("\n---\n", 1); head = "---\n" + fm + "\n---\n"
    else: head, body = "", txt
    nb = link_body(body, p.parent)
    if nb != body: p.write_text(head + nb, encoding="utf-8")
r1 = okf_check("after"); print(r1.stdout.strip().split("\n")[-2:])
after = count("after")
print("generated:", len(VL), "verification entries,", len(MR), "requirements")
