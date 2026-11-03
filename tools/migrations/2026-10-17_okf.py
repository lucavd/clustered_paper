#!/usr/bin/env python3
"""One-off migration of the MMM documentation to Open Knowledge Format v0.2.

2026-10-17, session protocol v3. Run from the repo root on a clean tree:
  python3 tools/migrations/2026-10-17_okf.py --meta tools/migrations/2026-10-17_decision_meta.json

What it does (text of existing documents is unchanged except where listed):
  1. DECISIONS.md (single table) -> docs/decisions/D-nn.md, one concept per decision,
     text and rationale verbatim; titles and descriptions from the metadata file,
     reviewed and approved by the PI on 2026-10-17 (verified). DECISIONS.md becomes a pointer.
  2. Two new decisions: D-21 (plain-language summary at the top of every report) and
     D-22 (session protocol v3); D-09 becomes deprecated, superseded by D-22.
  3. docs/steps/S0..S4.md: one small Step concept per closed step (report, session note,
     PI approval with the commit that records it). Session notes stay in docs/sessions/.
  4. OKF frontmatter on roadmap, protocol, verification log, session notes, audit report.
  5. Roadmap: sections 1-2 (rules and ritual of protocol v1) replaced by a pointer to the
     protocol v3; the stale decision status line of section 3 replaced by a pointer.
  6. Decision codes become links with the decision title (everywhere); step codes become
     links (not in the roadmap, which defines the steps, nor in headings).
  7. docs/log.md, indexes, conformance check (tools/okf.py check).
"""
import argparse
import datetime
import json
import pathlib
import re
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import okf  # noqa: E402

AGENT, LEGACY, PI = "mmm-agent/skill-v3", "mmm-agent/skill-v1", "human:luca-vedovelli"
BUNDLE = pathlib.Path("docs")
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def git(*args, default=None):
    try:
        return subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout.strip() or default
    except Exception:
        return default


def last_commit_time(path):
    return git("log", "-1", "--format=%aI", "--", str(path), default=NOW)


def yq(s):
    return json.dumps(s, ensure_ascii=False)


def slug(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def split_row(line):
    s = line.strip()
    s = s[1:-1] if s.endswith("|") else s[1:]
    cells, cur, in_code = [], [], False
    for i, ch in enumerate(s):
        if ch == "`":
            in_code = not in_code
        if ch == "|" and not in_code and (i == 0 or s[i - 1] != "\\"):
            cells.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    cells.append("".join(cur).strip())
    return cells


def add_frontmatter(path, lines):
    path = pathlib.Path(path)
    txt = path.read_text(encoding="utf-8")
    assert not txt.startswith("---\n"), f"{path} already has frontmatter"
    path.write_text("---\n" + "\n".join(lines) + "\n---\n" + txt, encoding="utf-8")


def write_concept(path, fm_lines, body):
    pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(path).write_text("---\n" + "\n".join(fm_lines) + "\n---\n" + body, encoding="utf-8")


MORRIS = ['  - id: morris2019', '    resource: https://doi.org/10.1002/sim.8086',
          '    title: "Morris, White, Crowther (2019). Using simulation studies to evaluate statistical methods. Statistics in Medicine 38:2074-2102"']

NEW_DECISIONS = {
    "D-21": dict(
        group="Organizzazione del lavoro e dei file",
        title="Riassunto in parole semplici in testa ai report",
        description="Ogni report di step si apre con un riquadro In breve (TL;DR) in parole semplici, tipo abstract, perche' il PI non debba rileggere tutte le tabelle per capire l'esito.",
        sources=['  - id: commit-tldr', '    resource: "git:morris-v4@b035448"',
                 '    title: "Commit che ha aggiunto il riquadro ai report dei primi quattro step"',
                 '  - id: riga-originale', '    resource: "git:morris-v4@b035448:DECISIONS.md"',
                 '    title: "Riga originale nel registro unico, prima della sovrascrittura"'],
        body=(
            "# Decisione\n\n"
            "Ogni report di step si apre con un riquadro \"In breve (TL;DR)\", in parole semplici e tipo abstract, "
            "prima delle sezioni tecniche: domanda, cosa abbiamo fatto, risultato, sorprese o attese non superate, "
            "cosa cambia per il paper, prossimo passo. Al massimo una schermata; si scrive per ultimo, a esito noto. "
            "Il riquadro e' stato aggiunto a posteriori ai report dei primi quattro step (S0-S3) il 2026-10-17, "
            "senza toccare il resto.[^commit-tldr]\n\n"
            "# Motivazione\n\n"
            "Richiesta del PI in chat, 2026-10-17: non dover rileggere tutte le tabelle a ogni report. La riga scritta "
            "allora nel registro unico e' stata sovrascritta dal commit `bae4915`, che ha assegnato lo stesso numero a "
            "un'altra decisione (D-18); questa decisione la ripristina con un numero nuovo.[^riga-originale]\n\n"
            "[^commit-tldr]: Commit che ha aggiunto il riquadro ai report dei primi quattro step\n"
            "[^riga-originale]: Riga originale nel registro unico, prima della sovrascrittura\n"
        ),
    ),
    "D-22": dict(
        group="Organizzazione del lavoro e dei file",
        title="Protocollo di sessione v3",
        description="Previsioni giudicate con un intervallo di rilevanza, step chiuso solo senza dubbi aperti (rinvio ammesso), linguaggio semplice e preciso per il PI e documentazione in formato OKF; sostituisce la regola sulle riserve.",
        sources=['  - id: skill-v3', '    resource: "skill mmm-session v3, catalogo delle skill del PI"',
                 '    title: "Protocollo di sessione v3"',
                 '  - id: okf', '    resource: https://github.com/GoogleCloudPlatform/open-knowledge-format/blob/main/SPEC.md',
                 '    title: "Open Knowledge Format v0.2"',
                 '  - id: bradley1978', '    resource: https://doi.org/10.1111/j.2044-8317.1978.tb00581.x',
                 '    title: "Bradley (1978). Robustness? British Journal of Mathematical and Statistical Psychology 31:144-152"',
                 '  - id: q9', '    resource: "git:morris-v4@47d4572:tables/V4_S3_q9_h0_paired.csv"',
                 '    title: "Differenze appaiate usate per provare la regola di rilevanza"'],
        body=(
            "# Decisione\n\n"
            "Il protocollo di sessione v3 (skill `mmm-session`, pubblicata il 2026-10-17) sostituisce D-09 e cambia "
            "cinque cose.[^skill-v3]\n\n"
            "1. **Attese in due classi.** Cancelli: al massimo 5 controlli binari che, se falliscono, invalidano lo step. "
            "Previsioni: ciascuna con un intervallo di rilevanza, cioe' lo scarto oltre il quale cambierebbe una frase "
            "del paper o una scelta della roadmap; si giudicano con la stima +/- 2 MCSE: confermata se tutta dentro, "
            "smentita se tutta fuori, indeterminata se a cavallo. Margini di default dal criterio liberale di "
            "Bradley:[^bradley1978] errore di I tipo da 2,5% a 7,5%; copertura da 92,5% a 97,5%; errore relativo "
            "dell'SE del modello da -9% a +14%; bias entro 0,46 SE empirici; differenze fra versioni o metodi "
            "+/- 2,5 punti per cella e +/- 1 punto sulla media delle celle.\n"
            "2. **Chiusura senza riserve.** Verdetto chiuso o non chiuso. L'agente prova per primo a chiudere ogni "
            "dubbio; un dubbio e' chiuso quando e' spiegato, corretto, deciso dal PI o rinviato (tollerato, come "
            "ultima opzione, e scritto nel file dello step che lo riceve).\n"
            "3. **Chi decide.** Il PI decide cio' che cambia estimando, generatori dei dati, metodi, griglia, testo "
            "del paper, o costa piu' di circa 4 ore di calcolo; il resto e' una decisione tecnica dell'agente.\n"
            "4. **Linguaggio.** Semplice ma preciso, per un lettore che sul dettaglio del progetto e' un passo "
            "indietro rispetto a chi scrive; nessun codice interno senza la sua spiegazione; sigle per esteso alla "
            "prima occorrenza.\n"
            "5. **Documentazione in formato Open Knowledge Format v0.2**[^okf] nella cartella `docs/`: un concetto "
            "per file con frontmatter, una decisione per file, indici generati, controllo con `tools/okf.py` prima "
            "di ogni commit.\n\n"
            "# Motivazione\n\n"
            "Richiesta del PI in chat, 2026-10-17: il protocollo precedente creava piu' problemi di quanti ne "
            "risolvesse. Nei primi quattro step 8 attese erano risultate non superate, circa meta' per come erano "
            "formulate, e la regola sulle riserve chiedeva di risolverle una per una prima di andare avanti, spesso "
            "con una decisione del PI (per esempio D-17).\n\n"
            "Prova della regola di rilevanza su un'attesa reale del rerun correttivo: effetto della rigenerazione "
            "dei dati sotto H0 di GenBinaryY sull'errore di I tipo, scenari M3-M4. Con la regola precedente (entro "
            "3 MCSE) l'attesa falliva: 1.190 combinazioni cella per metodo su 1.347 entro la soglia. Con i margini "
            "nuovi (+/- 2,5 punti per cella), sui 7 metodi diversi dal LRT 1.165 combinazioni sono confermate, 20 "
            "indeterminate e nessuna smentita; il LRT e' smentito in 88 combinazioni su 162, come atteso dal cambio "
            "del suo modello (D-08).[^q9]\n\n"
            "[^skill-v3]: Protocollo di sessione v3\n"
            "[^bradley1978]: Bradley (1978), British Journal of Mathematical and Statistical Psychology\n"
            "[^okf]: Open Knowledge Format v0.2\n"
            "[^q9]: Differenze appaiate usate per provare la regola di rilevanza\n"
        ),
    ),
}

SESSIONS = {
    "S0": ("2026-10-15_S0.md", "Nota di sessione: baseline sul server e protocollo"),
    "S1": ("2026-10-15_S1.md", "Nota di sessione: protocollo ADEMP di revisione"),
    "S2": ("2026-10-16_S2.md", "Nota di sessione: diagnostica sui risultati esistenti"),
    "S3": ("2026-10-16_S3.md", "Nota di sessione: rerun correttivo e scenario ancorato alla coorte"),
    "S4": ("2026-10-17_S4.md", "Nota di sessione: riaggregazione unica con MCSE"),
}
STEP_CONCEPTS = {
    "S0": ("Baseline sul server e protocollo di sessione",
           "Prima di correggere qualcosa: registrate le impronte md5 dei 28.800 risultati grezzi, delle 34 tabelle e del manoscritto; nessun file alterato, manoscritto identico a quello del Mac.",
           ("2026-10-15T21:55:06+02:00", "b0826c5")),
    "S1": ("Protocollo ADEMP datato prima dei rerun",
           "Ricostruiti dal codice i parametri che hanno generato i dati, verificati con una rigenerazione identica (240 dataset su 240), e trovati tre difetti da correggere: generatore casuale, parita' fra H0 e H1, modello del LRT.",
           ("2026-10-15T22:36:53+02:00", "010028d")),
    "S2": ("Diagnostica sui risultati esistenti, senza nuove simulazioni",
           "I replicati mancanti sono disegni non stimabili (247 per ipotesi), i valori veri sono stati ricalcolati, l'esclusione dei fit con warning cambia l'errore di I tipo del glmer e la coorte reale e' lontana dal modello simulato.",
           ("2026-10-16T14:18:48+02:00", "b380452")),
    "S3": ("Rerun correttivo e scenario ancorato alla coorte",
           "18.720 nuovi chunk senza errori e con i semi intatti; il LRT con le pendenze casuali e' meno liberale e il Wald del glmer con le opzioni di default non sempre raggiunge il massimo.",
           ("2026-10-17T18:21:10+02:00", "bae4915")),
    "S4": ("Riaggregazione unica con MCSE",
           "Tutte le misure di prestazione ricalcolate con MCSE su 3.024 celle; 11 attese su 14 superate, la versione precedente riprodotta esattamente dove i dati non cambiano, il GLMM con pendenze casuali non stimabile nelle serie piu' corte.",
           None),
}
ROADMAP_RULES = (
    "## 1-2. Regole di lavoro e rituale\n\n"
    "Le regole di lavoro e il rituale di apertura e chiusura sono nel protocollo di sessione v3 (skill `mmm-session`), "
    "registrato come D-22. Dal 2026-10-17 sostituiscono le regole e il rituale scritti qui, che restano nella storia "
    "git (commit `47d4572`).\n\n"
)
ROADMAP_STATUS_LINE = (
    "Tabella storica delle opzioni proposte prima del rerun correttivo. Lo stato attuale di ogni decisione e' "
    "nell'[indice delle decisioni](decisions/index.md)."
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    a = ap.parse_args()
    meta = json.loads(pathlib.Path(a.meta).read_text(encoding="utf-8"))
    head = git("rev-parse", "--short", "HEAD", default="47d4572")

    # 1. decisions from the single-table register
    rows = [split_row(l) for l in pathlib.Path("DECISIONS.md").read_text(encoding="utf-8").splitlines() if l.startswith("| D-")]
    assert all(len(r) == 5 for r in rows) and len(rows) == len(meta), (len(rows), len(meta))
    for code, date, decision, rationale, state in rows:
        m = meta[code]
        fm = ["type: Decision", f"title: {yq(m['title'])}", f"description: {yq(m['description'])}",
              f"decision_id: {yq(code)}", f"decided_on: {date}",
              f"decided_by: {yq(AGENT if 'tecnica' in state else PI)}",
              f"status: {'deprecated' if code == 'D-09' else 'stable'}"]
        if code == "D-09":
            fm.append('superseded_by: "D-22.md"')
        fm += [f"group: {yq(m['group'])}", f"tags: [decisione, {slug(m['group'])}]",
               f"generated: {{ by: {AGENT}, at: {NOW} }}", f"verified: {{ by: {PI}, at: {NOW} }}",
               "sources:", "  - id: registro", f'    resource: "git:morris-v4@{head}:DECISIONS.md"',
               '    title: "Registro delle decisioni, tabella unica prima della conversione"']
        body = (f"# Decisione\n\n{decision}\n\n# Motivazione\n\n{rationale}\n\n"
                f"Stato nel registro originale: {state}.[^registro]\n\n"
                "[^registro]: Registro delle decisioni, tabella unica prima della conversione\n")
        if code == "D-09":
            body = body.replace("# Motivazione", "Superata dal protocollo di sessione v3 (D-22).\n\n# Motivazione", 1)
        write_concept(BUNDLE / "decisions" / f"{code}.md", fm, body)

    # 2. new decisions
    for code, d in NEW_DECISIONS.items():
        fm = ["type: Decision", f"title: {yq(d['title'])}", f"description: {yq(d['description'])}",
              f"decision_id: {yq(code)}", "decided_on: 2026-10-17", f"decided_by: {yq(PI)}", "status: stable",
              f"group: {yq(d['group'])}", f"tags: [decisione, {slug(d['group'])}]",
              f"generated: {{ by: {AGENT}, at: {NOW} }}", "sources:", *d["sources"]]
        write_concept(BUNDLE / "decisions" / f"{code}.md", fm, d["body"])

    # 3. session notes: frontmatter only (text unchanged except links, step 6)
    for sid, (fname, title) in SESSIONS.items():
        p = BUNDLE / "sessions" / fname
        add_frontmatter(p, ["type: Session Note", f"title: {yq(title)}",
                            'description: "Diario di lavoro dello step: attese scritte prima del calcolo, cosa e\' stato fatto, cosa e\' stato trovato e come andare avanti."',
                            f'step_id: "{sid}"', "status: stable",
                            f"generated: {{ by: {LEGACY}, at: {last_commit_time(p)} }}"])

    # 3b. step concepts
    for sid, (title, desc, appr) in STEP_CONCEPTS.items():
        fname = SESSIONS[sid][0]
        fm = ["type: Step", f"title: {yq(title)}", f"description: {yq(desc)}", f'step_id: "{sid}"', "status: stable",
              f'resource: "../../reports/morris/{sid}.html"']
        fm.append(f'pi_approval: {{ by: {PI}, at: {appr[0]}, evidence: "commit {appr[1]}" }}' if appr else "pi_approval: pending")
        fm += [f"generated: {{ by: {AGENT}, at: {NOW} }}", "sources:",
               "  - id: nota", f'    resource: "../sessions/{fname}"', '    title: "Nota di sessione"',
               "  - id: report", f'    resource: "../../reports/morris/{sid}.qmd"', '    title: "Sorgente del report"']
        appr_txt = (f"Approvato dal PI: il commit `{appr[1]}` lo registra." if appr else "In attesa dell'ok del PI.")
        body = ("# Esito\n\n"
                f"Il riassunto in parole semplici e' in testa al [report dello step](../../reports/morris/{sid}.html);[^report] "
                f"attese, lavoro svolto e riferimenti sono nella [nota di sessione](../sessions/{fname}).[^nota] {appr_txt}\n\n"
                "[^report]: Sorgente del report\n[^nota]: Nota di sessione\n")
        write_concept(BUNDLE / "steps" / f"{sid}.md", fm, body)

    # 4-5. roadmap, protocol, verification log, audit report
    rm = BUNDLE / "ROADMAP_MORRIS.md"
    t = rm.read_text(encoding="utf-8")
    i, j = t.index("## 1. Regole non negoziabili"), t.index("## 3. Decisioni richieste al PI")
    assert t.count("## 1. Regole non negoziabili") == 1 and i < j
    t = t[:i] + ROADMAP_RULES + t[j:]
    stale = [l for l in t.splitlines() if l.startswith("Stato al 2026-10-16:")]
    assert len(stale) == 1
    t = t.replace(stale[0], ROADMAP_STATUS_LINE)
    rm.write_text(t, encoding="utf-8")
    add_frontmatter(rm, ["type: Roadmap", 'title: "Roadmap verso la conformita\' a Morris 2019"',
                         'description: "I nove step per portare simulazione e manoscritto a conformita\' con Morris, White e Crowther (2019): stato, prodotti, dipendenze e attese di ciascuno."',
                         "status: stable", f"generated: {{ by: {AGENT}, at: {NOW} }}", "sources:", *MORRIS,
                         "  - id: audit", '    resource: "audit/MORRIS2019_conformita.csv"',
                         '    title: "Audit dei 28 requisiti, punto di partenza"'])
    pr = BUNDLE / "PROTOCOL_ADEMP.md"
    add_frontmatter(pr, ["type: Protocol", 'title: "Protocollo ADEMP della revisione dello studio di simulazione"',
                         'description: "Obiettivi, meccanismi generatori, estimandi, metodi e misure di prestazione secondo lo schema ADEMP di Morris 2019, datato e versionato prima di ogni rerun (versione 1.5)."',
                         "status: stable", f"generated: {{ by: {LEGACY}, at: {last_commit_time(pr)} }}", "sources:", *MORRIS,
                         "  - id: parametri", '    resource: "../tables/V4_dgm_parameters.csv"',
                         '    title: "Tabella dei parametri estratta dal codice eseguito"'])
    vl = BUNDLE / "VERIFICATION_LOG_v4.md"
    add_frontmatter(vl, ["type: Verification Log", 'title: "Log di verifica dei numeri della revisione"',
                         'description: "Ogni numero usato nei report e nel manoscritto della revisione, con il file che lo ha prodotto e il comando che lo rilegge (voci V-nn)."',
                         "status: stable", f"generated: {{ by: {LEGACY}, at: {last_commit_time(vl)} }}"])
    au = BUNDLE / "audit" / "MORRIS2019_conformita_report.md"
    add_frontmatter(au, ["type: Audit Report", 'title: "Audit di conformita\' a Morris 2019, punto di partenza"',
                         'description: "Esito dell\'audit iniziale dei 28 requisiti di Morris, White e Crowther (2019) sul manoscritto precedente, con il dettaglio per requisito."',
                         "status: stable", 'resource: "MORRIS2019_conformita.csv"', "sources:", *MORRIS])

    # 6. links
    titles = okf.decision_titles(BUNDLE)
    for p in sorted((BUNDLE / "decisions").glob("D-*.md")):
        okf.link_file(p, titles, BUNDLE)
    for fname, _ in SESSIONS.values():
        okf.link_file(BUNDLE / "sessions" / fname, titles, BUNDLE)
    okf.link_file(pr, titles, BUNDLE)
    okf.link_file(vl, titles, BUNDLE)
    okf.link_file(rm, titles, BUNDLE, link_steps=False)

    # DECISIONS.md pointer and log
    pathlib.Path("DECISIONS.md").write_text(
        "# DECISIONS\n\n"
        "Le decisioni del progetto sono in `docs/decisions/`, una per file, con titolo e descrizione in parole semplici. "
        "Indice: [docs/decisions/index.md](docs/decisions/index.md).\n\n"
        f"Fino al commit `{head}` (2026-10-17) questo file conteneva la tabella unica delle prime 24 decisioni; "
        "la conversione e' registrata in [docs/log.md](docs/log.md).\n", encoding="utf-8")
    (BUNDLE / "log.md").write_text(
        "# Registro delle modifiche\n\n"
        f"## {NOW[:10]}\n"
        "* **Creazione**: documentazione convertita in Open Knowledge Format v0.2 con il protocollo di sessione v3 "
        f"([{NEW_DECISIONS['D-22']['title']} (D-22)](decisions/D-22.md)). Le 24 decisioni del registro unico sono "
        "ora un file ciascuna in `decisions/`, con titolo e descrizione riletti e approvati dal PI; due decisioni "
        f"nuove ([{NEW_DECISIONS['D-21']['title']} (D-21)](decisions/D-21.md) e D-22, non ancora rilette). "
        "Un file per step in `steps/`. Frontmatter su roadmap, protocollo, log di verifica, note di sessione e audit; "
        "nel protocollo, nel log di verifica e nelle note il testo e' invariato salvo i link ai codici.\n"
        "* **Aggiornamento**: nella roadmap le regole e il rituale del protocollo precedente sono sostituiti da un "
        "rimando al protocollo v3, e la riga sullo stato delle decisioni da un rimando all'indice.\n",
        encoding="utf-8")
    okf.link_file(BUNDLE / "log.md", titles, BUNDLE, link_steps=False)

    sys.exit(1 if okf.check(BUNDLE) else 0)


if __name__ == "__main__":
    main()
