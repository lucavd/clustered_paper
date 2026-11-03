---
type: Step
title: "Attuazione della revisione finale: correzioni, revisione di stile e accorciamento"
description: "Applicare al manoscritto v5 il piano di implementazione della revisione finale, con le decisioni del PI del 2026-11-01, poi la revisione di stile e l'accorciamento."
step_id: "S8b"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-11-02T09:10:00+01:00, evidence: "decisioni in chat del 2026-11-02 e commit di chiusura dello step" }
generated: { by: mmm-agent/skill-v4, at: 2026-11-01T16:15:00+01:00 }
sources:
  - id: piano
    resource: "../../reports/morris/S8_logs/PIANO_IMPLEMENTAZIONE_S8.md"
    title: "Piano di implementazione S8 (64 voci)"
---
# In breve

## In breve: la revisione finale è stata applicata al manoscritto, e il testo è pronto per la tua rilettura?

**A cosa serviva.** Applicare al manoscritto tutte le correzioni chieste dalla revisione finale, rifare le due applicazioni cliniche con il modello marginalizzato multilivello che ha una variabile latente, accorciare il testo principale entro il limite della rivista, e far rileggere da un revisore indipendente le parti cambiate. Serve perché il manoscritto si possa inviare senza numeri sbagliati né affermazioni senza sostegno.

**Cosa abbiamo trovato.**

1. **Tutte le correzioni chieste sono state applicate, salvo tre motivate.**
   - Due non servivano più; una (i collegamenti interni della documentazione) è corretta in parte, il resto è spiegato.
   - Le sei fonti che il revisore aveva letto solo nel riassunto sono state lette per intero: confermano le citazioni; un esempio della sezione sui limiti è stato allineato alle fonti.
2. **Tre errori di calcolo nella catena dei numeri sono corretti, e nessun confronto fra metodi cambia verso.**
   - Riguardavano la potenza aggiustata per la taratura, i denominatori del confronto diretto fra le due strategie raccomandate e la tabella dei parametri di simulazione.
3. **Nelle applicazioni cliniche l'errore standard del modello è troppo piccolo sia per il modello misto sia per il modello marginalizzato, in misura simile.**
   - Il modello misto qui è il modello logistico a effetti casuali per paziente.
   - L'errore standard robusto del modello marginalizzato è vicino alla dispersione delle stime nella coorte neonatale, non nel trial respiratorio con pochi pazienti.
4. **Il testo principale è sceso a 8.111 parole contate, dentro il limite.**
   - Una tabella del testo è passata nel supplemento; il supplemento ha due tabelle e una figura nuove, e l'errore Monte Carlo in ogni didascalia.
5. **La rilettura indipendente ha trovato cinque affermazioni della discussione e delle conclusioni rimaste senza sostegno nei risultati dopo l'accorciamento: i valori sono tornati nei risultati.**
   - Gli altri rilievi, tutti minori o di formulazione, sono stati corretti uno per uno; tutti i requisiti di Morris per la conduzione e il resoconto degli studi di simulazione risultano soddisfatti.

**Cosa devi decidere.** Le due decisioni chieste (abstract più prudente sul conteggio zero, numeri dell'intervallo del test a inversione di segno solo nei risultati) sono state prese il 2 novembre e applicate.

**Comunque.** Prima dell'invio restano fuori dal testo scientifico: le dichiarazioni (finanziamento, conflitti di interesse, approvazione etica, uso di strumenti software), il codice pubblico con DOI e il commit citato, la riconferma del formato della rivista e la lettera di accompagnamento.

Non c'è nulla da decidere: se non hai dubbi, lo step è chiuso.

# Da leggere per primo

- Piano: `reports/morris/S8_logs/PIANO_IMPLEMENTAZIONE_S8.md`, voci strutturate in `ledger_S8.csv`, generatore `s8_mkplan.py`, note dei revisori del piano in `arbitrati_S8.csv` (26 note da leggere prima delle voci a cui si riferiscono).
- Decisioni del PI che il piano chiama DEC-1, DEC-2, DEC-3: [D-70](../decisions/D-70.md), [D-71](../decisions/D-71.md), [D-72](../decisions/D-72.md).
- Ordine di lavoro (indicazione del PI e [D-69](../decisions/D-69.md)): correzioni della catena e del registro, poi testo, poi accorciamento, poi secondo passaggio del revisore n.2, poi revisione di stile (sezione sotto).

# Indicazione del PI (2026-11-01): secondo passaggio del revisore n.2 prima dello stile

Ordine di lavoro: (1) correzioni della catena e del registro; (2) correzioni del testo; (3) accorciamento; (4) **secondo passaggio del revisore n.2** (skill `revisione-e-piano-implementazione`) sul testo corretto e accorciato; (5) revisione di stile (`manuscript-style-review`).

Il secondo passaggio copre:

- tutto cio' che e' cambiato: numeri nuovi o ricalcolati (potenza aggiustata, applicazioni cliniche rifatte, id nuovi del registro), frasi nuove (abstract, 4.1, 4.3, 4.5, conclusioni), frasi spostate nel supplemento (devono restare sostenute e attestate), tabelle e figure nuove;
- il controllo che ogni affermazione dell'abstract e delle conclusioni resti sostenuta da una frase dei risultati dopo l'accorciamento;
- la catena automatica su tutti i numeri (registro, cifre a mano, etichette, conteggio delle parole, formato Word).

Le parti non toccate dalle correzioni non si rifanno da capo: per quelle vale la revisione di S8, salvo che un numero da cui dipendono sia cambiato.

# Dubbi ricevuti

Rinviati il 2026-11-01 da [Revisione avversariale e ri-audit Morris (S8)](S8.md).

1. **Codici di altri registri ancora senza spiegazione** (ricevuto da [manoscritto (S7)](S7.md) e prima da [D-22](../decisions/D-22.md)). Cosa: un documento per requisito Morris e per voce del log di verifica, poi `tools/okf.py link`; 509 codici senza link in 34 documenti. Per chiuderlo: lo script prototipo in `~/MMM_rev_S8/Dubbi/` sul server.
2. **Figura dei falsi positivi in M6-M7**. Cosa: Figura [protocollo datato (S1)](S1.md) dalla sorgente fewclust (`tables/V5_S7f_rule_dp.csv`), coerente con la Tabella 4, senza titoli e con i flag esatti (voce F02 del piano). Per chiuderlo: script in `py/s8b/`, riga nella pipeline, didascalia, valori nel registro.
3. **Bradley (1978)**. Cosa: confermare sul testo primario il criterio liberale 0,5-1,5 alfa (oggi confermato da due fonti secondarie). Per chiuderlo: il PDF dal PI o dalla biblioteca.

# Attese pre-registrate (2026-11-01, prima di ogni calcolo)

Le attese riguardano le correzioni della catena, le applicazioni cliniche rifatte secondo [D-70](../decisions/D-70.md), l'accorciamento secondo [D-71](../decisions/D-71.md) e il secondo passaggio del revisore n.2. Il manoscritto si modifica sul posto (`manuscript/v5/`); la versione di partenza resta nel commit di questa pre-registrazione.

## Cancelli

| Codice | Cancello | Come si controlla |
|---|---|---|
| GB1 | Integrita' | md5 di `results_full/` prima e dopo uguali alla lista di S7b (28.800 file, 0 differenze; `R/s7g/s7g_md5.sh` in copia in `R/s8b/`) |
| GB2 | Catena dei numeri | alla chiusura: `tools/s7g_verify_registry.py` (in copia) rilegge tutte le voci del registro dai file sorgente con 0 discordanze; il controllo dei numeri del render (GF2) non trova cifre senza fonte; il controllo delle etichette delle frasi (GF4) passa; nessuna cifra nuova scritta a mano nel qmd (ogni numero nuovo ha un id con file, filtro e campo) |
| GB3 | Controprova indipendente delle correzioni | (a) i 79 id della potenza aggiustata prodotti dalla catena corretta coincidono con `corrected_value` di `adjpower_fix_S8.csv` (ricalcolo del revisore dalle righe) entro 1e-8 e con il controllo indipendente `R/s7f/s7f_table5_recompute.R` corretto; (b) dopo la correzione dei denominatori le quattro celle della Tabella S31 valgono 19,15, 38,58, 31,00 e 69,55 entro 0,01 (ricalcolo del revisore); (c) il modello marginalizzato con variabile latente sul dato completo della coorte neonatale stima 0,352 (tre decimali, misura del revisore) |
| GB4 | Copertura del piano | ognuna delle 64 voci di `ledger_S8.csv` ha un esito nel registro di attuazione `reports/morris/S8b_logs/attuazione_S8b.csv` (fatta, non applicata con motivo, rinviata con lo step che la riceve); per le voci di testo un controllo da script verifica che l'ancora non compaia piu' (o resti per un motivo scritto) e che il testo nuovo compaia nel render; le voci legate a rilievi bloccanti (abstract, 4.1, conclusioni: P23, P40, P40b) sono fatte prima del render finale |
| GB5 | Formato e bibliografia | abstract al massimo 200 parole; controllo del formato Word superato per manoscritto e supplemento; controllo delle citazioni: ogni numero citato ha la sua voce, ogni voce e' citata, numerazione in ordine di prima citazione; render docx e html di manoscritto e supplemento senza errori |

## Previsioni

| Codice | Quantita' | Atteso | Intervallo di rilevanza | Perche' conta |
|---|---|---|---|---|
| PB1 | Confronti di potenza aggiustata fra metodi (verso della differenza) nelle frasi di 3.8-3.9 e nelle Tabelle S10 e S23 che cambiano verso dopo la correzione della regola | 0 (cambiano 22 cifre stampate, stima del revisore) | 0 | un confronto che cambia verso cambia una frase dei risultati |
| PB2 | Potenza aggiustata del modello misto (test di Wald, approssimazione di Laplace) dovuta a fit con errore standard collassato: quota di dataset usabili sotto l'alternativa con \|stima/SE\| > 37 che superano la soglia, media sulle celle della griglia principale, intervallo bootstrap per blocchi `(model, i)`, B = 2000, seed 20261015 | circa 1,5 punti (sotto l'ipotesi nulla questi fit sono circa il 2,2% dei dataset usabili) | da 0 a 1 punto (margine aggregato) | dentro: trascurabile, nessuna frase; fuori: una frase in 3.8 con il valore (attesa da questa previsione) |
| PB3 | Celle della tabella dei parametri (`tables/V4_dgm_parameters.csv`) che cambiano rigenerandola da `R/v4/rundir` | 1 (D[3,3] di GLMM-AR1 M5 sotto H1, da 0,01 a 0,001) | esattamente 1 | ogni altra cella cambiata e' un errore della Tabella S1 da spiegare |
| PB4 | Cifre stampate che cambiano dopo la correzione dei denominatori del confronto diretto | 4 (celle della Tabella S31); quote di 3.4, pannello a della Figura 1 e confronto con due occasioni invariati alla cifra stampata | esattamente 4 | una cifra in piu' tocca una frase di 3.4 o la Figura 1 |
| PB5 | Applicazioni cliniche con errore standard del modello per entrambi i modelli: media sulle cinque dimensioni dei sottocampioni di (rapporto SE/dispersione delle stime del modello marginalizzato meno quello del modello misto), coorte e trial separati; intervallo bootstrap sui sottocampioni, B = 2000 | coorte +0,05; trial 0 | da -0,10 a +0,10 | dentro: il testo dice che i due modelli hanno errori standard ugualmente sotto la dispersione; fuori: il testo nomina la differenza |
| PB6 | Lunghezza del testo principale dopo correzioni e accorciamento (`tools/s7_wordcount.py`, 200 parole per tabella o figura) | circa 8.000 (7.824 della proposta piu' circa 200 di aggiunte) | da 7.500 a 8.200 | sopra 8.200 servono altri tagli o una lettera all'editor (decisione del PI) |
| PB7 | Secondo passaggio del revisore n.2: rilievi bloccanti (`MISMATCH` o `UNSUPPORTED` che cambiano abstract, 4.1 o conclusioni) e affermazioni di abstract e conclusioni senza una frase di sostegno nei risultati del testo principale | 0 e 0 | 0 e 0 | ogni rilievo bloccante torna nel testo prima della revisione di stile |
| PB8 | Requisiti di Morris con esito `Conforme` o `Dichiarato` dopo l'attuazione (ri-audit dei 5 non conformi o parziali di S8) | 28 su 28 | 28 su 28 | e' il criterio di invio della roadmap |

Margini. PB2 usa il margine aggregato di default per differenze fra versioni (1 punto); con circa 1.500 celle la MCSE della media e' sotto 0,05 punti, quindi l'esito non sara' indeterminato per precisione. PB5: la dispersione delle stime in ciascuna dimensione viene da circa 400-500 sottocampioni usabili, quindi l'errore del rapporto e' dell'ordine di 0,03-0,04 per dimensione e di circa 0,02 sulla media; un margine di 0,10 e' circa cinque volte quell'errore. PB6 segue il limite di SMMR (circa 8.000 parole) con la tolleranza che l'editor applica di solito; il valore e' un conteggio deterministico.

# Esiti (2026-11-01)

| Codice | Esito | Dettaglio |
|---|---|---|
| GB1 | superato | md5 di `results_full/` prima e dopo: 28.800 file, 0 differenze dalla lista di S7b (`reports/morris/S8b_logs/S8b_GB1_before.log`, `S8b_GB1_after.log`) |
| GB2 | superato | registro: 3.150 numeri da 53 file, 0 discordanze; GF2 0 cifre senza fonte; GF4 92 etichette nel testo principale e 41 nel supplemento, 0 non conformi, 0 frasi comparative senza etichetta |
| GB3 | superato | (a) 78 dei 79 id coincidono entro 1e-8, `x:cc_p11` differisce per costruzione (potenza grezza, testo corretto con P02); ricalcolo indipendente della Tabella S10 con differenza massima 5e-14; (b) 19,15, 38,58, 31,00, 69,55; (c) 0,352 |
| GB4 | superato | `reports/morris/S8b_logs/attuazione_S8b.csv`: 64 voci con esito; controllo delle ancore da script (26 sparite, 18 presenti con il testo nuovo verificato o motivo scritto, 20 senza ancora di testo) |
| GB5 | superato | abstract 200 parole; formato Word superato per manoscritto e supplemento; citazioni 1-33 in ordine di prima citazione, nessuna mancante o orfana; render docx, html e gfm senza errori |
| PB1 | confermata | 0 confronti su 84 cambiano verso; 22 cifre stampate cambiate |
| PB2 | smentita (atteso fuori) | 2,28 punti (MCSE 0,047; 2,19-2,37), fuori dall'intervallo 0-1: frase aggiunta (ora nella sezione S3 del supplemento, con la Tabella S10) |
| PB3 | smentita, spiegata | 9 righe cambiate, 1 cella stampata: la tabella vecchia era estratta prima delle correzioni [D-01](../decisions/D-01.md) e [D-02](../decisions/D-02.md); controprova 240 su 240 |
| PB4 | confermata | esattamente 4 cifre stampate cambiate |
| PB5 | confermata | coorte +0,043 (MCSE 0,003), trial +0,055 (MCSE 0,011), entrambi dentro ±0,10: il testo dice "by similar amounts". Il revisore n.2 chiedeva di nominare la differenza; ridimensionato, perche' la differenza e' dentro il margine pre-registrato |
| PB6 | confermata | 8.111 parole contate (6.911 + 6 x 200) alla chiusura, dopo le decisioni del PI sull'abstract e su 4.1; prima di quelle erano 8.180, dopo 91 parole tagliate o spostate nella revisione di stile per compensare le frasi di sostegno aggiunte |
| PB7 | smentita, corretta | 5 rilievi bloccanti (affermazioni di 4.1 e delle conclusioni senza frase di sostegno nei risultati: criterio precedente fino a 25 cluster, potenza sotto il 2% nei disegni piccoli, potenza del test sign-flip, intervallo invertito) e 0 sull'abstract; tutti corretti con frasi nelle sezioni 3.3, 3.4 e 3.7 prima della revisione di stile. Altri rilievi: 4 maggiori e 43 minori, tutti con un'azione nel log (`revisione_n2_S8b2_log.csv`) |
| PB8 | confermata | 28 su 28 (26 Conforme, 2 Dichiarato), colonna `esito_S8b` di `docs/audit/MORRIS2019_conformita.csv` |

Secondo passaggio del revisore n.2: due sotto-agenti, numeri e tabelle (413 verifiche, 400 PASS, 12 MISMATCH minori, 1 non verificabile) e testo e citazioni (115 verifiche, 70 PASS, 35 MISMATCH, 5 UNSUPPORTED, 5 non verificabili). Limiti dichiarati: bootstrap per blocchi non rieseguiti (controllati sui file di uscita); fit clinici del modello misto e delle GEE non rifatti. Bradley, Pepe e Anderson, Schildcrout e Heagerty, Morel, Ford e Westgate e Gosho, letti dal revisore solo nel riassunto o non letti, sono stati poi letti sul testo completo (PDF del PI in Zotero): confermano le citazioni; per Pepe e Anderson e Schildcrout e Heagerty l'esempio della frase di 4.5 e' cambiato da "covariate che rispondono a esiti passati" (non trattato: entrambi considerano covariate esogene) a "effetti ritardati di un'esposizione che varia nel tempo".

Revisione di stile: modifiche in `reports/morris/S8b_logs/style_S8b_changes.csv`, rilievi non modificati in `style_S8b_flags.md`.

# Dubbi rinviati

Nessuno. Le due scelte del PI sono registrate in [Abstract piu' prudente sul conteggio zero e numeri dell'intervallo sign-flip solo nei risultati (D-81)](../decisions/D-81.md); le celle che stampavano "-0.0" (Tabelle S6 e S30) ora stampano 0.0, e la frase della sezione 2.4 non attribuisce piu' al primo modo di fallimento di `nlm` una stima spostata, in linea con la correzione gia' fatta nel supplemento. I codici di altri registri non collegati restano come spiegato in [D-80](../decisions/D-80.md).
