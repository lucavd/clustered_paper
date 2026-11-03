---
type: Step
title: "Rafforzare il manoscritto v5 prima dell'invio: taratura, denominatori, modello misto, controllo di convergenza"
description: "Integrare nel testo del manoscritto v5, come parte del discorso generale, i punti fondati emersi da due letture critiche interne del 2026-10-30, con una sola analisi nuova (quadratura adattativa per il modello misto con sola intercetta casuale)."
step_id: "S7g"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-11-01T09:31:45+01:00, evidence: "ok del PI in chat sul riassunto In breve ('chiudi lo step'); commit S7g: approved by PI" }
resource: "../../reports/morris/S7g.html"
generated: { by: mmm-agent/skill-v4, at: 2026-10-31T17:50:00+01:00 }
sources:
  - id: rev_a
    resource: "artifact b551afef-5b75-4d10-acff-b75204db6987 (archivio delle sessioni)"
    title: "Lettura critica interna A del manoscritto v5 (2026-10-30)"
  - id: rev_b
    resource: "artifact 3ea46a71-c060-4a5f-9cae-27057b3728c0 (archivio delle sessioni)"
    title: "Lettura critica interna B del manoscritto v5 (2026-10-30, testo troncato in fondo)"
  - id: d61
    resource: "../decisions/D-61.md"
    title: "Figura 1 sulle strategie della regola; ordine della raccomandazione confermato (GEE prima, MMM alternativa)"
  - id: d63
    resource: "../decisions/D-63.md"
    title: "Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson"
  - id: d66
    resource: "../decisions/D-66.md"
    title: "La quadratura adattativa per il modello misto con sola intercetta casuale si valuta come analisi di sensibilita'"
---
# Apertura (2026-10-31)

1. **Lo step che ha scritto il manoscritto v5 e' chiuso e approvato.** S7f era stato committato alle 14:12 (b16c0b6) con `pi_approval: pending`. Il PI ha letto il riassunto e l'ha approvato all'apertura di questa sessione (commit 99b11ca). Nello stesso commit sono entrati i file non tracciati rimasti: log piccoli di S7d e S7e, script di verifica di S7f. Erano datati 23-27 ottobre e sul server non girava nessun processo R; il PI ha confermato che nessun'altra chat era attiva. I due log di run (231 e 63 MB) restano sul server, esclusi da git, con una ricevuta md5 in `reports/morris/large_logs_receipt.txt`.
2. **Testo di partenza invariato.** `manuscript/v5/MMM_v5.qmd` ha md5 `9a7cf1b25d20defef5c628c7603fb4f0` e il supplemento `d3f001213a8f111b840655e9136421f9`, come nel brief: i numeri di riga sotto valgono.
3. **Dove si modifica.** Il manoscritto si modifica sul posto ([Il manoscritto v5 si modifica sul posto nello step S7g; non si crea una v6 (D-67)](../decisions/D-67.md)).
4. **Due correzioni al brief, decise dal PI all'apertura.**
   - Vale la regola dell'intervallo esatto di Clopper-Pearson per "chiaramente sopra 7,5%" e "chiaramente sotto 2,5%" ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)). Il brief la contraddiceva in tre punti (punto 1, "Non si fa", appendice), corretti sotto e segnati come tali.
   - La quadratura adattativa si fa, benche' [Nessuna nuova simulazione dopo la revisione avversariale; i metodi non valutati si dichiarano fra i limiti (D-64)](../decisions/D-64.md) l'avesse esclusa: la nuova decisione la supera solo per questo punto ([La quadratura adattativa per il modello misto con sola intercetta casuale si valuta come analisi di sensibilita' (supera in parte la decisione del 26 ottobre) (D-66)](../decisions/D-66.md)).

# Regola di scrittura (indicazione del PI, 2026-10-31)

Il paper non è ancora stato inviato. **Tutto entra nel discorso generale del manoscritto, non come risposta a qualcuno.**

- Nessuna formula da lettera di risposta: niente "as suggested", "in response", "we now report", "additionally", "to address concerns", né sezioni o tabelle presentate come "richieste".
- Ogni aggiunta sta nel punto in cui un lettore la cercherebbe (metodi, risultati, limiti), con lo stesso tono del resto del testo e senza tono difensivo.
- **Integrare non vuol dire nascondere.** Ciò che è stato deciso dopo aver visto dei risultati resta dichiarato come tale: il controllo di convergenza, la regola che sostituisce il criterio della versione precedente, il confronto fra le due strategie finaliste, e la nuova analisi con quadratura adattativa, che si presenta come analisi di sensibilità aggiunta dopo il run principale.
- L'origine dei punti (due letture critiche interne) si registra solo in questo file e nelle decisioni, mai nel manoscritto né nel supplemento.
- Vincoli invariati: numeri solo dal registro (`n("...")`), frasi comparative con etichetta, abstract non strutturato di al massimo 200 parole, limite SMMR di circa 8.000 parole con tabelle e figure a 200 parole.

# Dubbi ricevuti

Ricevuti il 2026-10-31 dalla verifica delle due letture critiche interne contro `MMM_v5_resolved.md` (2026-10-26) e il docx v5. I numeri di riga si riferiscono a `MMM_v5.qmd` con l'md5 sopra.

## Interventi decisi dal PI (2026-10-31: "mi va bene cosa hai deciso")

1. **Cosa garantisce davvero la regola sui falsi positivi.**
   - *Perché.* "Never clearly above 7.5%" significa nessun eccesso rilevato, non taratura dimostrata. Con il criterio dell'intervallo esatto di Clopper-Pearson ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)), con 1.000 repliche usabili un punto viene segnalato solo se il tasso osservato e' almeno 9,3% (93 rifiuti); con 500 repliche, almeno 10,0% (50 rifiuti); sotto le 20 repliche usabili il punto non viene valutato. *(Corretto all'apertura: il brief dava 9,2% e 9,9%, calcolati con tasso - 2 MCSE di Wald, la regola sostituita da [Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md).)*
   - *Cosa.* Una tabella nel supplemento, per ogni strategia della Tabella 4, con:
     - tasso massimo osservato, con le sue repliche usabili e la MCSE;
     - punti chiaramente sopra 7,5% (estremo inferiore dell'intervallo esatto oltre 7,5%) e punti in cui un eccesso non e' escluso (estremo superiore oltre 7,5%, ma non chiaramente sopra);
     - punti non valutati (meno di 20 repliche usabili);
     - per Fay-Graubard, la differenza massima e la distribuzione per punto fra le tre correlazioni di lavoro (oggi il testo dà solo le medie, −0,29 e −0,31 punti).
   - *Testo.* Una frase nei metodi su cosa significa "clearly above" in pratica (riga 97, sezione 2.6); una frase nei risultati che rimanda alla tabella (riga 133, sezione 3.4); una frase sulla dispersione fra correlazioni di lavoro (sezione 3.5).
   - *Fonti.* `tables/V5_S7f_cp_S7d_dp_rates.csv`, `tables/V5_S7f_cp_S7e_dp_rates.csv` e `tables/V5_S7f_cp_S7b_rule_dp.csv` (colonne `n`, `nd`, `rate`, `mc`, `assess`, `hi`, `lo`, `pu`; flag gia' secondo [Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)). *(Corretto all'apertura: il brief citava le tabelle V4 con i flag di Wald.)*
2. **Un solo denominatore per "risultato usabile".**
   - *Il problema.* La sezione 2.5 (riga 89) dice "every attempted replicate in the denominator"; la didascalia della Tabella 4 (riga 283) dice "denominator the replicates with an estimable design". Le repliche con un braccio mai osservato dopo la prima occasione (247 per braccio) sono escluse da ogni denominatore.
   - *Cosa.*
     - Verificare nel codice quale denominatore usano Tabella 3, Tabella 4 e Tabella [manoscritto (S7)](S7.md), e allineare il testo al codice.
     - Nelle due applicazioni cliniche (riga 165 e paragrafo sul trial respiratorio; Tabelle S15 e S30), dare la quota usabile su tutti i sottocampioni, compresi quelli in cui il modello della media non è identificabile. Oggi il testo dà "Among the fits that were returned" (per esempio Fay-Graubard 45%, modello marginalizzato 100% nella coorte).
3. **Selezione sui risultati usabili.**
   - *Cosa.* Rendere visibili nel testo vicino alla sezione 2.5 le tabelle che già esistono: la Tabella [tabelle e figure (S6)](S6.md) (errore di I tipo sulle repliche usabili del metodo e su quelle comuni a tutti i metodi) e la Tabella S9 (per strato di stimabilità).
   - *Calcolo nuovo dalle righe esistenti.* Sensibilità della Tabella 4 ai limiti del filtro (|stima| < 10, SE < 10), con limiti 5 e 20. Una frase con l'esito.
4. **La conclusione sul modello misto vale per il test di Wald con approssimazione di Laplace.**
   - *Testo.* Abstract (riga 28: "The mixed model was anticonservative when the dependence was serial"), risultati sul quarto meccanismo, discussione (righe 185 e 193) e conclusioni (righe 197-199).
   - *Nel supplemento, sezione [revisione finale (S8)](S8.md), se non c'è già.* Su 20.000 soggetti sotto l'ipotesi nulla, nel meccanismo con copula il modello misto converge a 0,005 (SE 0,005) (`tables/V4_S7e_glmm_pseudotrue.csv`, step del quarto meccanismo S7e). L'eccesso quindi non dipende da un'ipotesi nulla del modello misto diversa da zero. Da rileggere dal file prima di scriverlo.
5. **Analisi nuova: quadratura adattativa per il modello misto con la sola intercetta casuale** ([La quadratura adattativa per il modello misto con sola intercetta casuale si valuta come analisi di sensibilita' (supera in parte la decisione del 26 ottobre) (D-66)](../decisions/D-66.md)). Vedi "Attese pre-registrate" sotto. L'esito decide la forma della frase del punto 4: rafforzata se l'anticonservativismo resta, riscritta se sparisce.
6. **Il controllo di convergenza dove il modello sostitutivo è sbagliato.**
   - *Cosa.* Negli scenari M6-M7 riportare bias/SE empirico, rapporto fra SE del modello e SE empirico, e copertura del modello marginalizzato con controllo.
   - *Fonti.* Parte c'è già in `tables/V4_S6b_ns_dp.csv` (colonne `bias`, `empse`, `modse`); la copertura va calcolata dalle righe se manca.
   - *Descrizione aggiuntiva.* Stessa sintesi separata per fit in cui il controllo è intervenuto o no, presentata come descrizione. Non è un criterio, perché l'intervento dipende dai dati.
   - *Testo.* Una frase nei risultati (riga 121, sezione 3.2) e colonne nella Tabella S14.
   - *Limiti.* Trend della deviazione standard latente fra 0,01 e 0,1 per occasione non simulato, con il rimando al consiglio, già nel testo (riga 181), di riportare quante volte il controllo interviene.
7. **Raccomandazione simmetrica, nell'ordine di [Figura 1 sulle strategie della regola; la GEE AR(1) con Fay e Graubard e t va prima, il modello marginalizzato con controllo e' un'alternativa comparabile (D-61)](../decisions/D-61.md).**
   - *Formula.* GEE AR(1) con Fay-Graubard e t preferibile nel dominio studiato; modello marginalizzato con controllo e t come alternativa competitiva. Abstract (riga 28), discussione 4.1, conclusioni (riga 199).
   - *Numeri per la frase.* Già decisi in [Figura 1 sulle strategie della regola; ordine della raccomandazione confermato (GEE prima, MMM alternativa) (D-61)](../decisions/D-61.md):
     - distanza media dal 5%: 1,06 contro 1,07;
     - GEE più vicina con T ≥ 3, soprattutto a 10 cluster;
     - modello marginalizzato più vicino con T = 2 e N ≥ 15;
     - usabilità GEE +1,23 punti.
   - Controllare che l'abstract non dica più di questo.
8. **Coorte neonatale: due numeri da spiegare.**
   - *Errore standard del modello misto.* È 0,51-0,58 volte la dispersione empirica. Ipotesi da verificare: il modello misto ha solo l'intercetta per paziente e ignora il livello paziente-occasione (varianza 3,8 nel modello completo). Calcolo leggero sugli stessi sottocampioni. Esito in una frase alla riga 165.
   - *Disegni non identificabili.* Attesi 26,8% (ipergeometrica: 182 su 207 pazienti senza ventilazione invasiva, 10 estratti), osservati 22% su 500 sottocampioni, circa 2,4 errori standard sotto l'atteso. Rileggere il valore esatto e la definizione da Tabella S15, quindi spiegare o correggere.
9. **Frasi da integrare (solo testo):**
   - **Osservazione non informativa** (limiti, riga 193): il numero di occasioni osservate non dipende dagli esiti (completo con probabilità 0,7, altrimenti uniforme da 1 a T−1); dropout informativo non studiato. I risultati riguardano un'interazione lineare trattamento per tempo.
   - **La regola è stata fissata dopo i risultati della versione precedente.** Verificare che i limiti lo dicano, oltre al controllo e al confronto fra le due strategie (la frase "replaces the criterion of an earlier version" è nei metodi).
   - **Potenza aggiustata per il livello** (riga 91): il valore critico è stimato sotto l'ipotesi nulla di ogni punto; serve a confrontare i metodi, non è una procedura applicabile a dati nuovi.
   - **Condizioni di validità del test a inversione di segno** (riga 85 o discussione): indipendenza fra cluster e simmetria dei contributi allo score sotto l'ipotesi nulla; eterogeneità fra cluster ammessa.
   - **Perché Bradley come criterio eliminatorio** (riga 91 o 97): l'analista conosce N e T ma non il meccanismo. Bradley liberale = da 0,5α a 1,5α, cioè dal 2,5% al 7,5%.
   - **Kauermann-Carroll** (riga 79 e sezione 4.3, riga 185): portare nel testo principale che `fewclust` concorda con il calcolo a mano ([riaggregazione con MCSE (S4)](S4.md)) e che i valori di letteratura calcolati con `geesmv` non sono direttamente confrontabili.
   - **Mole del calcolo** (sezione 2.2, riga 56): numero totale di dataset e di fit, e tempo complessivo, dai log dei run.
   - **Trial respiratorio** (riga 103): una frase in più sul perché la variabile latente ha solo l'intercetta con quattro visite.

## Non si fa (deciso 2026-10-31)

- **Simulare un trend intermedio della variabilità latente (0,02-0,05 per occasione).** Basta la frase di limite al punto 6; lo si fa solo se lo chiede un revisore della rivista.
- **Riscrivere Bradley come "±1 punto".** E' sbagliato: il criterio e' quello liberale di Bradley (da 2,5% a 7,5%) giudicato per punto con l'intervallo esatto di Clopper-Pearson ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)). *(Corretto all'apertura: il brief diceva che il paper usa ±2 MCSE e scartava Clopper-Pearson; vale [Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md), confermata dal PI il 2026-10-31.)*
- **La "lacuna" nella sezione 2.2.** È un artefatto di estrazione: il sorgente e il docx contengono "1 to T−1" come formula.
- **Riformulare la raccomandazione su exchangeable "se il pattern è simmetrico".** I dati non la sostengono: exchangeable con Fay-Graubard e' chiaramente sopra 7,5% in 11 punti (7 con la regola di Wald, prima di [Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)).

# Attese pre-registrate (2026-10-31, prima del run con quadratura e prima di modificare il testo)

## Cosa si fa

1. **Run con quadratura** (dubbio 5). Script `R/s7g/s7g_common.R`, `s7g_run.R`. Il modello misto dello studio (glmer, intercetta casuale per soggetto; formula `binY ~ time*trt (+ age in M2) + (1 | id)`, famiglia binomiale, `glmerControl` di default come nello studio) si rifitta sugli stessi dataset dello studio, rigenerati con i semi `20260805 + 1000*model + i` e L'Ecuyer-CMRG (M1c usa i semi di M1). Si fitta tre volte: `nAGQ = 1` (Laplace, deve riprodurre le righe salvate), `nAGQ = 10` e, per le repliche 1-100, `nAGQ = 25`. Cambia solo `nAGQ`.
   - Punti di disegno: M1 e M2 nei quattro meccanismi (dati generati sulla scala della popolazione, dati generati dal modello misto, lo stesso con correlazione AR(1) latente, copula con correlazione seriale) con N in {10, 25, 55, 145} e T in {2, 4, 10}; M1c nei tre meccanismi dello studio (non esiste per la copula) con N in {10, 25, 45} e T in {2, 4, 6}, perche' la griglia di M1c arriva a 45 soggetti e 6 occasioni. Totale 123 punti, 1.000 repliche, solo senza effetto del trattamento.
   - Prova di funzionamento (2 repliche x 5 punti, niente salvato): 0,39 s per fit con Laplace, 0,59 con 10 punti, 0,76 con 25 nel punto piu' grande. Stima 12-15 ore-core, meno di un'ora con 20-25 processi sotto la guardia termica (tetto 28).
2. **Analisi.** Script `R/s7g/s7g_agg.R`, scritto prima del run.
   - Risultato usabile: la definizione dello studio (`u0_flag`: stima e SE finiti, |stima| < 10, 0 < SE < 10, p finito; fit con avvisi di convergenza tenuti), disegni non stimabili esclusi da ogni denominatore come nello studio.
   - Errore di I tipo: Wald con riferimento t a N - p gradi di liberta' (p = 4 in M1 e M1c, 5 in M2), primario come nella Tabella 4; normale come descrizione. Tasso su usabili del fit; flag con l'intervallo esatto di Clopper-Pearson ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)), punto valutato con almeno 20 usabili.
   - Differenza quadratura (10 punti) meno Laplace per punto: sulle repliche usabili per entrambi i fit, appaiata; MCSE dalla deviazione standard delle differenze degli indicatori di rifiuto.
   - Aggregato per meccanismo: media sui punti della differenza; intervallo con bootstrap per blocchi di replicato `(model, i)`, B = 2000, seme 20261015 (regola 6), percentile 2,5-97,5.
3. **Testo** (dubbi 1-4 e 6-9). Si integra nel manoscritto e nel supplemento v5 con la catena di S7f (`R/s7f/s7f_pipeline.sh` e i suoi controlli). Ogni numero nuovo entra nel registro con il suo script.

## Cancelli

| Codice | Cancello | Come si controlla |
|---|---|---|
| GG1 | Integrita' e parita' | md5 di `results_full/` prima e dopo uguali alla lista di S7b (0 differenze); `R/testing/test_h0h1_parity.R` passa prima del run |
| GG2 | Riproduzione | per tutte le repliche dei 123 punti, stima e SE di `nAGQ = 1` uguali alle righe salvate (`results_v4/variant_cc/S4/rows_v4.rds` per i tre meccanismi, `results_v4/S7e/rows` per la copula) entro 1e-8; stesso insieme di repliche con o senza riga |
| GG3 | Completezza | 1.230 file di chunk, 0 falliti, 123.000 repliche con i tre o due fit registrati |
| GG4 | Stabilita' della quadratura | sulle repliche 1-100, fra 10 e 25 punti la decisione di rifiuto (t, N - p) coincide in almeno il 99,5% dei fit usabili per entrambi; se fallisce, i 10 punti non bastano e l'analisi si rifa' con 25 |
| GG5 | Numeri del manoscritto | catena di S7f: registro senza discordanze, cifre scritte a mano controllate, etichette sulle frasi comparative, abstract entro 200 parole, formato Word, `tools/okf.py check` |

## Previsioni

Margini di default del protocollo (differenza di errore di I tipo fra versioni o metodi: +/- 1 punto in media sui punti, +/- 2,5 punti per punto). Stima +/- 2 MCSE (per l'aggregato, intervallo del bootstrap per blocchi).

| Codice | Quantita' | Atteso | Intervallo di rilevanza | Perche' conta |
|---|---|---|---|---|
| PG1 | Differenza di errore di I tipo (t), quadratura a 10 punti meno Laplace, media sui punti, per ciascuno dei quattro meccanismi | 0 | da -1 a +1 punto | se fuori, la conclusione sul modello misto dipende dall'approssimazione e la frase dell'abstract va riscritta |
| PG2 | La stessa differenza punto per punto (123 punti) | 0 in tutti | da -2,5 a +2,5 punti | conteggio dei punti per esito; i punti smentiti si spiegano. Il punto piu' a rischio sono le 2 occasioni, dove Laplace e' meno accurata |
| PG3 | Fra i punti della copula chiaramente sopra 7,5% con Laplace, quota che resta chiaramente sopra con la quadratura | 100% | almeno 80% | e' la base della frase "anticonservativo quando la dipendenza e' seriale"; sotto 80% la frase cambia. E' un conteggio, senza intervallo |
| PG4 | Sensibilita' della Tabella 4 ai limiti del filtro (|stima| e SE sotto 5 e sotto 20 invece di 10) | stesse due strategie mai chiaramente sopra 7,5%, stessa scelta della regola | nessun cambiamento | se la scelta cambia, la raccomandazione dipende da un parametro arbitrario |
| PG5 | Coorte neonatale: rapporto fra SE del modello misto e dispersione empirica nei sottocampioni, con un effetto casuale in piu' per paziente-occasione | da 0,85 a 1,15 | da 0,85 a 1,15 | se resta vicino a 0,5-0,6 l'ipotesi del brief cade e la frase del testo dice che il motivo non e' identificato |

Esito della quadratura e testo (dubbio 4): se PG1 e PG3 sono confermate, la frase sul modello misto si rafforza (l'eccesso resta con la quadratura); se PG1 e' smentita nel meccanismo con copula o PG3 e' smentita, la frase si riscrive e il PI decide come.

# Cosa si chiude in questo step

Manoscritto e supplemento v5 aggiornati e resi sul server, con le stesse verifiche di S7f:

- numeri attestati dal registro;
- controllo delle cifre scritte a mano;
- etichette sulle frasi comparative;
- conteggio delle parole;
- abstract entro 200 parole;
- formato Word;
- `tools/okf.py check`.

Poi report, riassunto «In breve» e messaggio al PI. La revisione avversariale finale e il ri-audit Morris restano allo step successivo.

# Appendice: verifica delle due letture critiche (2026-10-30/31)

Fatta sul testo, non su codice o risultati grezzi.

| Punto | Lettura | Esito | Nota |
|---|---|---|---|
| Assenza di eccesso ≠ taratura dimostrata | A | regge | la regola descritta come "intervallo esatto": in realtà ±2 MCSE |
| Selezione sui risultati usabili, punti non valutati | A | regge | trovata in più l'incoerenza di denominatore (punto 2) |
| Controllo di convergenza come procedura a due stadi | A, B | regge in parte | frequenze, soglie e natura post hoc già nel testo; manca bias/copertura in M6-M7 |
| GEE contro modello marginalizzato, regola discreta | A | regge nella sostanza | conteggi "sotto 2,5%" sbagliati (116/24; veri 136/86, a 10 cluster 18/60) |
| Ipotesi nulle marginale e condizionale | A | già quasi coperto | quantificato nel testo; pseudo-vero del modello misto sotto copula vicino a 0 (S7e) |
| Modello misto: Laplace | A | regge | |
| Osservazione non informativa | A | regge | |
| Applicazioni cliniche non sono validazione; denominatori | A | regge | |
| Scelta AR(1) poggia su pochi punti | B | regge nella sostanza | numeri sbagliati (11 invece di 7; "sopravvivono" invece di "escluse") |
| Trend intermedio non simulato | B | regge, rischio basso | "ottimizzato sui parametri" non regge |
| Errore standard del modello misto 0,51-0,58 nella coorte | B | regge | sfuggito ad A |
| Lacuna in 2.2 | B | non regge | artefatto di estrazione |
| Bradley "±1 punto", Clopper-Pearson | B | non regge | |
| "Regola pre-specificata" (fra i punti di forza) | B | non regge | la regola sostituisce il criterio della versione precedente |
| 138 di 864 riportato come 149 | B | errore di citazione | |

**Nota all'apertura (2026-10-31).** Nelle righe "GEE contro modello marginalizzato" e "Scelta AR(1)" i valori indicati come veri (136 e 86 punti chiaramente sotto 2,5%; 7 punti sopra per la GEE exchangeable) sono quelli della regola di Wald. Con la regola in vigore ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)) i conteggi sono 116, 24 e 11. Il manoscritto v5 li prende dal registro dei numeri; il cancello sui numeri li ricontrolla.

# Esito delle attese (2026-11-01)

## Cancelli

| Codice | Esito | Dettaglio |
|---|---|---|
| GG1 | superato | md5 di `results_full/` prima e dopo uguali alla lista di S7b (28.800 file, 0 differenze); test di parita' superato |
| GG2 | superato | 122.904 righe confrontate con le righe salvate, 0 differenze oltre 1e-8 (`tables/V5_S7g_agq_gates.csv`) |
| GG3 | fallito come scritto, spiegato | 1.230 chunk, 0 falliti, 123.000 repliche; 96 repliche di GLMM-AR1 non hanno fit perche' il dataset non si genera (calibrazione singolare): sono esattamente le 96 non generate negli step dello studio, con gli stessi flag di disegno degenere. Il cancello come scritto chiedeva i fit su tutte le repliche, cosa che il disegno dello studio esclude in partenza |
| GG4 | superato | 9.844 fit usabili con 10 e 25 punti: decisione diversa in 10 (99,90% uguali), mediana di abs(z25 - z10) 2,3e-05 |
| GG5 | superato | registro 2.803 numeri da 51 file, 0 discordanze; 0 cifre a mano fuori classe; 114 etichette, 0 frasi comparative senza etichetta; abstract 199 parole; formato Word conforme; `tools/okf.py check` superato |

## Previsioni

| Codice | Atteso | Osservato | Esito |
|---|---|---|---|
| PG1 | 0, intervallo +/- 1 punto | GenBinaryY -3,4 [-3,7; -3,2]; GLMM-native -3,3 [-3,6; -3,0]; GLMM-AR1 -3,1 [-3,6; -2,7]; Copula-AR1 -4,0 [-4,6; -3,5] | smentita in tutti e quattro i meccanismi |
| PG2 | 0 in tutti i punti, +/- 2,5 | 27 confermati, 63 indeterminati, 33 smentiti su 123; i piu' grandi a T = 2 con N = 145 (fino a -17,8 punti) | smentita in 33 punti |
| PG3 | 100%, almeno 80% | 6 dei 12 punti della copula chiaramente sopra con Laplace restano chiaramente sopra (50%), tutti a T = 10 | smentita |
| PG4 | stessa scelta con limiti 5 e 20 | stesse due strategie mai chiaramente sopra 7,5%, stessa scelta; con limite 10 la Tabella 4 si riproduce per 17 strategie su 17 | confermata |
| PG5 | rapporto da 0,85 a 1,15 | 0,50-0,57 con l'intercetta per paziente-visita (varianza al bordo fino a 25 pazienti) | smentita |

Spiegazioni. PG1-PG3: l'approssimazione di Laplace e' imprecisa con pochi esiti binari per soggetto. Con la quadratura il test del modello misto sui suoi stessi dati non e' piu' chiaramente sopra 7,5% in nessun punto, e diventa chiaramente sotto 2,5% in 37 punti (3 con Laplace). L'eccesso resta solo a T = 10: 6 punti con la copula, 2 con GenBinaryY. Il PI ha deciso come riscrivere ([La conclusione sul test di Wald del modello misto si riscrive con l'esito della quadratura adattativa, senza altri calcoli (D-68)](../decisions/D-68.md)). PG5: l'ipotesi del brief cade. Una spiegazione cercata dopo (pendenza casuale per le zone posteriori, varianza 10,2 nella coorte completa) porta il rapporto a 0,66-0,73 a 25 e 50 pazienti: spiega una parte dello scarto, non tutto (`tables/V5_S7g_lus_glmm_slope.csv`). Le code pesanti non c'entrano: con una dispersione robusta il rapporto resta 0,57-0,62.

# Esito dei dubbi ricevuti

1. **Cosa garantisce la regola**: corretto. Frase nei metodi (soglia 9,3% con 1.000 usabili, 10% con 500), frase nei risultati 3.4, Tabella S32 (per strategia: chiaramente sopra, eccesso non escluso, non valutati, tasso massimo); dispersione di Fay-Graubard fra correlazioni di lavoro nella sezione 3.5 (90% delle differenze fra -1,4 e +0,4 punti, 10 punti su 1.896 oltre 2,5, massimo 4,2).
2. **Denominatori**: spiegato e corretto. Nel codice le Tabelle 3, 4 e [manoscritto (S7)](S7.md) usano lo stesso denominatore (repliche con disegno stimabile; identico su 1.896 punti fra le due catene). Cambiate le parole della sezione 2.5 e delle didascalie delle Tabelle 3, [manoscritto (S7)](S7.md) e S20. Coorte: quota usabile su tutti i 500 sottocampioni a 10 pazienti nel testo.
3. **Selezione sui risultati usabili**: corretto. Rimando alle Tabelle [tabelle e figure (S6)](S6.md) e S9 nella sezione 2.5; sensibilita' ai limiti del filtro (Tabella S33, PG4 confermata).
4. **Modello misto e Laplace**: corretto secondo [La conclusione sul test di Wald del modello misto si riscrive con l'esito della quadratura adattativa, senza altri calcoli (D-68)](../decisions/D-68.md). Il valore a cui converge il modello misto sulla copula, riletto dal file, va da 0,0046 a 0,060 con SE fino a 0,115 (il brief citava solo 0,005 a T = 10; a T = 4 e' 0,036 con SE 0,024); nel testo e in [revisione finale (S8)](S8.md) come intervallo di valori.
5. **Quadratura adattativa**: fatta; paragrafo 3.10 e Tabella S34.
6. **Controllo di convergenza in M6-M7**: corretto. Frase nella sezione 3.2 e Tabella S35 (nuova invece di colonne della Tabella S14: decisione tecnica sotto). Distorsione al massimo 0,40 errori standard; SE del modello in mediana dal 2,7 al 4,7% sotto quello empirico; copertura t media 94,8-95,6%, non copertura mai chiaramente sopra 7,5%; con e senza intervento 94,7-95,7%. Trend non simulato nei limiti.
7. **Raccomandazione simmetrica**: corretto. Conclusioni: GEE raccomandata "in the designs studied", modello marginalizzato "competitive alternative, slightly closer to nominal with two occasions and 15 or more clusters"; l'abstract non dice di piu'.
8. **Coorte neonatale**: SE del modello misto spiegato in parte (sopra). Disegni non identificabili spiegati: estrazioni casuali semplici (`rev2/T2_defs.R`), 109 su 500 contro 133,8 attesi, z = -2,5; agli altri tagli abs(z) al massimo 1,1; numero medio di ventilati 1,28 contro 1,21 atteso.
9. **Frasi da integrare**: tutte nel testo (osservazione non informativa e interazione lineare nei limiti; regola fissata dopo i risultati della versione precedente nei metodi e nei limiti; potenza aggiustata solo per confrontare i metodi; condizioni del test a inversione di segno; Bradley da 0,5 a 1,5 volte alpha; Kauermann-Carroll da fewclust e valori calcolati con geesmv non confrontabili; trial respiratorio). Mole del calcolo: nel testo il numero di dataset (4.752.000); tempo complessivo e numero di fit rinviati (sotto).

# Dubbi emersi durante lo step

1. **Previsioni sulla quadratura smentite**: deciso dal PI ([La conclusione sul test di Wald del modello misto si riscrive con l'esito della quadratura adattativa, senza altri calcoli (D-68)](../decisions/D-68.md)), testo riscritto.
2. **PG5 smentita**: spiegato in parte; il testo lo presenta come spiegazione cercata dopo.
3. **Pseudo-vero della copula diverso dal brief**: corretto nel testo (intervallo dei valori, non solo T = 10).
4. **Lunghezza**: 9.911 parole con 200 per tabella o figura (8.858 alla chiusura di S7f), contro circa 8.000 della rivista; gia' rinviata alla revisione finale da S7f, numero aggiornato nel file di [revisione finale (S8)](S8.md).
5. **Tempo di calcolo complessivo**: rinviato (sotto).
6. **Due arrotondamenti fuorvianti nel render** ("100%" per 99,90% di decisioni uguali; "95% to 95%"): corretto, frasi riscritte con conteggi o senza numeri.

# Dubbi rinviati

1. **Tempo complessivo di calcolo e numero di fit** (dubbio 9) -> [Revisione avversariale e ri-audit Morris (S8)](S8.md). Perche': i log dei run della versione 4 danno la durata solo dell'ultimo segmento (S6b ha un primo segmento senza riga DONE) e il run originale di `results_full/` non ha log nel repository; un numero senza fonte non entra. Per chiuderlo: ricostruire le durate dai tempi dei file di chunk, oppure lasciare nel testo solo il numero di dataset.

# Decisioni tecniche dell'agente

- I due log di run di S7d e S7e (231 e 63 MB) restano sul server, esclusi da git, con ricevuta md5 (`reports/morris/large_logs_receipt.txt`).
- Le tabelle nuove del supplemento vanno in coda (S32-S35, sezioni S11-S13) invece di colonne della Tabella S14, per non rinumerare le tabelle citate dal testo.
- La catena di S7f si riusa in copia (`R/s7g/s7g_pipeline.sh`, `py/s7g_build_spec.py`, `tools/s7g_check_*.py`) con i log in `reports/morris/S7g_logs/`, cosi' i log di S7f restano intatti; il registro resta `tables/V5_S7f_numbers.csv`, letto dal manoscritto.
- Le frasi aggiunte non usano formule da lettera di risposta; l'origine dei punti (letture critiche) non compare nel manoscritto.

# In breve

## In breve: il manoscritto regge ai punti deboli segnalati dalle due letture critiche?

**A cosa serviva.** Due letture critiche interne del manoscritto avevano segnalato alcuni punti deboli: cosa garantisce davvero la regola che sceglie la raccomandazione, come si contano i risultati usabili, il giudizio sul modello misto, e il controllo di convergenza del modello marginalizzato multilivello quando il modello che lo sostituisce è sbagliato. Li abbiamo portati nel testo come parte del discorso, con un solo calcolo nuovo: il modello misto stimato con la quadratura adattativa invece che con l'approssimazione di Laplace. Conta perché dopo questo passo il paper va alla revisione finale e all'invio.

**Cosa abbiamo trovato.**

1. **Buona parte dei falsi positivi in eccesso del modello misto viene dall'approssimazione di Laplace, non dal modello.**
   - Il modello misto qui è il modello logistico a effetti casuali per soggetto, con la sola intercetta. L'abbiamo ristimato con la quadratura adattativa sugli stessi dataset, in 123 situazioni senza effetto del trattamento.
   - In media i falsi positivi scendono di 3-4 su 100. Con 2 visite e 145 soggetti passano da 10-20 su 100 a 3-5.
   - Restano chiaramente sopra 7,5 su 100 solo con 10 visite: in 6 situazioni con i dati in cui la dipendenza è solo fra visite vicine (12-16 su 100) e in 2 con i dati generati sulla scala della popolazione (circa 9-10).
   - Con pochi soggetti il test diventa prudente: in 37 situazioni fa chiaramente meno di 2,5 falsi positivi su 100, contro 3 con Laplace.
   - Ci aspettavamo che non cambiasse nulla. Come hai deciso, abstract, risultati, discussione e conclusioni ora lo dicono; la raccomandazione non cambia.
2. **"Mai chiaramente sopra 7,5 falsi positivi su 100" vuol dire che non abbiamo visto un eccesso, non che il metodo è tarato: ora il paper lo scrive.**
   - Con 1.000 simulazioni utili una situazione viene segnalata solo da 9,3 falsi positivi su 100 in su.
   - Per la GEE raccomandata il valore più alto osservato è 8,6 su 100, e in 218 delle 2.376 situazioni un eccesso non si può escludere; per il modello marginalizzato le situazioni sono 438.
   - Se si sposta la soglia oltre la quale una stima o il suo errore standard non si considera usabile (5 o 20 invece di 10), la regola sceglie ancora le stesse due strategie, nello stesso ordine.
3. **I risultati usabili si contavano già allo stesso modo in tutte le tabelle: abbiamo corretto solo le parole.**
   - Il denominatore sono sempre le simulazioni con un disegno stimabile, comprese quelle in cui il modello non si riesce a stimare.
   - Nella coorte neonatale ora diamo anche la quota usabile su tutti i sottocampioni da 10 pazienti: modello misto 78%, modello marginalizzato 77%, GEE senza correzione 70%, GEE con la correzione di Fay e Graubard 34%, con quella di Mancl e DeRouen 15%.
4. **Il controllo di convergenza funziona anche dove il modello che sostituisce è sbagliato.**
   - Negli scenari in cui la variabilità latente cresce molto nel tempo, il modello marginalizzato con il controllo stima senza distorsione rilevante, il suo errore standard è più basso del vero di appena il 3-5%, e l'intervallo copre il valore vero in circa 95 casi su 100.
   - Vale sia quando il controllo interviene sia quando non interviene.
5. **Nella coorte neonatale l'errore standard troppo piccolo del modello misto si spiega solo in parte.**
   - L'ipotesi che mancasse il livello paziente-visita non regge: aggiungerlo non cambia nulla.
   - Regge in parte un'altra spiegazione: la differenza fra zone posteriori e anteriori varia molto da paziente a paziente, e il modello con la sola intercetta non lo vede. Tenendone conto, l'errore standard passa da poco più di metà a circa due terzi della dispersione vera delle stime, non a tutta. Il testo lo dice così.
   - I sottocampioni senza neonati ventilati sono il 22%, contro il 27% atteso: le estrazioni sono casuali e la differenza è del caso.

**Cosa devi decidere.** Nulla: come riscrivere la frase sul modello misto l'hai già deciso durante lo step.

**Comunque:**

- Il testo ora conta circa 9.900 parole, contando 200 parole per ogni tabella o figura, contro le circa 8.000 della rivista (alla fine del passo precedente erano circa 8.900). Lo accorceremo nella revisione finale, come già previsto.
- Del calcolo complessivo il paper dà il numero di dataset (4.752.000). Il tempo totale e il numero di modelli stimati non si ricostruiscono dai registri dei run: ci torneremo nella revisione finale.
- Le tabelle principali restano con l'approssimazione di Laplace; la quadratura sta in un paragrafo dei risultati e in una tabella del supplemento.

Non c'è nulla da decidere: se non hai dubbi, lo step è chiuso.

**Sui tempi di calcolo.** La stima per il calcolo con la quadratura l'ha fatta l'agente: 12-15 ore-core, dai tempi misurati su una prova di 10 repliche (0,39 secondi per stima con Laplace, 0,59 con la quadratura nel caso più grande). Il calcolo è durato 0,92 ore con 8-14 processi sul server; il brief prevedeva 0,5-1,5 ore.
