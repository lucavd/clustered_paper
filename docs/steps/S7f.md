---
type: Step
title: "Aggiornamento del manoscritto con le analisi aggiunte"
description: "Riportare nel manoscritto i risultati delle rianalisi, del pilota, del run con i metodi aggiunti e del terzo meccanismo, riscrivere le limitazioni e citare fewclust, con gli stessi controlli del manoscritto v4."
step_id: "S7f"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-10-31T14:40:00+01:00, evidence: "ok del PI in chat sul riassunto In breve (sessione S7g); commit S7f: approved by PI" }
resource: "../../reports/morris/S7f.html"
generated: { by: mmm-agent/skill-v4, at: 2026-10-21T21:46:43+02:00 }
sources:
  - id: d42
    resource: "../decisions/D-42.md"
    title: "Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale"
  - id: d43
    resource: "../decisions/D-43.md"
    title: "fewclust allineato alla regola del paper e citato nel manoscritto"
---
# Dubbi ricevuti

1. **Risultati nuovi da portare nel testo.** Ricevuto il 2026-10-21 dal piano approvato dal PI ([Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](../decisions/D-42.md)): risultati di [Rianalisi senza nuove simulazioni (S7b)](S7b.md), [Disegno e pilota dei test a livello di soggetto e di permutazione (S7c)](S7c.md), [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](S7d.md) e [Terzo meccanismo generatore (S7e)](S7e.md). Le frasi sui meccanismi come modelli di casa si scrivono qui.
2. **Citazione di fewclust.** Ricevuto il 2026-10-21 ([fewclust allineato alla regola del paper e citato nel manoscritto (D-43)](../decisions/D-43.md)). Il paper cita il pacchetto nella frase su codice e dati, nella discussione e, se le parole bastano, nell'abstract. Cancello: la versione citata, installata sul server da un commit registrato, applica la regola del paper (sui dati di prova sceglie il modello marginalizzato con il controllo di convergenza e la t, e come alternativa la GEE con Fay e Graubard e la t). Se il cancello fallisce, la citazione non entra e lo si dice al PI.
3. **Dal run con i metodi aggiunti (2026-10-23).** Ricevuto da [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](S7d.md). Cosa: (a) Figura 1 e testo secondo [La regola della Figura 1 si applica alle strategie con stima e intervallo: sceglie la GEE AR(1) con Fay e Graubard e la t; il test a inversione di segno e' complementare (D-50)](../decisions/D-50.md): GEE AR(1) con Fay e Graubard e t raccomandata, modello marginalizzato come alternativa equivalente, test a inversione di segno breve nel testo e in dettaglio nel supplemento (con l'intervallo di [Intervallo di confidenza dal test a inversione di segno (S7d2)](S7d2.md) se la copertura tiene); (b) Kauermann e Carroll da fewclust e nota sull'errore di geesmv nel supplemento ([La correzione di Kauermann e Carroll nel paper viene da fewclust; l'errore di geesmv si spiega nel supplemento (D-51)](../decisions/D-51.md)); (c) limiti della sezione 4.4 riscritti con i numeri di S7d (correlazione di lavoro: Fay e Graubard -0,3 punti con indipendenza e AR(1); fewclust contro geesmv: +0,02 punti e +0,4 punti di risposte; 2-3 visite e fino a 25 soggetti: il test risponde nel 99,9%); (d) fewclust citato nella versione registrata. Perche': cambia la raccomandazione principale. Cosa serve: tabelle `tables/V4_S7d_rule.csv`, `V4_S7d_boot_means.csv`, `V4_S7d_dp_rates.csv` e lo step S7d2 chiuso.

4. **Testo da scrivere dalla rianalisi senza nuove simulazioni.** Ricevuto il 2026-10-22 da [Rianalisi senza nuove simulazioni (S7b)](S7b.md) (step chiuso, non sono dubbi aperti): (a) gradi di liberta': N - p resta nel testo, N - 2 di Li e Redden come sensibilita' con una frase breve e una tabella nel supplemento ([Gradi di liberta' della t: N - p nel testo principale, N - 2 di Li e Redden come sensibilita' nel supplemento (D-45)](../decisions/D-45.md); `tables/V4_S7b_rule_summary.csv`, `tables/V4_S7b_df_shift.csv`, `tables/V4_S7b_df_flips.csv`); (b) la limitazione sul nullo marginale diventa quantitativa: rifiuto atteso in piu' al massimo 0,03 punti negli scenari principali e 0,26 in M6-M7, nessuna bandiera cambia (`tables/V4_S7b_null_summary.csv`); (c) seconda applicazione clinica su `respiratory` di geepack, paragrafo nel testo e tabella nel supplemento (`tables/V4_S7b_resp_ref.csv`, `V4_S7b_resp_usable.csv`, `V4_S7b_resp_bias.csv`, `V4_S7b_resp_compare.csv`), con la differenza dichiarata sulla specificazione del modello marginalizzato (`lv.formula = ~1` perche' le visite sono una serie temporale); (d) in entrambe le applicazioni il rapporto SE/SD con il fattore di popolazione finita sqrt(1 - N/P), formula nei metodi ([Il rapporto fra errore standard e dispersione nei sottocampioni si corregge per la popolazione finita (D-44)](../decisions/D-44.md); `tables/V4_S7b_seratio_fpc.csv`), che cambia i numeri della frase di sezione 3.7 (modello marginalizzato 0,85-0,99, modello misto 0,51-0,58 nella coorte neonatale); (e) la frase che inquadra GenBinaryY e GLMM-native come modello di casa, ciascuno, di un metodo; (f) nota su geesmv e la colonna `id` nel supplemento, accanto alla verifica numerica delle correzioni (`R/s7b/s7b_geesmv_check.R`).

5. **Dal terzo meccanismo generatore (2026-10-24).** Ricevuto da [Terzo meccanismo generatore (S7e)](S7e.md). Cosa: (a) il quarto meccanismo (copula gaussiana con margini logistici e correlazione solo seriale, [Il quarto meccanismo e' una copula gaussiana con margini logistici e correlazione solo seriale (rho = 0,6) (D-53)](../decisions/D-53.md)) nel testo principale come quarto meccanismo, con i dettagli nel supplemento ([Il quarto meccanismo entra nel testo principale: regola della Figura 1 sui quattro meccanismi, frase sul modello misto, dettagli nel supplemento (D-55)](../decisions/D-55.md)); (b) regola della Figura 1 sui 2.376 punti dei quattro meccanismi (GEE AR(1) con Fay e Graubard e t 92,1% di risposte, modello marginalizzato con controllo 90,9%; `tables/V4_S7e_rule.csv`); (c) frase sul modello misto: nel nuovo meccanismo 10,2% medio con la t, 219 punti su 480 chiaramente sopra 7,5%, errore standard della sola intercetta casuale fra 0,64 e 0,81 della dispersione (`tables/V4_S7e_summary.csv`, `V4_S7e_glmm_se_ratio.csv`); (d) GEE exchangeable con Fay e Graubard sopra 7,5% in 5 punti (T = 4, N 10-25); (e) limite sui meccanismi come modelli di casa riscritto, dichiarando che la correlazione AR(1) e' vicina al nuovo meccanismo. Cosa serve: tabelle `tables/V4_S7e_*.csv` e il report S7e.

6. **Dall'intervallo dal test a inversione di segno (2026-10-25).** Ricevuto da [Intervallo di confidenza dal test a inversione di segno (S7d2)](S7d2.md). Cosa: (a) una frase nella discussione e il dettaglio nel supplemento sull'intervallo per inversione ([L'intervallo per inversione del test a inversione di segno non si propone nel paper: una frase nella discussione e il dettaglio nel supplemento (D-58)](../decisions/D-58.md)): copertura 94,5% sotto [H1: Metodi in ordine ADEMP / tabella ADEMP](../audit/requirements/H1.md), mai fuori taratura; estremi finiti nel 79% dei dataset e nel 13,6% con 2-3 visite e fino a 25 soggetti (GEE AR(1) con Fay e Graubard e t usabile nel 47%); 17% piu' largo dove entrambi esistono (6% con almeno 75 soggetti); (b) numeri del test ricalcolati con le regole di [Nel paper il test a inversione di segno non e' calcolabile quando il fit ridotto diverge o lo score e' numericamente nullo (D-59)](../decisions/D-59.md) dalla colonna `p0_guard` di `results_v4/S7d2/rows/`: risposte 96,5% sotto H0 su 2.376 punti e 72,3% con 2-3 visite e fino a 25 soggetti (non 99,9%), errore di I tipo 5,4% e 5,2%, 0 punti fuori taratura; la frase "risponde quasi sempre" del dubbio 3 qui sopra va sostituita. Perche': cambia la frase sul limite delle 2-3 visite con pochi soggetti. Cosa serve: `tables/V4_S7d2_predictions.csv`, `V4_S7d2_flip_rules.csv`, `V4_S7d2_flip_rules_dp.csv`, `V4_S7d2_width_groups.csv`, `V4_S7d2_explain_finite.csv`, report `reports/morris/S7d2.html`.

# Piano

- Stessa pipeline del manoscritto v4: registro dei numeri esteso, controllo delle cifre scritte a mano, controllo delle frasi di confronto, conteggio delle parole (limite SMMR circa 8.000 con tabelle e figure a 200 parole; oggi 6.602), abstract entro 200 parole, controllo del formato Word.
- Nuova versione del manoscritto in `manuscript/v5/`, cosi' la v4 resta come riferimento.

Nota tecnica (2026-10-25): i dubbi ricevuti erano numerati 1, 2, 3, 3, 5, 6; il secondo 3 e' ora il 4, senza altre modifiche al testo.

# Attese pre-registrate (2026-10-25, prima di costruire il registro v5 e di scrivere il testo)

## Cosa si fa

1. **Manoscritto v5.** `manuscript/v5/MMM_v5.qmd` e `MMM_v5_supplement.qmd` partono dal testo v4 (commit 0ef16e7) e portano i sei dubbi ricevuti. La raccomandazione segue [La regola della Figura 1 si applica alle strategie con stima e intervallo: sceglie la GEE AR(1) con Fay e Graubard e la t; il test a inversione di segno e' complementare (D-50)](../decisions/D-50.md), con i numeri del test di [Nel paper il test a inversione di segno non e' calcolabile quando il fit ridotto diverge o lo score e' numericamente nullo (D-59)](../decisions/D-59.md).
2. **Kauermann e Carroll da fewclust in tutte le tabelle e figure** ([La correzione di Kauermann e Carroll nel paper viene da fewclust; l'errore di geesmv si spiega nel supplemento (D-51)](../decisions/D-51.md)). Tabelle 2, 3 e 5, le tabelle del supplemento e i pannelli b-d della Figura 1 vengono dalla catena S4-S6 (aggregazione, bootstrap per blocchi, tabelle), la stessa del v4. La si rifa' su un nuovo albero di variante, `kcfc`: le righe della variante con il controllo di convergenza (`results_v4/variant_cc/S4/rows_v4.rds`), con le righe `gee_kc` sostituite da errore standard, statistica e p della correzione di Kauermann e Carroll di fewclust 0.2.0 sulla GEE exchangeable (`results_v4/S7d/rows/`, `cor == "exch"`, `se_kc`), e l'usabilita' ricalcolata con la stessa funzione dello studio. Seme del bootstrap e B invariati (20261015, 2.000). Tutte le affermazioni si rigiudicano sui file della nuova variante.
3. **Regola della Figura 1 e Tabella 4 sui 2.376 punti dei quattro meccanismi**, ricostruite da `tables/V4_S7d_dp_rates.csv` (1.896 punti) e `tables/V4_S7e_dp_rates.csv` (480 punti), con il test a inversione di segno da `tables/V4_S7d2_flip_rules_dp.csv`.
4. **Registro dei numeri v5.** `R/s7f/s7f_extras.R` (valori nuovi) e `R/s7f/s7f_numbers.R` (registro `tables/V5_S7f_numbers.csv` dalla specifica `manuscript/v5/numbers_spec.csv`); controllo indipendente in Python come nel v4. I confronti nuovi si giudicano su `tables/V5_S7f_claims.csv`, che riporta stima e intervallo per blocchi gia' calcolati negli step precedenti (`tables/V4_S7d_boot_means.csv`, `V4_S7e_boot_means.csv`, `V4_S7d2_boot_means.csv`) con il margine di default di questa skill (differenze medie di errore di I tipo o copertura +/- 1 punto; bias / SE empirico +/- 0,46).

Le tabelle 2, 3 e 5 restano sugli otto metodi dello studio (Kauermann e Carroll da fewclust), con il riferimento normale di [Tabelle dei risultati con riferimento normale, raccomandazione con riferimento t (D-37)](../decisions/D-37.md); le strategie aggiunte (tre correlazioni di lavoro, CR2, test a inversione di segno) entrano nella Tabella 4 e nel supplemento. Decisione tecnica, registrata con questa pre-registrazione.

## Cancelli

- **GF1 (integrita').** md5 di `results_full/` (28.800 file) identici prima e dopo lo step; md5 delle cartelle di righe lette (`results_v4/S7d/rows`, `results_v4/S7e/rows`, `results_v4/S7d2/rows`, `results_v4/variant_cc/S4/rows_v4.rds`) identici prima e dopo; md5 di ogni file sorgente del registro identico alla costruzione e al render finale.
- **GF2 (numeri attestati).** Come G2 del v4 ([Riscrittura del manoscritto (S7)](S7.md)): zero cifre fuori dalle classi ammesse nel manoscritto e nel supplemento resi in modalita' traccia. Classi ammesse invariate, con due aggiunte fissate ora: il parametro del quarto meccanismo scritto come livello (rho = 0.6, la sua definizione in [Il quarto meccanismo e' una copula gaussiana con margini logistici e correlazione solo seriale (rho = 0,6) (D-53)](../decisions/D-53.md)) e le soglie del test non calcolabile (1e3, 1e-6) e del filtro di sensibilita' (10), che sono soglie di disegno.
- **GF3 (ricevute e ricostruzioni).** (a) Il controllo indipendente rilegge il 100% delle voci del registro: valore identico (|differenza| < 1e-9 sui numeri) e md5 identico. (b) La variante `kcfc` differisce da `variant_cc` solo nelle righe `gee_kc` (zero differenze nelle altre colonne e righe); in quelle righe la stima e' quella dello studio (|differenza| < 1e-8 dove entrambe esistono) e l'errore standard e' `se_kc` di fewclust. (c) La Tabella 4 ricostruita dai file per punto coincide con `tables/V4_S7e_rule.csv` (numero di punti chiaramente sopra 7,5% e sotto 2,5%, quota usabile media, potenza media, per ogni strategia in comune, |differenza| < 1e-9) e, per il test, con `tables/V4_S7d2_flip_rules.csv`.
- **GF4 (affermazioni).** Come G4 del v4: ogni frase comparativa porta un'etichetta che si risolve su `tables/variant_kcfc/V4_S5_claims.csv`, su `tables/V4_S6_m1c_claims.csv` (rigiudicato sulla variante se contiene Kauermann e Carroll), su `tables/V5_S7f_claims.csv` o su una voce del registro; zero frasi comparative senza etichetta; zero formulazioni in contrasto con la classe dopo la revisione riga per riga.
- **GF5 (completezza).** `quarto render` di manoscritto e supplemento in `.docx` e `.html` senza errori sul server; formato Word come [Formato del Word: Arial, tutto nero, testo a 11 punti (D-41)](../decisions/D-41.md); `tools/okf.py check` passa; abstract non strutturato di al massimo 200 parole; ognuno dei 6 dubbi ricevuti ha un esito (spiegato, corretto, deciso o rinviato).

## Condizione per la citazione di fewclust (dal dubbio 2, non e' un cancello dello step)

La versione citata come implementazione della raccomandazione, installata sul server da un commit registrato, applicata con `auto_cluster_fit()` ai dati di prova del pacchetto e a un dataset simulato con 25 soggetti e 4 visite, sceglie la strategia del paper di oggi (GEE con correlazione AR(1), Fay e Graubard e t; alternativa il modello marginalizzato con il controllo). La condizione del dubbio era scritta sulla regola del v4 (modello marginalizzato per primo), superata da [La regola della Figura 1 si applica alle strategie con stima e intervallo: sceglie la GEE AR(1) con Fay e Graubard e la t; il test a inversione di segno e' complementare (D-50)](../decisions/D-50.md): vale la regola in vigore. Attesa dichiarata ora: con fewclust 0.2.0 (commit 8538cad, l'unico registrato) la condizione **fallisce**, perche' il codice sceglie `"mmm"` per costruzione (`R/decision.R`, passo D5). Se fallisce, la citazione come implementazione della raccomandazione non entra (frase su codice e dati, discussione, abstract) e la scelta va al PI. Resta in ogni caso la menzione di fewclust 0.2.0 come software che ha calcolato le correzioni (nei metodi e nella Tabella S2, come geesmv), perche' e' un fatto del calcolo e non una raccomandazione.

## Previsioni

- **PF1 (lunghezza).** Parole del testo principale (dall'introduzione alle conclusioni, esclusi abstract, riferimenti, tabelle e didascalie) piu' 200 per tabella o figura del testo principale: attese 7.800; intervallo di rilevanza 7.000-9.000. Sopra 9.000 serve un taglio o una lettera all'editor; sotto 7.000 e' probabile che manchi una delle aggiunte dei sei dubbi. Conteggio esatto, senza MCSE.
- **PF2 (registro della voce).** Mediana della lunghezza delle frasi del testo principale (`sentence_stats`): attesa 19 parole; intervallo 14-24; zero lineette em o en; zero aperture di paragrafo con connettivi di riempimento.
- **PF3 (portata della sostituzione di Kauermann e Carroll).** Numero di sotto-affermazioni di `V4_S5_claims.csv` che non contengono `gee_kc` nel calcolo e cambiano classe passando da `variant_cc` a `kcfc`: atteso 0; intervallo [0, 0], perche' ogni cambio riscriverebbe una frase che non riguarda quella correzione. Le sotto-affermazioni che la contengono (per esempio la potenza grezza piu' alta e la copertura piu' lontana dal nominale) possono cambiare: si elencano e si riscrivono.
- **PF4 (livello di Kauermann e Carroll da fewclust sulla griglia principale, riferimento normale).** Tasso medio di rifiuto sotto H0, pooled, su GenBinaryY: atteso piu' basso di quello di geesmv nel v4 di almeno 1 punto (in S7d, con la t e sui 1.896 punti, la differenza media e' -4,2 punti); intervallo di rilevanza: differenza sotto -1 punto (la frase del v4 "il metodo piu' anticonservativo" cambia). Giudicata con l'intervallo per blocchi della variante.

# Esito delle attese (2026-10-31)

Cancelli (log in `reports/morris/S7f_logs/`, ultimo giro della pipeline `R/s7f/s7f_pipeline.sh` il 2026-10-27 alle 18:13):

- **GF1 superato.** `results_full/`: 28.800 file, md5 identici prima e dopo lo step e uguali alla lista di fine S7b (`results_v4/S7f/md5_results_full_{before,after}.txt`). Righe lette (`results_v4/S7d/rows`, `S7e/rows`, `S7d2/rows`, `variant_cc/S4/rows_v4.rds`): 95.041 file, md5 identici prima e dopo (digest 58b9089f2354441fc2150448b5c8362f). Gli md5 dei file sorgente del registro sono ricontrollati a ogni render da GF3a.
- **GF2 superato.** 0 cifre fuori dalle classi ammesse, nel manoscritto e nel supplemento resi in traccia.
- **GF3 superato.** (a) 2.297 voci del registro rilette da 48 file, 0 discordanze. (b) La variante `kcfc` differisce da `variant_cc` solo nelle righe `gee_kc`, con la stessa stima e l'errore standard di fewclust (`tables/V5_S7f_kcfc_gate.csv`). (c) La Tabella 4 ricostruita coincide con `tables/V4_S7e_rule.csv` su 17 strategie (usabilita', potenza, scelta) e, per il test, con `tables/V4_S7d2_flip_rules.csv` (`tables/V5_S7f_GF3c.csv`). Dopo il passaggio agli intervalli esatti ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)) i conteggi sopra e sotto cambiano per costruzione: per questi il controllo e' che la regola di Wald riproduca i flag salvati prima della sostituzione (25 controlli, `tables/V5_S7f_cp_gates.csv`).
- **GF4 superato.** 115 etichette di affermazione, 0 non risolte, 0 frasi comparative senza etichetta.
- **GF5 superato.** Render di manoscritto e supplemento in `.docx` e `.html` senza errori; formato Word conforme; abstract di 200 parole; ognuno dei 6 dubbi ricevuti ha un esito (sezione seguente); controllo OKF al commit.

Previsioni:

- **PF1 confermata.** 7.458 parole + 7 tabelle o figure x 200 = 8.858, dentro 7.000-9.000. E' sopra le circa 8.000 di SMMR scritte nel piano: il taglio o la lettera all'editor passano alla revisione finale.
- **PF2 confermata.** Mediana 20,0 parole per frase (344 frasi), dentro 14-24; 0 lineette em o en; la revisione di stile non ha trovato aperture con connettivi di riempimento.
- **PF3 confermata.** 0 sotto-affermazioni senza Kauermann e Carroll nel calcolo cambiano classe (`tables/V5_S7f_PF.csv`); le 29 righe che cambiano classe (`tables/V5_S7f_claims_diff.csv`, colonna `class_change`) contengono tutte quella correzione e sono state riscritte.
- **PF4 confermata.** Kauermann e Carroll da fewclust contro geesmv, GenBinaryY, riferimento normale: -6,2 punti (intervallo per blocchi da -6,4 a -6,0), sotto -1. La frase del v4 sul metodo "most anticonservative" non c'e' piu'.

Condizione per la citazione di fewclust: **fallita come atteso**. `auto_cluster_fit()` di fewclust 0.2.0 (commit 8538cad) sceglie `"mmm"` sui dati `fc_trial` e su un dataset simulato con 25 soggetti e 4 occasioni (`R/s7f/s7f_fewclust_condition.R`, `tables/V5_S7f_fewclust_condition.csv`). Deciso dal PI: [fewclust non si cita come implementazione della raccomandazione; resta nominato come software che ha calcolato le correzioni (D-65)](../decisions/D-65.md).

# Esito dei dubbi ricevuti

1. **Risultati nuovi nel testo: corretto.** Manoscritto e supplemento v5 (`manuscript/v5/`).
2. **Citazione di fewclust: deciso** ([fewclust non si cita come implementazione della raccomandazione; resta nominato come software che ha calcolato le correzioni (D-65)](../decisions/D-65.md)).
3. **Dal run con i metodi aggiunti: corretto.** (a) Figura 1 e testo sulle strategie della regola, con il modello marginalizzato come alternativa comparabile ([Figura 1 sulle strategie della regola; la GEE AR(1) con Fay e Graubard e t va prima, il modello marginalizzato con controllo e' un'alternativa comparabile (D-61)](../decisions/D-61.md)); (b) Kauermann e Carroll da fewclust e nota su geesmv nel supplemento ([La correzione di Kauermann e Carroll nel paper viene da fewclust; l'errore di geesmv si spiega nel supplemento (D-51)](../decisions/D-51.md)); (c) limiti riscritti con i numeri di S7d e del test secondo [Nel paper il test a inversione di segno non e' calcolabile quando il fit ridotto diverge o lo score e' numericamente nullo (D-59)](../decisions/D-59.md); (d) fewclust nominato nella versione registrata come software.
4. **Testo della rianalisi senza nuove simulazioni: corretto** (a-f).
5. **Dal terzo meccanismo: corretto** (a-e); i conteggi "chiaramente sopra/sotto" sono ora con l'intervallo esatto; la correlazione del quarto meccanismo sugli esiti binari e' dichiarata (0,36 fra occasioni adiacenti, 0,20 a distanza due).
6. **Dall'intervallo dal test a inversione di segno: corretto** (a, b).

# Dubbi emersi durante lo step

- **Figura 1: pannelli che contraddicevano il testo.** I pannelli b e d mostravano gli otto metodi della Tabella 2, non le due strategie della regola. Corretto ([Figura 1 sulle strategie della regola; la GEE AR(1) con Fay e Graubard e t va prima, il modello marginalizzato con controllo e' un'alternativa comparabile (D-61)](../decisions/D-61.md)).
- **Confronto diretto letto senza margini e con un errore Monte Carlo di indipendenza.** Rifatto con il bootstrap per blocchi (B = 2.000, seme 20261015) e i margini della Sez. 2.5 (`R/s7f/s7f_h2h_boot.R`, `tables/V5_S7f_h2h_boot.csv`, Tabella S31). Corretto; la presentazione e' decisa dal PI ([Figura 1 sulle strategie della regola; la GEE AR(1) con Fay e Graubard e t va prima, il modello marginalizzato con controllo e' un'alternativa comparabile (D-61)](../decisions/D-61.md)).
- **Revisione avversariale.** 2.249 numeri riletti identici; i rilievi e il loro esito sono in `tables/V5_S7f_review_outcomes.csv` (corretti, spiegati, decisi dal PI o rinviati; nessuno aperto). Report del revisore: artifact `revisione_n2_report.md`.
- **Flag "chiaramente sotto 2,5%" con zero rifiuti.** Deciso ([Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](../decisions/D-63.md)).
- **Tabelle 2, 3 e 5 non ricalcolate dal revisore (modello marginalizzato 5,65 contro 5,54).** Spiegato: ricostruite per intero dalle singole simulazioni (72 + 96 + 48 celle, `R/s7f/s7f_table{2,3,5}_recompute.R`), differenza massima sotto 1e-12; il 5,65 veniva dalla selezione delle righe del revisore; in Tabella 3 la differenza iniziale era nelle celle senza righe.
- **Estrazioni della coorte: 21,8% senza ventilati invasivi contro 26,8% atteso a 10 pazienti.** Spiegato: le 500 estrazioni sono campioni casuali semplici (`sample(ids, n)`, seme 20260807); test binomiale p = 0,012 su una sola delle cinque dimensioni, a 15, 20, 25 e 50 pazienti coerenti (`R/s7f/s7f_lus_draws_check.R`). Nel testo l'atteso e' scritto accanto all'osservato.
- **Nuove simulazioni proposte dal revisore.** Deciso ([Nessuna nuova simulazione dopo la revisione avversariale; i metodi non valutati si dichiarano fra i limiti (D-64)](../decisions/D-64.md)).
- **GEE penalizzata alla Firth.** Il PI ricordava di averla gia' provata; nel repo e nelle sessioni non c'e' traccia (solo il suggerimento di fewclust quando trova separazione e la proposta della revisione del 2026-08-09). Il 2026-10-31 il PI ha confermato di essersi sbagliato e che non si prova. Nei limiti resta la frase che la dichiara non valutata ([Nessuna nuova simulazione dopo la revisione avversariale; i metodi non valutati si dichiarano fra i limiti (D-64)](../decisions/D-64.md)).
- **Revisione di stile.** 75 modifiche approvate dal PI applicate, piu' una sezione 4.4 sulle applicazioni cliniche e sul quarto meccanismo; bibliografia rinumerata per prima citazione (64 citazioni, stesse 31 opere). Deciso ([Revisione di stile applicata; il supplemento v5 diventa la fonte (D-62)](../decisions/D-62.md)).
- **Pannello a della Figura 1 non leggibile.** Corretto con tre etichette ([Figura 1 sulle strategie della regola; la GEE AR(1) con Fay e Graubard e t va prima, il modello marginalizzato con controllo e' un'alternativa comparabile (D-61)](../decisions/D-61.md)).
- **fewclust stampa le prime righe della matrice del modello quando un fit fallisce.** Non tocca il paper; nota per il PI, che aggiorna il pacchetto in una chat separata.

# Dubbi rinviati

Passano alla [Revisione avversariale e ri-audit Morris (S8)](S8.md), dubbi ricevuti 5-10:

- confronto con Bie et al. [18] (punti 15 e 24 della revisione di stile; rilievo sulla stima fra cluster del modello marginalizzato), perche' la fonte non e' stata letta in questo step;
- regola "et al." della bibliografia secondo le istruzioni SAGE/SMMR (punto 31);
- stime della GEE nella sezione 3.10 senza errore standard e senza correlazione di lavoro (punto 22);
- media di Mancl-DeRouen e Kauermann-Carroll di Ford e Westgate [14], gia' calcolata dal revisore (esclusa dalla regola): decidere se una riga nel supplemento;
- testo di Bradley (1978), non disponibile al revisore;
- lunghezza: 8.858 parole contando 200 per tabella o figura, contro circa 8.000.

# Decisioni tecniche dell'agente

- Tabelle 2, 3 e 5 sugli otto metodi dello studio con Kauermann e Carroll da fewclust (registrata nella pre-registrazione).
- Confronto diretto delle due strategie con il bootstrap per blocchi e i margini della Sez. 2.5, come regola 6 del protocollo.
- [Revisione di stile applicata; il supplemento v5 diventa la fonte (D-62)](../decisions/D-62.md).
- Riassunto «In breve» con i numeri riletti da `R/s7f/s7f_inbreve_sources.R` (`tables/V5_S7f_inbreve_sources.csv`, 28 valori su 28 uguali al testo).

# In breve

## In breve: il manoscritto riporta bene le analisi aggiunte, e i suoi numeri tornano?

**A cosa serviva.** Abbiamo portato nel manoscritto i risultati delle quattro analisi fatte dopo la prima stesura: la rianalisi senza nuove simulazioni, il run con i metodi aggiunti, il quarto meccanismo generatore e l'intervallo ottenuto dal test a inversione di segno. Poi il testo è passato da un revisore avversariale, che ha ricalcolato i numeri, e da un revisore di stile. Conta perché da qui il paper va alla revisione finale: ogni numero e ogni confronto devono venire dai file.

**Cosa abbiamo trovato.**

1. **La raccomandazione resta la GEE con correlazione AR(1), la correzione di Fay e Graubard e la t; il modello marginalizzato multilivello con il controllo di convergenza è un'alternativa comparabile.**
   - Su 2.376 situazioni simulate con quattro modi diversi di generare i dati, nessuna delle due strategie è mai chiaramente sopra 7,5 falsi positivi su 100.
   - La GEE dà un risultato nel 92% dei dataset, il modello marginalizzato nel 91%: la differenza è sotto la soglia di importanza, quindi sul criterio della regola sono pari.
   - La GEE va prima per due ragioni. È l'unica differenza importante: con 10 soggetti e 3 o più visite il modello marginalizzato è conservativo (2,1 falsi positivi su 100 contro 4,0). E la GEE non richiede il controllo di convergenza.
2. **Il confronto diretto fra le due strategie l'avevo letto senza soglia di importanza; rifatto con le soglie del paper, le due sono quasi ovunque pari.**
   - In media la distanza dal 5% è la stessa: 1,1 punti per entrambe.
   - Con 2 visite e almeno 15 soggetti il modello marginalizzato è più vicino al 5% di circa un punto, cioè proprio la soglia: non lo scriviamo più "preferibile". Lì la GEE è conservativa (2,5 falsi positivi su 100).
   - Nella Figura 1 il pannello a dice "GEE" solo con 10 soggetti e 3 o più visite, e "GEE = MMM" altrove.
3. **Il revisore avversariale ha trovato i numeri giusti, ma alcune frasi dicevano più dei dati; le abbiamo corrette.**
   - "Mantengono il livello ovunque" è diventato "mai chiaramente sopra 7,5%".
   - Il test di Wald del modello misto (il modello logistico a effetti casuali per soggetto, stimato con l'approssimazione di Laplace di lme4) è fuori taratura anche sui dati generati dal suo stesso modello, in 55 situazioni su 696, non solo quando le visite vicine sono più correlate.
   - Il test a inversione di segno, con 2-3 visite e fino a 25 soggetti, trova un effetto vero in circa 4 dataset su 100: non scriviamo più che "dice se c'è un effetto".
   - "Chiaramente sotto 2,5%" ora si giudica con un intervallo esatto. Prima una situazione con zero falsi positivi contava come chiaramente sotto anche con pochi dataset utili. Le situazioni scendono da 136 a 116 per la GEE e da 86 a 24 per il modello marginalizzato; la scelta della regola non cambia.
4. **Le Tabelle 2, 3 e 5, ricostruite dalle singole simulazioni, coincidono con quelle stampate.**
   - Il valore diverso trovato dal revisore per il modello marginalizzato veniva da una sua selezione delle simulazioni.
   - Tutti i 2.297 numeri del manoscritto e del supplemento si rileggono dai file senza discordanze.
5. **Il testo conta 8.858 parole (200 per ogni tabella o figura), sopra le circa 8.000 della rivista.** L'abstract è di 200 parole.

**Cosa devi decidere.** Nulla.

**Comunque:**

- Quando faremo la revisione finale: il confronto con Bie e colleghi (va letta la fonte), la regola "et al." della bibliografia secondo le istruzioni dell'editore, l'errore standard delle stime della GEE nelle applicazioni cliniche, la media delle correzioni di Mancl e DeRouen e di Kauermann e Carroll proposta da Ford e Westgate (già calcolata, da mettere o no nel supplemento), il testo di Bradley e la lunghezza.
- fewclust, nella versione registrata, sceglie ancora per primo il modello marginalizzato: l'ho provato sui dati del pacchetto e su un dataset simulato con 25 soggetti e 4 visite. Nel paper resta quindi nominato solo come il software che ha calcolato le correzioni. Quando lo aggiornerai, ti segnalo anche che negli errori stampa le prime righe della matrice del modello.
- La GEE penalizzata alla Firth non si prova. Nei limiti resta una frase che dice che non è stata valutata, insieme alla quadratura adattiva per il modello misto e ai gradi di libertà proposti da Fay e Graubard.

Non c'è nulla da decidere: se non hai dubbi, lo step è chiuso.
