---
type: Step
title: "Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto"
description: "Rigenerare con gli stessi semi i dataset di tutti i 1.896 punti di disegno e applicarvi correlazioni di lavoro alternative, correzioni del sandwich da implementazioni indipendenti e i test scelti nel pilota."
step_id: "S7d"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-10-23T09:30:00+02:00, evidence: "chat del 2026-10-23: \"si va bene ma il test a inversione ne parlerei nel testo (brevemente) e poi nei supplementari nello specifico\", \"A\", \"Provarlo in uno step dedicato prima della riscrittura\" (D-50, D-51, D-52); commit S7d: closed, approved by PI" }
resource: "../../reports/morris/S7d.html"
generated: { by: mmm-agent/skill-v4, at: 2026-10-23T06:30:00+02:00 }
sources:
  - id: d42
    resource: "../decisions/D-42.md"
    title: "Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale"
  - id: s7c
    resource: "S7c.md"
    title: "Disegno e pilota dei test a livello di soggetto e di permutazione"
  - id: kc2001
    resource: https://doi.org/10.1198/016214501753381896
    title: "Kauermann, Carroll (2001). A note on the efficiency of sandwich covariance matrix estimation. JASA 96:1387-1396"
  - id: md2001
    resource: https://doi.org/10.1111/j.0006-341X.2001.00126.x
    title: "Mancl, DeRouen (2001). A covariance estimator for GEE with improved small-sample properties. Biometrics 57:126-134"
  - id: fg2001
    resource: https://doi.org/10.1111/j.0006-341X.2001.01198.x
    title: "Fay, Graubard (2001). Small-sample adjustments for Wald-type tests using sandwich estimators. Biometrics 57:1198-1206"
---
# Dubbi ricevuti

1. **Correlazione di lavoro, implementazione delle correzioni e test a livello di soggetto, dove contano.** Ricevuto il 2026-10-21 dal piano approvato dal PI ([Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](../decisions/D-42.md)). Limitazioni della sezione 4.4 del manoscritto v4: la GEE usa sempre la correlazione di lavoro exchangeable; le correzioni vengono solo da geesmv, che con pochi soggetti fallisce o va in time-out; con 2-3 visite e fino a 25 soggetti nessuna strategia e' insieme disponibile e tarata.
2. **Dal pilota dei test (2026-10-22).** Ricevuti da [Disegno e pilota dei test a livello di soggetto e di permutazione (S7c)](S7c.md), sezione "Dubbi rinviati": (a) geesmv con `corstr = "AR-M"` risponde in circa 0,02 s e non va mai in time-out: verificare che restituisca la correzione, con il cancello di coincidenza geesmv / fewclust; (b) il CR2 di clubSandwich fallisce sulla GEE AR(1) nel 59% dei dataset dei tempi: fissare prima del run come contarlo; (c) la versione di fewclust allineata al paper non e' in un commit registrato (nel pilota: 0.1.0, commit `5103bb3` di `~/MMM-package`); (d) geesmv e' il 93% del tempo stimato (circa 3.300 ore di processore su 1.896 punti): valutare con il PI fewclust su tutta la griglia e geesmv sul solo sottoinsieme del cancello; misurare il tempo di generazione dei dati; (e) scelta del PI: entra il solo test a inversione di segno ([Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](../decisions/D-46.md)), principale senza filtro sui coefficienti del modello ridotto e con il filtro come sensibilita' ([Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](../decisions/D-47.md)); la versione senza filtro e' nata dopo i risultati del pilota, quindi le sue attese si scrivono qui prima del run.

# Piano

- Griglia: tutti i 1.896 punti di disegno del paper (griglia principale, M1c, M6-M7); H0 e H1; 1.000 replicati per braccio; stessi semi (`20260805 + 1000m + i`, L'Ecuyer-CMRG). La griglia e' intera perche' la regola della Figura 1 chiede che una strategia non sia mai chiaramente sopra 7,5% in nessun punto: una strategia valutata solo su parte della griglia non si puo' confrontare con le altre (indicazione del PI del 2026-10-21: i tempi di calcolo non limitano la scelta).
- Cancello forte: sui dataset rigenerati la GEE exchangeable con le correzioni di geesmv riproduce le stime salvate replicato per replicato; parita' H0/H1 prima del run; md5 dei risultati grezzi.
- Metodi aggiunti: GEE con correlazione di lavoro indipendente e AR(1); correzioni di Fay e Graubard, Mancl e DeRouen e Kauermann e Carroll con fewclust (versione allineata al paper, commit registrato) e CR2 di clubSandwich; i test scelti nel pilota. Dove geesmv e fewclust funzionano entrambi, le correzioni devono coincidere (cancello).
- Analisi: tutte le misure con MCSE, regola della Figura 1 riapplicata con le strategie nuove, confronti per blocchi di replicato (B = 2000, seed 20261015).
- Tempo: misurato nel pilota solo per organizzare il run sul server (processi, limite di temperatura); non limita griglia ne' metodi.

# Preparazione prima del run (2026-10-22, nessuna misura entra nel paper)

Sonda `R/s7d/s7d_probe.R`: 990 dataset rigenerati con i semi dello studio (22 combinazioni scenario x meccanismo; N 10, 25, 145, in M1c 10, 25, 45; T 2, 3, 10, in M1c 2, 3, 6; replicati 1-5 sotto H0), un processo per combinazione, operazioni in serie. Controprova a mano `R/s7d/s7d_check_formulas.R` su 36 dataset di M1. Risultati in `results_v4/S7d/probe/` (fuori da git).

## Chiusura dei dubbi ricevuti 2(a)-2(d)

- **(a) geesmv con AR(1): spiegato, non utilizzabile.** La risposta in 0,02 s e' un errore immediato: geesmv riadatta il modello con `gee::gee(corstr = "AR-M")`, che si ferma con "M-dependence, M=1, but clustsize=1" appena un soggetto ha una sola visita (702 degli 858 fallimenti con AR(1) su 986 dataset; quasi tutti gli altri: valori attesi vicini a 1). Quota con una correzione: 13% (0,3% con T = 2; 0% con N = 45 e 145). Dove risponde, l'errore standard e' vicino a quello di fewclust (rapporto mediano 1,000 per Fay e Graubard) ma con dispersione, per il diverso stimatore della correlazione di gee. Con AR(1) le correzioni vengono quindi solo da fewclust (risponde nel 78% dei dataset della sonda, come con exchangeable).
- **(b) CR2 con AR(1): deciso (decisione tecnica sotto).** clubSandwich fallisce sulla GEE AR(1) nel 69% dei dataset della sonda ("leading minor not positive"); dove risponde, il suo CR3 non coincide con Mancl e DeRouen calcolato a mano sullo stesso fit (rapporto da 0,998 a 1,071 su 9 dataset), mentre con exchangeable e indipendente coincide esattamente (12/12 e 12/12). La sua matrice di lavoro AR(1) non e' quella di geepack; CR2 si calcola solo con exchangeable e indipendente.
- **(c) versione di fewclust: corretto.** fewclust 0.2.0 (la versione allineata al paper, che applica la formula pubblicata di Fay e Graubard) era solo una copia di lavoro sul portatile del PI. Registrata come commit `8538cad` del ramo `package` e installata in `~/R/s7d_lib` ([fewclust 0.2.0 registrato dall'agente nel ramo package cosi' com'era sul portatile del PI (D-49)](../decisions/D-49.md)).
- **(d) organizzazione del run: deciso dal PI.** Tempi misurati (secondi per dataset, media sui punti della sonda): generazione 0,011-0,017; geeglm 0,004-0,007 per correlazione; fewclust 0,005-0,011 per correzione; CR2 0,019-0,035; test a inversione di segno 0,005-0,007; geesmv 0,09-0,25 per correzione (time-out compresi). Riportati ai 1.896 punti x 2 bracci x 1.000 replicati (punto misurato piu' vicino): generazione 18 ore di processore, geeglm 21, fewclust circa 135, CR2 circa 145 (circa 95 senza AR(1)), test 8; geesmv circa 1.400-1.700. Il PI ha scelto geesmv sui soli replicati 1-50 ([geesmv solo sulle prime 50 simulazioni di ogni scenario nel run con i metodi aggiunti (D-48)](../decisions/D-48.md)): run di circa 380 ore di processore, circa 16 ore con 24 processi.

## Anomalia trovata: la correzione di Kauermann e Carroll di geesmv non e' quella dell'articolo

Kauermann e Carroll premoltiplicano i residui del soggetto per (I - H_ii)^(-1/2), con H_ii = D_i B^(-1) D_i' V_i^(-1), matrice non simmetrica[^kc2001]. geesmv 1.3 (`GEE.var.kc`, funzione interna `mat.sqrt.inv`) calcola la radice con `eigen()` e ricompone con la trasposta degli autovettori, operazione valida solo per matrici simmetriche; inoltre azzera gli autovalori negativi. Controprova su 36 dataset di M1 (`R/s7d/s7d_check_formulas.R`): fewclust coincide con la radice principale calcolata a mano (Q diag(l^(-1/2)) Q^(-1)) in tutti i 36 casi e nelle tre correlazioni; con indipendenza (stesso fit nei due pacchetti) la formula di geesmv riprodotta a mano coincide con l'output di geesmv entro lo 0,5% in 9 dataset su 12 (negli altri 3 lo scarto e' del 3%, 20% e 78%, non spiegato in questa sonda); geesmv coincide entro l'1% con la radice principale in 1 dataset su 24; il rapporto fra l'errore standard di geesmv e quello corretto va da 0,58 a 3,85. Mancl e DeRouen e Fay e Graubard non hanno radici di matrici: con indipendenza geesmv e fewclust coincidono in tutti i 731 e 732 dataset della sonda (scarto relativo massimo 2e-12). Il paper riporta Kauermann e Carroll da geesmv: con la t, 9,75% di falsi positivi in media e 951 punti chiaramente sopra 7,5% (`tables/V4_S7b_rule_dp.csv`). Il run calcola la versione corretta su tutta la griglia; come riportarla nel paper e' una scelta del PI da portare alla chiusura dello step, con i numeri.

## Decisioni tecniche (agente, prima del run)

- geesmv con AR(1) non si esegue (dubbio 2(a)): fallimento strutturale del pacchetto gee con soggetti a una sola visita.
- CR2 di clubSandwich solo con exchangeable e indipendente (dubbio 2(b)); con AR(1) non si calcola e non entra in nessuna tabella.
- Modello ridotto del test a inversione di segno: il modello medio della GEE dello scenario senza `time:trt`, cioe' `binY ~ time + trt` piu' `age` in M2 e `dose` in M5. Negli altri scenari coincide con quello di [Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](../decisions/D-46.md). Il nullo resta quello del paper (interazione nulla nel modello medio dello scenario).
- GEE: `geepack::geeglm` con il modello medio dello studio, `scale.fix = TRUE`, dati ordinati per soggetto e tempo; AR(1) sull'ordine delle visite nel soggetto (le serie sono troncate dall'inizio, quindi ordine e tempo coincidono).
- Usabilita' delle strategie di Wald nuove: come nello studio (stima, errore standard e p-value finiti, stima in valore assoluto sotto 10, errore standard fra 0 e 10); la stima e' quella di geeglm della stessa correlazione, l'errore standard quello della correzione. Test a inversione di segno: p-value finito. Denominatore: replicati con disegno stimabile ([Replicati con disegno non stimabile esclusi dal denominatore (D-11)](../decisions/D-11.md)).
- Riferimento dei test di Wald nuovi: t con N - p gradi di liberta', come nella regola della Figura 1 ([Gradi di liberta' della t: N - p nel testo principale, N - 2 di Li e Redden come sensibilita' nel supplemento (D-45)](../decisions/D-45.md)); N - 2 tabulato come sensibilita'.
- Strategia del paper di confronto per le correzioni: le righe salvate dello studio (`results_v4/variant_cc/S4/rows_v4.rds` e `results_v4/S6b/ns/rows_ns.rds`, metodi `gee_md`, `gee_kc`, `gee_fg`), su tutti i replicati. Nota di metodo: in quelle righe la stima viene da geeglm e la varianza dal riadattamento di gee dentro geesmv; con exchangeable i due fit differiscono poco (scarto mediano della stima 2e-4).
- Run: `R/s7d/s7d_run.R`, chunk di 100 replicati in `results_v4/S7d/rows/` (fuori da git), riprendibile; processi fra 20 e 28 regolati dalla temperatura (`R/s7d/s7d_thermal_guard.sh`, copia della guardia del run di S6b); `OMP_NUM_THREADS = 1`. Test a inversione di segno con `set.seed(20261022 + i)` dopo la generazione e dopo i fit, come nel pilota.

# Disegno del run (fissato prima del calcolo)

Per ogni replicato di ognuna delle 3.792 celle (1.896 punti x H0, H1): rigenerazione con i semi dello studio; tre GEE (exchangeable, indipendente, AR(1)); per ciascuna stima, sandwich non corretto di geepack, Mancl e DeRouen, Kauermann e Carroll e Fay e Graubard di fewclust 0.2.0, CR2 (exchangeable e indipendente); geesmv (tre correzioni, exchangeable e indipendente) sui replicati 1-50; test a inversione di segno, principale e con filtro.

Strategie confrontate nella regola della Figura 1 (tutte con la t a N - p): le 7 di Wald dello studio dalle righe salvate, piu' le nuove: GEE con indipendenza e AR(1) con il sandwich non corretto e le tre correzioni di fewclust; GEE exchangeable con le tre correzioni di fewclust; CR2 con exchangeable e indipendente; il test a inversione di segno (solo p-value). Regola come in [La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](../decisions/D-35.md): si escludono le strategie chiaramente sopra 7,5% (tasso - 2 MCSE > 7,5, con almeno 20 risultati usabili) in almeno un punto, fra le rimaste si sceglie quella che risponde piu' spesso.

# Attese pre-registrate (scritte prima di ogni calcolo del run)

## Cancelli

| codice | controllo | passa se |
|---|---|---|
| GE1 | integrita' e completezza | md5 di `results_full/` prima e dopo il run uguali fra loro e alla lista di S7b; `R/testing/test_h0h1_parity.R` passa prima del run; 37.920 file di chunk, 3.792.000 replicati, ognuno con la sua riga e, se generato, tre righe di GEE; ogni errore di generazione cade in un disegno non stimabile |
| GE2 | stessi dataset dello studio | GEE exchangeable: stima ed errore standard di geepack uguali a quelli salvati (metodo `gee`) entro 1e-8 in ogni replicato in cui la riga salvata e' finita; geesmv exchangeable, replicati 1-50: errori standard di Mancl e DeRouen, Kauermann e Carroll e Fay e Graubard uguali a quelli salvati entro 1e-8 dove entrambi sono finiti |
| GE3 | stesso test del pilota | nei 18 punti del pilota (M1, N 10 e 25, T 2, 3 e 10, tre meccanismi, H0 e H1) il p-value del test principale coincide con `p_flip_relaxed` di `results_v4/S7c/posthoc/` e quello con filtro con il test A di `results_v4/S7c/pilot/`, replicato per replicato (uguali, oppure mancanti insieme) |
| GE4 | controprova indipendente delle formule, prima del run | `R/s7d/s7d_gate_formulas.R` su 880 dataset (22 combinazioni scenario x meccanismo, N 10 e 25, T 3 e 10, replicati 1-5, H0 e H1): fewclust Mancl e DeRouen, Kauermann e Carroll (radice principale) e Fay e Graubard uguali al calcolo a mano dagli articoli[^md2001][^kc2001][^fg2001] entro 1e-8 relativo dove entrambi rispondono, con al piu' l'1% di dataset in cui risponde uno solo, nelle tre correlazioni; in M2 e M5 la somma dei contributi effettivi allo score coincide con la derivata numerica della log-verosimiglianza (numDeriv) entro 1e-6 relativo, su almeno 50 dataset |
| GE5 | coincidenza geesmv / fewclust | con correlazione indipendente (stesso fit nei due pacchetti), replicati 1-50: Mancl e DeRouen e Fay e Graubard di geesmv e di fewclust uguali entro 1e-8 relativo dove entrambi sono finiti. Kauermann e Carroll e' escluso dal cancello per l'anomalia descritta sopra; con exchangeable la coincidenza non e' esatta per costruzione (geesmv riadatta con gee) e se ne riporta la distribuzione |

## Previsioni

Margini della skill salvo dove motivato: errore di I tipo da 2,5% a 7,5%; differenze fra metodi +/- 1 punto in media sui punti di disegno e +/- 2,5 punti per punto. Aggregato = media sui 1.896 punti sotto H0 (o sui punti indicati), intervallo dal bootstrap per blocchi di replicato (B = 2000, seed 20261015, indici di `results_v4/S4/boot_idx.rds` e `results_v4/S6b/ns/boot_idx_67.rds`). Le bandiere "chiaramente sopra / sotto" contengono gia' la MCSE e si contano come interi. Le previsioni sul test a inversione di segno valgono per la versione principale (senza filtro, [Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](../decisions/D-47.md)), che nasce dopo il pilota e si giudica qui per la prima volta su dati nuovi per lei (1.878 punti su 1.896).

| codice | quantita' | atteso | intervallo di rilevanza | perche' conta |
|---|---|---|---|---|
| PE1 | errore di I tipo del test a inversione di segno, media sui 1.896 punti | 4,5% (pilota: 4,15%) | da 2,5% a 7,5% | entra nella regola solo se tarato |
| PE2 | punti in cui il test a inversione di segno e' chiaramente sopra 7,5% | 0 | 0 | un punto basta a escluderlo dalla regola della Figura 1 |
| PE3 | punti in cui il test a inversione di segno e' chiaramente sotto 2,5% | 20 (con 10-15 soggetti e 2 visite) | da 0 a 60 | oltre 60 (meta' della correzione di Fay e Graubard con la t, 116) cambierebbe la frase "prudente solo dove i dati non hanno quasi informazione" |
| PE4 | quota usabile del test a inversione di segno meno quella di Fay e Graubard con la t dello studio, per punto | >= 0 in ogni punto | nessun punto sotto -2,5 | il test dovrebbe rispondere almeno quanto le strategie del paper |
| PE5 | quota usabile del test a inversione di segno, media sui 176 punti con 2-3 visite e fino a 25 soggetti | 98% (Fay e Graubard 43,6%, modello marginalizzato 40,4%) | da 90% a 100% | sotto 90% resterebbe la frase "nessun metodo e' insieme disponibile e tarato" |
| PE6 | errore di I tipo di Fay e Graubard exchangeable con la t: fewclust meno lo studio (geesmv), media sui punti H0 | 0 punti | da -1 a +1 | se l'implementazione spostasse il tasso, i numeri del paper su Fay e Graubard andrebbero rifatti |
| PE7 | punti in cui la bandiera "chiaramente sopra 7,5%" di Fay e Graubard exchangeable con la t differisce fra fewclust e studio | 0 | da 0 a 2 | 3 o piu' cambierebbero la frase del paper (oggi 2 punti) |
| PE8 | quota usabile di Fay e Graubard exchangeable: fewclust meno studio, media sui punti | +1 punto (sonda: +2,8 nei punti piccoli) | da 0 a +2 | misura quanto la disponibilita' delle correzioni del paper dipende da geesmv; per punto mai sotto -2,5 |
| PE9 | errore di I tipo di Fay e Graubard con la t, indipendente meno exchangeable (fewclust, appaiato), media sui punti H0 | 0 punti | da -1 a +1 | se la correlazione di lavoro contasse, la limitazione andrebbe scritta come risultato |
| PE10 | come PE9, AR(1) meno exchangeable | 0 punti | da -1 a +1 | idem |
| PE11 | errore di I tipo di Kauermann e Carroll corretto (fewclust, exchangeable, t), media sui punti H0 | 5,8% (fra il sandwich non corretto, 6,94%, e Mancl e DeRouen, 4,78%) | da 4,8% a 6,8% | lo studio riporta 9,75% con geesmv: un valore in questo intervallo cambia la frase del paper su Kauermann e Carroll |
| PE12 | strategie di Wald nuove mai chiaramente sopra 7,5% in nessuno dei 1.896 punti | 0 | 0 | una strategia nuova mai sopra e con quota usabile oltre il 91,0% del modello marginalizzato cambierebbe la raccomandazione della Figura 1 |

Prima di fissare i margini: con 1.000 replicati per braccio la MCSE tipica dell'errore di I tipo per punto e' circa 0,7 punti; l'intervallo per blocchi della media su 1.896 punti e' atteso sotto +/- 0,2 punti, quindi i margini di +/- 1 punto non rendono le previsioni indeterminate per costruzione.

[^kc2001]: Kauermann, Carroll (2001), JASA 96:1387-1396.
[^md2001]: Mancl, DeRouen (2001), Biometrics 57:126-134.
[^fg2001]: Fay, Graubard (2001), Biometrics 57:1198-1206.

# Esito del calcolo (2026-10-23)

Run: `R/s7d/s7d_run.R`, 37.920 chunk (3.792.000 dataset), 20-28 processi regolati dalla temperatura (Tctl fra 75 e 86 gradi), 11,5 ore, 0 chunk falliti; pre-registrazione `47c032a`, script di aggregazione committato prima del run (`afb66a5`). Controlli: `R/s7d/s7d_gate_formulas.R`, `s7d_gate_reprod.R`, `s7d_ge5_detail.R`; aggregazione `s7d_agg.R`. Tabelle `tables/V4_S7d_*.csv`.

## Cancelli

| codice | esito | dettaglio |
|---|---|---|
| GE1 | passa | md5 dei 28.800 file di `results_full/` prima = dopo = lista di S7b; parita' H0/H1 passa; 37.920 chunk, 3.792.000 replicati unici, 580 errori di generazione tutti in disegni non stimabili, tre righe di GEE per ogni dataset generato |
| GE2 | passa | GEE exchangeable: 3.744.885 righe salvate finite, stima ed errore standard identici (scarto massimo 0); geesmv, replicati 1-50: 175.184-175.259 errori standard per correzione identici (scarto 0); in 1-3 casi per correzione geesmv ha risposto nello studio e non qui (causa non verificata; plausibilmente il limite di 5 s sotto carico) |
| GE3 | passa | 36.000 replicati del pilota, 0 differenze per la versione principale e per quella con filtro |
| GE4 | fallito come scritto, spiegato | (a) passa: fewclust coincide con il calcolo a mano entro 8,7e-13 relativo, al piu' 4 dataset su circa 670 con una sola risposta, nelle tre correlazioni (`tables/V4_S7d_GE4.csv`). (b) in 31 dataset su 240 di M2 e M5 lo score e la derivata numerica sono entrambi sotto 1e-10 in valore assoluto (modello ridotto quasi separato con 10 soggetti e 3 visite) e lo scarto relativo non ha senso; lo scarto assoluto massimo su tutti i 240 e' 1,0e-9, contro valori tipici fra 0,4 e 36. La formula e' giusta; la tolleranza relativa era mal posta vicino a zero |
| GE5 | fallito come scritto, spiegato | con indipendenza, 19 dataset su circa 175.000 hanno scarti grandi fra geesmv e fewclust: tutti con fit divergente (stima di `time:trt` fra 35 e 4e15 in valore assoluto), esclusi dal filtro di usabilita' dello studio; i due pacchetti si fermano in punti diversi di un'ottimizzazione che diverge. Fra i fit usabili lo scarto massimo e' 5,0e-10 per Mancl e DeRouen e 1,2e-12 per Fay e Graubard (`tables/V4_S7d_GE5_mismatches.csv`) |

## Previsioni

| codice | atteso | intervallo | osservato (intervallo per blocchi) | esito |
|---|---|---|---|---|
| PE1 | 4,5% | 2,5-7,5% | 5,29% (5,25-5,32) | confermata |
| PE2 | 0 | 0 | 0 | confermata |
| PE3 | 20 | 0-60 | 4 | confermata |
| PE4 | >= 0 | non sotto -2,5 | minimo 0 | confermata |
| PE5 | 98% | 90-100% | 99,9% | confermata |
| PE6 | 0 | -1/+1 | +0,02 (0,01-0,03) | confermata |
| PE7 | 0 | 0-2 | 0 | confermata |
| PE8 | +1 | 0/+2 | +0,42 (0,41-0,43); minimo per punto -0,20 | confermata |
| PE9 | 0 | -1/+1 | -0,29 (-0,30/-0,27) | confermata |
| PE10 | 0 | -1/+1 | -0,31 (-0,33/-0,29) | confermata |
| PE11 | 5,8% | 4,8-6,8% | 5,55% (5,52-5,59) | confermata |
| PE12 | 0 | 0 | 1 | smentita |

PE12 smentita, spiegazione: la GEE con correlazione AR(1), Fay e Graubard di fewclust e la t a N - p non e' chiaramente sopra 7,5% in nessuno dei 1.896 punti (il valore piu' alto: 8,38% con 25 soggetti, 4 visite, scenario con pendenza 0,25 della variabilita' latente, limite inferiore 6,54) e risponde nel 92,2% dei dataset contro il 91,0% del modello marginalizzato con controllo e t. Con la regola come scritta ([La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](../decisions/D-35.md)) sceglierebbe questa GEE. Ha meno potenza (19,9% contro 21,2% per tentativo) ed e' chiaramente sotto 2,5% in piu' punti (112 contro 77). Con la t a N - 2 resta a 0 punti sopra. La GEE con indipendenza e Fay e Graubard ha 1 punto sopra, con exchangeable 2 (come nello studio). Se il test a inversione di segno entra nella regola, la regola sceglie lui: 0 punti sopra, risponde nel 99,98% dei dataset (99,9% nei 176 punti con 2-3 visite e fino a 25 soggetti), potenza 21,6%.

# In breve

## In breve: le strategie aggiunte, provate su tutte le 1.896 situazioni, cambiano la raccomandazione del paper?

**A cosa serviva.** Il paper dichiara tre limiti: la GEE è provata con un solo tipo di correlazione di lavoro (exchangeable); le correzioni dell'errore standard vengono da un solo pacchetto (geesmv); con 2-3 visite e fino a 25 soggetti nessuna strategia è insieme disponibile e tarata. Abbiamo rigenerato tutti i 3,8 milioni di dataset dello studio con gli stessi semi. Su ciascuno abbiamo stimato la GEE anche con correlazione indipendente e AR(1), con le correzioni di Mancl e DeRouen, Kauermann e Carroll e Fay e Graubard calcolate da un secondo pacchetto (fewclust), e abbiamo applicato il test a inversione di segno sui contributi dei soggetti allo score scelto nella prova precedente. I dati rigenerati coincidono esattamente con quelli dello studio.

**Cosa abbiamo trovato.**

1. **Il test a inversione di segno è tarato ovunque e risponde quasi sempre, anche dove oggi nessun metodo risponde.**
   - In media 5,3 falsi positivi su 100; in nessuna delle 1.896 situazioni è fuori taratura verso l'alto; è prudente (sotto 2,5) solo in 4 situazioni, tutte con 10 soggetti e 2-3 visite.
   - Con 2-3 visite e fino a 25 soggetti dà un risultato nel 99,9% dei dataset; la correzione di Fay e Graubard nel 44% e il modello marginalizzato nel 40%.
   - Trova un effetto vero in 21,6 tentativi su 100, come il modello marginalizzato (21,2) e più di Fay e Graubard (20,5).
   - Dà solo un p-value: nessuna stima e nessun intervallo.
2. **Una strategia nuova supera la regola della Figura 1: la GEE con correlazione AR(1), la correzione di Fay e Graubard e la t.**
   - Non è mai chiaramente sopra 7,5 falsi positivi su 100 e risponde nel 92,2% dei dataset, contro il 91,0% del modello marginalizzato. La regola, così come l'abbiamo scritta, sceglie quella che risponde più spesso, quindi sceglierebbe questa GEE.
   - Il vantaggio è piccolo: poco più di un dataset su 100. Ha meno potenza (19,9 contro 21,2 tentativi su 100) ed è prudente in più situazioni (112 contro 77).
   - Le stesse correzioni con correlazione exchangeable o indipendente restano sopra in 1-2 situazioni, con 15-25 soggetti e variabilità latente che cresce molto nel tempo.
3. **La correlazione di lavoro conta poco, e le correzioni di Mancl e DeRouen e di Fay e Graubard del paper sono giuste.**
   - Passare da exchangeable a indipendente o AR(1) abbassa i falsi positivi di Fay e Graubard di circa 0,3 su 100 in media.
   - Calcolata con fewclust invece che con geesmv, la correzione di Fay e Graubard dà gli stessi falsi positivi (differenza 0,02 su 100) e risponde in 4 dataset su 1.000 in più; nessuna situazione passa sopra o sotto la soglia di 7,5 falsi positivi su 100.
4. **La correzione di Kauermann e Carroll del paper è sbagliata nel pacchetto geesmv.**
   - geesmv calcola la radice di una matrice non simmetrica come se fosse simmetrica. Un calcolo a mano dalla formula dell'articolo coincide con fewclust e non con geesmv.
   - Con la versione corretta i falsi positivi medi scendono da 9,8 a 5,6 su 100, e le situazioni fuori taratura da 951 a 133. Resta esclusa dalla regola, ma la frase del paper su Kauermann e Carroll va rifatta.
5. **Con la correlazione AR(1) geesmv non funziona in questi dati.** Il pacchetto che usa per riadattare il modello rifiuta la correlazione AR(1) appena un soggetto ha una sola visita: nella prova preliminare su 990 dataset ha risposto in circa 1 su 8. La GEE AR(1) del punto 2 è calcolata con fewclust.

**Cosa devi decidere.**

1. **Che cosa raccomanda la Figura 1? Consiglio di seguire la regola come l'avevamo scritta: per il test, il test a inversione di segno; per stima e intervallo, la GEE con correlazione AR(1), Fay e Graubard e la t, con il modello marginalizzato come alternativa equivalente.**
   - È la regola fissata prima di vedere questi risultati; cambiarla adesso per tenere il modello marginalizzato sarebbe una scelta fatta a posteriori.
   - Il test a inversione di segno è l'unico che chiude il limite sulle 2-3 visite con pochi soggetti.
   - Conseguenza: il messaggio principale del paper passa dal modello marginalizzato a una combinazione di test e stima; titolo, abstract e discussione vanno riscritti.
   - Alternativa: tenere il modello marginalizzato come raccomandazione per la stima e scrivere che la GEE AR(1) con Fay e Graubard lo eguaglia (differenza di un dataset su 100 nella quota di risposte, a fronte di meno potenza); il test a inversione di segno come raccomandazione per il test. Il paper cambia meno, ma va dichiarato che la regola è stata letta dopo i risultati.
   - Alternativa: lasciare la Figura 1 com'è e mettere le strategie nuove nel supplemento. È la più semplice, ma il paper raccomanderebbe una strategia che la sua stessa regola non sceglie più.
2. **Come riportiamo la correzione di Kauermann e Carroll? Consiglio di sostituire nel paper la versione di geesmv con quella corretta di fewclust e di spiegare l'errore di geesmv in una nota del supplemento.**
   - Alternativa: tenere entrambe, con la versione di geesmv dichiarata come "implementata in geesmv". Più trasparente sul confronto con altri studi che usano geesmv, ma appesantisce tabelle e figure.

**Comunque:**

- quando riscriveremo il paper, la frase sulle 2-3 visite con pochi soggetti e quelle sulla correlazione di lavoro e sulle correzioni useranno questi numeri;
- citeremo fewclust nella versione registrata per questo calcolo, o in una successiva se la aggiorni;
- i due controlli falliti alla lettera (una tolleranza relativa vicino a zero e 19 fit divergenti su 175.000) sono spiegati nel report e non toccano nessun risultato.

Appena decidi questi 2 punti, lo step è chiuso.

Il PI ha deciso i 2 punti il 2026-10-23 (il primo nella versione rivista dopo la sua domanda, sotto): step chiuso.

**Sui tempi di calcolo.** Il run è durato 11,5 ore sul server con 20-28 processi regolati dalla temperatura, contro le circa 16 ore stimate prima del run con 24 processi.

# Scelte del PI sul riassunto (2026-10-23)

- **Scelta 1, Figura 1.** Dopo il riassunto il PI ha osservato che il test a inversione di segno non da' una stima interpretabile clinicamente. L'agente ha rivisto la raccomandazione: la regola si applica alle sole strategie con stima e intervallo. Il PI l'ha accettata, con il test a inversione di segno breve nel testo e in dettaglio nel supplemento. Vedi [La regola della Figura 1 si applica alle strategie con stima e intervallo: sceglie la GEE AR(1) con Fay e Graubard e la t; il test a inversione di segno e' complementare (D-50)](../decisions/D-50.md), che supera [La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](../decisions/D-35.md).
- **Scelta 2, Kauermann e Carroll: opzione A.** Vedi [La correzione di Kauermann e Carroll nel paper viene da fewclust; l'errore di geesmv si spiega nel supplemento (D-51)](../decisions/D-51.md).
- **Proposta aggiunta: intervallo dal test.** Il PI ha scelto uno step dedicato prima della riscrittura. Vedi [Uno step dedicato prova l'intervallo di confidenza ottenuto invertendo il test a inversione di segno, prima della riscrittura del paper (D-52)](../decisions/D-52.md) e [Intervallo di confidenza dal test a inversione di segno (S7d2)](S7d2.md).

# Dubbi rinviati

Rinviato all'[Aggiornamento del manoscritto con le analisi aggiunte (S7f)](S7f.md), che lo riceve nel suo file: il testo da scrivere sui risultati di questo step.
