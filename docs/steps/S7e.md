---
type: Step
title: "Terzo meccanismo generatore"
description: "Un meccanismo che non e' il modello di casa di nessun metodo, simulato su tutta la griglia principale con tutte le strategie."
step_id: "S7e"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-10-24T22:05:00+02:00, evidence: "chat del 2026-10-24: \"ok\" sulla raccomandazione del riassunto In breve (D-55); commit S7e: closed, approved by PI" }
resource: "../../reports/morris/S7e.html"
generated: { by: mmm-agent/skill-v4, at: 2026-10-23T10:40:00+02:00 }
sources:
  - id: d42
    resource: "../decisions/D-42.md"
    title: "Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale"
  - id: d53
    resource: "../decisions/D-53.md"
    title: "Il quarto meccanismo e' una copula gaussiana con margini logistici e correlazione solo seriale (rho = 0,6)"
  - id: d54
    resource: "../decisions/D-54.md"
    title: "Organizzazione tecnica del run del quarto meccanismo: semi dello studio, correzioni da fewclust, geesmv sui primi 50 replicati"
  - id: d55
    resource: "../decisions/D-55.md"
    title: "Il quarto meccanismo entra nel testo principale: regola della Figura 1 sui quattro meccanismi, frase sul modello misto, dettagli nel supplemento"
  - id: s7d
    resource: "S7d.md"
    title: "Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto"
  - id: song2009
    resource: https://doi.org/10.1111/j.1541-0420.2008.01058.x
    title: "Song, Li, Yuan (2009). Joint regression analysis of correlated data using Gaussian copulas. Biometrics 65:60-68"
  - id: mv2012
    resource: https://doi.org/10.1214/12-EJS721
    title: "Masarotto, Varin (2012). Gaussian copula marginal regression. Electronic Journal of Statistics 6:1517-1549"
  - id: ep1991
    resource: https://doi.org/10.1080/00031305.1991.10475828
    title: "Emrich, Piedmonte (1991). A method for generating high-dimensional multivariate binary variates. The American Statistician 45:302-304"
---
# Dubbi ricevuti

1. **I due meccanismi distinti sono ciascuno il modello di casa di un metodo.** Ricevuto il 2026-10-21 dal piano approvato dal PI ([Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](../decisions/D-42.md)). Il meccanismo con il processo autoregressivo latente si e' comportato quasi come il modello misto, quindi lo studio poggia su due meccanismi, ciascuno di casa per un metodo. Il piano iniziale apriva questo step solo se il run con i test a livello di soggetto spostava le conclusioni; il PI ha poi indicato che i tempi non sono un vincolo, e la limitazione resta anche se le conclusioni non cambiano: lo step si fa comunque.

# Piano

- Meccanismo candidato: dati binari correlati generati con una copula sui margini del modello marginale, con correlazione fra visite che decresce con la distanza; da scegliere con la letteratura e calibrare sui coefficienti marginali bersaglio di M1-M5, con un valore vero marginale esatto per costruzione. Parametri e verifica della calibrazione nel file dello step prima del run.
- Griglia: tutta la griglia principale (16 valori di N, 6 di T, M1-M5: 480 punti di disegno), H0 e H1, 1.000 replicati per braccio, semi con uno schema nuovo e dichiarato.
- Metodi: tutte le strategie del paper (compreso il modello marginalizzato con il controllo di convergenza) e quelle aggiunte nel run con i test a livello di soggetto.
- Prodotti: `results_v4/S7e/`, `tables/V4_S7e_*.csv`, regola della Figura 1 riapplicata su tutti i punti di disegno.

# Scelta del meccanismo (2026-10-23, prima di ogni analisi su dati del nuovo meccanismo)

Nessun metodo e' stato applicato ai dati del nuovo meccanismo prima di questa scelta, salvo una sonda dei tempi (`R/s7e/s7e_timing.R`, 8 replicati in 7 celle, H0) di cui non si e' calcolata ne' letta alcuna misura di prestazione.

**Profilo di dipendenza dei meccanismi** (`R/s7e/s7e_profiles.R`; un dataset di 20.000 soggetti per meccanismo, M1, H0, 10 visite, seme 20261023; correlazione di Pearson dei residui standardizzati fra visite dello stesso soggetto; `tables/V4_S7e_dependence_summary.csv`):

| meccanismo | visite consecutive | distanza 3 | distanza 9 | media su 10 visite |
|---|---|---|---|---|
| GenBinaryY | 0,312 | 0,163 | 0,136 | 0,194 |
| GLMM-native | 0,132 | 0,133 | 0,129 | 0,134 |
| GLMM-AR1 | 0,150 | 0,135 | 0,087 | 0,131 |
| copula, solo seriale, rho 0,6 (scelta) | 0,361 | 0,111 | -0,002 | 0,141 |

I tre meccanismi dello studio hanno una componente costante per soggetto; la copula scelta no. Scelta del PI: [Il quarto meccanismo e' una copula gaussiana con margini logistici e correlazione solo seriale (rho = 0,6) (D-53)](../decisions/D-53.md). Organizzazione tecnica: [Organizzazione tecnica del run del quarto meccanismo: semi dello studio, correzioni da fewclust, geesmv sui primi 50 replicati (D-54)](../decisions/D-54.md), che sostituisce la riga del piano sui semi.

# Disegno del run (fissato prima del calcolo)

- **Meccanismo** (etichetta `Copula-AR1`): P(Y_ij = 1) = expit(X_ij' beta) con il modello medio e i coefficienti marginali dello studio (M1: -2, 0,1, 0,2, 0,2; M2: con age 0,01; M3 come M1; M4: -2, 0,2, 0,2, 0,05; M5: con dose 0,05); sotto H0 cambia solo `time:trt`, che vale 0. Y_ij = 1{Z_ij <= qnorm(expit(X_ij' beta))}, Z_i ~ N(0, R), R[j,k] = 0,6^|t_j - t_k|[^ep1991][^song2009][^mv2012]. Valore vero marginale = beta esattamente, su ogni disegno realizzato. Il modello misto non ha un valore vero nel generatore: distorsione, copertura e errore quadratico medio del modello misto non sono definiti, come per GenBinaryY ([Nessun valore vero condizionale per GenBinaryY (D-18)](../decisions/D-18.md)); il suo errore di I tipo si riporta.
- **Griglia**: M1-M5 x N in {10, 15, 20, 25, 35, ..., 145} x T in {2, 3, 4, 6, 8, 10} = 480 punti, H0 e H1, 1.000 replicati per braccio, 9.600 chunk di 100 replicati, 960.000 dataset.
- **Strategie per replicato** (`R/s7e/s7e_common.R`): modello marginalizzato (studio), modello marginalizzato con controllo di convergenza in M3-M5, GEE exchangeable di geepack, modello misto con Wald (opzioni di default) e riscalamento shrinkage, modello misto con LRT (bobyqa, pendenze casuali del Wald); GEE con correlazione exchangeable, indipendente e AR(1), ciascuna con sandwich non corretto e con Mancl e DeRouen, Kauermann e Carroll e Fay e Graubard di fewclust, CR2 di clubSandwich su exchangeable e indipendente, geesmv sui replicati 1-50; test a inversione di segno principale e con filtro.
- **Misure**: come nel protocollo (sez. 5), usabilita' U0, denominatore senza i disegni non stimabili, riferimento t a N - p per la regola della Figura 1 e normale per le tabelle, bandiere "chiaramente sopra 7,5%" = tasso - 2 MCSE > 7,5 con almeno 20 replicati usabili.
- **Confronto con gli altri meccanismi**: per ogni punto, tasso della copula meno la media dei tassi dei tre meccanismi dello studio nello stesso punto (stessi disegni, stessi semi); intervallo dal bootstrap per blocchi di replicato con gli indici di `results_v4/S4/boot_idx.rds` (B = 2000, seed 20261015), gli stessi per i quattro meccanismi.
- **Regola della Figura 1** ([La regola della Figura 1 si applica alle strategie con stima e intervallo: sceglie la GEE AR(1) con Fay e Graubard e la t; il test a inversione di segno e' complementare (D-50)](../decisions/D-50.md)): sulle strategie con stima e intervallo, sui 480 punti della copula da soli e sui 2.376 punti (1.896 + 480) insieme.

# Attese pre-registrate (scritte prima di ogni calcolo del run)

## Cancelli

| codice | controllo | passa se |
|---|---|---|
| GF1 | integrita' e completezza | md5 di `results_full/` prima e dopo il run uguali fra loro e alla lista di S7b; `R/testing/test_h0h1_parity.R` passa prima del run; 9.600 chunk, 960.000 replicati unici; ogni errore di generazione cade in un disegno non stimabile; ogni dataset generato ha le righe delle quattro tecniche dello studio e tre righe di GEE |
| GF2 | parita' H0/H1 e stesso disegno degli altri meccanismi | su 100 replicati per scenario a (N 10, T 2) e (N 25, T 6), prima del run: H0 e H1 hanno colonne del disegno identiche e gli stessi numeri normali (controllo: dove expit(X beta) e' piu' alto sotto H1, Y sotto H1 >= Y sotto H0, e viceversa, in ogni osservazione); le colonne del disegno (id, time, trt, age, dose) coincidono con quelle di GLMM-native e GLMM-AR1 allo stesso seme |
| GF3 | calibrazione del generatore, controprova indipendente delle formule | prima del run, un dataset di 100.000 soggetti per scenario (H1, 10 visite, seme 20261023): (a) in M1, M3 e M4 la prevalenza osservata in ogni cella (visita x trattamento) e, per ogni coppia di visite e braccio, la frequenza di Y = 1 in entrambe coincidono con expit(X beta) e con la probabilita' normale bivariata calcolata per integrazione numerica, con abs(z) < 4 in tutte le celle; (b) in M2 e M5 la GEE di indipendenza recupera ogni coefficiente marginale entro 4 errori standard robusti |
| GF4 | riproduzione e coincidenza fra percorsi di codice | (a) 3 chunk rigenerati dopo il run in un processo nuovo sono identici a quelli del run (stime, errori standard, p-value; esclusi i tempi); (b) in ogni replicato la stima della GEE exchangeable dello studio coincide con quella della GEE exchangeable aggiunta entro 1e-10 assoluto, e in M3-M5 il fit A del controllo di convergenza coincide con il modello marginalizzato dello studio entro 1e-10; (c) replicati 1-50, correlazione indipendente, fit usabili: Mancl e DeRouen e Fay e Graubard di geesmv e fewclust coincidono entro 1e-6 relativo o 1e-10 assoluto |

## Previsioni

Margini della skill: errore di I tipo da 2,5% a 7,5%; copertura da 92,5% a 97,5%; distorsione / SE empirico da -0,46 a +0,46; differenze fra meccanismi o metodi +/- 1 punto in media sui punti di disegno; quota usabile +/- 2,5 punti in media (dove una differenza di quota cambia la scelta della regola). Aggregato = media sui 480 punti (H0 per l'errore di I tipo e la quota usabile, H1 per distorsione e copertura); intervallo dal bootstrap per blocchi. Le bandiere si contano come interi. "Altri meccanismi" = media dei tre meccanismi dello studio nello stesso punto.

| codice | quantita' | atteso | intervallo di rilevanza | perche' conta |
|---|---|---|---|---|
| PF1 | punti in cui il test a inversione di segno e' chiaramente sopra 7,5% | 0 | 0 | il test e' esatto se i contributi dei soggetti sono simmetrici; un punto sopra toglierebbe la frase "tarato in tutti i meccanismi" |
| PF2 | errore di I tipo del test a inversione di segno, copula meno altri meccanismi | 0 punti (altri: 5,24-5,30%) | da -1 a +1 | idem, in media |
| PF3 | punti in cui la GEE AR(1) con Fay e Graubard e la t e' chiaramente sopra 7,5% | 0 | 0 | e' la strategia della Figura 1: un punto sopra la escluderebbe dalla regola |
| PF4 | punti in cui il modello marginalizzato con controllo di convergenza e la t e' chiaramente sopra 7,5% | 0 | 0 | e' l'alternativa equivalente della Figura 1; il suo modello di dipendenza qui e' sbagliato e gli errori standard sono basati sul modello |
| PF5 | errore di I tipo della GEE AR(1) con Fay e Graubard e la t, copula meno altri | 0 punti (altri: 4,23-4,50%) | da -1 a +1 | se il meccanismo spostasse la taratura, la raccomandazione andrebbe condizionata al tipo di dipendenza |
| PF6 | errore di I tipo del modello marginalizzato con controllo di convergenza e la t, copula meno altri | +0,5 punti (altri: 4,28-4,57%) | da -1 a +1 | idem; atteso un piccolo aumento perche' la transizione di primo ordine cattura solo in parte la dipendenza seriale |
| PF7 | distorsione / SE empirico del modello marginalizzato con controllo di convergenza, media sui punti H1 | 0 | da -0,46 a +0,46 | il modello medio e' giusto ma il modello di dipendenza no: una distorsione rilevante toglierebbe la frase che il modello marginalizzato stima l'effetto marginale anche fuori dal suo meccanismo |
| PF8 | copertura al 95% del modello marginalizzato con controllo di convergenza (riferimento normale), media sui punti H1 | 94% | da 92,5% a 97,5% | errori standard basati su un modello di dipendenza sbagliato potrebbero sottostimare la variabilita' |
| PF9 | errore di I tipo del modello misto con Wald e la t, copula meno altri | +1,5 punti (altri: 5,65-8,02%) | da +1 a +10 | il modello misto ignora la dipendenza seriale, quindi l'errore standard di un effetto che varia nel soggetto dovrebbe essere troppo piccolo; sotto +1 la frase sul modello misto non cambierebbe |
| PF10 | strategia scelta dalla regola della Figura 1 sui 2.376 punti (1.896 + 480 della copula) | la GEE AR(1) con Fay e Graubard e la t | stessa strategia | se cambiasse, cambierebbe la raccomandazione del paper |
| PF11 | quota usabile del modello marginalizzato con controllo di convergenza, copula meno altri, media sui punti H0 | 0 punti | da -2,5 a +2,5 | la regola sceglie per quota usabile; uno scarto oltre 2,5 punti cambierebbe il confronto con la GEE AR(1) |

Prima di fissare i margini: la MCSE dell'errore di I tipo per punto e' circa 0,7 punti; nel run precedente l'intervallo per blocchi delle medie sui punti era largo +/- 0,04 punti per un meccanismo, quindi per una differenza fra meccanismi su 480 punti e' atteso sotto +/- 0,2: i margini di +/- 1 punto non rendono le previsioni indeterminate per costruzione. La MCSE di distorsione / SE empirico per punto e' circa 0,03, molto sotto 0,46.

[^ep1991]: Emrich, Piedmonte (1991), The American Statistician 45:302-304.
[^song2009]: Song, Li, Yuan (2009), Biometrics 65:60-68.
[^mv2012]: Masarotto, Varin (2012), Electronic Journal of Statistics 6:1517-1549.

# Controlli prima del run (2026-10-23)

- **Integrita' e parita'** (parte di GF1): md5 dei 28.800 file di `results_full/` uguali alla lista di S7b (`results_v4/S7e/md5_results_full_before.txt`, 0 differenze); `R/testing/test_h0h1_parity.R` passa (griglia principale 15/15, M1c 3/3, controllo di [Varianza della pendenza su dose allineata in M5 (D-02)](../decisions/D-02.md) 8/8).
- **GF2 passa** (`R/s7e/s7e_gates_pre.R`, `tables/V4_S7e_GF2.csv`): 1.000 replicati (5 scenari x 2 disegni x 100); disegno identico fra H0 e H1 e stessi numeri normali in tutti; disegno identico a GLMM-native in tutti e a GLMM-AR1 in tutti quelli in cui GLMM-AR1 genera (10 errori di generazione di GLMM-AR1, tutti in disegni non stimabili, a N = 10 e T = 2).
- **GF3 passa** (`tables/V4_S7e_GF3.csv`): in M1, M3 e M4 la massima abs(z) e' 2,7 su 60 celle marginali e 270 celle di coppia; in M2 e M5 la GEE di indipendenza recupera ogni coefficiente entro 1,5 errori standard robusti.
- **Prova su 6 chunk** (`R/s7e/s7e_run.R 6`, 600 replicati, 0 chunk falliti): secondi per replicato da 1,4 (N = 10, T = 2) a 12,6 (M5, N = 145, T = 10), di cui 9,8 per le analisi dello studio (soprattutto il modello misto con tre effetti casuali e il suo LRT) e 1,8 per il controllo di convergenza del modello marginalizzato. Stima del run intero, fatta dall'agente da questi tempi e da quelli del run precedente per la parte GEE: circa 550 ore di processore, cioe' circa 23 ore con 24 processi. I 6 chunk restano come parte del run.

# Esito del calcolo (2026-10-24)

Run: `R/s7e/s7e_run.R`, 9.600 chunk (960.000 dataset), 20-25 processi regolati dalla temperatura (guardia `R/s7e/s7e_thermal_guard.sh`), 26,4 ore, 0 chunk falliti; pre-registrazione `9e4c8ae`, script di aggregazione committato prima del run (`8fbe6de`; dopo, due correzioni dello strumento descritte sotto). Aggregazione `R/s7e/s7e_agg.R`; riassunti e figure `R/s7e/s7e_summary.R`; controlli dopo il run `R/s7e/s7e_gate_reprod.R`, `R/s7e/s7e_explain.R`, `R/s7e/s7e_glmm_pseudotrue.R`, `R/s7e/s7e_glmm_se_ratio.R`, `R/s7e/s7e_glmm_singular.R`, `R/s7e/s7e_glmm_by_scen_T.R`. Tabelle `tables/V4_S7e_*.csv`.

**Correzioni dello strumento di aggregazione** (provate sui primi chunk, prima di leggere qualsiasi risultato aggregato): (1) il controllo di completezza contava come mancante una tecnica senza riga `time:trt`; si riportano ora i dataset mancanti e, fra questi, quelli con disegno stimabile; (2) i quantili del bootstrap ignorano i punti senza replicati usabili in una replica bootstrap (`na.rm`).

## Cancelli

| codice | esito | dettaglio |
|---|---|---|
| GF1 | fallito come scritto, spiegato | md5 dei 28.800 file di `results_full/` prima = dopo = lista di S7b; parita' H0/H1 passa; 9.600 chunk, 960.000 replicati unici, 0 errori di generazione, tre righe di GEE per ogni dataset. In 344 dataset su 960.000 manca la riga di `time:trt` del modello marginalizzato (253) o del modello misto (270): tutti in disegni non stimabili (un braccio mai osservato dopo il tempo 0), dove i fit restituiscono coefficienti senza `time:trt` e la riga e' tolta dal filtro sul termine; sono esclusi comunque dal denominatore ([Replicati con disegno non stimabile esclusi dal denominatore (D-11)](../decisions/D-11.md)). Nessun risultato cambia (`tables/V4_S7e_GF1.csv`, `V4_S7e_GF1_missing.csv`) |
| GF2 | passa | vedi i controlli prima del run |
| GF3 | passa | vedi i controlli prima del run |
| GF4 | passa | (a) 3 chunk rigenerati identici (`tables/V4_S7e_GF4a.csv`); (b) GEE exchangeable dello studio = GEE exchangeable aggiunta in 959.506 replicati, fit A del controllo di convergenza = modello marginalizzato dello studio in 559.769, scarto massimo 0; (c) geesmv = fewclust con indipendenza, replicati 1-50: Mancl e DeRouen 44.026 fit, Fay e Graubard 44.050, scarto massimo 2,5e-12 (`tables/V4_S7e_GF4.csv`) |

## Previsioni

Valori in `tables/V4_S7e_predictions.csv` (intervallo = percentili 2,5 e 97,5 del bootstrap per blocchi; per i conteggi il valore intero).

| codice | atteso | intervallo | osservato (intervallo per blocchi) | esito |
|---|---|---|---|---|
| PF1 | 0 | 0 | 0 | confermata |
| PF2 | 0 | da -1 a +1 | -0,04 (-0,11 a +0,03) | confermata |
| PF3 | 0 | 0 | 0 | confermata |
| PF4 | 0 | 0 | 0 | confermata |
| PF5 | 0 | da -1 a +1 | +0,23 (+0,16 a +0,30) | confermata |
| PF6 | +0,5 | da -1 a +1 | +0,71 (+0,64 a +0,79) | confermata |
| PF7 | 0 | da -0,46 a +0,46 | 0,027 (0,023 a 0,030); nessun punto fuori | confermata |
| PF8 | 94% | da 92,5% a 97,5% | 93,8% (93,7 a 93,9) | confermata |
| PF9 | +1,5 | da +1 a +10 | +3,2 (+3,0 a +3,5), 418 punti | confermata |
| PF10 | GEE AR(1) con Fay e Graubard e t | stessa | GEE AR(1) con Fay e Graubard e t (92,1% di risposte; modello marginalizzato con controllo 90,9%) | confermata |
| PF11 | 0 | da -2,5 a +2,5 | -0,58 (-0,66 a -0,50) | confermata |

PF5, PF6 e PF11 sono differenze che si distinguono dal caso ma restano dentro la soglia di importanza: si scrivono come trascurabili. PF9 su 418 punti perche' in 62 punti il modello misto con pendenze casuali non ha replicati usabili in almeno un meccanismo (serie di 2-3 visite, come in [GLMM con pendenze casuali non stimabile nelle serie piu' corte (D-20)](../decisions/D-20.md)); limitandosi ai 414 punti con almeno 20 replicati usabili in tutti i meccanismi lo scarto resta +3,2 punti.

## Osservazioni non previste, chiuse

1. **Modello misto molto liberale con la sola intercetta casuale.** Negli scenari M1 e M2 l'errore di I tipo del Wald cresce con N fino a 19,5% (con N >= 85 fra 8,9% e 19,5%). Non e' un valore limite diverso da zero: su 20.000 soggetti sotto H0 il modello misto converge a 0,005 (SE 0,005) con 10 visite e 0,036 (SE 0,024) con 4 (`tables/V4_S7e_glmm_pseudotrue.csv`), e le stime sono centrate (in M1 e M2, con N >= 85, distorsione / SE empirico sotto 0,04 in media per scenario e T e al massimo 0,09 per punto). E' l'errore standard: con N >= 85 l'errore standard del modello e' fra 0,64 e 0,81 volte la dispersione delle stime, che con la normale implica fra 11% e 21% di rifiuti, come osservato (`tables/V4_S7e_glmm_se_ratio.csv`). L'intercetta casuale assume la stessa correlazione fra tutte le visite e ignora quella, piu' forte, fra visite vicine, che pesa proprio sull'effetto che varia nel soggetto. Con le pendenze casuali (M3-M5) e 8-10 visite il rapporto e' circa 1 e l'errore torna fra 5,3% e 6,2%; con 3-6 visite resta fra 8% e 16% (`tables/V4_S7e_glmm_by_scen_T.csv`), lo stesso andamento degli altri meccanismi ma piu' marcato (negli altri: 4-15%). Spiegato.
2. **GEE exchangeable con Fay e Graubard chiaramente sopra 7,5% in 5 punti** (0 negli altri meccanismi): M1 con N 10-25 e M3 con N 15, tutti con T = 4, tassi fra 9,6% e 10,4% (`R/s7e/s7e_explain.R`). Negli stessi punti la GEE AR(1) con Fay e Graubard e' fra 7,0% e 8,3% (mai chiaramente sopra) e la correlazione indipendente fra 7,4% e 9,7% (chiaramente sopra in 1 punto della griglia). In media sui punti con T = 4 la exchangeable e' a 6,1% contro 5,4% dell'AR(1). Con dipendenza solo seriale la correlazione exchangeable e' la piu' lontana dal vero; non tocca la Figura 1, che usa la AR(1). Spiegato; la frase del paper sulla GEE exchangeable va aggiornata.
3. **La correlazione AR(1) e' vicina alla dipendenza di questo meccanismo.** Era stato detto al PI prima della scelta. La regola della Figura 1 sceglieva la GEE AR(1) gia' sui tre meccanismi dello studio (92,2% di risposte), quindi la conclusione non poggia sul nuovo meccanismo. Spiegato.
4. **Meno potenza per tutte le strategie marginali** (17-19 tentativi su 100 contro 20-22 negli altri meccanismi): con correlazione seriale l'informazione sulla pendenza nel soggetto e' minore. Il modello misto resta a 23-24 perche' liberale. Atteso, non tocca conclusioni.
5. **Modello marginalizzato senza il controllo di convergenza**: chiaramente sopra 7,5% in 59 punti (con il controllo 0); conferma la scelta di [Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente (D-32)](../decisions/D-32.md).

# In breve

## In breve: le conclusioni del paper reggono con dati che non vengono dal modello di nessun metodo?

**A cosa serviva.** I tre modi in cui lo studio genera i dati hanno tutti una componente che resta costante per ogni soggetto: è proprio ciò che assumono il modello misto (qui il modello logistico a effetti casuali per soggetto), il modello marginalizzato multilivello e la GEE con correlazione exchangeable. Abbiamo aggiunto un quarto modo in cui la correlazione fra visite dello stesso soggetto si spegne con la distanza, senza parte costante, mentre l'effetto del trattamento nel tempo sulla scala della popolazione è esattamente quello voluto. Abbiamo simulato tutte le 480 situazioni della griglia principale, 1.000 simulazioni per scenario con e senza effetto, con tutte le strategie del paper e quelle aggiunte nel passaggio precedente.

**Cosa abbiamo trovato.**

1. **La raccomandazione della Figura 1 non cambia: resta la GEE con correlazione AR(1), la correzione di Fay e Graubard e la t, con il modello marginalizzato come alternativa equivalente.**
   - Con i dati nuovi nessuna delle due è mai fuori taratura verso l'alto: in media 4,6 e 5,1 falsi positivi su 100.
   - Applicata a tutte le 2.376 situazioni dei quattro modi di generare i dati, la regola sceglie la stessa GEE: risponde nel 92,1% dei dataset, il modello marginalizzato nel 90,9%.
   - Un'avvertenza: la correlazione AR(1) somiglia alla dipendenza dei dati nuovi, che le sono quindi un po' favorevoli; la regola però la sceglieva già con i tre modi dello studio.
   - Con questa dipendenza tutte le strategie trovano un effetto vero un po' meno spesso: 17-19 tentativi su 100 invece di 20-22.
2. **Il modello marginalizzato stima bene l'effetto anche quando la sua descrizione della dipendenza è sbagliata, purché si usi il controllo di convergenza.**
   - La distorsione è trascurabile e l'intervallo al 95% contiene il valore vero nel 94% dei casi (95% con la t).
   - Rispetto agli altri modi ha meno di un falso positivo su 100 in più: sotto la soglia di importanza.
   - Senza il controllo di convergenza è fuori taratura verso l'alto in 59 situazioni.
3. **Il modello misto peggiora molto: 10 falsi positivi su 100 in media, 3 in più che con gli altri modi, ed è fuori taratura verso l'alto in 219 situazioni su 480 (con gli altri modi da 38 a 100).**
   - Quando ha solo l'intercetta casuale, il suo errore standard ignora la correlazione più forte fra visite vicine ed è circa fra due terzi e quattro quinti di quello giusto; con molti soggetti si arriva a quasi 20 falsi positivi su 100. La stima invece è centrata.
   - Con le pendenze casuali e 8-10 visite torna vicino a 5 su 100.
4. **Il test a inversione di segno resta tarato: 5,2 falsi positivi su 100 in media, mai fuori taratura verso l'alto, e risponde in almeno 98 dataset su 100 in ogni situazione.**
5. **La GEE usata finora nel paper, con correlazione exchangeable e Fay e Graubard, esce dalla taratura in 5 situazioni, tutte con 4 visite e 10-25 soggetti, con circa 10 falsi positivi su 100.** Con gli altri modi non succedeva mai.

**Cosa devi decidere.**

1. **Come entra il nuovo modo di generare i dati nel paper? Consiglio di presentarlo nel testo principale come quarto meccanismo, con la regola della Figura 1 calcolata su tutti e quattro, una frase sul modello misto e i dettagli nel supplemento.**
   - Risponde direttamente al limite "ogni meccanismo è il modello di un metodo", e porta un risultato nuovo sul modello misto.
   - Alternativa: solo nel supplemento, citato in una frase come verifica di robustezza. Il testo principale cambia meno, ma il limite resta scritto lì e il risultato sul modello misto passa in secondo piano.
   - Alternativa: in tutte le tabelle e figure principali accanto agli altri tre. È più completo, ma appesantisce un testo già vicino al limite di parole.

**Comunque:**

- quando riscriveremo il paper, la frase sui meccanismi come modelli di casa, quelle sul modello misto e sulla GEE exchangeable, e la regola della Figura 1 useranno questi numeri;
- il controllo di completezza fallito alla lettera (344 dataset su 960.000 senza la riga dell'effetto per il modello marginalizzato o per il modello misto, tutti in disegni in cui l'effetto non è stimabile e che sono comunque esclusi) è spiegato nel report e non tocca nessun risultato.

Appena decidi questo punto, lo step è chiuso.

Il PI ha deciso il punto il 2026-10-24, con la raccomandazione: step chiuso.

**Sui tempi di calcolo.** Il run è durato 26,4 ore sul server, con 20-25 processi regolati dalla temperatura. Prima del run l'agente aveva stimato circa 23 ore con 24 processi, dai tempi di una prova su 600 simulazioni e da quelli del run precedente per la GEE; il costo reale delle situazioni più grandi dello scenario con la dose è stato circa 1,4 volte quello stimato.

# Scelte del PI sul riassunto (2026-10-24)

- **Come entra il nuovo meccanismo nel paper: raccomandazione accolta.** Testo principale come quarto meccanismo, regola della Figura 1 sui quattro meccanismi, frase sul modello misto, dettagli nel supplemento. Vedi [Il quarto meccanismo entra nel testo principale: regola della Figura 1 sui quattro meccanismi, frase sul modello misto, dettagli nel supplemento (D-55)](../decisions/D-55.md).

# Dubbi rinviati

Nessuno. Il testo da scrivere passa all'[Aggiornamento del manoscritto con le analisi aggiunte (S7f)](S7f.md), dubbio ricevuto 5.
