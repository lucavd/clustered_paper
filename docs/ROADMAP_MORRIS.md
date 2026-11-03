---
type: Roadmap
title: "Roadmap verso la conformita' a Morris 2019"
description: "Gli step per portare simulazione e manoscritto a conformita' con Morris, White e Crowther (2019): stato, prodotti, dipendenze e attese di ciascuno."
status: stable
generated: { by: mmm-agent/skill-v4, at: 2026-10-21T21:46:43+02:00 }
sources:
  - id: morris2019
    resource: https://doi.org/10.1002/sim.8086
    title: "Morris, White, Crowther (2019). Using simulation studies to evaluate statistical methods. Statistics in Medicine 38:2074-2102"
  - id: audit
    resource: "audit/MORRIS2019_conformita.csv"
    title: "Audit dei 28 requisiti, punto di partenza"
---
# MMM: schema per allineare analisi e manoscritto a Morris et al. (2019)

Punto di partenza: audit `MORRIS2019_conformita.csv` (28 requisiti: 4 conformi, 10 parziali, 12 non conformi, 2 non verificabili).
Obiettivo: ogni requisito chiuso come *Conforme* oppure *Dichiarato* (limite scritto nel paper con il motivo). Nessun numero del manoscritto si cambia senza la riga di verifica corrispondente.

Il protocollo di sessione ricalca quelli di `notion-error-budget-session` e `geo-spatialtrans-session`, adattato a questo progetto: skill `mmm-session`.

## 0. Stato degli step

| step | stato | data | report |
|---|---|---|---|
| [baseline (S0)](steps/S0.md) | chiuso, approvato dal PI | 2026-10-15 | `reports/morris/S0.html` |
| [protocollo datato (S1)](steps/S1.md) | chiuso, riserve risolte ([LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md), [Valori del modello simulato confrontati con la coorte LUS (D-10)](decisions/D-10.md), P1 dichiarata), approvato dal PI | 2026-10-15 | `reports/morris/S1.html` |
| [diagnostica sui risultati esistenti (S2)](steps/S2.md) | chiuso, approvato dal PI ([Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md), [Replicati con disegno non stimabile esclusi dal denominatore (D-11)](decisions/D-11.md), [Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md) decise; [Valori veri ricalcolati, tabelle precedenti superate (D-12)](decisions/D-12.md) tecnica) | 2026-10-16 | `reports/morris/S2.html` |
| [rerun correttivo (S3)](steps/S3.md) | chiuso, approvato dal PI ([Parametri dello scenario ancorato alla coorte (D-14)](decisions/D-14.md), [Organizzazione dei file del rerun correttivo (D-15)](decisions/D-15.md); Q8b e Q9 decise: [Sensibilita' del Wald del GLMM con ottimizzatore bobyqa (D-16)](decisions/D-16.md), [Attesa sul LRT formulata male, nessun rerun (D-17)](decisions/D-17.md)) | 2026-10-17 | `reports/morris/S3.html` |
| [riaggregazione con MCSE (S4)](steps/S4.md) | chiuso, approvato dal PI ([Nessun valore vero condizionale per GenBinaryY (D-18)](decisions/D-18.md), [Regole tecniche dell'aggregazione (D-19)](decisions/D-19.md), [GLMM con pendenze casuali non stimabile nelle serie piu' corte (D-20)](decisions/D-20.md)) | 2026-10-17 | `reports/morris/S4.html` |
| [inferenza per blocchi (S5)](steps/S5.md) | chiuso, approvato dal PI ([Regole tecniche dell'inferenza fra celle (D-23)](decisions/D-23.md) tecnica; [Replicati per cella: restano 1.000 anche dopo l'inferenza per blocchi (D-24)](decisions/D-24.md), [Il conteggio delle celle vinte esce dal manoscritto (D-25)](decisions/D-25.md), [Lo scenario ancorato alla coorte entra nel testo principale (D-26)](decisions/D-26.md) decise dal PI) | 2026-10-18 | `reports/morris/S5.html` |
| [tabelle e figure (S6)](steps/S6.md) | chiuso, approvato dal PI; regola della Figura 1 e verifiche della regola di ricaduta rinviate ([Verifica della regola di ricaduta del modello marginalizzato in uno step dedicato (S6b)](decisions/D-33.md), decisa dal PI); [Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente (D-32)](decisions/D-32.md) decisa e fatta; [Figure senza titoli sui pannelli, spiegazione in legenda (D-31)](decisions/D-31.md) decisa dal PI; [Regole tecniche di tabelle e figure (D-28)](decisions/D-28.md) tecnica; [Nessuna simulazione in piu' nello scenario ancorato alla coorte (D-29)](decisions/D-29.md), [Il nested-loop plot dell'errore di I tipo e' la Figura 2 (D-30)](decisions/D-30.md) decise dal PI) | 2026-10-19 | `reports/morris/S6.html` |
| [verifica della regola di ricaduta (S6b)](steps/S6b.md) | chiuso, approvato dal PI | 2026-10-21 | `reports/morris/S6b.html` |
| [manoscritto (S7)](steps/S7.md) | chiuso, approvato dal PI ([Il manoscritto non data le regole scelte dagli autori (D-40)](decisions/D-40.md), [Formato del Word: Arial, tutto nero, testo a 11 punti (D-41)](decisions/D-41.md)) | 2026-10-21 | `reports/morris/S7.html` |
| [rianalisi senza nuove simulazioni (S7b)](steps/S7b.md) | chiuso, approvato dal PI ([Il rapporto fra errore standard e dispersione nei sottocampioni si corregge per la popolazione finita (D-44)](decisions/D-44.md), [Gradi di liberta' della t: N - p nel testo principale, N - 2 di Li e Redden come sensibilita' nel supplemento (D-45)](decisions/D-45.md)) | 2026-10-22 | `reports/morris/S7b.html` |
| [pilota dei test a livello di soggetto (S7c)](steps/S7c.md) | chiuso, approvato dal PI ([Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](decisions/D-46.md), [Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](decisions/D-47.md)) | 2026-10-22 | `reports/morris/S7c.html` |
| [run con i metodi aggiunti (S7d)](steps/S7d.md) | da fare | | |
| [terzo meccanismo (S7e)](steps/S7e.md) | da fare | | |
| [aggiornamento del manoscritto (S7f)](steps/S7f.md) | da fare | | |
| [revisione finale (S8)](steps/S8.md) | revisione, audit e piano fatti; attuazione nello step S8b ([D-69](decisions/D-69.md)) | | |
| [attuazione della revisione finale (S8b)](steps/S8b.md) | chiuso, approvato dal PI ([D-81](decisions/D-81.md)) | 2026-11-02 | `reports/morris/S8b.html` |
| [esito costante come nessun risultato (S9)](steps/S9.md) | da fare ([D-82](decisions/D-82.md)) | | |
| [revisione di numeri e stile (S10)](steps/S10.md) | da fare, dopo S9 | | |

## 1-2. Regole di lavoro e rituale

Le regole di lavoro e il rituale di apertura e chiusura sono nel protocollo di sessione v3 (skill `mmm-session`), registrato come [Protocollo di sessione v3 (D-22)](decisions/D-22.md). Dal 2026-10-17 sostituiscono le regole e il rituale scritti qui, che restano nella storia git (commit `47d4572`).

## 3. Decisioni richieste al PI (prima di S3)

Tabella storica delle opzioni proposte prima del rerun correttivo. Lo stato attuale di ogni decisione e' nell'[indice delle decisioni](decisions/index.md).

| id | questione | opzioni | raccomandazione |
|---|---|---|---|
| [Dati sotto H0 di GenBinaryY rigenerati in M3-M5 (D-01)](decisions/D-01.md) | H0 ≠ H1 nei parametri di associazione (GenBinaryY M3-M5, rilievo B3) | (a) rigenerare l'arm H0 con sigma/gamma uguali a H1; (b) dichiarare la differenza nel paper | **(a)**: 2.880 chunk; con lo stesso seme l'accoppiamento CRN fra gli arm resta |
| [Varianza della pendenza su dose allineata in M5 (D-02)](decisions/D-02.md) | Varianza della pendenza su dose in M5: 0.01 in AR1, 0.001 in GLMM-native (B4) | (a) allineare AR1 a 0.001 e rigenerare AR1 M5 (1.920 chunk); (b) dichiarare | **(a)**: ripristina la generalizzazione stretta che il paper afferma |
| [Replicati per cella: restano 1.000 (D-03)](decisions/D-03.md) | n_sim = 1000 (F1) | (a) tenerlo e giustificarlo con la MCSE ottenuta; (b) aggiungere replicati nelle celle da cui dipende una conclusione | **(a)**, con (b) solo se [inferenza per blocchi (S5)](steps/S5.md) mostra una conclusione entro 2 MCSE |
| [Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md) | Fit con warning di lme4/MMLB oggi contati come usabili (E2) | (a) primaria invariata più sensibilità che li esclude; (b) escluderli nella primaria | **(a)**: decidere dopo la quantificazione in [diagnostica sui risultati esistenti (S2)](steps/S2.md) |
| [Stesso seme fra celle, inferenza per blocchi di replicato (D-05)](decisions/D-05.md) | Seme comune fra arm, N, T e meccanismi ([F3: Dipendenza fra dataset simulati dichiarata](audit/requirements/F3.md)) | (a) tenerlo, dichiararlo e correggere l'inferenza; (b) rifare tutto con stream indipendenti (circa 50 h) | **(a)**: Morris §5.4 considera i CRN un vantaggio per i confronti |
| [Codice pubblico alla sottomissione, con DOI (D-06)](decisions/D-06.md) | Disponibilità del codice (H5) | repository pubblico al tag di submission con DOI Zenodo | sì |

## 4. Step

| step | titolo | dove | requisiti Morris | dipende da | prodotti | attese pre-registrate |
|---|---|---|---|---|---|---|
| [baseline (S0)](steps/S0.md) | Baseline sul server e protocollo di sessione | server | [H5: Disponibilita' di codice e dati](audit/requirements/H5.md), [H6: Protocollo / pre-specificazione](audit/requirements/H6.md) | | ramo `morris-v4`; `rev2/` tracciata (file grandi solo in manifest); skill `mmm-session`; `DECISIONS.md`; md5 baseline di `results_full/`, `tables/` e manoscritto | 28.800 chunk; md5 del manoscritto sul server identici al Mac |
| [protocollo datato (S1)](steps/S1.md) | Protocollo ADEMP di revisione, datato prima di ogni rerun | server | A1, B1, B2, B5, B6, C1, C2, D1, D3, E1, F1, F2, F3, H1, H6 | [baseline (S0)](steps/S0.md) | `docs/PROTOCOL_ADEMP.md`; tabella dei parametri estratta *dal codice* con uno script, non a mano | la tabella dei parametri coincide con i generatori (test automatico) |
| [diagnostica sui risultati esistenti (S2)](steps/S2.md) | Diagnostica sui chunk esistenti, senza nuove simulazioni | server | B2, C2, E2, E3, F3 | [protocollo datato (S1)](steps/S1.md) | (a) quota di fit usabili con warning per cella; (b) sim_id mancanti rieseguiti con il loro seme e motivo d'errore registrato; (c) valori veri per meccanismo × scenario × T; (d) distribuzione realizzata delle lunghezze delle serie; (e) correlazione delle stime fra arm e fra meccanismi a parità di `(model, i)`; (f) pilota H0 corretto su 12 celle; (g) confronto della coorte LUS con i valori del DGM: prevalenza, interazione, varianza dell'intercetta o ICC ([Valori del modello simulato confrontati con la coorte LUS (D-10)](decisions/D-10.md)) | (b) mancanti = 1.440 × 1.000 − attesi dai file, motivo noto per il 100%; (d) quota di serie piene ≈ 0,70 ± 0,01 |
| [rerun correttivo (S3)](steps/S3.md) | Rerun correttivo e scenario ancorato alla coorte | server, 50 worker | B3, B4, B6, E3, D1 | [diagnostica sui risultati esistenti (S2)](steps/S2.md) e [Dati sotto H0 di GenBinaryY rigenerati in M3-M5 (D-01)](decisions/D-01.md)/[Varianza della pendenza su dose allineata in M5 (D-02)](decisions/D-02.md)/[Generatore di numeri casuali fissato nei rerun (D-07)](decisions/D-07.md)/[LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md)/[Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md) | `results_v4/`: GenBinaryY M3-M5 H0 (2.880 chunk) e AR1 M5 (1.920 chunk); scenario M1c ancorato alla coorte (1.440 chunk, [Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md)), con calibrazione dei parametri verificata prima del run; `glmer_lrt` con pendenze casuali in M3-M5 ([LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md)); RNGkind L'Ecuyer-CMRG esplicito ([Generatore di numeri casuali fissato nei rerun (D-07)](decisions/D-07.md)); il generatore registra i fallimenti con il motivo; test di parità H0/H1 su tutti i 15 generatori | test di parità PASS 15/15; zero errori di chunk; stime H1 bit-identiche alle precedenti dove il DGM non cambia |
| [riaggregazione con MCSE (S4)](steps/S4.md) | Riaggregazione unica con MCSE | server | E1, E2, E4, C2 | [rerun correttivo (S3)](steps/S3.md) | `R/postrun/agg_v4.R` → `tables/V4_*.csv`; per cella e metodo: bias, SE empirico, SE del modello, errore relativo dello SE del modello, MSE/RMSE, copertura, copertura bias-eliminated, errore di I tipo, potenza raw e size-adjusted, quota usabile, tutti con MCSE; due definizioni di usabilità ([Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md)) | le celle non toccate da [rerun correttivo (S3)](steps/S3.md) riproducono `O_level_all.csv` entro 1e-12 |
| [inferenza per blocchi (S5)](steps/S5.md) | Inferenza fra celle sotto CRN e adeguatezza di n_sim | server | E4, F1, F3 | [riaggregazione con MCSE (S4)](steps/S4.md) | bootstrap per blocchi di replicato (B = 2000) per worst case, celle vinte, rango combinato, Spearman e gap per strato; BH + BY; matrice affermazione × MCSE | ogni affermazione comparativa del paper è classificata come risolta (> 2 MCSE) o non risolta |
| [tabelle e figure (S6)](steps/S6.md) | Tabelle e figure secondo Morris | server | [H2: MCSE in tabelle e figure](audit/requirements/H2.md), [H3: Precisione coerente con la MCSE](audit/requirements/H3.md), [H4: Esplorazione grafica delle stime](audit/requirements/H4.md) | [inferenza per blocchi (S5)](steps/S5.md) | Tabella 1 ADEMP; tabella dei parametri; tabelle con MCSE (o MCSE massima in didascalia); nested-loop plot dell'errore di I tipo; Figura 1 con intervalli MC e regola completa (rami 2-3 occasioni e 25-45 cluster); zip plot nel supplemento | §9 di figure-style senza sovrapposizioni |
| [verifica della regola di ricaduta (S6b)](steps/S6b.md) | Verifica della regola di ricaduta del modello marginalizzato | server | robustezza della strategia di MMM (aggiunto dal PI, [Verifica della regola di ricaduta del modello marginalizzato in uno step dedicato (S6b)](decisions/D-33.md)) | [tabelle e figure (S6)](steps/S6.md) | sensibilita' della soglia "al bordo" (0,01 e 0,0001) sui dataset esistenti; nuovo scenario con pendenza vera grande della deviazione standard latente; scelta della regola della Figura 1 | la regola di ricaduta resta tarata dove il modello senza pendenza e' sbagliato, oppure il paper la presenta come verificata solo con pendenza circa zero |
| [manoscritto (S7)](steps/S7.md) | Riscrittura del manoscritto (v4) | server | tutti i requisiti di reporting | [verifica della regola di ricaduta (S6b)](steps/S6b.md) | Metodi in ordine ADEMP; ogni numero rifonte su `V4_*` con un nuovo `V-nn`; arrotondamento alla MCSE; 7 commenti di Luca nel docx; code/data availability; supplemento rigenerato con sessionInfo | zero numeri senza puntatore; zero affermazioni "non risolte" scritte come differenze |
| [rianalisi senza nuove simulazioni (S7b)](steps/S7b.md) | Gradi di liberta' della t, valore nullo marginale, seconda applicazione clinica | server | da mappare nel ri-audit Morris | [manoscritto (S7)](steps/S7.md), [Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](decisions/D-42.md) | regola della Figura 1 con N - p e N - 2 su 1.896 punti; rifiuto atteso dovuto al nullo marginale; applicazione su `respiratory` di geepack | nel file dello step, all'apertura |
| [pilota dei test a livello di soggetto (S7c)](steps/S7c.md) | Disegno e pilota dei test a livello di soggetto e di permutazione | server | da mappare nel ri-audit Morris | [rianalisi senza nuove simulazioni (S7b)](steps/S7b.md) | test con nullo coincidente con l'interazione trattamento per tempo; pilota su poche celle; scelta del PI | nel file dello step, all'apertura |
| [run con i metodi aggiunti (S7d)](steps/S7d.md) | Correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto | server | da mappare nel ri-audit Morris | [pilota dei test a livello di soggetto (S7c)](steps/S7c.md) | 1.896 punti rigenerati con gli stessi semi; GEE indipendente e AR(1); fewclust e CR2; test scelti; regola riapplicata | riproduzione replicato per replicato della GEE exchangeable salvata |
| [terzo meccanismo (S7e)](steps/S7e.md) | Terzo meccanismo generatore | server | da mappare nel ri-audit Morris | [run con i metodi aggiunti (S7d)](steps/S7d.md) | meccanismo neutrale su 480 punti della griglia principale con tutte le strategie | calibrazione sui coefficienti marginali verificata prima del run |
| [aggiornamento del manoscritto (S7f)](steps/S7f.md) | Aggiornamento del manoscritto con le analisi aggiunte | server | tutti i requisiti di reporting | [terzo meccanismo (S7e)](steps/S7e.md), [fewclust allineato alla regola del paper e citato nel manoscritto (D-43)](decisions/D-43.md) | `manuscript/v5/`; limitazioni riscritte; citazione di fewclust | stessi controlli del manoscritto v4 |
| [revisione finale (S8)](steps/S8.md) | Revisione avversariale e ri-audit Morris | sotto-agenti | gate di submission | [aggiornamento del manoscritto (S7f)](steps/S7f.md) | `revisione-e-piano-implementazione` sul v4; `MORRIS2019_conformita.csv` con `esito_v4` | 28/28 Conforme o Dichiarato; nessun MISMATCH bloccante |

## 5. Mappa paper ← step

| sezione | cosa cambia | step |
|---|---|---|
| §2 (nuova apertura) | Aims; Tabella 1 ADEMP | [protocollo datato (S1)](steps/S1.md), [manoscritto (S7)](steps/S7.md) |
| §2.2 DGM | parametri numerici; troncamento delle serie; struttura effettiva della "random slope" in GenBinaryY; correzione della descrizione di GLMM-native; parità H0/H1 | [protocollo datato (S1)](steps/S1.md), [rerun correttivo (S3)](steps/S3.md), [manoscritto (S7)](steps/S7.md) |
| §2.3 Griglia | n_sim e MCSE attesa; schema dei semi e CRN; motivo dei replicati mancanti | [protocollo datato (S1)](steps/S1.md), [diagnostica sui risultati esistenti (S2)](steps/S2.md), [manoscritto (S7)](steps/S7.md) |
| §2.4 Misure | definizione di usabilità (errori contro warning); formule MCSE; inferenza per blocchi | [riaggregazione con MCSE (S4)](steps/S4.md), [inferenza per blocchi (S5)](steps/S5.md), [manoscritto (S7)](steps/S7.md) |
| §2.1 Software | sessionInfo; opzioni di stima; timeout di geesmv | [protocollo datato (S1)](steps/S1.md), [manoscritto (S7)](steps/S7.md) |
| §3.1-3.15 | tutti i valori da `V4_*`, con MCSE; affermazioni sotto la MCSE riformulate | S4-S7 |
| §3.8 | argomento di identificabilità con serie troncate | [diagnostica sui risultati esistenti (S2)](steps/S2.md), [manoscritto (S7)](steps/S7.md) |
| §4.2 | la regola della Figura 1 si descrive con i suoi effetti, senza dire quando e' stata scelta ([Il manoscritto non data le regole scelte dagli autori (D-40)](decisions/D-40.md)) | [manoscritto (S7)](steps/S7.md) |
| §4.4 | limitazioni mitigate o dichiarate ([Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](decisions/D-42.md)) | S7b-S7f |
| §4.4 | limiti: pendenza piccola, CRN, vantaggio "in casa" dei metodi | [protocollo datato (S1)](steps/S1.md), [manoscritto (S7)](steps/S7.md) |
| Figura 1 | intervalli MC; rami mancanti | [tabelle e figure (S6)](steps/S6.md) |
| Fine testo | code and data availability (GitHub + Zenodo) | [baseline (S0)](steps/S0.md), [manoscritto (S7)](steps/S7.md) |

## 6. Stima dei tempi

- [rerun correttivo (S3)](steps/S3.md): circa 4.800 chunk, piu' 1.440 dello scenario M1c ([Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md)). Alla velocità misurata (circa 19 chunk/min con 100 worker) servono 8-10 h con 50 worker. Consuntivo: 18.720 chunk (il blocco solo LRT di [LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md) aggiunge 12.480 chunk), 20,1 h con 50 worker, 1.000,7 ore-worker.
- [diagnostica sui risultati esistenti (S2)](steps/S2.md), [riaggregazione con MCSE (S4)](steps/S4.md) e [inferenza per blocchi (S5)](steps/S5.md) girano sui chunk esistenti e richiedono minuti o poche ore.
- I restanti step sono lavoro di scrittura.

## Aggiunta del 2026-10-23

- [Intervallo di confidenza dal test a inversione di segno (S7d2)](steps/S7d2.md): prima dell'[Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md), per [Uno step dedicato prova l'intervallo di confidenza ottenuto invertendo il test a inversione di segno, prima della riscrittura del paper (D-52)](decisions/D-52.md).
