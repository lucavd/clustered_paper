---
type: Protocol
title: "Protocollo ADEMP della revisione dello studio di simulazione"
description: "Obiettivi, meccanismi generatori, estimandi, metodi e misure di prestazione secondo lo schema ADEMP di Morris 2019, datato e versionato prima di ogni rerun (versione 1.6)."
status: stable
generated: { by: mmm-agent/skill-v1, at: 2026-10-16T14:18:48+02:00 }
sources:
  - id: morris2019
    resource: https://doi.org/10.1002/sim.8086
    title: "Morris, White, Crowther (2019). Using simulation studies to evaluate statistical methods. Statistics in Medicine 38:2074-2102"
  - id: parametri
    resource: "../tables/V4_dgm_parameters.csv"
    title: "Tabella dei parametri estratta dal codice eseguito"
---
# Protocollo ADEMP di revisione dello studio di simulazione MMM (v4)

- **Versione**: 1.8, 2026-10-22. La 1.0 (`cc534d6`, [protocollo datato (S1)](steps/S1.md)) e' stata emendata in 1.1 ([protocollo datato (S1)](steps/S1.md), [LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md) e [Valori del modello simulato confrontati con la coorte LUS (D-10)](decisions/D-10.md)), in 1.2 (apertura di [diagnostica sui risultati esistenti (S2)](steps/S2.md), prima di ogni calcolo: [Confronto con la coorte attraverso grandezze analoghe (D-10b)](decisions/D-10b.md) e verita' marginale GLMM-native per quadratura esatta) in 1.3 (chiusura di [diagnostica sui risultati esistenti (S2)](steps/S2.md), prima di [rerun correttivo (S3)](steps/S3.md) e [riaggregazione con MCSE (S4)](steps/S4.md): [Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md), [Replicati con disegno non stimabile esclusi dal denominatore (D-11)](decisions/D-11.md), [Valori veri ricalcolati, tabelle precedenti superate (D-12)](decisions/D-12.md)) in 1.4 ([Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md), scenario ancorato alla coorte, prima di [rerun correttivo (S3)](steps/S3.md)) e in 1.5 ([Parametri dello scenario ancorato alla coorte (D-14)](decisions/D-14.md) e [Organizzazione dei file del rerun correttivo (D-15)](decisions/D-15.md), parametri di M1c e struttura dei chunk, all'apertura di [rerun correttivo (S3)](steps/S3.md)) e in 1.6 (apertura di [rianalisi senza nuove simulazioni (S7b)](steps/S7b.md), prima di ogni calcolo dello step: gradi di liberta' della t, N - p e N - 2; rifiuto atteso dovuto al nullo marginale; seconda applicazione clinica, sez. 4b) e in 1.7 (chiusura di [pilota dei test a livello di soggetto (S7c)](steps/S7c.md): test a inversione di segno aggiunto per il run su tutti i punti di disegno, [Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](decisions/D-46.md), [Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](decisions/D-47.md)) e in 1.8 (apertura del [run con i metodi aggiunti (S7d)](steps/S7d.md), prima di ogni calcolo del run: GEE con correlazione indipendente e AR(1), correzioni di fewclust 0.2.0, CR2, geesmv sui replicati 1-50, [geesmv solo sulle prime 50 simulazioni di ogni scenario nel run con i metodi aggiunti (D-48)](decisions/D-48.md), [fewclust 0.2.0 registrato dall'agente nel ramo package cosi' com'era sul portatile del PI (D-49)](decisions/D-49.md)). Sez. 11.
- **Stato**: scritto e committato **prima** di qualsiasi rerun ([rerun correttivo (S3)](steps/S3.md)) e di qualsiasi riaggregazione ([riaggregazione con MCSE (S4)](steps/S4.md)). Le attese di [protocollo datato (S1)](steps/S1.md) sono state committate prima dei calcoli di [protocollo datato (S1)](steps/S1.md) (`b0826c5`, 2026-10-15 21:55).
- **Riferimento**: Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Stat Med* 2019;38:2074-2102.
- **Fonti primarie**: `tables/V4_dgm_parameters.csv` (estratta dal codice eseguito con `R/protocol/extract_dgm_params.R`), `R/run_full.R`, `.rundir/`, `DECISIONS.md`.
- **Che cosa non e'**: non e' un protocollo pre-registrato dello studio originale. Il run di produzione e' del 2026-08-05/09; questo documento fissa *in anticipo* che cosa la revisione v4 rigenera, misura e riporta (requisito [H6: Protocollo / pre-specificazione](audit/requirements/H6.md), vedi sez. 9).

## 0. Riepilogo ADEMP (bozza della Tabella 1 del v4)

| | contenuto |
|---|---|
| **A**ims | Valutare, per l'interazione `time:trt` su scala logit in dati binari longitudinali con pochi cluster e serie brevi: (1) la validita' dell'errore di I tipo dei test; (2) la potenza, grezza e size-adjusted; (3) distorsione, precisione e copertura delle stime rispetto al loro estimando; (4) la quota di fit utilizzabili. Il confronto riguarda modelli marginalizzati multilivello (MMM), GEE con sandwich standard e con tre correzioni small-sample, GLMM con test di Wald e di rapporto di verosimiglianza, e GLM con riscalamento shrinkage. La robustezza si valuta su tre meccanismi generatori, ciascuno "di casa" per un metodo diverso o per nessuno. |
| **D**ata-generating mechanisms | Fattoriale completo: N (cluster) in {10, 15, 20, 25, 35, ..., 145} (16 livelli) x T (occasioni massime) in {2, 3, 4, 6, 8, 10} x 5 scenari (M1-M5) x 3 meccanismi (GenBinaryY, GLMM-native, GLMM-AR1) = 1.440 punti di disegno, ciascuno sotto H1 e H0 (2.880 celle). Parametri: sez. 2 e `tables/V4_dgm_parameters.csv`. |
| **E**stimands | Coefficiente `time:trt`: marginale (MMM, GEE, shrinkage) e condizionale sugli effetti casuali (GLMM). Valori veri per meccanismo: sez. 3. |
| **M**ethods | MMM (`MMLB::mm`), GEE exchangeable (`geepack`) con sandwich standard e Mancl-DeRouen, Kauermann-Carroll, Fay-Graubard (`geesmv`), GLMM (`lme4::glmer`, Laplace) con Wald e LRT, GLM shrinkage (Zeger). Dettagli e opzioni: sez. 4. |
| **P**erformance measures | Quota usabile; distorsione; SE empirico; SE del modello ed errore relativo; MSE e RMSE; copertura e copertura bias-eliminated; errore di I tipo; potenza grezza e size-adjusted. Tutte con MCSE. Sez. 5. |

## 1. Aims (requisito A1)

Obiettivi della **sola simulazione**. La raccomandazione pratica e l'applicazione alla coorte LUS neonatale sono parti distinte del manoscritto e non sono obiettivi di questo protocollo.

1. Errore di I tipo nominale al 5% del test sull'interazione `time:trt`, per metodo, al variare di N e T.
2. Potenza al 5%, grezza e size-adjusted (valore critico preso dalla distribuzione nulla della stessa cella).
3. Distorsione, SE empirico, accuratezza dello SE del modello, MSE e copertura al 95% della stima di `time:trt`, ciascuno rispetto all'estimando del metodo (sez. 3).
4. Quota di replicati con fit utilizzabile, come esito di primo livello (Morris sez. 5.1).

## 2. Data-generating mechanisms (requisiti B1, B2, B5, B6)

### 2.1 Tabella dei parametri estratta dal codice

La tabella sotto e' generata da `R/protocol/s1_tables.R` a partire da `tables/V4_dgm_parameters.csv`, che `R/protocol/extract_dgm_params.R` costruisce intercettando gli argomenti passati ai motori di generazione durante `make_all_data()`, piu' una lettura statica del corpo dei generatori. La tabella e' verificata da `R/testing/test_dgm_parameters.R`: un rigeneratore che conosce solo la tabella riproduce bit per bit i dataset del codice eseguito (test P3, sez. 8).

<!-- PARAMS:BEGIN -->
| DGM | M | mean model (coefficient order) | coefficients H1 | dependence (H1) | H0 differs from H1 in |
|---|---|---|---|---|---|
| GenBinaryY | 1 | ~time * trt ((Intercept), time, trt, time:trt) | marginal -2, 0.1, 0.2, 0.2 | lv ~1, sigma (1); t ~1, gamma (1); q = 20 | tested coefficient only |
| GenBinaryY | 2 | ~time * trt + age ((Intercept), time, trt, age, time:trt) | marginal -2, 0.1, 0.2, 0.01, 0.2 | lv ~1, sigma (1); t ~1, gamma (1); q = 20 | tested coefficient only |
| GenBinaryY | 3 | ~time * trt ((Intercept), time, trt, time:trt) | marginal -2, 0.1, 0.2, 0.2 | lv ~time, sigma (1, 0.01); t ~1, gamma (1); q = 20 | time:trt and sigma (H0: sigma (1, 0), gamma (1)) |
| GenBinaryY | 4 | ~time * trt ((Intercept), time, trt, time:trt) | marginal -2, 0.2, 0.2, 0.05 | lv ~time, sigma (1, 0.01); t ~trt, gamma (1, 0.01); q = 20 | time:trt and sigma, gamma (H0: sigma (1, 0), gamma (1, 0)) |
| GenBinaryY | 5 | ~time * trt + dose ((Intercept), time, trt, dose, time:trt) | marginal -2, 0.1, 0.2, 0.05, 0.2 | lv ~time + dose, sigma (1, 0.01, 0.01); t ~trt + steroids, gamma (1, 0.1, 0.05); q = 20 | time:trt and sigma, gamma (H0: sigma (1, 0, 0), gamma (1, 0, 0)) |
| GLMM-native | 1 | ~time * trt ((Intercept), time, trt, time:trt) | conditional -2.32, 0.116, 0.232, 0.232 | b_i on ~1, D = [1] | tested coefficient only |
| GLMM-native | 2 | ~time * trt + age ((Intercept), time, trt, age, time:trt) | conditional -2.32, 0.116, 0.232, 0.012, 0.232 | b_i on ~1, D = [1] | tested coefficient only |
| GLMM-native | 3 | ~time * trt ((Intercept), time, trt, time:trt) | conditional -2.32, 0.116, 0.232, 0.232 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001] | tested coefficient only |
| GLMM-native | 4 | ~time * trt ((Intercept), time, trt, time:trt) | conditional -2.32, 0.232, 0.232, 0.058 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001] | tested coefficient only |
| GLMM-native | 5 | ~time * trt + dose ((Intercept), time, trt, dose, time:trt) | conditional -2.32, 0.116, 0.232, 0.058, 0.232 | b_i on ~1 + time + dose, D = [1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.001] | tested coefficient only |
| GLMM-AR1 | 1 | ~time * trt ((Intercept), time, trt, time:trt) | marginal target -2, 0.1, 0.2, 0.2 | b_i on ~1, D = [1]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 2 | ~time * trt + age ((Intercept), time, trt, age, time:trt) | marginal target -2, 0.1, 0.2, 0.01, 0.2 | b_i on ~1, D = [1]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 3 | ~time * trt ((Intercept), time, trt, time:trt) | marginal target -2, 0.1, 0.2, 0.2 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 4 | ~time * trt ((Intercept), time, trt, time:trt) | marginal target -2, 0.2, 0.2, 0.05 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 5 | ~time * trt + dose ((Intercept), time, trt, dose, time:trt) | marginal target -2, 0.1, 0.2, 0.05, 0.2 | b_i on ~1 + time + dose, D = [1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.01]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |

Design draws, in execution order (table rows `draw_*`). They are the same in the three mechanisms and in both arms (checked); the scenarios differ only in how many covariates they draw:

- M1: nclust, trt
- M2: nclust, trt, age
- M3: nclust, trt, age
- M4: nclust, trt, age
- M5: nclust, trt, age, sex, steroids, dose

Full list (M5):

- follow-up: `nclust = ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))` (p_full_followup = 0.7)
- `nclust` <- `runif(n = N)` (level=subject; standardised=FALSE; assigned_as=nclust)
- `nclust` <- `sample(x = seq_len(ntime - 1), size = N, replace = TRUE)` (level=subject; standardised=FALSE; assigned_as=nclust)
- `trt` <- `rbinom(n = N, size = 1, prob = 0.5)` (level=subject; standardised=FALSE; assigned_as=trt)
- `age` <- `rnorm(n = N, mean = 50, sd = 10)` (level=subject; standardised=TRUE; assigned_as=age_raw)
- `sex` <- `rbinom(n = N, size = 1, prob = 0.5)` (level=subject; standardised=FALSE; assigned_as=sex)
- `steroids` <- `sample(x = 1:3, size = N, replace = TRUE)` (level=subject; standardised=FALSE; assigned_as=steroids)
- `dose` <- `rnorm(n = length(id))` (level=observation; standardised=TRUE; assigned_as=dose)

Covariates drawn but entering no part of the DGM (they consume random numbers):

- GenBinaryY M3: age
- GenBinaryY M4: age
- GenBinaryY M5: age, sex
- GLMM-native M3: age
- GLMM-native M4: age
- GLMM-native M5: age, sex, steroids
- GLMM-AR1 M3: age
- GLMM-AR1 M4: age
- GLMM-AR1 M5: age, sex, steroids

Seed: `set.seed(20260805L + 1000L * model + i)` per replicate; RNGkind in the workers: `L'Ecuyer-CMRG, Inversion, Rejection` (set by crew 1.3.2 / mirai 2.7.2 in every worker; NOT set by run_full.R (plain R session: Mersenne-Twister, Inversion, Rejection)).
<!-- PARAMS:END -->

### 2.2 Struttura delle serie (B2)

T e' un **massimo**, non il numero di occasioni di ogni soggetto: con probabilita' 0,7 un soggetto e' osservato a tutte le T occasioni, altrimenti a un numero uniforme in 1..T-1. A T = 2 quindi circa il 30% dei soggetti ha una sola osservazione. Il v4 lo dichiara in Metodi e riformula l'argomento di identificabilita' di sez. 3.8; la distribuzione realizzata delle lunghezze si misura in [diagnostica sui risultati esistenti (S2)](steps/S2.md)(d).

Dettaglio che tocca la sequenza dei numeri casuali: `ifelse()` valuta `sample()` solo se almeno un soggetto e' troncato, quindi in un dataset con tutti i soggetti completi (probabilita' 0,7^N, non trascurabile solo a N = 10) lo stream avanza di N numeri in meno. L'effetto e' identico nei due arm e nei tre meccanismi, quindi l'accoppiamento CRN non ne risente.

### 2.3 Natura della dipendenza per meccanismo (B5, D3)

- **GenBinaryY** (modello marginalizzato con variabile latente e transizione, MMLB): il carico della variabile latente e' il predittore **lineare** `sigma_ij = x_lv' sigma` (letto nel sorgente di `MMLB::GenBinaryY`: `sig = x.lv %*% sigma`), quindi con `lv.formula = ~time` e `sigma = (1, 0.01)` vale `1 + 0,01 time` e moltiplica **un'unica** variabile latente z per soggetto. E' un'intercetta casuale con carico variabile nel tempo (pendenza di rango 1, perfettamente correlata con l'intercetta), non una pendenza casuale indipendente. Lo stesso vale per gamma (`gamma_ij = x_t' gamma`). E' il meccanismo "di casa" per MMM: le formule di `mm()` coincidono con quelle del generatore sotto H1, ma `mm()` parametrizza `log(sigma)`, quindi il modello di associazione e' congeniale solo in approssimazione (log(1 + 0,01 t) e' circa 0,01 t) e in M1-M2 il valore vero di `log(sigma):(Intercept)` e' 0.
- **GLMM-native**: GLMM logistico con effetti casuali normali `b_i ~ N(0, D)`. E' il meccanismo "di casa" per il GLMM (la struttura casuale del `glmer` di Wald coincide con quella del generatore). I coefficienti condizionali sono **costanti fisse**, pari a circa 1,16 volte il marginale (fattore di Zeger con D[1,1] = 1, arrotondato a tre decimali), non calibrati per replicato.
- **GLMM-AR1**: GLMM con processo latente AR(1) sui residui; i coefficienti condizionali sono calibrati **su ogni replicato** perche' il marginale implicato coincida col target. Non e' il meccanismo di nessun metodo.

### 2.4 Motivazione dei valori (B6)

Motivazioni documentate nel codice: rho = 0,6 per una dipendenza seriale rilevabile senza dominare l'effetto di cluster; sigma2_ar = 0,5 per una varianza latente totale confrontabile fra meccanismi; M4 con interazione 0,05 come regime vicino al nullo. **Effetto 0,2 sull'interazione, intercetta -2 (prevalenza di circa il 12%) e varianza dell'intercetta 1 ([Valori del modello simulato confrontati con la coorte LUS (D-10)](decisions/D-10.md))**: si ancorano alla coorte LUS neonatale. In [diagnostica sui risultati esistenti (S2)](steps/S2.md) si confrontano prevalenza, interazione `time:trt` e varianza dell'intercetta (o ICC) stimate sulla coorte con i valori del DGM; il v4 riporta il confronto e, se i valori non coincidono, dichiara lo scarto e il motivo della scelta. Nella coorte non esistono trattamento ne' una vera dimensione temporale (22 pazienti su 207 con piu' di un timepoint), quindi `time:trt` non e' stimabile: il confronto avviene per analoghi ([Confronto con la coorte attraverso grandezze analoghe (D-10b)](decisions/D-10b.md)), cioe' prevalenza di `lus_severe`, varianze dell'intercetta di paziente e di paziente-timepoint con ICC latente, e i contrasti di regione, disease e invasive sulla scala logit come riferimento di grandezza dell'effetto; il v4 dichiara che 0,2 per occasione non ha un ancoraggio diretto nella coorte. Da dichiarare come limite (sez. 4.4 del paper): la varianza della pendenza 0,001 (SD 0,032 per occasione) rende la pendenza casuale quasi nulla, quindi le conclusioni sulle "random slopes" hanno portata limitata.

### 2.5 Modifiche del DGM previste per il v4 (decise il 2026-10-15)

- **[Dati sotto H0 di GenBinaryY rigenerati in M3-M5 (D-01)](decisions/D-01.md) (a)**: in GenBinaryY M3-M5 l'arm H0 viene rigenerato con sigma e gamma uguali a H1; resta diverso solo il coefficiente testato. 2.880 chunk in [rerun correttivo (S3)](steps/S3.md).
- **[Varianza della pendenza su dose allineata in M5 (D-02)](decisions/D-02.md) (a)**: in GLMM-AR1 M5 D[3,3] passa da 0,01 a 0,001, come in GLMM-native. 1.920 chunk in [rerun correttivo (S3)](steps/S3.md).
- **[LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md) (b)**: in M3-M5 (tre meccanismi, entrambi gli arm) si riesegue la sola tecnica `glmer_lrt` con la struttura casuale del Wald.
- **[Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md)**: scenario ancorato alla coorte M1c (modello di M1; prevalenza marginale 0,5945 e varianza latente dell'intercetta 3,545 dalla coorte LUS), tre meccanismi, N in {10, 15, 20, 25, 35, 45}, T in {2, 3, 4, 6}, 1.440 chunk in [rerun correttivo (S3)](steps/S3.md), seme di M1. Parametri dei generatori pre-registrati all'apertura di [rerun correttivo (S3)](steps/S3.md). Nel paper la catena e': griglia principale, scenario ancorato alla coorte, plasmode sui dati reali.
- **[Parametri dello scenario ancorato alla coorte (D-14)](decisions/D-14.md)** (parametri di M1c, 2026-10-16): intercetta marginale beta0(T) calibrata sull'arm H0 perche' la prevalenza realizzata attesa sul disegno troncato valga 0,5945 (0,2427; 0,1960; 0,1495; 0,0569 per T = 2, 3, 4, 6), stessa beta0(T) in H1; pendenze di M1. GenBinaryY sigma = sqrt(3,545), gamma = 1; GLMM-native D = 3,545 e beta_c = 1,492 x marginale (Zeger, prevalenza realizzata circa 0,589); GLMM-AR1 D = 3,545, sigma2_ar = 0,5, rho = 0,6, calibrazione per replicato. Verifiche: `docs/sessions/2026-10-16_S3.md`, Q10.
- Tutti gli altri generatori e tecniche restano invariati, e i loro chunk sono riusati senza modifiche.

## 3. Estimands (requisiti C1, C2)

- Estimando: coefficiente di `time:trt` sulla scala logit.
- **GenBinaryY**: valore vero marginale = `beta` del generatore (0,2; 0,05 in M4; 0 sotto H0). Per il GLMM su questi dati il target condizionale non e' un parametro del DGM: il v4 lo dichiara come approssimazione.
- **GLMM-native**: valore vero condizionale = `beta_c` (0,232; 0,058 in M4; 0 sotto H0). Il valore vero marginale e' derivato come limite della GEE di indipendenza con il modello medio dello scenario, pesato con la distribuzione realizzata del disegno (troncamento delle serie, trt ~ Bern(0,5)); la versione a griglia bilanciata si riporta accanto per confronto con `tables/F_truth_glmmnative.csv`. Le probabilita' marginali si ottengono con quadratura di Gauss-Hermite esatta in tutti gli scenari: dato il vettore delle covariate, z'b e' normale con varianza scalare z'Dz (v1.2; in M2 e M5 si integra anche su age e dose ~ N(0,1)). Si tabula per scenario e T in [diagnostica sui risultati esistenti (S2)](steps/S2.md)(c).
- **GLMM-AR1**: valore vero marginale = target (0,2; 0,05; 0) per costruzione sul disegno di ogni replicato; il valore condizionale cambia da replicato a replicato (sd fino a 0,0048 in [diagnostica sui risultati esistenti (S2)](steps/S2.md)(c)) e in [riaggregazione con MCSE (S4)](steps/S4.md) si ricalcola per replicato ([Valori veri ricalcolati, tabelle precedenti superate (D-12)](decisions/D-12.md)).
- **Fonti dei valori veri** ([Valori veri ricalcolati, tabelle precedenti superate (D-12)](decisions/D-12.md)): `tables/V4_S2_truth_values.csv`, `tables/V4_S2_truth_glmmnative.csv`. Sotto H0 la verita' marginale GLMM-native vale da -0,0011 a -0,0030, non 0.

## 4. Methods (requisiti D1, D3)

Codice eseguito: `.rundir/` (M3-M5) e `_superseded/rundir_m1m2_job23035d84/` (M1-M2), uguali per M1-M2 (test P4). Ogni replicato e' analizzato da tutti i metodi sullo stesso dataset (D2).

| metodo | implementazione e opzioni |
|---|---|
| MMM | `MMLB::mm`, `q = 20`; `lv.formula` e `t.formula` uguali a quelle del generatore GenBinaryY H1 dello scenario |
| GEE | `geepack::geeglm`, `family = binomial`, `corstr = "exchangeable"`, `scale.fix = TRUE`, sandwich robusto |
| GEE small-sample | `geesmv::GEE.var.md`, `GEE.var.kc`, `GEE.var.fg` (exchangeable), ciascuno sotto timeout rigido di 5 s via fork e SIGKILL; un timeout conta come fit non utilizzabile |
| GLMM Wald | `lme4::glmer`, Laplace (`nAGQ = 1`), ottimizzatore di default; effetti casuali `(1|id)` in M1-M2, `(1 + time|id)` in M3-M4, `(1 + time + dose|id)` in M5 |
| GLMM LRT | `glmer` con `bobyqa`, `maxfun = 1e5`; LRT con `anova()` contro il modello senza `time:trt`. Nei chunk attuali gli effetti casuali sono `(1|id)` in tutti gli scenari, anche in M3-M5 (rilievo di [protocollo datato (S1)](steps/S1.md)). **v4 ([LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md) (b))**: in M3-M5 il LRT usa la stessa struttura casuale del Wald, `(1 + time|id)` in M3-M4 e `(1 + time + dose|id)` in M5, e `glmer_lrt` si riesegue in [rerun correttivo (S3)](steps/S3.md) |
| GLM shrinkage | stima del GLMM riscalata di `1/sqrt(1 + c^2 D)`, `c = 16 sqrt(3)/(15 pi)`; SE per metodo delta con D fisso; p-value uguale al Wald del GLMM per invarianza |
| Test a inversione di segno (v1.7, dal run di [run con i metodi aggiunti (S7d)](steps/S7d.md)) | nullo `time:trt = 0` imposto dal modello ridotto `binY ~ time + trt` (correlazione di lavoro indipendente, `glm.fit`, `epsilon = 1e-12`, `maxit = 100`); contributi effettivi allo score per soggetto (Hemerik, Goeman, Finos 2020); statistica abs(sum_i g_i u_i), tutti i 2^N segni se N <= 12, altrimenti 1.999 segni casuali piu' l'identita', seme `20261022 + i` dopo la generazione; solo p-value. Principale: modello ridotto convergente con valori finiti; sensibilita': in piu' coefficienti sotto 10 in valore assoluto ([Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](decisions/D-47.md)) |
| GEE con correlazione di lavoro indipendente e AR(1) (v1.8, dal run di [run con i metodi aggiunti (S7d)](steps/S7d.md)) | `geepack::geeglm`, modello medio dello studio, `scale.fix = TRUE`, dati ordinati per soggetto e tempo; AR(1) sull'ordine delle visite; sandwich non corretto di geepack |
| Correzioni da implementazione indipendente (v1.8) | `fewclust::gee_vcov` (0.2.0, commit `8538cad` di lucavd/MMM, ramo package): Mancl-DeRouen, Kauermann-Carroll (radice principale di I - H_ii), Fay-Graubard (b = 0,75), sui fit exchangeable, indipendente e AR(1); `clubSandwich::vcovCR(type = "CR2")` sui fit exchangeable e indipendente. geesmv (exchangeable e indipendente, timeout 5 s) sui soli replicati 1-50, per la riproduzione delle righe dello studio e la coincidenza con fewclust; con AR(1) geesmv non e' utilizzabile (il pacchetto gee rifiuta AR-M con soggetti a una sola visita). Nota: il Kauermann-Carroll di geesmv 1.3 ricompone la radice di I - H_ii con la trasposta degli autovettori, valida solo per matrici simmetriche ([run con i metodi aggiunti (S7d)](steps/S7d.md)) |
| Modello ridotto del test a inversione di segno (v1.8) | modello medio della GEE dello scenario senza `time:trt`: `binY ~ time + trt`, piu' `age` in M2 e `dose` in M5 |

Software del run: R 4.6.0 (oggi 4.6.1: P5 riproduce i chunk bit per bit anche sotto 4.6.1), MMLB 1.0, lme4 2.0-6, geepack 1.3.13, geesmv 1.3, crew 1.3.2, mirai 2.7.2. Il `sessionInfo()` completo si allega al supplemento in [manoscritto (S7)](steps/S7.md).

Riferimento dei test di Wald: normale standard nelle tabelle primarie; distribuzione t con N - p gradi di liberta' (p = numero di coefficienti del modello medio: 4 in M1, M1c, M3, M4, M6, M7; 5 in M2 e M5) nella regola della Figura 1. **v1.6** ([rianalisi senza nuove simulazioni (S7b)](steps/S7b.md)): si aggiunge la t con N - 2 gradi di liberta' (Li e Redden 2015), ricalcolata da stima ed errore standard salvati, sulle stesse righe usabili; la regola della Figura 1 si riapplica con entrambe le scelte su tutti i 1.896 punti di disegno. Il test del rapporto di verosimiglianza resta sul chi quadrato a 1 grado di liberta'.

Congenialita' (D3), da dichiarare nel v4: MMM e' congeniale a GenBinaryY, il GLMM di Wald a GLMM-native, nessun metodo a GLMM-AR1; la GEE exchangeable non e' congeniale a nessuno.

## 4b. Applicazioni cliniche (fuori dagli obiettivi della simulazione; v1.6)

Le applicazioni non sono obiettivi del protocollo di simulazione (sez. 1) ma il loro disegno si fissa qui prima dei calcoli.

- **Coorte neonatale di ecografia polmonare** (manoscritto v4, sez. 2.7 e 3.7): codice `rev2/T3_rerun.R`, estrazioni `rev2/T3_draws.rds` (seme 20260807), misure `rev2/T5_analyse.R` e `rev2/T7_bias.R`.
- **Seconda applicazione, dataset `respiratory` di geepack** ([rianalisi senza nuove simulazioni (S7b)](steps/S7b.md)): 111 pazienti (centro x id), 4 visite ciascuno, esito binario (stato respiratorio buono), trattamento attivo contro placebo. Contrasto: `treatA` (riferimento placebo). Specificazione piena `outcome ~ treat + baseline + visit + age + sex + center`, ridotta `outcome ~ treat + baseline`. Strategie e opzioni come nella coorte neonatale (GEE exchangeable, correzioni Mancl-DeRouen, Kauermann-Carroll e Fay-Graubard di geesmv sulla stima della GEE con timeout di 5 s, GLMM a intercetta casuale con bobyqa, timeout rigidi), con una differenza dichiarata: il modello marginalizzato usa la specificazione dello scenario M1 (`lv.formula = ~1`, `t.formula = ~1`) perche' le visite sono una serie temporale. Sottocampioni di pazienti interi a N = 10, 15, 20, 25, 50, 500 estrazioni per N, seme 20260807 e stessa procedura di estrazione; riferimento = stima di ogni metodo sul dataset completo (le correzioni usano quella della GEE). Misure: quota di fit restituiti e usabili (stima e SE finiti, |stima| <= 10, SE < 10), rapporto SE del modello / deviazione standard delle stime, distorsione rispetto al riferimento, quota di estrazioni con una covariata binaria a un solo livello.

## 5. Performance measures (requisiti E1, E2)

Per cella (meccanismo x scenario x N x T x arm) e metodo, su `n_obs` replicati utilizzabili, con stima `θ̂_i`, SE del modello `s_i` e valore vero `θ`. Formule e MCSE di Morris (2019), Tabella 6:

| misura | stima | MCSE |
|---|---|---|
| distorsione | mean(θ̂_i) - θ | sqrt(empSE^2 / n) |
| SE empirico | empSE = sd(θ̂_i) | empSE / sqrt(2(n - 1)) |
| SE del modello | modSE = sqrt(mean(s_i^2)) | sqrt(var(s_i^2) / (4 n modSE^2)) |
| errore relativo di modSE | 100 (modSE/empSE - 1) | formula di Morris Tab. 6 |
| MSE | mean((θ̂_i - θ)^2) | sqrt(var((θ̂_i - θ)^2) / n) |
| copertura (95%) | mean(abs(θ̂_i - θ) <= 1,96 s_i) | sqrt(p(1 - p)/n) |
| copertura bias-eliminated | come sopra con mean(θ̂_i) al posto di θ | sqrt(p(1 - p)/n) |
| errore di I tipo (H0) | mean(p_i < 0,05) | sqrt(p(1 - p)/n) |
| potenza grezza (H1) | mean(p_i < 0,05) | sqrt(p(1 - p)/n) |
| potenza size-adjusted | rifiuto con il quantile 95% di abs(stat) sotto H0 della stessa cella | bootstrap per blocchi (sez. 7) |
| quota usabile | n_obs / replicati generati | sqrt(p(1 - p)/n) |

- **Usabilita'** (E2): `converged = TRUE` significa solo "nessun errore R"; i warning di lme4/MMLB stanno nella colonna `status`. [Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md) (a): primaria U0 (filtro di `agg_all5.R`, warning ammessi), sensibilita' U1 (esclusi tutti i warning); U2 (esclusi i soli warning di convergenza) tabulato. Quantificazione in [diagnostica sui risultati esistenti (S2)](steps/S2.md)(a), `tables/V4_S2_usable_pooled.csv`.
- **Denominatore** (E3, [Replicati con disegno non stimabile esclusi dal denominatore (D-11)](decisions/D-11.md)): replicati con disegno stimabile (entrambi gli arm osservati oltre il tempo 0). Le 247 chiavi per arm con disegno degenere si escludono e si dichiarano; sono le stesse nei tre meccanismi.
- **Scarto del valore nullo marginale sotto GLMM-native** (v1.6, [rianalisi senza nuove simulazioni (S7b)](steps/S7b.md)): per ogni punto di disegno GLMM-native, arm H0 e metodo marginale, con scarto delta = valore vero marginale della cella sotto H0 e s = SE empirico delle stime usabili, il rifiuto atteso in piu' e' Phi(-z + delta/s) + Phi(-z - delta/s) - 0,05, con z il quantile 0,975 della normale (e della t a N - p per il riferimento t). Si tabula in punti percentuali e si ricalcola la bandiera "chiaramente sopra 7,5%" dopo averlo sottratto dal tasso osservato.
- **Regole di reporting**: ogni numero con la sua MCSE (o la MCSE massima in didascalia); una differenza minore di 2 MCSE non si scrive come differenza; arrotondamento alla cifra giustificata dalla MCSE (Morris sez. 6.2).

## 6. Numero di replicati (requisito F1; [Replicati per cella: restano 1.000 (D-03)](decisions/D-03.md) (a))

n_sim = 1.000 per cella, mantenuto ([Replicati per cella: restano 1.000 (D-03)](decisions/D-03.md) (a)) e giustificato dalla MCSE che ne deriva (sqrt(p(1 - p)/n)):

- errore di I tipo e copertura al livello nominale: 0,69 punti percentuali con 1.000 replicati utilizzabili, 0,97 con 500 (quota usabile del 50%);
- potenza: al massimo 1,58 punti (p = 0,5), 1,26 a p = 0,8;
- servirebbero 1.900 replicati per una MCSE di 0,5 punti al 5%.

Una cella risolve quindi una deviazione dal 5% nominale di circa 1,4 punti (2 MCSE). Se [inferenza per blocchi (S5)](steps/S5.md) trova una conclusione del paper che si regge su una differenza entro 2 MCSE, si aggiungono replicati solo nelle celle interessate ([Replicati per cella: restano 1.000 (D-03)](decisions/D-03.md) (b), con lo stesso schema dei semi).

## 7. Numeri casuali e dipendenza fra dataset (requisiti F2, F3, F4; [Stesso seme fra celle, inferenza per blocchi di replicato (D-05)](decisions/D-05.md) (a))

- **Seme**: `set.seed(20260805L + 1000L * model + i)` per il replicato `i` dello scenario `model`, uguale per arm, N, T e meccanismo.
- **Generatore**: `RNGkind` = **L'Ecuyer-CMRG, Inversion, Rejection**. `run_full.R` non lo imposta: lo impostano i daemon di mirai (2.7.2, via crew 1.3.2) in ogni worker. Verificato in [protocollo datato (S1)](steps/S1.md) leggendo `RNGkind()` dentro un worker crew e riproducendo i chunk (P5: scarto 0 con L'Ecuyer-CMRG, 0/30 con Mersenne-Twister). La voce dell'audit [F2: Seme e generatore riportati](audit/requirements/F2.md) va corretta di conseguenza.
- **Conseguenza per [rerun correttivo (S3)](steps/S3.md) ([Generatore di numeri casuali fissato nei rerun (D-07)](decisions/D-07.md))**: i rerun devono usare lo stesso RNGkind, altrimenti i chunk rigenerati perdono l'accoppiamento CRN anche a parita' di seme. In [rerun correttivo (S3)](steps/S3.md) il worker chiama esplicitamente `RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")` prima di `set.seed()`, cosi' il risultato non dipende dal default di mirai.
- **Disegno CRN** ([Stesso seme fra celle, inferenza per blocchi di replicato (D-05)](decisions/D-05.md) (a)): stesso seme per H0 e [H1: Metodi in ordine ADEMP / tabella ADEMP](audit/requirements/H1.md), per tutte le N e T e per i tre meccanismi. In H0 e [H1: Metodi in ordine ADEMP / tabella ADEMP](audit/requirements/H1.md) il disegno e le estrazioni coincidono; GLMM-native e GLMM-AR1 condividono disegno ed effetti casuali in M1-M4; celle con N e T diverse condividono il prefisso dello stream. Il disegno e' tenuto (Morris sez. 5.4) e dichiarato.
- **Inferenza fra celle**: le procedure che trattano le celle come indipendenti si sostituiscono con un bootstrap per blocchi di replicato, che ricampiona gli indici `(model, i)` congiuntamente su tutte le celle (B = 2000, seed 20261015); il BH si affianca con il BY ([inferenza per blocchi (S5)](steps/S5.md)).
- **Parallelizzazione** (F4): chunk di 100 replicati con semi deterministici per replicato; ogni replicato si riproduce da solo (P5), gli stati RNG non sono salvati. Da dichiarare.

## 8. Provenienza e verifiche di S1 (requisiti H5, G1)

| chunk di `results_full/` | data | codice |
|---|---|---|
| M1-M2, tre meccanismi (11.520) | 2026-08-05/06 | job `23035d84`, copia in `_superseded/rundir_m1m2_job23035d84/` |
| M3-M5, tre meccanismi (17.280) | 2026-08-08/09 | commit `5f84297`, copia eseguita in `.rundir/` (= `R/` attuale, md5) |

Verifiche (report `reports/morris/S1.html`): la tabella coincide con il codice (P1, P3), la parita' H0/H1 si legge dalla tabella (P2), il codice M1-M2 e' equivalente (P4), i chunk si riproducono bit per bit (P5), i chunk grezzi non sono stati toccati (P6).

## 9. Pre-specificazione (requisito H6)

La regola decisionale di sez. 4.2 del paper ("unconditional validity") compare nel git per la prima volta il 2026-08-10 (`55a6363`), dopo l'avvio del run (2026-08-05). Salvo evidenza datata contraria fornita dal PI entro [manoscritto (S7)](steps/S7.md), il v4 la dichiara **definita dopo aver visto i risultati**. Da questo protocollo in poi la pre-specificazione e' documentata: attese datate e committate prima di ogni calcolo (`docs/sessions/`), decisioni datate (`DECISIONS.md`).

## 10. Decisioni (stato al 2026-10-15)

[Dati sotto H0 di GenBinaryY rigenerati in M3-M5 (D-01)](decisions/D-01.md) (a), [Varianza della pendenza su dose allineata in M5 (D-02)](decisions/D-02.md) (a), [Replicati per cella: restano 1.000 (D-03)](decisions/D-03.md) (a), [Stesso seme fra celle, inferenza per blocchi di replicato (D-05)](decisions/D-05.md) (a), [Codice pubblico alla sottomissione, con DOI (D-06)](decisions/D-06.md) si, [LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md) (b), [Valori del modello simulato confrontati con la coorte LUS (D-10)](decisions/D-10.md): decise dal PI il 2026-10-15. [Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md): dopo [diagnostica sui risultati esistenti (S2)](steps/S2.md) (per roadmap). [Generatore di numeri casuali fissato nei rerun (D-07)](decisions/D-07.md) (RNGkind esplicito in [rerun correttivo (S3)](steps/S3.md)): conseguenza della regola "schema dei semi invariato". [Nessuno step nuovo finche' restano riserve aperte (D-09)](decisions/D-09.md): regola di sessione (nessuna nuova chat con riserve aperte). Dettagli in `DECISIONS.md`.

## 11. Modifiche al protocollo

| versione | data | commit | modifica |
|---|---|---|---|
| 1.0 | 2026-10-15 | `cc534d6` | prima versione |
| 1.1 | 2026-10-15 | `c14dee4` | sez. 4: LRT allineato alle pendenze casuali in M3-M5 ([LRT del GLMM con le pendenze casuali del Wald (D-08)](decisions/D-08.md) (b)); sez. 2.4: valori di effetto, intercetta e varianza ancorati alla coorte LUS ([Valori del modello simulato confrontati con la coorte LUS (D-10)](decisions/D-10.md)); sez. 10 aggiornata. Nessun calcolo e' stato eseguito fra 1.0 e 1.1. |
| 1.2 | 2026-10-16 | `aa1c4c0` | sez. 3: verita' marginale GLMM-native per quadratura di Gauss-Hermite esatta (z'b normale con varianza z'Dz) al posto del Monte Carlo a `colMeans(Z)`, definita come limite della GEE di indipendenza pesato col disegno realizzato; sez. 2.4: [Confronto con la coorte attraverso grandezze analoghe (D-10b)](decisions/D-10b.md), confronto con la coorte per analoghi. Scritto all'apertura di [diagnostica sui risultati esistenti (S2)](steps/S2.md), prima di ogni calcolo dello step. |
| 1.3 | 2026-10-16 | `a4805f2` | sez. 5: [Fit con warning tenuti nell'analisi principale (D-04)](decisions/D-04.md) (a) e [Replicati con disegno non stimabile esclusi dal denominatore (D-11)](decisions/D-11.md) (denominatore); sez. 3: fonti dei valori veri e verita' condizionale AR1 per replicato ([Valori veri ricalcolati, tabelle precedenti superate (D-12)](decisions/D-12.md)). Scritto alla chiusura di [diagnostica sui risultati esistenti (S2)](steps/S2.md), sulla base dei suoi risultati e prima di [rerun correttivo (S3)](steps/S3.md) e [riaggregazione con MCSE (S4)](steps/S4.md). |
| 1.4 | 2026-10-16 | `5f009ad` | sez. 2.5: scenario ancorato alla coorte M1c ([Nuovo scenario ancorato alla coorte (M1c) (D-13)](decisions/D-13.md)). Decisione presa dopo i risultati di [diagnostica sui risultati esistenti (S2)](steps/S2.md)(g) e prima di qualsiasi calcolo di [rerun correttivo (S3)](steps/S3.md). |
| 1.5 | 2026-10-16 | questo commit | sez. 2.5: parametri di M1c ([Parametri dello scenario ancorato alla coorte (D-14)](decisions/D-14.md)); struttura dei chunk di [rerun correttivo (S3)](steps/S3.md) ([Organizzazione dei file del rerun correttivo (D-15)](decisions/D-15.md)). Scritto all'apertura di [rerun correttivo (S3)](steps/S3.md), prima di ogni calcolo dello step. |
| 1.6 | 2026-10-21 | questo commit | sez. 4: t con N - 2 gradi di liberta' accanto a N - p; sez. 4b (nuova): disegno delle due applicazioni cliniche, con la seconda su `respiratory` di geepack; sez. 5: rifiuto atteso dovuto allo scarto del nullo marginale sotto GLMM-native. Scritto all'apertura di [rianalisi senza nuove simulazioni (S7b)](steps/S7b.md), prima di ogni calcolo dello step. |
| 1.7 | 2026-10-22 | questo commit | sez. 4: test a inversione di segno sui contributi allo score per il run su tutti i punti di disegno, scelto dal PI alla chiusura di [pilota dei test a livello di soggetto (S7c)](steps/S7c.md) ([Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](decisions/D-46.md), [Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](decisions/D-47.md)). Scritto prima di ogni calcolo di [run con i metodi aggiunti (S7d)](steps/S7d.md). |
| 1.8 | 2026-10-22 | questo commit | sez. 4: GEE con correlazione indipendente e AR(1), correzioni di fewclust 0.2.0 e CR2, geesmv sui replicati 1-50, modello ridotto del test a inversione di segno con le covariate di M2 e M5 ([geesmv solo sulle prime 50 simulazioni di ogni scenario nel run con i metodi aggiunti (D-48)](decisions/D-48.md), [fewclust 0.2.0 registrato dall'agente nel ramo package cosi' com'era sul portatile del PI (D-49)](decisions/D-49.md)). Scritto prima di ogni calcolo del run di [run con i metodi aggiunti (S7d)](steps/S7d.md). |
