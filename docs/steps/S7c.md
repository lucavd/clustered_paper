---
type: Step
title: "Disegno e pilota dei test a livello di soggetto e di permutazione"
description: "Definire analisi a livello di soggetto e test di permutazione o bootstrap il cui nullo coincida con quello dell'interazione trattamento per tempo, provarle su poche celle e misurarne i tempi per organizzare il run."
step_id: "S7c"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-10-22T14:03:10+02:00, evidence: "chat del 2026-10-22: \"1 e 2 come consigli tu\" (D-46, D-47); commit S7c: closed, approved by PI" }
resource: "../../reports/morris/S7c.html"
generated: { by: mmm-agent/skill-v4, at: 2026-10-21T21:46:43+02:00 }
sources:
  - id: d42
    resource: "../decisions/D-42.md"
    title: "Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale"
  - id: anderson2001
    resource: https://doi.org/10.1139/f01-004
    title: "Anderson (2001). Permutation tests for univariate or multivariate analysis of variance and regression. Canadian Journal of Fisheries and Aquatic Sciences 58:626-639"
  - id: buzkova2016
    resource: https://doi.org/10.1515/em-2015-0010
    title: "Buzkova (2016). Interaction testing: residuals-based permutations and parametric bootstrap in continuous, count, and binary data. Epidemiologic Methods"
  - id: klinesantos2012
    resource: https://doi.org/10.1515/2156-6674.1006
    title: "Kline, Santos (2012). A score based approach to wild bootstrap inference. Journal of Econometric Methods 1:23-41"
  - id: hemerik2020
    resource: https://doi.org/10.1111/rssb.12369
    title: "Hemerik, Goeman, Finos (2020). Robust testing in generalized linear models by sign flipping score contributions. JRSS-B 82:841-864"
  - id: andreella2025
    resource: https://doi.org/10.1017/psy.2024.22
    title: "Andreella, Goeman, Hemerik, Finos (2025). Robust inference for generalized linear mixed models: a two-stage summary statistics approach based on score sign flipping. Psychometrika 90:531-553"
  - id: leyrat2018
    resource: https://doi.org/10.1093/ije/dyx169
    title: "Leyrat et al. (2018). Cluster randomized trials with a small number of clusters: which analyses should be used? International Journal of Epidemiology"
---
# Dubbi ricevuti

1. **Il nullo dei test a livello di soggetto e di permutazione non e' quello del paper.** Ricevuto il 2026-10-21, dal piano approvato dal PI ([Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](../decisions/D-42.md)). L'estimando del paper e' l'interazione trattamento per tempo sulla scala logit, e sotto H0 i due bracci differiscono ancora per l'effetto principale del trattamento (0,2 sulla scala logit). Quindi:
   - permutare l'assegnazione al trattamento non e' un test esatto per l'interazione, perche' sotto H0 i bracci non sono scambiabili;
   - un riassunto per soggetto (per esempio la pendenza dell'esito sul tempo) seguito da un confronto fra bracci verifica un nullo sulla scala della probabilita', che sotto H0 logit non e' esattamente zero; inoltre i soggetti con una sola osservazione (29-31% con T = 2) non danno una pendenza.
   Lo step deve scegliere test il cui nullo coincida con quello del paper, oppure misurare quanto se ne discosta e dichiararlo. Candidati da verificare in letteratura (Zotero, PubMed) prima di scegliere: test di permutazione su residui sotto il modello nullo, bootstrap selvaggio per cluster con il nullo imposto, analisi a livello di soggetto come in Leyrat e colleghi (riferimento [21] del manoscritto v4).

# Piano

- Pilota su poche celle (per esempio N 10 e 25, T 2, 3 e 10, tre meccanismi, M1, H0 e H1): scarto del nullo, falsi positivi, quota di risultati usabili, tempo per dataset. Nessuna misura entra nel paper.
- Tempi, sugli stessi dataset, solo per organizzare il run sul server (non limitano le scelte): GEE con correlazione di lavoro indipendente e AR(1), correzioni con geesmv, con fewclust e CR2 di clubSandwich, test scelti qui.
- Al PI: al massimo 3 opzioni su quali test entrano nel run, con cosa cambierebbe nel paper.
- Prodotti: report `reports/morris/S7c.qmd`; nuova versione datata di `docs/PROTOCOL_ADEMP.md` con i metodi scelti.

# Disegno dei test (fissato prima di ogni calcolo, 2026-10-22)

Letteratura consultata (Zotero e ricerca bibliografica, 2026-10-22). Un test di permutazione esatto per un'interazione in generale non esiste, perche' sotto il nullo dell'interazione le unita' non sono scambiabili se gli effetti principali non sono nulli[^anderson2001][^buzkova2016]; la via standard e' imporre il nullo con un modello ridotto e ricampionare cio' che sotto il nullo ha media zero (residui o contributi allo score)[^buzkova2016][^klinesantos2012][^hemerik2020]. Il test a sign-flip dei contributi allo score per soggetto e' asintoticamente esatto, robusto a una varianza sbagliata e gestisce i dati longitudinali passando l'identificativo del soggetto[^hemerik2020][^andreella2025]; il bootstrap selvaggio degli score con il nullo imposto ne e' l'equivalente econometrico per dati a gruppi[^klinesantos2012]. L'analisi a livello di gruppo in due stadi (modello ridotto sulle osservazioni, riassunto dei residui per gruppo, confronto fra bracci con una t) mantiene l'errore di I tipo anche con pochi gruppi[^leyrat2018].

Il nullo e' quello del paper: interazione `time:trt` uguale a zero sulla scala logit marginale, con intercetta, `time` e `trt` liberi. Tutti i test candidati lo impongono con lo stesso modello ridotto: regressione logistica `binY ~ time + trt` con correlazione di lavoro indipendente (stimata con `glm.fit`), cioe' la GEE di indipendenza sotto H0. Da questo fit:

- contributo effettivo allo score del soggetto i per il coefficiente di `time:trt`: u_i = sum_j x_ij (y_ij - mu_ij) - I_{x,nu} I_{nu,nu}^{-1} sum_j z_ij (y_ij - mu_ij), con x = time * trt, z = (1, time, trt) e I l'informazione di Fisher del modello logistico[^hemerik2020];
- riassunto per soggetto dei residui: r_i = sum_j time_ij (y_ij - mu_ij), la tendenza nel tempo dei residui del modello ridotto. Poiche' il modello ridotto contiene `time`, sum_i r_i = 0, e la differenza fra le medie dei bracci e' proporzionale allo score di `time:trt`: S = n1 n0 / N * (media_1 - media_0). Un soggetto con una sola osservazione al tempo 0 ha r_i = 0 e resta nell'analisi.

| test | nome nel codice | statistica e riferimento | nullo |
|---|---|---|---|
| A | `flip` | sign-flip dei contributi effettivi u_i per soggetto (Hemerik e colleghi), statistica abs(sum_i g_i u_i); tutti i 2^N segni se N <= 12, altrimenti 1.999 segni casuali piu' l'identita' | quello del paper, imposto |
| B | `subj_t` | analisi a livello di soggetto in due stadi: t di Student a varianza comune sui riassunti r_i fra i bracci, N - 2 gradi di liberta' (la scelta di Li e Redden) | quello del paper, imposto |
| C | `subj_perm` | permutazione dei bracci fra i riassunti r_i, statistica t di Welch (permutazione studentizzata, valida anche con varianze diverse fra bracci); tutte le combinazioni se sono al piu' 1.999, altrimenti 1.999 permutazioni casuali piu' l'identita' | quello del paper, imposto |
| D | `naive_perm` (confronto) | permutazione dell'assegnazione al trattamento per soggetto, rifit del modello completo con correlazione indipendente, statistica abs(z) con sandwich per soggetto; 999 permutazioni piu' l'identita' | non esatto: sotto H0 i bracci differiscono per l'effetto principale (0,2) |
| E | `naive_slope` (confronto) | pendenza OLS di y sul tempo per ogni soggetto con almeno due osservazioni, t di Student fra i bracci con n_s - 2 gradi di liberta' | scala della probabilita': sotto H0 logit non e' esattamente zero |

Risultato usabile (decisione tecnica, stessa logica del filtro dei fit del paper): A, B, C richiedono un fit ridotto convergente con coefficienti finiti e in valore assoluto sotto 10, almeno due soggetti per braccio e varianza dei contributi positiva; D richiede la statistica osservata finita (le permutazioni con statistica non finita si escludono dal riferimento e si contano); E richiede almeno due soggetti con pendenza per braccio e varianza positiva. I test danno solo un p-value: nessuna stima, quindi niente distorsione ne' copertura.

Semi: il dataset di ogni replicato si rigenera esattamente come nel run (seme `20260805 + 1000 * 1 + i`, L'Ecuyer-CMRG, indicatore di disegno degenere prima della generazione); i ricampionamenti usano poi `set.seed(20261022 + i)`, uguale fra celle, dopo la generazione.

# Piano operativo

- Celle del pilota: scenario M1, N in {10, 25}, T in {2, 3, 10}, tre meccanismi, H0 e H1: 18 punti di disegno, 36 celle, 1.000 replicati ciascuna (36.000 dataset). Denominatore: replicati con disegno stimabile ([Replicati con disegno non stimabile esclusi dal denominatore (D-11)](../decisions/D-11.md)).
- Confronto appaiato con le strategie del paper sugli stessi replicati, dai risultati salvati in `results_full/`: GEE exchangeable con Fay e Graubard e la t a N - p, modello marginalizzato con la t a N - p (in M1 il controllo di convergenza non cambia il fit, perche' il modello ha gia' `lv.formula = ~1`).
- Tempi per organizzare il run successivo (solo organizzazione, nessuna misura entra nel paper): 20 replicati per ciascuno dei 18 punti di disegno sotto H0 piu' N = 145, T = 10 nei tre meccanismi; tempo seriale per dataset di GEE exchangeable, indipendente e AR(1); correzioni di geesmv (con il limite di 5 s) per le tre correlazioni; correzioni di fewclust (`gee_vcov`, versione 0.1.0 dal commit `5103bb3` di `~/MMM-package`, libreria separata: non e' ancora la versione allineata al paper, ma le formule delle correzioni servono qui solo per i tempi); CR2 di clubSandwich; i test A-E.
- Codice in `R/s7c/`, risultati in `results_v4/S7c/`, tabelle `tables/V4_S7c_*.csv`. Al piu' 20 processi (limite di temperatura dell'host).

# Attese pre-registrate (scritte prima di ogni calcolo)

## Cancelli

| codice | controllo | passa se |
|---|---|---|
| GD1 | integrita' dei risultati grezzi | i md5 di `results_full/` prima e dopo lo step coincidono fra loro e con `results_v4/S7b/md5_results_full_after.txt` |
| GD2 | parita' H0/H1 | `R/testing/test_h0h1_parity.R` passa prima del run |
| GD3 | stessi dataset del run | sui dataset rigenerati la GEE exchangeable riproduce stima ed errore standard salvati in `results_full/` entro 1e-8 in ogni replicato del pilota in cui la riga salvata e' finita |
| GD4 | controprova delle formule | su 50 dataset: la somma dei contributi effettivi coincide con la derivata numerica (numDeriv) della log-verosimiglianza rispetto al coefficiente di `time:trt` nel fit ridotto (errore relativo sotto 1e-6); S = n1 n0 / N * (media_1 - media_0) entro 1e-10; la t di B coincide con `t.test(var.equal = TRUE)` e la statistica di C con `t.test()` di Welch entro 1e-10 |
| GD5 | completezza | 36.000 dataset, ciascuno con una riga per ognuno dei 5 test (p-value oppure motivo dichiarato di non usabilita') |

Controllo in piu', riportato ma non cancello: sullo stesso insieme di segni, il p-value di A coincide con `flipscores::flipscores(score_type = "effective", id = )`.

## Previsioni

Margini della skill: errore di I tipo da 2,5% a 7,5%; differenze fra test +/- 1 punto in media sulle celle e +/- 2,5 punti per cella. Per la potenza il margine e' +/- 5 punti in media sulle celle: una differenza piu' grande cambierebbe la frase del paper su quale strategia conviene con pochi soggetti. L'aggregato e' la media sulle 18 celle sotto H0 (o H1), con intervallo dal bootstrap per blocchi di replicato (B = 2000, seed 20261015). Le bandiere "chiaramente sopra / sotto" contengono gia' la MCSE e si contano come interi.

| codice | quantita' | atteso | intervallo di rilevanza | perche' conta |
|---|---|---|---|---|
| PD1 | errore di I tipo di A, media sulle 18 celle H0 | 4,5% (con 10 soggetti il riferimento e' discreto e un po' conservativo) | da 2,5% a 7,5% | A entra nel run solo se tarato |
| PD2 | errore di I tipo di B, media sulle 18 celle H0 | 5% | da 2,5% a 7,5% | idem per B |
| PD3 | errore di I tipo di C, media sulle 18 celle H0 | 4,5% | da 2,5% a 7,5% | idem per C |
| PD4 | celle H0 in cui A, B o C sono chiaramente sopra 7,5% | 0 per ciascun test | 0 (una cella basta a escludere un test dalla regola della Figura 1) | la regola chiede che una strategia non sia mai chiaramente sopra |
| PD5 | quota usabile di A, B, C meno quella di Fay e Graubard, per cella | sempre >= 0; da +5 a +20 punti con 10 soggetti | per cella non sotto -2,5 punti | i test che non stimano il modello completo dovrebbero rispondere piu' spesso con pochi soggetti |
| PD6 | errore di I tipo di D meno quello di C, appaiato, media sulle 18 celle H0 | 0 punti (l'effetto principale 0,2 e' piccolo, la non scambiabilita' sposta poco) | da -1 a +1 punto | misura quanto pesa la non esattezza della permutazione dell'assegnazione |
| PD7 | quota di soggetti senza pendenza in E con T = 2 (sui 6 punti con T = 2) | 30% | da 25% a 35% | E perde quasi un terzo dei soggetti proprio dove servirebbe |
| PD8 | rifiuto atteso in piu' di E dovuto al nullo sulla scala della probabilita', Phi(-z + delta/s) + Phi(-z - delta/s) - 0,05 con delta = media sui replicati della differenza fra bracci delle pendenze medie e s la sua deviazione standard, massimo sulle 18 celle H0 | sotto 0,5 punti | da 0 a 1 punto | se resta sotto 1 punto lo scarto del nullo di E e' trascurabile e conta solo la perdita di soggetti |
| PD9 | potenza di A, B, C meno quella della GEE con Fay e Graubard e la t, appaiata sui replicati usabili per entrambi, media sulle 18 celle H1 | -3 punti (correlazione di lavoro indipendente e test dello score contro Wald exchangeable) | da -5 a +5 punti | se un test tarato avesse molta meno potenza non varrebbe la pena raccomandarlo |

[^anderson2001]: Anderson (2001), Canadian Journal of Fisheries and Aquatic Sciences 58:626-639.
[^buzkova2016]: Bůžková (2016), Interaction testing: residuals-based permutations and parametric bootstrap in continuous, count, and binary data, Epidemiologic Methods, doi:10.1515/em-2015-0010.
[^klinesantos2012]: Kline, Santos (2012), A score based approach to wild bootstrap inference, Journal of Econometric Methods 1:23-41, doi:10.1515/2156-6674.1006.
[^hemerik2020]: Hemerik, Goeman, Finos (2020), Robust testing in generalized linear models by sign flipping score contributions, JRSS-B 82:841-864, doi:10.1111/rssb.12369.
[^andreella2025]: Andreella, Goeman, Hemerik, Finos (2025), Robust inference for generalized linear mixed models: a two-stage summary statistics approach based on score sign flipping, Psychometrika 90:531-553, doi:10.1017/psy.2024.22.
[^leyrat2018]: Leyrat e colleghi (2018), Cluster randomized trials with a small number of clusters: which analyses should be used?, International Journal of Epidemiology 47:321-331, doi:10.1093/ije/dyx169.

# Esito del calcolo (2026-10-22)

Run: `R/s7c/s7c_run.R` (36.000 dataset, 20 processi, circa 41 minuti), `R/s7c/s7c_agg.R`; analisi aggiunte `R/s7c/s7c_posthoc.R` e `R/s7c/s7c_posthoc_agg.R`; tempi `R/s7c/s7c_timing.R` (417 dataset, un processo). Tabelle `tables/V4_S7c_*.csv`.

## Cancelli

| codice | esito | dettaglio |
|---|---|---|
| GD1 | passa | 28.800 file di `results_full/`, md5 prima = dopo = lista di S7b |
| GD2 | passa | parita' H0/H1: griglia principale 15/15, M1c 3/3 |
| GD3 | passa | 33.448 replicati con GEE salvata finita: stima ed errore standard identici (scarto massimo 0) |
| GD4 | fallito come scritto, spiegato | score effettivo contro numDeriv entro 1,9e-10 e t contro `t.test` entro 1,3e-15 (passano); l'identita' S = n1 n0 / N (media_1 - media_0) ha scarto massimo 1,97e-10 contro 1e-10 richiesto. L'identita' vale solo al massimo esatto: al massimo calcolato S - n1/N sum_i r_i = n1 n0 / N (media_1 - media_0) tiene entro 1,8e-15, e sum_i r_i e' il residuo dell'equazione dello score del tempo a convergenza (al piu' 3,3e-10). La tolleranza era piu' stretta della precisione del fit; la formula e' giusta e nessun risultato ne dipende |
| GD5 | passa | 36.000 dataset; 72 errori di generazione, tutti con disegno degenere (esclusi dal denominatore come nel run); 5 righe di test per ognuno dei 35.928 generati |

Controllo in piu' (non cancello): `flipscores` 1.3.2 con `score_type = "effective"` e `id` non coincide con A (scarto massimo dei p-value 0,04 su 19 dataset con tutti i 1.024 segni). Non e' la stessa statistica: nel sorgente il pacchetto proietta senza i pesi di Fisher e usa i residui di Pearson; A segue la formula dell'articolo di Hemerik e colleghi ed e' controllata con numDeriv.

## Previsioni

| codice | atteso | intervallo | osservato (intervallo per blocchi, o +/- 2 MCSE) | esito | spiegazione |
|---|---|---|---|---|---|
| PD1 | 4,5% | 2,5-7,5% | A: 5,19% (4,80-5,63) | confermata | |
| PD2 | 5% | 2,5-7,5% | B: 0,51% (0,34-0,70) | smentita | il riassunto r_i porta il livello di ogni soggetto: la varianza della t e' da 2,1 a 8,1 volte quella vera (`V4_S7c_variance_ratio.csv`) |
| PD3 | 4,5% | 2,5-7,5% | C: 1,93% (1,63-2,24) | smentita | stessa causa di PD2 |
| PD4 | 0 celle | 0 | A 0, B 0, C 0 | confermata | |
| PD5 | >= 0 | non sotto -2,5 punti per cella | differenza minima 0 punti (A, B, C) | confermata | A risponde in 81-86% dei dataset con 25 soggetti e 2 visite, Fay e Graubard in 36-40% |
| PD6 | 0 punti | -1/+1 | D - C: +3,71 (3,21-4,21) | smentita | conseguenza di PD3: C e' sotto il nominale; D da solo e' 4,55% (4,24-4,87), mai fuori taratura in nessuna cella (da 2,9% a 5,5%) |
| PD7 | 30% | 25-35% | 29,8% | confermata | |
| PD8 | < 0,5 punti | 0-1 punto | 0,16 (da -0,30 a 0,62) | indeterminata | a cavallo solo dell'estremo inferiore (lo scarto vero non puo' essere negativo); l'estremo superiore e' sotto 1 punto, quindi la conclusione "trascurabile" non cambia |
| PD9 | -3 punti | -5/+5 | A: +0,81 (0,43-1,20); B: -8,08; C: -8,08 | A confermata; B, C smentite | B e C non hanno potenza per la stessa causa di PD2 |

## Analisi aggiunte dopo i risultati (dichiarate, esplorative)

Aggiunte per spiegare e correggere PD2-PD3 e per misurare il peso del filtro sui coefficienti; si confermano solo nel run successivo, con attese scritte prima. Stessi 36.000 dataset, stessi semi.

| variante | errore di I tipo medio (18 celle) | celle chiaramente sopra / sotto | quota usabile | potenza meno Fay e Graubard (appaiata) |
|---|---|---|---|---|
| A senza filtro sui coefficienti (fit ridotto convergente con valori finiti) | 4,15% (3,87-4,47) | 0 / 1 (1,5% con 10 soggetti, 2 visite, dati dal modello misto) | 99-100% | +0,81 (0,43-1,20) |
| t robusta dello score, N - 2 | 1,41% (1,23-1,57) | 0 / 12 | 99-100% | -3,38 |
| t a livello di soggetto con il tempo centrato nel soggetto, N - 2 | 4,89% (4,56-5,23) | 0 / 0; massimo 9,45% (limite inferiore 7,45) con 10 soggetti, 2 visite, dati sulla scala della popolazione | 87-100% | -0,30 (-0,73 a 0,18) |
| permutazione della stessa quantita', Welch | 6,07% (5,70-6,44) | 2 / 0 (11,2% e 9,9% con 10 soggetti e 3 visite) | 87-100% | +0,83 |

Spiegazione di PD2-PD3 con il calcolo: il rapporto fra la varianza che la t attribuisce alla differenza fra bracci e la varianza empirica di quella differenza sotto H0 va da 2,10 a 8,06 per r_i e da 1,04 a 1,58 per il riassunto centrato nel soggetto (`tables/V4_S7c_variance_ratio.csv`). Il livello del soggetto entra in r_i con peso pari alla somma dei suoi tempi; nella somma sul braccio si annulla, perche' il modello ridotto impone sum_{trt=1} e_ij = 0, ma resta nella varianza fra soggetti.

## Decisioni tecniche (agente)

- Modello ridotto stimato con `glm.fit` e `epsilon = 1e-12`, `maxit = 100`, fissato prima del run dopo il primo tentativo di GD4 (con la tolleranza di default lo scarto dell'identita' era 1,2e-7): migliora la precisione del fit, non tocca alcuna misura.
- `flipscores` 1.3.2 e `fewclust` 0.1.0 (commit `5103bb3` di `~/MMM-package`) in una libreria separata `~/R/s7c_lib`, fuori dalla libreria del progetto.
- Tempi misurati in serie (un processo) perche' le misure dentro processi biforcati sono risultate inaffidabili su questo host.

# In breve

## In breve: quali test, con il valore nullo del paper, vale la pena aggiungere al run su tutte le situazioni?

**A cosa serviva.** Il paper confronta strategie che stimano un modello e poi lo usano per il test. Volevamo provare test che non richiedono di stimare il modello completo e il cui valore nullo è esattamente quello del paper: nessuna differenza fra i bracci nell'andamento nel tempo, sulla scala logit, con livello ed effetto principale del trattamento liberi. Li abbiamo provati in 18 situazioni (10 e 25 soggetti, 2, 3 e 10 visite, i tre meccanismi generatori), con 1.000 simulazioni per scenario sia senza sia con effetto del trattamento, sugli stessi dataset del paper, così il confronto con le sue strategie è appaiato. Nessun numero di questa prova entra nel paper: serve a scegliere che cosa calcolare su tutte le 1.896 situazioni.

**Cosa abbiamo trovato.**

1. **Il test a segni alterni sui contributi dei soggetti allo score (Hemerik, Goeman e Finos) è tarato, risponde più spesso della GEE con la correzione di Fay e Graubard e ha almeno la stessa potenza.**
   - In media 5,2 falsi positivi su 100; in nessuna delle 18 situazioni è fuori taratura (il massimo è 6,7).
   - Con 25 soggetti e 2 visite dà un risultato in circa 8 dataset su 10, la correzione di Fay e Graubard in circa 4 su 10; con 10 soggetti e 2 visite in 3 su 10 contro meno di 1 su 10.
   - Sugli stessi dataset trova un effetto vero un po' più spesso: 0,8 test su 100 in più, in media.
   - Senza il filtro sui coefficienti del modello ridotto (lo stesso filtro che il paper applica ai fit) risponde in quasi tutti i dataset, con 4,2 falsi positivi su 100 in media; è prudente solo dove i dati non contengono quasi informazione (1,5 su 100 con 10 soggetti e 2 visite).
   - Dà solo un p-value: nessuna stima e nessun intervallo.
2. **Le due analisi a livello di soggetto scritte prima del calcolo non funzionano, e abbiamo capito perché.**
   - Il confronto fra bracci della tendenza dei residui di ciascun soggetto, con la t o con una permutazione, dà 0,5 e 1,9 falsi positivi su 100 e non trova quasi mai un effetto vero (meno di 2 volte su 100).
   - Quella tendenza porta con sé il livello del soggetto: nel test complessivo si annulla, ma gonfia la varianza fra soggetti da 2 a 8 volte.
   - Centrando il tempo dentro ogni soggetto la varianza torna giusta. Con questa correzione, introdotta dopo aver visto i risultati, la t dà 4,9 falsi positivi su 100 in media e la stessa potenza di Fay e Graubard, ma con 10 soggetti e 2 visite arriva a 9,5 su 100: non è fuori taratura per meno di un decimo di punto. La stessa quantità con la permutazione è fuori taratura in 2 situazioni su 18.
3. **Permutare l'assegnazione al trattamento non è esatto per l'interazione, ma qui sbaglia poco: 4,6 falsi positivi su 100 in media, mai fuori taratura.** Richiede però mille stime del modello per dataset: su tutta la griglia costerebbe circa un terzo del tempo dell'intero run.
4. **La pendenza di ciascun soggetto sulla scala della probabilità ha un valore nullo diverso da quello del paper, ma lo scarto pesa al più 0,6 falsi positivi su 100; il suo problema è che con 2 visite il 30% dei soggetti non ha una pendenza.**
5. **Un controllo delle formule è fallito alla lettera, per una tolleranza scritta più stretta della precisione del fit: la formula è giusta e nessun risultato ne è toccato.**

**Cosa devi decidere.**

1. **Quali test aggiungiamo al run su tutte le 1.896 situazioni? Consiglio il solo test a segni alterni.**
   - È l'unico con il valore nullo del paper, tarato in tutte le situazioni provate, che risponde più spesso delle strategie del paper. Su tutta la griglia costa circa 4 ore di calcolo.
   - Alternativa: aggiungere anche la t a livello di soggetto con il tempo centrato. È l'analisi in due stadi che il paper cita (Leyrat e colleghi) e costa poco, ma è nata dopo aver visto questa prova ed è al limite con 10 soggetti e 2 visite: la regola della Figura 1 potrebbe escluderla.
   - Alternativa: aggiungere anche la permutazione dell'assegnazione al trattamento, per documentare su tutta la griglia quanto pesa la sua non esattezza. Circa 1.070 ore di calcolo in più, quasi due giorni con 24 processi; nel paper darebbe una frase nelle limitazioni.
2. **Quale versione del test a segni alterni è la principale? Consiglio quella senza filtro sui coefficienti del modello ridotto.**
   - Risponde in quasi tutti i dataset, e la quota di risultati usabili conta nella regola della Figura 1; non è mai fuori taratura verso l'alto ed è prudente solo dove i dati non hanno quasi informazione.
   - Alternativa: la versione con il filtro, scritta prima del calcolo. È più vicina a 5 su 100 dove risponde, ma con 10 soggetti e 2 visite risponde in 3 dataset su 10.
   - In ogni caso si calcolano entrambe: l'altra entra come analisi di sensibilità.

**Comunque:**

- nel run su tutte le situazioni scriveremo prima del calcolo le attese per i test scelti, così le correzioni introdotte qui dopo i risultati vengono confermate o smentite su dati nuovi per loro;
- prima di quel run verificheremo che geesmv con la correlazione AR(1) restituisca davvero la correzione: nei dataset dei tempi risponde in 2 centesimi di secondo e non va mai al limite di 5 secondi, a differenza delle altre correlazioni;
- nello stesso run terremo conto che il CR2 di clubSandwich fallisce sulla GEE con correlazione AR(1) in circa 6 dataset su 10;
- la versione di fewclust usata qui per i tempi è la 0.1.0 presente sul server; quella allineata al paper non è ancora in un commit registrato, e servirà prima di usarla come implementazione di confronto e di citarla;
- quando aggiorneremo il protocollo, la versione nuova conterrà i test scelti qui;
- se un test a segni alterni entra nel paper, andrà affiancato al modello marginalizzato multilivello, che fornisce stima e intervallo.

Appena decidi questi 2 punti, lo step è chiuso.

Il PI ha deciso i 2 punti il 2026-10-22: step chiuso.

**Sui tempi di calcolo.** Stima dell'agente da tempi misurati in serie su 417 dataset (10, 25 e 145 soggetti), con un processore, riportati al numero di osservazioni delle 1.896 situazioni con 1.000 simulazioni per scenario sia senza sia con effetto del trattamento. Il run successivo come è pianificato (GEE con tre correlazioni di lavoro, le tre correzioni di geesmv e di fewclust per ciascuna, il CR2 di clubSandwich, il test a segni alterni) costa circa 3.300 ore di processore, circa 6 giorni con 24 processi. Il 93% di quel tempo è nelle correzioni di geesmv, che con pochi soggetti arrivano spesso al limite di 5 secondi; le stesse correzioni di fewclust costano 34 volte meno. Il tempo per generare i dati non è misurato qui.


# Scelte del PI sul riassunto (2026-10-22)

- **Scelta 1, test da aggiungere al run: il solo test a inversione di segno.** Vedi [Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](../decisions/D-46.md).
- **Scelta 2, versione principale: senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita'.** Vedi [Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](../decisions/D-47.md).

# Dubbi rinviati

Rinviati al [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](S7d.md), che li riceve nel suo file:

1. geesmv con `corstr = "AR-M"` risponde in circa 0,02 s e non va mai in time-out (exchangeable e indipendente: fino a 0,57 s e time-out): verificare che restituisca la correzione e non un valore degenere, con il cancello di coincidenza geesmv / fewclust gia' previsto.
2. `clubSandwich::vcovCR(type = "CR2")` fallisce sulla GEE AR(1) nel 59% dei dataset dei tempi ("leading minor not positive"): decidere come contarlo (non usabile) prima del run.
3. Versione di fewclust: qui 0.1.0 dal commit `5103bb3`; la versione allineata al paper ([fewclust allineato alla regola del paper e citato nel manoscritto (D-43)](../decisions/D-43.md)) non e' in un commit registrato.
4. Organizzazione del run: geesmv e' il 93% del tempo stimato; valutare (decisione del PI, oltre 4 ore) se usare fewclust su tutta la griglia e geesmv solo sul sottoinsieme del cancello di coincidenza. Il tempo di generazione dei dati va misurato.
