---
type: Audit Report
title: "Audit di conformita' a Morris 2019, punto di partenza"
description: "Esito dell'audit iniziale dei 28 requisiti di Morris, White e Crowther (2019) sul manoscritto precedente, con il dettaglio per requisito."
status: stable
resource: "MORRIS2019_conformita.csv"
sources:
  - id: morris2019
    resource: https://doi.org/10.1002/sim.8086
    title: "Morris, White, Crowther (2019). Using simulation studies to evaluate statistical methods. Statistics in Medicine 38:2074-2102"
---
# Conformità a Morris, White & Crowther (2019), *Stat Med* 38:2074–2102

**Oggetto:** studio di simulazione MMM (codice in `~/Documents/Projects/MMM`, commit HEAD `cccef18`; run eseguito da `5f84297` con `R/run_full.R`) e manoscritto finale `ARTICOLO 2026/FINAL/MMM_manuscript_SMMR.md` (17 ago 2026). Il `.docx` del 27 ago coincide con il testo .md al 99,7% e in più contiene 7 commenti di L. Vedovelli.
**Metodo:** ogni requisito di Morris 2019 (ADEMP, §4 RNG ed esecuzione, §5 analisi e MCSE, §6 reporting, §8.1 protocollo e codice) è stato confrontato con il testo e con il codice effettivamente eseguito. La tabella completa con l'evidenza file per file è in `MORRIS2019_conformita.csv`.

## Esito

| esito | n |
|---|---|
| Conforme | 4 (B7, C1, D2, G1 con riserve) |
| Parziale | 10 |
| Non conforme | 12 |
| Non verificabile / assente | 2 |

**Lo studio** (il codice) rispetta i punti di Morris che contano di più per la validità del confronto: tutti i metodi sono applicati agli stessi dataset (§4.3), i fallimenti sono trattati come esito (§5.1) con un'analisi sull'insieme comune, ogni replicato è riproducibile e la verifica del codice è ampia. Ci sono però **tre discrepanze fra codice e testo che toccano il DGM**: [B3: Arm nullo e alternativo differiscono solo per il parametro testato](requirements/B3.md), [B4: Descrizione dei meccanismi coerente con il codice](requirements/B4.md) e B2.

**Il paper** non è conforme alle richieste di reporting di Morris: mancano i valori dei parametri dei DGM, la giustificazione di n_sim, il seme e la dipendenza fra dataset, le MCSE nelle tabelle e nelle figure, la disponibilità del codice e la struttura ADEMP.

## Rilievi da risolvere prima della submission (in ordine di gravità)

1. **B3: H0 e H1 non differiscono solo per il parametro testato (GenBinaryY, M3–M5).** In `R/generation_model3-5.R` sigma e/o gamma sono diversi fra i due arm. Per esempio in M5 gamma vale c(1, 0.1, 0.05) sotto H1 e c(1, 0, 0) sotto H0, cioè la transizione dipende da trt e steroids solo sotto H1. Questo ha due conseguenze. L'errore di I tipo di questi scenari è stimato sotto un DGM privo della struttura di dipendenza che il testo descrive. La potenza size-adjusted, inoltre, applica il valore critico del nullo a dati con una dipendenza diversa. I commenti del codice dicono il contrario ("tested parameter ONLY"). GLMM-native e AR1 sono corretti.
2. **B4: §2.2 descrive in modo errato GLMM-native.** I coefficienti condizionali non sono "calibrated per replicate": sono costanti fisse, pari a 1,16 volte il marginale (approssimazione di Zeger). La calibrazione per replicato riguarda AR1. Inoltre in M5 la varianza della pendenza su dose è 0.01 in AR1 e 0.001 in GLMM-native, quindi per M5 AR1 non è una generalizzazione stretta di GLMM-native.
3. **B2: T è un massimo, non il numero di occasioni.** Il 30% dei soggetti ha da 1 a T−1 occasioni, e a T=2 ha una sola osservazione (media 1,7 per soggetto). L'argomento di §3.8 ("two observations… saturate") va riformulato. Il troncamento va dichiarato in §2.3.
4. **[F3: Dipendenza fra dataset simulati dichiarata](requirements/F3.md): la dipendenza fra dataset non è dichiarata.** Il seme `20260805 + 1000·model + i` è lo stesso per H0/H1, per ogni N e T e per i tre meccanismi. In particolare GLMM-native e GLMM-AR1 condividono disegno ed effetti casuali (M1–M4). Per Morris §5.4 i numeri casuali comuni sono un vantaggio nei confronti fra metodi, ma vanno riportati (§6.1). Ci sono due conseguenze: (i) la "near-duplicate" di §2.2 e §4.4 e i ranking identici sono in parte costruiti dal rumore condiviso; (ii) l'inferenza che tratta come indipendenti i 1.440 punti o le 11.163 celle (BH, Spearman p = 4,9×10⁻²¹, Wilcoxon, bootstrap delle celle vinte) va motivata o rifatta con un bootstrap per blocchi di replicato.
5. **E1/H2/H3: misure di prestazione e MCSE.** §2.4 dichiara bias, SE empirico, SE del modello, MSE e copertura con MCSE. Nei Risultati il SE empirico e quello del modello non compaiono mai, il bias compare solo per il GLMM e la MCSE solo due volte. Alcune affermazioni sono sotto la MCSE: per livello di N la MCSE è di circa 0,125 pp (`tables/G_level_fixed.csv`), mentre §3.1 riporta deviazioni di "0.001" e "0.009" e §4.1 "below 0.1 points". Andrebbero riformulate come "entro la MCSE" e arrotondate di conseguenza.
6. **[F1: n_sim giustificato](requirements/F1.md): n_sim non è giustificato.** Con 1.000 replicati la MCSE per cella al 5% è 0,69 pp, e 0,98 pp quando metà dei fit è usabile. Va dichiarata la precisione attesa (Morris §5.3).
7. **E2/E3: definizione di "usabile".** `converged` vuol dire soltanto "nessun errore R". I warning di convergenza o singolarità di lme4 restano fra i fit usabili. I fallimenti di generazione vengono scartati senza registrarne il motivo, e da qui vengono le 132 celle con 971–999 replicati.
8. **D1/H5: software e codice.** Il supplemento non è aggiornato. `renv.lock` non contiene geesmv e indica R 4.5.2 contro l'R 4.6 del run. Le opzioni di stima non sono dichiarate (Laplace, ottimizzatori, q=20, timeout di 5 s su geesmv fino al 16,8% per cella). Manca un code availability statement, e `_targets.R` non è il driver effettivo.
9. **H6: la pre-specificazione della regola decisionale (§4.2) non è verificabile.** "Unconditional validity" compare nel git per la prima volta il 2026-08-10, dopo il run del 2026-08-05.

## Interventi minori

- A1/H1: aggiungere un paragrafo Aims e una tabella ADEMP, e riordinare §2 in A-D-E-M-P.
- B5/B6: in GenBinaryY la "random slope" è un carico di una sola variabile latente che varia nel tempo. La varianza della pendenza di 0.001 (DS 0,032) va motivata e dichiarata fra i limiti.
- C2: tabulare i valori veri per meccanismo, scenario e T.
- D3: dichiarare che MMM è specificato correttamente su GenBinaryY e il GLMM su GLMM-native.
- E4: dare l'incertezza Monte Carlo anche per la potenza size-adjusted e per il rango combinato.
- H4: zip plot o istogrammi delle stime nel supplemento.
- Figura 1, pannello a (non è un requisito di Morris): la regola non copre il ramo "2–3 occasioni", che è la prima raccomandazione dell'abstract, né l'intervallo 25–45 cluster.

## Limiti di questa verifica

- I chunk grezzi (`~/MMM/results_full` sul server) e `O_level_all.csv`, citato da §2.3 come fonte dei numeri, non sono disponibili in locale. Per questo non ho quantificato quanti fit GLMM usabili portano warning (E2), né quanti replicati falliscono in generazione e perché (E3).
- L'unica tabella di livello disponibile in locale (`tables/G_level_fixed.csv`) dà per Fay-Graubard su GLMM-native una deviazione di 0,306 pp a N=125. È in tensione con "fluctuating within 0.2 points" (§3.1). Non l'ho riconciliata perché il testo dichiara un'altra fonte.
- Non ho modificato nessun numero del manoscritto: i rilievi sono riportati, non corretti.
