# Accorciamento del testo principale di MMM v5: proposta per lo step successivo

## In breve: si può portare il testo principale sotto le 7.850 parole senza perdere il filo?

**A cosa serviva.** Il testo principale contava 9.908 parole (8.508 di testo più 200 per ciascuna delle 7 tabelle e figure). Il tetto è 7.850, con circa 150 parole già riservate alle aggiunte decise dalla revisione. Qui si propone cosa tenere e cosa spostare nel supplemento, senza cancellare nulla e senza toccare le frasi che reggono abstract, raccomandazione e conclusioni.

**Cosa abbiamo trovato.**

1. **Sì: il testo accorciato conta 7.824 parole, 26 sotto il tetto.** Il conteggio viene dal contatore del progetto, applicato a un render simulato in cui le modifiche sono state fatte davvero: 6.624 parole di testo più 200 per ciascuna delle 6 tabelle e figure rimaste.
2. **La Tabella 5, quella della potenza, va nel supplemento: è identica alla Tabella S10.** Ho confrontato i 48 valori sul registro dei numeri del server: coincidono tutti. Lo spostamento da solo vale 200 parole.
3. **Si spostano 83 frasi (2.103 parole) e se ne aggiungono 12 brevi di collegamento.** Le frasi spostate sono dettagli di implementazione, analisi di sensibilità, diagnosi e ripetizioni di tabelle. Vanno in 13 sezioni del supplemento che esistono già e in 4 sezioni nuove; 18 di loro ripetono cose che il supplemento dice già.
4. **Nessuna frase che regge l'abstract, la raccomandazione o le conclusioni è stata spostata.** Per ognuna delle 15 frasi di abstract e conclusioni c'è una frase dei Risultati (le prime due, di contesto, sono rette dall'Introduzione) che resta nel testo; l'elenco è più sotto.
5. **Restano nel testo alcune cose che si potevano spostare, perché servono al filo o a leggere onestamente la regola.** Sono: la giustificazione dei valori dei parametri, la scelta di 1.000 simulazioni per scenario, il limite della regola (con 1.000 simulazioni uno scenario risulta fuori soglia solo da un 9,3% osservato) e la definizione di "il test non risponde" per il test a inversione di segno, da cui dipende il 72% dell'abstract.

**Cosa devi decidere.** Nulla di nuovo: la proposta segue le indicazioni del 1 novembre. Quando le modifiche verranno applicate restano da fare solo le operazioni redazionali elencate in fondo.

**Comunque**: quando riscriveremo il manoscritto, le frasi spostate entrano nel supplemento parola per parola, con i rimandi corretti ("of the main text", "Table S10"), e la didascalia della Tabella S10 riceve le informazioni della didascalia della Tabella 5.

Non c'è nulla da decidere: se non hai dubbi, questa parte della revisione è chiusa.

## 1. Conteggio

Contatore del progetto (`s7_wordcount.py`, Python 3.13 locale, necessario per l'f-string con backslash), stessa regola: testo da "1. Introduction" a "Code and data availability", senza commenti, tabelle, didascalie e formule, più 200 per ogni tabella o figura.

| render | testo | tabelle e figure | totale |
|---|---|---|---|
| `MMM_v5_render.md` (partenza, verificato) | 8.508 | 7 × 200 = 1.400 | 9.908 |
| `MMM_v5_render_short.md` (simulato) | 6.624 | 6 × 200 = 1.200 | **7.824** |

Output del contatore sul render simulato:

```
P1: words 6624 + 200 x 6 display items = 7824
P2: sentences 302, median 20.0, mean 22.0, p90 35, max 66
abstract words 199
em dashes 0, en dashes 0
```

Il render simulato è il render originale con: i 31 paragrafi toccati ricostruiti dalle frasi del qmd (valori `r n()` sostituiti con quelli letti dal render stesso, 339 id allineati senza conflitti), i paragrafi non toccati copiati verbatim, la Tabella 5 tolta. La ricostruzione del qmd non modificato restituisce esattamente 8.508 parole, quindi il metodo non introduce scarti. Abstract, bibliografia, altre tabelle e figure sono invariati.

## 2. Parole per sezione

| sezione | prima | dopo | risparmio |
|---|---|---|---|
| 1. Introduction | 574 | 462 | 112 |
| 2. Methods | 53 | 53 | 0 |
| 2.1 Aims | 103 | 103 | 0 |
| 2.2 Data-generating mechanisms | 839 | 673 | 166 |
| 2.3 Estimands and true values | 141 | 141 | 0 |
| 2.4 Methods compared | 812 | 545 | 267 |
| 2.5 Performance measures | 712 | 531 | 181 |
| 2.6 Decision rule | 205 | 132 | 73 |
| 2.7 Clinical applications | 352 | 249 | 103 |
| 3.1 Type I error on the main grid | 341 | 264 | 77 |
| 3.2 The convergence check of the marginalized model | 355 | 174 | 181 |
| 3.3 Availability of a result | 148 | 148 | 0 |
| 3.4 Reference distribution and the decision rule | 576 | 531 | 45 |
| 3.5 Working correlation and implementation of the corrections | 144 | 50 | 94 |
| 3.6 The fourth mechanism | 238 | 192 | 46 |
| 3.7 Sign-flip score test | 147 | 85 | 62 |
| 3.8 Power and coverage | 117 | 39 | 78 |
| 3.9 Scenario anchored to the cohort | 139 | 71 | 68 |
| 3.10 Adaptive quadrature for the mixed model | 228 | 147 | 81 |
| 3.11 Clinical applications | 449 | 265 | 184 |
| 4.1 The recommendation | 584 | 584 | 0 |
| 4.2 The convergence check | 190 | 124 | 66 |
| 4.3 Relation to previous work | 240 | 240 | 0 |
| 4.4 The clinical applications and the fourth mechanism | 135 | 135 | 0 |
| 4.5 Limitations | 444 | 444 | 0 |
| 5. Conclusions | 242 | 242 | 0 |
| Tabella 5 (oggetto grafico) | 200 | 0 | 200 |
| **Totale** | **9908** | **7824** | **2084** |

Il dettaglio per paragrafo (frasi tenute verbatim dal qmd, frasi ponte, destinazione, frasi spostate con identificativo) è in `accorciamento_S8.csv`; le 83 frasi spostate, verbatim e con la sezione di destinazione, sono in `accorciamento_S8_moved.csv`. Gli identificativi delle frasi sono "riga del qmd.indice della frase" contando da 0 (per l'elenco dei meccanismi alle righe 62-65: "62.voce.indice", con la voce numerata 1-4 come nell'elenco).

## 3. Abstract e conclusioni: la frase che le regge dopo gli spostamenti

Ogni frase di sostegno indicata è nel testo accorciato (controllo automatico: nessuna di esse è fra le frasi spostate). Testo come nel render.

**Abstract, frase 1.** "Binary outcomes measured repeatedly on few subjects are analysed with generalized estimating equations (GEE), mixed models or marginalized multilevel models." *(sfondo: regge l'Introduzione)*
- 34.1 (1. Introduction): "A generalized estimating equation (GEE) targets a population-averaged contrast and treats the dependence as a nuisance."
- 34.2 (1. Introduction): "A generalized linear mixed model (GLMM) targets a cluster-specific contrast and models the dependence explicitly."
- 34.3 (1. Introduction): "A marginalized multilevel model (MMM) specifies both: a marginal mean model of scientific interest and a conditional model that generates the dependence [1,2]."

**Abstract, frase 2.** "Small-sample corrections to the sandwich variance have been evaluated mainly by number of clusters, rarely by number of occasions." *(sfondo: regge l'Introduzione)*
- 42.4 (1. Introduction): "These corrections have been evaluated on subjects followed over time, with as few as ten subjects and four occasions [11-15]."
- 42.5 (1. Introduction): "They have also been evaluated on parallel cluster-randomized trials, where the simulation evidence has been reviewed systematically [16,17]."
- 42.10 (1. Introduction): "Whether their recommendations hold with two or three occasions per subject has not been established."

**Abstract, frase 3.** "We compared analysis strategies over 2,376 simulated design points (10-145 clusters, 2-10 occasions, four mechanisms, one purely serial) and in subsamples of two clinical datasets."
- 133.0 (3.4 Reference distribution and the decision rule): "Table 4 applies the rule of Section 2.6 to the 2,376 design points of the four mechanisms."
- 147.0 (3.6 The fourth mechanism): "On Copula-AR1 (480 design points) neither the GEE strategy nor the marginalized model with the check was clearly above 7.5% at any design point (Table 4)."
- 169.3 (3.11 Clinical applications): "Over all subsamples of this size, including those, a usable result was obtained by the mixed model in 78%, by the marginalized model in 77%, by the uncorrected GEE in 70%, by Fay-Graubard in 34% and by Mancl-DeRouen in 15%."
- 171.1 (3.11 Clinical applications): "In subsamples of 10 patients, with the full mean model, a usable result was returned by the mixed model in 53% of draws, by the marginalized model in 51%, by the uncorrected GEE in 47%, by Fay-Graubard in 23% and by Mancl-DeRouen in 23%; as in the cohort, the corrections mostly failed with an estimate but no standard error (Table S30)."

**Abstract, frase 4.** "GEE with a first-order autoregressive working correlation, the Fay-Graubard correction and a t reference, and the marginalized multilevel model with a convergence check and a t reference, were never clearly above 7.5% type I error at any design point and returned a usable result in 92% and 91% of datasets."
- 133.1 (3.4 Reference distribution and the decision rule): "Among the strategies with an interval, only two were never clearly above 7.5%: GEE with AR(1) working correlation, Fay-Graubard and a t reference; and the marginalized model with the check and a t reference."
- 133.2 (3.4 Reference distribution and the decision rule): "The rule selected the GEE strategy, which returned a usable result in 92% of datasets, against 91% for the marginalized model, a difference within the margin for usable shares (Section 2.5)."

**Abstract, frase 5.** "Their calibration differed materially only at ten clusters, where GEE was closer to nominal."
- 139.0 (3.4 Reference distribution and the decision rule): "Compared directly, with the block bootstrap and the margins of Section 2.5 (Figure 1, Table S31), the two strategies had similar mean absolute departures from 5% over all design points (both 1.1 percentage points), similar usable shares (difference 1.2 points) and similar power (difference 1.4 points)."
- 139.1 (3.4 Reference distribution and the decision rule): "The GEE strategy was closer to 5% at ten clusters with three or more occasions (1.6 against 2.9 points), where the marginalized model was conservative (mean type I error 2.1% against 4.0%)."
- 139.2 (3.4 Reference distribution and the decision rule): "With two occasions and 15 or more clusters the GEE strategy was conservative (2.5%) and the marginalized model was closer to 5% by about the margin (difference 1.0 points); with 15 or more clusters and three or more occasions the two were similar."

**Abstract, frase 6.** "The mixed-model Wald test was anticonservative under the Laplace approximation; with adaptive quadrature its excess remained only with ten occasions."
- 135.1 (3.4 Reference distribution and the decision rule): "Mancl-DeRouen with a t reference was clearly above 7.5% at 16 to 57 design points across working correlations, Kauermann-Carroll at 132 to 199, and the Wald test of the mixed model at 515, 55 of them on data generated by its own model (GLMM-native)."
- 147.4 (3.6 The fourth mechanism): "The Wald test of the mixed model (Laplace approximation) with a t reference rejected a true null in 10% of usable replicates on average, 3.2 points more than on the three study mechanisms, and was clearly above 7.5% at 223 of the 480 design points."
- 165.2 (3.10 Adaptive quadrature for the mixed model): "With a t reference, quadrature lowered the rejection rate under the null by 3.1 to 4.0 percentage points on average, depending on the mechanism."
- 165.4 (3.10 Adaptive quadrature for the mixed model): "The Wald test was clearly above 7.5% at 33 of these design points under the Laplace approximation and at 8 with quadrature, all with T = 10: 6 on Copula-AR1 (12% to 16%) and 2 on GenBinaryY (9.3% and 9.6%)."

**Abstract, frase 7.** "With two or three occasions and up to 25 clusters both strategies answered in fewer than half of datasets; a sign-flip score test answered in 72% and kept its level, but gives no estimate." *("gives no estimate": Metodi 2.1 (52.2) e frase ponte di 2.4)*
- 137.0 (3.4 Reference distribution and the decision rule): "With two or three occasions and at most 25 clusters, the GEE strategy and the marginalized model returned a usable result in fewer than half of the datasets (46% and 39%)."
- 151.0 (3.7 Sign-flip score test): "Over the 2,376 design points the sign-flip test rejected a true null in 5.4% of datasets on average, and it was clearly above 7.5% or clearly below 2.5% at no design point."
- 151.2 (3.7 Sign-flip score test): "It answered in 97% of datasets with an estimable design and, with two or three occasions and at most 25 clusters, in 72%, against 46% for the GEE strategy (Table 4)."
- 151.4 (3.7 Sign-flip score test): "Its type I error in those small designs was 5.2%, and its power per attempted dataset over all design points was 21%."

**Abstract, frase 8.** "We recommend the GEE strategy, which needs no convergence check; the marginalized model is a comparable alternative." *(raccomandazione)*
- 133.2 (3.4 Reference distribution and the decision rule): "The rule selected the GEE strategy, which returned a usable result in 92% of datasets, against 91% for the marginalized model, a difference within the margin for usable shares (Section 2.5)."
- 139.0 (3.4 Reference distribution and the decision rule): "Compared directly, with the block bootstrap and the margins of Section 2.5 (Figure 1, Table S31), the two strategies had similar mean absolute departures from 5% over all design points (both 1.1 percentage points), similar usable shares (difference 1.2 points) and similar power (difference 1.4 points)."
- 119.1 (3.2 The convergence check of the marginalized model): "Without the check, the marginalized model was clearly above 7.5% at 149 of the 864 null design points of these scenarios, 90 of them with eight or ten occasions."
- 119.2 (3.2 The convergence check of the marginalized model): "With the check it was clearly above 7.5% at 0."

**Conclusioni, frase 1.** "With few clusters and short series, the choice of analysis for a clustered binary outcome decides whether the test keeps its nominal level and whether it returns an answer at all."
- 135.0 (3.4 Reference distribution and the decision rule): "With exchangeable or independence working correlation, Fay-Graubard with a t reference was clearly above 7.5% at 11 and 2 design points, none of them in the main grid or M1c."
- 135.1 (3.4 Reference distribution and the decision rule): "Mancl-DeRouen with a t reference was clearly above 7.5% at 16 to 57 design points across working correlations, Kauermann-Carroll at 132 to 199, and the Wald test of the mixed model at 515, 55 of them on data generated by its own model (GLMM-native)."
- 135.2 (3.4 Reference distribution and the decision rule): "The likelihood-ratio test was clearly above 7.5% at 412 design points."
- 125.1 (3.3 Availability of a result): "With random slopes and two occasions, the mixed model returned a usable Wald test in 0.002% of attempts."
- 125.4 (3.3 Availability of a result): "From four occasions every method returned a usable result in more than 90% of attempts."
- 137.0 (3.4 Reference distribution and the decision rule): "With two or three occasions and at most 25 clusters, the GEE strategy and the marginalized model returned a usable result in fewer than half of the datasets (46% and 39%)."

**Conclusioni, frase 2.** "GEE with an AR(1) working correlation, the Fay-Graubard correction and a t reference, and the marginalized multilevel model fitted with a convergence check and a t reference, were never clearly above 7.5% at any simulated design point of four mechanisms, and their calibration differed beyond the relevance margin only at ten clusters."
- 133.1 (3.4 Reference distribution and the decision rule): "Among the strategies with an interval, only two were never clearly above 7.5%: GEE with AR(1) working correlation, Fay-Graubard and a t reference; and the marginalized model with the check and a t reference."
- 139.1 (3.4 Reference distribution and the decision rule): "The GEE strategy was closer to 5% at ten clusters with three or more occasions (1.6 against 2.9 points), where the marginalized model was conservative (mean type I error 2.1% against 4.0%)."
- 139.2 (3.4 Reference distribution and the decision rule): "With two occasions and 15 or more clusters the GEE strategy was conservative (2.5%) and the marginalized model was closer to 5% by about the margin (difference 1.0 points); with 15 or more clusters and three or more occasions the two were similar."

**Conclusioni, frase 3.** "The Wald test of the mixed model fitted by the Laplace approximation was anticonservative; with adaptive quadrature, examined for a random intercept, the excess remained only with ten occasions, mainly under serial dependence, and with few clusters the test became conservative."
- 135.1 (3.4 Reference distribution and the decision rule): "Mancl-DeRouen with a t reference was clearly above 7.5% at 16 to 57 design points across working correlations, Kauermann-Carroll at 132 to 199, and the Wald test of the mixed model at 515, 55 of them on data generated by its own model (GLMM-native)."
- 147.4 (3.6 The fourth mechanism): "The Wald test of the mixed model (Laplace approximation) with a t reference rejected a true null in 10% of usable replicates on average, 3.2 points more than on the three study mechanisms, and was clearly above 7.5% at 223 of the 480 design points."
- 165.2 (3.10 Adaptive quadrature for the mixed model): "With a t reference, quadrature lowered the rejection rate under the null by 3.1 to 4.0 percentage points on average, depending on the mechanism."
- 165.4 (3.10 Adaptive quadrature for the mixed model): "The Wald test was clearly above 7.5% at 33 of these design points under the Laplace approximation and at 8 with quadrature, all with T = 10: 6 on Copula-AR1 (12% to 16%) and 2 on GenBinaryY (9.3% and 9.6%)."
- 165.6 (3.10 Adaptive quadrature for the mixed model): "Quadrature also made the test conservative, mostly with few clusters: it was clearly below 2.5% at 37 design points, against 3 under the Laplace approximation (Table S34)."

**Conclusioni, frase 4.** "With two or three occasions and up to 25 clusters, neither strategy returned a usable result in half of the datasets."
- 137.0 (3.4 Reference distribution and the decision rule): "With two or three occasions and at most 25 clusters, the GEE strategy and the marginalized model returned a usable result in fewer than half of the datasets (46% and 39%)."

**Conclusioni, frase 5.** "A sign-flip score test can still test for an effect there, with little power, and no analysis estimates the effect reliably; in that regime the design limits what any analysis can show."
- 151.2 (3.7 Sign-flip score test): "It answered in 97% of datasets with an estimable design and, with two or three occasions and at most 25 clusters, in 72%, against 46% for the GEE strategy (Table 4)."
- 151.4 (3.7 Sign-flip score test): "Its type I error in those small designs was 5.2%, and its power per attempted dataset over all design points was 21%."
- 137.0 (3.4 Reference distribution and the decision rule): "With two or three occasions and at most 25 clusters, the GEE strategy and the marginalized model returned a usable result in fewer than half of the datasets (46% and 39%)."

**Conclusioni, frase 6.** "In the designs studied we recommend GEE with an AR(1) working correlation, the Fay-Graubard correction and a t reference with N - p degrees of freedom, which needs no convergence check and was closer to nominal at ten clusters." *(raccomandazione)*
- 133.2 (3.4 Reference distribution and the decision rule): "The rule selected the GEE strategy, which returned a usable result in 92% of datasets, against 91% for the marginalized model, a difference within the margin for usable shares (Section 2.5)."
- 131.2 (3.4 Reference distribution and the decision rule): "Under the t reference the marginalized model became conservative at ten clusters."
- 131.3 (3.4 Reference distribution and the decision rule): "There it was clearly below 2.5% at 16 of 102 design points, against 4 for Fay-Graubard."
- 139.1 (3.4 Reference distribution and the decision rule): "The GEE strategy was closer to 5% at ten clusters with three or more occasions (1.6 against 2.9 points), where the marginalized model was conservative (mean type I error 2.1% against 4.0%)."

**Conclusioni, frase 7.** "The marginalized model with the convergence check and a t reference is a competitive alternative, slightly closer to nominal with two occasions and 15 or more clusters." *(raccomandazione)*
- 139.2 (3.4 Reference distribution and the decision rule): "With two occasions and 15 or more clusters the GEE strategy was conservative (2.5%) and the marginalized model was closer to 5% by about the margin (difference 1.0 points); with 15 or more clusters and three or more occasions the two were similar."

Anche le frasi dei Risultati che reggono la raccomandazione in 4.1 restano: 133.1, 133.2, 139.0-139.3 (strategie e confronto diretto), 143.0 e 137.1 (correlazione di lavoro e stessa scelta sui tre meccanismi di studio), 127.2 e 135.2 (criterio precedente e test del rapporto di verosimiglianza), 135.0 e 147.7 (Fay-Graubard con correlazione exchangeable), 137.0 e 151.2 (piccoli disegni), 125.1 (pendenza casuale con due occasioni). La sezione 4.1 non è toccata.

## 4. Destinazioni nel supplemento

| sezione | frasi | parole | frasi spostate |
|---|---|---|---|
| S1. Parameters of the data-generating mechanisms | 4 | 74 | 62.3.2, 67.5, 69.4, 71.1 |
| S2. Software | 1 | 20 | 71.5 |
| S3. Tables of the main grid | 16 | 317 | 81.8, 89.1, 91.11, 91.12, 93.9, 113.0, 113.2, 115.0, 115.1, 131.4, 143.1, 155.0, 155.1, 155.2, 155.3, 155.4 |
| S4. Numerical behaviour of the small-sample corrections | 5 | 131 | 79.4, 79.5, 143.2, 143.3, 143.4 |
| S5. Tables of the scenario anchored to the cohort (M1c) | 4 | 68 | 161.4, 161.5, 161.6, 161.7 |
| S6. Degrees of freedom of the t reference | 1 | 18 | 83.1 |
| S7. Sign-flip score test and the interval obtained by inverting it | 6 | 166 | 85.2, 85.3, 85.5, 85.6, 151.1, 151.3 |
| S8. The fourth mechanism | 4 | 118 | 62.4.2, 62.4.5, 67.6, 147.5 |
| S9. Second clinical application: the respiratory trial | 4 | 127 | 103.1, 103.2, 171.2, 171.3 |
| S10. The two strategies of the rule compared | 2 | 59 | 97.7, 139.4 |
| S11. What the decision rule guarantees | 3 | 85 | 89.6, 91.7, 97.6 |
| S12. Adaptive quadrature for the mixed model | 2 | 98 | 165.0, 165.1 |
| S13. Estimation by the marginalized model with the check in scenarios M6-M7 | 4 | 138 | 121.3, 121.5, 121.7, 121.4 |
| S14 (new). Previous evaluations of the small-sample corrections | 6 | 154 | 42.1, 42.2, 42.6, 42.7, 42.8, 44.7 |
| S15 (new). Further details of the methods and performance measures | 5 | 121 | 79.0, 83.3, 89.2, 93.10, 93.11 |
| S16 (new). The convergence check of the marginalized model: diagnosis and sensitivity analyses | 9 | 193 | 81.6, 81.7, 119.4, 119.5, 119.6, 185.1, 185.2, 185.3, 185.4 |
| S17 (new). The neonatal cohort: subsampling details | 7 | 216 | 101.9, 101.10, 169.2, 169.4, 169.5, 169.7, 169.8 |

Delle 83 frasi spostate, 18 ripetono cose che il supplemento dice già: la riga del supplemento che lo dimostra è nella colonna `repeats_supplement` del file delle frasi spostate e nella colonna `note` del CSV principale. Le altre portano nel supplemento un contenuto nuovo per il supplemento, sempre con i loro numeri del registro.

**Sezioni nuove del supplemento.**

- **S14. Previous evaluations of the small-sample corrections.** Dall'Introduzione: i numeri delle revisioni sui trial a cluster (42.1-42.2), le unità dentro i cluster e le taglie nei lavori precedenti (42.6-42.8), i tre risultati citati e non rivendicati (44.7).
- **S15. Further details of the methods and performance measures.** Da 2.4 e 2.5: la frase introduttiva con la Tabella 5 (79.0), il riferimento chi quadrato e gli ottimizzatori del test del rapporto di verosimiglianza e del modello misto con bobyqa (83.3), l'utilizzabilità del test del rapporto di verosimiglianza (89.2), le convenzioni di arrotondamento (93.10-93.11).
- **S16. The convergence check of the marginalized model: diagnosis and sensitivity analyses.** Da 2.4, 3.2 e 4.2: quando è stato fissato il controllo e su quante stime (81.6), le soglie alternative 0,0001 e 0,01 (81.7, 119.5-119.6), l'effetto sulla potenza aggiustata (119.4), il meccanismo dei due modi di fallire di `nlm` (185.1-185.4).
- **S17. The neonatal cohort: subsampling details.** Da 2.7 e 3.11: cosa preserva il sottocampionamento (101.9), la probabilità attesa di nessun ventilato invasivamente (101.10, 169.2), le percentuali sui fit restituiti (169.4-169.5), l'intercetta casuale per paziente-tempo e la pendenza casuale per le zone posteriori (169.7-169.8).

**Sezioni esistenti che ricevono testo.** S1 (valori dei parametri, M6-M7, giustificazione della frase su GLMM-AR1, dataset non generati), S2 (numero di dataset rigenerati), S3 (Tabella 5 dentro la Tabella S10, potenza per ottimizzatore, correzione di Benjamini-Hochberg, correlazioni di Spearman, avvisi dei fit, spostamento con il riferimento t), S4 (geesmv e fewclust), S5 (potenza e celle fuori taratura di M1c), S6 (N - 2), S7 (calcolo del test a inversione di segno e filtro sui coefficienti), S8 (descrizione del quarto meccanismo e controllo su 20.000 soggetti), S9 (modello e risultati del trial respiratorio, interpretazione dello scarto degli errori standard), S10 (confronto non pre-specificato), S11 (criterio precedente, intervallo esatto, limiti di utilizzabilità), S12 (disegno del rifit con quadratura), S13 (stima del modello marginalizzato in M6-M7, Fay-Graubard in M7).

## 5. Frasi ponte nuove (12, al massimo 24 parole nel render)

- 42 (1. Introduction), al posto di 42.1, 19 parole nel render: `Cluster-randomized trials with ten to thirty clusters are common [6], and many stepped-wedge trials have fewer than ten [7].`
- 42 (1. Introduction), al posto di 42.6, 23 parole nel render: `In the original evaluations and in the cluster-randomized comparisons, the unit within a cluster is never a subject measured on several occasions [8,10,16,18,19].`
- 79 (2.4 Methods compared), al posto di 79.0, 18 parole nel render: `The eight methods of Table 2 and the further strategies of Table 1 were applied to every dataset.`
- 85 (2.4 Methods compared), al posto di 85.6, 20 parole nel render: `The test gives a p-value only; its computation and the interval obtained by inverting it are described in Section S7.`
- 103 (2.7 Clinical applications), al posto di 103.1, 17 parole nel render: `The contrast of interest was the treatment effect, adjusted for baseline status, visit, age, sex and centre.`
- 121 (3.2 The convergence check of the marginalized model), in coda al paragrafo, 19 parole nel render: `Bias, standard errors and coverage of the marginalized model with the check in these scenarios are in Section S13.`
- 143 (3.5 Working correlation and implementation of the corrections), in coda al paragrafo, 11 parole nel render: `The two implementations of the corrections are compared in Section S4.`
- 155 (3.8 Power and coverage), al posto di 155.0, 12 parole nel render: `Raw and size-adjusted power of the eight methods are in Table S10.`
- 165 (3.10 Adaptive quadrature for the mixed model), al posto di 165.0, 17 parole nel render: `The refit covered `r n("f:gq|npts")` null design points of scenarios M1, M2 and M1c (Section S12, Table S34).`
- 169 (3.11 Clinical applications), al posto di 169.2, 24 parole nel render: `In subsamples of `r n("x:lus_N10")` patients, `r n("x:lus_inv_none")`% of draws contained no invasively ventilated patient, so that the mean model could not be fitted (Table S15).`
- 169 (3.11 Clinical applications), al posto di 169.4, 15 parole nel render: `Most failures of the corrections returned a point estimate without a standard error (Table S15).`
- 185 (4.2 The convergence check), al posto di 185.1, 24 parole nel render: `Because the slope is exponentiated, a slope near zero lies on the boundary, where `nlm` fails in the two ways described in Section 2.4.`

## 6. Commenti claim toccati

In dieci gruppi il commento `<!-- claim -->` o `<!-- noclaim -->` era attaccato a una frase spostata ma copriva anche frasi tenute. Dove si applica anche alla frase tenuta è stato riattaccato verbatim; dove conteneva id che riguardano solo la frase spostata, nel testo restano solo gli id della frase tenuta e gli altri seguono la frase nel supplemento: 113.1 tiene `cl:C12a|GB` (i `cl:C12c` vanno in S3); 121.6 tiene `reg:ns_mmmt_hi` e `reg:ns_nocc_t_hi` (i due `reg:ns_fgt_*` vanno in S13). Per non spezzare `cl:C17b` la frase 131.1 resta nel testo. I casi sono elencati nella colonna `note`.

## 7. Da fare quando le modifiche verranno applicate

1. Tabella S10: aggiungere alla didascalia la definizione di size adjustment, "Kauermann-Carroll computed with `fewclust`" e "Monte Carlo standard errors at most `r n(\"x:t5_mcse\")` points" (0,27) dalla didascalia della Tabella 5.
2. Nelle frasi spostate, "Table 5" diventa "Table S10" (79.4, 143.4) e i rimandi alle sezioni del testo principale prendono "of the main text" (per esempio 119.5, "Section 2.5").
3. In S7, 85.5 (filtro sui coefficienti) e 151.1 (il suo risultato) vanno uno dopo l'altro; 85.2-85.3 precedono il paragrafo sull'intervallo.
4. Titolo della 3.8: con la sola copertura e il rimando alla Tabella S10 può diventare "Coverage and power" (i titoli non entrano nel conteggio).
5. Rieseguire il controllo dei claim sul qmd accorciato: i commenti del supplemento devono portare gli stessi id delle frasi spostate.

## 8. Previsione pre-registrata

| previsione | atteso | osservato | esito |
|---|---|---|---|
| Conteggio finale, con 200 per tabella o figura, senza spostare frasi che reggono abstract o conclusioni | al massimo 8.000 | 7.824 (6.624 + 6 × 200); 0 frasi di sostegno spostate | confermata |

Il conteggio è deterministico (nessuna incertezza Monte Carlo): l'esito non dipende da intervalli.

## 9. File

- `accorciamento_S8.csv`: una riga per paragrafo toccato (31) più una riga per la Tabella 5.
- `accorciamento_S8_moved.csv`: le 83 frasi spostate, verbatim dal qmd, con destinazione e prova di ripetizione.
- `MMM_v5_render_short.md`: il render simulato su cui è stato fatto il conteggio.