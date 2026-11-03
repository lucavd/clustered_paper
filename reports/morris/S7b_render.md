# S7b: rianalisi senza nuove simulazioni

2026-10-21

## In breve: le tre limitazioni che si potevano affrontare senza nuove simulazioni cambiano qualcosa nel paper?

**A cosa serviva.** Il paper riscritto dichiara tre limitazioni che si
potevano verificare sui risultati già salvati, senza generare nuovi
dati: i gradi di libertà della distribuzione t (noi usiamo il numero di
soggetti meno i coefficienti del modello, Li e Redden usano il numero di
soggetti meno 2); il fatto che, nei dati generati dal modello misto, il
valore nullo dell’effetto sulla scala della popolazione non è
esattamente zero; e il fatto che l’applicazione clinica usa una sola
coorte. Abbiamo ricalcolato i test con l’altra scelta di gradi di
libertà su tutte le 1.896 situazioni simulate, quantificato cella per
cella il rifiuto in più dovuto a quel valore nullo, e ripetuto
l’applicazione clinica su un dataset pubblico (lo studio respiratorio di
geepack: 111 pazienti, 4 visite, trattamento attivo contro placebo).

**Cosa abbiamo trovato.**

1.  **Con i gradi di libertà di Li e Redden la raccomandazione non
    cambia: il modello marginalizzato multilivello con il controllo di
    convergenza e la distribuzione t resta l’unica strategia mai
    chiaramente sopra 7,5 falsi positivi su 100, e la regola della
    Figura 1 lo sceglie di nuovo.**
    - Passare al numero di soggetti meno 2 alza i falsi positivi di meno
      di un quarto di punto in media sulle 1.896 situazioni, per ogni
      test di Wald; a 10 soggetti li alza di 0,6 punti per il modello
      marginalizzato e di 0,9 per la correzione di Fay e Graubard.
    - A 10 soggetti il modello marginalizzato resta conservativo, ma
      meno: chiaramente sotto 2,5 falsi positivi su 100 in 28 situazioni
      su 126, invece di 53.
    - La correzione di Fay e Graubard con la t diventa chiaramente sopra
      7,5 in 7 situazioni invece di 2: tutte con dati generati sulla
      scala della popolazione, variabilità latente che cresce nel tempo
      e al più 25 soggetti. La frase del paper che la descrive resta
      vera, con i numeri aggiornati.
2.  **Il valore nullo non esattamente zero del modello misto non può
    spostare nessuna cella: il rifiuto in più che produce è al massimo
    0,03 falsi positivi su 100 negli scenari principali e 0,26 negli
    scenari con variabilità crescente (145 soggetti, 10 visite).**
    - Sottraendolo dai tassi osservati, nessuna delle bandiere
      “chiaramente sopra 7,5” o “chiaramente sotto 2,5” cambia, per
      nessun metodo che stima l’effetto sulla scala della popolazione.
    - Un solo valore grande (1,5 punti) viene da una cella in cui lo
      stimatore con shrinkage ha 2 soli risultati usabili, entrambi
      collassati a zero: la cella non entra nelle valutazioni perché
      sotto i 20 risultati usabili.
    - La limitazione può diventare una frase quantitativa.
3.  **Nello studio respiratorio l’effetto del trattamento è ben
    identificato (log odds 1,26, errore standard 0,35, per il modello
    marginalizzato e la GEE; 2,06 sulla scala condizionale per il
    modello misto) e nei sottocampioni si ripetono i risultati della
    coorte neonatale, con una differenza.**
    - Con 10 pazienti e il modello con sei covariate, un risultato
      usabile arriva dal modello marginalizzato in 51 estrazioni su 100,
      dalla GEE in 47, dal modello misto in 53, da Fay e Graubard e da
      Mancl e DeRouen in 23: le correzioni falliscono per un errore
      numerico che non restituisce l’errore standard, come nella coorte
      neonatale. Con il modello ridotto (trattamento e stato basale)
      tutti i metodi rispondono in circa 8 estrazioni su 10.
    - La distorsione del modello marginalizzato rispetto alla stima
      sull’intero dataset scende da 0,15 con 10 pazienti a zero con 50,
      come nella coorte neonatale (da 0,16 a 0,02). La stima
      condizionale del modello misto resta invece sotto quella
      dell’intero dataset di 0,59 con 10 pazienti e di 0,15 con 50.
    - La differenza: il rapporto fra errore standard del modello e
      dispersione delle stime nei sottocampioni, per il modello
      marginalizzato con sei covariate, va da 0,57 con 10 pazienti a
      1,18 con 50 (nella coorte neonatale da 0,87 a 1,14); con 10-25
      pazienti e sette coefficienti le stime sono instabili e la loro
      dispersione è gonfiata da stime grandi che il filtro di usabilità
      tiene.
4.  **Un dettaglio di calcolo che tocca entrambe le applicazioni: la
    dispersione delle stime nei sottocampioni è più piccola del vero
    perché i sottocampioni sono estratti senza reinserimento da una
    coorte finita.**
    - Con 50 pazienti su 111 la dispersione va divisa per 0,74; con 50
      su 207 per 0,87. Il rapporto errore standard su dispersione, che
      il paper riporta per la coorte neonatale, è quindi gonfiato
      soprattutto a 50 pazienti: per il modello marginalizzato passa da
      1,14 a 0,99 nella coorte neonatale e da 1,18 a 0,87 nello studio
      respiratorio; per il modello misto nella coorte neonatale da
      0,55-0,65 a 0,51-0,58.
    - Le conclusioni del paper non cambiano (il modello misto sottostima
      l’errore standard in ogni caso), ma i numeri della frase sì.
5.  **Un errore di strumento evitato: il pacchetto geesmv usa come
    identificativo dei gruppi qualunque colonna si chiami `id`, se
    esiste, ignorando quella indicata.**
    - Nello studio respiratorio esiste una colonna `id` che numera i
      pazienti dentro ciascun centro: le correzioni risultavano quasi
      doppie. Rinominata la colonna, le tre correzioni coincidono entro
      la quarta cifra con un calcolo indipendente (il pacchetto
      clubSandwich e una formula scritta a mano).
    - Simulazioni e coorte neonatale non sono toccate: lì la colonna dei
      gruppi si chiama già `id`.

**Cosa devi decidere.**

1.  **Nelle due applicazioni riportiamo il rapporto errore standard su
    dispersione corretto per la popolazione finita? Consiglio di sì.**
    - L’errore standard del modello stima la variabilità di campioni
      indipendenti; la dispersione fra sottocampioni senza reinserimento
      la sottostima di un fattore noto (radice di 1 meno la frazione
      estratta). Correggere costa una frase nei metodi e cambia i numeri
      della frase sul rapporto nella sezione dei risultati.
    - Alternativa: tenere il rapporto grezzo e aggiungere una frase di
      cautela sul fattore. Non cambia nulla nel calcolo, ma il lettore
      vede un rapporto sopra 1 a 50 pazienti che non è un difetto dei
      metodi.
    - Alternativa: togliere il rapporto a 50 pazienti e tenerlo solo
      fino a 25, dove il fattore vale al più 0,94. Si perde il punto in
      cui il rapporto si avvicina a 1.
2.  **Quale scelta di gradi di libertà va nel testo principale?
    Consiglio di tenere la nostra (soggetti meno coefficienti) e di
    aggiungere quella di Li e Redden come analisi di sensibilità, con
    una frase e una tabella nel supplemento.**
    - La raccomandazione è la stessa con entrambe; con quella di Li e
      Redden il modello marginalizzato è meno conservativo a 10 soggetti
      e Fay e Graubard è chiaramente sopra la soglia in 7 situazioni
      invece di 2.
    - Alternativa: passare a Li e Redden come scelta principale. Allinea
      il paper al lavoro che citiamo, ma tutti i numeri della Tabella 4
      e della Figura 1 vanno rigenerati con quella scelta.

**Comunque:**

- quando aggiorneremo il manoscritto, le tre limitazioni diventeranno
  frasi con i numeri di questo step: la frase sui gradi di libertà, la
  frase quantitativa sul valore nullo del modello misto, e un paragrafo
  con una tabella nel supplemento per lo studio respiratorio (con la
  differenza dichiarata: lì il modello marginalizzato include la
  componente latente, come nella simulazione, perché le visite sono una
  serie temporale);
- nello stesso aggiornamento entrerà anche la frase che inquadra i due
  meccanismi generatori come modello di casa, ciascuno, di un metodo;
- la nota sul pacchetto geesmv e la colonna `id` andrà nella
  documentazione del codice e nel supplemento, accanto alla verifica
  numerica delle correzioni.

Appena decidi questi 2 punti, lo step è chiuso.

## Dubbi rinviati

Allo [Aggiornamento del manoscritto con le analisi aggiunte
(S7f)](../../docs/steps/S7f.md) passano, come testo da scrivere e non
come dubbi aperti: la frase sui gradi di liberta’ (N - p nel testo, N -
2 come sensibilita’, salvo diversa scelta del PI, vedi riassunto,
decisione 2); la frase quantitativa sul valore nullo marginale (massimo
0,03 punti negli scenari principali, 0,26 in M6-M7,
`tables/V4_S7b_null_summary.csv`); il paragrafo e la tabella
supplementare dello studio respiratorio (`tables/V4_S7b_resp_*.csv`); il
rapporto SE su dispersione con il fattore di popolazione finita in
entrambe le applicazioni se il PI lo approva
(`tables/V4_S7b_seratio_fpc.csv`); l’inquadramento dei meccanismi come
modello di casa di un metodo; la nota su geesmv e la colonna `id`.

## Cancelli e previsioni

Le attese sono nel file dello step ([Rianalisi senza nuove simulazioni
(S7b)](../../docs/steps/S7b.md)), committate prima dei calcoli
(`ed6fd75`). Le bandiere “chiaramente sopra / sotto” contengono gia’ la
MCSE: i conteggi di punti di disegno si giudicano come interi. Per PC6 e
PC7 la MCSE propagata e’ quella del SE empirico della cella (excess x
1/sqrt(2(n - 1))).

| cancello | esito | dettaglio |
|----|----|----|
| GC1 | superato | results_full/: tag=before files=28800 list_diff_vs_S7=1 sum_diff_vs_S6b=0 \| tag=after files=28800 list_diff_vs_S7=1 sum_diff_vs_S6b=0 (la riga di differenza con la lista dello step del manoscritto e’ il suo timestamp di testa; le somme coincidono con la lista verificata in S6b). Input dello step: md5 identici prima e dopo (results_v4/S7b/md5_inputs_before.txt, md5_inputs_after.txt) |
| GC2 | superato | 28440 design point x strategy rows compared with the stored tables (22680 study, 5760 new scenarios): n, nd, rate, hi, lo, pu, upow identical: TRUE; normal p recomputed from estimate and SE on 23682062 usable Wald rows, 0 beyond 1e-9 (max deviation 4.9e-15) |
| GC3 | superato | 7 Wald strategies x 1896 design points with the t at N - 2; usable rows with a non-finite statistic: 0 |
| GC4 | superato | 12 cells (4 with the largest excess, 8 at random, seed 20261021), 1e6 normal draws each: max deviation from the formula 0.045 points (normal), 0.049 (t) |
| GC5 | superato | respiratory: 444 rows, 111 patients, 4 visits each: TRUE; md5 of the CSV export 281e0d3d43ddca923b9debb6e2e996f1; the 2,500 cohort draws (seed 20260807) reproduced: TRUE |

| codice | quantita | osservato | intervallo | esito | nota |
|----|----|----|----|----|----|
| PC1 | design points with MMM (check) + t at N - 2 clearly above 7.5% | 0 +/- 0 | 0 a 0 | confermata |  |
| PC2 | strategy chosen by the rule with N - 2 (1 = mmm\|t2, 0 = other) | 1 +/- 0 | 1 a 1 | confermata | chosen: mmm\|t2 |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, gee (points) | 0.212 +/- 0.0077 | 0 a 1 | confermata |  |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, gee_fg (points) | 0.107 +/- 0.006 | 0 a 1 | confermata |  |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, gee_kc (points) | 0.126 +/- 0.006 | 0 a 1 | confermata |  |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, gee_md (points) | 0.127 +/- 0.0064 | 0 a 1 | confermata |  |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, glmer (points) | 0.0489 +/- 0.0038 | 0 a 1 | confermata |  |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, mmm (points) | 0.0955 +/- 0.0065 | 0 a 1 | confermata |  |
| PC3 | mean shift of the type I error rate, N - 2 minus N - p, all design points, shrink (points) | 0.0629 +/- 0.0047 | 0 a 1 | confermata |  |
| PC4 | mean shift at N = 10, mmm (points) | 0.633 +/- 0.08 | 0 a 2.5 | confermata |  |
| PC4 | mean shift at N = 10, gee_fg (points) | 0.851 +/- 0.075 | 0 a 2.5 | confermata |  |
| PC5 | design points with Fay-Graubard + t at N - 2 clearly above 7.5% | 7 +/- 0 | 2 a 10 | confermata |  |
| PC6 | largest expected excess rejection over all GLMM-native design points and marginal methods, normal reference (points) | 1.55 +/- 2.2 | 0 a 1.5 | indeterminata | M7 N=10 T=2 shrink; empirical SE 0.0077, rate 0.0 |
| PC6b | the same, shrinkage excluded (not pre-registered, explanation of PC6) | 0.26 +/- 0.012 | 0 a 1.5 | confermata | M7 N=145 T=10 mmm; empirical SE 0.0525, rate 6.8 |
| PC7 | largest expected excess rejection, main scenarios M1-M5 and M1c, normal reference (points) | 0.0292 +/- 0.0013 | 0 a 0.5 | confermata | M4 N=145 T=10 mmm; empirical SE 0.0457 |
| PC8 | hi flags on GLMM-native that change after subtracting the excess, mmm + t | 0 +/- 0 | 0 a 0 | confermata |  |
| PC8 | hi flags on GLMM-native that change after subtracting the excess, gee_fg + t | 0 +/- 0 | 0 a 0 | confermata |  |
| PC9 | share of draws of 10 patients with a binary covariate at one level (%) | 8.8 +/- 2.5 | 0 a 25 | confermata |  |
| PC10 | usable on returned, 10 patients, full specification, MMM (%) | 74.1 +/- 4.7 | 75 a 100 | indeterminata | returned 347, usable 257 of 500 |
| PC11 | mean model SE / SD of estimates, MMM, full specification, 50 patients | 1.18 +/- 0.074 | 0.91 a 1.14 | indeterminata |  |
| PC12 | bias / SD of estimates, MMM, full specification, 50 patients | -0.0111 +/- 0.089 | -0.46 a 0.46 | confermata |  |

Spiegazioni delle previsioni non confermate:

- **PC6 (indeterminata).** Il massimo del rifiuto atteso in piu’ (1,55
  punti) viene dalla cella M7, 10 soggetti, 2 occasioni, stimatore con
  shrinkage: 2 soli replicati usabili su 971, stime collassate a zero
  (SE empirico 0,0077 contro 1,2-1,6 degli altri metodi nella stessa
  cella), tasso osservato 0. La cella non e’ valutabile (sotto 20
  replicati) e lo shrinkage non e’ una strategia raccomandata; escluso
  quel caso il massimo e’ 0,26 punti (PC6b, M7, 145 soggetti, 10
  occasioni, modello marginalizzato) e nessuna bandiera cambia (PC8, e 0
  cambi su tutte le 696 x 6 celle per entrambi i riferimenti). La
  previsione sbagliava il SE empirico atteso a 145 soggetti e 10
  occasioni (0,05, non 0,03), per questo il valore atteso era 0,5-1,0.
  Non tocca alcuna conclusione.
- **PC10 (indeterminata).** Usabili sui restituiti a 10 pazienti,
  specificazione piena, modello marginalizzato: 74,1% (257 su 347
  restituiti), contro un intervallo di rilevanza da 75% in su. Dei 90
  restituiti ma non usabili, 77 hanno una stima oltre 10 in valore
  assoluto (quasi-separazione con 7 coefficienti e 40 righe) e 13 non
  hanno l’errore standard. La previsione riguardava il meccanismo dei
  fallimenti, non una frase del paper: la frase del paper (il modello
  marginalizzato risponde almeno quanto Fay e Graubard con pochi
  soggetti) vale a maggior ragione (51% contro 23% su 500).
- **PC11 (indeterminata).** Rapporto SE del modello su deviazione
  standard delle stime, modello marginalizzato, specificazione piena, 50
  pazienti: 1,18 +/- 0,07 contro l’intervallo 0,91-1,14. Spiegato
  dall’anomalia dei sottocampioni senza reinserimento (appendice A5):
  con il fattore di popolazione finita sqrt(1 - 50/111) = 0,74 il
  rapporto e’ 0,87. La stessa correzione porta il rapporto della coorte
  neonatale a 50 pazienti da 1,14 a 0,99. Decisione 1 del riassunto.

## Cosa e’ girato

- Commit di pre-registrazione `ed6fd75` (cancelli GC1-GC5, previsioni
  PC1-PC12, decisioni tecniche 1-4, protocollo 1.6). La decisione
  tecnica 5 (colonna `id` per geesmv) e’ stata aggiunta durante il
  lavoro, prima del run sui sottocampioni, dopo il controllo di
  appendice A6.
- Script, in ordine: `R/s7b/s7b_md5.sh before` (GC1) -\>
  `R/s7b/s7b_df.R` (parte 1: t a N - 2, regola, spostamento con
  bootstrap per blocchi; GC2, GC3; 1,2 minuti, 1 processo) -\>
  `R/s7b/s7b_null.R` (parte 2: rifiuto atteso dal nullo marginale, GC4;
  meno di un minuto) -\> `R/s7b/s7b_resp_run.R` (parte 3: dataset
  `respiratory`, GC5, 2 specificazioni x 2.500 estrazioni x 6 metodi;
  circa 4 minuti su 20 processi) -\> `R/s7b/s7b_resp_analyse.R` -\>
  `R/s7b/s7b_preds.R` (tabella delle previsioni, fattore di popolazione
  finita) -\> `R/s7b/s7b_md5.sh after`. Log in
  `reports/morris/S7b_logs/`.
- Nessun dataset simulato e’ stato generato; la parita’ H0/H1 non si
  riesegue (decisione tecnica 4).
- Input riletti: `results_v4/variant_cc/S4/rows_v4.rds` (23.781.670
  righe dei metodi primari), `results_v4/S6b/ns/rows_ns.rds`
  (6.011.274), `tables/V4_S2_degenerate_designs.csv`,
  `tables/V4_S6b_ns_dp.csv`, `tables/V4_S4_truth_by_cell.csv`,
  `tables/V4_S6b_ns_truth.csv`, `rev2/T3_draws.rds`,
  `rev2/T5_tab316.csv`, `rev2/T5_seratio.csv`, `rev2/T7_bias.csv`.
- Tabelle prodotte: `tables/V4_S7b_*.csv` (elenco in appendice A8).

## Appendice

### A1. Regola della Figura 1 con i tre riferimenti (tutti i 1.896 punti di disegno)

n_hi = punti chiaramente sopra 7,5%; n_lo = chiaramente sotto 2,5%; n_na
= non valutabili (meno di 20 replicati usabili); pu = quota usabile
media; upow = potenza per tentativo media. `ref` n = normale, t = t a
N - p, t2 = t a N - 2.

| strategia  | n_hi | n_lo | n_na | tasso_max | pu    | upow  | ammessa | scelta |
|------------|------|------|------|-----------|-------|-------|---------|--------|
| mmm\|t     | 0    | 77   | 0    | 8.442     | 90.98 | 21.21 | TRUE    | TRUE   |
| mmm\|t2    | 0    | 49   | 0    | 8.459     | 90.98 | 21.35 | TRUE    | TRUE   |
| gee_fg\|t  | 2    | 116  | 0    | 10.42     | 91.91 | 20.53 | FALSE   | FALSE  |
| gee_fg\|t2 | 7    | 113  | 0    | 10.65     | 91.91 | 20.66 | FALSE   | FALSE  |
| gee_kc\|t  | 951  | 48   | 4    | 24.9      | 91.86 | 29.12 | FALSE   | FALSE  |
| gee_kc\|t2 | 963  | 45   | 4    | 24.9      | 91.86 | 29.28 | FALSE   | FALSE  |
| gee_md\|t  | 32   | 124  | 0    | 14.14     | 91.85 | 20.64 | FALSE   | FALSE  |
| gee_md\|t2 | 55   | 118  | 0    | 15.43     | 91.85 | 20.78 | FALSE   | FALSE  |
| gee\|t     | 290  | 4    | 0    | 26.04     | 94.07 | 23.35 | FALSE   | FALSE  |
| gee\|t2    | 317  | 3    | 0    | 27.18     | 94.07 | 23.57 | FALSE   | FALSE  |
| glmer\|t   | 265  | 18   | 256  | 19.82     | 80.65 | 23.17 | FALSE   | FALSE  |
| glmer\|t2  | 268  | 15   | 256  | 19.82     | 80.65 | 23.28 | FALSE   | FALSE  |
| shrink\|t  | 340  | 11   | 256  | 24.24     | 81.84 | 23.65 | FALSE   | FALSE  |
| shrink\|t2 | 344  | 7    | 256  | 24.24     | 81.84 | 23.77 | FALSE   | FALSE  |

Per fonte (studio = griglia principale e M1c, 1.512 punti; nuovi =
M6-M7, 384 punti), riferimento t a N - 2:

| strategia  | fonte  | n_dp | n_hi | n_lo | tasso_max |
|------------|--------|------|------|------|-----------|
| gee_fg\|t2 | nuovi  | 384  | 7    | 22   | 10.65     |
| gee_fg\|t2 | studio | 1512 | 0    | 91   | 9.731     |
| gee_kc\|t2 | nuovi  | 384  | 228  | 8    | 20.3      |
| gee_kc\|t2 | studio | 1512 | 735  | 37   | 24.9      |
| gee_md\|t2 | nuovi  | 384  | 24   | 23   | 15.43     |
| gee_md\|t2 | studio | 1512 | 31   | 95   | 13.02     |
| gee\|t2    | nuovi  | 384  | 68   | 0    | 26.74     |
| gee\|t2    | studio | 1512 | 249  | 3    | 27.18     |
| glmer\|t2  | nuovi  | 384  | 78   | 1    | 14.7      |
| glmer\|t2  | studio | 1512 | 190  | 14   | 19.82     |
| mmm\|t2    | nuovi  | 384  | 0    | 8    | 7.967     |
| mmm\|t2    | studio | 1512 | 0    | 41   | 8.459     |
| shrink\|t2 | nuovi  | 384  | 86   | 0    | 16.08     |
| shrink\|t2 | studio | 1512 | 258  | 7    | 24.24     |

A 10 soggetti (126 punti) e fino a 25 soggetti (504 punti), modello
marginalizzato con il controllo e Fay-Graubard:

| strategia  | gruppo       | n_dp | n_hi | n_lo | tasso_max | pu    | upow  |
|------------|--------------|------|------|------|-----------|-------|-------|
| gee_fg\|t  | tutti N\<=25 | 504  | 2    | 79   | 10.42     | 74.17 | 6.476 |
| gee_fg\|t2 | tutti N\<=25 | 504  | 7    | 76   | 10.65     | 74.17 | 6.886 |
| mmm\|t     | tutti N\<=25 | 504  | 0    | 74   | 7.792     | 72.56 | 6.29  |
| mmm\|t2    | tutti N\<=25 | 504  | 0    | 46   | 8.028     | 72.56 | 6.744 |
| gee_fg\|t  | tutti N=10   | 126  | 0    | 21   | 9.312     | 57.96 | 3.116 |
| gee_fg\|t2 | tutti N=10   | 126  | 2    | 19   | 10.37     | 57.96 | 3.876 |
| mmm\|t     | tutti N=10   | 126  | 0    | 53   | 4.914     | 56.45 | 2.246 |
| mmm\|t2    | tutti N=10   | 126  | 0    | 28   | 5.651     | 56.45 | 3.087 |

### A2. Spostamento del tasso di errore di I tipo, N - 2 meno N - p

Media sui punti di disegno dei tassi di cella; MCSE dal bootstrap per
blocchi di replicato (B = 2000, indici di `results_v4/S4/boot_idx.rds` e
`results_v4/S6b/ns/boot_idx_67.rds`).

| metodo | tutti           | N = 10        | N \<= 25      | N \>= 35        |
|--------|-----------------|---------------|---------------|-----------------|
| gee    | 0,212 +/- 0,008 | 1,97 +/- 0,10 | 0,75 +/- 0,03 | 0,017 +/- 0,002 |
| gee_fg | 0,107 +/- 0,006 | 0,85 +/- 0,07 | 0,36 +/- 0,02 | 0,015 +/- 0,002 |
| gee_kc | 0,126 +/- 0,006 | 0,86 +/- 0,07 | 0,41 +/- 0,02 | 0,022 +/- 0,003 |
| gee_md | 0,127 +/- 0,006 | 1,08 +/- 0,08 | 0,43 +/- 0,02 | 0,015 +/- 0,002 |
| glmer  | 0,049 +/- 0,004 | 0,19 +/- 0,03 | 0,14 +/- 0,01 | 0,014 +/- 0,002 |
| mmm    | 0,096 +/- 0,006 | 0,63 +/- 0,08 | 0,31 +/- 0,02 | 0,017 +/- 0,002 |
| shrink | 0,063 +/- 0,005 | 0,29 +/- 0,04 | 0,19 +/- 0,02 | 0,014 +/- 0,002 |

Punti di disegno in cui Fay-Graubard con la t e’ chiaramente sopra 7,5%
(N - p: 2; N - 2: 7):

| strategia  | meccanismo | scenario | N   | T   | usabili | tasso | mcse   |
|------------|------------|----------|-----|-----|---------|-------|--------|
| gee_fg\|t  | GenBinaryY | M7       | 15  | 6   | 864     | 9.838 | 1.013  |
| gee_fg\|t  | GenBinaryY | M7       | 25  | 4   | 912     | 10.42 | 1.012  |
| gee_fg\|t2 | GenBinaryY | M6       | 10  | 6   | 704     | 10.37 | 1.149  |
| gee_fg\|t2 | GenBinaryY | M6       | 20  | 4   | 815     | 9.571 | 1.03   |
| gee_fg\|t2 | GenBinaryY | M7       | 10  | 10  | 815     | 10.31 | 1.065  |
| gee_fg\|t2 | GenBinaryY | M7       | 15  | 4   | 667     | 10.34 | 1.179  |
| gee_fg\|t2 | GenBinaryY | M7       | 15  | 6   | 864     | 10.65 | 1.049  |
| gee_fg\|t2 | GenBinaryY | M7       | 20  | 6   | 935     | 9.519 | 0.9598 |
| gee_fg\|t2 | GenBinaryY | M7       | 25  | 4   | 912     | 10.64 | 1.021  |

### A3. Rifiuto atteso in piu’ dovuto al valore nullo marginale (GLMM-native, arm H0, metodi marginali)

| gruppo | riferimento | n_celle | eccesso_max | eccesso_medio | max | delta/SE |  | n_hi | n_hi_dopo |
|----|----|----|----|----|----|----|----|----|----|
| M1-M5 | n | 2880 | 0.0292 | 2.254e-03 | 0.05048 | 440 | 440 | 0 | 0 |
| M1c | n | 144 | 0.0001923 | 3.243e-05 | 0.004097 | 18 | 18 | 0 | 0 |
| M6-M7 | n | 1152 | 1.545 | 2.226e-02 | 0.3657 | 186 | 186 | 0 | 0 |
| M1-M5 | t | 2880 | 0.02849 | 2.143e-03 | 0.05048 | 319 | 319 | 0 | 0 |
| M1c | t | 144 | 0.000176 | 2.721e-05 | 0.004097 | 6 | 6 | 0 | 0 |
| M6-M7 | t | 1152 | 0.6755 | 2.047e-02 | 0.3657 | 134 | 134 | 0 | 0 |

Per metodo, riferimento t, tutti i gruppi:

| metodo | n_celle | eccesso_max | eccesso_medio | n_hi | n_hi_dopo | cambi |
|--------|---------|-------------|---------------|------|-----------|-------|
| gee    | 696     | 0.25        | 0.006986      | 90   | 90        | 0     |
| gee_fg | 696     | 0.25        | 0.00699       | 0    | 0         | 0     |
| gee_kc | 696     | 0.25        | 0.00699       | 298  | 298       | 0     |
| gee_md | 696     | 0.25        | 0.00699       | 1    | 1         | 0     |
| mmm    | 696     | 0.2534      | 0.007044      | 0    | 0         | 0     |
| shrink | 696     | 0.6755      | 0.007692      | 70   | 70        | 0     |

Le 6 celle con eccesso sopra 0,25 punti (riferimento normale):

| scenario | N   | T   | metodo | usabili | tasso | delta     | SE_emp   | eccesso |
|----------|-----|-----|--------|---------|-------|-----------|----------|---------|
| M7       | 10  | 2   | shrink | 2       | 0     | -0.002818 | 0.007706 | 1.545   |
| M7       | 145 | 10  | mmm    | 1000    | 6.8   | -0.007893 | 0.05247  | 0.2597  |
| M7       | 145 | 10  | gee    | 1000    | 6.3   | -0.007893 | 0.05282  | 0.2563  |
| M7       | 145 | 10  | gee_fg | 1000    | 5.8   | -0.007893 | 0.05282  | 0.2563  |
| M7       | 145 | 10  | gee_kc | 1000    | 19    | -0.007893 | 0.05282  | 0.2563  |
| M7       | 145 | 10  | gee_md | 1000    | 5.8   | -0.007893 | 0.05282  | 0.2563  |

Controprova GC4 (12 celle, 10^6 estrazioni normali):

| scenario | N   | T   | metodo | delta/SE   | formula_n | MC_n    | formula_t | MC_t      |
|----------|-----|-----|--------|------------|-----------|---------|-----------|-----------|
| M7       | 10  | 2   | shrink | -0.3657    | 1.545     | 1.552   | 6.755e-01 | 0.6879    |
| M7       | 145 | 10  | mmm    | -0.1504    | 2.597e-01 | 0.2147  | 2.534e-01 | 0.2044    |
| M7       | 145 | 10  | gee    | -0.1495    | 2.563e-01 | 0.2527  | 2.500e-01 | 0.2435    |
| M7       | 145 | 10  | gee_fg | -0.1495    | 2.563e-01 | 0.2403  | 2.500e-01 | 0.234     |
| M1       | 115 | 2   | mmm    | -0.001358  | 2.112e-05 | -0.0093 | 2.047e-05 | -0.01398  |
| M7       | 10  | 6   | gee_md | -0.01054   | 1.273e-03 | 0.0213  | 5.437e-04 | 0.0005393 |
| M3       | 65  | 10  | gee_kc | -0.02542   | 7.401e-03 | 0.0029  | 6.981e-03 | -0.007593 |
| M1       | 25  | 2   | gee_md | -0.0009977 | 1.140e-05 | 0.0047  | 9.502e-06 | 0.017     |
| M1       | 35  | 8   | mmm    | -0.007146  | 5.849e-04 | 0.0132  | 5.191e-04 | 0.01002   |
| M7       | 135 | 3   | mmm    | -0.01014   | 1.178e-03 | 0.0223  | 1.147e-03 | 0.02292   |
| M4       | 135 | 4   | gee_kc | -0.01424   | 2.324e-03 | 0.02    | 2.263e-03 | 0.01662   |
| M3       | 75  | 3   | gee_md | -0.003387  | 1.314e-04 | -0.0086 | 1.250e-04 | -0.01092  |

### A4. Studio respiratorio (geepack): dataset completo e sottocampioni

Contrasto trattamento attivo contro placebo (log odds). Specificazione
piena: `outcome ~ treat + baseline + visit + age + sex + center`;
ridotta: `outcome ~ treat + baseline`.

| specificazione | metodo | stima | SE     | stato |
|----------------|--------|-------|--------|-------|
| pieno          | MMM    | 1.264 | 0.3453 | ok    |
| pieno          | GEE    | 1.256 | 0.3463 | ok    |
| pieno          | GEE-MD | 1.256 | 0.3679 | ok    |
| pieno          | GEE-KC | 1.256 | 0.3569 | ok    |
| pieno          | GEE-FG | 1.256 | 0.3499 | ok    |
| pieno          | GLMM   | 2.057 | 0.5377 | ok    |
| ridotto        | MMM    | 1.241 | 0.3278 | ok    |
| ridotto        | GEE    | 1.247 | 0.3287 | ok    |
| ridotto        | GEE-MD | 1.247 | 0.3377 | ok    |
| ridotto        | GEE-KC | 1.247 | 0.3331 | ok    |
| ridotto        | GEE-FG | 1.247 | 0.3335 | ok    |
| ridotto        | GLMM   | 2.053 | 0.5331 | ok    |

Estrazioni con una covariata binaria a un solo livello (su 500 per N):

| N   | estrazioni | % una sola modalita’ | % trattamento | % sesso | % centro | % basale |
|-----|------------|----------------------|---------------|---------|----------|----------|
| 10  | 500        | 8.8                  | 0             | 8.8     | 0        | 0.2      |
| 15  | 500        | 1.8                  | 0             | 1.8     | 0        | 0        |
| 20  | 500        | 0.8                  | 0             | 0.8     | 0        | 0        |
| 25  | 500        | 0                    | 0             | 0       | 0        | 0        |
| 50  | 500        | 0                    | 0             | 0       | 0        | 0        |

Fit usabili (stima e SE finiti, \|stima\| \<= 10, SE \< 10) su 500
estrazioni, con la MCSE binomiale, e sui restituiti:

| specificazione | metodo | 10   | 15   | 20   | 25    | 50    |
|----------------|--------|------|------|------|-------|-------|
| pieno          | GEE    | 46,8 | 83,0 | 94,0 | 99,2  | 100,0 |
| pieno          | GEE-FG | 23,4 | 67,8 | 88,2 | 97,6  | 100,0 |
| pieno          | GEE-KC | 27,0 | 69,4 | 88,2 | 97,6  | 100,0 |
| pieno          | GEE-MD | 23,0 | 67,6 | 87,6 | 97,4  | 100,0 |
| pieno          | GLMM   | 53,4 | 88,2 | 96,8 | 100,0 | 100,0 |
| pieno          | MMM    | 51,4 | 85,0 | 95,8 | 100,0 | 100,0 |
| ridotto        | GEE    | 79,2 | 94,4 | 98,4 | 100,0 | 100,0 |
| ridotto        | GEE-FG | 72,4 | 91,6 | 98,4 | 99,8  | 100,0 |
| ridotto        | GEE-KC | 72,4 | 91,6 | 98,4 | 99,8  | 100,0 |
| ridotto        | GEE-MD | 70,6 | 91,6 | 98,4 | 99,8  | 100,0 |
| ridotto        | GLMM   | 78,4 | 94,0 | 98,2 | 99,8  | 100,0 |
| ridotto        | MMM    | 78,4 | 93,2 | 98,2 | 100,0 | 100,0 |

MCSE massima della quota su 500: 2,24 punti.

| specificazione | metodo | 10   | 15   | 20   | 25    | 50    |
|----------------|--------|------|------|------|-------|-------|
| pieno          | GEE    | 52,1 | 84,5 | 94,8 | 99,2  | 100,0 |
| pieno          | GEE-FG | 26,1 | 69,0 | 88,9 | 97,6  | 100,0 |
| pieno          | GEE-KC | 30,1 | 70,7 | 88,9 | 97,6  | 100,0 |
| pieno          | GEE-MD | 25,6 | 68,8 | 88,3 | 97,4  | 100,0 |
| pieno          | GLMM   | 59,1 | 89,8 | 97,6 | 100,0 | 100,0 |
| pieno          | MMM    | 74,1 | 92,6 | 97,6 | 100,0 | 100,0 |
| ridotto        | GEE    | 79,5 | 94,4 | 98,4 | 100,0 | 100,0 |
| ridotto        | GEE-FG | 72,7 | 91,6 | 98,4 | 99,8  | 100,0 |
| ridotto        | GEE-KC | 72,7 | 91,6 | 98,4 | 99,8  | 100,0 |
| ridotto        | GEE-MD | 70,9 | 91,6 | 98,4 | 99,8  | 100,0 |
| ridotto        | GLMM   | 78,4 | 94,0 | 98,2 | 99,8  | 100,0 |
| ridotto        | MMM    | 80,2 | 93,4 | 98,2 | 100,0 | 100,0 |

Modi di fallimento a 10 pazienti:

| specificazione | metodo | restituiti | senza_SE | separati | usabili | errori | timeout |
|----------------|--------|------------|----------|----------|---------|--------|---------|
| pieno          | GLMM   | 452        | 0        | 185      | 267     | 48     | 0       |
| pieno          | MMM    | 347        | 13       | 77       | 257     | 54     | 0       |
| pieno          | GEE    | 449        | 73       | 142      | 234     | 51     | 0       |
| pieno          | GEE-KC | 449        | 309      | 5        | 135     | 309    | 0       |
| pieno          | GEE-FG | 449        | 309      | 23       | 117     | 309    | 0       |
| pieno          | GEE-MD | 449        | 317      | 17       | 115     | 314    | 0       |
| ridotto        | GEE    | 498        | 25       | 77       | 396     | 2      | 0       |
| ridotto        | MMM    | 489        | 10       | 87       | 392     | 4      | 0       |
| ridotto        | GLMM   | 500        | 0        | 108      | 392     | 0      | 0       |
| ridotto        | GEE-KC | 498        | 136      | 0        | 362     | 136    | 0       |
| ridotto        | GEE-FG | 498        | 136      | 0        | 362     | 136    | 0       |
| ridotto        | GEE-MD | 498        | 143      | 2        | 353     | 143    | 0       |

Rapporto SE medio del modello / deviazione standard delle stime (grezzo,
come nella coorte neonatale):

| specificazione | metodo | 10   | 15   | 20   | 25   | 50   |
|----------------|--------|------|------|------|------|------|
| pieno          | GEE    | 0,55 | 0,67 | 0,77 | 0,82 | 1,18 |
| pieno          | GEE-FG | 2,32 | 1,39 | 1,10 | 1,00 | 1,22 |
| pieno          | GEE-KC | 1,11 | 1,01 | 1,01 | 0,99 | 1,27 |
| pieno          | GEE-MD | 2,29 | 1,47 | 1,27 | 1,19 | 1,37 |
| pieno          | GLMM   | 0,71 | 0,70 | 0,71 | 0,78 | 0,94 |
| pieno          | MMM    | 0,57 | 0,64 | 0,75 | 0,82 | 1,18 |
| ridotto        | GEE    | 0,90 | 0,91 | 0,95 | 0,99 | 1,25 |
| ridotto        | GEE-FG | 1,15 | 1,06 | 1,06 | 1,07 | 1,29 |
| ridotto        | GEE-KC | 1,10 | 1,03 | 1,04 | 1,06 | 1,29 |
| ridotto        | GEE-MD | 1,43 | 1,19 | 1,15 | 1,14 | 1,33 |
| ridotto        | GLMM   | 1,12 | 1,10 | 1,04 | 1,06 | 1,19 |
| ridotto        | MMM    | 0,86 | 0,84 | 0,91 | 0,98 | 1,24 |

Distorsione rispetto alla stima sul dataset completo (le correzioni
contro la GEE), con n usabili:

| specificazione | metodo | 10     | 15     | 20     | 25     | 50     |
|----------------|--------|--------|--------|--------|--------|--------|
| pieno          | GEE    | 0,189  | 0,457  | 0,212  | 0,223  | 0,004  |
| pieno          | GEE-FG | 0,130  | 0,334  | 0,183  | 0,235  | 0,004  |
| pieno          | GEE-KC | 0,176  | 0,386  | 0,183  | 0,235  | 0,004  |
| pieno          | GEE-MD | 0,093  | 0,401  | 0,179  | 0,237  | 0,004  |
| pieno          | GLMM   | -0,591 | -0,082 | -0,148 | -0,001 | -0,146 |
| pieno          | MMM    | 0,154  | 0,394  | 0,161  | 0,206  | -0,005 |
| ridotto        | GEE    | 0,054  | 0,152  | 0,063  | 0,117  | -0,012 |
| ridotto        | GEE-FG | 0,064  | 0,139  | 0,063  | 0,119  | -0,012 |
| ridotto        | GEE-KC | 0,064  | 0,139  | 0,063  | 0,119  | -0,012 |
| ridotto        | GEE-MD | 0,040  | 0,139  | 0,063  | 0,119  | -0,012 |
| ridotto        | GLMM   | -0,300 | -0,060 | -0,084 | 0,088  | -0,041 |
| ridotto        | MMM    | 0,031  | 0,135  | 0,047  | 0,113  | -0,016 |

| metodo | N   | n   | stima_media | riferimento | distorsione | mcse    | distorsione/SD |
|--------|-----|-----|-------------|-------------|-------------|---------|----------------|
| GEE    | 10  | 234 | 1.445       | 1.256       | 0.189       | 0.1613  | 0.07661        |
| GEE    | 50  | 500 | 1.26        | 1.256       | 0.00421     | 0.02023 | 0.009308       |
| GEE-FG | 10  | 117 | 1.386       | 1.256       | 0.1302      | 0.1564  | 0.07698        |
| GEE-FG | 50  | 500 | 1.26        | 1.256       | 0.00421     | 0.02023 | 0.009308       |
| GEE-KC | 10  | 135 | 1.431       | 1.256       | 0.1759      | 0.1849  | 0.08189        |
| GEE-KC | 50  | 500 | 1.26        | 1.256       | 0.00421     | 0.02023 | 0.009308       |
| GEE-MD | 10  | 115 | 1.348       | 1.256       | 0.09291     | 0.1876  | 0.04619        |
| GEE-MD | 50  | 500 | 1.26        | 1.256       | 0.00421     | 0.02023 | 0.009308       |
| GLMM   | 10  | 267 | 1.466       | 2.057       | -0.5907     | 0.1539  | -0.2349        |
| GLMM   | 50  | 500 | 1.911       | 2.057       | -0.1457     | 0.03204 | -0.2035        |
| MMM    | 10  | 257 | 1.418       | 1.264       | 0.1541      | 0.1468  | 0.06548        |
| MMM    | 50  | 500 | 1.259       | 1.264       | -0.00504    | 0.02026 | -0.01112       |

Confronto con la coorte neonatale, specificazione piena, 10 e 50
pazienti (stesse definizioni):

| metodo | N | usabili/500 resp | usabili/500 LUS | SE/SD resp | SE/SD LUS | bias resp | bias LUS |
|----|----|----|----|----|----|----|----|
| GEE | 10 | 46.8 | 70 | 0.5515 | 0.8823 | 0.189 | 0.1179 |
| GEE | 50 | 100 | 100 | 1.185 | 1.147 | 0.00421 | 0.01052 |
| GEE-FG | 10 | 23.4 | 34.2 | 2.321 | 1.08 | 0.1302 | 0.2008 |
| GEE-FG | 50 | 100 | 99.8 | 1.225 | 1.16 | 0.00421 | 0.01076 |
| GEE-KC | 10 | 27 | 34 | 1.11 | 1.108 | 0.1759 | 0.2023 |
| GEE-KC | 50 | 100 | 99.8 | 1.274 | 1.172 | 0.00421 | 0.01076 |
| GEE-MD | 10 | 23 | 15 | 2.289 | 1.813 | 0.09291 | 0.1112 |
| GEE-MD | 50 | 100 | 98.8 | 1.371 | 1.173 | 0.00421 | 0.01175 |
| GLMM | 10 | 53.4 | 78.2 | 0.7124 | 0.5478 | -0.5907 | 0.09281 |
| GLMM | 50 | 100 | 100 | 0.9371 | 0.6539 | -0.1457 | 0.006431 |
| MMM | 10 | 51.4 | 77 | 0.5748 | 0.8659 | 0.1541 | 0.1618 |
| MMM | 50 | 100 | 100 | 1.176 | 1.138 | -0.00504 | 0.01525 |

### A5. Sottocampioni senza reinserimento: il fattore di popolazione finita

Le estrazioni di N pazienti su P senza reinserimento hanno una
variabilita’ fra estrazioni ridotta di circa (1 - N/P) rispetto a
campioni indipendenti dalla stessa popolazione (fattore di popolazione
finita della media campionaria, esatto per la media, approssimato per un
coefficiente di regressione). L’errore standard del modello stima la
variabilita’ di campioni indipendenti: il rapporto SE/SD grezzo e’
quindi gonfiato di 1/sqrt(1 - N/P), cioe’ 1,05 a 10 pazienti e 1,35 a 50
nello studio respiratorio (P = 111), 1,03 e 1,15 nella coorte neonatale
(P = 207). Rapporto corretto = rapporto grezzo x sqrt(1 - N/P),
specificazione piena:

| applicazione    | metodo | 10   | 15   | 20   | 25   | 50   |
|-----------------|--------|------|------|------|------|------|
| neonatal cohort | GEE    | 0,86 | 1,01 | 0,91 | 0,99 | 1,00 |
| neonatal cohort | GEE-FG | 1,05 | 1,02 | 0,93 | 1,03 | 1,01 |
| neonatal cohort | GEE-KC | 1,08 | 1,04 | 0,94 | 1,04 | 1,02 |
| neonatal cohort | GEE-MD | 1,77 | 1,23 | 1,00 | 1,06 | 1,02 |
| neonatal cohort | GLMM   | 0,53 | 0,58 | 0,51 | 0,58 | 0,57 |
| neonatal cohort | MMM    | 0,84 | 0,97 | 0,88 | 0,97 | 0,99 |
| respiratory     | GEE    | 0,53 | 0,62 | 0,70 | 0,72 | 0,88 |
| respiratory     | GEE-FG | 2,21 | 1,30 | 1,00 | 0,88 | 0,91 |
| respiratory     | GEE-KC | 1,06 | 0,94 | 0,92 | 0,87 | 0,94 |
| respiratory     | GEE-MD | 2,18 | 1,37 | 1,15 | 1,04 | 1,02 |
| respiratory     | GLMM   | 0,68 | 0,65 | 0,64 | 0,69 | 0,69 |
| respiratory     | MMM    | 0,55 | 0,59 | 0,68 | 0,72 | 0,87 |

Cosa guardare: il rapporto del modello marginalizzato a 50 pazienti
scende sotto 1 in entrambe le applicazioni; quello del modello misto
nella coorte neonatale resta fra 0,51 e 0,58.

### A6. Controllo indipendente delle correzioni di geesmv sul dataset respiratorio

Script `R/s7b/s7b_geesmv_check.R`, log
`reports/morris/S7b_logs/s7b_geesmv_check.log`. Dataset completo,
specificazione piena, SE del coefficiente del trattamento. Con la
colonna `id` originale di `respiratory` (paziente entro centro) geesmv
usava quella come identificativo dei gruppi: Fay-Graubard 0,6348,
Mancl-DeRouen 0,6738, Kauermann-Carroll 0,4509. Con la colonna dei
gruppi chiamata `id` (centro x paziente): Fay-Graubard 0,3499,
Mancl-DeRouen 0,3679, Kauermann-Carroll 0,3569. Calcolo indipendente
sullo stesso fit di geepack (alpha = 0,331): sandwich non corretto
0,3463 (= `clubSandwich` CR0); Mancl-DeRouen 0,3679 a mano = CR3 0,3679;
Kauermann-Carroll = CR2 0,3569; Fay-Graubard a mano (b = 0,75) 0,3499.
La causa e’ nel sorgente di `GEE.var.*`:
`if (is.null(data$id)) data$id <- data[, id]`, cioe’ l’argomento `id` e’
ignorato se una colonna `id` esiste. Nelle simulazioni e nella coorte
neonatale la colonna dei gruppi si chiama `id`.

### A7. Da dove vengono i numeri del riassunto

| frase | fonte | valore |
|----|----|----|
| 0 situazioni chiaramente sopra 7,5 con N - 2 (modello marginalizzato, t); regola: scelto di nuovo | tables/V4_S7b_rule_summary.csv, mmm\|t2, src tutti | n_hi 0; recommended TRUE |
| meno di un quarto di punto in media (spostamento massimo fra i metodi) | tables/V4_S7b_df_shift.csv, group all design points | max shift 0,212 |
| a 10 soggetti +0,6 (marginalizzato) e +0,9 (Fay-Graubard) | tables/V4_S7b_df_shift.csv, group N = 10 | 0,63, 0,85 |
| chiaramente sotto 2,5 a 10 soggetti: 28 su 126 invece di 53 | tables/V4_S7b_rule_summary.csv, src tutti N=10, mmm\|t e mmm\|t2, n_lo | 28, 53 |
| Fay-Graubard t: 7 situazioni invece di 2; tutte GenBinaryY M6-M7, N \<= 25 | tables/V4_S7b_rule_summary.csv gee_fg\|t2 n_hi; V4_S7b_rule_dp.csv hi == TRUE | 7, 2; N max 25; scenari GenBinaryY M6, GenBinaryY M7 |
| eccesso massimo 0,03 (scenari principali) e 0,26 (M6-M7, 145 soggetti, 10 visite) | tables/V4_S7b_null_summary.csv, method all, ref n, max_excess | 0,029; 0,260 |
| nessuna bandiera cambia | tables/V4_S7b_null_summary.csv, n_hi_flip, n_lo_flip | 0 cambi |
| cella shrinkage: 1,5 punti, 2 usabili, tasso 0 | tables/V4_S7b_null_shift.csv, M7 N=10 T=2 shrink | eccesso 1,55, n 2, tasso 0,00 |
| dataset completo: 1,26 (SE 0,35) marginalizzato e GEE; 2,06 (0,54) misto | tables/V4_S7b_resp_ref.csv, spec pieno | MMM 1,26 (0,35), GEE 1,26 (0,35), GLMM 2,06 (0,54) |
| 10 pazienti, piena: 51 / 47 / 53 / 23 / 23 su 100 | tables/V4_S7b_resp_usable.csv, N=10, pct_usable_500 | 51,4 / 46,8 / 53,4 / 23,4 / 23,0 |
| ridotta: circa 8 su 10 | tables/V4_S7b_resp_usable.csv, ridotto, N=10 | 78,4 / 79,2 / 78,4 / 72,4 |
| distorsione marginalizzato 0,15 -\> 0,00; coorte neonatale 0,16 -\> 0,02; misto -0,59 e -0,15 | tables/V4_S7b_resp_bias.csv pieno N=10, 50; tables/V4_S7b_resp_compare.csv lus_bias | 0,15 -\> -0,005; LUS 0,16 -\> 0,02; GLMM -0,59, -0,15 |
| rapporto SE/SD marginalizzato piena 0,58 -\> 1,18; coorte neonatale 0,87 -\> 1,14 | tables/V4_S7b_resp_seratio.csv pieno MMM; V4_S7b_resp_compare.csv lus_seratio | 0,57 -\> 1,18; LUS 0,87 -\> 1,14 |
| fattore 0,74 (50/111) e 0,87 (50/207); rapporti corretti 0,99, 0,87; misto 0,51-0,58 | tables/V4_S7b_seratio_fpc.csv, spec pieno | fpc 0,74, 0,87; MMM corretto LUS 0,99, resp 0,87; GLMM LUS 0,51-0,58 |
| geesmv: 0,635 contro 0,350 | R/s7b/s7b_geesmv_check.R, log reports/morris/S7b_logs/s7b_geesmv_check.log (appendice A6) | 0,6348 / 0,3499 (A6) |
| 8,8% estrazioni a 10 pazienti con sesso a un solo livello | tables/V4_S7b_resp_degeneracy.csv, N=10 | 8,8% |

### A8. Documenti prodotti

`tables/V4_S7b_gates.csv`, `V4_S7b_preds.csv`, `V4_S7b_dp_rates.csv`,
`V4_S7b_rule_dp.csv`, `V4_S7b_rule_summary.csv`, `V4_S7b_df_shift.csv`,
`V4_S7b_df_flips.csv`, `V4_S7b_null_shift.csv`,
`V4_S7b_null_summary.csv`, `V4_S7b_null_flips.csv`,
`V4_S7b_null_gc4.csv`, `V4_S7b_resp_ref.csv`, `V4_S7b_resp_fits.csv`,
`V4_S7b_resp_degeneracy.csv`, `V4_S7b_resp_usable.csv`,
`V4_S7b_resp_seratio.csv`, `V4_S7b_resp_bias.csv`,
`V4_S7b_resp_compare.csv`, `V4_S7b_seratio_fpc.csv`; `results_v4/S7b/`
(md5, `rows_h0_marg.rds`, `resp_draws.rds`); script `R/s7b/`; log
`reports/morris/S7b_logs/`.
