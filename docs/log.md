# Registro delle modifiche

## 2026-11-03 (codice pubblico)

* **Decisione del PI**: [D-92](decisions/D-92.md). Pubblicati `lucavd/clustered_paper` (commit d15be8d, dal commit privato 9237bcc) e `lucavd/fewclust` (commit 2bfaa6e, da fewclust 0.4.2, f05c914) con `tools/public_sync.sh`; controlli superati (nessun file della coorte, nessuna credenziale, nessun file oltre 50 MB).

## 2026-11-03 (S11 chiuso)

* **Step**: [Correzioni dopo la revisione indipendente (S11)](steps/S11.md): approvato dal PI, che conferma i contributi degli autori della title page (opzione raccomandata). Restano rinviati alla preparazione dell'invio il DOI di Zenodo e l'intestazione di provenienza delle costanti di fewclust.

## 2026-11-03 (S11: calcolo, correzioni e riassunto)

* **Step**: [Correzioni dopo la revisione indipendente (S11)](steps/S11.md): regole [D-89](decisions/D-89.md) e [D-90](decisions/D-90.md) applicate a tabelle, figure e testo senza simulazioni nuove (variante delle tabelle cv, [D-91](decisions/D-91.md)); 32 rilievi della revisione con un esito; pipeline di chiusura e passaggi testo (registro 3.376 numeri, 0 discordanze; 7.969 parole). Cancelli GC2, GC4, GC5 superati; GC1 e GC3 non superati come scritti e spiegati con un calcolo; previsione PC3 smentita (perdita di risultati del modello marginalizzato con pendenza casuale anche a 6 occasioni), scritta nella sezione 3.3. fewclust 0.4.1 (regola [Nel controllo di convergenza del modello marginalizzato il fit senza pendenza conta solo se nlm dichiara convergenza; altrimenti nessun risultato (D-90)](decisions/D-90.md)) e 0.4.2 (costanti dalle tabelle finali, commit f05c914, ramo `package`). Approvazione del PI in attesa (contributi degli autori).

## 2026-11-03 (S11: pre-registrazione)

* **Step**: [Correzioni dopo la revisione indipendente (S11)](steps/S11.md): cancelli GC1-GC5 e previsioni PC1-PC8 scritti prima dei calcoli, dopo una sonda in sola lettura (la regola sul fit di riserva toglie il risultato al modello marginalizzato con il controllo in 402.676 dataset; confine ricalcolato uguale a quello salvato in tutti i fit).
* **Decisione tecnica**: [D-91](decisions/D-91.md).

## 2026-11-03 (S10 chiuso)

* **Step**: [Revisione di numeri e stile del manoscritto corretto (S10)](steps/S10.md): approvato dal PI con tre decisioni; report reso.
* **Decisioni del PI**: [D-89](decisions/D-89.md) (dataset separati restano risultati; fallimento silenzioso della GEE nel paper; regola sugli esiti costanti confermata), [D-90](decisions/D-90.md) (fit di riserva del controllo solo se converge; sostituisce [Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente (D-32)](decisions/D-32.md)).
* **Step che riceve dubbi**: [Correzioni dopo la revisione indipendente (S11)](steps/S11.md).

## 2026-11-02 (S10: revisione indipendente e riassunto)

* **Step**: [Revisione di numeri e stile del manoscritto corretto (S10)](steps/S10.md): inventario dei dataset degeneri su tutti i 4.752.000 dataset rigenerati e controprova indipendente di 524 numeri del testo principale, 367 del supplemento e 151 frasi assolute; cancelli GD1-GD4 superati, GD5 superato sul sorgente invariato; previsioni PD1, PD2, PD4 confermate, PD3, PD5-PD8 smentite. Due rilievi cambiano la raccomandazione (dataset separati; fit di riserva non convergenti nel controllo del modello marginalizzato): portati al PI. Rilievi in `tables/V7_S10_findings.csv`.
* **Step che riceve dubbi**: [Correzioni dopo la revisione indipendente (S11)](steps/S11.md).

## 2026-11-02 (S9 chiuso)

* **Step**: [Esito costante come nessun risultato (S9)](steps/S9.md): approvato dal PI con le tre scelte raccomandate. Testo su fewclust aggiunto (sezione 4.1, codice e dati), titolo Word a 18 punti allineato a sinistra; pipeline di chiusura rifatta (run6, run7: registro 3.184 numeri, 0 discordanze). fewclust 0.3.1 chiuso (ramo `package`, commit 626393b).
* **Decisioni del PI**: [D-85](decisions/D-85.md), [D-86](decisions/D-86.md), [D-87](decisions/D-87.md), [D-88](decisions/D-88.md).
* **Step che riceve dubbi**: [Revisione di numeri e stile (S10)](steps/S10.md), con l'indicazione del PI sulle revisioni (inventario dei dataset degeneri per metodo prima della revisione di stile).

## 2026-11-02 (S9: calcolo e riassunto)

* **Step**: [Esito costante come nessun risultato (S9)](steps/S9.md): regola applicata a tutte le strategie con la variante nc della catena delle tabelle e una riga in ogni script per replicato; cancelli GS1-GS5 superati; previsioni PS1-PS4, PS6, PS7 confermate, PS5 smentita in una regione su 19 e PS8 smentita, spiegate. La regola sceglie ora il modello marginalizzato con il controllo; la GEE AR(1) con Fay-Graubard e' fuori taratura in 2 punti con 10 cluster e 2 occasioni. Testo v6 portato nel sorgente e riscritto; esportato come v7. fewclust 0.3.1 (ramo `package`, commit c688ecc). Tre decisioni chieste al PI; in attesa di approvazione.
* **Decisione del PI**: [D-84](decisions/D-84.md).
* **Step che riceve dubbi**: [Revisione di numeri e stile (S10)](steps/S10.md).

## 2026-11-02 (S9: pre-registrazione)

* **Step**: [Esito costante come nessun risultato (S9)](steps/S9.md): cancelli e previsioni scritti prima dei calcoli, dopo un controllo di portata in sola lettura (solo GEE senza correzione e Kauermann-Carroll di fewclust contavano dataset costanti fra le strategie dello studio; nessun sottocampione clinico con esito costante).
* **Decisione tecnica**: [D-83](decisions/D-83.md).

## 2026-11-02 (preparazione all'invio a SMMR e controllo degli esiti costanti)

* **Fuori dagli step**: manoscritto, supplemento, title page e lettera v6 preparati per SMMR (revisione double-anonymized, riferimenti Sage Vancouver, nuova Tabella 2, dichiarazioni e uso dell'AI); artifact del progetto elencati nel file dello step S9.
* **Controllo**: la GEE contava come usabili 9.646 dataset con esito costante; senza di essi la GEE AR(1) con Fay-Graubard e' chiaramente sopra 7,5% in 2 punti. Script in `R/s9/provenance/`, risultato `tables/V6_S9_const_check_table5.csv`.
* **Decisione del PI**: [D-82](decisions/D-82.md).
* **Step aperti**: [Esito costante come nessun risultato (S9)](steps/S9.md), [Revisione di numeri e stile (S10)](steps/S10.md).

## 2026-11-02 (S8b chiuso)

* **Step**: [Attuazione della revisione finale (S8b)](steps/S8b.md): approvato dal PI; le sei fonti lette solo nel riassunto sono state lette per intero dai PDF del PI; due decisioni del PI applicate all'abstract e alla sezione 4.1.
* **Decisione del PI**: [D-81](decisions/D-81.md).

## 2026-11-01 (revisione finale, S8b: attuazione)

* **Step**: [Attuazione della revisione finale (S8b)](steps/S8b.md): 64 voci del piano con esito (60 fatte), correzioni della catena, applicazioni cliniche rifatte, accorciamento a 8.179 parole contate, Tabelle S36-S37, Figura [protocollo datato (S1)](steps/S1.md), errori Monte Carlo nelle didascalie, secondo passaggio del revisore n.2 (5 rilievi bloccanti corretti), revisione di stile. GB1-GB5 superati; in attesa di approvazione del PI.
* **Decisioni tecniche**: [D-73](decisions/D-73.md), [D-74](decisions/D-74.md), [D-75](decisions/D-75.md), [D-76](decisions/D-76.md), [D-77](decisions/D-77.md), [D-78](decisions/D-78.md), [D-79](decisions/D-79.md), [D-80](decisions/D-80.md).
* **Audit**: `docs/audit/MORRIS2019_conformita.csv` con le colonne `esito_S8b`, `evidenza_S8b` (28 su 28 conformi o dichiarati); schede generate in `docs/verification/` e `docs/audit/requirements/`.

## 2026-11-01 (revisione finale, S8 chiuso)

* **Step**: [Revisione avversariale e ri-audit Morris (S8)](steps/S8.md): approvato dal PI sul riassunto In breve; tre dubbi rinviati allo step di attuazione.
* **Decisioni del PI**: [D-70](decisions/D-70.md), [D-71](decisions/D-71.md), [D-72](decisions/D-72.md).
* **Step aperto**: [Attuazione della revisione finale (S8b)](steps/S8b.md).

## 2026-11-01 (revisione finale, S8: revisione, audit e piano)

* **Step**: [Revisione avversariale e ri-audit Morris (S8)](steps/S8.md): revisione avversariale di manoscritto e supplemento v5 (3.517 affermazioni, nove verifiche in parallelo, fino a tre giri), ri-audit dei 28 requisiti di Morris (23 conformi o dichiarati), indagine sui 12 dubbi ricevuti, piano di implementazione ancorato (64 voci) e proposta di accorciamento (9.908 -> 7.824 parole). GS1, GS2, GS3, GS5 superati, GS4 fallito come scritto per un numero e spiegato; PS1, PS2, PS4 smentite, PS3, PS5, PS6 confermate. Manoscritto non modificato. In attesa delle tre decisioni del PI.
* **Decisione tecnica**: [La revisione finale si fa in due step (D-69)](decisions/D-69.md).
* **Audit**: `docs/audit/MORRIS2019_conformita.csv` con le colonne `esito_v8`, `evidenza_v8`, `nota_v8`.

## 2026-11-01 (S7g chiuso)

* **Step**: [Rafforzare il manoscritto v5 prima dell'invio: taratura, denominatori, modello misto, controllo di convergenza (S7g)](steps/S7g.md): approvato dal PI sul riassunto In breve; due dubbi rinviati alla revisione finale (tempo complessivo di calcolo, lunghezza).

## 2026-11-01 (rafforzamento del manoscritto v5, S7g: risultati e testo)

* **Step**: [Rafforzare il manoscritto v5 prima dell'invio: taratura, denominatori, modello misto, controllo di convergenza (S7g)](steps/S7g.md): run con quadratura adattativa (123 punti, 0,92 ore, 0 chunk falliti), GG1, GG2, GG4, GG5 superati, GG3 fallito come scritto e spiegato; PG1-PG3 smentite (la quadratura abbassa l'errore di I tipo del modello misto di 3,1-4,0 punti), PG4 confermata, PG5 smentita e spiegata in parte; manoscritto e supplemento v5 aggiornati (paragrafo 3.10, Tabelle S32-S35), 2.803 numeri attestati, abstract 199 parole; due dubbi rinviati alla revisione finale.
* **Decisione del PI**: [La conclusione sul test di Wald del modello misto si riscrive con l'esito della quadratura adattativa, senza altri calcoli (D-68)](decisions/D-68.md).

## 2026-10-31 (rafforzamento del manoscritto v5, S7g: apertura e pre-registrazione)

* **Step**: [Rafforzare il manoscritto v5 prima dell'invio: taratura, denominatori, modello misto, controllo di convergenza (S7g)](steps/S7g.md): brief corretto su due punti all'apertura (vale l'intervallo esatto di Clopper-Pearson; la quadratura adattativa si fa); disegno del run con quadratura (123 punti di disegno, 1.000 repliche, solo senza effetto del trattamento), 5 cancelli e 5 previsioni scritti prima di ogni calcolo.
* **Decisione del PI**: [La quadratura adattativa per il modello misto con sola intercetta casuale si valuta come analisi di sensibilita' (supera in parte la decisione del 26 ottobre) (D-66)](decisions/D-66.md).
* **Decisione tecnica**: [Il manoscritto v5 si modifica sul posto nello step S7g; non si crea una v6 (D-67)](decisions/D-67.md).

## 2026-10-31 (aggiornamento del manoscritto, S7f chiuso)

* **Step**: [Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md): manoscritto e supplemento v5 con i risultati di S7b, S7d, S7e e S7d2; revisione avversariale (2.249 numeri riletti identici, rilievi corretti, spiegati, decisi o rinviati) e revisione di stile (75 modifiche); Tabelle 2, 3 e 5 ricostruite dalle singole simulazioni; GF1-GF5 superati, PF1-PF4 confermate; condizione per la citazione di fewclust fallita come atteso. 2.297 numeri attestati, 0 discordanze; 8.858 parole. Report `reports/morris/S7f.html`. Sei punti rinviati alla [Revisione avversariale e ri-audit Morris (S8)](steps/S8.md).
* **Decisioni del PI**: [Figura 1 sulle strategie della regola; la GEE AR(1) con Fay e Graubard e t va prima, il modello marginalizzato con controllo e' un'alternativa comparabile (D-61)](decisions/D-61.md), [Chiaramente sopra 7,5% e chiaramente sotto 2,5% si giudicano con l'intervallo esatto di Clopper-Pearson (D-63)](decisions/D-63.md), [Nessuna nuova simulazione dopo la revisione avversariale; i metodi non valutati si dichiarano fra i limiti (D-64)](decisions/D-64.md), [fewclust non si cita come implementazione della raccomandazione; resta nominato come software che ha calcolato le correzioni (D-65)](decisions/D-65.md).
* **Decisione tecnica**: [Revisione di stile applicata; il supplemento v5 diventa la fonte (D-62)](decisions/D-62.md).

## 2026-10-25 (intervallo dal test a inversione di segno, S7d2 chiuso)

* **Decisioni del PI**: [L'intervallo per inversione del test a inversione di segno non si propone nel paper: una frase nella discussione e il dettaglio nel supplemento (D-58)](decisions/D-58.md), [Nel paper il test a inversione di segno non e' calcolabile quando il fit ridotto diverge o lo score e' numericamente nullo (D-59)](decisions/D-59.md). [Intervallo di confidenza dal test a inversione di segno (S7d2)](steps/S7d2.md) chiuso e approvato; il testo da scrivere passa all'[Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md), dubbio ricevuto 6.

## 2026-10-25 (intervallo dal test a inversione di segno, S7d2: risultati)

* **Step**: [Intervallo di confidenza dal test a inversione di segno (S7d2)](steps/S7d2.md): 4.752.000 dataset rigenerati con gli stessi semi (8,0 ore, 0 chunk falliti), dati e p-value del run precedente riprodotti esattamente. GG1 e GG2 passano; GG3 (b) e GG4 falliti come scritti e spiegati. 5 previsioni confermate, 3 smentite e spiegate: copertura del test 94,5%, mai fuori taratura; intervallo finito nel 79% dei dataset e nel 13,6% con 2-3 visite e fino a 25 soggetti; 17% piu' largo della GEE AR(1) con Fay e Graubard e t. Il test del run precedente restituisce p-value casuali dove l'interazione non e' identificata (28% dei dataset con 2-3 visite e fino a 25 soggetti). Report `reports/morris/S7d2.html`; tabelle `tables/V4_S7d2_*.csv`. In attesa di 2 scelte del PI.

## 2026-10-24 (intervallo dal test a inversione di segno, S7d2: controlli prima del run)

* **Step**: [Intervallo di confidenza dal test a inversione di segno (S7d2)](steps/S7d2.md): controprova su 1.080 dataset e prova a secco su 800; quattro correzioni dello strumento prima del run (attese invariate); tempi misurati 0,14 s per dataset. Dubbio aperto: il test del run precedente restituisce p-value casuali dove l'interazione non e' identificata.
* **Decisione tecnica**: [Correzioni dello strumento dell'intervallo prima del run: camminata a passi fissi, partenza da zero, test non calcolabile con fit divergente o score nullo (D-57)](decisions/D-57.md).

## 2026-10-24 (intervallo dal test a inversione di segno, S7d2: pre-registrazione)

* **Step**: [Intervallo di confidenza dal test a inversione di segno (S7d2)](steps/S7d2.md): disegno dell'intervallo (offset nel modello ridotto, segni fissi per dataset, ricerca per raddoppi e bisezione, stima della GEE di indipendenza come centro), valore vero per meccanismo, 4 cancelli e 8 previsioni scritti prima di ogni calcolo.
* **Decisione del PI**: [L'intervallo dal test a inversione di segno si prova su tutte le 2.376 situazioni dei quattro meccanismi (D-56)](decisions/D-56.md).

## 2026-10-24 (terzo meccanismo generatore, S7e chiuso)

* **Decisione del PI**: [Il quarto meccanismo entra nel testo principale: regola della Figura 1 sui quattro meccanismi, frase sul modello misto, dettagli nel supplemento (D-55)](decisions/D-55.md). [Terzo meccanismo generatore (S7e)](steps/S7e.md) chiuso e approvato; il testo da scrivere passa all'[Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md), dubbio ricevuto 5.

## 2026-10-24 (terzo meccanismo generatore, S7e: risultati)

* **Step**: [Terzo meccanismo generatore (S7e)](steps/S7e.md): 960.000 dataset dalla copula gaussiana con correlazione solo seriale (26,4 ore, 0 chunk falliti). GF2-GF4 passano, GF1 fallito come scritto e spiegato (344 dataset con disegno non stimabile senza la riga dell'effetto). 11 previsioni su 11 confermate: la regola della Figura 1 sui 2.376 punti dei quattro meccanismi sceglie ancora la GEE AR(1) con Fay-Graubard e t (92,1%); modello marginalizzato con controllo mai sopra 7,5% e centrato; modello misto 10,2% medio, 219 punti su 480 sopra (errore standard dell'intercetta casuale 0,64-0,81 della dispersione); GEE exchangeable con Fay-Graubard sopra in 5 punti. Report `reports/morris/S7e.html`; tabelle `tables/V4_S7e_*.csv`. In attesa di 1 scelta del PI.

## 2026-10-23 (terzo meccanismo generatore, S7e: pre-registrazione)

* **Step**: [Terzo meccanismo generatore (S7e)](steps/S7e.md): profilo di dipendenza dei tre meccanismi dello studio misurato (tutti con una componente costante per soggetto, correlazione fra visite lontane 0,09-0,14) e generatore a copula gaussiana con margini logistici scritto in `R/v4/rundir/copula/`. Pre-registrazione: 4 cancelli, 11 previsioni, prima di ogni analisi sui dati del nuovo meccanismo.
* **Decisioni**: del PI, [Il quarto meccanismo e' una copula gaussiana con margini logistici e correlazione solo seriale (rho = 0,6) (D-53)](decisions/D-53.md); tecnica, [Organizzazione tecnica del run del quarto meccanismo: semi dello studio, correzioni da fewclust, geesmv sui primi 50 replicati (D-54)](decisions/D-54.md).

## 2026-10-23 (run con i metodi aggiunti, S7d chiuso)

* **Decisioni del PI**: [La regola della Figura 1 si applica alle strategie con stima e intervallo: sceglie la GEE AR(1) con Fay e Graubard e la t; il test a inversione di segno e' complementare (D-50)](decisions/D-50.md) (supera [La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](decisions/D-35.md)), [La correzione di Kauermann e Carroll nel paper viene da fewclust; l'errore di geesmv si spiega nel supplemento (D-51)](decisions/D-51.md), [Uno step dedicato prova l'intervallo di confidenza ottenuto invertendo il test a inversione di segno, prima della riscrittura del paper (D-52)](decisions/D-52.md). [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](steps/S7d.md) chiuso e approvato; nuovo step [Intervallo di confidenza dal test a inversione di segno (S7d2)](steps/S7d2.md); il testo da scrivere passa all'[Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md).

## 2026-10-23 (run con i metodi aggiunti, S7d: risultati)

* **Step**: [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](steps/S7d.md): 3.792.000 dataset rigenerati con gli stessi semi (11,5 ore, 0 chunk falliti); GEE salvata riprodotta esattamente; GE1-GE3 passano, GE4 e GE5 falliti come scritti e spiegati (tolleranza relativa vicino a zero; 19 fit divergenti). 11 previsioni confermate, 1 smentita: la GEE AR(1) con Fay-Graubard di fewclust e la t non e' mai chiaramente sopra 7,5% e risponde piu' del modello marginalizzato (92,2% contro 91,0%). Test a inversione di segno: 5,29%, mai sopra, 99,9% di risposte con 2-3 visite e fino a 25 soggetti. Kauermann-Carroll corretto: 5,55% contro 9,75% di geesmv. Report `reports/morris/S7d.html`; tabelle `tables/V4_S7d_*.csv`. In attesa di 2 scelte del PI.

## 2026-10-22 (run con i metodi aggiunti, S7d: pre-registrazione)

* **Step**: [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](steps/S7d.md): dubbi ricevuti 2(a)-2(d) chiusi con una sonda su 990 dataset (geesmv con AR(1) fallisce per un limite del pacchetto gee; CR2 di clubSandwich non affidabile con AR(1); fewclust 0.2.0 registrato; tempi misurati, generazione compresa). Anomalia: il Kauermann-Carroll di geesmv non usa la radice principale di I - H (controprova a mano). Pre-registrazione: 5 cancelli, 12 previsioni. Protocollo alla versione 1.8.
* **Decisioni del PI**: [geesmv solo sulle prime 50 simulazioni di ogni scenario nel run con i metodi aggiunti (D-48)](decisions/D-48.md), [fewclust 0.2.0 registrato dall'agente nel ramo package cosi' com'era sul portatile del PI (D-49)](decisions/D-49.md).

## 2026-10-22 (pilota dei test con il valore nullo del paper, S7c chiuso)

* **Decisioni del PI**: [Nel run su tutte le situazioni si aggiunge il solo test a inversione di segno sui contributi allo score (D-46)](decisions/D-46.md) e [Test a inversione di segno: principale senza filtro sui coefficienti del modello ridotto, con il filtro come sensibilita' (D-47)](decisions/D-47.md). [Disegno e pilota dei test a livello di soggetto e di permutazione (S7c)](steps/S7c.md) chiuso e approvato; protocollo alla versione 1.7 con il test a inversione di segno; la scelta passa al [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](steps/S7d.md) nel dubbio ricevuto 2(e).

## 2026-10-22 (pilota dei test con il valore nullo del paper, S7c)

* **Step**: [Disegno e pilota dei test a livello di soggetto e di permutazione (S7c)](steps/S7c.md): pre-registrazione (5 cancelli, 9 previsioni) committata prima dei calcoli (`dd8a3e6`); 36.000 dataset rigenerati con gli stessi semi (GEE salvata riprodotta in 33.448 replicati su 33.448). Il test a segni alterni sui contributi allo score e' tarato (5,19% in media, nessuna cella fuori taratura) e risponde piu' spesso di Fay-Graubard; le due analisi a livello di soggetto pre-registrate sono troppo conservative (riassunto non centrato nel soggetto, spiegato con il rapporto di varianza); varianti aggiunte dopo i risultati dichiarate come esplorative. Un cancello (identita' dello score) fallito alla lettera per tolleranza, spiegato. Tempi per il run successivo: circa 3.300 ore di processore, 93% in geesmv. Report `reports/morris/S7c.html`; tabelle `tables/V4_S7c_*.csv`; script `R/s7c/`. In attesa di 2 scelte del PI. Cinque dubbi rinviati a [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](steps/S7d.md).

## 2026-10-22 (rianalisi senza nuove simulazioni, S7b chiuso)

* **Decisioni del PI**: [Il rapporto fra errore standard e dispersione nei sottocampioni si corregge per la popolazione finita (D-44)](decisions/D-44.md) e [Gradi di liberta' della t: N - p nel testo principale, N - 2 di Li e Redden come sensibilita' nel supplemento (D-45)](decisions/D-45.md). [Rianalisi senza nuove simulazioni (S7b)](steps/S7b.md) chiuso e approvato; il testo da scrivere passa allo [Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md) come dubbio ricevuto 3.

## 2026-10-21 (rianalisi senza nuove simulazioni, S7b)

* **Step**: [Rianalisi senza nuove simulazioni (S7b)](steps/S7b.md): protocollo alla versione 1.6 (t a N - 2 accanto a N - p, seconda applicazione clinica su `respiratory` di geepack, rifiuto atteso dal nullo marginale), pre-registrazione (5 cancelli, 12 previsioni) committata prima dei calcoli (`ed6fd75`). Tutti i cancelli superati. Con N - 2 la regola della Figura 1 sceglie di nuovo il modello marginalizzato con il controllo e la t (0 punti chiaramente sopra 7,5% su 1.896); Fay-Graubard con la t chiaramente sopra in 7 punti invece di 2. Il nullo marginale non esattamente zero produce al massimo 0,03 punti di rifiuto in piu' negli scenari principali e 0,26 in M6-M7, nessuna bandiera cambia. Studio respiratorio: risultati in linea con la coorte neonatale; trovato e corretto l'uso da parte di geesmv di qualunque colonna `id` presente; rilevato il fattore di popolazione finita nei rapporti SE/SD dei sottocampioni (decisione per il PI). Report `reports/morris/S7b.html`; tabelle `tables/V4_S7b_*.csv`; script `R/s7b/`. In attesa dell'approvazione del PI.

## 2026-10-21 (piano per le limitazioni)

* **Decisioni**: [Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](decisions/D-42.md) e [fewclust allineato alla regola del paper e citato nel manoscritto (D-43)](decisions/D-43.md). Creati i file degli step [Rianalisi senza nuove simulazioni (S7b)](steps/S7b.md), [Disegno e pilota dei test a livello di soggetto e di permutazione (S7c)](steps/S7c.md), [Run con correlazioni di lavoro, correzioni indipendenti e test a livello di soggetto (S7d)](steps/S7d.md), [Terzo meccanismo generatore (S7e)](steps/S7e.md) e [Aggiornamento del manoscritto con le analisi aggiunte (S7f)](steps/S7f.md); roadmap aggiornata.

## 2026-10-21 (manoscritto v4, commenti del PI nel Word)

* **Aggiornamento**: riportati nel sorgente Quarto i commenti del PI nel Word del manoscritto v4. Abstract e conclusioni ora dicono "GEE with the Fay-Graubard correction with a t reference" (modifica del PI nell'abstract, estesa alla frase delle conclusioni). Formato in [Formato del Word: Arial, tutto nero, testo a 11 punti (D-41)](decisions/D-41.md), con un controllo nella pipeline. Tolti gli spazi doppi lasciati dai commenti di controllo fra una frase e l'altra. Nuovo render: G2 0 violazioni, G3 1.713 numeri con 0 discordanze, G4 86 confronti in regola, abstract 190 parole, testo 6.602 parole, formato 0 scostamenti. Il commento sul pacchetto fewclust e la mitigazione delle limitazioni passano a [Revisione avversariale e ri-audit Morris (S8)](steps/S8.md) come dubbi ricevuti.

## 2026-10-21 (manoscritto, S7 chiuso)

* **Chiusura**: [Riscrittura del manoscritto (S7)](steps/S7.md) approvato dal PI. Scelta 1: il testo v4 e' la base della revisione finale. Scelta 2: [Il manoscritto non data le regole scelte dagli autori (D-40)](decisions/D-40.md). Due dubbi passano a [Revisione avversariale e ri-audit Morris (S8)](steps/S8.md).

## 2026-10-21 (manoscritto, S7: scelta 2 del PI)

* **Decisione**: [Il manoscritto non data le regole scelte dagli autori (D-40)](decisions/D-40.md). Tolte dal manoscritto v4 le frasi su quando sono state scelte la regola della Figura 1 e il criterio della versione precedente (2.6, 4.1, 4.4); la datazione resta nei documenti del progetto come controllo sull'agente. Nuovo render, controlli G2-G4 superati, 6.598 parole.

## 2026-10-21 (manoscritto, S7)

* **Aggiornamento**: [manoscritto (S7)](steps/S7.md) calcolato. Manoscritto e supplemento v4 in `manuscript/v4/` (Quarto, resi in Word e html), con i numeri letti dal registro `tables/V4_S7_numbers.csv` (1.713 voci, voce [V-62](verification/V-62.md) del log di verifica). Cancelli G1-G5 superati; P1 smentita e spiegata (6.684 parole contro 7.000-9.000, due mancanze dei metodi corrette), P2 confermata. Riassunto con due scelte per il PI. Due dubbi rinviati a [revisione finale (S8)](steps/S8.md), file dello step creato.

## 2026-10-21 (apertura del manoscritto, S7)

* **Aggiornamento**: [manoscritto (S7)](steps/S7.md) aperto. Quattro scelte del PI prima della pre-registrazione: [Il manoscritto si riscrive attorno al nuovo messaggio, nei limiti di SMMR (D-36)](decisions/D-36.md), [Tabelle dei risultati con riferimento normale, raccomandazione con riferimento t (D-37)](decisions/D-37.md), [Manoscritto in Quarto con i numeri letti da un registro, reso in Word (D-38)](decisions/D-38.md), [Testo a due cifre significative, tabelle alla cifra giustificata dalla MCSE (D-39)](decisions/D-39.md). Cancelli G1-G5 e previsioni P1-P2 scritti prima di costruire il registro e il testo.

## 2026-10-21 (procedura per il pacchetto)

* **Creazione**: specifica del controllo di convergenza di MMM per il pacchetto fewclust ("procedura Vedovelli", nome scelto dal PI; da implementare dopo il paper), sul ramo `package` in `fewclust/dev/vedovelli_procedure_mmm.md`, con sei dataset di prova e risultati attesi in `fewclust/dev/vedovelli_fixtures/` (script `R/s6b/s6b_vedovelli_fixtures.R`; i fit rifatti coincidono con quelli salvati). Regola e evidenza: [La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](decisions/D-35.md).

## 2026-10-21 (chiusura di S6b)

* **Chiusura di S6b**: approvato dal PI con la regola della Figura 1 ([La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](decisions/D-35.md), PI). Pannello a della Figura 1 rifatto (`reports/morris/variant_cc/S6_figures/`). [manoscritto (S7)](steps/S7.md) riceve un dubbio con cosa scrivere su regola, controllo, soglia e scenari nuovi.

## 2026-10-20 (verifica della regola di ricaduta, S6b)

* **Aggiornamento**: [verifica della regola di ricaduta (S6b)](steps/S6b.md) calcolato. Parte A: la soglia "al bordo" (0,0001 e 0,01 oltre a 0,001) non cambia nulla (PA1-PA6 confermate, `tables/V4_S6b_thr_*.csv`, varianti in `tables/variant_t4/`, `tables/variant_t2/`). Parte B: quattro scenari nuovi M6-M7 ([Nuovo scenario con pendenza grande della deviazione standard latente: griglia intera, scala della popolazione e modello misto (D-34)](decisions/D-34.md), PI), 768.000 dataset; MMM con la regola mai chiaramente sopra 7,5% (PB2), Fay e Graubard con t sopra in 2 punti (PB3 smentita), la regola interviene nel 36-51% dei fit (PB1 smentita, spiegata). Riassunto con una decisione per il PI (regola della Figura 1). Run limitato a 16-28 processi per la temperatura della CPU (fino a 94 C).

## 2026-10-19 (rianalisi di MMM, S6)

* **Chiusura di [tabelle e figure (S6)](steps/S6.md)**: approvato dal PI. La regola della Figura 1 e le due verifiche della regola di ricaduta di MMM (sensibilita' della soglia sui dataset esistenti; nuovo scenario con pendenza vera grande della deviazione standard latente) passano al nuovo step [verifica della regola di ricaduta (S6b)](steps/S6b.md) ([Verifica della regola di ricaduta del modello marginalizzato in uno step dedicato (S6b)](decisions/D-33.md), PI). Il manoscritto ([manoscritto (S7)](steps/S7.md)) dipende ora da S6b.
* **Rianalisi**: MMM in M3-M5 rianalizzato su 1,69 milioni di dataset esistenti con il controllo di convergenza ([Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente (D-32)](decisions/D-32.md), PI). Cancelli G6-G10 superati; P9-P11 confermate, P12 indeterminata e spiegata. Tabelle e figure principali in `tables/variant_cc/` e `reports/morris/variant_cc/`; cambiano classe 30 sotto-affermazioni (`tables/V4_S6cc_claims_diff.csv`), rinvio al manoscritto. Crash del server il 2026-10-18 (errori PCIe sul disco di sistema sotto carico; ventola SYS_FAN7 sotto soglia critica); limite di 50 core.

## 2026-10-18 (apertura di tabelle e figure, S6)
* **Aggiornamento**: causa dell'eccesso di MMM con serie lunghe trovata su 48.000 dataset esistenti rianalizzati (attese P5-P8 registrate prima, commit 20c60aa): pendenza della deviazione standard latente vicino al bordo della parametrizzazione di `MMLB::mm`, due modi di guasto di `nlm`; prova per intervento e regola di ricaduta (P8 confermata). Il PI ha escluso di presentarlo come limite del software; proposta la rianalisi di MMM in M3-M5 con controllo di convergenza.
* **Aggiornamento**: figure senza titoli, spiegazioni in didascalia ([Figure senza titoli sui pannelli, spiegazione in legenda (D-31)](decisions/D-31.md), PI). Domanda del PI su MMM con serie lunghe: diagnostiche post hoc (`tables/V4_S6_mmm_diag_summary.csv`, `V4_S6_mmm_tail.csv`, `V4_S6_mmm_se_probe_summary.csv`, `V4_S6_mmm_outlier.csv`); eccesso di rifiuti da fit di `mm()` non convergenti, l'SE robusto non lo corregge. Nuova decisione aperta per il PI.
* **Aggiornamento**: decisioni del PI [Nessuna simulazione in piu' nello scenario ancorato alla coorte (D-29)](decisions/D-29.md) e [Il nested-loop plot dell'errore di I tipo e' la Figura 2 (D-30)](decisions/D-30.md). Regola della Figura 1 rivista dopo la prima lettura del PI: strategie giudicate per punto di disegno (`tables/V4_S6_rule_safety.csv`), pannello c rifatto; riassunto riscritto con una sola decisione aperta.
* **Aggiornamento**: [tabelle e figure (S6)](steps/S6.md): tabelle v4 con MCSE (griglia principale e M1c), Figura 1 con intervalli Monte Carlo, nested-loop plot e zip plot; affermazioni di M1c classificate, con la sensibilita' bobyqa; riassunto con tre scelte per il PI. Decisione tecnica [Regole tecniche di tabelle e figure (D-28)](decisions/D-28.md). [manoscritto (S7)](steps/S7.md) riceve tre dubbi rinviati.
* **Aggiornamento**: file dello step [tabelle e figure (S6)](steps/S6.md) con cancelli G1-G5, previsioni P1-P4 e regole tecniche di arrotondamento, scritti prima di ogni calcolo.
* **Aggiornamento**: il PI ha approvato lo step [inferenza per blocchi (S5)](steps/S5.md), sul riassunto «In breve» del protocollo v4 con due precisioni (i 2.300 e 1.600-2.000 replicati non riguardano la frase dell'abstract; nello scenario M1c le affermazioni principali che non reggono sono tre, non due).
* **Creazione**: decisione [Riassunto «In breve» scritto per chi non ha seguito il lavoro (D-27)](decisions/D-27.md), dal protocollo di sessione v4; testo approvato dal PI in chat come modello di forma, file non ancora riletto. [Riassunto in parole semplici in testa ai report (D-21)](decisions/D-21.md) passa a superata.

## 2026-10-18 (S5)
* **Aggiornamento**: [inferenza per blocchi (S5)](steps/S5.md) chiuso con le tre scelte del PI: replicati invariati ([Replicati per cella: restano 1.000 anche dopo l'inferenza per blocchi (D-24)](decisions/D-24.md)), celle vinte fuori dal manoscritto ([Il conteggio delle celle vinte esce dal manoscritto (D-25)](decisions/D-25.md)), scenario M1c nel testo principale ([Lo scenario ancorato alla coorte entra nel testo principale (D-26)](decisions/D-26.md)). Log di verifica [V-53](verification/V-53.md) ... [V-61](verification/V-61.md); audit Morris [E4: MCSE delle misure complesse](audit/requirements/E4.md), [F1: n_sim giustificato](audit/requirements/F1.md), [F3: Dipendenza fra dataset simulati dichiarata](audit/requirements/F3.md); il dubbio su M1c in [tabelle e figure (S6)](steps/S6.md) diventa attivo; [manoscritto (S7)](steps/S7.md) riceve la rimozione delle celle vinte.
* **Creazione**: file dello step [inferenza per blocchi (S5)](steps/S5.md) con attese pre-registrate, TL;DR e tre decisioni richieste al PI; decisione tecnica [Regole tecniche dell'inferenza fra celle (D-23)](decisions/D-23.md); file dello step [tabelle e figure (S6)](steps/S6.md) (bozza) con il dubbio rinviato su M1c.
* **Aggiornamento**: [manoscritto (S7)](steps/S7.md) riceve due dubbi rinviati da [inferenza per blocchi (S5)](steps/S5.md) (affermazioni cambiate o trascurabili; sensibilita' bobyqa). Nella roadmap lo stato di [inferenza per blocchi (S5)](steps/S5.md).

## 2026-10-18
* **Aggiornamento**: il PI ha approvato lo step [riaggregazione con MCSE (S4)](steps/S4.md) e le decisioni [Riassunto in parole semplici in testa ai report (D-21)](decisions/D-21.md) e [Protocollo di sessione v3 (D-22)](decisions/D-22.md), che ora risultano verificate.
* **Creazione**: file dello step [manoscritto (S7)](steps/S7.md) (bozza) con il dubbio rinviato sui codici di altri registri ancora senza spiegazione.

## 2026-10-17
* **Creazione**: documentazione convertita in Open Knowledge Format v0.2 con il protocollo di sessione v3 ([Protocollo di sessione v3 (D-22)](decisions/D-22.md)). Le 24 decisioni del registro unico sono ora un file ciascuna in `decisions/`, con titolo e descrizione riletti e approvati dal PI; due decisioni nuove ([Riassunto in parole semplici in testa ai report (D-21)](decisions/D-21.md) e [Protocollo di sessione v3 (D-22)](decisions/D-22.md), non ancora rilette). Un file per step in `steps/`. Frontmatter su roadmap, protocollo, log di verifica, note di sessione e audit; nel protocollo, nel log di verifica e nelle note il testo e' invariato salvo i link ai codici.
* **Aggiornamento**: nella roadmap le regole e il rituale del protocollo precedente sono sostituiti da un rimando al protocollo v3, e la riga sullo stato delle decisioni da un rimando all'indice.
