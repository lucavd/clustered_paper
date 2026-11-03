---
type: Step
title: "Rianalisi senza nuove simulazioni"
description: "Gradi di liberta' della distribuzione t, scarto del valore nullo marginale sotto il modello misto e seconda applicazione clinica su dati pubblici, senza generare nuovi dataset."
step_id: "S7b"
status: stable
pi_approval: { by: human:luca-vedovelli, at: 2026-10-22T07:20:19+02:00, evidence: "chat del 2026-10-22: decisione 1 \"ok\" (D-44), decisione 2 \"no teniamo N-p come principale, e Li redden come sens\" (D-45); commit S7b: closed, approved by PI" }
resource: "../../reports/morris/S7b.html"
generated: { by: mmm-agent/skill-v4, at: 2026-10-21T21:57:34+02:00 }
sources:
  - id: d42
    resource: "../decisions/D-42.md"
    title: "Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale"
  - id: protocol
    resource: "../PROTOCOL_ADEMP.md"
    title: "Protocollo ADEMP della revisione dello studio di simulazione"
  - id: liredden2015
    resource: https://doi.org/10.1002/bimj.201300191
    title: "Li, Redden (2015). Small sample performance of bias-corrected sandwich estimators for cluster-randomized trials with binary outcomes. Statistics in Medicine 34:281-296"
---
# Dubbi ricevuti

1. **Mitigare tre limitazioni senza nuove simulazioni.** Ricevuto il 2026-10-21 dal PI, che ha approvato il piano di [Le limitazioni del paper si mitigano in cinque step dedicati prima della revisione finale (D-42)](../decisions/D-42.md). Le limitazioni sono quelle della sezione 4.4 del manoscritto v4.
   - **Gradi di liberta' della t.** Il paper usa N - p, Li e Redden usano N - 2 (con N soggetti). Stima ed errore standard di ogni replicato sono salvati (per esempio `results_v4/S6cc/cc_all.rds`, colonne `cc_est`, `cc_se`), quindi i p-value con N - 2 si ricalcolano senza rifare fit. Si riapplica la regola della Figura 1 con entrambe le scelte, su tutti i 1.896 punti di disegno.
   - **Valore nullo marginale non esattamente zero sotto il modello misto** (da -0,0030 a -0,0011 negli scenari principali, -0,0079 in M6-M7). Per ogni cella si calcola il rifiuto atteso in piu' dovuto a quello scarto, dall'errore standard empirico, e si dice se puo' spostare una cella oltre 7,5%.
   - **Seconda applicazione clinica.** Dataset `respiratory` del pacchetto geepack (111 pazienti, 4 visite, esito binario, trattamento attivo contro placebo), con lo stesso schema di sottocampionamento e le stesse strategie dell'applicazione sulla coorte neonatale.
   - **Inquadramento dei meccanismi** (ciascuno e' il modello di casa di un metodo): solo testo, si scrive nello step di aggiornamento del manoscritto.

# Piano

- Apertura: aggiornare `docs/PROTOCOL_ADEMP.md` (nuova versione datata) con le due scelte di gradi di liberta' e con la seconda applicazione, prima di ogni calcolo.
- Calcolo: sul server, nessun nuovo dataset simulato.
- Prodotti: `tables/V4_S7b_*.csv`, report `reports/morris/S7b.qmd`.

# Dove stanno i dati che questo step rilegge

- Griglia principale e scenario ancorato alla coorte (1.512 punti di disegno): `results_v4/variant_cc/S4/rows_v4.rds`, una riga per replicato e metodo con stima, errore standard, p-value normale e usabilita' U0; il modello marginalizzato multilivello vi compare gia' con il controllo di convergenza ([Il modello marginalizzato con controllo di convergenza e' la strategia del paper (D-32)](../decisions/D-32.md)). Tabella di riferimento: `tables/variant_cc/V4_S6_rule_safety_dp.csv`.
- Scenari con variabilita' latente crescente M6-M7 (384 punti di disegno): `results_v4/S6b/ns/rows_ns.rds` (metodo `mmm` = con il controllo). Tabella di riferimento: `tables/V4_S6b_ns_rule_safety_dp.csv`.
- Replicati con disegno degenere: `tables/V4_S2_degenerate_designs.csv` (griglia e M1c), colonna `design_degenerate` (M6-M7).
- Valori nulli marginali sotto GLMM-native: `tables/V4_S4_truth_by_cell.csv` (griglia, M1c) e `tables/V4_S6b_ns_truth.csv` (M6-M7), colonna del valore vero marginale con `label = manual`.
- Applicazione sulla coorte neonatale: `rev2/T3_rerun.R` (fit), `rev2/T3_draws.rds` (estrazioni, seme 20260807), `rev2/T5_analyse.R` e `rev2/T7_bias.R` (misure).

# Decisioni tecniche dell'agente (registrate qui, prima dei calcoli)

1. **Gradi di liberta'.** Per ogni strategia di Wald (modello marginalizzato con il controllo, GEE non corretta, Mancl-DeRouen, Kauermann-Carroll, Fay-Graubard, Wald del modello misto, GLM con shrinkage) il p-value con N - 2 gradi di liberta' e' `2 * pt(-abs(stima/SE), df = N - 2)`, sulle stesse righe usabili (U0, disegno non degenere) usate per N - p. Il test del rapporto di verosimiglianza non cambia. Le bandiere della regola (chiaramente sopra 7,5%: tasso - 2 MCSE > 7,5; chiaramente sotto 2,5%: tasso + 2 MCSE < 2,5; non valutabile sotto 20 replicati usabili) sono quelle di `R/s6/s6_rule_safety.R`.
2. **Rifiuto atteso dovuto al nullo marginale.** Per ogni punto di disegno GLMM-native, arm H0 e metodo marginale (modello marginalizzato, le quattro GEE, shrinkage), con scarto delta = valore vero marginale sotto H0 della cella e SE empirico s = deviazione standard delle stime usabili: rifiuto atteso in piu' = Phi(-z + delta/s) + Phi(-z - delta/s) - 0,05 con z = qnorm(0,975) (riferimento normale) e la stessa formula con il quantile della t a N - p gradi di liberta' (riferimento t). E' il rifiuto in piu' che un test esatto con SE pari a quello empirico avrebbe per il solo effetto dello scarto. Si riporta in punti percentuali; si ricalcola la bandiera "chiaramente sopra 7,5%" dopo aver sottratto il rifiuto in piu' dal tasso osservato.
3. **Seconda applicazione clinica.** Dataset `respiratory` di geepack 1.3.13 (identificativo del paziente = centro x id, 111 pazienti, 444 righe, 4 visite ciascuno). Contrasto di interesse: trattamento attivo contro placebo (coefficiente `treatA`, livello di riferimento placebo). Specificazione piena `outcome ~ treat + baseline + visit + age + sex + center`, specificazione ridotta `outcome ~ treat + baseline`, in analogia con le specificazioni "pieno" e "ridotto" della coorte neonatale. Metodi e opzioni come in `rev2/T3_rerun.R`: GEE exchangeable (`geeglm`), tre correzioni geesmv sulla stima puntuale della GEE con timeout di 5 s, modello misto a intercetta casuale con bobyqa, timeout rigidi di 120 s (MMM e GLMM) e 60 s (GEE). Una sola differenza dichiarata: il modello marginalizzato si specifica come nello scenario M1 della simulazione (`lv.formula = ~1`, `t.formula = ~1`), perche' le quattro visite del dataset sono una serie temporale come nella simulazione, mentre nella coorte neonatale le dieci zone polmonari non lo erano e si era usato il solo termine di transizione. Sottocampionamento di pazienti interi a N = 10, 15, 20, 25, 50, 500 estrazioni per N, seme 20260807 e stesso ordine di estrazione della coorte neonatale; riferimento = stima di ciascun metodo sul dataset completo (le correzioni usano il riferimento della GEE). Misure: quota di fit restituiti e usabili (stima e SE finiti, |stima| <= 10, SE < 10) su 500 e sui restituiti; rapporto fra SE medio del modello e deviazione standard delle stime; distorsione rispetto al riferimento; quota di estrazioni con una covariata binaria a un solo livello.
4. **Nessun run di simulazione**: la parita' H0/H1 non si riesegue perche' non si generano dati; i md5 dei risultati grezzi si confrontano con quelli registrati alla chiusura del manoscritto (`results_v4/S7/md5_after.txt`).

5. **Colonna dei gruppi per geesmv** (aggiunta durante il lavoro, prima del run sui sottocampioni). `geesmv::GEE.var.*` usa `data$id` come identificativo dei gruppi ogni volta che quella colonna esiste, qualunque sia l'argomento `id`. Il dataset `respiratory` ha una colonna `id` che numera i pazienti entro centro (1-56): con essa le correzioni uscivano quasi doppie (Fay-Graubard 0,635 contro 0,350 del calcolo indipendente). Nello script la colonna dei gruppi si chiama `id` (= centro x paziente) e i tre errori standard corretti coincidono con `clubSandwich` (CR2 = Kauermann-Carroll, CR3 = Mancl-DeRouen) e con una formula scritta a mano (Fay-Graubard) alla quarta cifra. Simulazioni e coorte neonatale usano gia' `id` come colonna dei gruppi e non sono toccate.

# Attese pre-registrate (scritte prima di ogni calcolo)

## Cancelli

| codice | controllo | passa se |
|---|---|---|
| GC1 | integrita' dei risultati grezzi | i md5 di `results_full/` coincidono con `results_v4/S7/md5_after.txt`; i md5 di `results_v4/variant_cc/S4/rows_v4.rds`, `results_v4/S6b/ns/rows_ns.rds` e delle tabelle rilette sono uguali prima e dopo lo step |
| GC2 | riproduzione dei conteggi con N - p | i conteggi di rifiuto con la t a N - p gradi di liberta', ricalcolati da stima ed errore standard, riproducono esattamente `n`, `nd`, `rate`, `hi`, `lo` di `tables/variant_cc/V4_S6_rule_safety_dp.csv` e di `tables/V4_S6b_ns_rule_safety_dp.csv` per ogni strategia e punto di disegno (1.896 punti); il p-value normale ricalcolato da stima ed errore standard riproduce quello salvato entro 1e-9 su ogni riga usabile di Wald |
| GC3 | completezza | ogni strategia di Wald ha 1.896 punti di disegno con il riferimento a N - 2; nessuna riga usabile perde il p-value nel ricalcolo |
| GC4 | controprova della formula del rifiuto atteso | su 12 celle (le quattro con lo scarto piu' grande e otto scelte a caso con seme 20261021) la formula analitica concorda entro 0,05 punti percentuali con una simulazione di 10^6 statistiche normali con media delta/s |
| GC5 | dataset e schema di estrazione | `respiratory` ha 111 pazienti, 444 righe, 4 visite per paziente, md5 `281e0d3d43ddca923b9debb6e2e996f1` del suo export CSV; le 2.500 estrazioni di pazienti rigenerate con il seme 20260807 e la stessa procedura, applicate agli identificativi della coorte neonatale, riproducono `rev2/T3_draws.rds` |

## Previsioni

Margini: regola della skill (aggregato +/- 1 punto, cella +/- 2,5 punti) dove non motivato altrimenti. Le bandiere "chiaramente sopra" contengono gia' la MCSE: i conteggi di punti di disegno si giudicano come numeri interi, senza intervallo di confidenza.

| codice | quantita' | atteso | intervallo di rilevanza | perche' conta |
|---|---|---|---|---|
| PC1 | punti di disegno (su 1.896) in cui il modello marginalizzato con il controllo e la t a N - 2 e' chiaramente sopra 7,5% | 0 | 0 (un solo punto cambia la frase "mai chiaramente sopra 7,5%" del paper) | la raccomandazione regge con la scelta di Li e Redden |
| PC2 | esito della regola della Figura 1 con N - 2 | la stessa strategia (modello marginalizzato con il controllo, t) | qualsiasi altra strategia scelta e' smentita | Figura 1 |
| PC3 | spostamento medio del tasso di errore di I tipo, N - 2 meno N - p, media sui 1.896 punti, per ogni strategia di Wald | da +0,1 a +0,4 punti (N - 2 ha piu' gradi di liberta' di N - p, quindi rifiuta di piu'; con p = 4 o 5 la differenza conta solo sotto 35 soggetti) | da 0 a +1 punto | se superasse 1 punto la scelta dei gradi di liberta' entrerebbe nei metodi come scelta sensibile |
| PC4 | spostamento medio ai 10 soggetti (media sui punti con N = 10), per il modello marginalizzato con il controllo e per Fay-Graubard | da +0,5 a +1,5 punti | da 0 a +2,5 punti | e' dove il paper dice che il modello marginalizzato e' conservativo con la t |
| PC5 | punti in cui Fay-Graubard con la t a N - 2 e' chiaramente sopra 7,5% | da 2 a 6 (con N - p sono 2, entrambi in M7 su GenBinaryY) | da 2 a 10 | la frase "chiaramente sopra solo con pochi soggetti e variabilita' che cresce" resta vera se restano pochi e in M6-M7 |
| PC6 | rifiuto atteso in piu' dovuto al nullo marginale, massimo su tutti i punti GLMM-native (griglia, M1c, M6-M7) e i metodi marginali | da 0,5 a 1,0 punti (in M7 con 145 soggetti e 10 visite: scarto 0,008, SE empirico circa 0,03) | da 0 a 1,5 punti (1,5 punti e' quanto serve per portare una cella al 6% oltre il 7,5%) | lo scarto non puo' da solo produrre una bandiera |
| PC7 | rifiuto atteso in piu', massimo negli scenari principali M1-M5 e M1c | sotto 0,15 punti | da 0 a 0,5 punti | la limitazione si puo' dichiarare trascurabile sulla griglia principale |
| PC8 | bandiere "chiaramente sopra 7,5%" su GLMM-native che cambiano sottraendo il rifiuto in piu', per il modello marginalizzato con il controllo (t) e per Fay-Graubard (t) | 0 | 0 | Tabella 4 e Figura 1 |
| PC9 | quota di estrazioni a 10 pazienti di `respiratory` con una covariata binaria (trattamento, sesso, centro, baseline) a un solo livello | da 2% a 10% (ipergeometrica con 54/57 trattati, 56/55 centri, baseline 50/61; il sesso e' sbilanciato) | da 0% a 25% | nella coorte neonatale la stessa quota era il motivo principale dei fit mancanti |
| PC10 | quota di fit usabili sui restituiti, a 10 pazienti, specificazione piena, modello marginalizzato | da 85% a 100% (nella coorte neonatale i fallimenti venivano dalle covariate degeneri, qui rare) | da 75% a 100% | il paper dice che con pochi soggetti il modello marginalizzato risponde quanto Fay-Graubard |
| PC11 | rapporto SE del modello / deviazione standard delle stime, modello marginalizzato, specificazione piena, a 50 pazienti | da 0,91 a 1,14 | da 0,91 a 1,14 (margine della skill per il SE) | nella coorte neonatale il rapporto era sotto 1 a tutti gli N |
| PC12 | distorsione del modello marginalizzato rispetto alla stima sul dataset completo, specificazione piena, a 50 pazienti, divisa per la deviazione standard delle stime | da -0,2 a +0,2 | da -0,46 a +0,46 | la distorsione a 50 pazienti e' una frase del paper |

# In breve

## In breve: le tre limitazioni che si potevano affrontare senza nuove simulazioni cambiano qualcosa nel paper?

**A cosa serviva.** Il paper riscritto dichiara tre limitazioni che si potevano verificare sui risultati già salvati, senza generare nuovi dati: i gradi di libertà della distribuzione t (noi usiamo il numero di soggetti meno i coefficienti del modello, Li e Redden usano il numero di soggetti meno 2); il fatto che, nei dati generati dal modello misto, il valore nullo dell'effetto sulla scala della popolazione non è esattamente zero; e il fatto che l'applicazione clinica usa una sola coorte. Abbiamo ricalcolato i test con l'altra scelta di gradi di libertà su tutte le 1.896 situazioni simulate, quantificato cella per cella il rifiuto in più dovuto a quel valore nullo, e ripetuto l'applicazione clinica su un dataset pubblico (lo studio respiratorio di geepack: 111 pazienti, 4 visite, trattamento attivo contro placebo).

**Cosa abbiamo trovato.**

1. **Con i gradi di libertà di Li e Redden la raccomandazione non cambia: il modello marginalizzato multilivello con il controllo di convergenza e la distribuzione t resta l'unica strategia mai chiaramente sopra 7,5 falsi positivi su 100, e la regola della Figura 1 lo sceglie di nuovo.**
   - Passare al numero di soggetti meno 2 alza i falsi positivi di meno di un quarto di punto in media sulle 1.896 situazioni, per ogni test di Wald; a 10 soggetti li alza di 0,6 punti per il modello marginalizzato e di 0,9 per la correzione di Fay e Graubard.
   - A 10 soggetti il modello marginalizzato resta conservativo, ma meno: chiaramente sotto 2,5 falsi positivi su 100 in 28 situazioni su 126, invece di 53.
   - La correzione di Fay e Graubard con la t diventa chiaramente sopra 7,5 in 7 situazioni invece di 2: tutte con dati generati sulla scala della popolazione, variabilità latente che cresce nel tempo e al più 25 soggetti. La frase del paper che la descrive resta vera, con i numeri aggiornati.
2. **Il valore nullo non esattamente zero del modello misto non può spostare nessuna cella: il rifiuto in più che produce è al massimo 0,03 falsi positivi su 100 negli scenari principali e 0,26 negli scenari con variabilità crescente (145 soggetti, 10 visite).**
   - Sottraendolo dai tassi osservati, nessuna delle bandiere "chiaramente sopra 7,5" o "chiaramente sotto 2,5" cambia, per nessun metodo che stima l'effetto sulla scala della popolazione.
   - Un solo valore grande (1,5 punti) viene da una cella in cui lo stimatore con shrinkage ha 2 soli risultati usabili, entrambi collassati a zero: la cella non entra nelle valutazioni perché sotto i 20 risultati usabili.
   - La limitazione può diventare una frase quantitativa.
3. **Nello studio respiratorio l'effetto del trattamento è ben identificato (log odds 1,26, errore standard 0,35, per il modello marginalizzato e la GEE; 2,06 sulla scala condizionale per il modello misto) e nei sottocampioni si ripetono i risultati della coorte neonatale, con una differenza.**
   - Con 10 pazienti e il modello con sei covariate, un risultato usabile arriva dal modello marginalizzato in 51 estrazioni su 100, dalla GEE in 47, dal modello misto in 53, da Fay e Graubard e da Mancl e DeRouen in 23: le correzioni falliscono per un errore numerico che non restituisce l'errore standard, come nella coorte neonatale. Con il modello ridotto (trattamento e stato basale) tutti i metodi rispondono in circa 8 estrazioni su 10.
   - La distorsione del modello marginalizzato rispetto alla stima sull'intero dataset scende da 0,15 con 10 pazienti a zero con 50, come nella coorte neonatale (da 0,16 a 0,02). La stima condizionale del modello misto resta invece sotto quella dell'intero dataset di 0,59 con 10 pazienti e di 0,15 con 50.
   - La differenza: il rapporto fra errore standard del modello e dispersione delle stime nei sottocampioni, per il modello marginalizzato con sei covariate, va da 0,57 con 10 pazienti a 1,18 con 50 (nella coorte neonatale da 0,87 a 1,14); con 10-25 pazienti e sette coefficienti le stime sono instabili e la loro dispersione è gonfiata da stime grandi che il filtro di usabilità tiene.
4. **Un dettaglio di calcolo che tocca entrambe le applicazioni: la dispersione delle stime nei sottocampioni è più piccola del vero perché i sottocampioni sono estratti senza reinserimento da una coorte finita.**
   - Con 50 pazienti su 111 la dispersione va divisa per 0,74; con 50 su 207 per 0,87. Il rapporto errore standard su dispersione, che il paper riporta per la coorte neonatale, è quindi gonfiato soprattutto a 50 pazienti: per il modello marginalizzato passa da 1,14 a 0,99 nella coorte neonatale e da 1,18 a 0,87 nello studio respiratorio; per il modello misto nella coorte neonatale da 0,55-0,65 a 0,51-0,58.
   - Le conclusioni del paper non cambiano (il modello misto sottostima l'errore standard in ogni caso), ma i numeri della frase sì.
5. **Un errore di strumento evitato: il pacchetto geesmv usa come identificativo dei gruppi qualunque colonna si chiami `id`, se esiste, ignorando quella indicata.**
   - Nello studio respiratorio esiste una colonna `id` che numera i pazienti dentro ciascun centro: le correzioni risultavano quasi doppie. Rinominata la colonna, le tre correzioni coincidono entro la quarta cifra con un calcolo indipendente (il pacchetto clubSandwich e una formula scritta a mano).
   - Simulazioni e coorte neonatale non sono toccate: lì la colonna dei gruppi si chiama già `id`.

**Cosa devi decidere.**

1. **Nelle due applicazioni riportiamo il rapporto errore standard su dispersione corretto per la popolazione finita? Consiglio di sì.**
   - L'errore standard del modello stima la variabilità di campioni indipendenti; la dispersione fra sottocampioni senza reinserimento la sottostima di un fattore noto (radice di 1 meno la frazione estratta). Correggere costa una frase nei metodi e cambia i numeri della frase sul rapporto nella sezione dei risultati.
   - Alternativa: tenere il rapporto grezzo e aggiungere una frase di cautela sul fattore. Non cambia nulla nel calcolo, ma il lettore vede un rapporto sopra 1 a 50 pazienti che non è un difetto dei metodi.
   - Alternativa: togliere il rapporto a 50 pazienti e tenerlo solo fino a 25, dove il fattore vale al più 0,94. Si perde il punto in cui il rapporto si avvicina a 1.
2. **Quale scelta di gradi di libertà va nel testo principale? Consiglio di tenere la nostra (soggetti meno coefficienti) e di aggiungere quella di Li e Redden come analisi di sensibilità, con una frase e una tabella nel supplemento.**
   - La raccomandazione è la stessa con entrambe; con quella di Li e Redden il modello marginalizzato è meno conservativo a 10 soggetti e Fay e Graubard è chiaramente sopra la soglia in 7 situazioni invece di 2.
   - Alternativa: passare a Li e Redden come scelta principale. Allinea il paper al lavoro che citiamo, ma tutti i numeri della Tabella 4 e della Figura 1 vanno rigenerati con quella scelta.

**Comunque:**

- quando aggiorneremo il manoscritto, le tre limitazioni diventeranno frasi con i numeri di questo step: la frase sui gradi di libertà, la frase quantitativa sul valore nullo del modello misto, e un paragrafo con una tabella nel supplemento per lo studio respiratorio (con la differenza dichiarata: lì il modello marginalizzato include la componente latente, come nella simulazione, perché le visite sono una serie temporale);
- nello stesso aggiornamento entrerà anche la frase che inquadra i due meccanismi generatori come modello di casa, ciascuno, di un metodo;
- la nota sul pacchetto geesmv e la colonna `id` andrà nella documentazione del codice e nel supplemento, accanto alla verifica numerica delle correzioni.

Appena decidi questi 2 punti, lo step è chiuso.

Il PI ha deciso i 2 punti il 2026-10-22: step chiuso.

# Scelte del PI sul riassunto (2026-10-22)

- **Scelta 1, correzione per la popolazione finita: si'.** Vedi [Il rapporto fra errore standard e dispersione nei sottocampioni si corregge per la popolazione finita (D-44)](../decisions/D-44.md).
- **Scelta 2, gradi di liberta': N - p resta nel testo principale, N - 2 di Li e Redden come sensibilita', breve nel testo e con una tabella nel supplemento.** Vedi [Gradi di liberta' della t: N - p nel testo principale, N - 2 di Li e Redden come sensibilita' nel supplemento (D-45)](../decisions/D-45.md). Alla domanda del PI se con Li e Redden il modello marginalizzato "va meglio per la regola": per i due criteri della regola e' uguale (0 punti chiaramente sopra, stessa quota usabile), ma e' meno conservativo a 10 soggetti e i concorrenti peggiorano.

# Dubbi rinviati

Allo [Aggiornamento del manoscritto con le analisi aggiunte (S7f)](S7f.md) passano, come testo da scrivere e non come dubbi aperti: la frase sui gradi di liberta' (N - p nel testo, N - 2 come sensibilita' breve nel testo e con tabella nel supplemento, [Gradi di liberta' della t: N - p nel testo principale, N - 2 di Li e Redden come sensibilita' nel supplemento (D-45)](../decisions/D-45.md)); la frase quantitativa sul valore nullo marginale (massimo 0,03 punti negli scenari principali, 0,26 in M6-M7, `tables/V4_S7b_null_summary.csv`); il paragrafo e la tabella supplementare dello studio respiratorio (`tables/V4_S7b_resp_*.csv`); il rapporto SE su dispersione con il fattore di popolazione finita in entrambe le applicazioni ([Il rapporto fra errore standard e dispersione nei sottocampioni si corregge per la popolazione finita (D-44)](../decisions/D-44.md), `tables/V4_S7b_seratio_fpc.csv`); l'inquadramento dei meccanismi come modello di casa di un metodo; la nota su geesmv e la colonna `id`.
