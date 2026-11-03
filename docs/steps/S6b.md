---
type: Step
title: "Verifica della regola di ricaduta del modello marginalizzato"
description: "Sensibilita' della soglia della regola di ricaduta di MMM sui dataset esistenti, nuovo scenario con pendenza vera grande della deviazione standard latente, e scelta della regola della Figura 1."
step_id: "S6b"
status: stable
resource: "../../reports/morris/S6b.html"
pi_approval: { by: human:luca-vedovelli, at: 2026-10-21T10:00:00+02:00, evidence: "scelta del PI in chat (prima alternativa, \"ok va bene\"); commit S6b: closed, approved by PI" }
generated: { by: mmm-agent/skill-v4, at: 2026-10-19T18:00:00+02:00 }
sources:
  - id: s6
    resource: "S6.md"
    title: "Tabelle e figure secondo Morris"
  - id: d32
    resource: "../decisions/D-32.md"
    title: "Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente"
  - id: d33
    resource: "../decisions/D-33.md"
    title: "Verifica della regola di ricaduta del modello marginalizzato in uno step dedicato (S6b)"
---
# Dubbi ricevuti

1. **Sensibilita' della soglia "al bordo".** Rinviato il 2026-10-19 da [tabelle e figure (S6)](S6.md) ([Verifica della regola di ricaduta del modello marginalizzato in uno step dedicato (S6b)](../decisions/D-33.md)). La regola di [Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente (D-32)](../decisions/D-32.md) dichiara un fit "al bordo" se `nlm` non converge (codice diverso da 1 e 2) oppure se una pendenza stimata della deviazione standard latente, exp(a_k), e' sotto 0,001. La soglia 0,001 e' stata scelta guardando i fit della diagnostica di [tabelle e figure (S6)](S6.md). Cosa serve, contato su `results_v4/S6cc/cc_all.rds` (1.685.204 fit usabili di M3-M5):
   - soglia 0,0001: nessun fit nuovo; 21.879 fit che convergono e hanno 0,0001 <= exp(a_k) < 0,001 tornano al fit dello studio (colonne `A_*` gia' salvate);
   - soglia 0,01: 19.962 fit che convergono e hanno 0,001 <= exp(a_k) < 0,01 richiedono il fit senza pendenza, che non e' stato calcolato (gli altri 1.178.907 fit al bordo ce l'hanno gia'). Stessi dataset, rigenerati con il loro seme come in `R/s6/s6_mmm_cc_run.R`.
   Per chiuderlo: rifare con le due soglie le misure della variante (errore di I tipo e potenza di MMM in M3-M5, punti chiaramente sopra 7,5%, classi delle sotto-affermazioni), con lo stesso motore, e dire se qualche conclusione di [tabelle e figure (S6)](S6.md) dipende dalla soglia.

2. **Nuovo scenario con pendenza vera grande della deviazione standard latente.** Rinviato il 2026-10-19 da [tabelle e figure (S6)](S6.md) ([Verifica della regola di ricaduta del modello marginalizzato in uno step dedicato (S6b)](../decisions/D-33.md)). In tutti i generatori dello studio la variabilita' latente cambia poco nel tempo: GenBinaryY con sigma = (1; 0,01), cioe' deviazione standard latente 1 + 0,01 t (`R/v4/rundir/generation_model3.R`, M3-M4; in M5 anche 0,01 per la dose); GLMM-native e GLMM-AR1 con pendenza casuale nel tempo di varianza 0,001, deviazione standard 0,03 (`R/v4/rundir/manual/generation_model3.R`, `R/v4/rundir/ar1/generation_models_ar1.R`). Il modello di ricaduta senza pendenza e' quindi quasi giusto e la regola e' favorita per costruzione. Per chiuderlo: simulare uno scenario in cui il modello senza pendenza e' sbagliato e verificare se la regola resta tarata. Le scelte di disegno spettano al PI e si portano all'apertura (sezione seguente).

3. **Regola della Figura 1.** Rinviata il 2026-10-19 da [tabelle e figure (S6)](S6.md) su scelta del PI. Le opzioni e la loro evidenza sono nel riassunto e nell'appendice A3 di `reports/morris/S6.html` (`tables/variant_cc/V4_S6_rule_safety.csv`): Fay-Graubard con riferimento t come regola, con MMM corretto citato come alternativa ugualmente tarata (raccomandazione di [tabelle e figure (S6)](S6.md)); MMM corretto come prima scelta; una via di mezzo (FG|t in generale, MMM corretto con N >= 35 e serie corte); il criterio del v3 alla lettera. Per chiuderlo: il PI sceglie dopo l'esito dei punti 1 e 2; si registra la decisione e si rifa' il pannello a della Figura 1.

# Scelte di disegno da portare al PI all'apertura

Per il punto 2, prima della pre-registrazione (al massimo tre opzioni, con una raccomandazione e il costo stimato da tempi misurati):

- quale struttura: copia di M3 (pendenza solo nel tempo) con la sola pendenza cambiata, oppure anche M5;
- quale pendenza vera e su quali generatori: per esempio GenBinaryY con sigma = (1; 0,1) e (1; 0,25), GLMM-native con varianza della pendenza 0,01 e 0,04; GLMM-AR1 si', no;
- quale griglia: intera (16 N x 6 T) oppure un sottoinsieme (per esempio N = 10, 25, 55, 145 x T = 4, 6, 8, 10), bracci H0 e H1, 1.000 replicati;
- quali metodi: almeno MMM come nello studio, MMM con la regola e Fay-Graubard con riferimento t; tutti gli otto metodi primari se il nuovo scenario deve entrare nelle tabelle del paper;
- semi: un nuovo indice di scenario (per esempio `20260805 + 1000*6 + i`) con L'Ecuyer-CMRG, registrato come decisione, perche' cambia lo schema dei semi;
- parita' H0/H1 e test dei parametri del generatore estesi al nuovo scenario prima del run; al massimo 50 core, un thread per processo (note del server).

# Attese pre-registrate, parte A: soglia della regola di ricaduta (2026-10-19, prima di ogni calcolo)

Cosa si fa. La regola di [Il modello marginalizzato si analizza con un controllo di convergenza negli scenari con pendenza della deviazione standard latente (D-32)](../decisions/D-32.md) dichiara un fit "al bordo" se `nlm` non converge (codice diverso da 1 e 2) oppure se una pendenza stimata della deviazione standard latente, exp(a_k), e' sotto 0,001. Si rifanno le misure con due soglie: 0,0001 (variante `t4`) e 0,01 (variante `t2`); il riferimento e' la variante con 0,001 (`cc`), oggi principale. Stessi dataset di M3-M5 (seme `20260805 + 1000*model + i`, L'Ecuyer-CMRG), stesso motore.

- Soglia 0,0001: nessun fit nuovo; i fit che convergono con 0,0001 <= exp(a_k) < 0,001 tornano al fit dello studio (colonne `A_*` di `results_v4/S6cc/cc_all.rds`).
- Soglia 0,01: per i fit che convergono con 0,001 <= exp(a_k) < 0,01 si calcola il fit senza pendenza (`lv.formula = ~1`), con lo script `R/s6b/s6b_thr_run.R`, copia della parte di fit di `R/s6/s6_mmm_cc_run.R`. Se il fit senza pendenza fallisce, resta il fit dello studio, come nella regola.
- Le due varianti passano dalla stessa catena di S4-S6 (`MMM_VARIANT=t4`, `MMM_VARIANT=t2`) fino alla tabella di sicurezza della regola (`R/s6/s6_rule_safety.R`); il confronto con `cc` e' in `R/s6b/s6b_thr_compare.R`.

## Cancelli

- **GA1 (parita').** `R/testing/test_h0h1_parity.R` passa prima del run.
- **GA2 (integrita').** md5 di `results_full/` (28.800 file), di `results_v4/S6cc/cc_all.rds` e di `results_v4/variant_cc/S4/rows_v4.rds` identici prima e dopo lo step.
- **GA3 (stesso motore).** (a) Per i fit da rifare con la soglia 0,01, il fit dello studio rifatto ha stima, SE e p identici a `A_est`, `A_se`, `A_p` salvati (|differenza| < 1e-8). (b) Su 2.000 fit gia' al bordo estratti a caso (seme 20261019), il fit senza pendenza rifatto ha stima e SE identici a `B_est`, `B_se` salvati. Zero discordanze.
- **GA4 (costruzione).** Lo script che costruisce le varianti, con soglia 0,001, riproduce stima, SE e p della strategia salvata (`cc_est`, `cc_se`, `cc_p`) in tutti i 1.685.204 fit usabili e un insieme di righe identico a `results_v4/variant_cc/S4/rows_v4.rds`. I conteggi dell'inventario si ritrovano esatti: 21.879 fit tornano al fit dello studio con 0,0001; 19.962 richiedono il fit senza pendenza con 0,01.
- **GA5 (completezza e confini).** Ogni fit da rifare ha un esito registrato; la catena gira senza errori per `t4` e `t2`; fuori dalle righe di MMM in M3-M5 l'insieme di righe di ciascuna variante e' identico a quello dello studio.

## Previsioni

Ragionamento: le due soglie spostano circa l'1% dei fit di M3-M5 (21.879 e 19.962 su 1.685.204). Abbassare la soglia rimette in gioco fit con pendenza piccola ma non nulla, vicini al bordo e quindi in parte con stima spostata: ci si aspetta un piccolo aumento dei falsi positivi con serie lunghe. Alzarla manda al fit senza pendenza fit in cui la pendenza vera (circa 0,01 per occasione) e' quasi nulla: atteso un effetto quasi nullo. Tutte le differenze sono variante meno `cc`, appaiate sugli stessi dataset; MCSE dal bootstrap per blocchi sui replicati (B = 2000, seed 20261015, matrice `results_v4/S4/boot_idx.rds`).

- **PA1.** Differenza di errore di I tipo di MMM, media sulle 864 celle H0 di M3-M5: attesa +0,1 punti (0,0001) e 0,0 (0,01); intervallo di rilevanza da -1 a +1 punto (margine di default per la media sulle celle).
- **PA2.** La stessa sulle 288 celle con T = 8-10: attesa +0,4 punti (0,0001) e -0,1 (0,01); intervallo da -1 a +1.
- **PA3.** Differenza di potenza di MMM, media sulle 864 celle H1 di M3-M5: attesa +0,1 punti (0,0001) e -0,1 (0,01); intervallo da -1 a +1 (la frase del paper e' che la potenza resta praticamente la stessa).
- **PA4.** Celle H0 di M3-M5 in cui la differenza appaiata supera +/- 2,5 punti oltre 2 MCSE (margine cella per cella): attese 0 con entrambe le soglie; intervallo 0-8, come per i punti chiaramente sopra 7,5% della rianalisi di [tabelle e figure (S6)](S6.md). Le celle smentite si spiegano.
- **PA5 (decisiva per la Figura 1).** Punti di disegno (1.512, griglia principale e scenario della coorte) in cui MMM con riferimento t e' chiaramente sopra 7,5% (tasso - 2 MCSE > 7,5): con 0,001 sono 0; attesi 0 con entrambe le soglie; intervallo 0-0, perche' anche un solo punto escluderebbe la strategia dalla regola "prima la taratura". Si riporta anche MMM con riferimento normale (1 punto con 0,001).
- **PA6.** Sotto-affermazioni del paper (156) che cambiano classe rispetto a `cc`: attese 0 con entrambe le soglie; intervallo 0-0. Ogni cambiamento si spiega uno per uno (soglia di rilevanza attraversata o intervallo che passa da un lato all'altro).

Controllo dei margini contro la MCSE: nella rianalisi di [tabelle e figure (S6)](S6.md) le differenze appaiate medie avevano MCSE 0,03-0,04 punti (`tables/V4_S6cc_preds.csv`), con 70% dei fit cambiati; qui ne cambia l'1%, quindi 2 MCSE attese sotto 0,1 punti, molto sotto il margine di 1 punto. Nessuna previsione e' indeterminata per costruzione.

## Parte B: nuovo scenario con pendenza vera grande

Le scelte di disegno si portano al PI con i tempi misurati (prova di sola durata `R/s6b/s6b_timing_pilot.R`, nessuna misura di prestazione); le attese della parte B si scrivono qui dopo la decisione e prima del run.

## Parte B: attese pre-registrate (2026-10-19, dopo la scelta del PI e prima del calcolo)

Disegno: [Nuovo scenario con pendenza grande della deviazione standard latente: griglia intera, scala della popolazione e modello misto (D-34)](../decisions/D-34.md). Scenari M6 e M7 (copie di M3 con la sola pendenza cambiata), GenBinaryY e GLMM-native, 768 celle per 1.000 replicati. Script: `R/s6b/s6b_ns_common.R` (generatori, un replicato, fit della regola), `R/s6b/s6b_ns_run.R` (run), `R/s6b/s6b_ns_agg_lib.R` e `R/s6b/s6b_ns_agg.R` (misure, con le definizioni di `R/s5/s5_build.R`, `R/s6/s6_dp_rates.R`, `R/s6/s6_rule_safety.R`). Nelle misure "MMM" e' la strategia con la regola; il fit dello studio resta come "MMM senza controllo".

### Cancelli

- **GB1 (parita' e generatori).** Per ogni scenario nuovo, con N = 10, 55, T = 3, 10 e due semi: (a) H1 e H0 passano al motore gli stessi argomenti salvo il coefficiente testato (H0 = 0); (b) il generatore nuovo passa gli argomenti del generatore v4 di M3 con lo stesso seme, salvo la pendenza; (c) la pendenza vale quella di [Nuovo scenario con pendenza grande della deviazione standard latente: griglia intera, scala della popolazione e modello misto (D-34)](../decisions/D-34.md). `R/testing/test_h0h1_parity.R` passa (gia' verificato per la parte A, si ripete).
- **GB2 (stesso motore).** Il motore nuovo, con indice di seme 3 e il generatore di M3 non modificato, riproduce sui replicati 1-20 di quattro celle di M3 (GenBinaryY N = 25, T = 4; GLMM-native N = 55, T = 8; due bracci) ogni riga dello studio (stima, SE, p, utilizzabilita', nove metodi) e i fit della regola salvati in `results_v4/S6cc/cc_all.rds`; il fit A della regola coincide con la riga MMM. Zero discordanze.
- **GB3 (completezza).** 7.680 blocchi da 100 replicati; ogni replicato ha un esito o una riga di traccia per ogni metodo.
- **GB4 (controprova delle misure).** Il codice di aggregazione di S6b, applicato alle righe di M3 della variante principale, riproduce `tables/variant_cc/V4_S6_dp_rates.csv` (n, nd, k, kt) e `V4_S6_rule_safety_dp.csv` (tasso, esiti, quota utilizzabile, potenza) in tutti i punti di M3 di GenBinaryY e GLMM-native, e i valori veri marginali di GLMM-native M3 di `tables/V4_S4_truth_by_cell.csv` (|differenza| < 1e-10).
- **GB5 (bootstrap e integrita').** Le matrici di indici per k = 1-5 rigenerate con il seme 20261015 coincidono con `results_v4/S4/boot_idx.rds`; md5 di `results_full/` invariati alla fine dello step (GA2).

### Previsioni

Ragionamento: nello studio la regola interveniva in circa il 70% dei fit (1.178.907 su 1.685.204 al bordo), perche' la pendenza vera, circa 0,01, sta vicino al bordo 0 della parametrizzazione di `MMLB::mm`. Con pendenza vera grande la stima resta lontana dal bordo: la regola dovrebbe intervenire solo quando `nlm` non converge, soprattutto con N piccolo. Il modello senza pendenza, usato in quei casi, e' sbagliato solo nella struttura di dipendenza; la media marginale resta quella giusta. Ci si aspetta quindi una regola quasi neutra e un MMM tarato. Per Fay e Graubard con riferimento t non c'e' motivo di attendersi cambiamenti: la correzione non dipende dal modello di dipendenza.

- **PB1.** Quota di fit MMM in cui la regola usa il fit senza pendenza, braccio H0, tutte le celle, per scenario: attesa 15% (GenBinaryY 0,1), 5% (GenBinaryY 0,25), 35% (GLMM-native 0,01), 15% (GLMM-native 0,04); intervallo di rilevanza 0-40%. Oltre il 40% la regola sostituirebbe un modello giusto con uno sbagliato nella maggior parte dei dataset, e la frase del paper su cosa fa la regola cambierebbe.
- **PB2 (decisiva per la Figura 1).** Punti di disegno dei quattro scenari nuovi (384) in cui MMM con la regola e riferimento t e' chiaramente sopra 7,5% (tasso - 2 MCSE > 7,5): attesi 0; intervallo 0-0. Anche un solo punto lo escluderebbe dalla regola "prima la taratura".
- **PB3 (decisiva per la Figura 1).** Lo stesso per Fay e Graubard con riferimento t: attesi 0; intervallo 0-0.
- **PB4.** Errore di I tipo di MMM con la regola, media sulle 96 celle H0 di ciascuno scenario: atteso 5,3%; intervallo 4-6% (come nella rianalisi di [tabelle e figure (S6)](S6.md)).
- **PB5.** Differenza di errore di I tipo MMM con la regola meno MMM senza controllo, media sulle 96 celle H0 di ciascuno scenario: attesa -0,3 punti (GenBinaryY) e -0,5 (GLMM-native); intervallo da -1 a +1 punto.
- **PB6.** Differenza di potenza MMM con la regola meno MMM senza controllo, media sulle 96 celle H1 di ciascuno scenario: attesa 0; intervallo da -1 a +1 punto.

Controllo dei margini contro la MCSE: la MCSE di una cella e' circa 0,7 punti (tasso 5%, 1.000 replicati); la media su 96 celle con numeri casuali comuni ha MCSE stimata 0,1-0,2 punti, quindi 2 MCSE sotto 0,4, contro un margine di 1 punto. Le differenze appaiate hanno MCSE piu' piccole (0,03 nella rianalisi di [tabelle e figure (S6)](S6.md)). PB2 e PB3 sono conteggi su giudizi gia' fatti con la MCSE del punto.

Riportati senza previsione: tassi, copertura, bias e rapporto fra SE del modello e SE empirico di tutti i metodi; il tasso di rifiuto di MMM nei soli fit in cui interviene la regola; il valore vero marginale sotto H0 nel modello misto.
# Esito delle attese

Parte A (soglia "al bordo"):

- Cancelli GA1-GA5 superati: parita' H0/H1 (8/8); md5 di `results_full/` (28.800 file) e dei due ingressi di S6cc identici prima e dopo (`tables/V4_S6b_md5_check.csv`); 19.962 fit dello studio e 2.000 fit senza pendenza rifatti con 0 discordanze; con soglia 0,001 la costruzione riproduce `cc_est`, `cc_se`, `cc_p` in 1.685.204 fit e l'insieme di righe della variante principale; conteggi 21.879 e 19.962 ritrovati; fuori da MMM in M3-M5 le righe delle due varianti sono identiche a quelle dello studio (`tables/V4_S6b_thr_gates.csv`, `V4_S6b_thr_integrity.csv`).
- PA1-PA6 confermate con entrambe le soglie: differenze medie di errore di I tipo e potenza fra -0,0012 e +0,0007 punti (MCSE <= 0,0011); 0 celle oltre +/- 2,5 punti; MMM con riferimento t chiaramente sopra 7,5% in 0 punti di disegno su 1.512; 0 sotto-affermazioni su 156 cambiano classe (`tables/V4_S6b_thr_preds.csv`). Spiegazione: nei fit che cambiano campo la stima con e senza pendenza coincide quasi (mediana |differenza| / SE 2e-5 e 8e-4; decisioni del test cambiate in 23-38 fit su circa 10.000 per braccio, nei due versi; `tables/V4_S6b_thr_flips.csv`).

Parte B (scenari M6-M7, [Nuovo scenario con pendenza grande della deviazione standard latente: griglia intera, scala della popolazione e modello misto (D-34)](../decisions/D-34.md)):

- GB1, GB2, GB4, GB5 superati (`tables/V4_S6b_ns_gates_pre.csv`), dopo la correzione di due errori degli strumenti di controllo: GB1 confrontava formule con ambienti diversi (sempre "diverse"); GB4 ometteva le combinazioni metodo-punto senza righe (n = 0). Nessun numero dello studio e' cambiato.
- GB3 fallito come scritto e superato con il controllo corretto (`tables/V4_S6b_ns_gates_post.csv`): 365 replicati su 768.000 non avevano la riga di un'analisi; rifatti, tutte le 516 analisi mancanti restituiscono un fit senza il termine time:trt (disegno collineare, quasi tutti con N = 10), che lo script, come lo studio, non salva; contano come non utilizzabili, come nello studio. Errore dello strumento di controllo, non dei dati.
- PB2 confermata: MMM con la regola e riferimento t chiaramente sopra 7,5% in 0 punti su 384.
- PB3 smentita: Fay e Graubard con riferimento t chiaramente sopra 7,5% in 2 punti su 384, entrambi GenBinaryY con pendenza 0,25 (N = 15, T = 6: 9,8%, MCSE 1,0; N = 25, T = 4: 10,4%, MCSE 1,0). Spiegata: in quello scenario con N <= 25 e T >= 3 Fay e Graubard con t sta fra 6,7% e 10,4% in tutti i 20 punti, con pendenza 0,1 fra 5,2% e 9,1% (`tables/V4_S6b_ns_diag_fg.csv`): la correzione per piccoli campioni non basta quando la varianza latente cresce molto nel tempo. E' un risultato e informa la scelta della regola.
- PB1 smentita in tre scenari su quattro: la regola usa il fit senza pendenza nel 44,3%, 35,6%, 51,4%, 43,6% dei fit H0 (atteso 5-35%). Spiegata (`tables/V4_S6b_ns_diag_fallback.csv`): soprattutto `nlm` non converge (codici 3-5: 26-31% dei fit, fino al 58% con T = 2-3 e N <= 55); in piu', con T = 4-6, il 13-36% dei fit converge ma con pendenza stimata sotto 0,001. Con N >= 65 e T = 8-10 la quota scende al 6-37%. Nei fit in cui la regola interviene il tasso di rifiuto sotto H0 e' 5,0-5,2% (per gruppi di N e T fra 1,4% e 6,8%): la regola sostituisce spesso un modello giusto con uno sbagliato nella dipendenza, ma senza perdere la taratura.
- PB4 confermata in tre scenari (5,27%, 5,21%, 5,65%), indeterminata in GenBinaryY con pendenza 0,25 (5,84% +/- 0,16). Spiegata: l'eccesso viene da T = 8-10 con N <= 55 (media 5,1-6,2%), nessun punto chiaramente sopra 7,5% (`tables/V4_S6b_ns_diag_mmm.csv`); non tocca una conclusione, perche' la regola della Figura 1 giudica punto per punto.
- PB5 confermata (da -0,59 a -0,89 punti), PB6 confermata (da -0,07 a -0,30 punti) (`tables/V4_S6b_ns_preds.csv`).
- Riportato senza previsione: MMM senza controllo chiaramente sopra 7,5% in 8 punti su 384 (6 in GLMM-native con varianza 0,04 e T = 10); valore vero marginale sotto H0 in GLMM-native da -0,002 a -0,008 (`tables/V4_S6b_ns_truth.csv`).
- Controllo del codice di aggregazione su 300 blocchi prima della fine del run (`S6B_TEST=1`, uscita in una cartella temporanea): l'ultima riga stampata mostrava valori parziali delle previsioni. Le attese erano gia' committate; nessuna e' stata cambiata.
- Calcolo: 768.000 dataset in 23,1 ore (16,7 + 6,4, un rilancio automatico dopo il limite di tempo del lavoro), con 16-28 processi invece di 48: la CPU del server arrivava a 94 C con circa 33 processi; un controllo automatico ha tenuto la temperatura sotto 90 C (`reports/morris/S6b_logs/workers_changes.log`).
- Evidenza per la regola su tutte le 1.896 situazioni (1.512 dello studio, 384 nuove; `tables/V4_S6b_rule_all.csv`): MMM con la regola e t: 0 chiaramente sopra 7,5%, 77 chiaramente sotto 2,5%, risponde nel 91,0%, potenza per tentativo 21,2%; Fay e Graubard con t: 2, 116, 91,9%, 20,5%.

# In breve

## In breve: la regola che corregge il modello marginalizzato regge, e quale metodo raccomandiamo nella Figura 1?

**A cosa serviva.** Nello step precedente abbiamo corretto il modello marginalizzato multilivello con una regola: se l'algoritmo non converge, o se la stima di quanto cresce nel tempo la variabilità tra soggetti finisce sul limite (zero), si usa il modello senza quella crescita. La regola era stata costruita su dati in cui la variabilità cresce pochissimo, dove il modello di ricaduta è quasi giusto. Abbiamo controllato due cose: se il risultato dipende dalla soglia che decide "sul limite", e se la regola resta tarata quando la variabilità cresce davvero e il modello di ricaduta è sbagliato. Da qui dipende quale metodo raccomandare nella Figura 1.

**Cosa abbiamo trovato.**

1. **La soglia della regola non conta.**
   - Con una soglia dieci volte più piccola o dieci volte più grande cambia il modello usato in circa 20.000 dei 1,69 milioni di dataset già simulati; lì i due modelli danno quasi la stessa stima.
   - Falsi positivi e capacità di trovare un effetto vero si spostano al più di circa un millesimo di test su 100, cioè in modo trascurabile; nessuna delle 156 frasi di confronto del paper cambia.
2. **Quando la variabilità tra soggetti cresce davvero, il modello marginalizzato corretto resta tarato in tutte le 384 situazioni nuove.**
   - Abbiamo simulato 768.000 dataset in quattro scenari nuovi: dati generati sulla scala della popolazione, con deviazione standard latente che in 10 visite passa da 1 a 1,9 oppure a 3,25 (nello studio arriva a 1,09); dati generati dal modello misto, con una pendenza casuale nel tempo 10 o 40 volte più variabile che nello studio.
   - Con la distribuzione t non supera mai chiaramente 7,5 falsi positivi su 100; in media ne fa da 5,2 a 5,8 su 100, a seconda dello scenario.
   - Senza la correzione sbaglierebbe: supera chiaramente la soglia in 8 situazioni, quasi tutte con 10 visite nei dati del modello misto.
3. **La regola interviene molto più spesso di quanto mi aspettassi, ma senza danno.**
   - Usa il modello senza crescita in un dataset su tre fino a uno su due, soprattutto perché l'algoritmo non converge con 2-3 visite o con pochi soggetti; con molti soggetti e 8-10 visite succede in 6-37 dataset su 100.
   - Dove interviene, i falsi positivi restano 5 su 100, e la capacità di trovare un effetto vero cala di meno di mezzo test su 100.
4. **La correzione di Fay e Graubard con la distribuzione t, che oggi è la regola consigliata, esce di taratura in uno scenario nuovo.**
   - Sulla scala della popolazione, con la deviazione standard latente che arriva a 3,25 e fino a 25 soggetti, fa da 7 a 10 falsi positivi su 100; in due situazioni (15 soggetti e 6 visite, 25 soggetti e 4 visite) supera chiaramente 7,5, con 10 falsi positivi su 100.
   - Mettendo insieme le simulazioni dello studio e quelle nuove (1.896 situazioni), il modello marginalizzato corretto con la t è l'unica strategia che non supera mai chiaramente la soglia.

**Cosa devi decidere.**

1. **Quale metodo raccomandiamo nella Figura 1? Consiglio: il modello marginalizzato multilivello con il controllo di convergenza e la distribuzione t, a ogni numero di soggetti e di visite, con la correzione di Fay e Graubard con la t citata come alternativa.**
   - Il principio è quello che hai già accettato: chi analizza i dati non sa come sono stati generati, quindi la raccomandazione deve reggere in tutte le situazioni simulate. Si escludono i metodi che anche una sola volta superano chiaramente 7,5 falsi positivi su 100; fra quelli rimasti si sceglie quello che risponde più spesso. Con gli scenari nuovi ne resta uno solo.
   - Il prezzo. Risponde un po' meno spesso (91 dataset su 100 contro 92) e la regola è più lunga, perché va spiegato il controllo. In cambio è troppo prudente in meno situazioni (77 contro 116) e trova un effetto vero un po' più spesso (21 volte su 100 contro 20,5). Con 2-3 visite e fino a 25 soggetti non funziona bene nessuno dei due: rispondono in 4 dataset su 10 e trovano un effetto vero meno di 2 volte su 100.
   - Cosa cambia nel paper. La raccomandazione si rovescia rispetto a quella proposta nello step precedente. Fay e Graubard si cita come alternativa ugualmente semplice, con l'avvertimento che con pochi soggetti e variabilità che cresce molto nel tempo può dare fino a 10 falsi positivi su 100.
   - Alternativa: tenere Fay e Graubard con la t come regola e presentare gli scenari nuovi come prova di resistenza nel materiale supplementare. La regola resta in una frase senza controllo di convergenza, ma non vale più "in tutte le situazioni simulate": serve l'avvertimento per pochi soggetti e variabilità crescente, e il paper deve dire che in quel caso il modello marginalizzato corretto è più sicuro.

**Comunque, quando riscriveremo il paper:**
- nei metodi una frase dirà che la soglia della regola non cambia i risultati, con le due soglie provate;
- i quattro scenari nuovi andranno nel materiale supplementare, con le stesse tabelle degli altri e la figura dei falsi positivi in ogni situazione; nel testo diremo cosa fa la regola, quanto spesso interviene e perché, e che senza di essa il modello marginalizzato sbaglia anche qui;
- nei dati del modello misto con pendenza casuale grande, l'assenza di effetto sulla scala della popolazione vale solo in modo approssimato (valore vero fino a -0,008, contro 0,2 dell'effetto cercato): lo diremo, come già per lo studio;
- rifaremo il pannello della Figura 1 con la regola che scegli.

Appena decidi questo punto, lo step è chiuso.

**Sui tempi di calcolo.** La simulazione è durata 23 ore (tempo misurato), invece delle 8-11 che avevo stimato con una prova di cronometraggio su 48 processori. Il motivo: la CPU del server arrivava a 94 gradi con circa 33 processi, e un controllo automatico ha tenuto il run fra 16 e 28 processi, quasi sempre 20-22. Il server sembra raffreddarsi poco: con 128 processori, oltre i 25 processi circa scalda troppo. Al crash del 18 ottobre era stata segnalata una ventola sotto la soglia critica: vale la pena farla controllare.

# Decisioni del PI

1. **Regola della Figura 1** (2026-10-21): il modello marginalizzato multilivello con controllo di convergenza e riferimento t, a ogni numero di soggetti e di visite; Fay e Graubard con t come alternativa, con l'avvertimento per pochi soggetti e variabilita' crescente ([La regola della Figura 1 raccomanda il modello marginalizzato con controllo di convergenza e riferimento t (D-35)](../decisions/D-35.md)). Pannello a rifatto in `reports/morris/variant_cc/S6_figures/S6_fig1.png` (`fig1_rule.json`; nel pannello c la linea in evidenza e' ora MMM con t); i dati del pannello sono gli stessi (`S6_figures_data.csv` identico a meno dell'ordine delle righe).

# Dubbi rinviati

Nessun dubbio rinviato da questo step. La riscrittura delle frasi del paper e del pannello della Figura 1 spetta al manoscritto ([il manoscritto (S7)](S7.md)), come già previsto.
