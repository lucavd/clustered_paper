---
type: Report
title: Revisione sulla letteratura del manoscritto v5 (S8, track Letteratura)
description: Confronto con la letteratura letta in sessione, punti nascosti, analisi proposte ed eseguite sulle righe esistenti, frasi da riformulare in 4.3 e 4.5.
status: draft
generated: {by: mmm-agent/skill-v4, at: 2026-11-01T12:30:00+01:00}
version: 2 (giro 2 di contestazione)
---

# Revisione sulla letteratura (S8, track Letteratura)

Manoscritto revisionato: `manuscript/v5/MMM_v5.qmd` (md5 2f0c4bb9a2e2b61f2da4ed3cacc11a70), commit cc192d7, ramo `morris-v4`. Sessione del track: frame bb9fafd3-5b98-4bbc-a29e-5d7f9bb8c332. Lavoro sul server in `~/MMM_rev_S8/Letteratura/` (nessuna scrittura in `~/MMM`).

## Sintesi

1. **Due frasi di sintesi vanno oltre cio' che lo studio mostra.**
   - L'abstract raccomanda le due strategie senza dire che il contrasto e' un'interazione trattamento x tempo, che varia dentro il soggetto. Bie et al. 2021 mostrano che il modello marginalizzato si comporta diversamente con un trattamento assegnato al cluster (LIT-003, bloccante).
   - Le conclusioni dicono che con 2-3 occasioni e al piu' 25 cluster "no analysis estimates the effect reliably". Ma il 54% dei dataset stimabili fallisce con le firme numeriche della separazione, e la GEE penalizzata non e' stata valutata, come dice 4.5 stesso (LIT-007, bloccante).
2. **Il verdetto "mai chiaramente sopra 7,5%" e' incerto in una replica, la graduatoria no.**
   - Ricampionando i blocchi di replicato con gli indici del progetto, la GEE AR(1) con Fay-Graubard resta senza flag in 38 ricampioni su 100, il modello marginalizzato in 18, entrambi in 7.
   - Il calcolo e' pessimista. Con tassi contratti, che e' un calcolo ottimista, le probabilita' salgono a 93 e 80 su 100.
   - In tutti i calcoli le due strategie hanno il minor numero atteso di flag. Rilievo ridimensionato a minor (LIT-004).
3. **La scelta di AR(1) rispetto all'indipendenza e' coerente con la regola, ma lo scarto e' piccolo.** L'indipendenza ha 2 punti con flag e la stessa quota usabile, e Thompson et al. 2021 (SMMR) la preferiscono. Rilievo minor (LIT-005).
4. **La sezione 4.3 non usa il suo argomento piu' forte e caratterizza male due fonti.**
   - Argomento non usato: Fay-Graubard con riferimento t e' anche la scelta della letteratura recente (Thompson 2021 negli stepped-wedge con esiti binari; Hemming 2025 per esiti continui).
   - Prima fonte: McNeish [31] e' uno studio di simulazione, non una rassegna, e conclude a favore di un metodo che il paper non valuta.
   - Seconda fonte: la frase sul modello marginalizzato ignora Bie et al. [18], che lo studiano proprio in piccoli campioni.
5. **I limiti (4.5) sono completi rispetto alla decisione registrata (D-64), non rispetto alla letteratura.** Nel primo giro avevo scritto che mancavano analisi richieste da D-64: e' sbagliato e l'ho ritirato (sezione 10). Mancano gli stimatori del modello misto raccomandati con pochi cluster (McNeish 2016; Hemming 2025, letti per intero) e la correzione di Morel.

Conteggi del log (giro 2): 51 righe, di cui 38 dell'inventario (4.3 e 4.5) e 13 righe LIT. Esiti: 29 PASS, 17 UNSUPPORTED (2 bloccanti, 5 major, 10 minor), 4 MISMATCH (2 major, 2 minor), 1 UNVERIFIABLE.

## 1. Tabella di confronto

Lettura: T = testo completo letto in sessione (Zotero `fulltext()` o PMC); A = solo abstract (PubMed/OpenAlex); M = solo metadati (Crossref). Ogni DOI e' verificato su Crossref o PubMed in questa sessione.

| Studio (DOI) | Lettura | Domanda | Disegno della simulazione | Metodi | Criterio | Conclusioni | Rispetto al manoscritto |
|---|---|---|---|---|---|---|---|
| Mancl & DeRouen 2001 (10.1111/j.0006-341X.2001.00126.x) | T | correzione del sandwich GEE | 10-40 soggetti, 16-64 siti per bocca, binario | robusto, jackknife, MD | dimensione vicina al nominale | MD vicino al nominale con 10 soggetti | allineato; cluster molto piu' grandi (fino a 64 contro 10) |
| Kauermann & Carroll 2001 (10.1198/016214501753382309) | T | efficienza del sandwich | teoria e simulazione | sandwich, KC | copertura | sandwich con t troppo corto; nessuna sottocopertura con disegno bilanciato per individuo | allineato |
| Fay & Graubard 2001 (10.1111/j.0006-341X.2001.01198.x) | T | test di Wald con sandwich | K = 20 cluster, Poisson e binomiale; pazienti in 38 ambulatori | FG con b = 0,75 e gradi di liberta' propri | dimensione | FG corregge la dimensione | il manoscritto usa FG con t a N - p, non i gradi di liberta' di FG (dichiarato) |
| Pan & Wall 2002 (10.1002/sim.1142) | T | t/F per il sandwich | K = 10, 20, 30 soggetti, 20 osservazioni ciascuno, indipendenza | gradi di liberta' di Pan-Wall | dimensione | t/F con gradi di liberta' stimati | non valutato; contraddice "never a subject" in senso lato (LIT-002) |
| Lu et al. 2007 (10.1111/j.1541-0420.2007.00764.x) | T | MD contro KC | ampia gamma di cluster; trial di comunita' con disegno trasversale ripetuto | MD, KC | copertura | entrambi riducono la distorsione | allineato |
| Li & Redden 2015 (10.1002/sim.6344) | T | correzioni del sandwich nei CRT binari | 10, 20, 30 cluster; 10-150 osservazioni; ICC 0,001-0,05; CV 0-1; 3.000 repliche | z e t (K - 2) con MD, KC, FG; scambiabile; SAS | 4,2-5,8% (IC binomiale) | KC con t nominale da 10 cluster se CV < 0,6, altrimenti FG | KC qui anticonservativo: le differenze di disegno sono piu' d'una (LIT/I03481) |
| Li & Redden 2015b (10.1186/s12874-015-0026-x) | A | gradi di liberta' del GLMM nei CRT | CRT binari, pochi cluster | Residual, Containment, BW, Satterthwaite, KR | errore di I tipo | BW nominale da 10 cluster | GLMM qui solo con t N - p e LRT chi-quadrato (PoW-4) |
| Bie et al. 2021 (10.1002/sim.9126) | T | MMM contro GEE in piccoli campioni | K e N in 6-100; cluster trasversali; effetti casuali normali, gamma, varianza dipendente dal braccio; 120 punti di quadratura; 10.000 repliche | GEE con MD, MMM (MMLB) con correzioni ad hoc, LRT | rapporto SE, errore di I tipo | SE del MMM poco sottostimati dentro il cluster, molto fra cluster; correzioni ad hoc troppo conservative; preferiscono GEE con MD | allineato sul contrasto dentro il cluster; non citato per l'MMM (I03486); 20 punti qui (LIT-012); effetto fra cluster non studiato qui (LIT-003) |
| Hemming et al. 2025 (10.1016/j.jclinepi.2025.111838) | T | rassegna sistematica delle correzioni nei CRT paralleli | 20 studi, 4-200 cluster | GLMM, GEE, analisi per cluster | errore di I tipo | binario: GLMM con BW o MD; FG raccomandato per continui; oltre 40 cluster per garanzia; esclude permutazione e stepped-wedge | allineato su FG e MD; non usato in 4.3 (LIT-009) |
| Thompson et al. 2021, SMMR (10.1177/0962280220958735) | T | correzioni GEE negli stepped-wedge binari | 6-54 cluster, 3 o 6 sequenze, 24-300 osservazioni per cluster, 690 scenari, 1.000 repliche | FG, MD, KC (geesmv e approssimazione di Gallis), MBN, MW; gradi di liberta' C-P, CP-P, FG, PW | copertura 93,6-96,4% | FG e KC-approx con indipendenza quasi non distorti da 12 cluster; DF_FG nominale da 12, C-P da 18; conservativi sotto; KC di geesmv instabile; media MD/KC troppo grande | non citato (LIT-001); allineato su FG con t e su KC di geesmv (I03482); preferisce l'indipendenza (LIT-005) |
| Thompson et al. 2022 (10.1186/s12874-022-01699-2) | A | analisi di CRT binari con pochi cluster | 8-30 cluster | analisi per cluster, GLMM (AGQ, pseudo-verosimiglianza ristretta), GEE KC e FG | errore di I tipo | FG GEE controlla l'errore in >97% degli scenari; GEE la meno potente | allineato su FG (LIT-009) |
| Ford & Westgate 2017 (10.1002/bimj.201600182) | A | correzione migliorata nei CRT | CRT con pochi cluster | media MD/KC | errore di I tipo | la media e' nominale con piu' coerenza | qui la media e' anticonservativa (PoW-5) |
| Ford & Westgate 2018 [14] (10.1002/sim.7917) | A | correzioni in studi longitudinali piccoli | non letto | Westgate-Burchett; media MD/KC con N - p | errore di I tipo | i due metodi piu' coerenti | non valutati nel paper (PoW-5 lo fa sulle righe) |
| Ford & Westgate 2020 (10.1002/sim.8575) | A | GEE negli stepped-wedge binari | pochi cluster, molti periodi | combinazioni di correzioni e gradi di liberta' | errore di I tipo | GEE valida con correzione e df adatti | non citato (LIT-001) |
| Westgate & Burchett 2016 [13] (10.1002/sim.6967) | A | potenza in studi longitudinali piccoli | spaziatura uguale | sandwich modificato e df | potenza, validita' | aumenta i gradi di liberta' | non valutato |
| Gosho et al. 2023 [15] (10.1002/sim.9744) | A | GEE corrette per dati binari sparsi | 20-30 soggetti, prevalenza 5% | GEE, BCGEE, PGEE; molte correzioni del sandwich | copertura, convergenza | PGEE converge quasi sempre; Morel accettabile | non valutato; e' il caso del piccolo regime (PoW-3, LIT-007) |
| Mondol & Rahman 2019 (10.1002/sim.8126) | A | separazione nella GEE | dati longitudinali piccoli o sparsi | PGEE (Firth) con sandwich corretto | convergenza, errore di I tipo | PGEE consigliata per N <= 50 | non valutato (LIT-007) |
| Ishii et al. 2024 (10.1186/s12874-024-02368-2) | A | software | | pacchetto R geessbin: GEE, BCGEE, PGEE con 11 correzioni | | | rende PoW-6 eseguibile |
| Ouyang et al. 2026 (10.1002/sim.70782) | A | correzioni con effetti variabili nel tempo negli stepped-wedge | continui e binari | sandwich, KC, MD, MBN, FG | copertura | per esiti binari solo MBN sempre affidabile | Morel non valutato qui (I03505, PoW-7) |
| Austin 2024 (10.1002/sim.10260) | A | correzioni per covariate di cluster in studi osservazionali | GEE logistica | LZ, KC, FG, MD | SE e IC | KC sotto 15 cluster, FG 15-40, MD con t per gli IC | parzialmente allineato |
| Huang et al. 2022 (10.1080/19345747.2022.2100301) | A | CR2 per esiti binari con pochi cluster | CRT | CR2 con df di Bell-McCaffrey | SE, copertura | quasi non distorto da 10 cluster | qui CR2 con N - p (LIT-013) |
| McNeish 2016 [31] (10.1080/00273171.2016.1236237) | T | metodi di stima del GLMM logistico con pochi cluster | cluster trasversali, piu' effetti casuali, 2.500 repliche, SAS | AGQ (10 punti), Laplace, RPQL con e senza KR | distorsione, copertura | RPQL preferibile ad AGQ e Laplace con pochi cluster | caratterizzato come rassegna (I03487); RPQL non valutato |
| McNeish & Stapleton 2016 (10.1080/00273171.2016.1167008) | A | 12 metodi con pochissimi cluster | CRT educativo | GEE, multilivello, bayesiano, effetti fissi | errore di I tipo, potenza | GEE va male; effetti fissi bene | non citato |
| Qiu et al. 2024 (10.1002/sim.9950) | A | test GLMM nei CRT con 10-20 cluster | covariate | LRT, Wald, df | errore di I tipo | LRT con df between-within vicino al nominale | qui LRT chi-quadrato (PoW-4) |
| Capanu et al. 2013 (10.1002/sim.5866); Ju et al. 2020 (10.1186/s12874-020-01035-6) | A | Laplace contro AGQ | binari | LA, AGQ, PQL | distorsione | accurati con molte osservazioni per effetto casuale; AGQ non superiore con eventi rari | coerente con 3.10 |
| Schildcrout & Heagerty 2007 [22] (10.1111/j.1541-0420.2006.00680.x) | T (titolo, abstract) | MMM per serie lunghe | | modello marginalizzato con transizione | | | citato correttamente |
| Schildcrout & Heagerty 2005 (10.1093/biostatistics/kxi033); Pepe & Anderson 1994 (10.1080/03610919408813210) | A / M | GEE con covariate tempo-varianti | | correlazioni di lavoro | distorsione | non indipendenza puo' distorcere | requisito della AR(1) non detto (I03497) |
| Rokicki et al. 2018 (10.1097/MLR.0000000000000830) | A | inferenza con pochi gruppi | difference-in-differences longitudinale | GEE corretta, permutazione, wild cluster bootstrap, aggregazione | copertura | permutazione, wild bootstrap e GEE corretta vicini al nominale | permutazione e wild bootstrap non valutati |
| Cameron et al. 2008 (10.1162/rest.90.3.414); Canay et al. 2017 (10.3982/ECTA13081); Gail et al. 1996 (10.1002/(SICI)1097-0258(19960615)15:11<1069::AID-SIM220>3.0.CO;2-Q); Rabideau & Wang 2021 (10.1002/sim.9040) | M / A | wild bootstrap e test di randomizzazione con pochi cluster | | | | | alternative al test a inversione di segno non discusse |
| Hemerik et al. 2020 [27] (10.1111/rssb.12369) | M | test a inversione di segno | | | | | citato; contenuto verificato dal track Citazioni |
| Kahan et al. 2016 [6]; Barker et al. 2017 [7] | T | numero di cluster nei CRT | rassegne | | | mediana 25 (15-44) su 99; 45% degli stepped-wedge sotto 10 cluster | citati correttamente |
| Tong et al. 2025 (10.1177/17407745241276137) | A | stepped-wedge con 2-6 cluster | rassegna di 61 trial | | | 4 su 61 con correzione per pochi cluster; nessuno con permutazione | rafforza la motivazione dell'introduzione |

**Dove il manoscritto e' allineato.** Il sandwich non corretto e' anticonservativo. Fay-Graubard con riferimento t e' la correzione piu' stabile, conservativa con pochissimi cluster. Mancl-DeRouen e' quasi nominale in media ma anticonservativo in alcuni scenari. La Laplace e' imprecisa con pochi dati binari per cluster. Il modello marginalizzato ha SE quasi corretti per un contrasto dentro il cluster, e il riferimento t lo rende conservativo con 10 cluster (coerente con Bie 2021).

**Dove se ne discosta senza spiegarlo.**
- Kauermann-Carroll e' anticonservativo qui e nominale in Li & Redden: le cause possibili non sono discusse tutte.
- Viene scelta la AR(1), mentre Thompson 2021 preferisce l'indipendenza.
- Il modello misto e' valutato solo con stimatori di massima verosimiglianza, mentre la letteratura con pochi cluster raccomanda RPQL con Kenward-Roger o gradi di liberta' between-within.
- Il modello marginalizzato usa 20 punti di quadratura (Bie et al. ne usano 120).

## 2. Punti nascosti che un revisore chiederebbe

1. **Stima-obiettivo e generalizzabilita'.** Il contrasto e' un'interazione dentro il soggetto, con il trattamento assegnato al soggetto. I trial a cluster dell'introduzione stimano invece l'effetto di un trattamento assegnato al cluster. Per quel contrasto Bie 2021 trovano SE del modello marginalizzato molto sottostimati con pochi cluster (LIT-003).
2. **Stabilita' della regola.** La regola chiede zero flag su 2.376 punti, ciascuno con circa 1.000 repliche. In una replica il verdetto zero si ripeterebbe con probabilita' fra circa 0,4 e 0,9 per la strategia raccomandata, mentre la graduatoria non cambia (PoW-1, LIT-004). La regola e' stata scelta dopo aver visto i risultati (dichiarato), e la stabilita' non e' riportata.
3. **Scelta della correlazione di lavoro.** La scelta di AR(1) rispetto all'indipendenza e' decisa da 2 punti, uno nel quarto meccanismo che favorisce la AR(1) per costruzione. Con covariate tempo-varianti non esogene la AR(1) puo' distorcere la stima (Pepe & Anderson; Schildcrout & Heagerty 2005).
4. **Correzioni non valutate.**
   - La correzione di Morel-Bokossa-Neerchal: Ouyang 2026 la trovano l'unica affidabile per esiti binari con effetto variabile nel tempo.
   - La media MD/KC di Ford & Westgate (calcolata qui: PoW-5).
   - I gradi di liberta' di Fay-Graubard e di Pan-Wall.
   - CR2 con i gradi di liberta' di Bell-McCaffrey.
   - L'approssimazione KC di Gallis et al.
5. **Separazione nei piccoli disegni.** Il 54% dei dataset stimabili con T <= 3 e N <= 25 fallisce con stima divergente o informazione singolare, firme numeriche della separazione (PoW-3). La GEE penalizzata (Firth) e' il rimedio pubblicato ed e' disponibile in R (geessbin).
6. **Modello misto con gli stimatori raccomandati per pochi cluster.** Mancano RPQL con Kenward-Roger (McNeish 2016) e i gradi di liberta' between-within (Li & Redden 2015b; Hemming 2025; Qiu 2024). Il test del rapporto di verosimiglianza con riferimento F e' calcolato qui (PoW-4).
7. **Quadratura del modello marginalizzato.** Bie et al. usano 120 punti, qui 20. In M1c (varianza 3,5) il modello marginalizzato con riferimento normale rifiuta nel 6,8-7,0%.
8. **Alternative al test a inversione di segno.** Mancano la permutazione delle etichette di trattamento, il wild cluster bootstrap e gli intervalli per inversione di test di randomizzazione (Gail 1996; Cameron 2008; Canay 2017; Rabideau & Wang 2021; Rokicki 2018).
9. **Meccanismi.** Nessun meccanismo ha effetti casuali non normali o una varianza degli effetti casuali diversa fra bracci, i casi in cui Bie 2021 trovano il modello marginalizzato sensibile (LIT-008).
10. **Criterio di Bradley rispetto alla prassi del campo.** Li & Redden usano 4,2-5,8%, Thompson 2021 una copertura di 93,6-96,4%. Il paper usa 2,5-7,5% (convenzione registrata, non un errore). Nei confronti di 4.3 conviene dire che "nominal" in [16] e' un criterio piu' stretto.
11. **CR2 e Kauermann-Carroll.** CR2 e' lo stimatore di Bell-McCaffrey, che per la GEE coincide in teoria con KC. Sulle righe i due SE hanno rapporto mediano 1,000-1,004, eppure la Tabella 4 li tratta come metodi distinti (254 contro 199 punti). Va spiegato. Il confronto e' provvisorio: campione di 12 chunk.
12. **Potenza aggiustata.** Il valore critico per punto di disegno usa la distribuzione nulla vera, ignota in pratica (gia' dichiarato in 2.5). Nessun rilievo ulteriore dalla letteratura.

## 3. Analisi proposte (PoW) e loro esito

Base dei tempi: 0,77 secondi-core per replica alla cella piu' grande (generazione, modello marginalizzato, modello misto e GEE; note del server). Il solo percorso GEE vale circa 0,21 secondi-core: e' una stima derivata, 0,77 meno 0,32 del modello marginalizzato e 0,24 del modello misto, misurati alla stessa cella. Il tetto e' di 25 processi. Le celle piccole costano meno della piu' grande, quindi le stime sono tetti.

### Eseguite sulle righe esistenti (nessuna nuova simulazione)

**PoW-1. Stabilita' della regola (alta priorita'; rifatta nel giro 2 con i blocchi del protocollo).**

- **Cosa si ricampiona.**
  - Sono gli indici di replicato i = 1..1000 dentro ciascun modello di generazione (blocchi (model, i)).
  - Si usano gli indici salvati del progetto: `results_v4/S4/boot_idx.rds` per i modelli 1-5 e `results_v4/S6b/ns/boot_idx_67.rds` per 6-7, cioe' 2.000 righe di 1.000 indici per modello, generate con il seme 20261015.
  - M1c usa gli indici del modello 1. Nel quarto meccanismo lo scenario Mk usa quelli del modello k.
  - Dentro un ricampione b, lo stesso vettore di indici si applica a tutti i punti del modello, in ogni meccanismo, N e T (numeri casuali comuni), come in `R/s7d/s7d_agg.R`.
  - Per ogni punto si ricalcolano n (usabili), k (rifiuti) e nd (disegni stimabili), poi i flag con Clopper-Pearson e almeno 20 repliche usabili, poi la regola: ammessa la strategia senza flag; scelta quella con la quota usabile media piu' alta.
  - La v1 del giro 1 usava un solo vettore di pesi per tutti i modelli: le cifre cambiano al piu' di 0,5 punti percentuali.
- **Cosa rappresenta il 37,9%.** Approssima la probabilita' che una replica indipendente dello studio, con tassi veri uguali agli osservati in ogni punto, dia zero punti chiaramente sopra 7,5% per la GEE AR(1) con FG.
- **Perche' e' pessimista.** I tassi osservati vicini a 7,5% contengono gia' rumore: trattarli come veri e aggiungere la variabilita' del ricampionamento sovrastima la probabilita' di un flag. Lo stesso vale per il plug-in grezzo, che somma su punti indipendenti le probabilita' binomiali di flag con p = tasso osservato. Le due stime coincidono: 0,379 e 0,373, flag attesi 0,95 e 0,94.
- **Contrappeso ottimista.** Plug-in con tassi contratti: shrinkage empirical Bayes verso la media della strategia nel meccanismo, con varianza fra punti stimata come varianza osservata meno varianza binomiale. E' ottimista perche' riduce anche gli scarti veri fra punti.
- **Cancello superato.** n, k e nd per punto coincidono con le tabelle sorgente per tutte le 14 strategie (0 differenze); i flag coincidono con la Tabella 4.
- **Esito:**

| strategia | flag osservati | flag attesi: bootstrap | flag attesi: plug-in grezzo | flag attesi: plug-in contratto | P(zero flag): bootstrap | P(zero flag): plug-in grezzo | P(zero flag): plug-in contratto | usabili (%) |
|---|---|---|---|---|---|---|---|---|
| GEE AR(1), FG | 0 | 0.95 | 0.94 | 0.08 | 0.379 | 0.373 | 0.927 | 92.12 |
| Modello marginalizzato con controllo | 0 | 1.70 | 1.72 | 0.22 | 0.182 | 0.165 | 0.803 | 90.89 |
| GEE indipendenza, FG | 2 | 3.39 | 3.35 | 0.74 | 0.017 | 0.014 | 0.448 | 92.33 |
| GEE scambiabile, FG | 11 | 13.32 | 13.25 | 4.16 | 0.000 | 0.000 | 0.009 | 92.21 |
| GEE AR(1), MD | 16 | 24.70 | 24.76 | 9.80 | 0.000 | 0.000 | 0.000 | 92.01 |
| GEE AR(1), media MD/KC | 51 | 60.57 | 60.45 | 34.50 | 0.000 | 0.000 | 0.000 | 91.88 |
| GEE indipendenza, media MD/KC | 102 | 108.52 | 108.66 | 82.22 | 0.000 | 0.000 | 0.000 | 92.20 |
| GEE scambiabile, media MD/KC | 113 | 121.69 | 122.00 | 93.79 | 0.000 | 0.000 | 0.000 | 91.97 |
| GEE AR(1), KC | 132 | 135.77 | 135.90 | 104.39 | 0.000 | 0.000 | 0.000 | 91.90 |

- **Scelta della regola nei ricampioni:** nessuna strategia 50,6%, GEE AR(1) FG 36,9%, modello marginalizzato 10,9%, GEE indipendenza FG 1,65%. Entrambe le strategie della regola ammesse nel 7,15%. La quota usabile della GEE con indipendenza supera quella della AR(1) nel 100% dei ricampioni (scarto 0,21 punti, entro il margine di 2,5).
- **Interpretazione.**
  - La probabilita' che una replica dia di nuovo zero flag sta fra circa 0,4 e 0,9 per la GEE AR(1) e fra 0,2 e 0,8 per il modello marginalizzato: non si sa quanto dell'eccesso osservato vicino a 7,5% sia rumore.
  - La graduatoria per flag attesi e' identica nei tre calcoli: GEE AR(1) FG, poi modello marginalizzato, poi GEE con indipendenza FG, poi GEE scambiabile FG, poi tutte le altre con almeno 9,8 flag attesi anche nel calcolo ottimista.
  - Quindi la raccomandazione regge. La frase "never clearly above" descrive la realizzazione e conviene affiancarla da una frase sulla replica (LIT-004).
- **Frase proposta per 4.5 (+63 parole):** "The zero count is not a stable property of these strategies: treating the observed rates as true, which is pessimistic near 7.5%, a replication would leave the GEE strategy and the marginalized model without a flagged design point with probability about 0.4 and 0.2 (0.9 and 0.8 after shrinking the rates), yet in every calculation they had the fewest expected flags (Section S11)."
- **Versione minima (+24):** "The zero count of flagged design points is uncertain in a replication; the ranking of the strategies by expected flags is not (Section S11)." Tabella e metodo in S11.
- **Riferimento:** Morris, White & Crowther 2019 (10.1002/sim.8086).

**PoW-2. CR2 contro Kauermann-Carroll (media priorita', campione).**
- Protocollo: 12 chunk di GLMM-native M1 (N = 10, 25, 65; T = 4, 10; repliche 1-200).
- Esito: rapporto mediano se_cr2/se_kc tra 1,000 e 1,004. Percentili 1-99% tra 0,85 e 1,20 (con un estremo a 39 a N = 10, T = 4, scambiabile); a N = 65 entro +/-0,2%.
- Interpretazione: CR2 e' di fatto una seconda implementazione di KC, vicina ma non identica con 10 cluster. Riguarda I03507 e LIT-013. Provvisorio (campione).

**PoW-3. Perche' la GEE AR(1) con FG non risponde nei piccoli disegni (alta priorita'; definizioni precisate nel giro 2).**

- **Regione:** T <= 3 e N <= 25, i 216 punti dei quattro meccanismi (righe S7d e S7e).
- **Denominatore:** dataset generati con disegno stimabile (`design_degenerate` FALSE), 215.012 per braccio. Il denominatore della Tabella 4 include anche i pochi dataset non generati, da cui 45,9% contro 46,0%.
- **Unita':** il dataset. Per questa strategia c'e' un fit per dataset (geeglm con AR(1) e SE di FG da fewclust), quindi quota di dataset e quota di fit coincidono.
- **Classificazione, in quest'ordine:**
  1. errore di geeglm (`gee_ok` FALSE);
  2. stima dell'interazione non finita o |stima| >= 10 (colonna `est`);
  3. SE di FG non calcolabile con il messaggio di fewclust "model-based information matrix is singular (aliased or degenerate design, possibly separation)" (colonna `msg_fg`);
  4. altri fallimenti dell'SE;
  5. usabile.
- **Esito (H0):** usabile 98.992 (46,0%); stima non finita o |stima| >= 10: 85.797 (39,9%); informazione singolare: 29.681 (13,8%); altri fallimenti dell'SE: 542 (0,3%); errori di geeglm: 0. Le classi 2 e 3 insieme fanno il 53,7% (H1: 52,6%). Per T: con T = 2 usabile 30,1%, con T = 3 61,9%.
- **Cosa non e'.** Non e' un controllo diretto della separazione: le righe non contengono gli esiti per braccio e occasione. Sono firme numeriche compatibili con separazione o quasi-separazione. Nel testo va scritto "con stime divergenti o informazione singolare, firme della separazione", non "per separazione".
- **Interpretazione:** se la GEE penalizzata risponde e mantiene il livello, la frase delle conclusioni cade; se fallisce anche lei, la frase regge, ma limitata ai metodi valutati (LIT-007).
- **Riferimenti (solo abstract):** Mondol & Rahman 2019 (10.1002/sim.8126); Gosho et al. 2023 (10.1002/sim.9744).

**PoW-4. Test del rapporto di verosimiglianza con riferimento F(1, N - p) (media priorita').**
- Protocollo: righe `variant_cc/S4/rows_v4.rds`, `S6b/ns/rows_ns.rds` e S7e; statistica ricavata dal p-value chi-quadrato e riferita a F(1, N - p).
- Cancello superato: con il chi-quadrato i flag sono 223 + 0 + 189 = 412, identici alla Tabella 4.
- Esito:
  - Con F(1, N - p): 293 punti chiaramente sopra 7,5% (124 sulla griglia principale e M1c, 0 in M6-M7, 169 sul quarto meccanismo).
  - Tasso medio 5.9% contro 6.9%.
  - Con N <= 25 e T >= 4 sulla griglia principale i flag scendono da 49 a 0. L'eccesso residuo sta a T = 2 con N > 25 (74 punti) e sul quarto meccanismo.
- Interpretazione: l'argomento di 4.1 contro il criterio precedente regge anche con il riferimento F (LIT-006). L'eccesso dell'LRT non e' un problema di riferimento con pochi cluster.
- Riferimenti: Li & Redden 2015b (10.1186/s12874-015-0026-x); Qiu et al. 2024 (10.1002/sim.9950).

**PoW-5. Media di Mancl-DeRouen e Kauermann-Carroll (controprova della previsione PS5, che appartiene al track Dubbi).**
- Protocollo: varianza media (se_md^2 + se_kc^2)/2, correzioni fewclust, t a N - p, stessi flag della Tabella 4.
- Esito: chiaramente sopra 7,5% in 51 punti (AR(1)), 102 (indipendenza) e 113 (scambiabile). Quote usabili 91,9-92,2%. Previsione PS5 ("almeno 1 punto") confermata.
- Interpretazione: qui la media e' anticonservativa perche' KC lo e'. Il risultato e' coerente con Thompson 2021, che la trovano inadatta negli stepped-wedge, e contrario a Ford & Westgate 2017/2018. Se ne puo' dare il risultato in 3.5 invece di dichiararla non valutata.
- Riferimenti: Ford & Westgate 2017 (10.1002/bimj.201600182), 2018 (10.1002/sim.7917).

### Non eseguite: richiedono di rigenerare i dataset (decisione PI, D-64)

| PoW | Punto nascosto | Protocollo | Costo stimato | Se va in un verso | Se va nell'altro | Riferimento | Priorita' |
|---|---|---|---|---|---|---|---|
| 6. GEE penalizzata nei piccoli disegni | separazione (PoW-3) | 216 punti con T <= 3 e N <= 25, entrambi i bracci, 1.000 dataset; PGEE (geessbin) con Morel e FG, t a N - p, stessa regola | 432.000 repliche: 25 ore-core (percorso GEE), tetto 92; 1,0-3,7 h con 25 processi | risponde e mantiene il livello: cadono LIT-007 e la frase "no strategy ... available and calibrated" di 4.1 | non risponde o non tiene il livello: la frase regge, ma limitata ai metodi valutati | 10.1002/sim.8126; 10.1002/sim.9744; 10.1186/s12874-024-02368-2 | alta |
| 7. Correzioni non valutate | Morel, gradi di liberta' FG e Pan-Wall, KC-approx, CR2 con Bell-McCaffrey | 2.376 punti sotto H0, una sola passata GEE | 139 ore-core (percorso GEE), tetto 508; 5,5-20 h; entrambi i bracci 11 h | una di esse senza flag con quota usabile maggiore: cambia la strategia raccomandata (abstract, 4.1, Figura 1) | tutte con flag: la raccomandazione si rafforza | 10.1002/sim.70782; 10.1177/0962280220958735; 10.1177/1536867X20931001 | alta |
| 8. Punti di quadratura del modello marginalizzato | 20 contro 120 punti | M1c e M6-M7 (456 punti) sotto H0, 50 o 120 punti | 128-270 ore-core; 5-11 h (ipotesi: costo lineare nei punti) | il tasso di M1c scende verso 5%: il modello marginalizzato migliora | nessun cambiamento: il limite si chiude | 10.1002/sim.9126 | media |
| 9. Modello misto con stimatori per pochi cluster | RPQL con Kenward-Roger; AGQ con pendenze casuali (GLMMadaptive); gradi di liberta' between-within | sottoinsieme M1-M5 con N <= 25 | da misurare (pbkrtest non supporta glmer binomiale su questo host); AGQ con pendenze casuali circa 1-2 s per fit | il modello misto diventa ammesso: cambia la frase delle conclusioni sul modello misto | resta escluso: la frase si rafforza | 10.1080/00273171.2016.1236237; 10.1186/s12874-015-0026-x | media |
| 10. Permutazione e wild bootstrap | alternative al test a inversione di segno | piccoli disegni sotto H0, 1.999 permutazioni | 252 ore-core, circa 10 h (ipotesi 2 ms per fit); versione sullo score molto meno | piu' potente del test a inversione di segno: va citata come complemento | uguale: la scelta e' confermata | 10.1162/rest.90.3.414; 10.3982/ECTA13081; 10.1002/sim.9040 | bassa |
| 11. Nuovo meccanismo o nuova stima-obiettivo | effetto del trattamento assegnato al cluster; effetti casuali non normali o con varianza diversa fra bracci | griglia principale completa, entrambi i bracci | 616 ore-core, circa 25 h (tetto) | il modello marginalizzato fallisce fra cluster: la restrizione di LIT-003 diventa un risultato | regge: la raccomandazione si generalizza | 10.1002/sim.9126 | media (alternativa a costo zero: restringere il testo, LIT-003) |

## 4. Frasi della discussione (e di sintesi) da riformulare

Le versioni ponderate sono nella colonna `note` del log; qui solo l'elenco, per gravita'.

| iid | Dove | Perche' | Esito |
|---|---|---|---|
| LIT-003 | Abstract, "We recommend the GEE strategy ... comparable alternative" | la stima-obiettivo (interazione dentro il soggetto) non e' detta; Bie 2021 mostra un comportamento diverso fra cluster | UNSUPPORTED, bloccante |
| LIT-007 | Conclusioni, "no analysis estimates the effect reliably" | la separazione spiega il 54% dei fallimenti e la GEE penalizzata non e' valutata | UNSUPPORTED, bloccante |
| LIT-004 | 4.1/4.5, "neither was ever clearly above 7.5%" | in una replica il verdetto zero e' incerto (probabilita' 0,4-0,9 e 0,2-0,8); la graduatoria e' stabile | UNSUPPORTED, minor |
| LIT-005 | 4.1, AR(1) | l'indipendenza ha 2 flag e la stessa quota usabile; la letteratura la preferisce | UNSUPPORTED, minor |
| I03481 | 4.3, differenza con Li & Redden | spiegazione incompleta (ICC, dimensione, stima-obiettivo, gradi di liberta') | UNSUPPORTED, major |
| I03486 | 4.3, modello marginalizzato "moderate to long series" | ignora Bie 2021 | UNSUPPORTED, major |
| I03487, I00656 | 4.3, McNeish "reviewed" | studio di simulazione con conclusione a favore di RPQL | MISMATCH, major |
| LIT-009 | 4.3, sezione intera | manca l'accordo con Hemming 2025 e Thompson 2021/2022 | UNSUPPORTED, major |
| I03505 | 4.5, metodi non valutati | mancano gli stimatori del modello misto per pochi cluster (McNeish, Hemming) e la correzione di Morel; la parte su D-64 e' ritirata | UNSUPPORTED, major |
| LIT-001 | Introduzione, letteratura sugli stepped-wedge | Thompson 2021 (SMMR) e Ford & Westgate 2020 assenti | UNSUPPORTED, major |
| I03482, I03497, I03504, I03507, LIT-002, LIT-008, LIT-012, LIT-013 | vari | completamenti | UNSUPPORTED, minor |
| I00671, LIT-010 | [31] in 4.5; voce [18] | citazione che non sostiene il limite di lme4; ultimo autore di [18] sbagliato (McGee G) | MISMATCH, minor |

## 5. Riproducibilita'

- Script (artifact e server): `~/MMM_rev_S8/Letteratura/lit_pow_rule_stability.R` (PoW-1, PoW-3, PoW-5) e `~/MMM_rev_S8/Letteratura/lit_pow_lrt_f.R` (PoW-4). Esecuzione: `cd ~/MMM_rev_S8/Letteratura && OMP_NUM_THREADS=1 Rscript --vanilla <script>`, 2 processi R (`mclapply`, `mc.cores = 2`). Durate: 1,4 minuti per la lettura delle righe; il job PoW-1 e' durato circa 10 minuti, il job PoW-4 circa 5.
- Output in `~/MMM_rev_S8/Letteratura/`: `lit_pow1_gate.csv`, `lit_pow1_pass.csv`, `lit_pow1_choice.csv`, `lit_pow1_joint.csv`, `lit_pow1_flagged_points.csv`, `lit_pow3_reasons.csv`, `lit_pow3_by_T.csv`, `lit_pow4_lrt_F.csv`, `lit_pow4_lrt_F_byNT.csv`, `lit_pow_sessionInfo.txt`.
- Temperatura alla partenza: Tctl 86,6 C (PoW-1) e 70,2 C (PoW-4), sotto il limite di 90 C.
- Letteratura: Zotero `build_index()` della sessione (7.459 voci), poi `fts()` e `fulltext()`; PubMed E-utilities; OpenAlex con chiave; Crossref per i DOI.

## 6. Deviazioni

- Righe assegnate: il task non elenca le righe dell'inventario del track. Le ho prese come le 38 righe di tipo `discussion_sentence` e `citation` delle sezioni 4.3 e 4.5 (punto 4 del task). Le altre affermazioni contestate sono righe LIT.
- Fonti lette sul solo abstract (A nella tabella): Ford & Westgate 2017/2018/2020, Westgate & Burchett 2016, Gosho 2023, Mondol & Rahman 2019, Ouyang 2026, Thompson 2022, Qiu 2024, Huang 2022, Austin 2024, Rokicki 2018, Tong 2025, Rabideau & Wang 2021, Capanu 2013, Ju 2020, Li & Redden 2015b, McNeish & Stapleton 2016, Schildcrout & Heagerty 2005. Solo metadati (M): Pepe & Anderson 1994, Cameron 2008, Canay 2017, Gail 1996, Hemerik 2020, Pustejovsky & Tipton 2018. Il testo completo di Ford & Westgate 2018 e Gosho 2023 non e' accessibile (closed / 403).
- PoW-2 su un campione di 12 chunk di un solo scenario: risultato provvisorio.
- Gallis et al. 2020 non letta: e' citata solo attraverso Thompson 2021.
- PoW-1 del giro 1 (pesi comuni a tutti i modelli, seme proprio) e' sostituita dalla v2 con gli indici (model, i) del progetto; i file v1 restano per confronto.
- I numeri del registro dentro le frasi di 4.3 e 4.5 non sono ricalcolati in questo track, salvo quelli che la PoW-1 riottiene dalle righe (132, 199, 57, 2.376).

## 7. Proposta di abstract (giro 2)

La proposta soddisfa insieme due richieste. La prima e' nominare la stima-obiettivo (LIT-003). La seconda e' la riserva sul modello misto del track Discussione: "in random-intercept models" nella frase 6, togliendo il doppio "and a t reference" nella frase 4. Conteggio con il metodo del progetto: parole separate da spazi dopo aver tolto i commenti HTML, ogni numero in linea vale una parola. Il conteggio e' fatto sul render, e lo stesso codice sul testo qmd (ogni `r n()` = una parola) da' lo stesso valore.

- **Abstract attuale: 199 parole. Proposta: 200 parole.**
- **Modifiche:**
  - frase 2: "Small-sample corrections to the sandwich variance" diventa "Small-sample sandwich corrections" (-3);
  - frase 3: "We compared analysis strategies over 2,376 simulated design points" diventa "We compared strategies for a within-subject treatment-by-time interaction over 2,376 design points" (+3);
  - frase 4: "...the Fay-Graubard correction and a t reference, and the marginalized multilevel model with a convergence check and a t reference, were..." diventa "...and the Fay-Graubard correction, and the marginalized multilevel model with a convergence check, both with a t reference, were..." (-2);
  - frase 6: "with adaptive quadrature in random-intercept models its excess..." (+3);
  - frasi 1, 5, 7 e 8 invariate.

Testo (render):

> Binary outcomes measured repeatedly on few subjects are analysed with generalized estimating equations (GEE), mixed models or marginalized multilevel models. Small-sample sandwich corrections have been evaluated mainly by number of clusters, rarely by number of occasions. We compared strategies for a within-subject treatment-by-time interaction over 2,376 design points (10-145 clusters, 2-10 occasions, four mechanisms, one purely serial) and in subsamples of two clinical datasets. GEE with a first-order autoregressive working correlation and the Fay-Graubard correction, and the marginalized multilevel model with a convergence check, both with a t reference, were never clearly above 7.5% type I error at any design point and returned a usable result in 92% and 91% of datasets. Their calibration differed materially only at ten clusters, where GEE was closer to nominal. The mixed-model Wald test was anticonservative under the Laplace approximation; with adaptive quadrature in random-intercept models its excess remained only with ten occasions. With two or three occasions and up to 25 clusters both strategies answered in fewer than half of datasets; a sign-flip score test answered in 72% and kept its level, but gives no estimate. We recommend the GEE strategy, which needs no convergence check; the marginalized model is a comparable alternative.

Testo (qmd, con i numeri del registro):

```
Binary outcomes measured repeatedly on few subjects are analysed with generalized estimating equations (GEE), mixed models or marginalized multilevel models. Small-sample sandwich corrections have been evaluated mainly by number of clusters, rarely by number of occasions. We compared strategies for a within-subject treatment-by-time interaction over `r n("f:x|ndp_all5")` design points (10-145 clusters, 2-10 occasions, four mechanisms, one purely serial) and in subsamples of two clinical datasets. GEE with a first-order autoregressive working correlation and the Fay-Graubard correction, and the marginalized multilevel model with a convergence check, both with a t reference, were never clearly above 7.5% type I error at any design point and returned a usable result in `r n("f:t4|gee_ar1_fg|pu")`% and `r n("f:t4|mmm_cc|pu")`% of datasets.<!-- claim: reg:t4 as=descr --> Their calibration differed materially only at ten clusters, where GEE was closer to nominal.<!-- claim: reg:hb as=descr --> The mixed-model Wald test was anticonservative under the Laplace approximation; with adaptive quadrature in random-intercept models its excess remained only with ten occasions.<!-- noclaim: summary of Sections 3.6 and 3.10 --> With two or three occasions and up to 25 clusters both strategies answered in fewer than half of datasets; a sign-flip score test answered in `r n("f:t4|flip|pu_small")`% and kept its level, but gives no estimate.<!-- claim: reg:t4 as=descr --> We recommend the GEE strategy, which needs no convergence check; the marginalized model is a comparable alternative.<!-- noclaim: recommendation -->
```

Nelle conclusioni basta "In the designs studied, for a treatment-by-time interaction, we recommend ..." (+4 parole).

## 8. Livello di lettura delle fonti dei rilievi major e bloccanti (giro 2)

| Fonte | Lettura | Frase usata (parafrasi) e posizione | Rilievi che ne dipendono | Esito dopo il controllo |
|---|---|---|---|---|
| Bie et al. 2021 [18] | testo completo (Zotero itemID 3965) | con pochi cluster gli SE del MMM sono poco sottostimati per le associazioni dentro il cluster e molto per quelle fra cluster (abstract p. 5298; risultati p. 5303); X2 e' una covariata binaria di cluster, l'indicatore di braccio (sezione 3.1, p. 5301); 120 punti di quadratura, pochi punti peggiorano le stime con sigma grande (sezione 3.2, p. 5302); le correzioni ad hoc dei gradi di liberta' per il MMM sono troppo conservative e gli autori preferiscono GEE con MD (sezione 4, p. 5310) | LIT-003 (bloccante), I03486 (major), LIT-008, LIT-012 | confermati |
| Bie et al. 2021, autori | Crossref, giro 2 | api.crossref.org/works/10.1002/sim.9126: Ruofan Bie (first), Sebastien Haneuse, Nathan Huey, Jonathan Schildcrout, Glen McGee; Statistics in Medicine 40(24):5298-5312, 2021. Il frontespizio del PDF in Zotero da' gli stessi cinque autori | LIT-010 | confermato: l'ultimo autore e' McGee G, non Shepherd BE |
| McNeish 2016 [31] | testo completo (Zotero itemID 17950) | con pochi cluster la distorsione nota di RPQL e' trascurabile e preferibile a quella di AGQ e Laplace (Conclusions, p. 13 del PDF, paginazione 1-15 della versione online); la Laplace soffre con pochi cluster o cluster piccoli (p. 6); AGQ e' disponibile in lme4 (p. 5) e ha difficolta' con effetti casuali correlati (p. 12) | I03487 e I00656 (major), I03505 (major), LIT-011, I03504 | confermati |
| Thompson et al. 2021, SMMR | testo completo (PMC) | KC di geesmv instabile, ripetuto con l'approssimazione di Gallis et al. (sezione 4.1.6); copertura bassa sotto 93,6%, alta sopra 96,4% (sezione 4.2); FG e KC buoni con indipendenza, DF_FG nominale da 12 cluster, C-P da 18 (sezione 7); la media KC/MD darebbe SE troppo grandi (sezione 7); MD, MW e MBN sovracorretti (sezione 7) | LIT-001 (major), LIT-009 (major), LIT-005, I03482, I03505 | confermati |
| Hemming et al. 2025 [17] | testo completo (Zotero itemID 17931) | con esiti binari, GLMM con between-within o GEE con MD; FG raccomandato con esiti continui (abstract e "What is new", p. 1); molte correzioni sviluppate e valutate prima su dati longitudinali (sezione 1.2, p. 2); esclusi permutazione, bootstrap e stepped-wedge (sezione 3.2, p. 3) | LIT-009 (major), I03505 (major), LIT-001, LIT-002 | confermati |
| Ford & Westgate 2018 [14] | solo abstract (PubMed 30073684; testo chiuso, 403) | la media MD/KC con gradi di liberta' soggetti meno parametri e il metodo di Westgate-Burchett raggiungono il livello nominale piu' spesso | PoW-5 (controprova, non un rilievo), parte di I03505 | nessun rilievo major poggia solo su questa fonte; in I03505 e' marcata come abstract |
| Gallis et al. 2020 (xtgeebcv) | NON letta (solo metadati Crossref) | citata solo attraverso Thompson 2021, sezione 4.1.6 | I03482, I03507 (entrambi minor) | il riferimento a Gallis passa per Thompson; nessun rilievo major |
| Li & Redden 2015 [16] | testo completo (Zotero itemID 17897) | KC con t nominale da 10 cluster (abstract, p. 281); 10-150 osservazioni per cluster, ICC 0,001-0,05, SAS 9.3 (sezione 2.1, p. 283); K - 2 gradi di liberta', criterio 4,2-5,8% (sezione 2.2, p. 284); KC se CV < 0,6, altrimenti FG (p. 287) | I03481 (major), I03479 | confermati |
| Mondol & Rahman 2019; Gosho et al. 2023 [15] | solo abstract | la GEE penalizzata alla Firth converge e da' stime finite con separazione | LIT-007 (bloccante) | il bloccante poggia sul testo (4.5 dichiara la GEE penalizzata non valutata) e sulle righe (PoW-3); la letteratura e' di supporto, al livello abstract |

## 9. Versioni minime delle aggiunte e cosa va nel supplemento (giro 2)

Ogni voce sostituisce una frase o una sua parte che esiste una volta sola nel render (ancora verificata). Il saldo e' in parole.

| iid | sezione | parole prima | parole dopo | saldo | versione minima | nel supplemento |
|---|---|---|---|---|---|---|
| I03481 | 4.3 | 29 | 30 | +1 | Their trials had low intracluster correlation, large clusters, a cluster-level treatment and K - 2 degrees of freedom [16], so the difference may reflect the design rather than the estimator. | elenco completo delle differenze (ICC 0,001-0,05, 10-150 per cluster, CV < 0,6, criterio 4,2-5,8%, SAS) in S4 |
| I03486 | 4.3 | 21 | 28 | +7 | With few clusters, Bie et al. found the standard error of the marginalized model nearly correct for within-cluster contrasts, as here, but too small for between-cluster ones [18]. | il resto di Bie (120 punti di quadratura, correzioni ad hoc conservative, effetti casuali non normali) in S2 o S13 |
| I03487 | 4.3 | 14 | 17 | +3 | With few clusters McNeish preferred penalized quasi-likelihood with Kenward-Roger to quadrature [31]; it was not evaluated here. | nessuno |
| LIT-009 | 4.3 | 12 | 31 | +19 | We reproduce the anticonservatism of the uncorrected sandwich with few clusters [8,9] and, as in stepped-wedge trials with binary outcomes [T21], find Fay-Graubard with a t reference the most stable correction. | confronto per studio (Thompson 2021 e 2022, Li & Redden, Austin 2024) in una tabella del supplemento |
| LIT-001 | 1 | 18 | 23 | +5 | They have also been evaluated on parallel cluster-randomized trials, where the simulation evidence has been reviewed systematically [16,17], and on stepped-wedge trials [T21,FW20]. | nessuno (citazioni) |
| LIT-002 | 1 | 18 | 11 | -7 | In the evaluations cited here, the unit within a cluster is | nessuno (accorcia) |
| LIT-005 | 4.1 | 16 | 29 | +13 | The working correlation changed the level of Fay-Graubard by less than half a point on average, and with independence it was clearly above 7.5% at only two design points, | i due punti e i flag attesi in S10 o S11 |
| LIT-007 | 5 | 18 | 20 | +2 | and none of the analyses evaluated estimates the effect reliably; in that regime the design limits what they can show. | quota per causa (PoW-3) in S11 |
| I03505 | 4.5 | 17 | 28 | +11 | Penalized GEE [15], the Morel correction [12], the degrees of freedom proposed with the Fay-Graubard correction [10] and penalized quasi-likelihood for the mixed model [31] were not evaluated. | elenco completo (CR2 con Bell-McCaffrey, permutazione, wild bootstrap, punti di quadratura del modello marginalizzato) e risultato della media MD/KC in S4 |
| I03504 | 4.5 | 8 | 8 | +0 | could not be examined with random slopes [26]. | nessuno |
| I03482 | 4.3 | 17 | 24 | +7 | geesmv computes a different quantity, as Thompson et al. also found [T21], so published values obtained with it are not directly comparable with ours. | l'intera frase su fewclust/geesmv puo' andare in S4 (risparmio 30 parole) |
| LIT-003c | 5 | 12 | 16 | +4 | In the designs studied, for a treatment-by-time interaction, we recommend GEE with an AR(1) working correlation, | nessuno |

- **Saldo di tutte le versioni minime: +65 parole.**
- Se la frase su fewclust e geesmv (I03482, 32 parole) va intera in S4 invece di essere estesa, il saldo scende a +26.
- La frase su LIT-004 aggiunge +24 nella versione minima, oppure +63 nella versione completa.
- Restano solo nel supplemento: I03497 (Pepe-Anderson), I03507 (CR2 come Bell-McCaffrey), LIT-008, LIT-012, LIT-013.

## 10. Diario dei cambiamenti del giro 2

1. **Ritirato:** "in 4.5 mancano tre analisi che D-64 chiede di dichiarare" (I03505, LIT-012). Il corpo di D-64 chiede solo Laplace, quadratura adattiva, GEE penalizzata e gradi di liberta' di FG, e 4.5 li dichiara tutti. L'errore veniva dalla riga dell'indice delle decisioni, che elenca le sei analisi proposte. I03505 passa da MISMATCH a UNSUPPORTED major per la sola parte di letteratura.
2. **Ridimensionati a minor:** LIT-004 e LIT-005. Con i blocchi (model, i) del protocollo le cifre del bootstrap non cambiano (37,9% contro 37,5%), ma il plug-in ottimista mostra che il verdetto zero puo' replicarsi con probabilita' alta.
3. **Precisato:** LIT-007 si basa sul testo e sulle righe, non sulla letteratura letta in abstract; "per separazione" diventa "con le firme numeriche della separazione".
4. **Ricontrollato su Crossref in giro 2:** l'ultimo autore di [18] e' Glen McGee.
