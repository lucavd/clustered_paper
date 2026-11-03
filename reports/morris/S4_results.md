### Tabella E. Attese pre-registrate e osservato

| id | attesa | tolleranza | osservato | superato |
|---|---|---|---|---|
| Q0 | md5 di results_full/ (28.800), results_v4 S2 (126) e S3 (18.720) e di rows_timetrt.rds invariati prima e dopo | zero falliti | prima: results_full: 28800 checked, 0 FAIL lines; results_v4 S2: 126 checked, 0 FAIL lines; results_v4 S3: 18720 checked, 0 FAIL lines; rows_timetrt.rds 092a5ad60b60e0ea28cb6b78ba2adbeb (expected 092a5ad60b60e0ea28cb6b78ba2adbeb) \|\| dopo: results_full: 28800 checked, 0 FAIL lines; results_v4 S2: 126 checked, 0 FAIL lines; results_v4 S3: 18720 checked, 0 FAIL lines; rows_timetrt.rds 092a5ad60b60e0ea28cb6b78ba2adbeb (expected 092a5ad60b60e0ea28cb6b78ba2adbeb) | si |
| Q1a | nessuna riga glmer_lrt v3 in M3-M5; re_structure come D-08 in M3-M5 e (1\|id) in M1c | esatta | righe v3: 0; re_structure inattese: 0 | si |
| Q1b | celle F1 e F2 interamente da S3 | esatta | righe non da S3: 0 | si |
| Q1c | 3.024 celle, 1.512 design point, ogni cella con righe per tutti gli 8 metodi primari | esatta | celle 3.024, design point 1.512; celle con meno di 8 metodi: 342 (qualunque riga), 342 (righe time:trt); duplicati 0 | no |
| Q1d | tabella lunga 3.024 x 9 x 3 x 2 = 163.296 righe | esatta | 163.296 | si |
| Q2 | n_den: 2.878.518 (griglia) e 143.742 (M1c); zero righe U0 su chiavi degeneri; ogni chiave non degenere con una riga di analisi | esatta | 2.878.518; 143.742; U0 su degeneri 0; chiavi senza righe 0 | si |
| Q3 | righe non toccate da S3 riproducono O_level_all (n, k_rej) e O_coverage_all (n, cov_gen) | 100%, 1e-12 | O_level_all (H0): 12.395/12.395 n uguali, diff. max 5.3e-16, mancanti 0, inverse mancanti 0; O_coverage_all (H1): 13.144/13.144 n uguali, diff. max 1.1e-15, mancanti 0 | si |
| Q4 | n_U0, n_U1, n_U2 = V4_S2_usable_by_cell nelle combinazioni non toccate | esatta | 17.725 righe; U0 17.725, U1 17.725, U2 17.725 uguali; mancanti 0 | si |
| Q5a | rsimsum 0.13.1 su 300 combinazioni: 17 misure e MCSE entro 1e-10 | 1e-10 | 5.070 confronti, oltre tolleranza 0, diff. max 4.4e-13 | si |
| Q5b | MSE = bias^2 + SE empirico^2 (n - 1)/n in tutte le combinazioni con n_obs >= 2 | 1e-12 | marginal: 0/100.893 oltre, diff. max 2.9e-14; conditional per replicate: 10.052/10.580 oltre, diff. max 4.5e-03; conditional: 0/10.579 oltre, diff. max 2.1e-13 | no |
| Q6 | verita' AR1 per replicato (riproduzione S2, biiezione con i degeneri, replay = generatore); quadratura M1 | 1e-12; < 1e-8 | S2 reproduction M1-M4: ok; failures = degenerate keys, main: ok; failures = degenerate keys, M1c: ok; replay = generator, 200 keys: ok; quadrature reproduces S2 M1: ok | si |
| Q7 | SE bootstrap della potenza grezza entro il 10% dalla MCSE binomiale | >= 95% | 58.198/58.198 (1,0000); rapporto mediano 1,000 | si |
| Q8 | D-08: ogni cella del blocco L con almeno un replicato usabile in entrambe le versioni | 100% | 1.021/1.248 | no |
| Q9 | glmer_bobyqa: p = 2 Phi(-abs(z)); n_obs <= n_obs(glmer_lrt) per cella | esatta | p diversi 0 su 2.634.153 righe; celle 3.024/3.024 (eccesso massimo 0) | si |

### E1. Q1c: celle senza righe time:trt per un metodo primario

| metodo | scenario | T | coppie cella x metodo | con righe | solo righe di fallimento | glmer mai convergente nella cella | motivo piu' frequente |
|---|---|---|---|---|---|---|---|
| glmer | M3 | 2 | 84 | 84 | 84 | 84 | number of observations (=34) < number of random effects (=40) for term |
| glmer | M4 | 2 | 84 | 84 | 84 | 84 | number of observations (=59) < number of random effects (=70) for term |
| glmer | M5 | 2 | 96 | 96 | 96 | 96 | number of observations (=34) < number of random effects (=60) for term |
| glmer | M5 | 3 | 78 | 78 | 78 | 78 | number of observations (=65) < number of random effects (=75) for term |
| glmer_lrt | M3 | 2 | 84 | 84 | 84 | 84 | fit returned NULL |
| glmer_lrt | M4 | 2 | 84 | 84 | 84 | 84 | fit returned NULL |
| glmer_lrt | M5 | 2 | 96 | 96 | 96 | 96 | fit returned NULL |
| glmer_lrt | M5 | 3 | 78 | 78 | 78 | 78 | fit returned NULL |
| shrink | M3 | 2 | 84 | 0 | 0 | 84 |  |
| shrink | M4 | 2 | 84 | 0 | 0 | 84 |  |
| shrink | M5 | 2 | 96 | 0 | 0 | 96 |  |
| shrink | M5 | 3 | 78 | 0 | 0 | 78 |  |

### E2. Quota usabile per definizione (griglia principale, entrambi gli arm, denominatore D-11)

| method | U0 | U1 | U2 |
|---|---|---|---|
| gee | 0,942 | 0,938 | 0,942 |
| gee_fg | 0,922 | 0,922 | 0,922 |
| gee_kc | 0,922 | 0,922 | 0,922 |
| gee_md | 0,922 | 0,922 | 0,922 |
| glmer | 0,811 | 0,635 | 0,635 |
| glmer_bobyqa | 0,813 | 0,769 | 0,769 |
| glmer_lrt | 0,865 | 0,798 | 0,798 |
| mmm | 0,907 | 0,783 | 0,907 |
| shrink | 0,821 | 0,640 | 0,640 |

### E3. MCSE nell'analisi primaria (U0, own; celle con n_obs >= 100)

Scala: proporzioni (I tipo, copertura, potenza, quota usabile), unita' del coefficiente (bias, SE empirico), punti percentuali (errore relativo dello SE).

| misura | arm | combinazioni | MCSE mediana | MCSE massima |
|---|---|---|---|---|
| rej | H0 | 12.738 | 0.0077 | 0.0324 |
| rej | H1 | 12.752 | 0.0105 | 0.0314 |
| cover | H1 | 10.558 | 0.0077 | 0.0314 |
| becover | H1 | 11.432 | 0.0078 | 0.0313 |
| bias | H1 | 10.558 | 0.0093 | 0.2214 |
| empse | H1 | 11.432 | 0.0068 | 0.1781 |
| relerr_modse | H1 | 11.432 | 2.2916 | 51.4197 |
| power_sa | H1 | 12.752 | 0.0126 | 0.4026 |
| p_usable | H0 | 12.738 | 0.0000 | 0.0160 |

### E4. Valori veri di M1c in GLMM-native (quadratura, Q6c descrittivo)

| arm | T | condizionale | marginale | errore 40 vs 80 nodi |
|---|---|---|---|---|
| H0 | 2 | 0,000000 | 0,000569 | 1.1e-08 |
| H0 | 3 | 0,000000 | 0,000568 | 9.5e-09 |
| H0 | 4 | 0,000000 | 0,000569 | 8.3e-09 |
| H0 | 6 | 0,000000 | 0,000572 | 5.3e-09 |
| H1 | 2 | 0,298393 | 0,190843 | 4.1e-08 |
| H1 | 3 | 0,298393 | 0,192294 | 2.6e-08 |
| H1 | 4 | 0,298393 | 0,193990 | 1.1e-08 |
| H1 | 6 | 0,298393 | 0,197849 | 7.3e-09 |

| controllo | chiavi | valore | discordanze NA | superato |
|---|---|---|---|---|
| Q6a S2 reproduction M1-M4 | 9.600 | 5.6e-16 | 0 | si |
| Q6a failures = degenerate keys, main | 960.000 | 1.0e-08 | 0 | si |
| Q6a failures = degenerate keys, M1c | 48.000 | 1.0e-08 | 0 | si |
| Q6b replay = generator, 200 keys | 200 | 0.0e+00 | 0 | si |
| Q6c quadrature reproduces S2 M1 | 12 | 4.4e-16 | 0 | si |

### E5. Controprova rsimsum (Q5a) e identita' (Q5b)

| controllo | combinazioni | confronti | diff. max | oltre tolleranza | tolleranza |
|---|---|---|---|---|---|
| rsimsum n | 300 | 300 | 0.0e+00 | 0 | 1.0e-10 |
| rsimsum bias | 300 | 300 | 4.2e-16 | 0 | 1.0e-10 |
| rsimsum empse | 300 | 300 | 4.9e-15 | 0 | 1.0e-10 |
| rsimsum modse | 300 | 300 | 5.1e-15 | 0 | 1.0e-10 |
| rsimsum relerr | 300 | 300 | 4.4e-13 | 0 | 1.0e-10 |
| rsimsum mse | 300 | 300 | 5.3e-15 | 0 | 1.0e-10 |
| rsimsum cover | 300 | 300 | 6.7e-16 | 0 | 1.0e-10 |
| rsimsum becover | 300 | 300 | 6.7e-16 | 0 | 1.0e-10 |
| rsimsum power | 300 | 300 | 5.6e-16 | 0 | 1.0e-10 |
| rsimsum mcse_bias | 300 | 270 | 4.7e-16 | 0 | 1.0e-10 |
| rsimsum mcse_empse | 300 | 300 | 3.9e-16 | 0 | 1.0e-10 |
| rsimsum mcse_modse | 300 | 300 | 3.3e-16 | 0 | 1.0e-10 |
| rsimsum mcse_relerr | 300 | 300 | 5.0e-14 | 0 | 1.0e-10 |
| rsimsum mcse_mse | 300 | 300 | 4.7e-16 | 0 | 1.0e-10 |
| rsimsum mcse_cover | 300 | 300 | 6.2e-17 | 0 | 1.0e-10 |
| rsimsum mcse_becover | 300 | 300 | 6.2e-17 | 0 | 1.0e-10 |
| rsimsum mcse_power | 300 | 300 | 5.2e-17 | 0 | 1.0e-10 |
| identity MSE, marginal | 100.893 | 100.893 | 2.9e-14 | 0 | 1.0e-12 |
| identity MSE, conditional per replicate | 10.580 | 10.580 | 4.5e-03 | 10.052 | 1.0e-12 |
| identity MSE, conditional | 10.579 | 10.579 | 2.1e-13 | 0 | 1.0e-12 |

### E6. Effetto di D-08 sul LRT: differenza appaiata v4 - v3 (U0), media sulle celle

d I tipo/potenza: differenza del tasso di rifiuto sui replicati usabili in entrambe le versioni; d usabile: differenza della quota usabile sui replicati con disegno stimabile. Celle oltre 2 MCSE appaiate.

| gruppo | meccanismo | scenario | arm | celle | d rifiuto medio | celle oltre 2 MCSE (rifiuto) | d usabile medio | celle oltre 2 MCSE (usabile) | n appaiati minimo |
|---|---|---|---|---|---|---|---|---|---|
| D-01 + D-08 | GenBinaryY | M3 | H0 | 96 | -0,0455 | 80 | -0,1644 | 16 | 0 |
| D-01 + D-08 | GenBinaryY | M4 | H0 | 96 | -0,0455 | 80 | -0,1647 | 18 | 0 |
| D-01 + D-08 | GenBinaryY | M5 | H0 | 96 | -0,0357 | 53 | -0,3317 | 38 | 0 |
| D-02 + D-08 | GLMM-AR1 | M5 | H0 | 96 | -0,0151 | 28 | -0,3329 | 39 | 0 |
| D-02 + D-08 | GLMM-AR1 | M5 | H1 | 96 | -0,0380 | 56 | -0,3331 | 39 | 0 |
| D-08 only | GLMM-AR1 | M3 | H0 | 96 | -0,0248 | 78 | -0,1648 | 16 | 0 |
| D-08 only | GLMM-AR1 | M3 | H1 | 96 | -0,0425 | 78 | -0,1650 | 17 | 0 |
| D-08 only | GLMM-AR1 | M4 | H0 | 96 | -0,0222 | 70 | -0,1651 | 17 | 0 |
| D-08 only | GLMM-AR1 | M4 | H1 | 96 | -0,0246 | 72 | -0,1652 | 17 | 0 |
| D-08 only | GLMM-native | M3 | H0 | 96 | -0,0196 | 64 | -0,1651 | 16 | 0 |
| D-08 only | GLMM-native | M3 | H1 | 96 | -0,0376 | 79 | -0,1652 | 17 | 0 |
| D-08 only | GLMM-native | M4 | H0 | 96 | -0,0165 | 71 | -0,1653 | 17 | 0 |
| D-08 only | GLMM-native | M4 | H1 | 96 | -0,0165 | 60 | -0,1654 | 17 | 0 |
| D-08 only | GLMM-native | M5 | H0 | 96 | -0,0124 | 33 | -0,3329 | 38 | 0 |
| D-08 only | GLMM-native | M5 | H1 | 96 | -0,0290 | 54 | -0,3333 | 38 | 0 |
| D-08 only | GenBinaryY | M3 | H1 | 96 | -0,0661 | 79 | -0,1646 | 16 | 0 |
| D-08 only | GenBinaryY | M4 | H1 | 96 | -0,0490 | 80 | -0,1648 | 18 | 0 |
| D-08 only | GenBinaryY | M5 | H1 | 96 | -0,0656 | 61 | -0,3318 | 39 | 0 |

### E7. Wald glmer (default) e Wald bobyqa (D-16), media sulle celle con n_obs >= 100 per entrambi

| meccanismo | arm | celle | rifiuto glmer | rifiuto bobyqa | celle con abs(z) > 2 (z non appaiato, conservativo) | quota usabile glmer | quota usabile bobyqa |
|---|---|---|---|---|---|---|---|
| GLMM-AR1 | H0 | 412 | 0,0648 | 0,0649 | 88 | 0,945 | 0,947 |
| GLMM-AR1 | H1 | 414 | 0,2789 | 0,2800 | 59 | 0,946 | 0,948 |
| GLMM-native | H0 | 414 | 0,0621 | 0,0616 | 80 | 0,944 | 0,946 |
| GLMM-native | H1 | 414 | 0,2682 | 0,2697 | 47 | 0,949 | 0,951 |
| GenBinaryY | H0 | 412 | 0,0857 | 0,0825 | 117 | 0,930 | 0,934 |
| GenBinaryY | H1 | 413 | 0,3036 | 0,3036 | 98 | 0,934 | 0,938 |

