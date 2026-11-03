### Tabella E. Attese pre-registrate e osservato

| id | attesa | tolleranza | osservato | superato |
|---|---|---|---|---|
| Q0 | md5 di results_full/ e results_v4 (S2) invariati prima e dopo | zero falliti | prima: Q0 results_full OK 28800 files; prima: Q0 results_v4 (S2) OK 126 files; dopo: Q0 results_full OK 28800 files; dopo: Q0 results_v4 (S2) OK 126 files | si |
| Q1 | parita' H0/H1: 15/15 + 3/3 (M1c); controprova sui generatori originali 12/15, falliti GenBinaryY M3-M5 | esatta | v4 15/15 + 3/3; originali 12/15, falliti: GenBinaryY M3, GenBinaryY M4, GenBinaryY M5 | si |
| Q2 | AR1 M5 v4 differisce dall'originale solo in D[3,3] (0,01 -> 0,001) e D = D di GLMM-native M5 | esatta | 8/8 | si |
| Q3 | banco con codice originale riproduce results_full/ bit per bit, 300 replicati | 300/300 | 300/300 (righe di traccia escluse: 4) | si |
| Q4 | disegno identico al generatore di riferimento (F1, F2: originale; C: M1) | 100% | F1 288.000/288.000; F2 191.888/191.888; C 143.914/143.914 | si |
| Q5 | blocco L: gee rieseguita identica a results_full/ | 100% | 12.480/12.480 chunk identici; sim_id uguali 12.480/12.480; righe 5.277.866 contro 5.277.866 | si |
| Q6 | 18.720 chunk, zero errori, ogni sim_id con una traccia | zero | 18.720/18.720 chunk; sim completi 18.720; errori di crew: vedi run.log | si |
| Q7 | design_degenerate in biiezione con le chiavi di S2; errori di generazione solo su chiavi degeneri | esatta | main 1.728/1.728 celle, C 144/144 celle; errori su chiavi non degeneri: 0 | si |
| Q8a | re_structure del LRT: pendenze in M3-M5, (1\|id) in M1c | 100% | F1 2.880/2.880; F2 1.920/1.920; L 12.480/12.480; C 1.440/1.440 | si |
| Q8b | abs(LRT - Wald) < 1e-3 in almeno il 95% delle coppie senza warning, per scenario x meccanismo | >= 0,95 in ogni riga | 6/12 righe; quota minima 0,608 (GenBinaryY M5) | no |
| Q9 | GenBinaryY M3-M4: differenza appaiata H0 v4 - H0 originale entro 3 MCSE | >= 99% | usable 469/474 (0,9895); type I 1.190/1.347 (0,8834) | no |
| Q10.1 | beta0(T) entro 5e-4 dai valori pre-registrati, residuo < 1e-10 | 4/4 | 4/4; beta0 = 0,2427; 0,1960; 0,1495; 0,0569 | si |
| Q10.2 | prevalenza attesa H0: 0,5945 (GenBinaryY), 0,589 +/- 0,001 (GLMM-native) | 4/4 | 4/4; GLMM-native 0,5889; 0,5889; 0,5889; 0,5890 | si |
| Q10.3 | prevalenza realizzata a N = 45 entro 3 MCSE dall'attesa | 24/24 | 24/24; \|z\| massimo 1,91 | si |
| Q10.4 | AR1: errore massimo di calibrazione < 1e-8 sui replicati non degeneri | < 1e-8 | 9.51e-09 | si |
| Q10.5 | costanti e argomenti dei motori come D-14 | tutti | 18/18 | si |

### E1. Parita' H0/H1 per generatore (Q1)

| dgp | scen | v4 | v4_differisce_in | originale | originale_differisce_in |
|---|---|---|---|---|---|
| GLMM-AR1 | M1 | 12/12 | target_beta_m | 12/12 | target_beta_m |
| GLMM-AR1 | M1c | 9/9 | target_beta_m | NA | NA |
| GLMM-AR1 | M2 | 12/12 | target_beta_m | 12/12 | target_beta_m |
| GLMM-AR1 | M3 | 12/12 | target_beta_m | 12/12 | target_beta_m |
| GLMM-AR1 | M4 | 12/12 | target_beta_m | 12/12 | target_beta_m |
| GLMM-AR1 | M5 | 12/12 | target_beta_m | 12/12 | target_beta_m |
| GLMM-native | M1 | 12/12 | beta_c | 12/12 | beta_c |
| GLMM-native | M1c | 9/9 | beta_c | NA | NA |
| GLMM-native | M2 | 12/12 | beta_c | 12/12 | beta_c |
| GLMM-native | M3 | 12/12 | beta_c | 12/12 | beta_c |
| GLMM-native | M4 | 12/12 | beta_c | 12/12 | beta_c |
| GLMM-native | M5 | 12/12 | beta_c | 12/12 | beta_c |
| GenBinaryY | M1 | 12/12 | beta | 12/12 | beta |
| GenBinaryY | M1c | 9/9 | beta | NA | NA |
| GenBinaryY | M2 | 12/12 | beta | 12/12 | beta |
| GenBinaryY | M3 | 12/12 | beta | 0/12 | beta,sigma |
| GenBinaryY | M4 | 12/12 | beta | 0/12 | beta,sigma,gamma |
| GenBinaryY | M5 | 12/12 | beta | 0/12 | beta,sigma,gamma |

### E2. Calibrazione di M1c (Q10.1, Q10.2)

| T | beta0 | preregistrato | residuo | prev. H1 marginale | GLMM-native H0 | GLMM-native H1 |
|---|---|---|---|---|---|---|
| 2 | 0,2427 | 0,2427 | 0.0e+00 | 0,6038 | 0,5889 | 0,5978 |
| 3 | 0,1960 | 0,1960 | 0.0e+00 | 0,6139 | 0,5889 | 0,6076 |
| 4 | 0,1495 | 0,1495 | -1.1e-16 | 0,6233 | 0,5889 | 0,6170 |
| 6 | 0,0569 | 0,0569 | -1.1e-16 | 0,6399 | 0,5890 | 0,6339 |

### E3. Prevalenza realizzata di M1c a N = 45, replicati 1-1.000 (Q10.3)

| meccanismo | T | arm | generati | prevalenza | MCSE | attesa | z |
|---|---|---|---|---|---|---|---|
| GLMM-AR1 | 2 | H0 | 1000 | 0,5926 | 0,0020 | 0,5945 | -0,95 |
| GLMM-AR1 | 2 | H1 | 1000 | 0,6019 | 0,0020 | 0,6038 | -0,94 |
| GLMM-AR1 | 3 | H0 | 1000 | 0,5938 | 0,0018 | 0,5945 | -0,39 |
| GLMM-AR1 | 3 | H1 | 1000 | 0,6139 | 0,0018 | 0,6139 | 0,02 |
| GLMM-AR1 | 4 | H0 | 1000 | 0,5921 | 0,0017 | 0,5945 | -1,42 |
| GLMM-AR1 | 4 | H1 | 1000 | 0,6207 | 0,0017 | 0,6233 | -1,57 |
| GLMM-AR1 | 6 | H0 | 1000 | 0,5929 | 0,0016 | 0,5945 | -1,01 |
| GLMM-AR1 | 6 | H1 | 1000 | 0,6385 | 0,0016 | 0,6399 | -0,87 |
| GLMM-native | 2 | H0 | 1000 | 0,5900 | 0,0020 | 0,5889 | 0,60 |
| GLMM-native | 2 | H1 | 1000 | 0,6001 | 0,0019 | 0,5978 | 1,17 |
| GLMM-native | 3 | H0 | 1000 | 0,5912 | 0,0018 | 0,5889 | 1,25 |
| GLMM-native | 3 | H1 | 1000 | 0,6111 | 0,0018 | 0,6076 | 1,91 |
| GLMM-native | 4 | H0 | 1000 | 0,5884 | 0,0018 | 0,5889 | -0,30 |
| GLMM-native | 4 | H1 | 1000 | 0,6167 | 0,0017 | 0,6170 | -0,15 |
| GLMM-native | 6 | H0 | 1000 | 0,5892 | 0,0017 | 0,5890 | 0,09 |
| GLMM-native | 6 | H1 | 1000 | 0,6341 | 0,0016 | 0,6339 | 0,12 |
| GenBinaryY | 2 | H0 | 1000 | 0,5955 | 0,0021 | 0,5945 | 0,48 |
| GenBinaryY | 2 | H1 | 1000 | 0,6052 | 0,0021 | 0,6038 | 0,64 |
| GenBinaryY | 3 | H0 | 1000 | 0,5952 | 0,0020 | 0,5945 | 0,33 |
| GenBinaryY | 3 | H1 | 1000 | 0,6153 | 0,0020 | 0,6139 | 0,74 |
| GenBinaryY | 4 | H0 | 1000 | 0,5913 | 0,0020 | 0,5945 | -1,61 |
| GenBinaryY | 4 | H1 | 1000 | 0,6200 | 0,0019 | 0,6233 | -1,70 |
| GenBinaryY | 6 | H0 | 1000 | 0,5929 | 0,0018 | 0,5945 | -0,89 |
| GenBinaryY | 6 | H1 | 1000 | 0,6376 | 0,0017 | 0,6399 | -1,28 |

### E4. Completezza, disegno e identita' per blocco (Q4-Q7)

| blocco | attesi | scritti | sim completi | errori di generazione | righe di traccia di analisi | errori su chiavi non degeneri |
|---|---|---|---|---|---|---|
| C | 1.440 | 1.440 | 1.440 | 86 | 0 | 0 |
| F1 | 2.880 | 2.880 | 2.880 | 0 | 0 | 0 |
| F2 | 1.920 | 1.920 | 1.920 | 112 | 0 | 0 |
| L | 12.480 | 12.480 | 12.480 | 172 | 0 | 0 |

### E5. LRT contro Wald senza warning (Q8b)

| meccanismo | scenario | coppie | entro | quota | superato |
|---|---|---|---|---|---|
| GenBinaryY | M3 | 124.154 | 117.155 | 0,9436 | no |
| GenBinaryY | M4 | 123.290 | 116.848 | 0,9477 | no |
| GenBinaryY | M5 | 8.368 | 5.089 | 0,6082 | no |
| GLMM-AR1 | M5 | 14.935 | 11.358 | 0,7605 | no |
| GLMM-native | M3 | 124.397 | 118.686 | 0,9541 | si |
| GLMM-native | M4 | 125.264 | 120.151 | 0,9592 | si |
| GLMM-native | M5 | 16.942 | 13.338 | 0,7873 | no |
| GLMM-AR1 | M3 | 125.661 | 119.948 | 0,9545 | si |
| GLMM-AR1 | M4 | 124.701 | 119.618 | 0,9592 | si |
| GenBinaryY | M1c | 41.861 | 39.767 | 0,9500 | no |
| GLMM-native | M1c | 42.864 | 41.196 | 0,9611 | si |
| GLMM-AR1 | M1c | 42.865 | 41.162 | 0,9603 | si |

### E6. Effetto di D-01 su GenBinaryY M3-M5 H0 (Q9)

M3-M4: 1.480 combinazioni cella x metodo; superamenti di 3 MCSE elencati sotto (162 righe). M5 esplorativo: I tipo 54 su 663 oltre 3 MCSE, quota usabile 43 su 331.

| M | N | T | metodo | d usabile | z_u | I tipo orig | I tipo v4 | z_t1 |
|---|---|---|---|---|---|---|---|---|
| 3 | 10 | 2 | glmer_lrt | -0,9656 | -156,53 | 0,0424 | 0,0000 | -1,00 |
| 4 | 10 | 2 | glmer_lrt | -0,9724 | -174,82 | 0,0449 | 0,0417 | -1,45 |
| 3 | 15 | 2 | glmer_lrt | -0,9968 | -540,98 | 0,0586 | 0,0000 |  NA |
| 4 | 15 | 2 | glmer_lrt | -0,9948 | -428,88 | 0,0728 | 0,0000 |  NA |
| 3 | 10 | 3 | glmer_lrt | -0,0042 | -1,27 | 0,0920 | 0,0202 | -8,18 |
| 4 | 10 | 3 | glmer_lrt | -0,0137 | -2,99 | 0,1092 | 0,0247 | -8,64 |
| 3 | 15 | 3 | glmer_lrt | -0,0010 | -0,45 | 0,1239 | 0,0325 | -9,19 |
| 4 | 15 | 3 | glmer_lrt | -0,0051 | -2,24 | 0,1204 | 0,0498 | -6,89 |
| 3 | 20 | 3 | glmer_lrt | -0,0010 | -1,00 | 0,1473 | 0,0545 | -8,71 |
| 4 | 20 | 3 | glmer_lrt | -0,0010 | -0,58 | 0,1424 | 0,0552 | -7,96 |
| 3 | 25 | 3 | glmer_lrt | 0,0010 | 1,00 | 0,1542 | 0,0540 | -8,98 |
| 4 | 25 | 3 | glmer_lrt | 0,0000 |  NA | 0,1393 | 0,0601 | -7,19 |
| 3 | 35 | 3 | glmer_lrt | 0,0000 |  NA | 0,1500 | 0,0640 | -7,25 |
| 4 | 35 | 3 | glmer_lrt | -0,0010 | -1,00 | 0,1280 | 0,0631 | -5,86 |
| 3 | 45 | 3 | glmer_lrt | 0,0000 |  NA | 0,1250 | 0,0770 | -4,35 |
| 4 | 45 | 3 | glmer_lrt | 0,0000 |  NA | 0,1200 | 0,0720 | -4,28 |
| 3 | 55 | 3 | glmer_lrt | 0,0000 |  NA | 0,1160 | 0,0610 | -4,90 |
| 4 | 55 | 3 | glmer_lrt | 0,0000 |  NA | 0,1350 | 0,0660 | -6,40 |
| 3 | 65 | 3 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0550 | -3,01 |
| 4 | 65 | 3 | glmer_lrt | 0,0000 |  NA | 0,1100 | 0,0650 | -4,12 |
| 3 | 75 | 3 | glmer_lrt | 0,0000 |  NA | 0,1000 | 0,0540 | -4,35 |
| 4 | 75 | 3 | glmer_lrt | 0,0000 |  NA | 0,0840 | 0,0480 | -3,66 |
| 3 | 85 | 3 | glmer_lrt | 0,0000 |  NA | 0,0970 | 0,0550 | -3,86 |
| 3 | 95 | 3 | glmer_lrt | 0,0000 |  NA | 0,0950 | 0,0480 | -4,50 |
| 4 | 95 | 3 | glmer_lrt | 0,0000 |  NA | 0,0870 | 0,0560 | -3,01 |
| 3 | 105 | 3 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0410 | -4,39 |
| 3 | 115 | 3 | glmer_lrt | 0,0000 |  NA | 0,0890 | 0,0460 | -4,15 |
| 4 | 115 | 3 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0420 | -4,61 |
| 4 | 125 | 3 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0360 | -5,06 |
| 3 | 135 | 3 | glmer_lrt | 0,0000 |  NA | 0,0740 | 0,0430 | -3,07 |
| 4 | 135 | 3 | glmer_lrt | 0,0000 |  NA | 0,0910 | 0,0390 | -5,27 |
| 4 | 145 | 3 | glmer_lrt | 0,0000 |  NA | 0,0760 | 0,0430 | -3,40 |
| 3 | 10 | 4 | glmer_lrt | 0,0000 | 0,00 | 0,0985 | 0,0359 | -7,80 |
| 4 | 10 | 4 | glmer_lrt | -0,0010 | -0,45 | 0,1172 | 0,0327 | -8,85 |
| 3 | 15 | 4 | glmer_lrt | 0,0010 | 1,00 | 0,1158 | 0,0312 | -9,01 |
| 4 | 15 | 4 | glmer_lrt | 0,0020 | 1,41 | 0,1216 | 0,0371 | -9,01 |
| 3 | 20 | 4 | glmer_lrt | 0,0000 |  NA | 0,1230 | 0,0370 | -9,24 |
| 4 | 20 | 4 | glmer_lrt | 0,0010 | 1,00 | 0,1381 | 0,0550 | -8,94 |
| 3 | 25 | 4 | glmer_lrt | 0,0000 |  NA | 0,1121 | 0,0330 | -8,90 |
| 4 | 25 | 4 | glmer_lrt | 0,0000 |  NA | 0,1130 | 0,0400 | -8,07 |
| 3 | 35 | 4 | glmer_lrt | 0,0000 |  NA | 0,0970 | 0,0480 | -6,29 |
| 4 | 35 | 4 | glmer_lrt | 0,0000 |  NA | 0,1000 | 0,0370 | -7,47 |
| 3 | 45 | 4 | glmer_lrt | 0,0000 |  NA | 0,0880 | 0,0310 | -6,92 |
| 4 | 45 | 4 | glmer_lrt | 0,0000 |  NA | 0,1020 | 0,0390 | -7,81 |
| 3 | 55 | 4 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0300 | -7,22 |
| 4 | 55 | 4 | glmer_lrt | 0,0000 |  NA | 0,0910 | 0,0280 | -8,06 |
| 3 | 65 | 4 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0240 | -6,79 |
| 4 | 65 | 4 | glmer_lrt | 0,0000 |  NA | 0,0710 | 0,0260 | -6,30 |
| 3 | 75 | 4 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0240 | -7,52 |
| 4 | 75 | 4 | glmer_lrt | 0,0000 |  NA | 0,1000 | 0,0320 | -8,28 |
| 3 | 85 | 4 | glmer_lrt | 0,0000 |  NA | 0,0900 | 0,0250 | -7,83 |
| 4 | 85 | 4 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0210 | -7,41 |
| 3 | 95 | 4 | glmer_lrt | 0,0000 |  NA | 0,0700 | 0,0220 | -6,32 |
| 4 | 95 | 4 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0270 | -6,21 |
| 3 | 105 | 4 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0280 | -6,59 |
| 4 | 105 | 4 | glmer_lrt | 0,0000 |  NA | 0,0760 | 0,0270 | -6,89 |
| 3 | 115 | 4 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0290 | -7,27 |
| 4 | 115 | 4 | glmer_lrt | 0,0000 |  NA | 0,0530 | 0,0300 | -4,16 |
| 3 | 125 | 4 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0190 | -8,00 |
| 4 | 125 | 4 | glmer_lrt | 0,0000 |  NA | 0,0740 | 0,0280 | -6,79 |
| 3 | 135 | 4 | glmer_lrt | 0,0000 |  NA | 0,0700 | 0,0140 | -7,42 |
| 4 | 135 | 4 | glmer_lrt | 0,0000 |  NA | 0,0870 | 0,0360 | -6,79 |
| 3 | 145 | 4 | glmer_lrt | 0,0000 |  NA | 0,0700 | 0,0210 | -6,63 |
| 4 | 145 | 4 | glmer_lrt | 0,0000 |  NA | 0,0880 | 0,0310 | -7,37 |
| 3 | 10 | 6 | glmer_lrt | -0,0010 | -1,00 | 0,1181 | 0,0571 | -7,55 |
| 4 | 10 | 6 | glmer_lrt | 0,0010 | 1,00 | 0,1256 | 0,0602 | -8,07 |
| 3 | 15 | 6 | glmer_lrt | 0,0000 |  NA | 0,1140 | 0,0580 | -6,95 |
| 4 | 15 | 6 | glmer_lrt | 0,0000 |  NA | 0,1160 | 0,0510 | -7,83 |
| 3 | 20 | 6 | glmer_lrt | 0,0000 |  NA | 0,0680 | 0,0280 | -6,29 |
| 4 | 20 | 6 | glmer_lrt | 0,0000 |  NA | 0,0990 | 0,0440 | -7,63 |
| 3 | 25 | 6 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0450 | -6,59 |
| 4 | 25 | 6 | glmer_lrt | 0,0000 |  NA | 0,0880 | 0,0450 | -6,40 |
| 3 | 35 | 6 | glmer_lrt | 0,0000 |  NA | 0,1020 | 0,0480 | -7,27 |
| 4 | 35 | 6 | glmer_lrt | 0,0000 |  NA | 0,0940 | 0,0600 | -5,32 |
| 3 | 45 | 6 | glmer_lrt | 0,0000 |  NA | 0,0650 | 0,0350 | -4,80 |
| 4 | 45 | 6 | glmer_lrt | 0,0000 |  NA | 0,0900 | 0,0460 | -6,63 |
| 3 | 55 | 6 | glmer_lrt | 0,0000 |  NA | 0,0760 | 0,0360 | -5,75 |
| 4 | 55 | 6 | glmer_lrt | 0,0000 |  NA | 0,0840 | 0,0560 | -4,72 |
| 3 | 65 | 6 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0490 | -5,07 |
| 4 | 65 | 6 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0450 | -5,06 |
| 3 | 75 | 6 | glmer_lrt | 0,0000 |  NA | 0,1010 | 0,0500 | -7,04 |
| 4 | 75 | 6 | glmer_lrt | 0,0000 |  NA | 0,0820 | 0,0570 | -3,52 |
| 3 | 85 | 6 | glmer_lrt | 0,0000 |  NA | 0,0790 | 0,0470 | -4,77 |
| 4 | 85 | 6 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0510 | -5,96 |
| 3 | 95 | 6 | glmer_lrt | 0,0000 |  NA | 0,0720 | 0,0520 | -4,11 |
| 4 | 95 | 6 | glmer_lrt | 0,0000 |  NA | 0,0790 | 0,0540 | -4,26 |
| 3 | 105 | 6 | glmer_lrt | 0,0000 |  NA | 0,0760 | 0,0530 | -4,31 |
| 4 | 105 | 6 | glmer_lrt | 0,0000 |  NA | 0,0890 | 0,0600 | -4,96 |
| 3 | 115 | 6 | glmer_lrt | 0,0000 |  NA | 0,0840 | 0,0490 | -5,55 |
| 4 | 115 | 6 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0450 | -5,69 |
| 3 | 125 | 6 | glmer_lrt | 0,0000 |  NA | 0,0890 | 0,0480 | -6,38 |
| 4 | 125 | 6 | glmer_lrt | 0,0000 |  NA | 0,0790 | 0,0480 | -5,31 |
| 3 | 135 | 6 | glmer_lrt | 0,0000 |  NA | 0,0870 | 0,0440 | -6,40 |
| 4 | 135 | 6 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0460 | -5,28 |
| 3 | 145 | 6 | glmer_lrt | 0,0000 |  NA | 0,0730 | 0,0450 | -5,01 |
| 4 | 145 | 6 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0560 | -5,64 |
| 3 | 10 | 8 | glmer_lrt | -0,0010 | -1,00 | 0,1301 | 0,0711 | -7,78 |
| 4 | 10 | 8 | glmer_lrt | 0,0000 |  NA | 0,1142 | 0,0581 | -7,30 |
| 3 | 15 | 8 | glmer_lrt | 0,0000 |  NA | 0,1230 | 0,0490 | -8,81 |
| 4 | 15 | 8 | glmer_lrt | 0,0000 |  NA | 0,1050 | 0,0560 | -6,51 |
| 3 | 20 | 8 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0400 | -5,96 |
| 4 | 20 | 8 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0500 | -5,81 |
| 3 | 25 | 8 | glmer_lrt | 0,0000 |  NA | 0,0980 | 0,0400 | -7,70 |
| 4 | 25 | 8 | glmer_lrt | 0,0000 |  NA | 0,0770 | 0,0440 | -5,36 |
| 3 | 35 | 8 | glmer_lrt | 0,0000 |  NA | 0,1020 | 0,0580 | -6,22 |
| 4 | 35 | 8 | glmer_lrt | 0,0000 |  NA | 0,0970 | 0,0490 | -6,55 |
| 3 | 45 | 8 | glmer_lrt | 0,0000 |  NA | 0,0800 | 0,0490 | -5,02 |
| 4 | 45 | 8 | glmer_lrt | 0,0000 |  NA | 0,0760 | 0,0460 | -5,06 |
| 3 | 55 | 8 | glmer_lrt | 0,0000 |  NA | 0,0760 | 0,0510 | -4,03 |
| 4 | 55 | 8 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0550 | -4,36 |
| 3 | 65 | 8 | glmer_lrt | 0,0000 |  NA | 0,0880 | 0,0540 | -5,19 |
| 4 | 65 | 8 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0560 | -4,92 |
| 3 | 75 | 8 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0530 | -5,09 |
| 4 | 75 | 8 | glmer_lrt | 0,0000 |  NA | 0,0770 | 0,0520 | -4,39 |
| 3 | 85 | 8 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0560 | -5,21 |
| 4 | 85 | 8 | glmer_lrt | 0,0000 |  NA | 0,0980 | 0,0650 | -4,98 |
| 3 | 95 | 8 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0410 | -5,19 |
| 4 | 95 | 8 | glmer_lrt | 0,0000 |  NA | 0,0960 | 0,0660 | -4,92 |
| 3 | 105 | 8 | glmer_lrt | 0,0000 |  NA | 0,0750 | 0,0470 | -4,86 |
| 4 | 105 | 8 | glmer_lrt | 0,0000 |  NA | 0,0900 | 0,0510 | -5,78 |
| 3 | 115 | 8 | glmer_lrt | 0,0000 |  NA | 0,0820 | 0,0550 | -4,48 |
| 4 | 115 | 8 | glmer_lrt | 0,0000 |  NA | 0,0940 | 0,0600 | -5,32 |
| 3 | 125 | 8 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0590 | -4,80 |
| 4 | 125 | 8 | glmer_lrt | 0,0000 |  NA | 0,0780 | 0,0500 | -5,01 |
| 3 | 135 | 8 | glmer_lrt | 0,0000 |  NA | 0,0800 | 0,0490 | -5,31 |
| 4 | 135 | 8 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0570 | -5,17 |
| 3 | 145 | 8 | glmer_lrt | 0,0000 |  NA | 0,0720 | 0,0490 | -3,92 |
| 4 | 145 | 8 | glmer_lrt | 0,0000 |  NA | 0,0910 | 0,0640 | -4,90 |
| 3 | 10 | 10 | glmer_lrt | 0,0000 |  NA | 0,1050 | 0,0540 | -6,91 |
| 3 | 10 | 10 | mmm | -0,0172 | -3,07 | 0,0992 | 0,0953 | 0,00 |
| 4 | 10 | 10 | glmer_lrt | 0,0000 |  NA | 0,1255 | 0,0622 | -8,06 |
| 3 | 15 | 10 | glmer_lrt | 0,0000 |  NA | 0,1040 | 0,0500 | -7,27 |
| 4 | 15 | 10 | glmer_lrt | 0,0000 |  NA | 0,1080 | 0,0600 | -6,43 |
| 3 | 20 | 10 | glmer_lrt | 0,0000 |  NA | 0,1070 | 0,0660 | -5,96 |
| 4 | 20 | 10 | glmer_lrt | 0,0000 |  NA | 0,0990 | 0,0460 | -6,83 |
| 3 | 25 | 10 | gee_fg | -0,0010 | -1,00 | 0,0782 | 0,0662 | -3,01 |
| 3 | 25 | 10 | glmer_lrt | 0,0000 |  NA | 0,0900 | 0,0530 | -5,87 |
| 4 | 25 | 10 | glmer_lrt | 0,0000 |  NA | 0,1020 | 0,0430 | -7,40 |
| 3 | 35 | 10 | glmer_lrt | 0,0000 |  NA | 0,1000 | 0,0670 | -5,09 |
| 4 | 35 | 10 | glmer_lrt | 0,0000 |  NA | 0,1100 | 0,0760 | -5,07 |
| 3 | 45 | 10 | glmer_lrt | 0,0000 |  NA | 0,1020 | 0,0560 | -6,26 |
| 4 | 45 | 10 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0490 | -5,36 |
| 3 | 55 | 10 | glmer_lrt | 0,0000 |  NA | 0,0910 | 0,0580 | -5,36 |
| 4 | 55 | 10 | glmer_lrt | 0,0000 |  NA | 0,0980 | 0,0660 | -5,00 |
| 3 | 65 | 10 | glmer_lrt | 0,0000 |  NA | 0,0870 | 0,0610 | -4,25 |
| 4 | 65 | 10 | glmer_lrt | 0,0000 |  NA | 0,0960 | 0,0560 | -5,63 |
| 3 | 75 | 10 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0600 | -5,12 |
| 4 | 75 | 10 | glmer_lrt | 0,0000 |  NA | 0,0970 | 0,0510 | -6,64 |
| 3 | 85 | 10 | glmer_lrt | 0,0000 |  NA | 0,0860 | 0,0480 | -5,96 |
| 4 | 85 | 10 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0490 | -6,11 |
| 3 | 95 | 10 | glmer_lrt | 0,0000 |  NA | 0,0880 | 0,0530 | -5,41 |
| 4 | 95 | 10 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0480 | -5,73 |
| 3 | 105 | 10 | glmer_lrt | 0,0000 |  NA | 0,0980 | 0,0520 | -6,38 |
| 4 | 105 | 10 | glmer_lrt | 0,0000 |  NA | 0,0960 | 0,0610 | -5,41 |
| 3 | 115 | 10 | glmer_lrt | 0,0000 |  NA | 0,0770 | 0,0520 | -4,26 |
| 4 | 115 | 10 | glmer_lrt | 0,0000 |  NA | 0,0850 | 0,0510 | -5,60 |
| 3 | 125 | 10 | glmer_lrt | 0,0000 |  NA | 0,1200 | 0,0770 | -6,26 |
| 4 | 125 | 10 | glmer_lrt | 0,0000 |  NA | 0,0920 | 0,0650 | -4,48 |
| 3 | 135 | 10 | glmer_lrt | 0,0000 |  NA | 0,0950 | 0,0470 | -6,67 |
| 4 | 135 | 10 | glmer_lrt | 0,0000 |  NA | 0,0900 | 0,0520 | -5,69 |
| 3 | 145 | 10 | glmer_lrt | 0,0000 |  NA | 0,0830 | 0,0590 | -4,15 |
| 4 | 145 | 10 | glmer_lrt | 0,0000 |  NA | 0,1000 | 0,0690 | -5,16 |

Media della differenza appaiata dell'errore di I tipo in M3-M4, per metodo (cella pesata 1): glmer_lrt -0,0460; gee 0,0002; gee_kc -0,0002; mmm -0,0001; gee_fg 0,0003; gee_md 0,0006; shrink 0,0004; glmer 0,0005.

### E6b. Q9 per metodo (scomposizione diagnostica, non pre-registrata)

In F1 `glmer_lrt` cambia anche per D-08 (pendenze casuali nel modello del LRT), quindi la sua differenza appaiata non misura solo D-01.

| scenari | metodo | I tipo valutabili | I tipo oltre 3 MCSE | media d I tipo | usabile valutabili | usabile oltre 3 MCSE |
|---|---|---|---|---|---|---|
| M3-M4 | gee | 177 | 0 | 0,0002 | 57 | 0 |
| M3-M4 | gee_fg | 171 | 1 | 0,0003 | 51 | 0 |
| M3-M4 | gee_kc | 178 | 0 | -0,0002 | 52 | 0 |
| M3-M4 | gee_md | 171 | 0 | 0,0006 | 51 | 0 |
| M3-M4 | glmer | 156 | 0 | 0,0005 | 56 | 0 |
| M3-M4 | glmer_lrt | 162 | 156 | -0,0460 | 20 | 4 |
| M3-M4 | mmm | 173 | 0 | -0,0001 | 152 | 1 |
| M3-M4 | shrink | 159 | 0 | 0,0004 | 35 | 0 |
| M5 (esplorativo) | gee | 96 | 1 | 0,0013 | 40 | 0 |
| M5 (esplorativo) | gee_fg | 93 | 0 | 0,0012 | 41 | 0 |
| M5 (esplorativo) | gee_kc | 93 | 1 | 0,0010 | 41 | 0 |
| M5 (esplorativo) | gee_md | 93 | 0 | 0,0011 | 41 | 0 |
| M5 (esplorativo) | glmer | 63 | 0 | 0,0058 | 30 | 3 |
| M5 (esplorativo) | glmer_lrt | 65 | 52 | -0,0368 | 14 | 8 |
| M5 (esplorativo) | mmm | 95 | 0 | -0,0014 | 95 | 2 |
| M5 (esplorativo) | shrink | 65 | 0 | 0,0084 | 29 | 1 |

Senza `glmer_lrt`, M3-M4: I tipo 1184/1185 entro 3 MCSE, quota usabile 453/454.

### E5b. Indagine su Q8b: rifit di 400 replicati campionati (esplorativa, non pre-registrata)

Campione: per GenBinaryY M5, GLMM-AR1 M5, GLMM-native M5, GenBinaryY M3 e GenBinaryY M1c, 60 coppie discordanti e 20 concordanti da 60 chunk estratti a caso (seed 20261017). d logLik = logLik(LRT completo, bobyqa) - logLik(Wald, ottimizzatore di default).

| meccanismo | scenario | concordante | n | riproduce Wald | riproduce LRT | stessi theta | singolare | abs(d logLik) < 1e-3 | LRT piu' alta | Wald piu' alta | mediana abs(d stima) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GenBinaryY | M5 | no | 60 | 60 | 60 | 60 | 60 | 29 | 30 | 1 | 0.13 |
| GenBinaryY | M5 | si | 20 | 20 | 20 | 20 | 20 | 20 | 0 | 0 | 2.5e-05 |
| GLMM-AR1 | M5 | no | 60 | 60 | 60 | 60 | 60 | 21 | 35 | 4 | 0.022 |
| GLMM-AR1 | M5 | si | 20 | 20 | 20 | 20 | 20 | 20 | 0 | 0 | 2e-05 |
| GLMM-native | M5 | no | 60 | 60 | 60 | 60 | 60 | 29 | 28 | 3 | 0.013 |
| GLMM-native | M5 | si | 20 | 20 | 20 | 20 | 19 | 20 | 0 | 0 | 3e-05 |
| GenBinaryY | M3 | no | 60 | 60 | 60 | 60 | 51 | 25 | 24 | 11 | 0.15 |
| GenBinaryY | M3 | si | 20 | 20 | 20 | 20 | 5 | 20 | 0 | 0 | 2.4e-06 |
| GenBinaryY | M1c | no | 60 | 60 | 60 | 60 | 12 | 22 | 6 | 32 | 0.13 |
| GenBinaryY | M1c | si | 20 | 20 | 20 | 20 | 0 | 20 | 0 | 0 | 3.2e-06 |

### E7. Tempo di calcolo per blocco

| blocco | chunk | ore-worker | s/chunk mediana | s/chunk massimo |
|---|---|---|---|---|
| C | 1.440 | 38,3 | 89,7 | 279,6 |
| L | 12.480 | 383,5 | 84,5 | 968,6 |
| F2 | 1.920 | 270,3 | 432,1 | 1585,8 |
| F1 | 2.880 | 308,6 | 301,1 | 2026,9 |

