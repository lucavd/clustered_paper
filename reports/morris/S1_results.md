| # | attesa (pre-registrata, commit b0826c5) | osservato | esito |
|---|---|---|---|
| P1 | 30 combinazioni; valori = audit B1, B3, B4 (esatti) | 30 combinazioni; 99/101 controlli; falliti: GLMM-native M2 H1 beta_c == 1.16 x beta_m (beta_c = -2.32, 0.116, 0.232, 0.012, 0.232; 1.16 x beta_m = -2.32, 0.116, 0.232, 0.0116, 0.232); GLMM-native M2 H0 beta_c == 1.16 x beta_m (beta_c = -2.32, 0.116, 0.232, 0.012, 0; 1.16 x beta_m = -2.32, 0.116, 0.232, 0.0116, 0) | non superato alla lettera (vedi sez. 3) |
| P2 | 12/15 coppie con parita; eccezioni GenBinaryY M3, M4, M5 | 12/15; eccezioni: GenBinaryY M3 (also differ: sigma); GenBinaryY M4 (also differ: sigma, gamma); GenBinaryY M5 (also differ: sigma, gamma) | superato |
| P3 | 240/240 dataset bit-identici dal solo contenuto della tabella | 240/240 (errori di generazione in entrambi i percorsi: 0) | superato |
| P4 | 96/96 dataset M1-M2 identici fra codice del job 23035d84 e `.rundir/` | 96/96 | superato |
| P5 | 30 rifit i = 1, N = 25, T = 4: max abs diff <= 1e-8 | 30/30; 606 righe confrontate; max abs diff stima 0, SE 0; status identici 30/30 | superato (secondo run, vedi sez. 3) |
| P5' | (non pre-registrato) stesso rifit sotto Mersenne-Twister | 0/30; max abs diff stima 1.59e+07 | controprova del RNGkind |
| P6 | `md5sum -c` su `results_full/` prima e dopo: 28.800/28.800 | before (2026-10-15 21:55:21): 28800 OK, 0 falliti; after (2026-10-15 22:07:01): 28800 OK, 0 falliti; final (2026-10-15 22:15:56): 28800 OK, 0 falliti | superato |
