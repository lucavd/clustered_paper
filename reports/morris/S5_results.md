### Tabella E. Cancelli e previsioni

| id | atteso | osservato | esito |
|---|---|---|---|
| G1 | md5 invariati prima e dopo | prima: results_full: 28800 checked, 0 FAIL lines; results_v4 S2: 126 checked, 0 FAIL lines; results_v4 S3: 18720 checked, 0 FAIL lines; results_v4 S4: 5 checked, 0 FAIL lines \|\| dopo: results_full: 28800 checked, 0 FAIL lines; results_v4 S2: 126 checked, 0 FAIL lines; results_v4 S3: 18720 checked, 0 FAIL lines; results_v4 S4: 5 checked, 0 FAIL lines | superato |
| G2 | b = 0 e MCSE di S4 riprodotte entro 1e-12 | 7 of 7 checks pass (max 6.7e-16) | superato (dopo la correzione dello strumento, sez. 3) |
| G3 | SE bootstrap = SE per linearizzazione entro 7% | 27 of 27 pass; rel. diff. -0.033 to 0.020 | superato |
| G4 | ogni riga con un esito; >= 1.990 replicati finiti | missing 0; sub-claims with < B - 10 finite: 0 | superato |
| P1 | numeri v3 riprodotti alla cifra stampata | 275 of 288 reproduced, 0 within one unit, 13 not reproduced | smentita, spiegata |
| P2 | v4 - v3 entro +/- 1 punto (63 valori) | mean diff 0.020; 63 of 63 confirmed | confermata |
| P3 | rapporto SE blocchi / SE indipendenti circa 2, sopra 1,25 | ratio 0.97 to 1.04; 0 of 27 above 1.25 | smentita, spiegata |
| P4a C01 | GB risolta e rilevante; NAT e AR1 trascurabile o indeterminata | C01\|AR1: trascurabile; C01\|GB: risolta e rilevante; C01\|NAT: trascurabile | confermata |
| P4b C02 | trascurabile sui tre | C02\|AR1: trascurabile; C02\|GB: trascurabile; C02\|NAT: trascurabile | confermata |
| P4c C03 | risolta, rilevanza indeterminata sui tre | C03\|AR1: risolta, rilevanza indeterminata; C03\|GB: trascurabile; C03\|NAT: trascurabile | smentita su GB e NAT, spiegata |
| P4d C04 | uguaglianza confermata (margine 72) | MMM - FG = 135,4 [113,2; 182,1]: smentita | smentita, spiegata |
| P4e C20 | C20b rilevante, C20a trascurabile | C20a\|AR1: trascurabile; C20a\|GB: trascurabile; C20a\|NAT: trascurabile; C20b\|AR1: risolta e rilevante; C20b\|GB: risolta e rilevante; C20b\|NAT: risolta e rilevante | confermata |
| P4f C14 | vantaggio del LRT ridotto di oltre 5 punti | usabilita': scarto da 20,0 a 7,1 [6,9; 7,2]; validita': da 24,1 a 14,7 [14,5; 14,9] | confermata |

### Tabella A. Matrice affermazione x MCSE (griglia principale, U0, insieme proprio)

Stima e intervallo sono orientati nel verso del manoscritto: un valore positivo va nella direzione affermata. Intervallo percentile del bootstrap per blocchi (B = 2.000). Le righe con il nome della affermazione sono descritte in `reports/morris/S5_claims_inventory.csv`.

| sotto-affermazione | sezione | valore v3 | stima v4 | intervallo 95% | margine | esito | esito, intervallo basic | esito con bobyqa | replicati per risolvere |
|---|---|---|---|---|---|---|---|---|---|
| C01\|AR1 | Abstract; 3.1; 4.1 | 2.06/1.50/1.55 (runner-up GLMM 4.43/1.83/2.18) | 0,240 | [-0,189; 0,438] | 1,00 | trascurabile |  |  |  |
| C01\|GB | Abstract; 3.1; 4.1 | 2.06/1.50/1.55 (runner-up GLMM 4.43/1.83/2.18) | 1,492 | [1,008; 1,905] | 1,00 | risolta e rilevante |  | risolta, rilevanza indeterminata |  |
| C01\|NAT | Abstract; 3.1; 4.1 | 2.06/1.50/1.55 (runner-up GLMM 4.43/1.83/2.18) | 0,260 | [-0,196; 0,480] | 1,00 | trascurabile |  |  |  |
| C02\|AR1 | Abstract; 3.1 | 1.00/0.52/0.55 (runner-up MD 1.41/0.88/0.97) | 0,404 | [0,360; 0,430] | 1,00 | trascurabile |  |  |  |
| C02\|GB | Abstract; 3.1 | 1.00/0.52/0.55 (runner-up MD 1.41/0.88/0.97) | 0,423 | [0,399; 0,475] | 1,00 | trascurabile |  |  |  |
| C02\|NAT | Abstract; 3.1 | 1.00/0.52/0.55 (runner-up MD 1.41/0.88/0.97) | 0,253 | [0,131; 0,310] | 1,00 | trascurabile |  |  |  |
| C03\|AR1 | Abstract; 4.1 | 0.797/0.848/1.012 pp | 1,034 | [0,975; 1,092] | 1,00 | risolta, rilevanza indeterminata |  |  | 3.253 |
| C03\|GB | Abstract; 4.1 | 0.797/0.848/1.012 pp | 0,729 | [0,665; 0,792] | 1,00 | trascurabile |  |  |  |
| C03\|NAT | Abstract; 4.1 | 0.797/0.848/1.012 pp | 0,848 | [0,787; 0,910] | 1,00 | trascurabile |  |  |  |
| C03b\|AR1 | 4.1 | 0.544/0.577/0.782 pp | 0,819 | [0,754; 0,882] | 1,00 | trascurabile |  |  |  |
| C03b\|GB | 4.1 | 0.544/0.577/0.782 pp | 0,472 | [0,407; 0,540] | 1,00 | trascurabile |  |  |  |
| C03b\|NAT | 4.1 | 0.544/0.577/0.782 pp | 0,577 | [0,513; 0,640] | 1,00 | trascurabile |  |  |  |
| C04\|all | 3.1; 4.1 | 321 vs 308 | 135,417 | [113,210; 182,070] | 72,00 | smentita |  |  |  |
| C06a\|all | 3.1 | 1.80 (runner-up GEE 2.17) | 0,355 | [0,120; 0,426] | 1,00 | trascurabile |  |  |  |
| C06b\|all | 3.1 | 1.02 (runner-up MD 1.14) | 0,116 | [0,063; 0,160] | 1,00 | trascurabile |  |  |  |
| C06c\|all | 3.1 | 1.18 (runner-up MD 1.60) | 0,232 | [0,163; 0,264] | 1,00 | trascurabile |  |  |  |
| C06d\|all | 3.1 | 1.80 vs 1.74 | 0,077 | [0,032; 0,253] | 1,00 | confermata |  |  |  |
| C07\|all | Abstract; Figure 1 | Spearman -0.357 | 0,333 | [0,190; 0,548] | 0,30 | risolta, rilevanza indeterminata |  | non risolta | 37.486 |
| C08\|mad | 1; Abstract | min pairwise Spearman WC 0.976, MAD 0.833, pooled rate 0.857 | 0,905 | [0,905; 0,905] |  | confermata |  |  |  |
| C08\|madcell | 1; Abstract | min pairwise Spearman WC 0.976, MAD 0.833, pooled rate 0.857 | 0,952 | [0,857; 0,952] |  | confermata |  |  |  |
| C08\|rate | 1; Abstract | min pairwise Spearman WC 0.976, MAD 0.833, pooled rate 0.857 | 0,905 | [0,905; 0,905] |  | confermata |  |  |  |
| C08\|wc | 1; Abstract | min pairwise Spearman WC 0.976, MAD 0.833, pooled rate 0.857 | 0,976 | [0,881; 1,000] |  | confermata |  |  |  |
| C09\|all | 2.2; 4.4 | 0.217 pp | 0,196 | [0,134; 0,280] | 1,00 | confermata |  |  |  |
| C10b\|BH | Abstract; 3.2 | 23.2% (runner-up MD 28.1%) | 5,069 | [2,917; 5,972] | 5,00 | risolta, rilevanza indeterminata |  |  | 505.285 |
| C10b\|BY | Abstract; 3.2 | 23.2% (runner-up MD 28.1%) | 6,111 | [4,373; 6,667] | 5,00 | risolta, rilevanza indeterminata | risolta e rilevante |  | 1.166 |
| C10c\|BH | Abstract; 3.2 | 87.2% (next LRT 60.6%) | 33,853 | [28,480; 32,356] | 5,00 | risolta e rilevante |  |  |  |
| C10c\|BY | Abstract; 3.2 | 87.2% (next LRT 60.6%) | 34,732 | [32,196; 35,683] | 5,00 | risolta e rilevante |  |  |  |
| C11a\|AR1 | 3.3 | -2.13/-1.78/-1.84 | 0,811 | [0,781; 0,842] | 1,00 | trascurabile |  |  |  |
| C11a\|GB | 3.3 | -2.13/-1.78/-1.84 | 0,971 | [0,935; 1,007] | 1,00 | risolta, rilevanza indeterminata |  |  | 1.638 |
| C11a\|NAT | 3.3 | -2.13/-1.78/-1.84 | 0,806 | [0,773; 0,838] | 1,00 | trascurabile |  |  |  |
| C11b\|AR1 | 3.3 | +0.45/+0.51/+0.48 | 0,472 | [0,442; 0,504] | 1,00 | trascurabile |  |  |  |
| C11b\|GB | 3.3 | +0.45/+0.51/+0.48 | 0,433 | [0,396; 0,468] | 1,00 | trascurabile |  |  |  |
| C11b\|NAT | 3.3 | +0.45/+0.51/+0.48 | 0,514 | [0,485; 0,544] | 1,00 | trascurabile |  |  |  |
| C11c\|AR1 | 3.3 | max 0.69 | 0,354 | [0,332; 0,375] | 1,00 | confermata |  |  |  |
| C11c\|GB | 3.3 | max 0.69 | 0,699 | [0,670; 0,730] | 1,00 | confermata |  |  |  |
| C11c\|NAT | 3.3 | max 0.69 | 0,352 | [0,331; 0,373] | 1,00 | confermata |  |  |  |
| C11d\|GB | 3.3 | 8.93 vs 8.48 | -1,756 | [-1,826; -1,689] | 1,00 | smentita e rilevante |  | smentita |  |
| C12a\|GB | 3.5 | 25.90% | 18,518 | [17,894; 19,094] |  | risolta |  |  |  |
| C12b\|GB | 3.5 | 5.51% (z = 3.9) | 0,484 | [0,228; 0,758] | 2,50 | trascurabile |  |  |  |
| C12c\|AR1 | 3.5 | min 5.19/5.24 on the mixed mechanisms | 0,258 | [0,004; 0,361] | 2,50 | trascurabile |  |  |  |
| C12c\|GB | 3.5 | min 5.19/5.24 on the mixed mechanisms | 0,185 | [-0,055; 0,395] | 2,50 | trascurabile |  |  |  |
| C12c\|NAT | 3.5 | min 5.19/5.24 on the mixed mechanisms | 0,378 | [0,069; 0,407] | 2,50 | trascurabile |  |  |  |
| C13a\|all | 3.8; Abstract | 0.0% (3 of 144,000) | 0,002 | [0,000; 0,005] | 1,00 | confermata |  |  |  |
| C13c\|all | 3.8 | 98.8 vs 68.6 | 30,368 | [30,128; 30,633] | 5,00 | risolta e rilevante |  |  |  |
| C13d\|all | 3.8; 4.1 | 5.6 vs 1.8 | 3,771 | [3,250; 4,528] | 1,00 | risolta e rilevante |  |  |  |
| C14a\|all | Abstract; 3.9; 4.1 | 98.7 (next GEE 78.7) | 7,098 | [6,943; 7,245] | 5,00 | risolta e rilevante |  |  |  |
| C14b\|all | Abstract; 3.9; 4.1 | 89.3 (next GEE 65.2) | 14,697 | [14,505; 14,879] | 5,00 | risolta e rilevante |  |  |  |
| C14c\|all | 3.9 | 64.5 (next MMM 66.6) | 2,020 | [1,925; 2,118] | 5,00 | trascurabile |  |  |  |
| C14d\|all | 3.9 | 66.6 vs 68.7-69.0 | 2,466 | [2,391; 2,540] | 5,00 | confermata |  |  |  |
| C15a\|s1 | Abstract; 3.10; 4.1 | 2.62/1.56/2.00 (runner-up 3.06/1.95/2.47) | 0,393 | [-0,037; 0,609] | 1,00 | trascurabile |  |  |  |
| C15a\|s2 | Abstract; 3.10; 4.1 | 2.62/1.56/2.00 (runner-up 3.06/1.95/2.47) | 0,352 | [0,068; 0,501] | 1,00 | trascurabile |  |  |  |
| C15a\|s3 | Abstract; 3.10; 4.1 | 2.62/1.56/2.00 (runner-up 3.06/1.95/2.47) | 0,469 | [0,152; 0,619] | 1,00 | trascurabile |  |  |  |
| C15b\|s5 | 3.10; 4.1 | 0.61 (next FG 0.63) | 0,018 | [0,005; 0,024] | 1,00 | trascurabile |  |  |  |
| C15c\|s4 | 3.10 | 1.13 (next MD 1.32) | 0,195 | [0,174; 0,236] | 1,00 | trascurabile |  |  |  |
| C15d\|s1 | Abstract; 3.10 | -1.27/-1.27/-0.75/+1.21/+0.40 | 1,253 | [0,757; 1,393] | 1,00 | risolta, rilevanza indeterminata | risolta e rilevante |  | 1.671 |
| C15d\|s2 | Abstract; 3.10 | -1.27/-1.27/-0.75/+1.21/+0.40 | 1,247 | [0,972; 1,435] | 1,00 | risolta, rilevanza indeterminata | risolta e rilevante |  | 883 |
| C15d\|s3 | Abstract; 3.10 | -1.27/-1.27/-0.75/+1.21/+0.40 | 0,752 | [0,372; 0,749] | 1,00 | trascurabile | risolta, rilevanza indeterminata |  |  |
| C15d\|s4 | Abstract; 3.10 | -1.27/-1.27/-0.75/+1.21/+0.40 | 1,166 | [0,994; 1,188] | 1,00 | risolta, rilevanza indeterminata | risolta e rilevante |  | 372 |
| C15d\|s5 | Abstract; 3.10 | -1.27/-1.27/-0.75/+1.21/+0.40 | 0,414 | [0,333; 0,466] | 1,00 | trascurabile |  |  |  |
| C15e\|all | Abstract; 3.10 | Spearman 0.245 | 0,243 | [0,171; 0,229] | 0,10 | risolta e rilevante |  |  |  |
| C16a\|AR1 | 3.11 | 3.79/0.72/2.37 pp | 1,909 | [1,380; 2,130] | 2,00 | risolta, rilevanza indeterminata |  | trascurabile | 18.096 |
| C16a\|GB | 3.11 | 3.79/0.72/2.37 pp | 2,216 | [2,019; 2,499] | 2,00 | risolta e rilevante | risolta, rilevanza indeterminata |  |  |
| C16a\|NAT | 3.11 | 3.79/0.72/2.37 pp | 0,321 | [-0,036; 0,990] | 2,00 | trascurabile |  | smentita |  |
| C16b\|all | 3.11 | 3.79 vs 0.72 | 1,895 | [1,225; 2,369] | 2,00 | risolta, rilevanza indeterminata |  | risolta e rilevante | 30.250 |
| C16c\|GB | 3.11 | 34.57 (next shrink 32.15) | 2,423 | [2,309; 2,538] | 2,00 | risolta e rilevante |  |  |  |
| C16d\|GB | 3.11 | 22.96 vs 25.78 | 3,406 | [3,275; 3,555] | 2,00 | risolta e rilevante |  |  |  |
| C17a\|AR1 | Abstract; 3.12 | shift -0.53 to -1.15 | 0,568 | [0,544; 0,590] | 1,00 | trascurabile |  |  |  |
| C17a\|GB | Abstract; 3.12 | shift -0.53 to -1.15 | 0,653 | [0,627; 0,678] | 1,00 | trascurabile |  |  |  |
| C17a\|NAT | Abstract; 3.12 | shift -0.53 to -1.15 | 0,526 | [0,503; 0,548] | 1,00 | trascurabile |  |  |  |
| C17b\|all | Abstract; 3.12 | 1.38-3.71 | 1,389 | [1,342; 1,437] | 1,00 | risolta e rilevante |  |  |  |
| C17c\|AR1 | Abstract; 3.12; 4.1 | FG 0.80/1.08/0.95 (next 2.02/1.16/1.39) | 0,388 | [-0,012; 0,975] | 1,00 | trascurabile |  |  |  |
| C17c\|GB | Abstract; 3.12; 4.1 | FG 0.80/1.08/0.95 (next 2.02/1.16/1.39) | 0,927 | [0,491; 1,284] | 1,00 | risolta, rilevanza indeterminata |  |  | 30.347 |
| C17c\|NAT | Abstract; 3.12; 4.1 | FG 0.80/1.08/0.95 (next 2.02/1.16/1.39) | 0,085 | [-0,393; 0,651] | 1,00 | trascurabile |  |  |  |
| C17d\|AR1 | 3.12 | 1.22/1.42/1.67 | 1,453 | [1,084; 1,798] | 1,00 | risolta e rilevante |  |  |  |
| C17d\|GB | 3.12 | 1.22/1.42/1.67 | 0,927 | [0,491; 1,284] | 1,00 | risolta, rilevanza indeterminata |  |  | 30.347 |
| C17d\|NAT | 3.12 | 1.22/1.42/1.67 | 1,420 | [1,076; 1,743] | 1,00 | risolta e rilevante |  |  |  |
| C17e\|AR1 | 3.12 | MD 3.05/1.16/1.39 vs MMM 2.02/2.50/2.63 | 1,065 | [0,602; 1,393] | 1,00 | risolta, rilevanza indeterminata |  |  | 39.589 |
| C17e\|GB | 3.12 | MD 3.05/1.16/1.39 vs MMM 2.02/2.50/2.63 | 1,346 | [0,836; 1,851] | 1,00 | risolta, rilevanza indeterminata |  |  | 2.255 |
| C17e\|NAT | 3.12 | MD 3.05/1.16/1.39 vs MMM 2.02/2.50/2.63 | 1,335 | [0,878; 1,714] | 1,00 | risolta, rilevanza indeterminata |  |  | 1.573 |
| C17f\|AR1 | Abstract; 3.12; 4.1 | 0.76 vs 0.26 (GenBinaryY) | 0,491 | [0,451; 0,530] | 1,00 | trascurabile |  |  |  |
| C17f\|GB | Abstract; 3.12; 4.1 | 0.76 vs 0.26 (GenBinaryY) | 0,506 | [0,469; 0,545] | 1,00 | trascurabile |  |  |  |
| C17f\|NAT | Abstract; 3.12; 4.1 | 0.76 vs 0.26 (GenBinaryY) | 0,561 | [0,523; 0,601] | 1,00 | trascurabile |  |  |  |
| C19a\|AR1 | 3.14 | 0.409/0.369/0.392 | -0,042 | [-0,053; -0,028] | 0,10 | smentita |  |  |  |
| C19a\|GB | 3.14 | 0.409/0.369/0.392 | -0,092 | [-0,105; -0,063] | 0,10 | smentita |  |  |  |
| C19a\|NAT | 3.14 | 0.409/0.369/0.392 | -0,040 | [-0,047; -0,024] | 0,10 | smentita |  |  |  |
| C19b\|AR1 | 3.14 | factor 1.33-2.00 | 0,684 | [0,650; 0,707] | 0,10 | risolta e rilevante |  |  |  |
| C19b\|NAT | 3.14 | factor 1.33-2.00 | 0,611 | [0,585; 0,634] | 0,10 | risolta e rilevante |  |  |  |
| C19c\|AR1 | 3.14 | 2.3/2.8/3.3% | 0,029 | [0,026; 0,037] | 0,10 | confermata |  |  |  |
| C19c\|GB | 3.14 | 2.3/2.8/3.3% | 0,023 | [0,018; 0,028] | 0,10 | confermata |  |  |  |
| C19c\|NAT | 3.14 | 2.3/2.8/3.3% | 0,027 | [0,022; 0,032] | 0,10 | confermata |  |  |  |
| C20a\|AR1 | Abstract; 3.15 | 94.55/94.83/94.91 (next MD 94.24/94.59/94.64) | 0,269 | [0,245; 0,293] | 1,00 | trascurabile |  |  |  |
| C20a\|GB | Abstract; 3.15 | 94.55/94.83/94.91 (next MD 94.24/94.59/94.64) | 0,317 | [0,292; 0,341] | 1,00 | trascurabile |  |  |  |
| C20a\|NAT | Abstract; 3.15 | 94.55/94.83/94.91 (next MD 94.24/94.59/94.64) | 0,244 | [0,221; 0,268] | 1,00 | trascurabile |  |  |  |
| C20b\|AR1 | Abstract; 3.15 | 86.36/88.94/88.84 | 3,640 | [3,563; 3,721] | 1,00 | risolta e rilevante |  |  |  |
| C20b\|GB | Abstract; 3.15 | 86.36/88.94/88.84 | 3,703 | [3,595; 3,811] | 1,00 | risolta e rilevante |  |  |  |
| C20b\|NAT | Abstract; 3.15 | 86.36/88.94/88.84 | 3,422 | [3,350; 3,502] | 1,00 | risolta e rilevante |  |  |  |
| C21a\|all | Abstract; 4.1 | 2.33 (next MD 3.67) | 1,222 | [1,222; 1,333] |  | risolta |  |  |  |
| C23a\|AR1 | 4.2 | 35/45/45 | 25,000 | [25,000; 25,000] |  | vedi tabella dei punti di cambio |  |  |  |
| C23a\|GB | 4.2 | 35/45/45 | 25,000 | [25,000; 25,000] |  | vedi tabella dei punti di cambio |  |  |  |
| C23a\|NAT | 4.2 | 35/45/45 | 25,000 | [25,000; 25,000] |  | vedi tabella dei punti di cambio |  |  |  |
| C23b\|AR1 | 4.2 | 45/55/55 | 5,000 | [5,000; 5,000] |  | vedi tabella dei punti di cambio |  |  |  |
| C23b\|GB | 4.2 | 45/55/55 | 5,000 | [3,000; 5,000] |  | vedi tabella dei punti di cambio |  |  |  |
| C23b\|NAT | 4.2 | 45/55/55 | 3,000 | [3,000; 5,000] |  | vedi tabella dei punti di cambio |  |  |  |
| C23c\|AR1 | 4.2 | 32.6/32.2/32.5 | 23,703 | [23,137; 24,263] | 5,00 | risolta e rilevante |  |  |  |
| C23c\|GB | 4.2 | 32.6/32.2/32.5 | 24,985 | [24,362; 25,557] | 5,00 | risolta e rilevante |  |  |  |
| C23c\|NAT | 4.2 | 32.6/32.2/32.5 | 23,053 | [22,462; 23,635] | 5,00 | risolta e rilevante |  |  |  |
| C23d\|AR1\|N25 | 4.2 | 16-33; 10-11; 3-5 | -0,910 | [-1,214; -0,585] | 5,00 | smentita |  |  |  |
| C23d\|AR1\|N35 | 4.2 | 16-33; 10-11; 3-5 | -5,983 | [-6,237; -5,713] | 5,00 | smentita e rilevante |  |  |  |
| C23d\|GB\|N25 | 4.2 | 16-33; 10-11; 3-5 | -0,673 | [-1,008; -0,367] | 5,00 | smentita |  |  |  |
| C23d\|GB\|N35 | 4.2 | 16-33; 10-11; 3-5 | -5,430 | [-5,690; -5,170] | 5,00 | smentita e rilevante |  |  |  |
| C23d\|NAT\|N25 | 4.2 | 16-33; 10-11; 3-5 | -1,437 | [-1,757; -1,125] | 5,00 | smentita |  |  |  |
| C23d\|NAT\|N35 | 4.2 | 16-33; 10-11; 3-5 | -6,203 | [-6,463; -5,927] | 5,00 | smentita e rilevante |  |  |  |
| C23e\|AR1 | 4.2 |  | 1,383 | [0,802; 1,690] | 1,00 | risolta, rilevanza indeterminata | risolta e rilevante |  | 1.411 |
| C23e\|GB | 4.2 |  | 2,565 | [2,013; 3,056] | 1,00 | risolta e rilevante |  |  |  |
| C23e\|NAT | 4.2 |  | 1,337 | [0,862; 1,789] | 1,00 | risolta, rilevanza indeterminata |  |  | 1.978 |
| C24a\|AR1 | 4.3 | 5.42/4.06/3.93 | 3,924 | [3,642; 4,231] |  | confermata |  |  |  |
| C24a\|GB | 4.3 | 5.42/4.06/3.93 | 5,472 | [5,124; 5,806] |  | confermata |  |  |  |
| C24a\|NAT | 4.3 | 5.42/4.06/3.93 | 4,061 | [3,761; 4,363] |  | confermata |  |  |  |
| C24b\|all | 4.3 | 13.85% | 4,455 | [4,241; 4,664] |  | risolta |  |  |  |
| C24c\|AR1 | 4.3 | 5.24/4.73/4.82 | 4,805 | [4,744; 4,869] |  | confermata |  |  |  |
| C24c\|GB | 4.3 | 5.24/4.73/4.82 | 5,264 | [5,198; 5,334] |  | confermata |  |  |  |
| C24c\|NAT | 4.3 | 5.24/4.73/4.82 | 4,729 | [4,667; 4,791] |  | confermata |  |  |  |
| C25a\|AR1 | 3.1; 4.1 | <= 0.2 pp | 0,189 | [0,163; 0,445] | 1,00 | confermata |  |  |  |
| C25a\|NAT | 3.1; 4.1 | <= 0.2 pp | 0,203 | [0,176; 0,487] | 1,00 | confermata |  |  |  |
| C25b\|AR1 | 3.1 | 0.9-1.7 | 1,434 | [1,346; 1,522] | 1,00 | risolta e rilevante |  |  |  |
| C25b\|GB | 3.1 | 0.9-1.7 | 1,486 | [1,396; 1,577] | 1,00 | risolta e rilevante |  |  |  |
| C25b\|NAT | 3.1 | 0.9-1.7 | 1,173 | [1,087; 1,259] | 1,00 | risolta e rilevante |  |  |  |

### Tabella B. Esiti per classe (tutte le righe calcolate e non calcolate)

| esito | N |
|---|---|
| trascurabile | 42 |
| risolta e rilevante | 26 |
| confermata | 24 |
| descrittiva | 20 |
| risolta, rilevanza indeterminata | 17 |
| smentita |  7 |
| vedi tabella dei punti di cambio |  6 |
| smentita e rilevante |  4 |
| risolta |  3 |
| superseded |  3 |
| out_of_scope |  3 |
| not_operationalized |  1 |

### Tabella C. Punti di cambio della validita' incondizionata (C23a, C23b) e instabilita' del caso peggiore (C22)

| meccanismo | cambio v4 (N) | cambio v3 (N) | P(cambio = v3) | metodo dopo il cambio | P(metodo dopo = MD) | N distinti dell'argmax del caso peggiore di MMM |
|---|---|---|---|---|---|---|
| GB | 25 | 35 | 0,000 | gee_fg | 0,049 | 13 |
| NAT | 25 | 45 | 0,000 | gee_md | 0,973 | 14 |
| AR1 | 25 | 45 | 0,000 | gee_fg | 0,007 | 14 |

### Tabella D. Celle vinte contro la seconda (C05): punti di disegno con differenza significativa

Su 1.440 punti di disegno con almeno due metodi stimabili.

| conteggio | valore | esito |
|---|---|---|
| C05\|uncorrected | 33 | descrittiva |
| C05\|BH | 0 | confermata |
| C05\|BY | 0 | confermata |

## Appendice

### A1. Numeri del v3 non riprodotti (P1)

Riprodotti 275 di 288 numeri.

| numero | stampato | ricalcolato | fonte |
|---|---|---|---|
| C08\|mad | 0.833 | 0,929 | V-92 |
| crank\|gee_fg | 2.330 | 1,667 | V-67 |
| crank\|gee_md | 3.670 | 2,889 | V-67 |
| crank\|gee | 4.000 | 4,333 | V-67 |
| crank\|shrink | 4.330 | 5,667 | V-67 |
| crank\|glmer | 4.670 | 5,889 | V-67 |
| crank\|mmm | 4.670 | 4,000 | V-67 |
| crank\|glmer_lrt | 5.000 | 4,778 | V-67 |
| crank\|gee_kc | 7.330 | 6,778 | V-67 |
| crank25\|glmer | 2.670 | 3,778 | V-67 |
| crank25\|glmer_lrt | 2.670 | 3,000 | V-67 |
| crank25\|shrink | 3.330 | 4,778 | V-67 |
| crank25\|gee_fg | 3.670 | 3,000 | V-67 |

Pesature della deviazione media per la correlazione fra meccanismi (C08, da `O_level_all.csv`):

| pesatura | coppia | rho |
|---|---|---|
| perN | GenBinaryY vs GLMM-native | 0,929 |
| perN | GenBinaryY vs GLMM-AR1 | 0,929 |
| perN | GLMM-native vs GLMM-AR1 | 1,000 |
| perCell | GenBinaryY vs GLMM-native | 0,833 |
| perCell | GenBinaryY vs GLMM-AR1 | 0,833 |
| perCell | GLMM-native vs GLMM-AR1 | 1,000 |
| perNT | GenBinaryY vs GLMM-native | 0,881 |
| perNT | GenBinaryY vs GLMM-AR1 | 0,905 |
| perNT | GLMM-native vs GLMM-AR1 | 0,976 |

### A2. Stabilita' dal v3 al v4 (P2)

| n | media v4 - v3 | massimo abs | confermati |
|---|---|---|---|
| 63 | 0,020 | 0,220 | 63 |

### A3. Peso dei numeri casuali comuni (P3)

| quantita | SE blocchi | SE indipendente | rapporto |
|---|---|---|---|
| rate\|AR1\|gee | 0,0386 | 0,0391 | 0,989 |
| rate\|AR1\|gee_fg | 0,0328 | 0,0338 | 0,971 |
| rate\|AR1\|gee_kc | 0,0471 | 0,0462 | 1,021 |
| rate\|AR1\|gee_md | 0,0341 | 0,0344 | 0,994 |
| rate\|AR1\|glmer | 0,0391 | 0,0395 | 0,990 |
| rate\|AR1\|glmer_lrt | 0,0368 | 0,0369 | 0,996 |
| rate\|AR1\|mmm | 0,0370 | 0,0370 | 1,000 |
| rate\|AR1\|shrink | 0,0393 | 0,0402 | 0,978 |
| rate\|GB\|gee | 0,0413 | 0,0404 | 1,022 |
| rate\|GB\|gee_fg | 0,0356 | 0,0353 | 1,010 |
| rate\|GB\|gee_kc | 0,0523 | 0,0506 | 1,032 |
| rate\|GB\|gee_md | 0,0365 | 0,0359 | 1,017 |
| rate\|GB\|glmer | 0,0453 | 0,0451 | 1,006 |
| rate\|GB\|glmer_lrt | 0,0415 | 0,0411 | 1,010 |
| rate\|GB\|mmm | 0,0383 | 0,0376 | 1,019 |
| rate\|GB\|shrink | 0,0472 | 0,0460 | 1,026 |
| rate\|NAT\|gee | 0,0391 | 0,0389 | 1,007 |
| rate\|NAT\|gee_fg | 0,0330 | 0,0336 | 0,982 |
| rate\|NAT\|gee_kc | 0,0462 | 0,0455 | 1,016 |
| rate\|NAT\|gee_md | 0,0338 | 0,0341 | 0,990 |
| rate\|NAT\|glmer | 0,0395 | 0,0387 | 1,022 |
| rate\|NAT\|glmer_lrt | 0,0363 | 0,0359 | 1,010 |
| rate\|NAT\|mmm | 0,0373 | 0,0364 | 1,025 |
| rate\|NAT\|shrink | 0,0402 | 0,0394 | 1,020 |
| C03\|AR1 | 0,0310 | 0,0312 | 0,995 |
| C03\|GB | 0,0334 | 0,0332 | 1,005 |
| C03\|NAT | 0,0315 | 0,0304 | 1,035 |

Correlazione degli indicatori di rifiuto dello stesso replicato fra celle (H0, MMM e FG, scenari M1 e M3):

| relazione | coppie di celle | correlazione media |
|---|---|---|
| same N and T, mechanisms GenBinaryY~GLMM-native | 384 | 0,0021 |
| same mechanism, other N, same T | 8.640 | 0,0001 |
| same N and T, mechanisms GenBinaryY~GLMM-AR1 | 384 | -0,0000 |
| same mechanism, other N and T | 43.200 | -0,0001 |
| other mechanism, other N or T | 109.440 | -0,0001 |
| same N and T, mechanisms GLMM-AR1~GLMM-native | 384 | -0,0008 |
| same mechanism, same N, other T | 2.880 | -0,0008 |

Per confronto, la correlazione mediana delle stime fra H0 e H1 della stessa cella (S2, `V4_S2_corr_arms.csv`, celle con almeno 100 coppie) e' 0,872.

### A4. Controprova del bootstrap (G3)

| quantita | SE linearizzazione | SE bootstrap | scarto relativo |
|---|---|---|---|
| rate\|GB\|mmm | 0,0377 | 0,0383 | 0,015 |
| rate\|GB\|gee | 0,0408 | 0,0413 | 0,012 |
| rate\|GB\|gee_md | 0,0362 | 0,0365 | 0,008 |
| rate\|GB\|gee_kc | 0,0514 | 0,0523 | 0,018 |
| rate\|GB\|gee_fg | 0,0353 | 0,0356 | 0,010 |
| rate\|GB\|glmer | 0,0454 | 0,0453 | -0,001 |
| rate\|GB\|shrink | 0,0463 | 0,0472 | 0,020 |
| rate\|GB\|glmer_lrt | 0,0414 | 0,0415 | 0,001 |
| rate\|NAT\|mmm | 0,0366 | 0,0373 | 0,019 |
| rate\|NAT\|gee | 0,0386 | 0,0391 | 0,015 |
| rate\|NAT\|gee_md | 0,0340 | 0,0338 | -0,005 |
| rate\|NAT\|gee_kc | 0,0462 | 0,0462 | 0,001 |
| rate\|NAT\|gee_fg | 0,0337 | 0,0330 | -0,019 |
| rate\|NAT\|glmer | 0,0393 | 0,0395 | 0,007 |
| rate\|NAT\|shrink | 0,0399 | 0,0402 | 0,008 |
| rate\|NAT\|glmer_lrt | 0,0361 | 0,0363 | 0,004 |
| rate\|AR1\|mmm | 0,0372 | 0,0370 | -0,005 |
| rate\|AR1\|gee | 0,0395 | 0,0386 | -0,023 |
| rate\|AR1\|gee_md | 0,0349 | 0,0341 | -0,022 |
| rate\|AR1\|gee_kc | 0,0472 | 0,0471 | -0,002 |
| rate\|AR1\|gee_fg | 0,0339 | 0,0328 | -0,033 |
| rate\|AR1\|glmer | 0,0394 | 0,0391 | -0,008 |
| rate\|AR1\|shrink | 0,0401 | 0,0393 | -0,020 |
| rate\|AR1\|glmer_lrt | 0,0377 | 0,0368 | -0,025 |
| C03\|GB | 0,0331 | 0,0334 | 0,010 |
| C03\|NAT | 0,0310 | 0,0315 | 0,016 |
| C03\|AR1 | 0,0309 | 0,0310 | 0,004 |

### A5. Riproduzione di S4 (G2)

| controllo | n | scarto massimo |
|---|---|---|
| n_obs | 54.432 | 0.0e+00 |
| rej | 49.737 | 5.6e-16 |
| p_usable | 54.432 | 5.6e-16 |
| cover | 40.736 | 6.7e-16 |
| power_sa (H1 own) | 12.879 | 5.6e-16 |
| boot_se_power_raw | 24.801 | 5.0e-16 |
| mcse_power_sa (H1 own) | 12.879 | 5.6e-16 |

### A6. Statistiche non lisce: stima fuori dall'intervallo percentile

| nome | stima | media bootstrap | percentile | basic |
|---|---|---|---|---|
| C10a\|raw | 49,982 | 53,515 | [52,665; 54,400] | [45,564; 47,298] |
| C10a\|BH | 44,323 | 48,737 | [47,025; 49,754] | [38,892; 41,621] |
| C10a\|BY | 35,362 | 37,552 | [36,879; 38,211] | [32,513; 33,845] |
| C10c\|BH | 33,853 | 30,443 | [28,480; 32,356] | [35,350; 39,226] |
| C15d\|s3 | 0,752 | 0,561 | [0,372; 0,749] | [0,756; 1,133] |
| C15e\|all | 0,243 | 0,200 | [0,171; 0,229] | [0,257; 0,316] |

### A7. Celle vinte e regola sui pareggi (C04)

| metodo | v3, primo in ordine | v3, pareggi divisi | v4, primo in ordine | v4, pareggi divisi |
|---|---|---|---|---|
| mmm | 321 | 337,3 | 306 | 323,3 |
| gee | 172 | 154,7 | 183 | 166,3 |
| gee_md | 210 | 283,5 | 185 | 251,4 |
| gee_kc | 57 | 57,0 | 55 | 56,0 |
| gee_fg | 308 | 213,9 | 272 | 187,9 |
| glmer | 216 | 147,8 | 168 | 113,9 |
| shrink | 39 | 116,8 | 27 | 83,9 |
| glmer_lrt | 117 | 129,0 | 244 | 257,4 |

| versione | metodi a pari merito | punti di disegno |
|---|---|---|
| v3 | gee_fg+gee_md | 158 |
| v3 | glmer+shrink | 136 |
| v3 | gee_md+mmm |   9 |
| v3 | glmer_lrt+mmm |   9 |
| v4 | gee_fg+gee_md | 138 |
| v4 | glmer+shrink |  95 |
| v4 | glmer+glmer_lrt+shrink |  14 |
| v4 | gee_fg+gee_md+mmm |  11 |

### A8. RMSE su celle proprie e su celle comuni (C19a)

| meccanismo | metodo | RMSE, celle proprie | celle proprie | RMSE, celle comuni | celle comuni |
|---|---|---|---|---|---|
| GLMM-AR1 | gee_kc | 0,4371 | 480 | 0,3584 | 418 |
| GLMM-AR1 | gee_fg | 0,4372 | 480 | 0,3595 | 418 |
| GLMM-AR1 | gee_md | 0,4389 | 480 | 0,3606 | 418 |
| GLMM-AR1 | mmm | 0,4499 | 480 | 0,3700 | 418 |
| GLMM-AR1 | shrink | 0,3784 | 419 | 0,3743 | 418 |
| GLMM-AR1 | gee | 0,4762 | 480 | 0,3881 | 418 |
| GLMM-AR1 | glmer | 0,6233 | 418 | 0,6233 | 418 |
| GLMM-native | gee_kc | 0,4373 | 480 | 0,3542 | 417 |
| GLMM-native | gee_fg | 0,4380 | 480 | 0,3546 | 417 |
| GLMM-native | gee_md | 0,4378 | 480 | 0,3549 | 417 |
| GLMM-native | mmm | 0,4491 | 480 | 0,3648 | 417 |
| GLMM-native | shrink | 0,3689 | 418 | 0,3689 | 417 |
| GLMM-native | gee | 0,4771 | 480 | 0,3814 | 417 |
| GLMM-native | glmer | 0,5877 | 417 | 0,5877 | 417 |
| GenBinaryY | gee_kc | 0,4372 | 480 | 0,3700 | 419 |
| GenBinaryY | gee_md | 0,4384 | 480 | 0,3714 | 419 |
| GenBinaryY | gee_fg | 0,4386 | 480 | 0,3717 | 419 |
| GenBinaryY | mmm | 0,4474 | 480 | 0,3752 | 419 |
| GenBinaryY | gee | 0,4764 | 480 | 0,3999 | 419 |
| GenBinaryY | shrink | 0,4075 | 419 | 0,4075 | 419 |

### A9. Usabilita' a N = 35, v3 e v4 (C23)

| metodo | N | T | scenari | usabile v3 (%) | usabile v4 (%) |
|---|---|---|---|---|---|
| gee_fg | 35 | T<=3 | M1-M2 | 74,5 | 74,5 |
| gee_fg | 35 | T<=3 | M3-M5 | 75,6 | 75,5 |
| gee_fg | 35 | T>=4 | M1-M2 | 99,5 | 99,5 |
| gee_fg | 35 | T>=4 | M3-M5 | 99,6 | 99,6 |
| gee_md | 35 | T<=3 | M1-M2 | 74,5 | 74,5 |
| gee_md | 35 | T<=3 | M3-M5 | 75,6 | 75,5 |
| gee_md | 35 | T>=4 | M1-M2 | 99,5 | 99,5 |
| gee_md | 35 | T>=4 | M3-M5 | 99,6 | 99,6 |
| glmer_lrt | 35 | T<=3 | M1-M2 | 100,0 | 100,0 |
| glmer_lrt | 35 | T<=3 | M3-M5 | 99,9 | 33,3 |
| glmer_lrt | 35 | T>=4 | M1-M2 | 100,0 | 100,0 |
| glmer_lrt | 35 | T>=4 | M3-M5 | 100,0 | 99,9 |

### A10. Sensibilita' bobyqa: sotto-affermazioni che cambiano esito

| sotto-affermazione | stima | esito | stima bobyqa | esito bobyqa |
|---|---|---|---|---|
| C01\|GB | 1,492 | risolta e rilevante | 1,385 | risolta, rilevanza indeterminata |
| C07\|all | 0,333 | risolta, rilevanza indeterminata | 0,119 | non risolta |
| C11d\|GB | -1,756 | smentita e rilevante | -0,974 | smentita |
| C16a\|AR1 | 1,909 | risolta, rilevanza indeterminata | -0,728 | trascurabile |
| C16a\|NAT | 0,321 | trascurabile | -1,719 | smentita |
| C16b\|all | 1,895 | risolta, rilevanza indeterminata | 4,847 | risolta e rilevante |

### A11. Scenario ancorato alla coorte (M1c), descrittivo

Stesse definizioni della griglia principale su 72 punti di disegno (6 N, 4 T); caso peggiore sui 6 N.

| sotto-affermazione | stima | intervallo 95% |
|---|---|---|
| C01\|GB | -1,224 | [-2,216; -0,574] |
| C02\|GB | 0,247 | [-0,004; 0,404] |
| C01\|NAT | -1,244 | [-2,384; -0,771] |
| C02\|NAT | 0,175 | [-0,022; 0,347] |
| C01\|AR1 | -1,729 | [-2,729; -0,955] |
| C02\|AR1 | 0,259 | [-0,043; 0,345] |
| C03\|GB | 1,830 | [1,510; 2,160] |
| C03b\|GB | 1,641 | [1,330; 1,965] |
| C03\|NAT | 1,777 | [1,480; 2,084] |
| C03b\|NAT | 1,502 | [1,210; 1,796] |
| C03\|AR1 | 1,929 | [1,619; 2,230] |
| C03b\|AR1 | 1,738 | [1,441; 2,040] |
| C20a\|GB | -0,791 | [-0,923; -0,653] |
| C20b\|GB | -4,326 | [-4,742; -3,917] |
| C20a\|NAT | -0,372 | [-0,768; 0,103] |
| C20b\|NAT | -1,422 | [-1,726; -1,124] |
| C20a\|AR1 | -0,462 | [-0,728; -0,111] |
| C20b\|AR1 | -1,117 | [-1,438; -0,793] |

