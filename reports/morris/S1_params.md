| DGM | M | mean model (coefficient order) | coefficients H1 | dependence (H1) | H0 differs from H1 in |
|---|---|---|---|---|---|
| GenBinaryY | 1 | ~time * trt ((Intercept), time, trt, time:trt) | marginal -2, 0.1, 0.2, 0.2 | lv ~1, sigma (1); t ~1, gamma (1); q = 20 | tested coefficient only |
| GenBinaryY | 2 | ~time * trt + age ((Intercept), time, trt, age, time:trt) | marginal -2, 0.1, 0.2, 0.01, 0.2 | lv ~1, sigma (1); t ~1, gamma (1); q = 20 | tested coefficient only |
| GenBinaryY | 3 | ~time * trt ((Intercept), time, trt, time:trt) | marginal -2, 0.1, 0.2, 0.2 | lv ~time, sigma (1, 0.01); t ~1, gamma (1); q = 20 | time:trt and sigma (H0: sigma (1, 0), gamma (1)) |
| GenBinaryY | 4 | ~time * trt ((Intercept), time, trt, time:trt) | marginal -2, 0.2, 0.2, 0.05 | lv ~time, sigma (1, 0.01); t ~trt, gamma (1, 0.01); q = 20 | time:trt and sigma, gamma (H0: sigma (1, 0), gamma (1, 0)) |
| GenBinaryY | 5 | ~time * trt + dose ((Intercept), time, trt, dose, time:trt) | marginal -2, 0.1, 0.2, 0.05, 0.2 | lv ~time + dose, sigma (1, 0.01, 0.01); t ~trt + steroids, gamma (1, 0.1, 0.05); q = 20 | time:trt and sigma, gamma (H0: sigma (1, 0, 0), gamma (1, 0, 0)) |
| GLMM-native | 1 | ~time * trt ((Intercept), time, trt, time:trt) | conditional -2.32, 0.116, 0.232, 0.232 | b_i on ~1, D = [1] | tested coefficient only |
| GLMM-native | 2 | ~time * trt + age ((Intercept), time, trt, age, time:trt) | conditional -2.32, 0.116, 0.232, 0.012, 0.232 | b_i on ~1, D = [1] | tested coefficient only |
| GLMM-native | 3 | ~time * trt ((Intercept), time, trt, time:trt) | conditional -2.32, 0.116, 0.232, 0.232 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001] | tested coefficient only |
| GLMM-native | 4 | ~time * trt ((Intercept), time, trt, time:trt) | conditional -2.32, 0.232, 0.232, 0.058 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001] | tested coefficient only |
| GLMM-native | 5 | ~time * trt + dose ((Intercept), time, trt, dose, time:trt) | conditional -2.32, 0.116, 0.232, 0.058, 0.232 | b_i on ~1 + time + dose, D = [1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.001] | tested coefficient only |
| GLMM-AR1 | 1 | ~time * trt ((Intercept), time, trt, time:trt) | marginal target -2, 0.1, 0.2, 0.2 | b_i on ~1, D = [1]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 2 | ~time * trt + age ((Intercept), time, trt, age, time:trt) | marginal target -2, 0.1, 0.2, 0.01, 0.2 | b_i on ~1, D = [1]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 3 | ~time * trt ((Intercept), time, trt, time:trt) | marginal target -2, 0.1, 0.2, 0.2 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 4 | ~time * trt ((Intercept), time, trt, time:trt) | marginal target -2, 0.2, 0.2, 0.05 | b_i on ~1 + time, D = [1, 0.01; 0.01, 0.001]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |
| GLMM-AR1 | 5 | ~time * trt + dose ((Intercept), time, trt, dose, time:trt) | marginal target -2, 0.1, 0.2, 0.05, 0.2 | b_i on ~1 + time + dose, D = [1, 0.01, 0.01; 0.01, 0.001, 0; 0.01, 0, 0.01]; AR(1) rho = 0.6, sigma2_ar = 0.5 | tested coefficient only |

Design draws, in execution order (table rows `draw_*`). They are the same in the three mechanisms and in both arms (checked); the scenarios differ only in how many covariates they draw:

- M1: nclust, trt
- M2: nclust, trt, age
- M3: nclust, trt, age
- M4: nclust, trt, age
- M5: nclust, trt, age, sex, steroids, dose

Full list (M5):

- follow-up: `nclust = ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))` (p_full_followup = 0.7)
- `nclust` <- `runif(n = N)` (level=subject; standardised=FALSE; assigned_as=nclust)
- `nclust` <- `sample(x = seq_len(ntime - 1), size = N, replace = TRUE)` (level=subject; standardised=FALSE; assigned_as=nclust)
- `trt` <- `rbinom(n = N, size = 1, prob = 0.5)` (level=subject; standardised=FALSE; assigned_as=trt)
- `age` <- `rnorm(n = N, mean = 50, sd = 10)` (level=subject; standardised=TRUE; assigned_as=age_raw)
- `sex` <- `rbinom(n = N, size = 1, prob = 0.5)` (level=subject; standardised=FALSE; assigned_as=sex)
- `steroids` <- `sample(x = 1:3, size = N, replace = TRUE)` (level=subject; standardised=FALSE; assigned_as=steroids)
- `dose` <- `rnorm(n = length(id))` (level=observation; standardised=TRUE; assigned_as=dose)

Covariates drawn but entering no part of the DGM (they consume random numbers):

- GenBinaryY M3: age
- GenBinaryY M4: age
- GenBinaryY M5: age, sex
- GLMM-native M3: age
- GLMM-native M4: age
- GLMM-native M5: age, sex, steroids
- GLMM-AR1 M3: age
- GLMM-AR1 M4: age
- GLMM-AR1 M5: age, sex, steroids

Seed: `set.seed(20260805L + 1000L * model + i)` per replicate; RNGkind in the workers: `L'Ecuyer-CMRG, Inversion, Rejection` (set by crew 1.3.2 / mirai 2.7.2 in every worker; NOT set by run_full.R (plain R session: Mersenne-Twister, Inversion, Rejection)).
