# s7e_glmm_pseudotrue.R -- S7e explanation (after the run, not a prediction): the value of time:trt that the study's GLMM
# converges to on copula data under H0 (marginal time:trt = 0). One dataset of 20,000 subjects per (scenario, T), H0, seed 20261024;
# glmer with the study's random structure (Laplace, bobyqa). A non-zero limit makes the Wald rejection rate grow with N.
source(path.expand("~/MMM/R/s7e/s7e_common.R"))
env <- load_env_cop()
G <- CJ(model = c(1L, 3L), ntime = c(2L, 4L, 10L))[!(model == 3L & ntime == 2L)]
res <- rbindlist(lapply(seq_len(nrow(G)), function(k) {
  m <- G$model[k]; Tt <- G$ntime[k]; RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(20261024L)
  d <- as.data.frame(env$make_all_data(20000L, Tt, 1L, m, "manual")); d$id <- factor(d$id)
  f <- if (m == 1L) binY ~ time * trt + (1 | id) else binY ~ time * trt + (1 + time | id)
  t0 <- Sys.time(); g <- suppressWarnings(lme4::glmer(f, data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa")))
  s <- summary(g)$coefficients
  data.table(model = m, ntime = Tt, n_obs = nrow(d), est_timetrt = s["time:trt", 1], se_timetrt = s["time:trt", 2],
             re_var_intercept = as.numeric(lme4::VarCorr(g)$id[1, 1]), secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}))
print(res, digits = 4); fwrite(res, file.path(TABDIR, "V4_S7e_glmm_pseudotrue.csv"))
