# =============================================================================
# s2_g_lus.R -- S2(g), Q7 (D-10, D-10b): compare the LUS cohort with the DGM by
# analogues. Data read-only: ARTICOLO 2026/LURID_bis/lurid_long.rds (built from
# data_raw/LURID.xlsx by 01_prepare_data.R). Only aggregates are written.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
d0 <- readRDS(file.path(REPO, "ARTICOLO 2026/LURID_bis/lurid_long.rds"))
# patient-timepoint blocks of ten zones, as rev2/T2_defs.R
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along)
d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
d0$pt <- factor(paste(d0$id, d0$tp, sep = "_"))
d <- d0 |> filter(!is.na(lus_severe)) |>
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior",
                                   grepl("P$", as.character(section)) ~ "posterior", TRUE ~ "lateral"),
                         levels = c("anterior", "lateral", "posterior")),
         id = droplevels(factor(id)), pt = droplevels(factor(pt)), hemi = factor(paste(pt, side, sep = "_"))) |>
  arrange(id) |> as.data.frame()
OUT <- list(); add <- function(item, value, se = NA_real_, lo = NA_real_, hi = NA_real_, dgm = NA_real_, note = "")
  OUT[[length(OUT) + 1L]] <<- data.table(item, value, se, lo, hi, dgm, note)
add("n_patients", n_distinct(d$id)); add("n_patient_timepoints", nlevels(d$pt)); add("n_zone_outcomes", nrow(d))
cat(sprintf("Q7.1 rows %d (expected 2296), patients %d, patient-timepoints %d\n", nrow(d), n_distinct(d$id), nlevels(d$pt)))

ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
mA <- glmer(lus_severe ~ region + (1 | id), family = binomial, data = d, control = ctl)
mB <- glmer(lus_severe ~ region + (1 | id) + (1 | pt), family = binomial, data = d, control = ctl)
mC <- glmer(lus_severe ~ region + (1 | id) + (1 | pt) + (1 | hemi), family = binomial, data = d, control = ctl)
vc <- function(m) { v <- as.data.frame(VarCorr(m)); setNames(v$vcov, v$grp) }
vA <- vc(mA); vB <- vc(mB); vC <- vc(mC)
cat(sprintf("Q7.1 patient variance, patient level only: %.4f (expected 3.545 +/- 0.001: %s)\n", vA["id"], abs(vA["id"] - 3.545) <= 0.001))
print(list(A = vA, B = vB, C = vC))
# profile CIs on the SD scale, from the glmer deviance function (lme4's
# confint(profile) fails here: a lower deviance than bobyqa's in A, boundary
# components in B/C). pars = c(theta, beta); theta = SD of each intercept term.
prof_ci <- function(m, f) {
  dv0 <- glmer(f, family = binomial, data = d, control = ctl, devFunOnly = TRUE)
  th <- getME(m, "theta"); be <- fixef(m); nt <- length(th)
  dv <- function(x) { r <- tryCatch(dv0(x), error = function(e) NA_real_)
                      if (!is.finite(r)) { try(dv0(c(th, be)), silent = TRUE); r <- 1e10 }; r }
  lo <- c(rep(0, nt), rep(-Inf, length(be)))
  opt <- optim(c(th, be), dv, method = "L-BFGS-B", lower = lo, control = list(maxit = 5000, factr = 1e5))
  d_glmer <- -2 * as.numeric(logLik(m))   # NB deviance() of a GLMM is the conditional deviance, not -2 logLik
  dmin <- min(opt$value, d_glmer); start <- if (opt$value < d_glmer) opt$par else c(th, be)
  pdev <- function(k, s) {               # profile deviance with theta_k fixed at s
    g <- function(q) { x <- numeric(length(start)); x[k] <- s; x[-k] <- q; dv(x) }
    o <- optim(start[-k], g, method = "L-BFGS-B", lower = lo[-k], control = list(maxit = 5000, factr = 1e5))
    o$value
  }
  crit <- qchisq(0.95, 1)
  rbindlist(lapply(seq_len(nt), function(k) {
    sk <- start[k]; f0 <- function(s) pdev(k, s) - dmin - crit
    lwr <- if (sk < 1e-3 || f0(0) <= 0) 0 else uniroot(f0, c(0, sk))$root
    upr <- uniroot(f0, c(max(sk, 1e-3), max(5 * sk, 5)), extendInt = "upX")$root
    data.table(grp = sub("[.].*$", "", names(th)[k]), sd_hat = sk, sd_lo = lwr, sd_hi = upr,
               lr_at_sd1 = pdev(k, 1) - dmin, dev_gain_vs_glmer = d_glmer - dmin)
  }))
}
ciA <- prof_ci(mA, lus_severe ~ region + (1 | id))
ciB <- prof_ci(mB, lus_severe ~ region + (1 | id) + (1 | pt))
ciC <- prof_ci(mC, lus_severe ~ region + (1 | id) + (1 | pt) + (1 | hemi))
print(list(ciA = ciA, ciB = ciB, ciC = ciC))
lat <- pi^2 / 3
icc <- function(v) v / (sum(v) + lat)
put_var <- function(tag, v, ci) for (g in names(v)) {
  r <- ci[grp == g]
  lo <- if (nrow(r)) r$sd_lo^2 else NA_real_; hi <- if (nrow(r)) r$sd_hi^2 else NA_real_
  lr1 <- if (nrow(r)) r$lr_at_sd1 else NA_real_
  add(sprintf("var_%s_model%s", g, tag), unname(v[g]), lo = lo, hi = hi, dgm = 1,
      note = sprintf("%s; profile LR at variance 1 = %.2f", if (is.na(lo)) "profile CI unavailable" else if (1 >= lo && 1 <= hi) "coerente" else "scarto dichiarato", lr1))
  add(sprintf("icc_latent_%s_model%s", g, tag), unname(v[g]) / (sum(v) + lat), dgm = 1 / (1 + lat))
}
put_var("A", vA, ciA); put_var("B", vB, ciB); put_var("C", vC, ciC)
add("LRT_pt_level_B_vs_A_chi2", as.numeric(2 * (logLik(mB) - logLik(mA))), note = "boundary mixture p = 0.5 * pchisq(chi2, 1, lower = FALSE)")

# prevalence with patient-cluster bootstrap
prev <- mean(d$lus_severe)
ids <- levels(d$id); sp <- split(d$lus_severe, d$id)
set.seed(20261015); RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261015)
bt <- replicate(2000, { s <- sample(ids, replace = TRUE); x <- unlist(sp[s], use.names = FALSE); mean(x) })
ci_prev <- quantile(bt, c(0.025, 0.975))
add("prevalence_lus_severe", prev, se = sd(bt), lo = ci_prev[1], hi = ci_prev[2], dgm = plogis(-2),
    note = if (plogis(-2) >= ci_prev[1] && plogis(-2) <= ci_prev[2]) "coerente" else "scarto dichiarato")
add("logit_prevalence", qlogis(prev), lo = qlogis(ci_prev[1]), hi = qlogis(ci_prev[2]), dgm = -2)
cat(sprintf("Q7.2 prevalence %.4f [%.4f, %.4f] vs DGM %.4f\n", prev, ci_prev[1], ci_prev[2], plogis(-2)))
# DGM marginal prevalence over the realised design (GenBinaryY M1 marginal mean), for orientation
for (T in c(2, 10)) {
  g <- CJ(time = 0:(T - 1), trt = 0:1)[, w := fifelse(time == 0, 1, 0.7 + 0.3 * (T - 1 - time) / (T - 1))]
  add(sprintf("dgm_marginal_prevalence_M1_H1_T%d", T), g[, sum(w * plogis(-2 + 0.1 * time + 0.2 * trt + 0.2 * time * trt)) / sum(w)])
}

# effect-size references (D-10b): full-cohort contrasts, manuscript mean model
f <- lus_severe ~ region + disease + invasive + age_weeks
mm_fit <- mm(f, t.formula = ~ 1, data = d, id = id, verbose = FALSE)
V <- if (!is.null(mm_fit$rob.cov)) mm_fit$rob.cov else mm_fit$mod.cov
gee_fit <- geeglm(f, family = binomial, id = id, corstr = "exchangeable", data = d)
glmm_fit <- glmer(update(f, . ~ . + (1 | id)), family = binomial, data = d, control = ctl)
sg <- summary(gee_fit)$coefficients; sq <- summary(glmm_fit)$coefficients
terms <- setdiff(names(mm_fit$beta), "(Intercept)")
terms <- terms[terms %in% rownames(sg)]
for (tm in terms) {
  add(sprintf("contrast_MMM_%s", tm), unname(mm_fit$beta[tm]), se = sqrt(V[tm, tm]), note = "marginal")
  add(sprintf("contrast_GEE_%s", tm), sg[tm, 1], se = sg[tm, 2], note = "marginal")
  if (tm %in% rownames(sq)) add(sprintf("contrast_GLMM_%s", tm), sq[tm, 1], se = sq[tm, 2], note = "conditional")
}
cat(sprintf("Q7.1 posterior vs anterior: MMM %.4f (SE %.4f), expected 0.378 (0.168): %s; GEE %.4f (%.4f); GLMM %.4f (%.4f)\n",
            mm_fit$beta["regionposterior"], sqrt(V["regionposterior", "regionposterior"]),
            abs(mm_fit$beta["regionposterior"] - 0.378) <= 0.001 && abs(sqrt(V["regionposterior", "regionposterior"]) - 0.168) <= 0.001,
            sg["regionposterior", 1], sg["regionposterior", 2], sq["regionposterior", 1], sq["regionposterior", 2]))
for (T in T_LEVELS) add(sprintf("dgm_cumulative_interaction_T%d", T), 0.2 * (T - 1), note = "0.2 per occasion x (T - 1); no coherence criterion (D-10b)")
R <- rbindlist(OUT)
fwrite(R, file.path(TABDIR, "V4_S2_lus_vs_dgm.csv"))
options(width = 200); print(R, digits = 4)
