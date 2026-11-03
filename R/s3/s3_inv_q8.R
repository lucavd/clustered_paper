# S3, exploratory investigation of Q8b (not a pre-registered expectation).
# For sampled replicates, regenerate the data (v4 code, seed, RNGkind), refit the
# Wald glmer exactly as analysis_model*.R (default control) and the LRT full model
# exactly as analyze_glmer_lrt (bobyqa, maxfun 1e5), and compare: number of
# variance parameters (D-08 check), logLik, singularity, |d time:trt|.
source(path.expand("~/MMM/R/s3/s3_common.R")); library(data.table)
set.seed(20261017)
tk <- as.data.table(s3_tasks())
G <- list(list("GenBinaryY", "M5", "L"), list("GLMM-AR1", "M5", "F2"), list("GLMM-native", "M5", "L"),
          list("GenBinaryY", "M3", "L"), list("GenBinaryY", "M1c", "C"))
pairs_of <- function(t) {
  x <- as.data.table(readRDS(t$file)); w <- if (t$block == "L") as.data.table(readRDS(t$orig_file)) else x
  lr <- x[method == "glmer_lrt" & term == "time:trt"]; wa <- w[method == "glmer" & term == "time:trt"]
  m <- merge(lr[, .(sim_id = as.integer(sim_id), e_lrt = estimate, s_lrt = status, c_lrt = converged)],
             wa[, .(sim_id = as.integer(sim_id), e_w = estimate, s_w = status, c_w = converged)], by = "sim_id")
  m <- m[c_lrt & c_w & is.na(s_lrt) & is.na(s_w)]
  m[, `:=`(dgp = t$dgp, scen = t$scen, model = t$model, N = t$N, ntime = t$ntime, label = t$label, close = abs(e_lrt - e_w) < 1e-3)]
}
S <- rbindlist(lapply(G, function(g) {
  tt <- tk[dgp == g[[1]] & scen == g[[2]] & block == g[[3]]]; tt <- tt[sample(.N, 60)]
  P <- rbindlist(lapply(seq_len(nrow(tt)), function(k) pairs_of(tt[k])))
  rbind(P[close == FALSE][sample(.N, min(.N, 60))], P[close == TRUE][sample(.N, min(.N, 20))])
}))
cat("sampled replicates:", nrow(S), "\n")
re_of <- function(model) c("1" = "(1 | id)", "3" = "(1 + time | id)", "4" = "(1 + time | id)", "5" = "(1 + time + dose | id)")[[as.character(model)]]
rhs_of <- function(model) if (model == 5) "time*trt + dose" else "time*trt"
R <- rbindlist(parallel::mclapply(seq_len(nrow(S)), function(k) {
  s <- S[k]; env <- load_env(s$dgp, RUNDIR_V4); set_rep_seed(s$model, s$sim_id)
  d <- gen_data(env, s$scen, s$model, s$N, s$ntime, s$sim_id, s$label); d <- as.data.frame(d); d$id <- factor(d$id)
  f <- as.formula(paste("binY ~", rhs_of(s$model), "+", re_of(s$model)))
  w <- suppressWarnings(suppressMessages(lme4::glmer(f, data = d, family = binomial)))
  l <- suppressWarnings(suppressMessages(lme4::glmer(f, data = d, family = binomial,
         control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))))
  data.table(s[, .(dgp, scen, N, ntime, label, sim_id, close)],
    reproduces_wald = abs(lme4::fixef(w)[["time:trt"]] - s$e_w) < 1e-8, reproduces_lrt = abs(lme4::fixef(l)[["time:trt"]] - s$e_lrt) < 1e-8,
    n_theta_w = length(lme4::getME(w, "theta")), n_theta_l = length(lme4::getME(l, "theta")),
    ll_w = as.numeric(logLik(w)), ll_l = as.numeric(logLik(l)), sing_w = lme4::isSingular(w), sing_l = lme4::isSingular(l),
    d_est = abs(lme4::fixef(w)[["time:trt"]] - lme4::fixef(l)[["time:trt"]]))
}, mc.cores = 40, mc.preschedule = FALSE))
R[, d_ll := ll_l - ll_w]
fwrite(R, file.path(REPO, "tables/V4_S3_q8_investigation.csv"))
options(width = 220)
print(R[, .(n = .N, reproduces_w = sum(reproduces_wald), reproduces_l = sum(reproduces_lrt), same_n_theta = sum(n_theta_w == n_theta_l),
            any_singular = sum(sing_w | sing_l), d_ll_abs_lt_1e3 = sum(abs(d_ll) < 1e-3), lrt_higher_ll = sum(d_ll > 1e-3), wald_higher_ll = sum(d_ll < -1e-3),
            median_d_est = median(d_est)), by = .(dgp, scen, close)])
