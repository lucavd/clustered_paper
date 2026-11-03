
# s10tab_common.R -- S10 track "Tabelle": independent recomputation of Tables 3-5 (v7) from per-replicate rows.
# Definitions taken from the rendered Methods (manuscript/v5/MMM_v5_render.md, Sections 2.2-2.6), NOT from the
# project's aggregation code. Read-only on every project file; writes only under results_v4/S10/tabelle/.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
REPO <- "/home/user/MMM"; OUT <- file.path(REPO, "results_v4/S10/tabelle"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NCORES <- as.integer(Sys.getenv("S10TAB_CORES", "2"))
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
KD <- c("dgp", "scen", "N", "ntime", "label")
# number of mean-model coefficients p (Section 2.2 / Table 2): M1 (4), M2 adds age (5), M3, M4 as M1 (4),
# M5 adds dose (5), M1c as M1 (4), M6, M7 as M3 (4); Copula-AR1 takes the mean model of each scenario.
P_OF <- c(M1 = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M1c = 4L, M6 = 4L, M7 = 4L)
norm_scen <- function(s) ifelse(tolower(s) == "m1c", "M1c", toupper(s))
# usable Wald result (Section 2.5): fit returned without error, estimate, SE, p finite, |est| < 10, 0 < SE < 10;
# a dataset with a constant outcome is never usable; denominators = replicates with an estimable design.
wald_ok <- function(ok, est, se) { o <- (ok %in% TRUE) & is.finite(est) & is.finite(se) & abs(est) < 10 & se > 0 & se < 10; o[is.na(o)] <- FALSE; o }
zrej <- function(use, est, se, df) { z <- abs(est / se); list(rn = use & z > qnorm(0.975), rt = use & z > qt(0.975, df)) }
# convergence check of the marginalized model (Section 2.4): boundary if nlm code not in {1, 2} or a slope
# exp(a_k) < 0.001; such a fit is replaced by the fit without the time trend, kept if that refit failed.
own_cc <- function(cc) {
  cc <- as.data.table(cc)
  bnd <- !(cc$A_code %in% c(1L, 2L)) | (is.finite(cc$A_slope_min) & cc$A_slope_min < 0.001)
  # literal reading: the refit "failed" = returned an error (B_ok FALSE); variant "fin": also failed when its
  # estimate or SE is not finite (sensitivity, to locate any difference with the project's rows)
  useB <- bnd & (cc$B_ok %in% TRUE)
  useBf <- useB & is.finite(cc$B_est) & is.finite(cc$B_se)
  hasA <- if ("A_has_term" %in% names(cc)) cc$A_has_term %in% TRUE else rep(TRUE, nrow(cc))
  okA <- (cc$A_ok %in% TRUE) & hasA
  data.table(sim_id = cc$sim_id, bnd_own = bnd, used_own = ifelse(useB, "B", "A"),
             ok = ifelse(useB, TRUE, okA), est = ifelse(useB, cc$B_est, cc$A_est), se = ifelse(useB, cc$B_se, cc$A_se),
             ok_f = ifelse(useBf, TRUE, okA), est_f = ifelse(useBf, cc$B_est, cc$A_est), se_f = ifelse(useBf, cc$B_se, cc$A_se),
             used_proj = cc$used, bnd_proj = cc$bnd)
}
recs_std <- function(r) { r <- as.data.table(r)[, c(K6, "generated", "design_degenerate", "n_obs", "sum_y"), with = FALSE]
  r[, scen := norm_scen(scen)]; r[, estim := (generated %in% TRUE) & !(design_degenerate %in% TRUE)]
  r[, const := estim & (sum_y == 0 | sum_y == n_obs)]; r[] }

# per-design-point counts for study-like rows (method, term, estimate, se, p, converged)
study_counts <- function(S, R, methods, lrt = "glmer_lrt") {
  S <- as.data.table(S); S[, scen := norm_scen(scen)]
  S <- S[method %in% methods & term %in% "time:trt"]
  dup <- S[, .N, by = c("method", "sim_id")][N > 1L, .N]
  out <- list()
  for (m in methods) {
    X <- merge(R[, c("sim_id", "N", "scen", "estim", "const"), with = FALSE], S[method == m, .(sim_id, estimate, se, p, converged)], by = "sim_id", all.x = TRUE)
    base <- X$estim & !X$const
    if (m == lrt) { use <- base & (X$converged %in% TRUE) & is.finite(X$p); rn <- rt <- use & X$p < 0.05
    } else { use <- base & wald_ok(X$converged, X$estimate, X$se); rr <- zrej(use, X$estimate, X$se, X$N - P_OF[X$scen]); rn <- rr$rn; rt <- rr$rt }
    out[[m]] <- data.table(strategy = m, n = sum(use), kn = sum(rn %in% TRUE), kt = sum(rt %in% TRUE), dup = dup)
  }
  rbindlist(out)
}
cc_counts <- function(CC, R) {
  o <- own_cc(CC); X <- merge(R[, c("sim_id", "N", "scen", "estim", "const"), with = FALSE], o, by = "sim_id", all.x = TRUE)
  base <- X$estim & !X$const
  u1 <- base & wald_ok(X$ok, X$est, X$se); r1 <- zrej(u1, X$est, X$se, X$N - P_OF[X$scen])
  u2 <- base & wald_ok(X$ok_f, X$est_f, X$se_f); r2 <- zrej(u2, X$est_f, X$se_f, X$N - P_OF[X$scen])
  list(dp = rbind(data.table(strategy = "mmm_cc", n = sum(u1), kn = sum(r1$rn %in% TRUE), kt = sum(r1$rt %in% TRUE), dup = 0L),
                  data.table(strategy = "mmm_cc_f", n = sum(u2), kn = sum(r2$rn %in% TRUE), kt = sum(r2$rt %in% TRUE), dup = 0L)),
       agree = X[, .(n_cc = .N, bnd_disagree = sum(xor(bnd_own %in% TRUE, bnd_proj %in% TRUE)), used_disagree = sum(used_own != used_proj, na.rm = TRUE),
                     used_proj_na = sum(is.na(used_proj)))])
}
