# s6_mmm_boundary_sum.R -- summary of results_v4/S6/mmm_boundary_<dgp>.rds (s6_mmm_boundary.R) and
# evaluation of P5-P8 exactly as pre-registered in docs/steps/S6.md.
# Usable fit (as in the study): fit ran, estimate and model-based SE finite.
# "At the boundary": some latent-SD slope exp(a_k) < 0.001, or nlm code not in {1, 2}.
# Outputs: tables/V4_S6_mmm_bnd_cells.csv (per cell), V4_S6_mmm_bnd_pooled.csv (T >= 8 pooled),
#          V4_S6_mmm_bnd_preds.csv (P5-P8), V4_S6_mmm_bnd_repro.csv (fit A vs stored study fit)
source(path.expand("~/MMM/R/s6/s6_common.R"))
X <- rbindlist(lapply(c("GenBinaryY", "GLMM-native", "GLMM-AR1"), function(g) readRDS(file.path(REPO, sprintf("results_v4/S6/mmm_boundary_%s.rds", g)))), fill = TRUE)
z <- qnorm(0.975)
# --- reproduction of the study fit --------------------------------------------------------------
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[method == "mmm" & label == "manual" & term == "time:trt" & scen %in% c("M1", "M3", "M4", "M5") & sim_id <= 200 &
       N %in% c(10, 25, 55, 145) & ntime %in% c(2, 4, 6, 8, 10), .(dgp, model = as.integer(substr(scen, 2, 2)), N, ntime, i = sim_id, est_st = estimate, se_st = se)]
X <- merge(X, V, by = c("dgp", "model", "N", "ntime", "i"), all.x = TRUE)
REP <- X[, .(n = .N, n_stored = sum(!is.na(est_st)), n_A = sum(A_ok %in% TRUE),
             same_est = sum(abs(A_est - est_st) < 1e-8, na.rm = TRUE), same_se = sum(abs(A_se_mod - se_st) < 1e-8, na.rm = TRUE)), by = dgp]
fwrite(REP, file.path(TABDIR, "V4_S6_mmm_bnd_repro.csv")); print(REP)
# --- flags --------------------------------------------------------------------------------------
slope_min <- function(a, nm) mapply(function(a, nm) {
  if (is.na(a)) return(NA_real_); v <- as.numeric(strsplit(a, ";")[[1]]); n <- strsplit(nm, ";")[[1]]
  k <- grepl("^log\\(sigma\\):", n) & !grepl("\\(Intercept\\)", n); if (!any(k)) return(NA_real_); min(exp(v[k])) }, a, nm)
X[, useA := A_ok %in% TRUE & is.finite(A_est) & is.finite(A_se_mod)]
X[, useB := B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se_mod)]
X[, smin := slope_min(A_alpha, A_alpha_names)]
X[, nonconv := !(A_code %in% c(1L, 2L))]
X[, bnd := useA & (nonconv | (!is.na(smin) & smin < 1e-3))]
X[, rA := abs(A_est / A_se_mod) > z][, rB := abs(B_est / B_se_mod) > z][, rR := abs(A_est / A_se_rob) > z]
X[, rule_est := fifelse(bnd & useB, B_est, A_est)][, rule_se := fifelse(bnd & useB, B_se_mod, A_se_mod)]
X[, rRule := abs(rule_est / rule_se) > z]
X[, Tb := fifelse(ntime >= 8, "T>=8", "T<=6")]
saveRDS(X, file.path(REPO, "results_v4/S6/mmm_boundary_all.rds"))
pr <- function(x) 100 * mean(x, na.rm = TRUE); mc <- function(p, n) 100 * sqrt(p / 100 * (1 - p / 100) / n)
summ <- function(D) D[useA == TRUE, .(n_use = .N, pct_bnd = pr(bnd), pct_nonconv = pr(nonconv), pct_slope0 = pr(!is.na(smin) & smin < 1e-3),
          rej_all = pr(rA), n_int = sum(!bnd), rej_int = pr(rA[!bnd]), n_bnd = sum(bnd), rej_bnd = pr(rA[bnd]),
          rej_rule = pr(rRule), rej_rob = pr(rR), rej_B = pr(rB),
          empse_all = sd(A_est), rmsse_all = sqrt(mean(A_se_mod^2)), empse_int = sd(A_est[!bnd]), rmsse_int = sqrt(mean(A_se_mod[!bnd]^2)),
          empse_rule = sd(rule_est), rmsse_rule = sqrt(mean(rule_se^2)))]
C <- X[, summ(.SD), by = .(dgp, model, ntime, N)]
C[, `:=`(mcse_all = mc(rej_all, n_use), mcse_int = mc(rej_int, n_int), mcse_rule = mc(rej_rule, n_use))]
fwrite(C, file.path(TABDIR, "V4_S6_mmm_bnd_cells.csv"))
P <- X[, summ(.SD), by = .(dgp, model, Tb)]
P[, `:=`(mcse_all = mc(rej_all, n_use), mcse_int = mc(rej_int, n_int), mcse_rule = mc(rej_rule, n_use))]
fwrite(P, file.path(TABDIR, "V4_S6_mmm_bnd_pooled.csv"))
inb <- function(r) r >= 2.5 & r <= 7.5
P8 <- P[model %in% 3:5 & Tb == "T>=8"]
p5 <- P8[, all(inb(rej_int)) & all(rej_bnd > rej_int)]
p6 <- C[model %in% 3:5, .(b = pct_bnd[N == 145] >= pct_bnd[N == 10] - 2 * 100 * sqrt(0.25 / 200)), by = .(dgp, model, ntime)][, all(b)]
c1 <- C[model == 1]
p7 <- X[model == 1 & useA == TRUE, pr(nonconv)] < 2 & c1[n_use >= 100, all(inb(rej_all))]
p8 <- P8[, all(inb(rej_rule))]
PRED <- data.table(id = c("P5", "P6", "P7", "P8"), confirmed = c(p5, p6, p7, p8),
  detail = c(sprintf("non-boundary rejection %s; boundary rejection %s", paste(sprintf("%.1f", P8$rej_int), collapse = "/"), paste(sprintf("%.1f", P8$rej_bnd), collapse = "/")),
             sprintf("cells where boundary share at N=145 is not below N=10 minus 2 MCSE: %d of %d", C[model %in% 3:5, .(b = pct_bnd[N == 145] >= pct_bnd[N == 10] - 2 * 100 * sqrt(0.25 / 200)), by = .(dgp, model, ntime)][, sum(b)], C[model %in% 3:5, uniqueN(paste(dgp, model, ntime))]),
             sprintf("M1 non-convergence %.2f%%; M1 cells (n>=100) outside band: %d of %d", X[model == 1 & useA == TRUE, pr(nonconv)], c1[n_use >= 100, sum(!inb(rej_all))], c1[n_use >= 100, .N]),
             sprintf("rule rejection %s", paste(sprintf("%.1f", P8$rej_rule), collapse = "/"))))
fwrite(PRED, file.path(TABDIR, "V4_S6_mmm_bnd_preds.csv")); print(PRED)
print(P[order(model, dgp, Tb), .(dgp, model, Tb, n_use, pct_bnd, pct_nonconv, pct_slope0, rej_all, rej_int, rej_bnd, rej_rule, rej_rob, rej_B)], digits = 3)
