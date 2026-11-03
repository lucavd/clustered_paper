# =============================================================================
# s4_checks.R -- S4 expectations Q3, Q4, Q5, Q7, Q8, Q9 (Q0 in s4_q0.sh,
# Q1-Q2 in agg_v4.R, Q6 in s4_truth.R). Reads files only.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
suppressMessages(library(rsimsum))
P  <- fread(file.path(TABDIR, "V4_perf_all.csv.gz"))
untouched_cell <- function(dgp, scen, label) !((dgp == "GenBinaryY" & scen %in% c("M3", "M4", "M5") & label == "manual") | (dgp == "GLMM-AR1" & scen == "M5"))
P[, model := as.integer(model)]

# ---- Q3: reproduction of O_level_all / O_coverage_all ---------------------------
OL <- fread(file.path(TABDIR, "O_level_all.csv")); OC <- fread(file.path(TABDIR, "O_coverage_all.csv"))
OL[, scen := paste0("M", model)]; OC[, scen := paste0("M", model)]
OL[, use := model %in% 1:2 | (set == "own" & method != LRT & untouched_cell(dgp, scen, "manual"))]
OC[, use := (set == "own" & !(dgp == "GLMM-AR1" & model == 5)) | (set == "intersection" & model %in% 1:2)]
PL <- P[arm == "H0" & usab == "U0" & scen != "M1c" & method %in% PRIM, .(dgp, scen, N, ntime, method, set, n_obs, rej)]
ML <- merge(OL[use == TRUE], PL, by = c("dgp", "scen", "N", "ntime", "method", "set"), all.x = TRUE)
q3l <- ML[, .(rows = .N, n_equal = sum(n == n_obs, na.rm = TRUE), max_rate_diff = max(abs(k_rej / n - rej), na.rm = TRUE), missing = sum(is.na(n_obs)))]
# reverse direction: untouched v4 rows with n_obs > 0 must exist in O_level_all
PLr <- PL[n_obs > 0 & (scen %in% c("M1", "M2") | (set == "own" & method != LRT & untouched_cell(dgp, scen, "manual")))]
q3l_rev <- PLr[!OL, on = c("dgp", "scen", "N", "ntime", "method", "set"), .N]
PC <- P[arm == "H1" & usab == "U0" & scen != "M1c" & method %in% WALD, .(dgp, scen, N, ntime, method, set, n_obs, cover_v3marg)]
MC <- merge(OC[use == TRUE], PC, by = c("dgp", "scen", "N", "ntime", "method", "set"), all.x = TRUE)
q3c <- MC[, .(rows = .N, n_equal = sum(n == n_obs, na.rm = TRUE), max_cov_diff = max(abs(cov_gen / 100 - cover_v3marg), na.rm = TRUE), missing = sum(is.na(n_obs)))]
Q3 <- rbind(q3l[, .(table = "O_level_all (H0)", rows, n_equal, max_diff = max_rate_diff, missing, reverse_missing = q3l_rev)],
            q3c[, .(table = "O_coverage_all (H1)", rows, n_equal, max_diff = max_cov_diff, missing, reverse_missing = NA_integer_)])
Q3[, pass := n_equal == rows & missing == 0 & max_diff <= 1e-12 & (is.na(reverse_missing) | reverse_missing == 0)]
fwrite(Q3, file.path(TABDIR, "V4_S4_q3_reproduction.csv")); cat("Q3\n"); print(Q3)

# ---- Q4: usability counts against S2 --------------------------------------------
US <- fread(file.path(TABDIR, "V4_S2_usable_by_cell.csv")); US[, scen := paste0("M", model)]
US <- US[untouched_cell(dgp, scen, label) & !(method == LRT & model %in% 3:5)]
PU <- dcast(P[set == "own" & scen != "M1c", .(dgp, scen, N, ntime, label, method, usab, n_obs)], dgp + scen + N + ntime + label + method ~ usab, value.var = "n_obs")
MU <- merge(US, PU, by = c("dgp", "scen", "N", "ntime", "label", "method"), all.x = TRUE)
Q4 <- MU[, .(rows = .N, U0_equal = sum(n_U0 == U0, na.rm = TRUE), U1_equal = sum(n_U1 == U1, na.rm = TRUE), U2_equal = sum(n_U2 == U2, na.rm = TRUE), missing = sum(is.na(U0)))]
Q4[, pass := U0_equal == rows & U1_equal == rows & U2_equal == rows & missing == 0]
fwrite(Q4, file.path(TABDIR, "V4_S4_q4_usable_vs_S2.csv")); cat("Q4\n"); print(Q4)

# ---- Q5: rsimsum on 300 random combinations, independent selection ----------------
V  <- readRDS(file.path(OUT_S4, "rows_v4.rds"))[!is.na(term)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
V <- V[!DEG, on = c("model", "N", "ntime", "sim_id")]
TR <- readRDS(file.path(OUT_S4, "truth_ar1_rep.rds")); TT <- fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))
set.seed(20261017L)
cand <- P[n_obs >= 10 & method != LRT & !truth_type %in% c("none", "none (D-18)")]
S5 <- cand[sample(.N, 300L)]
q5_one <- function(k) {
  s <- S5[k]; uc <- s$usab
  d <- V[dgp == s$dgp & scen == s$scen & N == s$N & ntime == s$ntime & label == s$label & method == s$method & get(uc) == TRUE]
  if (s$set == "intersection") {
    ids <- V[dgp == s$dgp & scen == s$scen & N == s$N & ntime == s$ntime & label == s$label & method %in% PRIM & get(uc) == TRUE,
             .(k = uniqueN(method)), by = sim_id][k == length(PRIM), sim_id]
    d <- d[sim_id %in% ids]
  }
  tt <- TT[dgp == s$dgp & scen == s$scen & ntime == s$ntime & label == s$label]
  if (s$truth_type == "marginal") d[, th := tt$truth_marginal]
  else if (s$truth_type == "conditional") d[, th := tt$truth_conditional]
  else d[TR[dgp == s$dgp & scen == s$scen & N == s$N & ntime == s$ntime & label == s$label], on = "sim_id", th := i.truth_cond]
  d <- as.data.frame(d[, .(estimate, se, th)])
  const <- length(unique(d$th)) == 1
  ss <- simsum(data = d, estvarname = "estimate", se = "se", true = if (const) d$th[1] else "th", control = list(mcse = TRUE))
  sm <- as.data.table(ss$summ)
  g <- function(st, col = "est") sm[stat == st][[col]]
  ref <- c(n = g("nsim"), bias = g("bias"), empse = g("empse"), modse = g("modelse"), relerr = g("relerror"), mse = g("mse"),
           cover = g("cover"), becover = g("becover"), power = g("power"),
           mcse_bias = g("bias", "mcse"), mcse_empse = g("empse", "mcse"), mcse_modse = g("modelse", "mcse"), mcse_relerr = g("relerror", "mcse"),
           mcse_mse = g("mse", "mcse"), mcse_cover = g("cover", "mcse"), mcse_becover = g("becover", "mcse"), mcse_power = g("power", "mcse"))
  mine <- c(n = s$n_obs, bias = s$bias, empse = s$empse, modse = s$modse, relerr = s$relerr_modse, mse = s$mse,
            cover = s$cover, becover = s$becover, power = s$rej, mcse_bias = s$mcse_bias, mcse_empse = s$mcse_empse,
            mcse_modse = s$mcse_modse, mcse_relerr = s$mcse_relerr, mcse_mse = s$mcse_mse, mcse_cover = s$mcse_cover,
            mcse_becover = s$mcse_becover, mcse_power = s$mcse_rej)
  dd <- abs(ref - mine); if (!const) dd["mcse_bias"] <- NA
  data.table(k = k, dgp = s$dgp, scen = s$scen, N = s$N, ntime = s$ntime, label = s$label, method = s$method, usab = uc, set = s$set,
             truth_type = s$truth_type, const_truth = const, measure = names(dd), abs_diff = as.numeric(dd),
             ref = as.numeric(ref), mine = as.numeric(mine))
}
Q5D <- rbindlist(parallel::mclapply(seq_len(nrow(S5)), q5_one, mc.cores = 30))
fwrite(Q5D, file.path(TABDIR, "V4_S4_q5_rsimsum_detail.csv.gz"))
Q5 <- Q5D[, .(combos = uniqueN(k), compared = sum(!is.na(abs_diff)), max_abs_diff = max(abs_diff, na.rm = TRUE),
              over_1e10 = sum(abs_diff > 1e-10, na.rm = TRUE)), by = measure]
cat("Q5 rsimsum\n"); print(Q5)
# identity MSE = bias^2 + empSE^2 (n - 1)/n
idn <- P[n_obs >= 2 & is.finite(mse), .(truth_type, d = abs(mse - (bias^2 + empse^2 * (n_obs - 1) / n_obs)))]
Q5I <- idn[, .(combos = .N, max_abs_diff = max(d), over_1e12 = sum(d > 1e-12)), by = truth_type]
cat("Q5 identity\n"); print(Q5I)
fwrite(rbind(Q5[, .(check = paste("rsimsum", measure), combos, compared, max_abs_diff, over_tol = over_1e10, tol = 1e-10)],
             Q5I[, .(check = paste("identity MSE,", truth_type), combos, compared = combos, max_abs_diff, over_tol = over_1e12, tol = 1e-12)]),
       file.path(TABDIR, "V4_S4_q5_formulas.csv"))

# ---- Q7: bootstrap SE of raw power vs binomial MCSE ---------------------------------
q7 <- P[arm == "H1" & rej >= 0.05 & rej <= 0.95 & n_obs >= 100 & is.finite(boot_se_power_raw)]
q7[, ratio := boot_se_power_raw / mcse_rej]
Q7 <- q7[, .(combos = .N, within_10pct = sum(abs(ratio - 1) <= 0.10), share = mean(abs(ratio - 1) <= 0.10),
             median_ratio = median(ratio), q01 = quantile(ratio, 0.01), q99 = quantile(ratio, 0.99))]
Q7[, pass := share >= 0.95]
fwrite(Q7, file.path(TABDIR, "V4_S4_q7_bootstrap.csv")); cat("Q7\n"); print(Q7)

# ---- Q8: effect of D-08 on the LRT (paired v4 - v3) -----------------------------------
LV3 <- readRDS(file.path(OUT_S4, "lrt_v3_m3m5.rds"))[!is.na(term)]
LV3[, model := as.integer(model)]; LV3 <- LV3[!DEG, on = c("model", "N", "ntime", "sim_id")]
L4 <- V[method == LRT & model %in% 3:5]
KE <- unique(rbind(LV3[, .(dgp, scen, N, ntime, label, model)], L4[, .(dgp, scen, N, ntime, label, model)]))
KE <- KE[, .(sim_id = 1:1000), by = .(dgp, scen, N, ntime, label, model)][!DEG, on = c("model", "N", "ntime", "sim_id")]
KE[LV3, on = c("dgp", "scen", "N", "ntime", "label", "sim_id"), `:=`(u_old = i.U0, p_old = i.p)]
KE[L4, on = c("dgp", "scen", "N", "ntime", "label", "sim_id"), `:=`(u_new = i.U0, p_new = i.p)]
KE[is.na(u_old), u_old := FALSE]; KE[is.na(u_new), u_new := FALSE]
KE[, group := fifelse(dgp == "GenBinaryY" & label == "manual", "D-01 + D-08", fifelse(dgp == "GLMM-AR1" & scen == "M5", "D-02 + D-08", "D-08 only"))]
D8 <- KE[, {
  b <- u_old & u_new; dr <- (p_new[b] < 0.05) - (p_old[b] < 0.05); du <- u_new - u_old
  .(n_den = .N, n_both = sum(b), rej_old = mean(p_old[b] < 0.05), rej_new = mean(p_new[b] < 0.05),
    d_rej = mean(dr), mcse_d_rej = sd(dr) / sqrt(sum(b)), d_usable = mean(du), mcse_d_usable = sd(du) / sqrt(.N))
}, by = .(group, dgp, scen, N, ntime, label)]
D8[, `:=`(arm = ifelse(label == "mmlb", "H1", "H0"), z_rej = d_rej / mcse_d_rej, z_usable = d_usable / mcse_d_usable)]
fwrite(D8, file.path(TABDIR, "V4_S4_d08_lrt_effect.csv"))
D8S <- D8[, .(cells = .N, mean_d_rej = mean(d_rej, na.rm = TRUE), cells_rej_over2 = sum(abs(z_rej) > 2, na.rm = TRUE),
              mean_d_usable = mean(d_usable), cells_usable_over2 = sum(abs(z_usable) > 2, na.rm = TRUE), min_n_both = min(n_both)),
          by = .(group, dgp, scen, arm)][order(group, dgp, scen, arm)]
fwrite(D8S, file.path(TABDIR, "V4_S4_d08_summary.csv")); cat("Q8\n"); print(D8S)
Q8 <- D8[group == "D-08 only", .(cells = .N, cells_n_both_ge1 = sum(n_both >= 1))][, pass := cells_n_both_ge1 == cells]
nb0 <- D8[group == "D-08 only" & n_both == 0]
nb0[KE[, .(u_new_any = any(u_new), u_old_any = any(u_old)), by = .(dgp, scen, N, ntime, label)], on = c("dgp", "scen", "N", "ntime", "label"), `:=`(u_new_any = i.u_new_any, u_old_any = i.u_old_any)]
INV <- fread(file.path(LOG_S4, "S4_q1c_investigation.csv"))[method == LRT]
nb0[, in_q1c := FALSE]; nb0[INV, on = c("dgp", "scen", "N", "ntime", "label"), in_q1c := TRUE]
Q8[, `:=`(cells_n_both_0 = nrow(nb0), of_which_v4_never_usable = sum(!nb0$u_new_any), of_which_in_q1c_list = sum(nb0$in_q1c))]
fwrite(nb0, file.path(LOG_S4, "S4_q8_cells_without_pairs.csv"))
print(Q8)

# ---- Q9: glmer_bobyqa ---------------------------------------------------------------
bq <- V[method == SENS]
q9p <- bq[, .(rows = .N, p_mismatch = sum(abs(p - 2 * pnorm(-abs(estimate / se))) > 0, na.rm = TRUE))]
nb <- dcast(P[usab == "U0" & set == "own" & method %in% c(LRT, SENS), .(dgp, scen, N, ntime, label, method, n_obs)],
            dgp + scen + N + ntime + label ~ method, value.var = "n_obs")
q9n <- nb[, .(cells = .N, cells_bobyqa_le_lrt = sum(glmer_bobyqa <= glmer_lrt), max_excess = max(glmer_bobyqa - glmer_lrt))]
Q9 <- cbind(q9p, q9n)[, pass := p_mismatch == 0 & cells_bobyqa_le_lrt == cells]
fwrite(Q9, file.path(TABDIR, "V4_S4_q9_bobyqa.csv")); cat("Q9\n"); print(Q9)
fwrite(rbindlist(list(Q8[, .(check = "Q8", pass)], Q9[, .(check = "Q9", pass)], Q7[, .(check = "Q7", pass)],
                      Q4[, .(check = "Q4", pass)], Q3[, .(check = paste("Q3", table), pass)])), file.path(LOG_S4, "S4_checks_pass.csv"))
