# =============================================================================
# agg_v4.R -- S4: single re-aggregation of the v4 row set with Monte Carlo SEs.
# Measures and MCSE: protocol sec. 5 (Morris et al. 2019, Table 6); technical
# rules D-19; truths D-12/D-18; usability D-04; denominator D-11; rows D-15.
# Inputs : results_v4/S4/rows_v4.rds, results_v4/S4/truth_ar1_rep.rds,
#          tables/V4_S4_truth_by_cell.csv, tables/V4_S2_degenerate_designs.csv
# Outputs: tables/V4_perf_all.csv.gz (every cell x method x U x set),
#          tables/V4_perf_primary.csv (U0, own), results_v4/S4/boot_idx.rds,
#          tables/V4_S4_q1_q2.csv (composition and denominator checks)
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
t0 <- Sys.time()
B <- 2000L; SEED_BOOT <- 20261015L
V  <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
TR <- readRDS(file.path(OUT_S4, "truth_ar1_rep.rds"))
TT <- fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]

# ---- Q1 / Q2 checks on the row set ----------------------------------------------
dup <- V[!is.na(term), .(cnt = .N), by = c(CELL, "method", "sim_id")][cnt > 1, .N]
lrt35 <- V[method == LRT & model %in% 3:5 & !is.na(term)]
re_exp <- function(scen) fifelse(scen == "M1c", "(1|id)", fifelse(scen == "M5", "(1 + time + dose|id)", "(1 + time|id)"))
q1a_src  <- lrt35[scen != "M1c", sum(src == "full")]
q1a_re   <- V[method == LRT & (model %in% 3:5 | scen == "M1c") & !is.na(term), sum(re_structure != re_exp(scen) | is.na(re_structure))]
f12 <- V[(dgp == "GenBinaryY" & scen %in% c("M3", "M4", "M5") & label == "manual") | (dgp == "GLMM-AR1" & scen == "M5")]
q1b <- f12[, sum(!src %in% c("F1", "F2"))]
cells <- unique(V[, ..CELL]); dps <- unique(V[, ..DP])
mpc <- V[!is.na(term) & method %in% PRIM, .(nm = uniqueN(method)), by = CELL]
mpc_any <- V[method %in% PRIM, .(nm = uniqueN(method)), by = CELL]
miss_m <- CJ_cells_methods <- unique(V[, ..CELL])[, .(method = PRIM), by = CELL][!unique(V[!is.na(term) & method %in% PRIM, c(CELL, "method"), with = FALSE]), on = c(CELL, "method")]
fwrite(miss_m, file.path(LOG_S4, "S4_q1c_cells_methods_without_timetrt.csv"))
INV <- V[miss_m, on = c(CELL, "method"), nomatch = NULL][, .(n_rows = .N, n_failure_rows = sum(is.na(term)),
         top_status = names(sort(table(substr(status, 1, 70)), decreasing = TRUE))[1]), by = c(CELL, "method")]
INV <- merge(miss_m, INV, by = c(CELL, "method"), all.x = TRUE)
INV[is.na(n_rows), `:=`(n_rows = 0L, n_failure_rows = 0L)]
GC <- V[method == "glmer", .(glmer_converged_any = any(converged & !is.na(term))), by = CELL]
INV <- merge(INV, GC, by = CELL, all.x = TRUE)
fwrite(INV, file.path(LOG_S4, "S4_q1c_investigation.csv"))
print(INV[, .(combos = .N, with_rows = sum(n_rows > 0), all_failure = sum(n_rows > 0 & n_rows == n_failure_rows),
              glmer_never_converged = sum(!glmer_converged_any)), by = .(method, scen, ntime)][order(method, scen, ntime)])
V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
q2_u0_deg <- V[deg == TRUE, sum(U0)]
KEYS <- unique(V[, c(DP, "model"), with = FALSE])
KEYS <- KEYS[, .(sim_id = 1:1000), by = c(DP, "model")]
KEYS[, deg := FALSE]; KEYS[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
ND <- KEYS[, .(n_den = sum(!deg)), by = DP]
has <- unique(V[!is.na(term), c(CELL, "sim_id"), with = FALSE])
nd_keys <- KEYS[deg == FALSE][, .(label = c("mmlb", "manual")), by = c(DP, "sim_id")]
q2_noanalysis <- nd_keys[!has, on = c(CELL, "sim_id"), .N]
Q12 <- data.table(
  check = c("Q1 duplicate (cell, method, sim_id) rows", "Q1a glmer_lrt rows from results_full in M3-M5",
            "Q1a glmer_lrt rows with unexpected re_structure (M3-M5, M1c)", "Q1b F1/F2 cell rows not from S3",
            "Q1c cells", "Q1c design points", "Q1c cells with fewer than 8 primary methods (any row)", "Q1c cells with fewer than 8 primary methods (time:trt rows; descriptive)",
            "Q2 sum n_den main grid", "Q2 sum n_den M1c", "Q2 U0 rows on degenerate keys",
            "Q2 non-degenerate keys without any analysis row"),
  observed = c(dup, q1a_src, q1a_re, q1b, nrow(cells), nrow(dps), sum(mpc_any$nm < 8) + (nrow(cells) - nrow(mpc_any)), sum(mpc$nm < 8) + (nrow(cells) - nrow(mpc)),
               2L * ND[scen != "M1c", sum(n_den)] , 2L * ND[scen == "M1c", sum(n_den)], q2_u0_deg, q2_noanalysis),
  expected = c(0, 0, 0, 0, 3024, 1512, 0, NA, 2878518, 143742, 0, 0))
Q12[, pass := observed == expected]
if (Sys.getenv("S4_DEV") != "1") fwrite(Q12, file.path(TABDIR, "V4_S4_q1_q2.csv")); print(Q12)

# ---- bootstrap index counts, one matrix per scenario model (D-19) ----------------
set.seed(SEED_BOOT, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
BIDX <- lapply(1:5, function(m) matrix(sample.int(1000L, B * 1000L, replace = TRUE), nrow = B))
if (Sys.getenv("S4_DEV") != "1") saveRDS(BIDX, file.path(OUT_S4, "boot_idx.rds"))
CNT <- lapply(BIDX, function(M) { C <- matrix(0L, B, 1000L); for (b in seq_len(B)) C[b, ] <- tabulate(M[b, ], 1000L); C })

# ---- measures --------------------------------------------------------------------
perf <- function(est, se, p, th, n_den, est_measures, truth_ok) {
  n <- length(p); out <- c(n_obs = n, p_usable = n / n_den, mcse_p_usable = mcse_p(n / n_den, n_den))
  rej <- if (n > 0) mean(p < 0.05) else NA_real_
  out <- c(out, rej = rej, mcse_rej = mcse_p(rej, n), uncond_valid = (n - sum(p < 0.05)) / n_den,
           mcse_uncond = mcse_p((n - sum(p < 0.05)) / n_den, n_den))
  na <- c(bias = NA, mcse_bias = NA, empse = NA, mcse_empse = NA, modse = NA, mcse_modse = NA, relerr_modse = NA,
          mcse_relerr = NA, mse = NA, mcse_mse = NA, rmse = NA, mcse_rmse = NA, cover = NA, mcse_cover = NA,
          becover = NA, mcse_becover = NA)
  if (!est_measures || n < 2) return(c(out, na))
  e <- est - th
  empse <- sd(est); s2 <- se^2; modse <- sqrt(mean(s2)); v2 <- var(s2)
  becov <- mean(abs(est - mean(est)) <= Z975 * se)
  r <- na
  r[c("empse", "mcse_empse", "modse", "mcse_modse", "relerr_modse", "mcse_relerr", "becover", "mcse_becover")] <-
    c(empse, empse / sqrt(2 * (n - 1)), modse, sqrt(v2 / (4 * n * modse^2)), 100 * (modse / empse - 1),
      100 * (modse / empse) * sqrt(v2 / (4 * n * modse^4) + 1 / (2 * (n - 1))), becov, mcse_p(becov, n))
  if (truth_ok) {
    mse <- mean(e^2); cov <- mean(abs(e) <= Z975 * se)
    r[c("bias", "mcse_bias", "mse", "mcse_mse", "rmse", "mcse_rmse", "cover", "mcse_cover")] <-
      c(mean(e), sd(e) / sqrt(n), mse, sqrt(var(e^2) / n), sqrt(mse), sqrt(var(e^2) / n) / (2 * sqrt(mse)), cov, mcse_p(cov, n))
  }
  c(out, r)
}
# type 7 quantile at 0.05 of a bootstrap sample given counts over sorted values
boot_sa <- function(C, p0, s0, p1, s1) {
  o0 <- s0[order(p0[s0])]; v0 <- p0[o0]; o1 <- s1[order(p1[s1])]; v1 <- p1[o1]
  if (length(o0) < 2 || length(o1) < 1) return(c(se_sa = NA, se_raw = NA))
  rc <- function(M) if (ncol(M) == 1L) M else t(apply(M, 1, cumsum))   # row cumsums, keeps B x n shape
  cs0 <- rc(C[, o0, drop = FALSE]); n0 <- cs0[, ncol(cs0)]
  h <- (n0 - 1) * 0.05 + 1; lo <- floor(h); fr <- h - lo
  j_lo <- rowSums(cs0 < lo) + 1L; j_hi <- rowSums(cs0 < pmin(lo + 1, n0)) + 1L
  q <- v0[j_lo] + fr * (v0[j_hi] - v0[j_lo]); q[n0 < 2] <- NA
  cs1 <- cbind(0L, rc(C[, o1, drop = FALSE])); n1 <- cs1[, ncol(cs1)]
  k_sa <- cs1[cbind(seq_len(nrow(C)), findInterval(q, v1, left.open = TRUE) + 1L)]
  k_raw <- cs1[, sum(v1 < 0.05) + 1L]
  pw <- k_sa / n1; pr <- k_raw / n1
  c(se_sa = sd(pw, na.rm = TRUE), se_raw = sd(pr, na.rm = TRUE))
}

TRU <- function(dgp_, scen_, ntime_, label_) TT[dgp == dgp_ & scen == scen_ & ntime == ntime_ & label == label_]
MARG_V3 <- TT[dgp == "GenBinaryY" & label == "mmlb" & scen %in% sprintf("M%d", 1:5), .(th = truth_marginal[1]), by = scen]
setkey(TR, dgp, scen, N, ntime, label, sim_id)
V <- V[!is.na(term)]
V[, c("term", "stat", "converged", "status", "re_structure", "block", "design_degenerate", "src") := NULL]
SPL <- split(V, by = DP, keep.by = TRUE); rm(V); gc(verbose = FALSE)
DEV <- Sys.getenv("S4_DEV") == "1"
if (DEV) { SPL <- SPL[c(1, 700, 1505)]; TABDIR <- tempdir() }

one_dp <- function(D) {
  k <- D[1]; m <- k$model; C <- CNT[[m]]
  degm <- KEYS[dgp == k$dgp & scen == k$scen & N == k$N & ntime == k$ntime, deg]
  ok <- !degm; n_den <- sum(ok)
  arr <- function(lab, meth, col) { x <- rep(NA, 1000); r <- D[label == lab & method == meth]; x[r$sim_id] <- r[[col]]; x }
  get <- function(lab) lapply(setNames(METHODS, METHODS), function(meth) {
    u <- lapply(c(U0 = "U0", U1 = "U1", U2 = "U2"), function(uc) { x <- arr(lab, meth, uc); x[is.na(x)] <- FALSE; x & ok })
    list(u = u, est = arr(lab, meth, "estimate"), se = arr(lab, meth, "se"), p = arr(lab, meth, "p"))
  })
  G <- list(mmlb = get("mmlb"), manual = get("manual"))
  inter <- lapply(G, function(g) lapply(c(U0 = "U0", U1 = "U1", U2 = "U2"), function(uc) Reduce(`&`, lapply(PRIM, function(mm) g[[mm]]$u[[uc]]))))
  res <- list(); r <- 0L
  for (lab in c("mmlb", "manual")) {
    tt <- TRU(k$dgp, k$scen, k$ntime, lab)
    trc <- if (k$dgp == "GLMM-AR1") { x <- rep(NA_real_, 1000); z <- TR[.(k$dgp, k$scen, k$N, k$ntime, lab)]; x[z$sim_id] <- z$truth_cond; x } else NULL
    for (meth in METHODS) for (uc in c("U0", "U1", "U2")) for (st in c("own", "intersection")) {
      g <- G[[lab]][[meth]]; sel <- g$u[[uc]]; if (st == "intersection") sel <- sel & inter[[lab]][[uc]]
      s <- which(sel)
      if (meth == LRT) { th <- NA_real_; ttype <- "none"; tok <- FALSE; em <- FALSE }
      else if (meth %in% MARGINAL_METHODS) { th <- rep(tt$truth_marginal, length(s)); ttype <- "marginal"; tok <- TRUE; em <- TRUE }
      else if (k$dgp == "GenBinaryY") { th <- NA_real_; ttype <- "none (D-18)"; tok <- FALSE; em <- TRUE }
      else if (k$dgp == "GLMM-native") { th <- rep(tt$truth_conditional, length(s)); ttype <- "conditional"; tok <- TRUE; em <- TRUE }
      else { th <- trc[s]; ttype <- "conditional per replicate"; tok <- TRUE; em <- TRUE }
      v <- perf(g$est[s], g$se[s], g$p[s], th, n_den, em, tok)
      cv3 <- NA_real_
      if (meth %in% WALD && k$scen != "M1c" && length(s) > 0) cv3 <- mean(abs(g$est[s] - MARG_V3[scen == k$scen]$th) <= Z975 * g$se[s])
      sa <- c(crit_p_sa = NA, power_sa = NA, mcse_power_sa = NA, boot_se_power_raw = NA)
      if (lab == "mmlb") {
        g0 <- G$manual[[meth]]; sel0 <- g0$u[[uc]]; if (st == "intersection") sel0 <- sel0 & inter$manual[[uc]]
        s0 <- which(sel0)
        if (length(s0) >= 2 && length(s) >= 1) {
          # Superseded for the paper by D-73 (S8b): the size-adjusted power printed in the manuscript comes from
          # R/s5/s5_build.R with the 95th percentile of |estimate/se|; crit_p_sa/power_sa here keep the old p rule.
          q <- quantile(g0$p[s0], 0.05, type = 7, names = FALSE)
          bs <- boot_sa(C, g0$p, s0, g$p, s)
          sa <- c(crit_p_sa = q, power_sa = mean(g$p[s] < q), mcse_power_sa = bs[["se_sa"]], boot_se_power_raw = bs[["se_raw"]])
        }
      }
      r <- r + 1L
      res[[r]] <- c(list(dgp = k$dgp, scen = k$scen, model = m, N = k$N, ntime = k$ntime, label = lab,
                         arm = ifelse(lab == "mmlb", "H1", "H0"), method = meth, usab = uc, set = st, n_den = n_den,
                         truth_type = ttype, truth_value = if (tok) mean(th) else NA_real_),
                    as.list(v), list(cover_v3marg = cv3), as.list(sa))
    }
  }
  rbindlist(res)
}
P <- parallel::mclapply(SPL, function(D) tryCatch(one_dp(D), error = function(e) conditionMessage(e)), mc.cores = 40, mc.preschedule = FALSE)
bad <- which(vapply(P, is.character, logical(1)))
if (length(bad)) { print(head(unlist(P[bad]))); stop("failed design points: ", length(bad)) }
P <- rbindlist(P)
P[arm == "H1", uncond_valid := NA_real_][arm == "H1", mcse_uncond := NA_real_]
setorder(P, dgp, scen, N, ntime, label, method, usab, set)
fwrite(P, file.path(TABDIR, "V4_perf_all.csv.gz"))
fwrite(P[usab == "U0" & set == "own"], file.path(TABDIR, "V4_perf_primary.csv"))
cat(sprintf("rows %d (expected 163296) | primary %d | elapsed %.1f min\n", nrow(P), P[usab == "U0" & set == "own", .N],
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
