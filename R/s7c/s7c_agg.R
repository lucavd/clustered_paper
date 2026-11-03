# s7c_agg.R -- S7c aggregation: gates GD3/GD5, per-cell measures, predictions PD1-PD9 with block bootstrap (B = 2000, seed 20261015).
source(path.expand("~/MMM/R/s7c/s7c_common.R"))
TAB <- file.path(REPO, "tables")
tk <- pilot_tasks()
X <- lapply(tk$file, readRDS)
recs  <- rbindlist(lapply(X, `[[`, "recs"), fill = TRUE)
tests <- rbindlist(lapply(X, `[[`, "tests"), fill = TRUE)
# ---- saved study rows (results_full) for the same replicates --------------------------------------------
u0 <- function(converged, est, se, p) { w <- converged & is.finite(est) & is.finite(se) & abs(est) < 10 & se > 0 & se < 10 & is.finite(p); w[is.na(w)] <- FALSE; w }
S <- rbindlist(lapply(seq_len(nrow(tk)), function(k) {
  o <- as.data.table(readRDS(tk$orig[k]))[term == "time:trt" & method %in% c("gee", "gee_fg", "mmm")]
  o[, .(dgp = tk$dgp[k], N = tk$N[k], ntime = tk$ntime[k], label = tk$label[k], sim_id = as.integer(sim_id), method, estimate, se, p, converged)]
}))
S[, usable := u0(converged, estimate, se, p)]
S[, p_t := fifelse(usable, 2 * pt(-abs(estimate / se), N - 4), NA_real_)]   # t with N - p, p = 4 in M1
# ---- GD3: regenerated GEE exchangeable equals saved -----------------------------------------------------
g3 <- merge(recs[, .(dgp, N, ntime, label, sim_id, gee_est, gee_se)], S[method == "gee", .(dgp, N, ntime, label, sim_id, s_est = estimate, s_se = se)],
            by = c("dgp", "N", "ntime", "label", "sim_id"), all = TRUE)
g3[, fin := is.finite(s_est) & is.finite(s_se)]
g3[, d := pmax(abs(gee_est - s_est), abs(gee_se - s_se))]
GD3 <- g3[fin == TRUE, .(n_compared = .N, n_bad = sum(!(d <= 1e-8 * pmax(1, abs(s_est), abs(s_se))) | is.na(d)), max_abs = max(d, na.rm = TRUE))]
# ---- GD5: completeness -------------------------------------------------------------------------------
TESTS <- c("flip", "subj_t", "subj_perm", "naive_perm", "naive_slope")
cnt <- tests[, .(k = uniqueN(test), has_all = all(TESTS %in% test), ok_row = all(is.finite(p) | !is.na(why))), by = .(dgp, N, ntime, label, sim_id)]
GD5 <- data.table(n_datasets = nrow(recs), n_generated = sum(recs$generated), n_gen_error = sum(!recs$generated),
                  gen_error_degenerate = sum(!recs$generated & recs$design_degenerate),
                  n_with_5 = sum(cnt$has_all & cnt$k == 5), n_rows_ok = sum(cnt$ok_row))
GD5[, pass := n_datasets == 36000 & n_with_5 == n_generated & n_rows_ok == n_generated]
# ---- per-replicate wide table ------------------------------------------------------------------------
W <- dcast(tests, dgp + N + ntime + label + sim_id ~ test, value.var = "p")
W <- merge(recs[, .(dgp, N, ntime, label, sim_id, generated, design_degenerate)], W, by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE)
SW <- dcast(S, dgp + N + ntime + label + sim_id ~ method, value.var = "p_t")
setnames(SW, c("gee", "gee_fg", "mmm"), c("gee_t", "gee_fg_t", "mmm_t"))
W <- merge(W, SW, by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE)
W <- W[design_degenerate == FALSE]          # D-11 denominator: estimable designs only
METH <- c(TESTS, "gee_fg_t", "mmm_t", "gee_t")
L <- melt(W, id.vars = c("dgp", "N", "ntime", "label", "sim_id"), measure.vars = METH, variable.name = "method", value.name = "p")
L[, method := as.character(method)]
C <- L[, .(nd = .N, n = sum(is.finite(p)), k = sum(p < 0.05, na.rm = TRUE)), by = .(dgp, N, ntime, label, method)]
C[, `:=`(rate = 100 * k / n, pu = 100 * n / nd)]
C[, mc := 100 * sqrt((rate / 100) * (1 - rate / 100) / n)]
C[, `:=`(assess = n >= 20, hi = n >= 20 & rate - 2 * mc > 7.5, lo = n >= 20 & rate + 2 * mc < 2.5)]
C[, rate_all := 100 * k / nd]
fwrite(C[order(label, method, dgp, N, ntime)], file.path(TAB, "V4_S7c_cells.csv"))
# ---- block bootstrap over replicate index i (model fixed = 1) ----------------------------------------
set.seed(20261015L); BI <- matrix(sample.int(1000L, 2000L * 1000L, replace = TRUE), 2000L)
Wt <- t(apply(BI, 1, tabulate, nbins = 1000L))             # 2000 x 1000 multiplicities
cellmat <- function(lab, expr_rej, expr_use) {             # cells x 1000 matrices for a per-replicate indicator pair
  D <- W[label == lab]; D[, cell := paste(dgp, N, ntime)]
  cells <- sort(unique(D$cell)); R <- U <- matrix(0, length(cells), 1000L)
  r <- eval(expr_rej, D); u <- eval(expr_use, D); r[!u | is.na(r)] <- 0; u[is.na(u)] <- FALSE
  R[cbind(match(D$cell, cells), D$sim_id)] <- r; U[cbind(match(D$cell, cells), D$sim_id)] <- as.numeric(u)
  list(R = R, U = U, cells = cells)
}
boot_mean <- function(M) {                                  # mean over cells of the cell rate; estimate and 95% percentile interval
  est <- mean(rowSums(M$R) / rowSums(M$U), na.rm = TRUE)
  bt <- rowMeans((Wt %*% t(M$R)) / (Wt %*% t(M$U)), na.rm = TRUE)
  c(est = 100 * est, lo = 100 * unname(quantile(bt, 0.025, na.rm = TRUE)), hi = 100 * unname(quantile(bt, 0.975, na.rm = TRUE)))
}
boot_diff <- function(M1, M2) {                             # mean over cells of (rate1 - rate2), same usable set per cell
  est <- mean(rowSums(M1$R) / rowSums(M1$U) - rowSums(M2$R) / rowSums(M2$U), na.rm = TRUE)
  bt <- rowMeans((Wt %*% t(M1$R)) / (Wt %*% t(M1$U)) - (Wt %*% t(M2$R)) / (Wt %*% t(M2$U)), na.rm = TRUE)
  c(est = 100 * est, lo = 100 * unname(quantile(bt, 0.025, na.rm = TRUE)), hi = 100 * unname(quantile(bt, 0.975, na.rm = TRUE)))
}
rej <- function(m) substitute(get(m) < 0.05, list(m = m)); use <- function(m) substitute(is.finite(get(m)), list(m = m))
both <- function(a, b) substitute(is.finite(get(a)) & is.finite(get(b)), list(a = a, b = b))
agg <- list()
for (m in METH) agg[[m]] <- boot_mean(cellmat("manual", rej(m), use(m)))
# PD6: D - C paired on replicates usable for both
pd6 <- boot_diff(cellmat("manual", rej("naive_perm"), both("naive_perm", "subj_perm")), cellmat("manual", rej("subj_perm"), both("naive_perm", "subj_perm")))
# PD9: power differences vs gee_fg_t, paired
pd9 <- lapply(c("flip", "subj_t", "subj_perm"), function(m)
  boot_diff(cellmat("mmlb", rej(m), both(m, "gee_fg_t")), cellmat("mmlb", rej("gee_fg_t"), both(m, "gee_fg_t"))))
names(pd9) <- c("flip", "subj_t", "subj_perm")
pwr <- list(); for (m in METH) pwr[[m]] <- boot_mean(cellmat("mmlb", rej(m), use(m)))
# PD5: usable share differences per cell (H0 and H1 cells)
U5 <- merge(C[method %in% c("flip", "subj_t", "subj_perm"), .(dgp, N, ntime, label, method, pu, n, nd)],
            C[method == "gee_fg_t", .(dgp, N, ntime, label, pu_fg = pu)], by = c("dgp", "N", "ntime", "label"))
U5[, dpu := pu - pu_fg]
# MCSE of a usable-share difference (paired replicates, conservative independent approximation)
U5[, mc_d := 100 * sqrt((pu / 100) * (1 - pu / 100) / nd + (pu_fg / 100) * (1 - pu_fg / 100) / nd)]
# PD7: subjects without a slope at T = 2 (H0 arm, design identical across arms)
S7 <- tests[test == "naive_slope" & ntime == 2 & label == "manual", .(noslope = sum(n_noslope, na.rm = TRUE), subj = sum(N_subj, na.rm = TRUE))]
ns_all <- tests[test == "naive_slope" & label == "manual" & is.finite(n_noslope), .(share = 100 * sum(n_noslope) / sum(N_subj)), by = .(ntime)]
# PD8: expected excess rejection of E from the probability-scale null
z <- qnorm(0.975)
E8 <- tests[test == "naive_slope" & label == "manual" & is.finite(slope_diff), .(delta = mean(slope_diff), s = sd(slope_diff), n = .N), by = .(dgp, N, ntime)]
E8[, excess := 100 * (pnorm(-z + delta / s) + pnorm(-z - delta / s) - 0.05)]
E8[, excess_mc := 100 * dnorm(z - abs(delta) / s) * (s / sqrt(n)) / s]   # delta-method MCSE of the excess
fwrite(E8, file.path(TAB, "V4_S7c_slope_null.csv"))
# ---- predictions table --------------------------------------------------------------------------------
judge <- function(lo, hi, a, b) if (lo >= a && hi <= b) "confermata" else if (hi < a || lo > b) "smentita" else "indeterminata"
P <- rbindlist(list(
  data.table(code = "PD1", what = "type I error A (flip), mean of 18 H0 cells", est = agg$flip[1], lo = agg$flip[2], hi = agg$flip[3], ra = 2.5, rb = 7.5),
  data.table(code = "PD2", what = "type I error B (subj_t), mean of 18 H0 cells", est = agg$subj_t[1], lo = agg$subj_t[2], hi = agg$subj_t[3], ra = 2.5, rb = 7.5),
  data.table(code = "PD3", what = "type I error C (subj_perm), mean of 18 H0 cells", est = agg$subj_perm[1], lo = agg$subj_perm[2], hi = agg$subj_perm[3], ra = 2.5, rb = 7.5),
  data.table(code = "PD6", what = "type I error D - C, paired, mean of 18 H0 cells (points)", est = pd6[1], lo = pd6[2], hi = pd6[3], ra = -1, rb = 1),
  data.table(code = "PD9_flip", what = "power A - GEE FG t, paired, mean of 18 H1 cells (points)", est = pd9$flip[1], lo = pd9$flip[2], hi = pd9$flip[3], ra = -5, rb = 5),
  data.table(code = "PD9_subj_t", what = "power B - GEE FG t, paired, mean of 18 H1 cells (points)", est = pd9$subj_t[1], lo = pd9$subj_t[2], hi = pd9$subj_t[3], ra = -5, rb = 5),
  data.table(code = "PD9_subj_perm", what = "power C - GEE FG t, paired, mean of 18 H1 cells (points)", est = pd9$subj_perm[1], lo = pd9$subj_perm[2], hi = pd9$subj_perm[3], ra = -5, rb = 5)))
P[, outcome := mapply(judge, lo, hi, ra, rb)]
hiC <- C[label == "manual" & method %in% c("flip", "subj_t", "subj_perm"), .(n_hi = sum(hi), n_lo = sum(lo), n_assess = sum(assess)), by = method]
P <- rbind(P, data.table(code = paste0("PD4_", hiC$method), what = "H0 cells clearly above 7.5% (count of 18)", est = hiC$n_hi, lo = NA, hi = NA, ra = 0, rb = 0,
                         outcome = ifelse(hiC$n_hi == 0, "confermata", "smentita")), fill = TRUE)
pd5 <- U5[, .(worst = min(dpu), worst_lo = min(dpu - 2 * mc_d), worst_hi = min(dpu + 2 * mc_d), n_below = sum(dpu + 2 * mc_d < -2.5), mean_d = mean(dpu)), by = method]
P <- rbind(P, data.table(code = paste0("PD5_", pd5$method), what = "usable share minus GEE FG, worst of 36 cells (points); n_below = cells clearly below -2.5",
                         est = pd5$worst, lo = pd5$worst_lo, hi = pd5$worst_hi, ra = -2.5, rb = Inf,
                         outcome = ifelse(pd5$n_below > 0, "smentita", ifelse(pd5$worst_lo >= -2.5, "confermata", "indeterminata"))), fill = TRUE)
sh7 <- 100 * S7$noslope / S7$subj
P <- rbind(P, data.table(code = "PD7", what = "share of subjects without slope in E at T = 2 (%)", est = sh7, lo = NA, hi = NA, ra = 25, rb = 35,
                         outcome = if (sh7 >= 25 && sh7 <= 35) "confermata" else "smentita"), fill = TRUE)
mx <- E8[which.max(excess)]
P <- rbind(P, data.table(code = "PD8", what = "max expected excess rejection of E over 18 H0 cells (points)", est = mx$excess, lo = mx$excess - 2 * mx$excess_mc, hi = mx$excess + 2 * mx$excess_mc,
                         ra = 0, rb = 1), fill = TRUE)
P[code == "PD8", outcome := judge(lo, hi, ra, rb)]
fwrite(P, file.path(TAB, "V4_S7c_predictions.csv"))
A <- rbindlist(lapply(names(agg), function(m) data.table(method = m, arm = "H0", est = agg[[m]][1], lo = agg[[m]][2], hi = agg[[m]][3])))
A <- rbind(A, rbindlist(lapply(names(pwr), function(m) data.table(method = m, arm = "H1", est = pwr[[m]][1], lo = pwr[[m]][2], hi = pwr[[m]][3]))))
fwrite(A, file.path(TAB, "V4_S7c_aggregate.csv")); fwrite(U5, file.path(TAB, "V4_S7c_usable_vs_fg.csv"))
fwrite(ns_all, file.path(TAB, "V4_S7c_noslope_share.csv"))
why <- tests[!is.finite(p), .N, by = .(test, why)][order(test, -N)]; fwrite(why, file.path(TAB, "V4_S7c_unusable_reasons.csv"))
G <- data.table(gate = c("GD3", "GD5"), pass = c(GD3$n_bad == 0, GD5$pass),
                detail = c(sprintf("compared %d, mismatches %d, max abs diff %.3g", GD3$n_compared, GD3$n_bad, GD3$max_abs),
                           sprintf("datasets %d, generated %d (generation errors %d, of which degenerate design %d), with 5 test rows %d, rows ok %d",
                                   GD5$n_datasets, GD5$n_generated, GD5$n_gen_error, GD5$gen_error_degenerate, GD5$n_with_5, GD5$n_rows_ok)))
fwrite(G, file.path(OUT_S7C, "gates_post.csv"))
nb <- tests[test == "naive_perm", .(n_bad_perm = sum(n_bad, na.rm = TRUE), n_perm = sum(n_ref - 1L, na.rm = TRUE))]
cat("GATES\n"); print(G)
cat("naive_perm non-finite permutation statistics:", nb$n_bad_perm, "of", nb$n_perm, "\n")
cat("PREDICTIONS\n"); print(P[, .(code, est = round(est, 2), lo = round(lo, 2), hi = round(hi, 2), outcome)])
cat("AGGREGATE\n"); print(A[, .(method, arm, est = round(est, 2), lo = round(lo, 2), hi = round(hi, 2))])
cat("CELLS H0 (rate, pu)\n"); print(dcast(C[label == "manual"], dgp + N + ntime ~ method, value.var = "rate", fun.aggregate = function(x) round(x, 1)))
print(dcast(C[label == "manual"], dgp + N + ntime ~ method, value.var = "pu", fun.aggregate = function(x) round(x, 0)))
cat("CELLS H1 power\n"); print(dcast(C[label == "mmlb"], dgp + N + ntime ~ method, value.var = "rate", fun.aggregate = function(x) round(x, 1)))
cat("noslope share by T\n"); print(ns_all); cat("E8\n"); print(E8[order(-excess)][1:6])
cat("unusable reasons\n"); print(why)
