# s7c_posthoc_agg.R -- aggregation of the post hoc variants (same denominators and bootstrap as s7c_agg.R)
source(path.expand("~/MMM/R/s7c/s7c_common.R"))
TAB <- file.path(REPO, "tables"); tk <- pilot_tasks()
PH <- rbindlist(lapply(sub("/pilot/", "/posthoc/", tk$file), readRDS), fill = TRUE)
X <- lapply(tk$file, readRDS); recs <- rbindlist(lapply(X, `[[`, "recs"), fill = TRUE); tests <- rbindlist(lapply(X, `[[`, "tests"), fill = TRUE)
Wp <- dcast(tests, dgp + N + ntime + label + sim_id ~ test, value.var = "p")
PH <- merge(recs[, .(dgp, N, ntime, label, sim_id, design_degenerate)], PH, by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE)
PH <- merge(PH, Wp[, .(dgp, N, ntime, label, sim_id, flip)], by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE)
PH <- PH[design_degenerate == FALSE]
# consistency: flip_relaxed equals flip wherever the strict filter passed
cons <- PH[is.finite(flip), .(n = .N, max_diff = max(abs(flip - p_flip_relaxed)))]
# saved gee_fg_t for power comparisons
u0 <- function(cv, est, se, p) { w <- cv & is.finite(est) & is.finite(se) & abs(est) < 10 & se > 0 & se < 10 & is.finite(p); w[is.na(w)] <- FALSE; w }
S <- rbindlist(lapply(seq_len(nrow(tk)), function(k) { o <- as.data.table(readRDS(tk$orig[k]))[term == "time:trt" & method == "gee_fg"]
  o[, .(dgp = tk$dgp[k], N = tk$N[k], ntime = tk$ntime[k], label = tk$label[k], sim_id = as.integer(sim_id), gee_fg_t = fifelse(u0(converged, estimate, se, p), 2 * pt(-abs(estimate / se), tk$N[k] - 4), NA_real_))] }))
PH <- merge(PH, S, by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE)
M <- c("p_flip_relaxed", "p_robust_t", "p_within_t", "p_within_perm", "flip", "gee_fg_t")
L <- melt(PH, id.vars = c("dgp", "N", "ntime", "label", "sim_id"), measure.vars = M, variable.name = "method", value.name = "p")
C <- L[, .(nd = .N, n = sum(is.finite(p)), k = sum(p < 0.05, na.rm = TRUE)), by = .(dgp, N, ntime, label, method)]
C[, `:=`(rate = 100 * k / n, pu = 100 * n / nd)]; C[, mc := 100 * sqrt((rate / 100) * (1 - rate / 100) / n)]
C[, `:=`(hi = n >= 20 & rate - 2 * mc > 7.5, lo = n >= 20 & rate + 2 * mc < 2.5)]
fwrite(C, file.path(TAB, "V4_S7c_posthoc_cells.csv"))
set.seed(20261015L); BI <- matrix(sample.int(1000L, 2000L * 1000L, replace = TRUE), 2000L); Wt <- t(apply(BI, 1, tabulate, nbins = 1000L))
cm <- function(lab, m, m2 = NULL) { D <- PH[label == lab]; D[, cell := paste(dgp, N, ntime)]; cells <- sort(unique(D$cell))
  p <- D[[m]]; u <- is.finite(p); if (!is.null(m2)) u <- u & is.finite(D[[m2]]); r <- u & p < 0.05
  R <- U <- matrix(0, length(cells), 1000L); R[cbind(match(D$cell, cells), D$sim_id)] <- r; U[cbind(match(D$cell, cells), D$sim_id)] <- u; list(R = R, U = U) }
bm <- function(M) { e <- mean(rowSums(M$R) / rowSums(M$U), na.rm = TRUE); b <- rowMeans((Wt %*% t(M$R)) / (Wt %*% t(M$U)), na.rm = TRUE)
  100 * c(est = e, lo = unname(quantile(b, .025, na.rm = TRUE)), hi = unname(quantile(b, .975, na.rm = TRUE))) }
bd <- function(m, ref) { A1 <- cm("mmlb", m, ref); A2 <- cm("mmlb", ref, m); e <- mean(rowSums(A1$R) / rowSums(A1$U) - rowSums(A2$R) / rowSums(A2$U), na.rm = TRUE)
  b <- rowMeans((Wt %*% t(A1$R)) / (Wt %*% t(A1$U)) - (Wt %*% t(A2$R)) / (Wt %*% t(A2$U)), na.rm = TRUE)
  100 * c(est = e, lo = unname(quantile(b, .025, na.rm = TRUE)), hi = unname(quantile(b, .975, na.rm = TRUE))) }
A <- rbindlist(lapply(M, function(m) data.table(method = m, h0 = t(bm(cm("manual", m))), h1 = t(bm(cm("mmlb", m))),
       pw_vs_fg = if (m != "gee_fg_t") t(bd(m, "gee_fg_t")) else t(c(NA, NA, NA)))))
fwrite(A, file.path(TAB, "V4_S7c_posthoc_aggregate.csv"))
# variance-ratio explanation for B (H0 cells): mean Student-t variance of m1 - m0 over its empirical variance, for r and for w
VR <- PH[label == "manual" & is.finite(dm_r), .(ratio_r = mean(v_r) / var(dm_r), ratio_w = mean(v_w) / var(dm_w), n = .N), by = .(dgp, N, ntime)]
fwrite(VR, file.path(TAB, "V4_S7c_variance_ratio.csv"))
cat("consistency flip_relaxed vs flip:", cons$n, "max diff", cons$max_diff, "\n")
print(A[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
cat("H0 rates\n"); print(dcast(C[label == "manual"], dgp + N + ntime ~ method, value.var = "rate", fun.aggregate = function(x) round(x, 1)))
cat("usable %\n"); print(dcast(C[label == "manual"], dgp + N + ntime ~ method, value.var = "pu", fun.aggregate = function(x) round(x, 0)))
cat("H1 power\n"); print(dcast(C[label == "mmlb"], dgp + N + ntime ~ method, value.var = "rate", fun.aggregate = function(x) round(x, 1)))
cat("flags hi/lo\n"); print(C[label == "manual", .(hi = sum(hi), lo = sum(lo)), by = method])
cat("variance ratio\n"); print(VR[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
