# =============================================================================
# s8b_pb2_collapsed.R -- S8b, prediction PB2: share of the size-adjusted power of the mixed-model Wald test (glmer,
# Laplace) that comes from fits with a collapsed standard error (|estimate/se| > 37, p numerically 0).
# Rows: results_v4/variant_kcfc/S4/rows_v4.rds, main grid M1-M5, method glmer, usable = U0 and estimable design (D-11),
# definitions of R/s5/s5_build.R. Threshold per cell: 95th percentile (type 7) of the usable null |z| (D-73).
# Per H1 cell: n1 usable, k = #(|z1| > q95), kbig = #(|z1| > q95 & |z1| > 37); Q_cell = 100 kbig / n1.
# Aggregate: mean over cells of Q_cell, and pooled 100 sum(kbig) / sum(n1). Interval: block bootstrap over replicate
# index within model (B = 2000, seed 20261015), jointly over cells, threshold held at its full-sample value (declared:
# in cells with q95 < 37 kbig does not depend on the threshold).
# Output: tables/V5_S8b_pb2_collapsed.csv (cells) and tables/V5_S8b_pb2_summary.csv.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L); setwd("/home/user/MMM")
V <- readRDS("results_v4/variant_kcfc/S4/rows_v4.rds")[!is.na(term) & method == "glmer" & scen %in% sprintf("M%d", 1:5),
       .(dgp, scen, N, ntime, label, sim_id, estimate, se, U0)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
DEG <- fread("tables/V4_S2_degenerate_designs.csv")[, .(model, N, ntime, sim_id)]
V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]; V[, u0 := U0 %in% TRUE & !deg]
V <- V[u0 == TRUE]; V[, za := abs(estimate / se)]
k <- c("dgp", "scen", "N", "ntime")
Q0 <- V[label == "manual", .(n0 = .N, q95 = if (.N >= 2) quantile(za, 0.95, type = 7, names = FALSE) else NA_real_, big0 = sum(za > 37)), by = k]
H <- merge(V[label == "mmlb"], Q0, by = k)
H <- H[is.finite(q95)]
H[, `:=`(rej = za > q95, big = za > q95 & za > 37)]
C <- H[, .(n1 = .N, k = sum(rej), kbig = sum(big), q95 = q95[1], n0 = n0[1], big0 = big0[1]), by = c(k, "model")]
C[, `:=`(pow_sa = 100 * k / n1, Q = 100 * kbig / n1)]
fwrite(C, "tables/V5_S8b_pb2_collapsed.csv")
# block bootstrap
set.seed(20261015L); B <- 2000L
C[, key := paste(dgp, scen, N, ntime)]; H[, key := paste(dgp, scen, N, ntime)]
est <- function(W) {   # W: list by model of 1000 x (B+1) weight matrices
  out <- matrix(0, 2, B + 1L); nc <- 0
  num_pool <- rep(0, B + 1L); den_pool <- rep(0, B + 1L); sumQ <- rep(0, B + 1L); ncb <- rep(0, B + 1L)
  for (m in 1:5) { Hm <- H[model == m]; keys <- C[model == m, key]
    Ib <- matrix(0, length(keys), 1000L); In <- Ib
    Ib[cbind(match(Hm$key, keys), Hm$sim_id)] <- as.numeric(Hm$big); In[cbind(match(Hm$key, keys), Hm$sim_id)] <- 1
    nb <- Ib %*% W[[m]]; nn <- In %*% W[[m]]
    Qm <- 100 * nb / nn; Qm[nn == 0] <- NA   # a cell with no usable H1 dataset in a bootstrap draw is left out of that draw's mean
    sumQ <- sumQ + colSums(Qm, na.rm = TRUE); ncb <- ncb + colSums(!is.na(Qm)); num_pool <- num_pool + colSums(nb); den_pool <- den_pool + colSums(nn) }
  rbind(mean = sumQ / ncb, pooled = 100 * num_pool / den_pool)
}
W <- lapply(1:5, function(m) { Cm <- matrix(0L, 1000L, B + 1L); Cm[, 1] <- 1L; for (b in 1:B) Cm[, b + 1] <- tabulate(sample.int(1000L, 1000L, TRUE), 1000L); Cm })
E <- est(W)
S <- data.table(quantity = rownames(E), est = E[, 1], mcse = apply(E[, -1], 1, sd), lo = apply(E[, -1], 1, quantile, 0.025), hi = apply(E[, -1], 1, quantile, 0.975),
                cells = nrow(C), cells_q95_above_37 = C[q95 > 37, .N], cells_Q_above_2.5 = C[Q > 2.5, .N], max_Q = max(C$Q),
                mean_pow_sa = mean(C$pow_sa), null_big_share = 100 * sum(C$big0) / sum(C$n0))
fwrite(S, "tables/V5_S8b_pb2_summary.csv"); print(S, digits = 4)
print(C[order(-Q)][1:10], digits = 4)
print(C[, .(cells = .N, meanQ = mean(Q), meanpow = mean(pow_sa)), by = dgp], digits = 4)
