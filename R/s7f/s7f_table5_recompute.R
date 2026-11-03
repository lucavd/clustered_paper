
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
# s7f_table5_recompute.R -- independent recomputation of Table 5 (raw and size-adjusted power, main grid) from the
# per-replicate rows, with the definitions of R/s5/s5_build.R (u0 = U0 and non-degenerate design; size adjustment:
# 95th percentile, type 7, of the usable null |estimate/se| of the same cell, D-73) and R/s5/s5_engine.R (pooled sum k / sum n).
V <- readRDS("results_v4/variant_kcfc/S4/rows_v4.rds")[!is.na(term) & scen %in% sprintf("M%d", 1:5)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
DEG <- fread("tables/V4_S2_degenerate_designs.csv")[, .(model, N, ntime, sim_id)]
V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]; V[, u0 := U0 & !deg]
PRIM8 <- c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt"); V <- V[method %in% PRIM8]
k <- c("dgp", "scen", "N", "ntime", "method")
# DEC-T1 (S8b, D-73): Wald tests use the declared rule, 95th percentile (type 7) of the usable null |estimate/se|,
# written as the 5th percentile of -|z| (reflection-symmetric); the likelihood-ratio test stays on p.
V[, pv := fifelse(method == "glmer_lrt", p, -abs(estimate / se))]
Q0 <- V[label == "manual" & u0 %in% TRUE & is.finite(pv), .(q0 = if (.N >= 2) quantile(pv, 0.05, type = 7, names = FALSE) else NA_real_), by = k]
H1 <- V[label == "mmlb", .(n = sum(u0 %in% TRUE), kk = sum(u0 %in% TRUE & is.finite(p) & p < 0.05), p = list(pv[u0 %in% TRUE & is.finite(pv)])), by = k]
H1 <- merge(H1, Q0, by = k, all.x = TRUE)
H1[, ksa := mapply(function(pp, q) if (is.na(q) || length(pp) < 1) NA_real_ else sum(pp < q), p, q0)]
PW <- H1[, .(raw = 100 * sum(kk) / sum(n), sa = 100 * sum(ksa[is.finite(ksa)]) / sum(n[is.finite(ksa)])), by = .(dgp, method)]
S <- c(GenBinaryY = "GB", "GLMM-native" = "NAT", "GLMM-AR1" = "AR1"); PW[, g := S[dgp]]
TC <- fread("tables/variant_kcfc/V4_S6_table_cells.csv")[variant == "main" & table_id == "T10"]
TC[, c("q", "g", "m") := tstrsplit(name, "|", fixed = TRUE)]
J <- merge(TC[, .(q, g, m, est)], melt(PW[, .(g, m = method, raw, sa)], id.vars = c("g", "m"), variable.name = "q", value.name = "recomp"), by = c("q", "g", "m"))
print(J[, .(cells = .N, max_abs_diff = max(abs(recomp - est))), by = q]); fwrite(J, "tables/V5_S7f_table5_recompute.csv")
