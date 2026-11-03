
# s10s_kccmp.R -- S10 "supplemento": explain small differences between the independent Kauermann-Carroll rows
# (fewclust SE from S7d, exchangeable) and the gee_kc rows of the variant used by the paper (variant_nc).
# Read-only use of results_v4/variant_nc/S4/rows_v4.rds only to explain differences (not to compute values).
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2)
RP <- "/home/user/MMM/results_v4"; D <- file.path(RP, "S10/supplemento"); K <- c("dgp","scen","N","ntime","label","sim_id")
V <- readRDS(file.path(RP, "variant_nc/S4/rows_v4.rds")); setDT(V)
V <- V[method %in% c("gee_kc", "mmm", "gee"), c(K, "method", "estimate", "se", "p", "converged", "U0"), with = FALSE]
G <- readRDS(file.path(D, "gee_s7d.rds"))[cor == "exch", c(K, "gee_ok", "est", "se_kc"), with = FALSE]
R <- readRDS(file.path(D, "recs_s7d.rds"))[, c(K, "n_obs", "sum_y", "design_degenerate"), with = FALSE]
M <- merge(V[method == "gee_kc"], G, by = K, all = TRUE); M <- merge(M, R, by = K, all.x = TRUE)
M[, const := sum_y == 0 | sum_y == n_obs]
M[, mine := gee_ok %in% TRUE & is.finite(est) & is.finite(se_kc) & abs(est) < 10 & se_kc > 0 & se_kc < 10 & !(const %in% TRUE)]
M[, theirs := U0 %in% TRUE]
cat("gee_kc rows in nc:", V[method == "gee_kc", .N], " S7d exch rows:", nrow(G), "\n")
print(M[, .N, by = .(mine, theirs)])
print(M[mine != theirs, .N, by = .(scen, label, const, gee_ok, conv = converged, est_fin = is.finite(estimate), se_eq = abs(se - se_kc) < 1e-10)][order(-N)][1:30])
print(M[mine & theirs, .(n = .N, se_diff = sum(abs(se - se_kc) > 1e-10), est_diff = sum(abs(estimate - est) > 1e-10))])
x <- M[mine != theirs][1:10]; print(x[, c(K, "estimate", "se", "est", "se_kc", "converged", "gee_ok", "U0"), with = FALSE])
