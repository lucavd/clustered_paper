# =============================================================================
# s9_const_keys.R -- S9 (D-82, D-83): keys of the non-degenerate datasets with a constant outcome, from the recs of
# the S7d run (results_v4/S7d/rows: GenBinaryY, GLMM-native, GLMM-AR1 on the 1,896 points; same seeds as the stored
# study rows) and of the S7e run (results_v4/S7e/rows: Copula-AR1, 480 points). Read-only on results_v4/*/rows.
# Outputs: results_v4/S9/const_keys.rds (K6 only), tables/V6_S9_const_keys_summary.csv (counts by dgp, label,
# T 2-3 and N <= 25), tables/V6_S9_const_scope_extra.csv (usable constant datasets in the stored rows that the step
# does not regenerate: variant_cc gee_md/gee_kc/gee_fg of geesmv, S7g adaptive-quadrature glmer).
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L)
REPO <- "/home/user/MMM"; setwd(REPO); dir.create("results_v4/S9", showWarnings = FALSE)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
RC0 <- c(K6, "generated", "design_degenerate", "n_obs", "sum_y")
f <- c(list.files("results_v4/S7d/rows", full.names = TRUE, pattern = "\\.rds$"), list.files("results_v4/S7e/rows", full.names = TRUE, pattern = "\\.rds$"))
R <- rbindlist(parallel::mclapply(f, function(x) as.data.table(readRDS(x)$recs)[, ..RC0], mc.cores = 6L))
stopifnot(nrow(R) == 4752000L, uniqueN(R[, ..K6]) == nrow(R))
R[, const := generated & !design_degenerate & (sum_y == 0 | sum_y == n_obs)]
K <- R[const == TRUE, ..K6]; setkeyv(K, K6)
cat("constant-outcome datasets:", nrow(K), "\n"); stopifnot(nrow(K) == 9711L)
saveRDS(K, "results_v4/S9/const_keys.rds")
SM <- R[, .(datasets = .N, nondegenerate = sum(!design_degenerate), constant = sum(const), all_zero = sum(const & sum_y == 0)),
        by = .(dgp, label, small = ntime <= 3 & N <= 25)][order(dgp, label, small)]
fwrite(SM, "tables/V6_S9_const_keys_summary.csv"); print(SM)
## stored rows the step does not regenerate
source("R/s9/s9_const.R")
V <- readRDS("results_v4/variant_cc/S4/rows_v4.rds")[method %in% c("gee", "gee_md", "gee_kc", "gee_fg") & term %in% "time:trt", c(K6, "method", "U0"), with = FALSE]
V[, u := U0 %in% TRUE]; V[, c0 := u]; s9_drop_const(V, "u")
EX <- V[, .(rows = .N, usable = sum(c0), usable_constant = sum(c0 & !u)), by = method][, source := "variant_cc rows_v4"]
rm(V); invisible(gc())
fg <- list.files("results_v4/S7g/rows", full.names = TRUE, pattern = "\\.rds$")
G <- rbindlist(lapply(fg, function(x) { y <- readRDS(x); r <- as.data.table(y$rows); r }), fill = TRUE)
cat("S7g row columns:", paste(names(G), collapse = ","), "\n")
u0 <- function(converged, estimate, se, p) { w <- converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p); w[is.na(w)] <- FALSE; w }
if (all(c("estimate", "se", "p", "converged") %in% names(G))) {
  if (!"method" %in% names(G)) G[, method := "glmer_agq"]
  G[, u := u0(converged, estimate, se, p)]; G[, c0 := u]; s9_drop_const(G, "u")
  mcol <- intersect(c("nAGQ", "agq", "method"), names(G))[1]
  EG <- G[, .(rows = .N, usable = sum(c0), usable_constant = sum(c0 & !u)), by = mcol]; setnames(EG, mcol, "method"); EG[, `:=`(method = as.character(method), source = "S7g rows")]
  EX <- rbind(EX, EG, fill = TRUE)
}
fwrite(EX, "tables/V6_S9_const_scope_extra.csv"); print(EX)
