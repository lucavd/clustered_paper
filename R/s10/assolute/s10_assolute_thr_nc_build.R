## S10 track 'assolute', follow-up on A0315: row sets t4nc and t2nc = the current row set (variant_nc: fewclust KC and
## the constant-outcome rule D-82) with the MMM rows of M3-M5 taken from variant_t4 / variant_t2 (threshold 1e-4 / 1e-2
## of the convergence check, built by R/s6b/s6b_thr_build.R from fits that already exist). No model is refitted.
## Read-only on results_v4/variant_{nc,t4,t2}; writes only under results_v4/S10/assolute/thr_shadow and the gate CSV.
## Usage: Rscript --vanilla s10_assolute_thr_nc_build.R
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
source("/home/user/MMM/R/s9/s9_const.R")
RV <- "/home/user/MMM/results_v4"; SH <- file.path(RV, "S10/assolute/thr_shadow/MMM")
t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
NC <- readRDS(file.path(RV, "variant_nc/S4/rows_v4.rds")); cat("nc loaded", nrow(NC), el(), "min\n")
K8 <- c("dgp", "scen", "N", "ntime", "label", "sim_id", "method", "term")
sel <- function(V) V[, which(method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term) & term == "time:trt")]
iN <- sel(NC)
G <- list()
for (tag in c("t4", "t2")) {
  TT <- readRDS(file.path(RV, sprintf("variant_%s/S4/rows_v4.rds", tag))); cat(tag, "loaded", nrow(TT), el(), "min\n")
  only_nc <- setdiff(names(NC), names(TT)); only_t <- setdiff(names(TT), names(NC))
  iT <- sel(TT)
  A <- NC[iN, ..K8][, iN := iN]; B <- TT[iT, ..K8][, iT := iT]
  stopifnot(!anyDuplicated(A, by = K8), !anyDuplicated(B, by = K8))
  m <- merge(A, B, by = K8); stopifnot(nrow(m) == length(iN), nrow(m) == length(iT))
  repl <- setdiff(intersect(names(NC), names(TT)), K8)
  ndiff <- sapply(repl, function(cl) { a <- NC[[cl]][m$iN]; b <- TT[[cl]][m$iT]; sum(!((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))) })
  W <- copy(NC)
  for (cl in repl) set(W, i = m$iN, j = cl, value = TT[[cl]][m$iT])
  # rows outside the MMM M3-M5 subset: how does the old threshold row set differ from the current one?
  same_order <- isTRUE(all.equal(NC[, ..K8], TT[, ..K8], check.attributes = FALSE))
  out_diff <- if (same_order) {
    o <- setdiff(seq_len(nrow(NC)), iN)
    d <- Reduce(`|`, lapply(intersect(c("estimate", "se", "p", "U0", "U1", "U2"), repl), function(cl) { a <- NC[[cl]][o]; b <- TT[[cl]][o]; !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)) }))
    x <- NC[o[d], .N, by = method]; paste(sprintf("%s:%d", x$method, x$N), collapse = " ")
  } else "row order differs"
  rm(TT, A, B); gc(verbose = FALSE)
  n0 <- s9_drop_const(W, "U0"); n1 <- s9_drop_const(W, "U1"); n2 <- s9_drop_const(W, "U2")   # D-82
  od <- file.path(SH, sprintf("results_v4/variant_%snc/S4", tag)); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  saveRDS(W, file.path(od, "rows_v4.rds"), compress = FALSE)
  G[[tag]] <- data.table(variant = paste0(tag, "nc"), rows = nrow(W), mmm_m3m5_rows = length(iN),
    mmm_rows_changed_by_threshold = paste(sprintf("%s:%d", names(ndiff)[ndiff > 0], ndiff[ndiff > 0]), collapse = " "),
    cols_only_in_nc = paste(only_nc, collapse = " "), cols_only_in_old_variant = paste(only_t, collapse = " "),
    rows_outside_subset_differing_old_vs_nc = out_diff, D82_switched_off = sprintf("U0 %d, U1 %d, U2 %d", n0, n1, n2), minutes = el())
  print(G[[tag]]); rm(W); gc(verbose = FALSE)
}
fwrite(rbindlist(G), file.path(RV, "S10/assolute/thr_nc_build_gates.csv"))
cat("done", el(), "min\n")
