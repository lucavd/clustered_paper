
# s10s_build.R -- S10 track "supplemento": compact per-replicate inputs for an independent recomputation.
# Reads only per-replicate row files (read-only); writes to results_v4/S10/supplemento/.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(parallel)}); setDTthreads(1)
RP <- "/home/user/MMM/results_v4"; OUT <- file.path(RP, "S10/supplemento"); dir.create(OUT, showWarnings = FALSE)
NC <- 3L; t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%6.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
K <- c("dgp","scen","N","ntime","label","sim_id")
rb <- function(L, nm) rbindlist(lapply(L, `[[`, nm), use.names = TRUE, fill = TRUE)
# ---- S7d: recs + GEE (three study mechanisms, 1,896 design points)
f <- list.files(file.path(RP, "S7d/rows"), full.names = TRUE); tm(paste("S7d files", length(f)))
L <- mclapply(f, function(x) { z <- readRDS(x)
  list(recs = z$recs[, c(K, "generated","design_degenerate","n_obs","sum_y","p_flip","p_flip_strict","flip_why","flip_strict_why"), with = FALSE],
       gee = z$rows[, c(K, "cor","gee_ok","est","se_unc","se_md","se_kc","se_fg","se_cr2","smv_md","smv_kc","smv_fg","smv_msg_md","smv_msg_fg"), with = FALSE]) }, mc.cores = NC)
saveRDS(rb(L, "recs"), file.path(OUT, "recs_s7d.rds")); saveRDS(rb(L, "gee"), file.path(OUT, "gee_s7d.rds")); rm(L); gc(); tm("S7d done")
# ---- S7e: Copula-AR1 (480 design points)
f <- list.files(file.path(RP, "S7e/rows"), full.names = TRUE); tm(paste("S7e files", length(f)))
L <- mclapply(f, function(x) { z <- readRDS(x)
  list(recs = z$recs[, c(K, "generated","design_degenerate","n_obs","sum_y","p_flip","p_flip_strict","flip_why","flip_strict_why"), with = FALSE],
       gee = z$gee[, c(K, "cor","gee_ok","est","se_unc","se_md","se_kc","se_fg","se_cr2","smv_md","smv_kc","smv_fg","smv_msg_md","smv_msg_fg"), with = FALSE],
       study = z$study, cc = if (nrow(z$cc)) z$cc else NULL) }, mc.cores = NC)
saveRDS(rb(L, "recs"), file.path(OUT, "recs_s7e.rds")); saveRDS(rb(L, "gee"), file.path(OUT, "gee_s7e.rds"))
saveRDS(rb(L, "study"), file.path(OUT, "study_s7e.rds")); saveRDS(rb(L, "cc"), file.path(OUT, "cc_s7e.rds")); rm(L); gc(); tm("S7e done")
# ---- S7d2: inverted sign-flip interval (all 2,376 design points)
f <- list.files(file.path(RP, "S7d2/rows"), full.names = TRUE); tm(paste("S7d2 files", length(f)))
L <- mclapply(f, function(x) { z <- as.data.table(readRDS(x))
  z[, c(K, "generated","design_degenerate","n_obs","sum_y","theta","p_theta","p0_guard","why0_guard","lo","hi","ar1_gee_ok","ar1_est","ar1_se_fg"), with = FALSE] }, mc.cores = NC)
saveRDS(rbindlist(L, fill = TRUE), file.path(OUT, "ci_s7d2.rds")); rm(L); gc(); tm("S7d2 done")
# ---- S7g: adaptive quadrature
f <- list.files(file.path(RP, "S7g/rows"), full.names = TRUE); tm(paste("S7g files", length(f)))
L <- mclapply(f, function(x) readRDS(x), mc.cores = NC)
saveRDS(rb(L, "recs"), file.path(OUT, "recs_s7g.rds")); saveRDS(rb(L, "rows"), file.path(OUT, "rows_s7g.rds")); rm(L); gc(); tm("S7g done")
# ---- study rows (variant_cc), main grid + M1c, and the study MMM without the check (results_v4/S4)
V <- readRDS(file.path(RP, "variant_cc/S4/rows_v4.rds")); setDT(V)
V <- V[, c(K, "method","estimate","se","p","converged","status","U0","U1"), with = FALSE]
saveRDS(V, file.path(OUT, "study_cc.rds"), compress = FALSE); rm(V); gc(); tm("variant_cc rows done")
S <- readRDS(file.path(RP, "S4/rows_v4.rds")); setDT(S)
saveRDS(S[method == "mmm", c(K, "method","estimate","se","p","converged","status","U0","U1"), with = FALSE], file.path(OUT, "study_mmm_nocc.rds"), compress = FALSE)
rm(S); gc(); tm("S4 mmm rows done")
cat("BUILD OK\n")
