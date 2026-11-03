# s10t_extract.R -- S10 track "Testo": read-only extraction of the per-replicate rows into compact tables.
# Inputs (read-only): results_v4/S7d/rows, results_v4/S7e/rows, results_v4/variant_cc/S4/rows_v4.rds,
#   results_v4/S6b/ns/rows_ns.rds, results_v4/S6cc/rows. No project aggregation function is used.
# Output: results_v4/S10/testo/x_*.qs2 (compact copies). At most 3 processes.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"; dir.create(OUT, showWarnings = FALSE)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
rd <- function(dir) {
  f <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  L <- parallel::mclapply(f, function(x) { y <- readRDS(x)
    rc <- y$recs[, c(K6, "generated", "design_degenerate", "n_obs", "sum_y", "p_flip"), with = FALSE]
    g <- if (!is.null(y$rows)) y$rows else y$gee
    g <- if (!is.null(g) && nrow(g)) g[, c(K6, "cor", "gee_ok", "est", "se_unc", "se_md", "se_kc", "se_fg", "se_cr2"), with = FALSE] else NULL
    st <- if (!is.null(y$study) && nrow(y$study)) y$study[, c(K6, "method", "term", "estimate", "se", "p", "converged", "status"), with = FALSE] else NULL
    cc <- if (!is.null(y$cc) && nrow(y$cc)) y$cc else NULL
    list(rc = rc, g = g, st = st, cc = cc) }, mc.cores = 3L)
  list(rc = rbindlist(lapply(L, `[[`, "rc")), g = rbindlist(lapply(L, `[[`, "g")),
       st = rbindlist(lapply(L, `[[`, "st"), fill = TRUE), cc = rbindlist(lapply(L, `[[`, "cc"), fill = TRUE))
}
A <- rd("results_v4/S7d/rows"); tm(sprintf("S7d: recs %d, gee %d", nrow(A$rc), nrow(A$g)))
B <- rd("results_v4/S7e/rows"); tm(sprintf("S7e: recs %d, gee %d, study %d, cc %d", nrow(B$rc), nrow(B$g), nrow(B$st), nrow(B$cc)))
RECS <- rbind(A$rc, B$rc); GEE <- rbind(A$g, B$g)
qs_save(RECS, file.path(OUT, "x_recs.qs2"), nthreads = 2L); qs_save(GEE, file.path(OUT, "x_gee.qs2"), nthreads = 2L)
qs_save(B$st, file.path(OUT, "x_cop_study.qs2"), nthreads = 2L); qs_save(B$cc, file.path(OUT, "x_cop_cc.qs2"), nthreads = 2L)
cat("copula cc columns:", paste(names(B$cc), collapse = ","), "\n")
rm(A, B, RECS, GEE); invisible(gc()); tm("saved recs/gee/copula")
keep <- c(K6, "method", "term", "estimate", "se", "p", "converged", "status", "U0")
SA <- readRDS("results_v4/variant_cc/S4/rows_v4.rds"); tm(sprintf("rows_v4 read: %d rows; cols %s", nrow(SA), paste(names(SA), collapse = ",")))
print(SA[, .N, by = .(method, has_term = !is.na(term))])
SA <- SA[method %in% c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt", "mmm_nocc") & (is.na(term) | term == "time:trt")]
for (nm in setdiff(keep, names(SA))) SA[, (nm) := NA]
SA <- SA[, ..keep]; invisible(gc())
SB <- readRDS("results_v4/S6b/ns/rows_ns.rds")
SB <- SB[method %in% c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt", "mmm_nocc") & (is.na(term) | term == "time:trt")]
for (nm in setdiff(keep, names(SB))) SB[, (nm) := NA]
ST <- rbind(SA[, ..keep], SB[, ..keep]); rm(SA, SB); invisible(gc())
print(ST[, .N, by = .(method, dgp)][order(method, dgp)])
qs_save(ST, file.path(OUT, "x_study.qs2"), nthreads = 2L); rm(ST); invisible(gc()); tm("saved study rows")
f <- list.files("results_v4/S6cc/rows", pattern = "\\.rds$", full.names = TRUE)
CC <- rbindlist(parallel::mclapply(f, readRDS, mc.cores = 3L), fill = TRUE)
CC[, scen := paste0("M", model)]
qs_save(CC, file.path(OUT, "x_s6cc.qs2"), nthreads = 2L); tm(sprintf("S6cc: %d rows; cols %s", nrow(CC), paste(names(CC), collapse = ",")))
f <- list.files("results_v4/S6b/ns/rows", full.names = TRUE); cat("S6b ns rows files:", length(f), " e.g.", head(basename(f), 2), "\n")
y <- readRDS(f[1]); cat("S6b ns chunk structure:\n"); str(y, max.level = 2, list.len = 40)
tm("done")
