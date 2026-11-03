# Provenance (2026-11-02): read-only scope check run from a scratch directory (writes to ./out); see docs/steps/S9.md.

# s9_const_scope.R -- S9, read-only scope check before pre-registration: for every source of saved per-replicate rows,
# how many results counted as usable come from datasets with a constant outcome (sum_y == 0 or sum_y == n_obs), by
# method and arm. Constant flag from the S7d recs (1,896 points, same seeds as the study rows) and the S7e recs (copula).
# Writes ./out/scope_*.csv in the working directory. Read-only on results_v4/.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L)
REPO <- "/home/user/MMM"; dir.create("out", showWarnings = FALSE)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
RC0 <- c("sim_id", "generated", "design_degenerate", "n_obs", "sum_y", "dgp", "scen", "N", "ntime", "label")
fd <- list.files(file.path(REPO, "results_v4/S7d/rows"), full.names = TRUE, pattern = "\\.rds$")
fe <- list.files(file.path(REPO, "results_v4/S7e/rows"), full.names = TRUE, pattern = "\\.rds$")
RECD <- rbindlist(parallel::mclapply(fd, function(f) as.data.table(readRDS(f)$recs)[, ..RC0], mc.cores = 6L))
LE <- parallel::mclapply(fe, function(f) { x <- readRDS(f)
  list(recs = as.data.table(x$recs)[, ..RC0], st = as.data.table(x$study)) }, mc.cores = 6L)
RECE <- rbindlist(lapply(LE, `[[`, "recs")); STE <- rbindlist(lapply(LE, `[[`, "st"), fill = TRUE); rm(LE); invisible(gc())
cat("S7d dgp:", unique(RECD$dgp), "| S7e dgp:", unique(RECE$dgp), "\n")
RECA <- rbind(RECD, RECE); RECA[, const := generated & !design_degenerate & (sum_y == 0 | sum_y == n_obs)]
CK <- RECA[, c(K6, "design_degenerate", "const", "sum_y", "n_obs"), with = FALSE]
cat(sprintf("recs %d unique %d const %d\n", nrow(RECA), uniqueN(RECA[, ..K6]), RECA[, sum(const)]))
fwrite(RECA[, .(recs = .N, const = sum(const), all0 = sum(const & sum_y == 0)), by = .(dgp, label)], "out/scope_recs.csv")
cnt <- function(X, src) { X <- merge(X, CK, by = K6, all.x = TRUE)
  cat(src, ": rows", nrow(X), "unmatched keys", X[is.na(const), .N], "\n")
  X[, .(rows = .N, usable = sum(use, na.rm = TRUE), usable_const = sum(use & const %in% TRUE, na.rm = TRUE),
        const_rows = sum(const %in% TRUE), unmatched = sum(is.na(const))), by = .(method, term, label)][, src := src] }
OUT <- list()
## ---- S4 chain, variant kcfc (main grid + M1c) ----
S <- readRDS(file.path(REPO, "results_v4/variant_kcfc/S4/rows_v4.rds"))
cat("rows_v4 cols:", names(S), "\n"); print(S[, .N, by = .(method, term)])
S <- as.data.table(S)[, intersect(c(K6, "method", "term", "U0"), names(S)), with = FALSE]
S[, use := U0 %in% TRUE]; S[, U0 := NULL]
OUT[["kcfc"]] <- cnt(S, "variant_kcfc/S4/rows_v4"); rm(S); invisible(gc())
## ---- S6b ns rows (M6-M7) ----
S <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds")); cat("rows_ns cols:", names(S), "\n")
S <- as.data.table(S)[, intersect(c(K6, "method", "term", "U0"), names(S)), with = FALSE]
S[, use := U0 %in% TRUE]; S[, U0 := NULL]
OUT[["ns"]] <- cnt(S, "S6b/ns/rows_ns"); rm(S); invisible(gc())
## ---- S7e study rows (copula) ----
cat("S7e study cols:", names(STE), "\n")
u0 <- function(converged, estimate, se, p) { w <- converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p); w[is.na(w)] <- FALSE; w }
STE[, use := u0(converged, estimate, se, p)]
OUT[["s7e"]] <- cnt(STE[, c(K6, "method", "term", "use"), with = FALSE], "S7e/rows$study"); rm(STE)
## ---- S7g adaptive quadrature rows ----
fg <- list.files(file.path(REPO, "results_v4/S7g/rows"), full.names = TRUE, pattern = "\\.rds$")
x <- readRDS(fg[1]); cat("S7g file names:", names(x), "\n"); str(x, max.level = 1)
fwrite(rbindlist(OUT, fill = TRUE), "out/scope_usable_const.csv")
print(rbindlist(OUT, fill = TRUE)[usable_const > 0 | unmatched > 0][order(src, method, term, label)], nrows = 300)
