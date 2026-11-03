## S10 track 'assolute': per-replicate facts needed by the absolute-sentence checks.
## Reads every chunk of results_v4/S7d/rows (37,920 files) and results_v4/S7e/rows
## (9,600 files); read-only. Writes <outdir>/recs_all.rds and fgpt_rows.rds.
## Usage: Rscript --vanilla s10_assolute_recs.R <outdir> <cores>
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
args <- commandArgs(TRUE); out <- args[1]; nc <- as.integer(args[2])
dir.create(out, showWarnings = FALSE, recursive = TRUE)
one <- function(f, run) {
  x <- readRDS(f)
  r <- as.data.table(x$recs)[, .(dgp, scen, N, ntime, label, sim_id, generated,
         design_degenerate, n_obs, sum_y, p_flip, flip_why)]
  g <- if (run == "S7d") x$rows else x$gee
  gc <- as.data.table(g)[, .(n_gee_cor = uniqueN(cor), gee_ok_any = any(gee_ok %in% TRUE)), by = sim_id]
  r <- merge(r, gc, by = "sim_id", all.x = TRUE)
  r[is.na(n_gee_cor), n_gee_cor := 0L]
  if (run == "S7e") {
    s <- as.data.table(x$study)
    sc <- dcast(s[, .N, by = .(sim_id, method)], sim_id ~ method, value.var = "N", fill = 0L)
    setnames(sc, setdiff(names(sc), "sim_id"), paste0("n_", setdiff(names(sc), "sim_id")))
    r <- merge(r, sc, by = "sim_id", all.x = TRUE)
  }
  r[, run := run]; r[, file := basename(f)]
  keep <- NULL
  if (grepl("_m2_N10_T2_manual_", basename(f)) && run == "S7d")
    keep <- as.data.table(g)[cor == "ar1"][, file := basename(f)]
  list(r = r, keep = keep)
}
fs <- c(list.files("/home/user/MMM/results_v4/S7d/rows", full.names = TRUE),
        list.files("/home/user/MMM/results_v4/S7e/rows", full.names = TRUE))
run <- ifelse(grepl("/S7d/", fs), "S7d", "S7e")
cat(length(fs), "files\n")
t0 <- Sys.time()
res <- parallel::mclapply(seq_along(fs), function(i) tryCatch(one(fs[i], run[i]),
         error = function(e) list(r = data.table(file = basename(fs[i]), err = conditionMessage(e)), keep = NULL)),
         mc.cores = nc, mc.preschedule = TRUE)
R <- rbindlist(lapply(res, `[[`, "r"), fill = TRUE)
K <- rbindlist(lapply(res, `[[`, "keep"), fill = TRUE)
saveRDS(R, file.path(out, "recs_all.rds")); saveRDS(K, file.path(out, "fgpt_rows.rds"))
cat("rows", nrow(R), "errors", if ("err" %in% names(R)) sum(!is.na(R$err)) else 0, "secs", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
