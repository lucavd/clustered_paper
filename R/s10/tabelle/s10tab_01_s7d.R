
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
# Script 1: S7d chunks (1,896 design points; GEE fewclust strategies + dataset universe)
FL <- list.files(file.path(REPO, "results_v4/S7d/rows"), pattern = "\\.rds$", full.names = TRUE)
cat("files", length(FL), "\n"); t0 <- Sys.time()
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
one <- function(f) {
  x <- readRDS(f); R <- recs_std(x$recs); W <- as.data.table(x$rows); W[, scen := norm_scen(scen)]
  W <- merge(W, R[, c(K6, "estim", "const"), with = FALSE], by = K6, all.x = TRUE)
  out <- list(); kc <- NULL
  for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) {
    if (ty == "cr2" && cn == "ar1") next
    X <- W[cor == cn]; se <- X[[SEC[[ty]]]]
    use <- (X$estim %in% TRUE) & !(X$const %in% TRUE) & wald_ok(X$gee_ok, X$est, se)
    rr <- zrej(use, X$est, se, X$N - P_OF[X$scen])
    s <- paste0("gee_", cn, "_", ty)
    out[[s]] <- data.table(strategy = s, n = sum(use), kn = sum(rr$rn, na.rm = TRUE), kt = sum(rr$rt, na.rm = TRUE), nrow = nrow(X))
    if (s == "gee_exch_kc") kc <- data.table(X[, K6, with = FALSE], use = use, rn = rr$rn %in% TRUE)
  }
  key <- R[1, KD, with = FALSE]
  list(dp = cbind(key, rbindlist(out), nd = sum(R$estim), nconst = sum(R$const), nrecs = nrow(R)),
       recs = R, kc = if (R$scen[1] %in% c("M1", "M2", "M3", "M4", "M5")) kc else NULL)
}
L <- parallel::mclapply(FL, function(f) tryCatch(one(f), error = function(e) list(err = paste(f, conditionMessage(e)))), mc.cores = NCORES)
err <- unlist(lapply(L, `[[`, "err")); cat("errors", length(err), "\n"); if (length(err)) print(head(err))
DP <- rbindlist(lapply(L, `[[`, "dp")); RECS <- rbindlist(lapply(L, `[[`, "recs")); KC <- rbindlist(lapply(L, `[[`, "kc"))
rm(L); invisible(gc())
cat(sprintf("read in %.1f min; recs %d unique %d; dp rows %d\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(RECS), uniqueN(RECS[, ..K6]), nrow(DP)))
saveRDS(DP, file.path(OUT, "s1_gee_dp.rds")); saveRDS(RECS, file.path(OUT, "s1_recs.rds")); saveRDS(KC, file.path(OUT, "s1_kc_rep.rds"))
print(RECS[, .(n = .N, estim = sum(estim), degenerate = sum(!estim), const = sum(const)), by = .(dgp, grp = fifelse(scen %in% c("M6","M7"), "M6M7", fifelse(scen == "M1c", "M1c", "main")), label)])
cat("done\n")
