
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
# Script 3: M6-M7 chunks (S6b/ns), Copula-AR1 chunks (S7e), sign-flip chunks (S7d2)
t0 <- Sys.time()
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
gee_counts <- function(W, R) {
  W <- as.data.table(W); W[, scen := norm_scen(scen)]
  W <- merge(W, R[, c("sim_id", "estim", "const"), with = FALSE], by = "sim_id", all.x = TRUE); out <- list()
  for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) { if (ty == "cr2" && cn == "ar1") next
    X <- W[cor == cn]; se <- X[[SEC[[ty]]]]
    use <- (X$estim %in% TRUE) & !(X$const %in% TRUE) & wald_ok(X$gee_ok, X$est, se); rr <- zrej(use, X$est, se, X$N - P_OF[X$scen])
    out[[paste(cn, ty)]] <- data.table(strategy = paste0("gee_", cn, "_", ty), n = sum(use), kn = sum(rr$rn, na.rm = TRUE), kt = sum(rr$rt, na.rm = TRUE), dup = 0L) }
  rbindlist(out) }
STM <- c("mmm", "glmer", "shrink", "glmer_lrt", "gee")
one_ns <- function(f) { x <- readRDS(f); R <- recs_std(cbind(as.data.table(x$recs), as.data.table(x$rows)[1, .(dgp, scen, N, ntime, label)]))
  cc <- cc_counts(x$cc, R); st <- study_counts(x$rows, R, c(STM, "mmm_nocc"))
  st[strategy == "mmm", strategy := "mmm_projrow"]
  list(dp = cbind(R[1, KD, with = FALSE], rbind(cc$dp, st), nd = sum(R$estim), nconst = sum(R$const), nrecs = nrow(R)), agree = cbind(R[1, KD, with = FALSE], cc$agree)) }
one_cop <- function(f) { x <- readRDS(f); R <- recs_std(x$recs); st <- study_counts(x$study, R, STM)
  if (R$scen[1] %in% c("M1", "M2")) { cc <- list(dp = copy(st[strategy == "mmm"])[, strategy := "mmm_cc"], agree = NULL)
    cc$dp <- rbind(cc$dp, copy(cc$dp)[, strategy := "mmm_cc_f"]) } else cc <- cc_counts(x$cc, R)
  st[strategy == "mmm", strategy := "mmm_projrow"]
  list(dp = cbind(R[1, KD, with = FALSE], rbind(cc$dp, st, gee_counts(x$gee, R)), nd = sum(R$estim), nconst = sum(R$const), nrecs = nrow(R)),
       agree = if (is.null(cc$agree)) NULL else cbind(R[1, KD, with = FALSE], cc$agree)) }
one_flip <- function(f) { x <- as.data.table(readRDS(f)); R <- recs_std(x)
  X <- merge(R, x[, .(sim_id, p0_guard, why0_guard)], by = "sim_id")
  use <- X$estim & !X$const & is.finite(X$p0_guard); rj <- use & X$p0_guard < 0.05
  list(dp = cbind(R[1, KD, with = FALSE], data.table(strategy = "flip", n = sum(use), kn = sum(rj), kt = sum(rj), dup = 0L), nd = sum(R$estim), nconst = sum(R$const), nrecs = nrow(R)),
       why = X[estim & !const, .N, by = .(why = substr(ifelse(is.na(why0_guard), "<NA>", why0_guard), 1, 40), finite_p = is.finite(p0_guard))],
       recs = R[, c(K6, "estim", "const", "n_obs", "sum_y"), with = FALSE]) }
run <- function(dir, fun) { FL <- list.files(file.path(REPO, dir), pattern = "\\.rds$", full.names = TRUE); cat(dir, length(FL), "files\n")
  L <- parallel::mclapply(FL, function(f) tryCatch(fun(f), error = function(e) list(err = paste(basename(f), conditionMessage(e)))), mc.cores = NCORES)
  err <- unlist(lapply(L, `[[`, "err")); cat("  errors", length(err), "\n"); if (length(err)) print(head(err)); L }
L <- run("results_v4/S6b/ns/rows", one_ns)
saveRDS(list(dp = rbindlist(lapply(L, `[[`, "dp")), agree = rbindlist(lapply(L, `[[`, "agree"))), file.path(OUT, "s3_ns.rds"))
cat(sprintf("ns done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
L <- run("results_v4/S7e/rows", one_cop)
saveRDS(list(dp = rbindlist(lapply(L, `[[`, "dp")), agree = rbindlist(lapply(L, `[[`, "agree"))), file.path(OUT, "s3_cop.rds"))
cat(sprintf("copula done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
L <- run("results_v4/S7d2/rows", one_flip)
FLIP <- list(dp = rbindlist(lapply(L, `[[`, "dp")), why = rbindlist(lapply(L, `[[`, "why"))[, .(N = sum(N)), by = .(why, finite_p)], recs = rbindlist(lapply(L, `[[`, "recs")))
saveRDS(FLIP, file.path(OUT, "s3_flip.rds")); print(FLIP$why[order(-N)])
cat(sprintf("flip done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
