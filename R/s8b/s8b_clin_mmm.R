# =============================================================================
# s8b_clin_mmm.R -- S8b, D-70: the marginalized multilevel model of the two clinical applications refitted with a latent
# variable (lv.formula = ~1, t.formula = ~1, q = 20) on the same subsamples, keeping both the model-based and the robust SE.
# Cohort: data preparation, draws (rev2/T3_draws.rds, seed 20260807), specs and hard timeout as rev2/T3_rerun.R; the only
# change is lv.formula = ~1, q = 20 (as the simulation and the trial). Trial: as R/s7b/s7b_resp_run.R (already with the
# latent variable), refitted to obtain the model-based SE, which the stored rows do not carry.
# Outputs: tables/V5_S8b_clin_mmm_fits.csv (app, spec, N, rep, estimate, se_mod, se_rob, status, code),
#          tables/V5_S8b_clin_mmm_ref.csv (full data), tables/V5_S8b_clin_mmm_checks.csv.
# Check: trial estimates and robust SE equal the stored rows (tables/V4_S7b_resp_fits.csv) to 1e-8.
# Run: S8B_CORES=2 Rscript --vanilla R/s8b/s8b_clin_mmm.R
# =============================================================================
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressMessages({ library(data.table); library(MMLB); library(dplyr); library(parallel) }); setDTthreads(1L)
setwd(path.expand("~/MMM")); NCORE <- as.integer(Sys.getenv("S8B_CORES", "2")); TO <- 120
hto <- function(expr, secs) { p <- mcparallel(expr); r <- mccollect(p, wait = FALSE, timeout = secs)
  if (is.null(r)) { tools::pskill(p$pid, tools::SIGKILL); mccollect(p, wait = FALSE); return(structure("timeout", class = "hard_timeout")) }; r[[1]] }
fit1 <- function(f, D, term) {
  m <- hto(try(mm(f, lv.formula = ~1, t.formula = ~1, data = D, id = id, q = 20, verbose = FALSE), silent = TRUE), TO)
  st <- if (inherits(m, "hard_timeout")) "timeout" else if (inherits(m, "try-error")) "error" else if (!term %in% names(m$beta)) "no_term" else "ok"
  if (st != "ok") return(data.table(estimate = NA_real_, se_mod = NA_real_, se_rob = NA_real_, status = st, code = NA_real_))
  sm <- try(sqrt(m$mod.cov[term, term]), silent = TRUE); sr <- try(sqrt(m$rob.cov[term, term]), silent = TRUE)
  data.table(estimate = unname(m$beta[term]), se_mod = if (inherits(sm, "try-error")) NA_real_ else as.numeric(sm),
             se_rob = if (inherits(sr, "try-error")) NA_real_ else as.numeric(sr), status = "ok",
             code = unname(m$control["convergence_code"]))
}
## ---- neonatal cohort (rev2/T3_rerun.R) ----
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along); d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
dL <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior", grepl("P$", as.character(section)) ~ "posterior", TRUE ~ "lateral"),
                         levels = c("anterior", "lateral", "posterior")), pt = factor(paste(id, tp, sep = "_")), id = droplevels(factor(id)))
FSL <- list(pieno = lus_severe ~ region + disease + invasive + age_weeks, ridotto = lus_severe ~ region, ridotto_age = lus_severe ~ region + age_weeks)
prepL <- function(s) as.data.frame(s %>% arrange(id) %>% mutate(id = droplevels(factor(id)), across(where(is.factor), droplevels)))
## ---- respiratory trial (R/s7b/s7b_resp_run.R) ----
data(respiratory, package = "geepack"); dR <- as.data.frame(respiratory)
dR$pid <- factor(paste(dR$center, dR$id, sep = "_")); dR$id <- as.integer(dR$pid)
dR$treat <- factor(as.character(dR$treat), levels = c("P", "A")); dR$sex <- factor(as.character(dR$sex)); dR$center <- factor(dR$center)
dR <- dR[order(dR$pid, dR$visit), ]
FSR <- list(pieno = outcome ~ treat + baseline + visit + age + sex + center, ridotto = outcome ~ treat + baseline)
prepR <- function(s) as.data.frame(s %>% arrange(id, visit) %>% mutate(id = as.integer(id), pid = droplevels(factor(pid)), across(where(is.factor), droplevels)))
t0 <- Sys.time()
REF <- rbind(rbindlist(lapply(names(FSL), function(fn) fit1(FSL[[fn]], prepL(dL), "regionposterior")[, `:=`(app = "neonatal cohort", spec = fn)])),
             rbindlist(lapply(names(FSR), function(fn) fit1(FSR[[fn]], prepR(dR), "treatA")[, `:=`(app = "respiratory", spec = fn)])))
fwrite(REF, "tables/V5_S8b_clin_mmm_ref.csv"); print(REF, digits = 6)
drL <- readRDS("rev2/T3_draws.rds"); drR <- readRDS("results_v4/S7b/resp_draws.rds")
JOBS <- c(unlist(lapply(names(FSL), function(fn) lapply(if (fn == "ridotto_age") which(sapply(drL, function(z) z$N == 10)) else seq_along(drL),
                                                        function(i) list(app = "neonatal cohort", spec = fn, i = i))), recursive = FALSE),
          unlist(lapply(names(FSR), function(fn) lapply(seq_along(drR), function(i) list(app = "respiratory", spec = fn, i = i))), recursive = FALSE))
cat("jobs:", length(JOBS), "\n")
res <- mclapply(JOBS, function(j) {
  if (j$app == "neonatal cohort") { z <- drL[[j$i]]; r <- fit1(FSL[[j$spec]], prepL(dL[dL$id %in% z$sel, ]), "regionposterior") }
  else { z <- drR[[j$i]]; r <- fit1(FSR[[j$spec]], prepR(dR[dR$pid %in% z$sel, ]), "treatA") }
  r[, `:=`(app = j$app, spec = j$spec, N = z$N, rep = z$rep)] }, mc.cores = NCORE, mc.preschedule = TRUE)
bad <- sapply(res, function(x) !is.data.table(x)); cat("failed children:", sum(bad), "\n"); stopifnot(!any(bad))
X <- rbindlist(res); setcolorder(X, c("app", "spec", "N", "rep")); fwrite(X, "tables/V5_S8b_clin_mmm_fits.csv")
cat(sprintf("fits: %d in %.1f min\n", nrow(X), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(X[, .(n = .N, ok = sum(status == "ok"), timeout = sum(status == "timeout")), by = .(app, spec)])
## ---- checks ----
ST <- fread("tables/V4_S7b_resp_fits.csv")[method == "MMM", .(spec, N, rep, est_st = estimate, se_st = se, status_st = status)]
J <- merge(X[app == "respiratory"], ST, by = c("spec", "N", "rep"))
both <- J[is.finite(estimate) & is.finite(est_st)]
RS <- fread("tables/V4_S7b_resp_ref.csv")[method == "MMM"]
CK <- data.table(check = c("trial: rows matched to stored", "trial: status equal to stored", "trial: max |estimate - stored|", "trial: max |robust SE - stored SE|",
                           "trial full data: |estimate - stored|", "cohort full data (pieno): estimate"),
  value = c(nrow(J), J[, sum(status == status_st)], both[, max(abs(estimate - est_st))], both[is.finite(se_rob) & is.finite(se_st), max(abs(se_rob - se_st))],
            abs(REF[app == "respiratory" & spec == "pieno", estimate] - RS[spec == "pieno", estimate]), REF[app == "neonatal cohort" & spec == "pieno", estimate]))
CK[, pass := c(value[1] == 5000, value[2] == 5000, value[3] < 1e-8, value[4] < 1e-8, value[5] < 1e-8, round(value[6], 3) == 0.352)]
fwrite(CK, "tables/V5_S8b_clin_mmm_checks.csv"); print(CK, digits = 6)
