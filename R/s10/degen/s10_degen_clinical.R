# =============================================================================
# s10_degen_clinical.R -- S10 "Degeneri", item 5: separation classes on the subsamples of the two clinical applications.
# Data preparation, draws, specifications and term exactly as R/s7f/s7f_kc_clinical.R and R/s8b/s8b_clin_mmm.R
# (neonatal cohort: rev2/T3_draws.rds, term regionposterior; respiratory trial: results_v4/S7b/resp_draws.rds, term treatA).
# Per draw: S = the ML estimate of the logistic mean model of the GEE (formula of the spec, independence) does not exist
# (lp_sep of R/s10/probe/s10_probe2.R); Si = detectseparation 0.4.0 declares the coefficient of the term infinite.
# Usable flags as in R/s8b/s8b_clinical_tables.R (KC rows from fewclust, MMM rows from the latent-variable refit with the
# model-based SE; usable = finite estimate, 0 < SE, |estimate| <= 10, SE < 10; the tables call "sep" |est| > 10 or SE > 10).
# Output: tables/V7_S10_degen_clinical.csv, tables/V7_S10_degen_clinical_draws.csv. Read-only on every input.
# =============================================================================
.libPaths(c("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu", path.expand("~/R/s10_lib")))
suppressMessages({ library(data.table); library(dplyr); library(detectseparation); library(lpSolveAPI); library(parallel) })
setDTthreads(1L); setwd("/home/user/MMM")
local({ src <- readLines("R/s10/probe/s10_probe2.R"); a <- grep("^lp_sep <- function", src)
  b <- a + which(src[(a + 1):length(src)] == "}")[1]; eval(parse(text = src[a:b]), envir = globalenv()) })
## ---- data and draws (verbatim preparation) ----
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along); d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
dL <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior", grepl("P$", as.character(section)) ~ "posterior", TRUE ~ "lateral"),
                         levels = c("anterior", "lateral", "posterior")), pt = factor(paste(id, tp, sep = "_")), id = droplevels(factor(id)))
FSL <- list(pieno = lus_severe ~ region + disease + invasive + age_weeks, ridotto = lus_severe ~ region, ridotto_age = lus_severe ~ region + age_weeks)
prepL <- function(s) as.data.frame(s %>% arrange(id) %>% mutate(id = droplevels(factor(id)), across(where(is.factor), droplevels)))
data(respiratory, package = "geepack"); dR <- as.data.frame(respiratory)
dR$pid <- factor(paste(dR$center, dR$id, sep = "_")); dR$id <- as.integer(dR$pid)
dR$treat <- factor(as.character(dR$treat), levels = c("P", "A")); dR$sex <- factor(as.character(dR$sex)); dR$center <- factor(dR$center)
dR <- dR[order(dR$pid, dR$visit), ]; stopifnot(nrow(dR) == 444L)
FSR <- list(pieno = outcome ~ treat + baseline + visit + age + sex + center, ridotto = outcome ~ treat + baseline)
prepR <- function(s) as.data.frame(s %>% arrange(id, visit) %>% mutate(id = as.integer(id), pid = droplevels(factor(pid)), across(where(is.factor), droplevels)))
drL <- readRDS("rev2/T3_draws.rds"); drR <- readRDS("results_v4/S7b/resp_draws.rds")
JOBS <- c(unlist(lapply(names(FSL), function(fn) lapply(if (fn == "ridotto_age") which(sapply(drL, function(z) z$N == 10)) else seq_along(drL),
                                                      function(i) list(app = "neonatal cohort", spec = fn, i = i))), recursive = FALSE),
          unlist(lapply(names(FSR), function(fn) lapply(seq_along(drR), function(i) list(app = "respiratory", spec = fn, i = i))), recursive = FALSE))
one <- function(j) {
  if (j$app == "neonatal cohort") { z <- drL[[j$i]]; D <- prepL(dL[dL$id %in% z$sel, ]); f <- FSL[[j$spec]]; term <- "regionposterior" }
  else { z <- drR[[j$i]]; D <- prepR(dR[dR$pid %in% z$sel, ]); f <- FSR[[j$spec]]; term <- "treatA" }
  y <- model.response(model.frame(f, D)); X <- model.matrix(f, D)
  y <- as.numeric(y); a <- if (length(unique(y)) < 2) TRUE else lp_sep(X, y)
  out <- data.table(app = j$app, spec = j$spec, N = z$N, rep = z$rep, nobs = nrow(X), ncoef = ncol(X), rank = qr(X)$rank, sum_y = sum(y),
                    term_in_X = term %in% colnames(X), S = isTRUE(a), lp_na = is.na(a), Si = FALSE, ds_outcome = NA, inf_terms = NA_character_)
  if (isTRUE(a) || is.na(a)) {
    s <- tryCatch(suppressWarnings(glm(f, data = D, family = binomial, method = "detect_separation")), error = function(e) e)
    if (!inherits(s, "error")) { cf <- s$coefficients; out$ds_outcome <- isTRUE(s$outcome)
      out$Si <- isTRUE(a) && term %in% names(cf) && is.infinite(cf[[term]]); out$inf_terms <- paste(names(cf)[is.infinite(cf)], collapse = ";") }
  }
  out
}
t0 <- Sys.time()
one_safe <- function(j) tryCatch(one(j), error = function(e) { z <- if (j$app == "neonatal cohort") drL[[j$i]] else drR[[j$i]]
  data.table(app = j$app, spec = j$spec, N = z$N, rep = z$rep, nobs = NA_integer_, ncoef = NA_integer_, rank = NA_integer_, sum_y = NA_real_,
             term_in_X = NA, S = NA, lp_na = NA, Si = NA, ds_outcome = NA, inf_terms = NA_character_, err = substr(conditionMessage(e), 1, 120)) })
CL <- rbindlist(mclapply(JOBS, one_safe, mc.cores = 4L), fill = TRUE)
if (!"err" %in% names(CL)) CL[, err := NA_character_]
cat("draws not classifiable:", CL[!is.na(err), .N], "\n"); print(CL[!is.na(err), .N, by = .(app, spec, err)])
cat(sprintf("draws classified: %d in %.1f min\n", nrow(CL), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
## ---- usable flags as in R/s8b/s8b_clinical_tables.R ----
swap_kc <- function(A, FC) { K <- c("spec", "N", "rep"); kc <- FC[rep > 0 & corr == "kc"]
  B <- A[method == "GEE-KC"]; B <- merge(B, kc[, c(K, "se", "status"), with = FALSE], by = K, all.x = TRUE, suffixes = c("", ".fc"))
  stopifnot(B[is.na(status.fc), .N] == 0L)
  B[, `:=`(se = fifelse(status.fc == "ok", se.fc, NA_real_), status = status.fc)]; B[, c("se.fc", "status.fc") := NULL]
  rbind(A[method != "GEE-KC"], B[, names(A), with = FALSE]) }
MF <- fread("tables/V5_S8b_clin_mmm_fits.csv")
swap_mmm <- function(A, app_) { K <- c("spec", "N", "rep"); X <- MF[app == app_]; B <- A[method == "MMM"]; n0 <- nrow(B)
  B <- merge(B, X[, .(spec, N, rep, est_new = estimate, se_new = se_mod, st_new = status)], by = K, all.x = TRUE)
  stopifnot(B[is.na(st_new), .N] == 0L, nrow(B) == n0)
  B[, `:=`(estimate = est_new, se = se_new, status = st_new)][, c("est_new", "se_new", "st_new") := NULL]
  rbind(A[method != "MMM"], B[, names(A), with = FALSE]) }
AL <- swap_mmm(swap_kc(fread("rev2/T3_subsample_rerun.csv"), fread("tables/V5_S7f_lus_fc.csv")), "neonatal cohort")[, app := "neonatal cohort"]
AR <- swap_mmm(swap_kc(fread("tables/V4_S7b_resp_fits.csv"), fread("tables/V5_S7f_resp_fc.csv")), "respiratory")[, app := "respiratory"]
A <- rbind(AL[, .(app, method, spec, N, rep, estimate, se)], AR[, .(app, method, spec, N, rep, estimate, se)])
A[, sep_paper := (is.finite(estimate) & abs(estimate) > 10) | (is.finite(se) & se > 10)]
A[, usable := is.finite(estimate) & is.finite(se) & se > 0 & !sep_paper & abs(estimate) <= 10 & se < 10]
# cross-check against the stored clinical usable tables
chkL <- merge(A[app == "neonatal cohort", .(u = sum(usable)), by = .(spec, method, N)], fread("tables/V5_S7f_lus_tab316.csv")[, .(spec, method, N, usable)], by = c("spec", "method", "N"))
chkR <- merge(A[app == "respiratory", .(u = sum(usable)), by = .(spec, method, N)], fread("tables/V5_S7f_resp_usable.csv")[, .(spec, method, N, usable)], by = c("spec", "method", "N"))
cat(sprintf("usable counts equal to stored tables: cohort %d/%d, trial %d/%d\n", chkL[, sum(u == usable)], nrow(chkL), chkR[, sum(u == usable)], nrow(chkR)))
J <- merge(A, CL, by = c("app", "spec", "N", "rep"), all.x = TRUE); stopifnot(nrow(J) == nrow(A))
J[, `:=`(S_na = is.na(S))]; J[is.na(S), `:=`(S = FALSE, Si = FALSE)]
TAB <- rbind(
  CL[, .(method = "(draws)", draws = .N, S = sum(S, na.rm = TRUE), Si = sum(Si, na.rm = TRUE), not_classifiable = sum(is.na(S)), lp_na = sum(lp_na, na.rm = TRUE), usable = NA_integer_, usable_on_S = NA_integer_, usable_on_Si = NA_integer_,
         sep_paper = NA_integer_, sep_paper_on_S = NA_integer_, usable_on_unclassifiable = NA_integer_), by = .(app, spec, N)],
  J[, .(draws = .N, S = sum(S), Si = sum(Si), not_classifiable = sum(S_na), lp_na = sum(lp_na, na.rm = TRUE), usable = sum(usable), usable_on_unclassifiable = sum(usable & S_na), usable_on_S = sum(usable & S), usable_on_Si = sum(usable & Si),
        sep_paper = sum(sep_paper), sep_paper_on_S = sum(sep_paper & S)), by = .(app, spec, N, method)])
setcolorder(TAB, c("app", "spec", "N", "method")); setorder(TAB, app, spec, N, method)
fwrite(TAB, "tables/V7_S10_degen_clinical.csv"); fwrite(CL, "tables/V7_S10_degen_clinical_draws.csv")
print(TAB[method == "(draws)"]); print(J[, .(usable = sum(usable), on_S = sum(usable & S), on_Si = sum(usable & Si)), by = .(app, method)])
fwrite(rbind(chkL[, .(app = "neonatal cohort", cells = .N, equal = sum(u == usable))], chkR[, .(app = "respiratory", cells = .N, equal = sum(u == usable))]),
       "results_v4/S10/degen/clinical_usable_check.csv")
