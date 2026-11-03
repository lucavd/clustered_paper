#!/usr/bin/env Rscript
# =============================================================================
# s7f_kc_clinical.R -- S7f: the three small-sample corrections of fewclust 0.2.0 on the subsamples of the two clinical
# applications, so that the Kauermann-Carroll rows of their tables use the published estimator (D-51), as in the
# simulation. Same data preparation, same draws, same GEE call as rev2/T3_rerun.R (neonatal cohort) and
# R/s7b/s7b_resp_run.R (respiratory trial); only the variance of the GEE estimate is recomputed.
# Mancl-DeRouen and Fay-Graubard are recomputed too, as a check against the stored geesmv values.
# Usage: Rscript --vanilla R/s7f/s7f_kc_clinical.R <lus|resp>
# Output: tables/V5_S7f_<app>_fc.csv (one row per spec, draw and correction; draw 0 = full data). Read-only on inputs.
# =============================================================================
.libPaths(c("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu", path.expand("~/R/s7d_lib")))
suppressMessages({ library(dplyr); library(geepack); library(data.table); library(parallel); library(fewclust) })
setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
setwd("/home/user/MMM")
APP <- commandArgs(TRUE)[1]; stopifnot(APP %in% c("lus", "resp"))
NCORE <- 20L
hto <- function(expr, secs) {
  p <- parallel::mcparallel(expr); r <- parallel::mccollect(p, wait = FALSE, timeout = secs)
  if (is.null(r)) { tools::pskill(p$pid, tools::SIGKILL); parallel::mccollect(p, wait = FALSE); return(structure("timeout", class = "hard_timeout")) }
  r[[1]]
}
if (APP == "lus") {
  TERM <- "regionposterior"
  d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
  d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along)
  d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
  d <- d0 %>% filter(!is.na(lus_severe)) %>%
    mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior",
                                     grepl("P$", as.character(section)) ~ "posterior", TRUE ~ "lateral"),
                           levels = c("anterior", "lateral", "posterior")),
           pt = factor(paste(id, tp, sep = "_")), id = droplevels(factor(id)))
  FS <- list(pieno = lus_severe ~ region + disease + invasive + age_weeks, ridotto = lus_severe ~ region,
             ridotto_age = lus_severe ~ region + age_weeks)
  draws <- readRDS("rev2/T3_draws.rds")
  prep <- function(dat) as.data.frame(dat %>% arrange(id) %>% mutate(id = droplevels(factor(id)), across(where(is.factor), droplevels)))
  sel_of <- function(z) d[d$id %in% z$sel, ]
  idx_of <- function(fn) if (fn == "ridotto_age") which(sapply(draws, function(z) z$N == 10)) else seq_along(draws)
} else {
  TERM <- "treatA"
  data(respiratory, package = "geepack")
  d <- as.data.frame(respiratory)
  d$pid <- factor(paste(d$center, d$id, sep = "_")); d$id <- as.integer(d$pid)
  d$treat <- factor(as.character(d$treat), levels = c("P", "A"))
  d$sex <- factor(as.character(d$sex)); d$center <- factor(d$center)
  d <- d[order(d$pid, d$visit), ]
  stopifnot(nrow(d) == 444L, nlevels(d$pid) == 111L)
  FS <- list(pieno = outcome ~ treat + baseline + visit + age + sex + center, ridotto = outcome ~ treat + baseline)
  draws <- readRDS("results_v4/S7b/resp_draws.rds")
  prep <- function(dat) as.data.frame(dat %>% arrange(id, visit) %>% mutate(id = as.integer(id), pid = droplevels(factor(pid)), across(where(is.factor), droplevels)))
  sel_of <- function(z) d[d$pid %in% z$sel, ]
  idx_of <- function(fn) seq_along(draws)
}
one <- function(dat, f) {
  D <- prep(dat)
  out <- data.table(corr = c("md", "kc", "fg"), gee_est = NA_real_, se = NA_real_, status = NA_character_)
  g <- hto(try(geeglm(f, family = binomial, id = id, corstr = "exchangeable", data = D), silent = TRUE), 60)
  if (inherits(g, "hard_timeout")) { out$status <- "no_base_gee"; return(out) }
  if (inherits(g, "try-error")) { out$status <- "no_base_gee"; return(out) }
  k <- match(TERM, names(coef(g)))
  if (is.na(k)) { out$status <- "no_base_gee"; return(out) }
  out$gee_est <- unname(coef(g)[k])
  for (j in seq_len(3)) {
    v <- hto(try(fewclust::gee_vcov(g, type = out$corr[j]), silent = TRUE), 30)
    if (inherits(v, "hard_timeout")) { out$status[j] <- "timeout"; next }
    if (inherits(v, "try-error")) { out$status[j] <- "error"; next }
    vv <- tryCatch(as.matrix(v)[k, k], error = function(e) NA_real_)
    if (!is.finite(vv) || vv <= 0) out$status[j] <- "bad_var" else { out$se[j] <- sqrt(vv); out$status[j] <- "ok" }
  }
  out
}
t0 <- Sys.time()
REF <- rbindlist(lapply(names(FS), function(fn) one(d, FS[[fn]])[, `:=`(spec = fn, N = NA_integer_, rep = 0L)]))
RES <- rbindlist(lapply(names(FS), function(fn) {
  res <- mclapply(idx_of(fn), function(i) { z <- draws[[i]]
    r <- try(one(sel_of(z), FS[[fn]]), silent = TRUE)
    if (inherits(r, "try-error")) r <- data.table(corr = c("md", "kc", "fg"), gee_est = NA_real_, se = NA_real_, status = "outer_error")
    r[, `:=`(spec = fn, N = z$N, rep = z$rep)] }, mc.cores = NCORE, mc.preschedule = TRUE)
  bad <- sapply(res, function(x) !is.data.table(x)); cat(sprintf("%s spec=%s: failed children=%d\n", APP, fn, sum(bad)))
  rbindlist(res[!bad]) }))
OUT <- rbind(REF, RES)
fwrite(OUT, sprintf("tables/V5_S7f_%s_fc.csv", APP))
cat(sprintf("%s: %d rows in %.1f min\n", APP, nrow(OUT), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(dcast(OUT[rep > 0, .N, by = .(spec, corr, status)], spec + corr ~ status, value.var = "N", fill = 0))
