# =============================================================================
# s7b_resp_run.R -- S7b part 3 (cluster column named `id` for geesmv, see below): second clinical application on geepack::respiratory, with the subsampling scheme and
# the strategies of the neonatal cohort (rev2/T3_rerun.R). Declared difference: the marginalized model is specified as in
# scenario M1 of the simulation (lv.formula = ~1, t.formula = ~1). Gate GC5: dataset and reproduction of rev2/T3_draws.rds.
# Outputs: results_v4/S7b/resp_draws.rds, tables/V4_S7b_resp_ref.csv, tables/V4_S7b_resp_fits.csv, tables/V4_S7b_resp_degeneracy.csv,
#          GC5 row in tables/V4_S7b_gates.csv
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R"))
suppressMessages({library(dplyr); library(lme4); library(geepack); library(MMLB); library(parallel)})
NCORE <- as.integer(Sys.getenv("S7B_NCORE", "20"))
TERM <- "treatA"
hto <- function(expr, secs) {   # hard timeout by fork + SIGKILL (as rev2/T3_rerun.R)
  p <- parallel::mcparallel(expr); r <- parallel::mccollect(p, wait = FALSE, timeout = secs)
  if (is.null(r)) { tools::pskill(p$pid, tools::SIGKILL); parallel::mccollect(p, wait = FALSE); return(structure("timeout", class = "hard_timeout")) }
  r[[1]]
}
TO <- c(mmm = 120, gee = 60, glmm = 120, smv = 5)
# ---- dataset -----------------------------------------------------------------------------------------------------------
data(respiratory, package = "geepack")
d <- as.data.frame(respiratory)
d$pid <- factor(paste(d$center, d$id, sep = "_"))
d$id <- as.integer(d$pid)   # the cluster column MUST be named `id`: geesmv keys on `data$id` whenever that column exists, and `respiratory` has a within-centre `id` (1..56) that would otherwise define the clusters
d$treat <- factor(as.character(d$treat), levels = c("P", "A"))
d$sex <- factor(as.character(d$sex)); d$center <- factor(d$center)
d <- d[order(d$pid, d$visit), ]
f <- tempfile(); write.csv(respiratory, f, row.names = FALSE); md5 <- unname(tools::md5sum(f))
cat("respiratory:", nrow(d), "rows,", nlevels(d$pid), "patients, visits per patient:", paste(names(table(table(d$pid))), collapse = "/"), "md5 csv", md5, "\n")
gc5_data <- nrow(d) == 444L && nlevels(d$pid) == 111L && all(table(d$pid) == 4L) && md5 == "281e0d3d43ddca923b9debb6e2e996f1"
# ---- GC5: the draw procedure reproduces rev2/T3_draws.rds on the cohort ids ------------------------------------------------
draw_ids <- function(ids, NLEV = c(10, 15, 20, 25, 50), REPS = 500, seed = 20260807) {
  set.seed(seed); draws <- list(); k <- 0L
  for (n in NLEV) for (r in seq_len(REPS)) { k <- k + 1L; draws[[k]] <- list(N = n, rep = r, sel = sample(ids, n)) }
  draws
}
d0 <- readRDS(file.path(REPO, "ARTICOLO 2026/LURID_bis/lurid_long.rds"))
dl <- d0 %>% filter(!is.na(lus_severe)) %>% mutate(id = droplevels(factor(id)))
old <- readRDS(file.path(REPO, "rev2/T3_draws.rds")); new <- draw_ids(levels(dl$id))
gc5_draws <- length(old) == length(new) && all(mapply(function(a, b) a$N == b$N && a$rep == b$rep && identical(as.character(a$sel), as.character(b$sel)), old, new))
cat("cohort draws reproduced:", gc5_draws, "\n")
GT <- fread(file.path(TABDIR, "V4_S7b_gates.csv"))[gate != "GC5"]
GT <- rbind(GT, data.table(gate = "GC5", pass = gc5_data && gc5_draws,
  detail = sprintf("respiratory: %d rows, %d patients, 4 visits each: %s; md5 of the CSV export %s; the 2,500 cohort draws (seed 20260807) reproduced: %s", nrow(d), nlevels(d$pid), all(table(d$pid) == 4L), md5, gc5_draws)))
setorder(GT, gate); fwrite(GT, file.path(TABDIR, "V4_S7b_gates.csv")); print(GT)
stopifnot(gc5_data, gc5_draws)
# ---- draws on respiratory -----------------------------------------------------------------------------------------------
draws <- draw_ids(levels(d$pid)); saveRDS(draws, file.path(OUT_S7B, "resp_draws.rds"))
one_level <- function(x) length(unique(as.character(x))) < 2
DG <- rbindlist(lapply(draws, function(z) { s <- d[d$pid %in% z$sel, ]
  data.table(N = z$N, rep = z$rep, nrows = nrow(s), nev = sum(s$outcome), one_treat = one_level(s$treat), one_sex = one_level(s$sex),
             one_center = one_level(s$center), one_baseline = one_level(s$baseline), n_treatA = length(unique(s$pid[s$treat == "A"]))) }))
DG[, any_one := one_treat | one_sex | one_center | one_baseline]
fwrite(DG, file.path(TABDIR, "V4_S7b_resp_degeneracy.csv"))
print(DG[, .(n = .N, pct_any_one = 100 * mean(any_one), pct_one_treat = 100 * mean(one_treat), pct_one_sex = 100 * mean(one_sex), pct_one_center = 100 * mean(one_center), pct_one_baseline = 100 * mean(one_baseline)), by = N])
# ---- fits (rev2/T3_rerun.R, MMM with lv.formula = ~1) ------------------------------------------------------------------
row1 <- function(meth, est, se, status) data.table(method = meth, estimate = est, se = se, status = status)
fit_set <- function(dat, f) {
  dat <- dat %>% arrange(id, visit) %>% mutate(id = as.integer(id), pid = droplevels(factor(pid)), across(where(is.factor), droplevels))
  D <- as.data.frame(dat); out <- list()
  m <- hto(try(mm(f, lv.formula = ~1, t.formula = ~1, data = D, id = id, q = 20, verbose = FALSE), silent = TRUE), TO["mmm"])
  if (inherits(m, "hard_timeout")) out$MMM <- row1("MMM", NA, NA, "timeout")
  else if (inherits(m, "try-error")) out$MMM <- row1("MMM", NA, NA, "error")
  else if (!TERM %in% names(m$beta)) out$MMM <- row1("MMM", NA, NA, "no_term")
  else { V <- if (!is.null(m$rob.cov)) m$rob.cov else m$mod.cov; se <- try(sqrt(V[TERM, TERM]), silent = TRUE)
         out$MMM <- row1("MMM", unname(m$beta[TERM]), if (inherits(se, "try-error")) NA_real_ else as.numeric(se), "ok") }
  g <- hto(try(geeglm(f, family = binomial, id = id, corstr = "exchangeable", data = D), silent = TRUE), TO["gee"])
  if (inherits(g, "hard_timeout")) { out$GEE <- row1("GEE", NA, NA, "timeout"); g <- NULL }
  else if (inherits(g, "try-error")) { out$GEE <- row1("GEE", NA, NA, "error"); g <- NULL }
  else { s <- try(summary(g)$coefficients, silent = TRUE)
    if (inherits(s, "try-error") || !TERM %in% rownames(s)) { out$GEE <- row1("GEE", NA, NA, "no_term"); g <- NULL }
    else out$GEE <- row1("GEE", s[TERM, 1], s[TERM, 2], "ok") }
  est_g <- if (!is.null(g)) unname(coef(g)[TERM]) else NA_real_
  for (nm in c("md", "kc", "fg")) {
    lab <- paste0("GEE-", toupper(nm))
    if (is.null(g)) { out[[lab]] <- row1(lab, NA, NA, "no_base_gee"); next }
    v <- hto(try(do.call(getExportedValue("geesmv", paste0("GEE.var.", nm)), list(f, id = "id", family = binomial, data = D, corstr = "exchangeable")), silent = TRUE), TO["smv"])
    if (inherits(v, "hard_timeout")) { out[[lab]] <- row1(lab, est_g, NA, "timeout"); next }
    if (inherits(v, "try-error")) { out[[lab]] <- row1(lab, est_g, NA, "error"); next }
    cn <- try(names(coef(g)), silent = TRUE); k <- match(TERM, cn); vv <- try(as.numeric(v$cov.beta)[k], silent = TRUE)
    if (inherits(vv, "try-error") || is.na(vv) || !is.finite(vv) || vv <= 0) out[[lab]] <- row1(lab, est_g, NA, "bad_var")
    else out[[lab]] <- row1(lab, est_g, sqrt(vv), "ok")
  }
  q <- hto(try(suppressWarnings(glmer(update(f, . ~ . + (1 | id)), family = binomial, data = dat, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))), silent = TRUE), TO["glmm"])
  if (inherits(q, "hard_timeout")) out$GLMM <- row1("GLMM", NA, NA, "timeout")
  else if (inherits(q, "try-error")) out$GLMM <- row1("GLMM", NA, NA, "error")
  else { s <- try(summary(q)$coefficients, silent = TRUE)
    if (inherits(s, "try-error") || !TERM %in% rownames(s)) out$GLMM <- row1("GLMM", NA, NA, "no_term")
    else out$GLMM <- row1("GLMM", s[TERM, 1], s[TERM, 2], "ok") }
  rbindlist(out)
}
FS <- list(pieno = outcome ~ treat + baseline + visit + age + sex + center, ridotto = outcome ~ treat + baseline)
REF <- rbindlist(lapply(names(FS), function(fn) fit_set(d, FS[[fn]])[, spec := fn]))
fwrite(REF, file.path(TABDIR, "V4_S7b_resp_ref.csv")); print(REF, digits = 6)
t0 <- Sys.time()
run_spec <- function(fn) {
  res <- mclapply(seq_along(draws), function(i) {
    z <- draws[[i]]; s <- d[d$pid %in% z$sel, ]
    r <- try(fit_set(s, FS[[fn]]), silent = TRUE)
    if (inherits(r, "try-error")) r <- rbindlist(lapply(c("MMM", "GEE", "GEE-MD", "GEE-KC", "GEE-FG", "GLMM"), function(mm) row1(mm, NA, NA, "outer_error")))
    r[, `:=`(spec = fn, N = z$N, rep = z$rep, nrows = nrow(s), nev = sum(s$outcome))]
    r
  }, mc.cores = NCORE, mc.preschedule = TRUE)
  bad <- sapply(res, function(x) !is.data.table(x)); cat(sprintf("spec=%s: failed children=%d\n", fn, sum(bad)))
  rbindlist(res[!bad])
}
ALL <- rbindlist(lapply(names(FS), run_spec))
cat(sprintf("fits done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
fwrite(ALL, file.path(TABDIR, "V4_S7b_resp_fits.csv"))
print(ALL[, .N, by = .(spec, method)][order(spec, method)])
print(dcast(ALL[, .(pct = round(100 * mean(status == "ok"), 2)), by = .(spec, method, N)], spec + method ~ N, value.var = "pct"))
