# s7c_timing.R [workers] -- serial per-dataset timing of the S7d candidate methods (organisation only; no measure enters the paper).
# 20 replicates (i = 1..20, H0) per design point: the 18 pilot points plus N = 145, T = 10 in the three mechanisms.
source(path.expand("~/MMM/R/s7c/s7c_common.R"))
.libPaths(c(.libPaths(), S7C_LIB))
suppressMessages({ library(clubSandwich); library(fewclust) })
a <- commandArgs(TRUE); W <- if (length(a)) as.integer(a[1]) else 10L
pts <- rbind(CJ(dgp = DGPS, N = PILOT_N, ntime = PILOT_T), CJ(dgp = DGPS, N = 145L, ntime = 10L))
tk <- pts[, .(i = 1:20), by = .(dgp, N, ntime)]
CORS <- c(exch = "exchangeable", ind = "independence", ar1 = "ar1")
SMV  <- c(exch = "exchangeable", ind = "independence", ar1 = "AR-M")
tim <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err"))
  list(t = proc.time()[["elapsed"]] - t0, ok = !(inherits(v, "err") || inherits(v, "hard_timeout") || is.null(v)),
       msg = if (inherits(v, "err")) substr(as.character(v), 1, 120) else if (inherits(v, "hard_timeout")) "timeout" else NA_character_) }
one <- function(k) {
  g <- tk[k]; env <- load_env(g$dgp, RUNDIR_V4)
  set_rep_seed(1L, g$i); d <- as.data.frame(gen_data(env, "M1", 1L, g$N, g$ntime, g$i, "manual")); d <- d[order(d$id, d$time), ]
  f <- binY ~ time * trt; rows <- list()
  add <- function(op, r) rows[[length(rows) + 1]] <<- data.table(op = op, sec = r$t, ok = r$ok, msg = r$msg)
  fits <- list()
  for (cn in names(CORS)) {
    r <- tim(fits[[cn]] <- geepack::geeglm(f, id = id, data = d, family = binomial, corstr = CORS[[cn]], scale.fix = TRUE)); add(paste0("geeglm_", cn), r)
    for (ty in c("md", "kc", "fg")) {
      add(paste0("geesmv_", ty, "_", cn), tim(env$.with_hard_timeout(do.call(getExportedValue("geesmv", paste0("GEE.var.", ty)),
            list(f, id = "id", family = binomial, data = d, corstr = SMV[[cn]])), secs = 5)))
      if (!is.null(fits[[cn]])) add(paste0("fewclust_", ty, "_", cn), tim(fewclust::gee_vcov(fits[[cn]], type = ty)))
    }
    if (!is.null(fits[[cn]])) add(paste0("cr2_", cn), tim(clubSandwich::vcovCR(fits[[cn]], cluster = d$id, type = "CR2")))
  }
  set.seed(RESAMPLE_SEED(g$i)); rt <- tryCatch(run_tests(d, timing = TRUE), error = function(e) NULL)
  if (!is.null(rt)) { tt <- attr(rt, "timing"); for (nm in names(tt)) add(paste0("test_", nm), list(t = tt[[nm]], ok = TRUE, msg = NA_character_)) }
  out <- rbindlist(rows); out[, `:=`(dgp = g$dgp, N = g$N, ntime = g$ntime, i = g$i, n_obs = nrow(d))]; out
}
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(tk)), function(k) tryCatch(one(k), error = function(e) data.table(op = "ERROR", msg = conditionMessage(e))),
                          mc.cores = W, mc.preschedule = FALSE)
R <- rbindlist(res, fill = TRUE); fwrite(R, file.path(OUT_S7C, "timing.csv"))
cat("timing rows", nrow(R), "errors", sum(R$op == "ERROR"), "wall", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
print(R[op != "ERROR", .(sec = round(mean(sec), 4), ok = round(mean(ok), 3)), by = .(op, N)][order(op, N)], nrows = 200)
