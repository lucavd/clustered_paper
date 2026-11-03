
# s10_probe.R -- S10 probe: regeneration timing, reproduction of n_obs/sum_y, separation check availability
.libPaths(c("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
S10LIB <- path.expand("~/R/s10_lib"); dir.create(S10LIB, showWarnings = FALSE)
if (!requireNamespace("detectseparation", lib.loc = c(S10LIB, .libPaths()), quietly = TRUE)) {
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
          HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)),
          install.packages.compile.from.source = "never")
  install.packages("detectseparation", lib = S10LIB, quiet = TRUE)
}
source(path.expand("~/MMM/R/s7d/s7d_common.R"))
.libPaths(c(.libPaths(), S10LIB))
suppressMessages(library(detectseparation)); cat("detectseparation", as.character(packageVersion("detectseparation")), "\n")
TK <- s7d_tasks()
pick <- TK[(dgp == "GenBinaryY" & scen == "M2" & N == 10 & ntime == 2 & label == "manual" & chunk == 1) |
           (dgp == "GLMM-native" & scen == "M1c" & N == 45 & ntime == 6 & label == "manual" & chunk == 1) |
           (dgp == "GLMM-AR1" & scen == "M5" & N == 25 & ntime == 3 & label == "mmlb" & chunk == 1)]
for (r in seq_len(nrow(pick))) { t <- pick[r]; if (!file.exists(t$file)) next
  rec <- readRDS(t$file)$recs; env <- load_env(t$dgp, RUNDIR_V4); gen <- s7d_gen(env, t)
  t0 <- proc.time()[[3]]; ok <- 0L; tsep <- 0; nsep <- 0L; nint <- 0L
  for (i in 1:100) { deg <- design_degenerate(t$model, t$N, t$ntime, i); set_rep_seed(t$model, i)
    d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) NULL); if (is.null(d)) next
    ok <- ok + (nrow(d) == rec$n_obs[i] && sum(d$binY) == rec$sum_y[i])
    t1 <- proc.time()[[3]]
    s <- tryCatch(glm(FORM_OF(t$scen), data = d, family = binomial, method = "detect_separation"), error = function(e) NULL)
    tsep <- tsep + proc.time()[[3]] - t1
    if (!is.null(s)) { nsep <- nsep + isTRUE(s$outcome); nint <- nint + isTRUE(is.infinite(s$coefficients[["time:trt"]])) } }
  cat(sprintf("%s %s N%d T%d %s: reproduced %d/100, total %.1fs (separation %.1fs), separated %d, time:trt infinite %d\n",
              t$dgp, t$scen, t$N, t$ntime, t$label, ok, proc.time()[[3]] - t0, tsep, nsep, nint)) }
