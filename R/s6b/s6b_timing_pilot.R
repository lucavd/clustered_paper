# s6b_timing_pilot.R -- S6b: wall-clock timing only (no performance measures are computed or printed),
# to size the new scenario with a large latent-SD slope. Seeds 77000000 + k, outside the study scheme.
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressMessages({ library(MMLB); library(lme4); library(geepack); library(geesmv); library(dplyr); library(tibble); library(MASS); library(fastGHQuad); library(data.table) })
data.table::setDTthreads(1L)
RD <- path.expand("~/MMM/R/v4/rundir")
for (f in c("analysis_model3.R", "analysis_smallsample.R", "do_analysis.R", "gen_glmm_data.R")) sys.source(file.path(RD, f), envir = globalenv())
gen_cov <- function(N, ntime) {
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))
  id <- rep(seq(N), nclust); trt <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1)); age <- scale(rep(rnorm(N, 50, 10), nclust))
  data.frame(id, time, trt, age) |> arrange(id, time)
}
gen_pop <- function(N, ntime, label, s1) {
  b <- if (label == "mmlb") c(-2, 0.1, 0.2, 0.2) else c(-2, 0.1, 0.2, 0)
  GenBinaryY(mean.formula = ~ time * trt, lv.formula = ~time, t.formula = ~1, beta = b, sigma = c(1, s1), gamma = 1,
             id = id, data = gen_cov(N, ntime), q = 20, Yname = "binY")
}
gen_glmm <- function(N, ntime, label, v1) {
  b <- if (label == "mmlb") c(-2.32, 0.116, 0.232, 0.232) else c(-2.32, 0.116, 0.232, 0)
  dd <- gen_cov(N, ntime)
  gen_glmm_data(mean.formula = ~ time * trt, re.formula = ~1 + time, beta_c = b, D = matrix(c(1, 0.01, 0.01, v1), 2),
                id = "id", data = dd, Yname = "binY")
}
G <- CJ(gen = c("pop_0.1", "pop_0.25", "glmm_0.04"), N = c(55L, 145L), ntime = c(6L, 10L), label = c("mmlb", "manual"))
REPS <- 6L
one <- function(g) {
  r <- G[g]; out <- vector("list", REPS)
  for (j in seq_len(REPS)) {
    RNGkind("L'Ecuyer-CMRG"); set.seed(77000000L + 100L * g + j)
    par <- as.numeric(sub(".*_", "", r$gen))
    t0 <- proc.time()[[3]]
    d <- as.data.frame(if (startsWith(r$gen, "pop")) gen_pop(r$N, r$ntime, r$label, par) else gen_glmm(r$N, r$ntime, r$label, par))
    tg <- proc.time()[[3]] - t0
    tm <- sapply(c("mmm", "gee", "glmer", "gee_ss", "glmer_lrt"), function(tc) { t0 <- proc.time()[[3]]
      try(do_analysis(d, r$N, r$ntime, r$label, 3, tc, j), silent = TRUE); proc.time()[[3]] - t0 })
    t0 <- proc.time()[[3]]; try(mm(binY ~ time * trt, lv.formula = ~1, t.formula = ~1, data = d, id = id, q = 20), silent = TRUE)
    tB <- proc.time()[[3]] - t0
    out[[j]] <- data.table(gen = r$gen, N = r$N, ntime = r$ntime, label = r$label, rep = j, t_gen = tg, t_B = tB, t(tm))
  }
  rbindlist(out)
}
t0 <- Sys.time()
L <- parallel::mclapply(seq_len(nrow(G)), function(g) tryCatch(one(g), error = function(e) { message("FAIL ", g, ": ", conditionMessage(e)); NULL }), mc.cores = 24, mc.preschedule = FALSE)
R <- rbindlist(Filter(Negate(is.null), L))
fwrite(R, "s6b_timing_pilot.csv")
S <- R[, lapply(.SD, mean), by = .(gen, N, ntime), .SDcols = c("t_gen", "mmm", "gee", "glmer", "gee_ss", "glmer_lrt", "t_B")]
S[, all_core_s := t_gen + mmm + gee + glmer + gee_ss + glmer_lrt + t_B]
print(S, digits = 3)
cat("wall", format(difftime(Sys.time(), t0, units = "mins")), "\n")
