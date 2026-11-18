# Mechanism illustration: lag correlation and prevalence by occasion, scenario M1, H0 ("manual"),
# large datasets per mechanism, generators of the study used unchanged (repo read-only).
.libPaths("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
source(path.expand("~/MMM/R/s3/s3_common.R"))   # load_env, set_rep_seed, RUNDIR_V4 (no writes at source time)
COP_DIR <- file.path(RUNDIR_V4, "copula")
load_env_cop <- function() {                      # file list of R/s7e/s7e_common.R::load_env_cop (generation part)
  env <- new.env(parent = globalenv())
  for (f in c("copula_params.R", "gen_copula_data.R", "generation_models_copula.R")) sys.source(file.path(COP_DIR, f), envir = env)
  for (f in c("manual/generate_data.R", "ar1/make_all_data.R")) sys.source(file.path(RUNDIR_V4, f), envir = env)
  env
}
N <- 2000L; TT <- 10L; NREP <- 10L; MODEL <- 1L; LABEL <- "manual"
MECH <- c("GenBinaryY", "GLMM-native", "GLMM-AR1", "Copula-AR1")
gen_one <- function(m, i) {
  env <- if (m == "Copula-AR1") load_env_cop() else load_env(m, RUNDIR_V4)
  set_rep_seed(MODEL, i)
  t0 <- proc.time()[[3]]
  d <- as.data.frame(env$make_all_data(N, TT, i, MODEL, LABEL))
  d <- d[order(d$id, d$time), c("id", "time", "trt", "binY")]
  par <- if (m == "Copula-AR1") c(rho = env$COPULA_RHO, floor = env$COPULA_FLOOR) else if (m == "GLMM-AR1") c(rho = env$AR1_RHO, sigma2 = env$AR1_SIGMA2) else c(none = NA)
  list(m = m, i = i, d = d, secs = proc.time()[[3]] - t0, par = par)
}
jobs <- expand.grid(i = seq_len(NREP), m = MECH, stringsAsFactors = FALSE)
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) tryCatch(gen_one(jobs$m[k], jobs$i[k]),
          error = function(e) list(m = jobs$m[k], i = jobs$i[k], err = conditionMessage(e))), mc.cores = 4L, mc.preschedule = FALSE)
for (x in res) if (!is.null(x$err)) cat("GEN ERROR", x$m, x$i, x$err, "\n")
res <- Filter(function(x) is.null(x$err), res)
cat("generated:", length(res), "datasets\n")
for (x in res) cat(sprintf("%s rep %d: %.1f s, n_subj %d\n", x$m, x$i, x$secs, length(unique(x$d$id))))
for (x in res[!duplicated(sapply(res, `[[`, "m"))]) cat("params", x$m, paste(names(x$par), x$par, collapse = " "), "\n")

pair_rows <- list(); prev_rows <- list(); pool_rows <- list(); info <- list()
for (x in res) {
  d <- x$d
  W <- matrix(NA_integer_, N, TT); W[cbind(d$id, d$time + 1L)] <- d$binY
  first <- !duplicated(d$id); arm <- integer(N); arm[d$id[first]] <- d$trt[first]
  info[[length(info) + 1]] <- data.frame(mechanism = x$m, rep = x$i, n_subj = length(unique(d$id)), n_obs = nrow(d), mean_y = mean(d$binY))
  for (j in 1:TT) for (a in 0:1) { y <- W[arm == a, j]; y <- y[!is.na(y)]
    prev_rows[[length(prev_rows) + 1]] <- data.frame(mechanism = x$m, rep = x$i, occasion = j - 1L, arm = a, n = length(y), y1 = sum(y)) }
  for (L in 1:(TT - 1)) {
    pa <- c(); pb <- c()
    for (j in 1:(TT - L)) { k <- j + L
      for (a in 0:1) { ok <- arm == a & !is.na(W[, j]) & !is.na(W[, k])
        ya <- W[ok, j]; yb <- W[ok, k]
        r <- if (length(ya) > 2 && sd(ya) > 0 && sd(yb) > 0) cor(ya, yb) else NA_real_
        pair_rows[[length(pair_rows) + 1]] <- data.frame(mechanism = x$m, rep = x$i, lag = L, occ_j = j - 1L, occ_k = k - 1L, arm = a, n = sum(ok), r = r) }
      ok <- !is.na(W[, j]) & !is.na(W[, k]); pa <- c(pa, W[ok, j]); pb <- c(pb, W[ok, k]) }
    pool_rows[[length(pool_rows) + 1]] <- data.frame(mechanism = x$m, rep = x$i, lag = L, n_pairs = length(pa), r_pooled = cor(pa, pb))
  }
}
pairs <- do.call(rbind, pair_rows); prev <- do.call(rbind, prev_rows); pool <- do.call(rbind, pool_rows); info <- do.call(rbind, info)
write.csv(pairs, "mech_pair_corr.csv", row.names = FALSE)
write.csv(prev, "mech_prev_counts.csv", row.names = FALSE)
write.csv(pool, "mech_pooled_corr.csv", row.names = FALSE)
write.csv(info, "mech_dataset_info.csv", row.names = FALSE)
cat("pairs rows", nrow(pairs), " NA r:", sum(is.na(pairs$r)), "\n")
print(aggregate(r ~ mechanism + lag, data = pairs, FUN = mean))
cat("MMLB version", as.character(packageVersion("MMLB")), " R", R.version.string, "\n")
