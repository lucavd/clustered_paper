# =============================================================================
# test_chunk_reproduction.R -- S1 test P5 (pre-registered, commit b0826c5).
#
# Refits replicate i = 1 of cell N = 25, T = 4 for all 15 generators x 2 arms
# with the code in .rundir/, for the first three techniques of run_full.R
# (mmm, gee, glmer -> also yields the derived 'shrink' rows), in the same order
# and after the same set.seed() as the production worker, and compares estimate
# and SE with the stored chunk in results_full/ (read only).
# Tolerance: max absolute difference <= 1e-8; NaN/NA must match NaN/NA.
# =============================================================================

LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressMessages({
  library(MMLB); library(lme4); library(geepack); library(geesmv)
  library(dplyr); library(tibble); library(MASS); library(fastGHQuad)
})
REPO   <- path.expand("~/MMM")
RUNDIR <- file.path(REPO, ".rundir")
source(file.path(REPO, "R/protocol/dgm_common.R"))
ANALYSIS <- c(sprintf("analysis_model%d.R", 1:5), "analysis_smallsample.R", "do_analysis.R")
N <- 25; Tt <- 4; i <- 1; TECHS <- c("mmm", "gee", "glmer"); TOL <- 1e-8
tab <- read.csv(file.path(REPO, "tables/V4_dgm_parameters.csv"), stringsAsFactors = FALSE)
rk <- unique(tab$value[tab$parameter == "RNGkind"]); stopifnot(length(rk) == 1L)
rk <- trimws(strsplit(rk, ",")[[1]])
if (nzchar(Sys.getenv("MMM_RNGKIND_OVERRIDE"))) rk <- trimws(strsplit(Sys.getenv("MMM_RNGKIND_OVERRIDE"), ",")[[1]])
RNGkind(kind = rk[1], normal.kind = rk[2], sample.kind = rk[3])
cat("RNGkind:", paste(RNGkind(), collapse = ", "), "\n")

out <- list()
for (g in DGPS) {
  env <- load_family(g, RUNDIR)
  for (f in ANALYSIS) sys.source(file.path(RUNDIR, f), envir = env)
  for (m in 1:5) for (arm in names(ARMS)) {
    lab  <- ARMS[[arm]]
    file <- file.path(REPO, "results_full",
                      sprintf("%s_m%d_N%d_T%d_%s_c01.rds", gsub("[^A-Za-z0-9]", "", g), m, N, Tt, lab))
    st <- readRDS(file)
    st <- as.data.frame(st[st$sim_id == i & st$method %in% c(TECHS, "shrink"), ])
    set.seed(seed_of(m, i))
    d <- env$make_all_data(N, Tt, i, m, lab)
    rf <- lapply(TECHS, function(tc)
      tryCatch(env$do_analysis(d, N, Tt, lab, m, tc, i), error = function(e) NULL))
    rf <- as.data.frame(dplyr::bind_rows(Filter(Negate(is.null), rf)))
    key <- function(x) paste(x$method, x$term, sep = "|")
    ks <- key(st); kr <- key(rf)
    common <- intersect(ks, kr)
    cmp <- function(col) {
      a <- st[[col]][match(common, ks)]; b <- rf[[col]][match(common, kr)]
      na_ok <- identical(is.na(a), is.na(b))
      dmax <- if (any(!is.na(a) & !is.na(b))) max(abs(a - b), na.rm = TRUE) else 0
      c(dmax = dmax, na_ok = na_ok)
    }
    ce <- cmp("estimate"); cs <- cmp("se")
    st_meth <- sort(unique(st$method)); rf_meth <- sort(unique(rf$method))
    status_same <- identical(st$status[match(common, ks)], rf$status[match(common, kr)])
    out[[length(out) + 1L]] <- data.frame(
      dgp = g, model = m, arm = arm, rows_stored = nrow(st), rows_refit = nrow(rf),
      rows_common = length(common), keys_identical = setequal(ks, kr),
      methods_stored = paste(st_meth, collapse = ","), methods_refit = paste(rf_meth, collapse = ","),
      max_abs_diff_estimate = ce[["dmax"]], max_abs_diff_se = cs[["dmax"]],
      na_pattern_same = as.logical(ce[["na_ok"]]) && as.logical(cs[["na_ok"]]),
      status_same = status_same,
      pass = setequal(ks, kr) && ce[["dmax"]] <= TOL && cs[["dmax"]] <= TOL &&
             as.logical(ce[["na_ok"]]) && as.logical(cs[["na_ok"]]),
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, out)
OUTF <- if (nzchar(Sys.getenv("MMM_RNGKIND_OVERRIDE"))) "tables/V4_S1_tests_P5_mersenne.csv" else "tables/V4_S1_tests_P5.csv"
write.csv(res, file.path(REPO, OUTF), row.names = FALSE)
cat(sprintf("P5: %d/%d generator x arm pass | max |diff| estimate = %.3g, se = %.3g | rows compared = %d\n",
            sum(res$pass), nrow(res), max(res$max_abs_diff_estimate), max(res$max_abs_diff_se),
            sum(res$rows_common)))
cat(sprintf("status strings identical: %d/%d\n", sum(res$status_same), nrow(res)))
if (any(!res$pass)) print(res[!res$pass, ], row.names = FALSE)
