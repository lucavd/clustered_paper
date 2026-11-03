# =============================================================================
# dgm_common.R -- helpers shared by the S1 parameter extraction and its tests.
#
# Loads each data-generating family from a staged copy of the code (.rundir/
# by default, i.e. the code that produced results_full/ for M3-M5) into its own
# environment, exactly as run_full.R does inside a worker, and provides
# recorders that capture (i) the arguments handed to the generating engines and
# (ii) every random-number call made while the design is built.
# =============================================================================

DGPS <- c("GenBinaryY", "GLMM-native", "GLMM-AR1")
ARMS <- c(H1 = "mmlb", H0 = "manual")
RNG_FNS <- c("rnorm", "rbinom", "runif", "sample")
rng_real <- function(fn) get(fn, envir = asNamespace(if (fn == "sample") "base" else "stats"))

FAMILY_FILES <- list(
  "GenBinaryY"  = c(sprintf("generation_model%d.R", 1:5),
                    "generate_data.R", "make_all_data.R"),
  "GLMM-native" = c("gen_glmm_data.R", sprintf("manual/generation_model%d.R", 1:5),
                    "manual/generate_data.R", "manual/make_all_data.R"),
  "GLMM-AR1"    = c("gen_ar1_data.R", "ar1/generation_models_ar1.R",
                    "manual/generate_data.R", "ar1/make_all_data.R")
)

ENGINE_OF <- c("GenBinaryY" = "GenBinaryY", "GLMM-native" = "gen_glmm_data",
               "GLMM-AR1" = "gen_ar1_data")

# Seed rule of R/run_full.R (checked against the file by the extractor).
seed_of <- function(model, i) 20260805L + 1000L * as.integer(model) + as.integer(i)

load_family <- function(dgp, rundir) {
  env <- new.env(parent = globalenv())
  for (f in FAMILY_FILES[[dgp]]) sys.source(file.path(rundir, f), envir = env)
  env
}

# Install recorders in a family environment. Engine wrappers record their
# arguments and then call the real engine; RNG wrappers record design-phase
# draws only (draws made inside an engine are the engine's business).
install_recorders <- function(env) {
  rec <- new.env()
  rec$engine_calls <- list(); rec$rng <- list(); rec$depth <- 0L
  wrap_engine <- function(name, real) {
    force(name); force(real)
    function(...) {
      rec$engine_calls[[length(rec$engine_calls) + 1L]] <- list(engine = name, args = list(...))
      rec$depth <- rec$depth + 1L
      on.exit(rec$depth <- rec$depth - 1L)
      real(...)
    }
  }
  wrap_rng <- function(fn) {
    force(fn)
    real <- rng_real(fn)
    function(...) {
      if (rec$depth == 0L) {
        mc <- match.call(real, as.call(c(as.name(fn), list(...))))
        rec$rng[[length(rec$rng) + 1L]] <- list(fn = fn, args = as.list(mc)[-1])
      }
      real(...)
    }
  }
  env$GenBinaryY <- wrap_engine("GenBinaryY", MMLB::GenBinaryY)
  for (nm in c("gen_glmm_data", "gen_ar1_data", "calibrate_beta_c_ar1"))
    if (exists(nm, envir = env, inherits = FALSE))
      assign(nm, wrap_engine(nm, get(nm, envir = env)), envir = env)
  for (fn in RNG_FNS) assign(fn, wrap_rng(fn), envir = env)
  rec
}

# Static walk of a function body: every RNG call in source order, with the
# variable it is assigned to, plus every assignment (target -> RHS text) and the
# follow-up truncation threshold (the literal compared with runif()).
static_design <- function(fun) {
  draws <- list(); assigns <- list(); thr <- NULL
  walk <- function(e, target) {
    if (!is.call(e)) return(invisible())
    fn <- paste(deparse(e[[1]]), collapse = "")
    if (fn %in% c("<-", "=") && length(e) == 3L) {
      target <- paste(deparse(e[[2]]), collapse = "")
      assigns[[length(assigns) + 1L]] <<- list(target = target,
                                              rhs = paste(deparse(e[[3]]), collapse = " "))
      walk(e[[3]], target); return(invisible())
    }
    if (fn == "<" && is.call(e[[2]]) && identical(e[[2]][[1]], as.name("runif")))
      thr <<- eval(e[[3]], baseenv())
    if (fn %in% RNG_FNS) {
      real <- rng_real(fn)
      mc <- match.call(real, e)
      draws[[length(draws) + 1L]] <<- list(fn = fn, target = target,
        call = paste(deparse(mc, width.cutoff = 500L), collapse = ""))
    }
    for (k in seq_along(e)[-1]) if (!identical(e[[k]], quote(expr = ))) walk(e[[k]], target)
  }
  walk(body(fun), NA_character_)
  list(draws = draws, assigns = assigns, threshold = thr)
}

# Text helpers for the CSV: numbers written with 15 significant digits, which
# reproduce every literal in the generators exactly; matrices row-major with ';'.
num_txt <- function(x) paste(sprintf("%.15g", as.numeric(x)), collapse = ", ")
mat_txt <- function(m) {
  m <- as.matrix(m)
  paste(apply(m, 1, num_txt), collapse = "; ")
}
txt_num <- function(s) as.numeric(strsplit(gsub(";", ",", s), ",")[[1]])
txt_mat <- function(s) {
  rows <- strsplit(s, ";")[[1]]
  do.call(rbind, lapply(rows, function(r) as.numeric(strsplit(r, ",")[[1]])))
}
f_txt <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = "")
