# =============================================================================
# test_dgm_parameters.R -- S1 tests P1-P4 (pre-registered in
# docs/sessions/2026-10-15_S1.md, commit b0826c5).
#
#  P1  table values == values stated in the Morris audit (B1, B3, B4)
#  P2  H0/H1 parity read off the table: which generators differ in more than
#      the tested coefficient
#  P3  counter-proof: a regenerator that knows ONLY the table (plus the
#      engines) reproduces make_all_data() bit for bit
#  P4  the staged code of the M1-M2 chunks (job 23035d84) and .rundir/ give
#      identical datasets for M1-M2
#
# Usage: Rscript --vanilla R/testing/test_dgm_parameters.R
# =============================================================================

LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(MMLB); library(dplyr); library(MASS); library(fastGHQuad)})
REPO   <- path.expand("~/MMM")
RUNDIR <- file.path(REPO, ".rundir")
OLDDIR <- file.path(REPO, "_superseded/rundir_m1m2_job23035d84")
source(file.path(REPO, "R/protocol/dgm_common.R"))
tab <- read.csv(file.path(REPO, "tables/V4_dgm_parameters.csv"), stringsAsFactors = FALSE)
res <- list()
# Production RNG kind, read from the table (row 'RNGkind'), set before any draw.
rk <- unique(tab$value[tab$parameter == "RNGkind"]); stopifnot(length(rk) == 1L)
rk <- trimws(strsplit(rk, ",")[[1]])
RNGkind(kind = rk[1], normal.kind = rk[2], sample.kind = rk[3])
cat("RNGkind set from table:", paste(RNGkind(), collapse = ", "), "\n")
log_res <- function(test, item, pass, detail = "")
  res[[length(res) + 1L]] <<- data.frame(test = test, item = item, pass = pass,
                                         detail = detail, stringsAsFactors = FALSE)
val <- function(dgp, model, arm, p) {
  v <- tab$value[tab$dgp == dgp & tab$model == model & tab$arm == arm & tab$parameter == p]
  if (length(v) != 1L) NA_character_ else v
}

# ---- P1 --------------------------------------------------------------------------
combos <- unique(tab[, c("dgp", "model", "arm")])
log_res("P1", "n_combinations", nrow(combos) == 30L, sprintf("%d combinations", nrow(combos)))
beta_m <- list("1" = c(-2, 0.1, 0.2, 0.2), "2" = c(-2, 0.1, 0.2, 0.01, 0.2),
               "3" = c(-2, 0.1, 0.2, 0.2), "4" = c(-2, 0.2, 0.2, 0.05),
               "5" = c(-2, 0.1, 0.2, 0.05, 0.2))
sig_H1 <- list(1, 1, c(1, 0.01), c(1, 0.01), c(1, 0.01, 0.01))
gam_H1 <- list(1, 1, 1, c(1, 0.01), c(1, 0.1, 0.05))
sig_H0 <- list(1, 1, c(1, 0), c(1, 0), c(1, 0, 0))
gam_H0 <- list(1, 1, 1, c(1, 0), c(1, 0, 0))
D_of <- function(m, dgp) {
  if (m <= 2) return(matrix(1))
  if (m <= 4) return(matrix(c(1, 0.01, 0.01, 0.001), 2, byrow = TRUE))
  matrix(c(1, 0.01, 0.01, 0.01, 0.001, 0, 0.01, 0, if (dgp == "GLMM-AR1") 0.01 else 0.001),
         3, byrow = TRUE)
}
eq <- function(a, b) length(a) == length(b) && all(a == b)
for (m in 1:5) for (arm in c("H1", "H0")) {
  bm <- beta_m[[m]]; if (arm == "H0") bm[length(bm)] <- 0
  g <- "GenBinaryY"
  log_res("P1", sprintf("%s M%d %s beta", g, m, arm), eq(txt_num(val(g, m, arm, "beta")), bm),
          val(g, m, arm, "beta"))
  s_exp <- if (arm == "H1") sig_H1[[m]] else sig_H0[[m]]
  g_exp <- if (arm == "H1") gam_H1[[m]] else gam_H0[[m]]
  log_res("P1", sprintf("%s M%d %s sigma", g, m, arm), eq(txt_num(val(g, m, arm, "sigma")), s_exp),
          val(g, m, arm, "sigma"))
  log_res("P1", sprintf("%s M%d %s gamma", g, m, arm), eq(txt_num(val(g, m, arm, "gamma")), g_exp),
          val(g, m, arm, "gamma"))
  log_res("P1", sprintf("%s M%d %s q", g, m, arm), eq(txt_num(val(g, m, arm, "q")), 20),
          val(g, m, arm, "q"))
  g <- "GLMM-native"
  bc <- txt_num(val(g, m, arm, "beta_c"))
  ex <- 1.16 * bm
  log_res("P1", sprintf("%s M%d %s beta_c == 1.16 x beta_m", g, m, arm),
          length(bc) == length(ex) && all(abs(bc - ex) < 1e-12),
          sprintf("beta_c = %s; 1.16 x beta_m = %s", val(g, m, arm, "beta_c"), num_txt(ex)))
  log_res("P1", sprintf("%s M%d %s beta_c == round(1.16 x beta_m, 3)", g, m, arm),
          eq(bc, round(ex, 3)), "")
  for (g in c("GLMM-native", "GLMM-AR1")) {
    Dm <- txt_mat(val(g, m, arm, "D"))
    log_res("P1", sprintf("%s M%d %s D", g, m, arm),
            all(dim(Dm) == dim(D_of(m, g))) && all(Dm == D_of(m, g)), val(g, m, arm, "D"))
  }
  g <- "GLMM-AR1"
  log_res("P1", sprintf("%s M%d %s target_beta_m", g, m, arm),
          eq(txt_num(val(g, m, arm, "target_beta_m")), bm), val(g, m, arm, "target_beta_m"))
  log_res("P1", sprintf("%s M%d %s rho, sigma2_ar", g, m, arm),
          eq(txt_num(val(g, m, arm, "rho")), 0.6) && eq(txt_num(val(g, m, arm, "sigma2_ar")), 0.5),
          paste(val(g, m, arm, "rho"), val(g, m, arm, "sigma2_ar")))
}

# ---- P2 --------------------------------------------------------------------------
for (g in DGPS) for (m in 1:5) {
  a <- tab[tab$dgp == g & tab$model == m & tab$arm == "H1" & tab$block != "derived", ]
  b <- tab[tab$dgp == g & tab$model == m & tab$arm == "H0" & tab$block != "derived", ]
  stopifnot(identical(a$parameter, b$parameter))
  diffs <- character(0)
  for (k in seq_len(nrow(a))) {
    if (identical(a$value[k], b$value[k]) && identical(a$terms[k], b$terms[k])) next
    if (a$parameter[k] %in% c("beta", "beta_c", "target_beta_m")) {
      terms <- trimws(strsplit(a$terms[k], ",")[[1]])
      d <- terms[txt_num(a$value[k]) != txt_num(b$value[k])]
      if (identical(d, "time:trt")) next
      diffs <- c(diffs, sprintf("%s[%s]", a$parameter[k], paste(d, collapse = "+")))
    } else diffs <- c(diffs, a$parameter[k])
  }
  log_res("P2", sprintf("%s M%d parity", g, m), length(diffs) == 0L,
          if (length(diffs)) paste("also differ:", paste(diffs, collapse = ", ")) else "tested coefficient only")
}

# ---- P3: regenerate from the table only ----------------------------------------
regen_from_table <- function(rows, env, N, ntime, i) {
  get1 <- function(p) rows$value[rows$parameter == p]
  dgp <- rows$dgp[1]; model <- rows$model[1]
  ev <- new.env(parent = globalenv()); ev$N <- N; ev$ntime <- ntime
  set.seed(seed_of(model, i))
  nclust <- eval(str2lang(get1("followup_rule")), ev)
  n_obs  <- sum(nclust)
  dat <- data.frame(id = rep(seq_len(N), nclust),
                    time = unlist(lapply(nclust, function(x) seq_len(x) - 1)))
  dr <- rows[grepl("^draw_", rows$parameter), ]
  dr <- dr[order(dr$parameter), ]
  for (k in seq_len(nrow(dr))) {
    if (dr$terms[k] == "nclust") next
    lvl <- sub(".*level=([a-z]+);.*", "\\1", dr$note[k])
    std <- grepl("standardised=TRUE", dr$note[k], fixed = TRUE)
    e <- str2lang(dr$value[k])
    if (lvl == "observation") { if (as.character(e[[1]]) == "sample") e$size <- n_obs else e$n <- n_obs }
    x <- eval(e, ev)
    if (lvl == "subject") x <- rep(x, nclust)
    if (std) x <- as.numeric(scale(x))
    dat[[dr$terms[k]]] <- x
  }
  cols <- trimws(strsplit(get1("data_columns"), ",")[[1]])
  dat <- dat[, intersect(cols, names(dat)), drop = FALSE]
  fml <- function(p) as.formula(get1(p), env = globalenv())
  if (dgp == "GenBinaryY") {
    # GenBinaryY resolves `id` through substitute(), so it must be passed as `id = id`
    id <- dat$id
    MMLB::GenBinaryY(mean.formula = fml("mean.formula"), lv.formula = fml("lv.formula"),
                     t.formula = fml("t.formula"), beta = txt_num(get1("beta")),
                     sigma = txt_num(get1("sigma")), gamma = txt_num(get1("gamma")),
                     id = id, data = dat, q = txt_num(get1("q")), Yname = "binY")
  } else if (dgp == "GLMM-native") {
    env$gen_glmm_data(mean.formula = fml("mean.formula"), re.formula = fml("re.formula"),
                      beta_c = txt_num(get1("beta_c")), D = txt_mat(get1("D")),
                      id = "id", data = dat, Yname = "binY")
  } else {
    D <- txt_mat(get1("D")); s2 <- txt_num(get1("sigma2_ar"))
    bc <- env$calibrate_beta_c_ar1(mean.formula = fml("mean.formula"),
                                   re.formula = fml("re.formula"),
                                   target_beta_m = txt_num(get1("target_beta_m")),
                                   D = D, sigma2_ar = s2, data = dat)
    env$gen_ar1_data(mean.formula = fml("mean.formula"), re.formula = fml("re.formula"),
                     beta_c = as.numeric(bc), D = D, rho = txt_num(get1("rho")),
                     sigma2_ar = s2, time = "time", id = "id", data = dat, Yname = "binY")
  }
}
same_data <- function(a, b) {
  if (inherits(a, "error") || inherits(b, "error"))
    return(inherits(a, "error") && inherits(b, "error") &&
           identical(conditionMessage(a), conditionMessage(b)))
  ca <- setdiff(names(a), "sim_id"); cb <- setdiff(names(b), "sim_id")
  identical(ca, cb) && nrow(a) == nrow(b) &&
    all(vapply(ca, function(k) identical(as.numeric(a[[k]]), as.numeric(b[[k]])), logical(1)))
}
safe <- function(expr) tryCatch(expr, error = function(e) e)
GRID <- expand.grid(N = c(10, 45), ntime = c(3, 6), i = c(1, 500))
n_err <- 0L
envs <- setNames(lapply(DGPS, load_family, rundir = RUNDIR), DGPS)
for (g in DGPS) for (m in 1:5) for (arm in c("H1", "H0")) {
  rows <- tab[tab$dgp == g & tab$model == m & tab$arm == arm, ]
  for (r in seq_len(nrow(GRID))) {
    N <- GRID$N[r]; Tt <- GRID$ntime[r]; i <- GRID$i[r]
    set.seed(seed_of(m, i))
    orig <- safe(envs[[g]]$make_all_data(N, Tt, i, m, ARMS[[arm]]))
    reg  <- safe(regen_from_table(rows, envs[[g]], N, Tt, i))
    if (inherits(orig, "error")) n_err <- n_err + 1L
    log_res("P3", sprintf("%s M%d %s N=%d T=%d i=%d", g, m, arm, N, Tt, i),
            same_data(orig, reg),
            paste(if (inherits(orig, "error")) paste("generation error:", conditionMessage(orig)) else "",
                  if (inherits(reg, "error")) paste("regenerator error:", conditionMessage(reg)) else ""))
  }
}

# ---- P4: M1-M2 staged code vs .rundir ---------------------------------------------
old_envs <- setNames(lapply(DGPS, load_family, rundir = OLDDIR), DGPS)
for (g in DGPS) for (m in 1:2) for (arm in c("H1", "H0")) for (r in seq_len(nrow(GRID))) {
  N <- GRID$N[r]; Tt <- GRID$ntime[r]; i <- GRID$i[r]
  set.seed(seed_of(m, i)); a <- safe(old_envs[[g]]$make_all_data(N, Tt, i, m, ARMS[[arm]]))
  set.seed(seed_of(m, i)); b <- safe(envs[[g]]$make_all_data(N, Tt, i, m, ARMS[[arm]]))
  log_res("P4", sprintf("%s M%d %s N=%d T=%d i=%d", g, m, arm, N, Tt, i), same_data(a, b))
}

out <- do.call(rbind, res)
write.csv(out, file.path(REPO, "tables/V4_S1_tests_P1_P4.csv"), row.names = FALSE)
summ <- aggregate(pass ~ test, out, function(x) sprintf("%d/%d", sum(x), length(x)))
print(summ, row.names = FALSE)
cat("P3 datasets whose generation raised an error (in both paths):", n_err, "\n")
fails <- out[!out$pass, ]
if (nrow(fails)) { cat("FAILURES:\n"); print(fails, row.names = FALSE) }
