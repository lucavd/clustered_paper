# =============================================================================
# s7b_common.R -- shared setup for S7b (reanalysis without new simulations).
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
OUT_S7B <- file.path(REPO, "results_v4/S7b"); dir.create(OUT_S7B, recursive = TRUE, showWarnings = FALSE)
LOG_S7B <- file.path(REPO, "reports/morris/S7b_logs"); dir.create(LOG_S7B, recursive = TRUE, showWarnings = FALSE)
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)   # mean-model coefficients (R/s5/s5_common.R, R/s6b/s6b_ns_agg_lib.R)
WALD7 <- c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink")
PRIM8 <- c(WALD7, "glmer_lrt")
flags <- function(H) {   # Bradley flags of R/s6/s6_rule_safety.R
  H[, rate := fifelse(n > 0, 100 * kk / n, NA_real_)]
  H[, mc := fifelse(n > 0, 100 * sqrt((rate / 100) * (1 - rate / 100) / n), NA_real_)]
  H[, `:=`(assess = n >= 20, pu = 100 * n / nd)]
  H[, `:=`(hi = assess & rate - 2 * mc > 7.5, lo = assess & rate + 2 * mc < 2.5)]
  H
}
