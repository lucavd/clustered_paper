# =============================================================================
# s5_common.R -- shared setup for S5 (inference across cells under CRN, D-23).
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
OUT_S5 <- file.path(REPO, "results_v4/S5"); dir.create(OUT_S5, recursive = TRUE, showWarnings = FALSE)
LOG_S5 <- file.path(REPO, "reports/morris/S5_logs"); dir.create(LOG_S5, recursive = TRUE, showWarnings = FALSE)
if (nzchar(VARIANT)) {                                    # D-32 variant tree (see s4_common.R)
  OUT_S5 <- file.path(REPO, "results_v4", paste0("variant_", VARIANT), "S5")
  LOG_S5 <- file.path(REPO, "reports/morris", paste0("variant_", VARIANT), "S5_logs")
  for (p in c(OUT_S5, LOG_S5)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
B_FULL <- 2000L
B_USE  <- as.integer(Sys.getenv("S5_B", B_FULL))          # S5_B < 2000 only for development runs
DEV    <- B_USE < B_FULL
WALD_ALL <- c(WALD, SENS)                                  # Wald routes incl. glmer_bobyqa
PRIM8 <- PRIM                                              # mmm gee gee_md gee_kc gee_fg glmer shrink glmer_lrt
WALD7 <- WALD
MARG5 <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc")
DGPS  <- c("GenBinaryY", "GLMM-native", "GLMM-AR1")
NP    <- c(4L, 5L, 4L, 4L, 5L)                             # fixed effects per scenario model (tref_new.R)
SETS  <- c("own", "intersection", "pair")
DPK   <- c("dgp", "scen", "model", "N", "ntime")
tag_out <- function(x) if (DEV) paste0("DEV_", x) else x
