# =============================================================================
# s7g_tabkeys.R -- S7g: registry source for the new supplement tables S32-S35 and for constants of the new sentences.
# Inputs : tables/V5_S7g_calib.csv, V5_S7g_filter_t4.csv, V5_S7g_m67_cc.csv, V5_S7g_m67_cc_byused.csv, V5_S7g_lus_glmm.csv,
#          V5_S7g_lus_glmm_slope.csv (if present), V5_S7f_table4.csv.
# Output : tables/V5_S7g_tab_keys.csv (key, value, mcse, kind, source, note), read by the token rule f:gt|<key>.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
EX <- list(); put <- function(key, value, kind = "sig2", src = "", note = "") {
  x <- data.table(k_ = key, value = as.character(value), mcse = "", kind = kind, source = src, note = note); setnames(x, "k_", "key"); EX[[length(EX) + 1]] <<- x }
T4 <- fread("tables/V5_S7f_table4.csv")
put("n_datasets", T4[strategy == "gee_ar1_fg", n_dp] * 2L * 1000L, "int", "tables/V5_S7f_table4.csv", "design points x 2 hypotheses x 1,000 replicates")
put("lim_lo", 5L, "int", "R/s7g/s7g_filter.R"); put("lim_hi", 20L, "int", "R/s7g/s7g_filter.R")
put("brad_lo", 0.5, "plain", "Bradley (1978)", "liberal criterion, lower multiple of alpha"); put("brad_hi", 1.5, "plain", "Bradley (1978)")
put("fgwc_central", 90L, "int", "tables/V5_S7g_fg_wc.csv", "5th to 95th percentile")
CAL <- fread("tables/V5_S7g_calib.csv")
for (i in seq_len(nrow(CAL))) { r <- CAL[i]; s <- r$strategy
  put(paste0("cal|", s, "|hi"), r$n_hi, "int", "tables/V5_S7g_calib.csv"); put(paste0("cal|", s, "|ne"), r$n_notexcl, "int", "tables/V5_S7g_calib.csv")
  put(paste0("cal|", s, "|na"), r$n_notassessed, "int", "tables/V5_S7g_calib.csv"); put(paste0("cal|", s, "|max"), r$max_rate, "tab1", "tables/V5_S7g_calib.csv")
  put(paste0("cal|", s, "|maxn"), r$max_n, "int", "tables/V5_S7g_calib.csv") }
F4 <- fread("tables/V5_S7g_filter_t4.csv")
for (i in seq_len(nrow(F4))) { r <- F4[i]
  put(sprintf("flt|%s|%d|hi", r$strategy, as.integer(r$L)), r$hi_all, "int", "tables/V5_S7g_filter_t4.csv")
  put(sprintf("flt|%s|%d|pu", r$strategy, as.integer(r$L)), r$pu, "tab1", "tables/V5_S7g_filter_t4.csv") }
put("flt_same_choice", as.integer(F4[chosen == TRUE, uniqueN(strategy)] == 1L && F4[eligible == TRUE, uniqueN(strategy)] == 2L), "int", "tables/V5_S7g_filter_t4.csv", "1 = same chosen strategy and same two eligible at every limit")
SH <- c("GenBinaryY" = "GB", "GLMM-native" = "NAT")
M <- fread("tables/V5_S7g_m67_cc.csv"); B <- fread("tables/V5_S7g_m67_cc_byused.csv")
for (i in seq_len(nrow(M))) { r <- M[i]; b <- sprintf("m67|%s|%s|%s|", SH[[r$dgp]], r$scen, r$arm)
  put(paste0(b, "npts"), r$n_points, "int", "tables/V5_S7g_m67_cc.csv"); put(paste0(b, "rb"), r$mean_rb, "sig2", "tables/V5_S7g_m67_cc.csv")
  put(paste0(b, "rbmax"), r$max_abs_rb, "sig2", "tables/V5_S7g_m67_cc.csv"); put(paste0(b, "ser"), 100 * r$median_ser, "tab1", "tables/V5_S7g_m67_cc.csv")
  put(paste0(b, "cov"), r$mean_cov, "tab1", "tables/V5_S7g_m67_cc.csv"); put(paste0(b, "covmin"), r$min_cov, "tab1", "tables/V5_S7g_m67_cc.csv")
  put(paste0(b, "covhi"), r$n_cov_clearly_low, "int", "tables/V5_S7g_m67_cc.csv"); put(paste0(b, "shB"), r$mean_share_B, "tab1", "tables/V5_S7g_m67_cc.csv")
  for (us in c("A", "B")) { u <- B[dgp == r$dgp & scen == r$scen & arm == r$arm & used == us]
    put(paste0(b, "cov", us), u$mean_cov, "tab1", "tables/V5_S7g_m67_cc_byused.csv"); put(paste0(b, "ser", us), 100 * u$median_ser, "tab1", "tables/V5_S7g_m67_cc_byused.csv") } }
L1 <- fread("tables/V5_S7g_lus_glmm.csv")
put("lus_pt_min", L1[re == "id_pt", min(ratio_fpc)], "sig2", "tables/V5_S7g_lus_glmm.csv"); put("lus_pt_max", L1[re == "id_pt", max(ratio_fpc)], "sig2", "tables/V5_S7g_lus_glmm.csv")
put("lus_pt_var_med_max", L1[re == "id_pt" & N <= 25, max(median_var_pt)], "sig2", "tables/V5_S7g_lus_glmm.csv", "largest median patient-timepoint variance, N <= 25")
put("lus_pt_Nmax0", L1[re == "id_pt" & median_var_pt < 1e-6, max(N)], "int", "tables/V5_S7g_lus_glmm.csv", "largest N with median variance at the boundary")
if (file.exists("tables/V5_S7g_lus_glmm_slope.csv")) { L2 <- fread("tables/V5_S7g_lus_glmm_slope.csv")
  put("lus_sl_min", L2[re == "id_slope", min(ratio_fpc)], "sig2", "tables/V5_S7g_lus_glmm_slope.csv"); put("lus_sl_max", L2[re == "id_slope", max(ratio_fpc)], "sig2", "tables/V5_S7g_lus_glmm_slope.csv")
  put("lus_sl_Nmin", L2[re == "id_slope", min(N)], "int", "tables/V5_S7g_lus_glmm_slope.csv")
  put("lus_sl_var_full", L2[re == "id_slope", full_var_pt[1]], "sig2", "tables/V5_S7g_lus_glmm_slope.csv", "random-slope variance of posterior zones, full cohort") }
X <- rbindlist(EX); stopifnot(!anyDuplicated(X$key)); fwrite(X, "tables/V5_S7g_tab_keys.csv"); cat(nrow(X), "keys\n")
