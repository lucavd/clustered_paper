# =============================================================================
# s7g_agq_keys.R -- S7g: registry source for the adaptive-quadrature numbers of the manuscript and supplement (decision D-68).
# Inputs : tables/V5_S7g_agq_diff.csv, V5_S7g_agq_dp.csv, V5_S7g_agq_mech.csv, V5_S7g_agq_gates.csv (R/s7g/s7g_agg.R).
# Output : tables/V5_S7g_agq_keys.csv (key, value, mcse, kind, source, note), read by the token rule f:gq|<key>.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
DF <- fread("tables/V5_S7g_agq_diff.csv"); DP <- fread("tables/V5_S7g_agq_dp.csv"); ME <- fread("tables/V5_S7g_agq_mech.csv")
EX <- list(); put <- function(key, value, kind = "sig2", note = "", src = "tables/V5_S7g_agq_diff.csv", mcse = NA_real_) {
  x <- data.table(k_ = key, value = as.character(value), mcse = if (is.na(mcse)) "" else sprintf("%.17g", mcse), kind = kind, source = src, note = note)
  setnames(x, "k_", "key"); EX[[length(EX) + 1]] <<- x }
put("npts", nrow(DF), "int", "design points"); put("nodes", 10L, "int", "nAGQ main", "R/s7g/s7g_common.R"); put("nodes_check", 25L, "int", "nAGQ check", "R/s7g/s7g_common.R")
put("check_reps", 100L, "int", "replicates with the 25-node check", "R/s7g/s7g_common.R")
W <- dcast(fread("tables/V5_S7g_agq_dp.csv")[nAGQ %in% c(10L, 25L)], dgp + scen + N + ntime ~ nAGQ, value.var = "n")
G <- fread("tables/V5_S7g_agq_gates.csv")
agree <- as.numeric(sub(".*same t decision in ([0-9.]+)%.*", "\\1", G[gate == "GG4", detail]))
put("check_agree", agree, "sig2", "same t decision, 25 vs 10 nodes (%)", "tables/V5_S7g_agq_gates.csv")
ncmp <- as.integer(sub("^([0-9]+) fits usable for both.*", "\\1", G[gate == "GG4", detail]))
put("check_n", ncmp, "int", "fits usable with 10 and 25 nodes", "tables/V5_S7g_agq_gates.csv")
put("check_disagree", as.integer(round(ncmp * (1 - agree / 100))), "int", "fits with a different t decision, 25 vs 10 nodes", "tables/V5_S7g_agq_gates.csv")
put("gg2_rows", as.integer(sub("^([0-9]+) replicate rows.*", "\\1", G[gate == "GG2", detail])), "int", "replicate rows compared with the stored Laplace fits", "tables/V5_S7g_agq_gates.csv")
put("diff_min", min(abs(ME$mean_diff_t)), "sig2", "smallest |mean difference| over mechanisms (points)", "tables/V5_S7g_agq_mech.csv")
put("diff_max", max(abs(ME$mean_diff_t)), "sig2", "largest |mean difference| over mechanisms (points)", "tables/V5_S7g_agq_mech.csv")
x <- DF[N == 145L & ntime == 2L & scen %in% c("M1", "M2")]
put("t2n145_lap_min", min(x$rate_t_lap)); put("t2n145_lap_max", max(x$rate_t_lap)); put("t2n145_agq_min", min(x$rate_t_agq)); put("t2n145_agq_max", max(x$rate_t_agq))
put("hi_lap", DF[, sum(hi_t_lap)], "int"); put("hi_agq", DF[, sum(hi_t_agq)], "int")
put("hi_agq_cop", DF[dgp == "Copula-AR1", sum(hi_t_agq)], "int"); put("hi_agq_gb", DF[dgp == "GenBinaryY", sum(hi_t_agq)], "int")
put("hi_agq_nat", DF[dgp == "GLMM-native", sum(hi_t_agq)], "int"); put("hi_lap_nat", DF[dgp == "GLMM-native", sum(hi_t_lap)], "int")
put("hi_agq_T", paste(sort(unique(DF[hi_t_agq == TRUE, ntime])), collapse = ","), "text", "occasions of the points clearly above with quadrature")
hc <- DF[hi_t_agq == TRUE & dgp == "Copula-AR1"]; hg <- DF[hi_t_agq == TRUE & dgp == "GenBinaryY"]
put("hi_cop_min", min(hc$rate_t_agq)); put("hi_cop_max", max(hc$rate_t_agq)); put("hi_gb_min", min(hg$rate_t_agq)); put("hi_gb_max", max(hg$rate_t_agq))
put("hi_cop_Nmin", min(hc$N), "int")
lo <- DP[, .(lo = sum(lo_t)), by = nAGQ]
put("lo_lap", lo[nAGQ == 1L, lo], "int", "points clearly below 2.5%, Laplace", "tables/V5_S7g_agq_dp.csv"); put("lo_agq", lo[nAGQ == 10L, lo], "int", "points clearly below 2.5%, quadrature", "tables/V5_S7g_agq_dp.csv")
put("lo_agq_Nmax", DP[nAGQ == 10L & lo_t == TRUE, max(N)], "int", "largest N among points clearly below 2.5% with quadrature", "tables/V5_S7g_agq_dp.csv")
put("lo_agq_nat", DP[nAGQ == 10L & dgp == "GLMM-native", sum(lo_t)], "int", "", "tables/V5_S7g_agq_dp.csv")
x <- DF[dgp == "Copula-AR1" & ntime == 10L & N >= 25L]
put("cop_t10_lap_min", min(x$rate_t_lap)); put("cop_t10_lap_max", max(x$rate_t_lap)); put("cop_t10_agq_min", min(x$rate_t_agq)); put("cop_t10_agq_max", max(x$rate_t_agq))
put("cop_t10_Nmin", 25L, "int", "")
## Table S33: mechanism x occasions
SH <- c("GenBinaryY" = "GB", "GLMM-native" = "NAT", "GLMM-AR1" = "AR1", "Copula-AR1" = "COP")
TB <- DF[, .(n = .N, rlap = mean(rate_t_lap), ragq = mean(rate_t_agq), hlap = sum(hi_t_lap), hagq = sum(hi_t_agq)), by = .(dgp, ntime)]
TB <- merge(TB, DP[nAGQ == 10L, .(lagq = sum(lo_t), pagq = mean(pu)), by = .(dgp, ntime)], by = c("dgp", "ntime"))
TB <- merge(TB, DP[nAGQ == 1L, .(llap = sum(lo_t), plap = mean(pu)), by = .(dgp, ntime)], by = c("dgp", "ntime"))
fwrite(TB, "tables/V5_S7g_agq_tab.csv")
for (i in seq_len(nrow(TB))) { r <- TB[i]; b <- sprintf("t_%s_%d_", SH[[r$dgp]], r$ntime)
  put(paste0(b, "n"), r$n, "int", "", "tables/V5_S7g_agq_tab.csv"); put(paste0(b, "rlap"), r$rlap, "tab1", "", "tables/V5_S7g_agq_tab.csv")
  put(paste0(b, "ragq"), r$ragq, "tab1", "", "tables/V5_S7g_agq_tab.csv"); put(paste0(b, "hlap"), r$hlap, "int", "", "tables/V5_S7g_agq_tab.csv")
  put(paste0(b, "hagq"), r$hagq, "int", "", "tables/V5_S7g_agq_tab.csv"); put(paste0(b, "llap"), r$llap, "int", "", "tables/V5_S7g_agq_tab.csv")
  put(paste0(b, "lagq"), r$lagq, "int", "", "tables/V5_S7g_agq_tab.csv"); put(paste0(b, "plap"), r$plap, "tab1", "", "tables/V5_S7g_agq_tab.csv")
  put(paste0(b, "pagq"), r$pagq, "tab1", "", "tables/V5_S7g_agq_tab.csv") }
for (g in names(SH)) { m <- ME[dgp == g]; b <- paste0("m_", SH[[g]], "_")
  put(paste0(b, "diff"), m$mean_diff_t, "tab1", "mean AGQ - Laplace (points)", "tables/V5_S7g_agq_mech.csv", m$mcse_t)
  put(paste0(b, "lo95"), m$lo95_t, "tab1", "", "tables/V5_S7g_agq_mech.csv"); put(paste0(b, "hi95"), m$hi95_t, "tab1", "", "tables/V5_S7g_agq_mech.csv")
  put(paste0(b, "npts"), m$n_points, "int", "", "tables/V5_S7g_agq_mech.csv") }
X <- rbindlist(EX); stopifnot(!anyDuplicated(X$key)); fwrite(X, "tables/V5_S7g_agq_keys.csv"); print(X[, .(key, value, kind)], nrows = 300)
