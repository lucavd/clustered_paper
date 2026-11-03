# =============================================================================
# s9_keys.R -- S9: registry source for numbers of the v6 text that are not a single cell of an existing CSV
# (token family f:n9|<key>, read by py/s8b_build_spec.py). Every value is computed here from a file.
# Output: tables/V6_S9_keys.csv (key, value, mcse, kind, source, how).
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
OUT <- list(); put <- function(key, value, kind, source, how) { stopifnot(length(value) == 1, !is.na(value)); OUT[[length(OUT) + 1]] <<- setnames(data.table(k_ = key, value = as.character(value), mcse = "", kind = kind, source = source, how = how), "k_", "key") }
## Table 2(a): design points per mechanism and scenario group (one strategy, one row per point)
RD <- unique(fread("tables/V5_S7f_rule_dp.csv")[, .(dgp, scen, N, ntime)])
stopifnot(nrow(RD) == 2376L)
RD[, sg := fifelse(scen %in% sprintf("M%d", 1:5), "M1-M5", fifelse(scen == "M1c", "M1c", "M6-M7"))]
C <- RD[, .N, by = .(dgp, sg)]
for (i in seq_len(nrow(C))) put(sprintf("ndp|%s|%s", C$dgp[i], C$sg[i]), C$N[i], "int", "tables/V5_S7f_rule_dp.csv", "design points of the mechanism in the scenario group")
put("grid_N_min", RD[sg == "M1-M5", min(N)], "int", "tables/V5_S7f_rule_dp.csv", "smallest N of the main grid")
put("grid_N_max", RD[sg == "M1-M5", max(N)], "int", "tables/V5_S7f_rule_dp.csv", "largest N of the main grid")
put("grid_T_min", RD[sg == "M1-M5", min(ntime)], "int", "tables/V5_S7f_rule_dp.csv", "smallest T of the main grid")
put("grid_T_max", RD[sg == "M1-M5", max(ntime)], "int", "tables/V5_S7f_rule_dp.csv", "largest T of the main grid")
put("m1c_N_max", RD[sg == "M1c", max(N)], "int", "tables/V5_S7f_rule_dp.csv", "largest N of M1c")
put("m1c_T_max", RD[sg == "M1c", max(ntime)], "int", "tables/V5_S7f_rule_dp.csv", "largest T of M1c")
## constant-outcome datasets (D-82)
K <- fread("tables/V6_S9_const_keys_summary.csv")
put("const_n", K[, sum(constant)], "int", "tables/V6_S9_const_keys_summary.csv", "non-degenerate datasets with a constant outcome")
put("const_all0_pct", 100 * K[, sum(all_zero)] / K[, sum(constant)], "sig2", "tables/V6_S9_const_keys_summary.csv", "of which all outcomes 0 (%)")
put("const_small_pct", 100 * K[small == TRUE, sum(constant)] / K[, sum(constant)], "sig2", "tables/V6_S9_const_keys_summary.csv", "of which with T 2-3 and N <= 25 (%)")
put("datasets_n", K[, sum(datasets)], "int", "tables/V6_S9_const_keys_summary.csv", "datasets of the 2,376 design points")
## the design points at which GEE AR(1) + FG + t is clearly above 7.5% (Table 5), with the marginalized model at the same points
A <- fread("tables/V5_S7f_cp_S7d_dp_rates.csv")[ref == "t"]; B <- fread("tables/V5_S7f_cp_S7e_dp_rates.csv")[ref == "t"]
AB <- rbind(A, B, fill = TRUE); HI <- AB[strategy == "gee_ar1_fg" & hi == TRUE][order(dgp)]
stopifnot(nrow(HI) == fread("tables/V5_S7f_table4.csv")[strategy == "gee_ar1_fg", hi_all])
PRE <- fread("results_v4/S9/pre/tables/V5_S7f_cp_S7d_dp_rates.csv")[ref == "t" & strategy == "gee_ar1_fg"]
put("fgpt|n", nrow(HI), "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", "points with GEE AR(1) FG t clearly above 7.5%")
put("fgpt|same_cell", as.integer(uniqueN(HI[, .(scen, N, ntime)]) == 1L), "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", "1 = all at the same scenario, N and T")
for (j in seq_len(nrow(HI))) { h <- HI[j]; tag <- sprintf("fgpt|%d", j)
  put(paste0(tag, "|N"), h$N, "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("N of point", j, h$dgp, h$scen))
  put(paste0(tag, "|T"), h$ntime, "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("T of point", j))
  put(paste0(tag, "|n"), h$n, "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("usable replicates under H0 at point", j))
  put(paste0(tag, "|k"), h$kk, "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("rejections at point", j))
  put(paste0(tag, "|rate"), h$rate, "sig2", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("rejection rate (%) at point", j))
  p <- PRE[dgp == h$dgp & scen == h$scen & N == h$N & ntime == h$ntime]; stopifnot(nrow(p) == 1)
  put(paste0(tag, "|n_pre"), p$n, "int", "results_v4/S9/pre/tables/V5_S7f_cp_S7d_dp_rates.csv", paste("usable replicates before the constant-outcome rule at point", j))
  put(paste0(tag, "|k_pre"), p$kk, "int", "results_v4/S9/pre/tables/V5_S7f_cp_S7d_dp_rates.csv", paste("rejections before the rule at point", j))
  m <- AB[strategy == "study_mmm" & dgp == h$dgp & scen == h$scen & N == h$N & ntime == h$ntime]; stopifnot(nrow(m) == 1)
  put(paste0(tag, "|mmm_n"), m$n, "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("usable replicates of the marginalized model at point", j))
  put(paste0(tag, "|mmm_k"), m$kk, "int", "tables/V5_S7f_cp_S7d_dp_rates.csv", paste("rejections of the marginalized model at point", j)) }
## fewclust cited as implementation (PI, 2026-11-02): check-3 counts of fewclust 0.3.1 (R/s9/s9_fewclust_coherence.R)
FC <- fread("tables/V6_S9_fewclust_coherence.csv")
put("fc|n_cmp", FC[item == "n_cmp", value], "int", "tables/V6_S9_fewclust_coherence.csv", "datasets regenerated from the study and compared with fewclust 0.3.1")
put("fc|n_coh", FC[item == "n_coh", value], "int", "tables/V6_S9_fewclust_coherence.csv", "of which the default coherence condition changed the reported fit")
X <- rbindlist(OUT); fwrite(X, "tables/V6_S9_keys.csv"); print(X)
