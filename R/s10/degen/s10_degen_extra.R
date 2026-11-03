# s10_degen_extra.R -- S10 "Degeneri": numbers quoted in the text of Sections 3.4, 4.1 and the abstract recomputed per version
# from results_v4/S10/degen/table5_points_versions.csv.gz; incremental reading of PD6 (version d against version a);
# combined gate table tables/V7_S10_degen_gates.csv (GD2, GD3, baseline reproduction).
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); suppressMessages(library(data.table)); setwd("/home/user/MMM")
P <- fread("results_v4/S10/degen/table5_points_versions.csv.gz")
fl <- P[scen == "M2" & N == 10 & ntime == 2 & dgp %in% c("GLMM-AR1", "GLMM-native") & strategy %in% c("gee_ar1_fg", "mmm_cc"), .(version, dgp, strategy, txt = sprintf("%d of %d", k, n))]
X <- rbind(
  fl[, .(version, quantity = sprintf("Section 3.4: H0 rejections of %s at M2 N=10 T=2 %s", strategy, dgp), value = txt)],
  P[strategy == "glmer", .(quantity = "Section 3.4: GLMM Wald clearly above 7.5%, of which GLMM-native", value = sprintf("%d, %d", sum(hi), sum(hi & dgp == "GLMM-native"))), by = version],
  P[strategy %in% c("gee_exch_fg", "gee_ind_fg"), .(quantity = sprintf("Section 3.4: %s clearly above 7.5%%, of which main grid (all N=10?)", strategy[1]), value = sprintf("%d, %d main, N=10 among them %d", sum(hi), sum(hi & grp == "main"), sum(hi & N == 10))), by = .(version, strategy)][, strategy := NULL],
  P[strategy %in% c("gee_ar1_md", "gee_exch_md", "gee_ind_md", "gee_ar1_kc", "gee_exch_kc", "gee_ind_kc", "glmer_lrt"), .(quantity = sprintf("Section 3.4: %s clearly above 7.5%%", strategy[1]), value = as.character(sum(hi))), by = .(version, strategy)][, strategy := NULL],
  P[strategy %in% c("gee_ar1_fg", "mmm_cc"), .(quantity = sprintf("Section 3.4: %s clearly below 2.5%% (all; N=10)", strategy[1]), value = sprintf("%d; %d", sum(lo), sum(lo & N == 10))), by = .(version, strategy)][, strategy := NULL],
  P[strategy %in% c("gee_ar1_fg", "mmm_cc", "flip"), .(quantity = sprintf("Abstract/3.4: %s usable %% (all; T2-3 & N<=25), integer as printed", strategy[1]), value = sprintf("%.0f%%; %.0f%%", mean(pu), mean(pu[small == TRUE]))), by = .(version, strategy)][, strategy := NULL],
  P[strategy %in% c("gee_ar1_fg", "mmm_cc", "flip"), .(quantity = sprintf("Section 3.4/4.1: %s power in T2-3 & N<=25 (%%)", strategy[1]), value = sprintf("%.2f", mean(upow[small == TRUE]))), by = .(version, strategy)][, strategy := NULL])
setorder(X, quantity, version); fwrite(X, "tables/V7_S10_degen_paper_numbers.csv"); print(X[version %in% c("base", "a", "b")], nrow = 200)
S <- fread("tables/V7_S10_degen_sensitivity.csv")[table == "Table 3"]
W <- dcast(S[version %in% c("a", "d")], row + column ~ version, value.var = "new"); W[, dd := d - a]; W[, what := sub(" .*$", "", column)]
cat(sprintf("PD6 incremental (d - a): rate/MAD %.3f to %.3f; worst case %.3f to %.3f\n", W[what %in% c("rate", "mad"), min(dd)], W[what %in% c("rate", "mad"), max(dd)], W[what == "wc", min(dd)], W[what == "wc", max(dd)]))
g2 <- fread("tables/V7_S10_degen_gates_gd2.csv"); g3 <- fread("tables/V7_S10_degen_gd3.csv"); gb <- fread("tables/V7_S10_degen_baseline.csv")
GA <- rbind(g2[, .(gate, check, pass, detail)], g3[, .(gate = "GD3", check, pass, detail = as.character(value))], gb[, .(gate = "baseline (step 2)", check, pass, detail)], fill = TRUE)
fwrite(GA, "tables/V7_S10_degen_gates.csv"); print(GA)

PD <- fread("tables/V7_S10_degen_predictions.csv")[!grepl("incremental", code)]
PD[code == "PD2", explanation := paste0(explanation, "; second part of the expectation not met: 53 Si usable results (27 and 26) carry all 23 H0 rejections (12 and 11) at the two points, against at most 10 expected")]
PD[code == "PD6", explanation := paste0(explanation, "; judged against the paper (baseline); the choice already changes under version a")]
PD <- rbind(PD, data.table(code = "PD6 (incremental reading)", quantity = "version d against version a: rule choice; Table 3 rate and MAD change",
  expected = "choice unchanged; Table 3 within +/-0.5", relevance = "categorical for the rule; -1..+1 for Table 3",
  observed = sprintf("chosen under a and d: gee_ar1_fg; rate/MAD change d - a %.3f to %.3f; worst case %.3f to %.3f", W[what %in% c("rate", "mad"), min(dd)], W[what %in% c("rate", "mad"), max(dd)], W[what == "wc", min(dd)], W[what == "wc", max(dd)]),
  verdict = "confermata", explanation = "alternative reading of 'rimozione in piu'': reported, not used for the verdict of PD6"))
fwrite(PD, "tables/V7_S10_degen_predictions.csv")
