# s7e_summary.R -- S7e: summary tables and figures for the report, from the tables written by s7e_agg.R (no new computation of
# performance). Outputs: tables/V4_S7e_summary.csv (per strategy: copula vs each study mechanism, main grid, t at N - p),
# tables/V4_S7e_bias_summary.csv, reports/morris/S7e_figures/*.png
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(data.table); library(ggplot2) })
REPO <- path.expand("~/MMM"); TAB <- file.path(REPO, "tables"); FIG <- file.path(REPO, "reports/morris/S7e_figures"); dir.create(FIG, showWarnings = FALSE)
NEW <- fread(file.path(TAB, "V4_S7e_dp_rates.csv"))[ref == "t"]
RS <- fread(file.path(TAB, "V4_S7d_dp_rates.csv"))[ref == "t" & scen %in% sprintf("M%d", 1:5)]; RS[strategy == "study_mmm", strategy := "mmm_cc"]
RB <- fread(file.path(TAB, "V4_S7b_rule_dp.csv"))[ref == "t" & scen %in% sprintf("M%d", 1:5) & strategy %in% c("glmer|t", "shrink|t")]; RB[, strategy := sub("\\|t$", "", strategy)]
KEY <- c("flip", "gee_ar1_fg", "mmm_cc", "gee_exch_fg", "gee_exch_md", "gee_exch_kc", "gee_exch_unc", "gee_ar1_unc", "gee_ind_fg", "glmer", "shrink", "glmer_lrt", "mmm_study")
LAB <- c(flip = "Sign-flip score test", gee_ar1_fg = "GEE AR(1) + Fay-Graubard + t", mmm_cc = "MMM with convergence check + t", gee_exch_fg = "GEE exch. + Fay-Graubard + t",
         gee_exch_md = "GEE exch. + Mancl-DeRouen + t", gee_exch_kc = "GEE exch. + Kauermann-Carroll + t", gee_exch_unc = "GEE exch., uncorrected sandwich + t",
         gee_ar1_unc = "GEE AR(1), uncorrected sandwich + t", gee_ind_fg = "GEE indep. + Fay-Graubard + t", glmer = "GLMM Wald + t", shrink = "Shrinkage GLM + t",
         glmer_lrt = "GLMM likelihood-ratio test", mmm_study = "MMM as in the study (no check) + t")
cols <- c("dgp", "scen", "N", "ntime", "strategy", "rate", "hi", "lo", "pu", "upow")
ALL <- rbind(NEW[strategy %in% KEY, ..cols], RS[strategy %in% KEY, ..cols], RB[, ..cols])
SM <- ALL[, .(points = .N, mean_type1 = mean(rate, na.rm = TRUE), clearly_above = sum(hi), clearly_below = sum(lo), mean_usable = mean(pu), mean_power = mean(upow)), by = .(strategy, dgp)]
SM[, label := LAB[strategy]]; setorder(SM, strategy, dgp)
fwrite(SM, file.path(TAB, "V4_S7e_summary.csv"))
BC <- fread(file.path(TAB, "V4_S7e_bias_cov.csv"))
BS <- BC[label == "mmlb", .(points = .N, mean_bias_rel = mean(bias_rel, na.rm = TRUE), max_abs_bias_rel = max(abs(bias_rel), na.rm = TRUE),
         mean_relerr_modse = 100 * mean(modse / empse - 1, na.rm = TRUE), mean_cov_n = mean(cov_n, na.rm = TRUE), mean_cov_t = mean(cov_t, na.rm = TRUE),
         cov_n_below_92.5 = sum(cov_n + 2 * mcse_cov_n < 92.5, na.rm = TRUE)), by = strategy]
BS[, label := LAB[strategy]]; fwrite(BS, file.path(TAB, "V4_S7e_bias_summary.csv"))
# ---- figures ------------------------------------------------------------------------------------------------------------
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
PR <- fread(file.path(TAB, "V4_S7e_dependence_profiles.csv"))[mech %in% c("GenBinaryY", "GLMM-native", "GLMM-AR1", "copula floor=0.00 rho=0.6")]
PR[mech == "copula floor=0.00 rho=0.6", mech := "Gaussian copula (serial, rho = 0.6)"]
g1 <- ggplot(PR, aes(lag, cor, colour = mech, shape = mech)) + geom_hline(yintercept = 0, colour = "grey60") + geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = 1:9) + labs(x = "Distance between visits (lag)", y = "Within-subject Pearson correlation", colour = NULL, shape = NULL) + th
ggsave(file.path(FIG, "fig1_dependence_profiles.png"), g1, width = 6.5, height = 4, dpi = 200)
F2 <- NEW[strategy %in% c("gee_ar1_fg", "mmm_cc", "glmer", "flip")]; F2[, lab := factor(LAB[strategy], levels = LAB[c("gee_ar1_fg", "mmm_cc", "glmer", "flip")])]
F2[, Tlab := factor(paste("T =", ntime), levels = paste("T =", c(2, 3, 4, 6, 8, 10)))]
g2 <- ggplot(F2[n >= 20], aes(N, rate, colour = scen)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 7.5, alpha = 0.12) +
  geom_hline(yintercept = 5, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.4) + geom_point(size = 0.8) +
  facet_grid(lab ~ Tlab) + labs(x = "Number of subjects (N)", y = "Type I error (%), nominal 5%", colour = "Scenario") + th + theme(strip.text.y = element_text(size = 7))
ggsave(file.path(FIG, "fig2_type1_copula.png"), g2, width = 9, height = 7.5, dpi = 200)
OLD <- rbind(RS[strategy %in% c("gee_ar1_fg", "mmm_cc", "flip"), ..cols], RB[strategy == "glmer", ..cols])
OM <- OLD[, .(old = mean(rate, na.rm = TRUE)), by = .(scen, N, ntime, strategy)]
F3 <- merge(F2[, .(scen, N, ntime, strategy, lab, Tlab, rate)], OM, by = c("scen", "N", "ntime", "strategy")); F3[, d := rate - old]
g3 <- ggplot(F3, aes(factor(ntime), d)) + geom_hline(yintercept = 0, colour = "grey40") + geom_hline(yintercept = c(-1, 1), linetype = 3, colour = "grey40") +
  geom_boxplot(outlier.size = 0.6, fill = "grey90") + facet_wrap(~ lab, nrow = 1) +
  labs(x = "Maximum number of visits (T)", y = "Type I error: copula minus mean of the study mechanisms (points)") + th
ggsave(file.path(FIG, "fig3_type1_difference.png"), g3, width = 9, height = 3.8, dpi = 200)
print(SM[strategy %in% c("flip", "gee_ar1_fg", "mmm_cc", "glmer")], digits = 3); print(BS, digits = 3)
