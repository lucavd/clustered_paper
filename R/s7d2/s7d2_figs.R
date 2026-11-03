# s7d2_figs.R -- figures of the S7d2 report from tables/V4_S7d2_dp.csv (reports/morris/S7d2_figures/).
# fig1: coverage at H1 of the test-inversion set and of GEE AR(1) + Fay-Graubard + t, per design point, by N and T;
# fig2: width ratio (geometric mean per point) interval of the test / GEE AR(1) FG t, H1, by N and T;
# fig3: small designs (T 2-3, N <= 25), H1: share of datasets with an answer: test at 0 as in S7d, test at 0 with the
#       D-57 rules, finite interval of the test, usable GEE AR(1) FG t.
source(path.expand("~/MMM/R/s3/s3_common.R")); suppressMessages({ library(data.table); library(ggplot2) })
D <- fread(Sys.getenv("S7D2_DP", file.path(REPO, "tables/V4_S7d2_dp.csv"))); OUT <- Sys.getenv("S7D2_FIGOUT", file.path(REPO, "reports/morris/S7d2_figures")); dir.create(OUT, showWarnings = FALSE)
D[, mech := factor(fifelse(dgp == "Copula-AR1", "Copula (serial only)", dgp), levels = c("GenBinaryY", "GLMM-native", "GLMM-AR1", "Copula (serial only)"))]
D[, T := factor(paste0("T = ", ntime), levels = paste0("T = ", c(2, 3, 4, 6, 8, 10)))]
th <- theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
H1 <- D[label == "mmlb"]
L <- rbind(H1[, .(mech, N, T, what = "Sign-flip test inversion", cov = cov_t, n = n_ut)], H1[, .(mech, N, T, what = "GEE AR(1) + Fay-Graubard + t", cov = cov_g, n = n_ug)])[n >= 20]
p1 <- ggplot(L, aes(N, cov, colour = mech)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 92.5, ymax = 97.5, alpha = 0.15) +
  geom_point(size = 0.7, alpha = 0.7) + facet_grid(what ~ T) + scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "Number of subjects (N)", y = "Coverage under H1 (%)") + coord_cartesian(ylim = c(85, 100)) + th
ggsave(file.path(OUT, "fig1_coverage.png"), p1, width = 9, height = 5, dpi = 150)
p2 <- ggplot(H1[is.finite(gm_ratio) & n_pr >= 20], aes(N, gm_ratio, colour = mech)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.9, ymax = 1.1, alpha = 0.15) +
  geom_hline(yintercept = 1, linewidth = 0.3) + geom_point(size = 0.7, alpha = 0.7) + facet_wrap(~ T, nrow = 1) + scale_y_log10() +
  scale_colour_brewer(palette = "Dark2", name = NULL) + labs(x = "Number of subjects (N)", y = "Width ratio test / GEE (log scale)") + th
ggsave(file.path(OUT, "fig2_width_ratio.png"), p2, width = 9, height = 3.6, dpi = 150)
S <- H1[small == TRUE, .(`Test at 0, as in the previous run` = mean(pu_flip_s7d), `Test at 0, score not zero (D-57)` = mean(pu_flip_guard),
                         `Finite interval of the test` = mean(pfin), `GEE AR(1) + Fay-Graubard + t usable` = mean(pu_g)), by = .(N, ntime)]
S <- melt(S, id.vars = c("N", "ntime"), variable.name = "what", value.name = "share")
S[, design := factor(sprintf("N = %d, T = %d", N, ntime), levels = unique(sprintf("N = %d, T = %d", N[order(ntime, N)], ntime[order(ntime, N)])))]
p3 <- ggplot(S, aes(design, share, fill = what)) + geom_col(position = position_dodge(0.8), width = 0.75) + scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Datasets with an answer (%), H1") + guides(fill = guide_legend(nrow = 2)) + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "fig3_small_designs.png"), p3, width = 9, height = 4.2, dpi = 150)
cat("figures written to", OUT, "\n")
