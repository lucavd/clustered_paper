# s10t_cov.R -- S10 "Testo": per design point x arm summaries of usable estimates and SEs (coverage, null excess, SE ratios)
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L); setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
CV <- qs_read(file.path(OUT, "cov_rows.qs2"))
S <- CV[, .(n = .N, mean_est = mean(estimate), sd_est = sd(estimate), mean_se = mean(se), rms_se = sqrt(mean(se^2)),
            m2 = mean(estimate^2)), by = c("strategy", DPK)]
# coverage needs the truth per replicate class: nominal for GenBinaryY, GLMM-AR1, Copula-AR1; GLMM-native from the truth file
tr <- fread("tables/V4_S2_truth_glmmnative.csv")[, .(dgp, scen = paste0("M", model), label, ntime, th_nat = truth_marginal_weighted)]
CV <- merge(CV, tr, by = c("dgp", "scen", "label", "ntime"), all.x = TRUE)
CV[, th := fifelse(label == "manual", 0, fifelse(scen == "M4", 0.05, 0.2))]
CV[dgp == "GLMM-native" & !is.na(th_nat), th := th_nat]
CV[, cov_n := abs(estimate - th) <= qnorm(0.975) * se]
C <- CV[, .(n = .N, k_cov = sum(cov_n), has_nat_truth = any(!is.na(th_nat))), by = c("strategy", DPK)]
fwrite(merge(S, C, by = c("strategy", DPK)), file.path(OUT, "cov_dp.csv"))
# GLMM Wald on the copula (model SE over the spread of the estimates)
R <- qs_read(file.path(OUT, "x_recs.qs2"))[dgp == "Copula-AR1"]; R[, est_design := generated & !design_degenerate]; R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
G <- qs_read(file.path(OUT, "x_cop_study.qs2"))[method == "glmer" & term %in% "time:trt"]
G <- merge(G, R[, c(K6, "est_design", "const"), with = FALSE], by = K6)
G[, use := est_design & !const & is.finite(estimate) & is.finite(se) & se > 0 & abs(estimate) < 10 & se < 10 & is.finite(p)]
fwrite(G[use == TRUE, .(n = .N, mean_est = mean(estimate), sd_est = sd(estimate), mean_se = mean(se), rms_se = sqrt(mean(se^2))), by = DPK], file.path(OUT, "cop_glmm_se.csv"))
