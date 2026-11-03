
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
S <- readRDS("results_v4/variant_kcfc/S4/rows_v4.rds")[term %in% "time:trt" & label == "manual" & scen %in% sprintf("M%d", 1:5)]
TC <- fread("tables/variant_kcfc/V4_S6_table_cells.csv")[table_id == "T3" & variant == "main"]
TC[, c("q", "g", "m") := tstrsplit(name, "|", fixed = TRUE)]
gm <- c(GB = "GenBinaryY", AR1 = "GLMM-AR1", NAT = "GLMM-native")
OUT <- list()
for (mm in unique(TC$m)) for (gg in names(gm)) {
  X <- S[method == mm & dgp == gm[[gg]]]; if (!nrow(X)) next
  P <- X[, .(n = sum(U0), k = sum(U0 & p < 0.05)), by = N][, d := abs(100 * k / n - 5)]
  OUT[[length(OUT) + 1]] <- data.table(m = mm, g = gg, rate = X[, 100 * sum(U0 & p < 0.05) / sum(U0)], wc = P[, max(d)], mad = P[, mean(d)])
}
R <- melt(rbindlist(OUT), id.vars = c("m", "g"), variable.name = "q", value.name = "recomp")
J <- merge(TC[, .(m, g, q, est)], R, by = c("m", "g", "q"), all.x = TRUE)
J[, diff := recomp - est]
print(J[, .(cells = .N, recomputed = sum(!is.na(recomp)), max_abs_diff = max(abs(diff), na.rm = TRUE)), by = m])
fwrite(J, "tables/V5_S7f_table2_recompute.csv")
