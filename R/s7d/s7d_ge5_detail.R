# s7d_ge5_detail.R -- explains the GE5 mismatches (independence, geesmv vs fewclust, MD and FG): lists the datasets
# with relative difference >= 1e-8 and whether the geeglm fit is usable by the study filter (|est| < 10, 0 < se < 10).
source(path.expand("~/MMM/R/s4/s4_common.R")); suppressMessages(library(data.table)); setDTthreads(8L)
f <- list.files(file.path(REPO, "results_v4/S7d/rows"), pattern = "_c01\\.rds$", full.names = TRUE)   # replicates 1-100 hold i <= 50
R <- rbindlist(parallel::mclapply(f, function(x) readRDS(x)$rows[cor == "ind" & sim_id <= 50], mc.cores = 16L), fill = TRUE)
rel <- function(x, y) abs(x - y) / pmax(abs(y), 1e-12)
R[, `:=`(r_md = rel(smv_md, se_md), r_fg = rel(smv_fg, se_fg))]
R[, usable_md := is.finite(est) & abs(est) < 10 & is.finite(se_md) & se_md > 0 & se_md < 10]
R[, usable_fg := is.finite(est) & abs(est) < 10 & is.finite(se_fg) & se_fg > 0 & se_fg < 10]
M <- R[(is.finite(r_md) & r_md >= 1e-8) | (is.finite(r_fg) & r_fg >= 1e-8)]
cat("mismatch datasets:", nrow(M), "\n")
print(M[, .(dgp, scen, N, ntime, label, sim_id, est = signif(est, 4), se_unc = signif(se_unc, 4), se_md = signif(se_md, 4), smv_md = signif(smv_md, 4),
            se_fg = signif(se_fg, 4), smv_fg = signif(smv_fg, 4), r_md = signif(r_md, 3), r_fg = signif(r_fg, 3), usable_md, usable_fg)][order(-abs(est))], nrows = 200)
cat("mismatches among usable (study filter): md", M[usable_md & is.finite(r_md) & r_md >= 1e-8, .N], " fg", M[usable_fg & is.finite(r_fg) & r_fg >= 1e-8, .N], "\n")
b <- R[usable_md & is.finite(r_md)]; cat("usable md both:", nrow(b), "max rel", format(max(b$r_md), digits = 3), "\n")
b <- R[usable_fg & is.finite(r_fg)]; cat("usable fg both:", nrow(b), "max rel", format(max(b$r_fg), digits = 3), "\n")
O <- R[, .(only_geesmv_md = sum(is.finite(smv_md) & !is.finite(se_md)), only_fc_md = sum(!is.finite(smv_md) & is.finite(se_md)),
           only_geesmv_md_usable = sum(is.finite(smv_md) & !is.finite(se_md) & is.finite(est) & abs(est) < 10 & smv_md < 10),
           only_fc_md_usable = sum(!is.finite(smv_md) & usable_md))]
print(O)
print(R[!is.finite(smv_md) & usable_md, .N, by = .(why = substr(smv_msg_md, 1, 60))][order(-N)][1:6])
print(R[is.finite(smv_md) & !is.finite(se_md), .N, by = .(why = substr(msg_md, 1, 60))][order(-N)][1:6])
fwrite(M, file.path(REPO, "tables/V4_S7d_GE5_mismatches.csv")); fwrite(O, file.path(REPO, "tables/V4_S7d_GE5_onesided.csv"))
