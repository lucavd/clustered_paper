
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
for (f in c("s3_ns.rds", "s3_cop.rds")) { a <- readRDS(file.path(OUT, f))$agree; print(f); print(a[, lapply(.SD, sum), .SDcols = c("n_cc", "bnd_disagree", "used_disagree", "used_proj_na"), by = .(dgp)]) }
