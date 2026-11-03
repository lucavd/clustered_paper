
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
ALL <- readRDS(file.path(OUT, "s4_all_dp.rds")); fwrite(ALL, file.path(OUT, "s4_all_dp.csv")); cat(nrow(ALL), "\n")
