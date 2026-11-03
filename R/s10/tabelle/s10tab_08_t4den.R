
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
# Follow-up: Table 4 denominators. Per design point x label: all replicates, generated, estimable (non-degenerate)
R <- readRDS(file.path(OUT, "s1_recs.rds"))
D <- R[, .(nrep = .N, ngen = sum(generated %in% TRUE), nd = sum(estim), ndeg = sum(!estim), ngen_deg = sum(generated %in% TRUE & !estim)), by = KD]
fwrite(D, file.path(OUT, "t4_den_counts.csv"))
print(D[!(scen %in% c("M6","M7","M1c")), .(nrep = sum(nrep), ngen = sum(ngen), nd = sum(nd), ndeg = sum(ndeg), ngen_deg = sum(ngen_deg)), by = .(dgp, label)])
