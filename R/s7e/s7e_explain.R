
source(path.expand("~/MMM/R/s7e/s7e_common.R")); data.table::setDTthreads(8L)
TK <- s7e_tasks(); L <- parallel::mclapply(TK$file, readRDS, mc.cores = 12L)
RECS <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); ST <- rbindlist(lapply(L, `[[`, "study"), fill = TRUE); rm(L)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
P <- CJ(k = seq_len(nrow(RECS)), method = TECHS_S7E); P <- cbind(RECS[P$k, c(K6, "design_degenerate"), with = FALSE], method = P$method)
have <- unique(ST[method %in% TECHS_S7E, c(K6, "method"), with = FALSE]); have[, h := TRUE]
P <- merge(P, have, by = c(K6, "method"), all.x = TRUE)
G <- P[is.na(h), .(missing = .N, degenerate = sum(design_degenerate)), by = .(method, scen, N, ntime)]
fwrite(G, file.path(TABDIR, "V4_S7e_GF1_missing.csv")); print(G[, .(missing = sum(missing), degenerate = sum(degenerate)), by = method])
print(P[is.na(h), .N, by = .(sim_id, scen, N, ntime, label)][, .N, by = .(n_methods_missing = N)])
D <- fread(file.path(TABDIR, "V4_S7e_dp_rates.csv"))[ref == "t" & strategy %in% c("gee_exch_fg", "gee_ind_fg", "gee_ar1_fg", "mmm_cc") ]
W <- dcast(D, scen + N + ntime ~ strategy, value.var = c("rate", "hi"))
print(W[hi_gee_exch_fg | hi_gee_ind_fg, .(scen, N, ntime, rate_gee_exch_fg, rate_gee_ind_fg, rate_gee_ar1_fg, rate_mmm_cc)], digits = 3)
print(D[strategy == "gee_exch_fg", .(mean_rate = mean(rate, na.rm=TRUE)), by = ntime][order(ntime)], digits = 3)
print(D[strategy == "gee_ar1_fg", .(mean_rate = mean(rate, na.rm=TRUE)), by = ntime][order(ntime)], digits = 3)
