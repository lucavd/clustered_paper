
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
V <- readRDS("results_v4/variant_kcfc/S4/rows_v4.rds")[!is.na(term) & scen %in% sprintf("M%d", 1:5) & label == "manual"]
print(V[method %in% c("glmer", "gee_fg") & scen == "M3", .(rows = .N, reps = uniqueN(paste(dgp, N, ntime, sim_id)), terms = paste(unique(term), collapse = "/")), by = method])
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
DEG <- fread("tables/V4_S2_degenerate_designs.csv")[, .(model, N, ntime, sim_id)]
V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]; V[, u0 := U0 & !deg]
NDt <- unique(V[method == "gee", .(dgp, scen, N, ntime, sim_id, deg)])[, .(nd = sum(!deg)), by = .(dgp, scen, N, ntime)]
G0 <- CJ(method = unique(V$method), k = NDt[, paste(dgp, scen, N, ntime)])   # every cell, also those where a method returned no row
NDt[, k := paste(dgp, scen, N, ntime)]; G0 <- merge(G0, NDt, by = "k")
P <- merge(G0, V[, .(n = sum(u0)), by = .(method, dgp, scen, N, ntime)], by = c("method", "dgp", "scen", "N", "ntime"), all.x = TRUE)[is.na(n), n := 0L]
P[, grp := fifelse(scen %in% c("M1", "M2"), "RI", "RS")]
Q <- P[, .(pu = mean(100 * n / nd)), by = .(grp, method, ntime)]
TC <- fread("tables/variant_kcfc/V4_S6_table_cells.csv")[variant == "main" & table_id == "T7"]
TC[, c("q", "grp", "m", "t") := tstrsplit(name, "|", fixed = TRUE)]; TC[, t := as.integer(t)]
J <- merge(TC[, .(m, grp, t, est)], Q, by.x = c("m", "grp", "t"), by.y = c("method", "grp", "ntime"))
print(J[, .(cells = .N, max_abs_diff = max(abs(pu - est))), by = m])
fwrite(J, "tables/V5_S7f_table3_recompute.csv")
