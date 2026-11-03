
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2)
D <- "/home/user/MMM/results_v4/S10/supplemento"
TOm <- fread(file.path(D, "s10s_timeouts.csv"))[method == "gee_md"][, run := "study rows (main grid, M1c)"]
NS <- readRDS("/home/user/MMM/results_v4/S6b/ns/rows_ns.rds")[method == "gee_md"]
NS[, to := grepl("timeout|time-out|timed out", status, ignore.case = TRUE)]
TOn <- NS[, .(n = .N, pct_to = 100 * mean(to)), by = .(dgp, scen, N, ntime, label)][, `:=`(method = "gee_md", run = "S6b rows (M6-M7)")]
E <- readRDS(file.path(D, "study_s7e.rds"))[method == "gee_md"]
E[, to := grepl("timeout|time-out|timed out", status, ignore.case = TRUE)]
TOe <- E[, .(n = .N, pct_to = 100 * mean(to)), by = .(dgp, scen, N, ntime, label)][, `:=`(method = "gee_md", run = "S7e rows (Copula-AR1)")]
A <- rbind(TOm, TOn, TOe, use.names = TRUE)
fwrite(A[order(-pct_to)], file.path(D, "s10s_timeouts_all.csv"))
print(A[, .(cells = .N, cells_with_to = sum(pct_to > 0), max_pct = max(pct_to), fits = sum(n), to = round(sum(n * pct_to / 100))), by = run])
print(A[order(-pct_to)][1:6])
print(table(unique(E$status[E$to])))
