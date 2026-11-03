
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
D <- fread(path.expand("~/MMM/rev2/DFREF_cells.csv"))
cat("colonne:", paste(names(D), collapse=","), "\n righe:", nrow(D), "\n")
NP <- c(4,5,4,4,5); NPW <- c(2,2,2,2,2)
D[, nobs_exp := N*(0.7*ntime + 0.3*ntime/2)]
D[, `:=`(f_bw = N - NP[model], f_win = nobs_exp - NP[model],
         f_bwin = nobs_exp - N - NP[model] + 1, f_bw_rule = nobs_exp - N - NPW[model])]
cat("\n=== V2-08: le formule dei gdl nel file ===\n")
cat("df_bw  == N - p                 :", all(D$df_bw==D$f_bw), "\n")
cat("df_win == nobs_exp - p          :", all(abs(D$df_win-round(D$f_win))<=1), "\n")
cat("df_bwin== nobs_exp - N - p + 1  :", all(abs(D$df_bwin-round(D$f_bwin))<=1), "\n")
cat("\ncaso N=10, T=4 (il caso citato 'a dieci cluster e quattro occasioni'):\n")
print(unique(D[N==10 & ntime==4, .(model, nobs_exp, df_bw, df_win, df_bwin,
     regola_bw = f_bw_rule)])[order(model)])
cat("\n=== worst case per CELLA (df_bwin) ===\n")
D[, rate_bwin := 100*k_bwin/n]
wc_cell <- D[method %in% c("mmm","gee_fg"), .(worst=max(abs(rate_bwin-5))), by=.(dgp,method)]
print(dcast(wc_cell, method~dgp, value.var="worst"), digits=5)
cat("\n=== worst case aggregati per N (poolando le occasioni) ===\n")
byN <- D[method %in% c("mmm","gee_fg"), .(n=sum(n), k=sum(k_bwin)), by=.(dgp,method,N)]
byN[, rate := 100*k/n]
wc_N <- byN[, .(worst=max(abs(rate-5))), by=.(dgp,method)]
print(dcast(wc_N, method~dgp, value.var="worst"), digits=5)
cat("\nMMM per meccanismo nell'ordine GenBinaryY / GLMM-native / GLMM-AR1:\n")
m <- wc_N[method=="mmm"]
for (d in c("GenBinaryY","GLMM-native","GLMM-AR1")) cat(sprintf("  %-12s %.4f\n", d, m[dgp==d, worst]))
cat("\nFG per meccanismo:\n")
fg <- wc_N[method=="gee_fg"]
for (d in c("GenBinaryY","GLMM-native","GLMM-AR1")) cat(sprintf("  %-12s %.4f\n", d, fg[dgp==d, worst]))
cat("\nDONE_JOB4\n")
