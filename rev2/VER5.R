
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
H1 <- R[label=="mmlb" & set=="own"]
cat("=== L3: si puo' fare §3.14 (RMSE, media per cella) da REV_cells.csv? ===\n")
H1[, rmse := sqrt(sd_est^2*(n-1)/n + (mean_est-truth)^2)]
S <- H1[method %in% c("mmm","gee_fg","gee_kc","gee_md","glmer","shrink"),
        .(rmse=mean(rmse), ncell=.N), by=.(method,dgp)]
print(dcast(S, method~dgp, value.var="rmse"))
q <- S[method %in% c("mmm","gee_fg","gee_kc","gee_md"), .(spread_pp=100*(max(rmse)-min(rmse))), by=dgp]
cat("spread da REV_cells.csv:\n"); print(q)
H <- fread(path.expand("~/MMM/rev2/REV_H1_targets.csv"))
H[, rmse2 := sqrt(sse_marg/n)]
S2 <- H[method %in% c("mmm","gee_fg","gee_kc","gee_md"), .(rmse=mean(rmse2)), by=.(method,dgp)]
cat("spread da REV_H1_targets.csv (fonte chunk):\n")
print(S2[, .(spread_pp=100*(max(rmse)-min(rmse))), by=dgp])
cat("\n=== L3: si puo' fare la coverage contro il condizionale da REV_cells.csv? ===\n")
cat("colonne di coverage disponibili:", paste(grep("cov", names(R), value=TRUE), collapse=","), "\n")
cat("valori distinti di truth per dgp/model (e' il marginale?):\n")
print(R[, .(truth=paste(sort(unique(truth)),collapse="/")), by=.(dgp,model)][order(dgp,model)])
cat("coverage marginale ricalcolabile:\n")
print(H1[method=="glmer", .(cov=round(100*sum(k_cov)/sum(n_cov),2)), by=dgp])
cat("\n=== L3: divario appaiato R2-01 da REV_cells.csv? ===\n")
cat("colonne per-replica presenti:", any(grepl("sim_id|estimate|^p$", names(R))), "\n")
cat("\nDONE5\n")
