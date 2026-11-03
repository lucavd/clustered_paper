
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table)})
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
H0 <- R[label=="manual" & set=="own"]
d <- H0[, .N, by=.(dgp,model,N,ntime,method)][N>1]
cat("duplicate keys:", nrow(d), "\n"); print(head(d,10))
cat("\nrows for one duplicate:\n")
if (nrow(d)) print(H0[dgp==d$dgp[1] & model==d$model[1] & N==d$N[1] & ntime==d$ntime[1] & method==d$method[1]])
cat("\nunique n_attempt per design point?\n")
print(H0[, .(k=uniqueN(n_attempt)), by=.(dgp,model,N,ntime)][k>1][1:5])
