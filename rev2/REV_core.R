
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table)})
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
cat("rows", nrow(R), "\n")
H0 <- R[label=="manual" & set=="own"]
cat("design points H0:", uniqueN(H0[,.(dgp,model,N,ntime)]), "\n")
cat("methods:", paste(sort(unique(H0$method)),collapse=","), "\n")

# ---- 3.1 : replicate-weighted rate within each N, worst case, MAD over the 16 N
byN <- H0[, .(n=sum(n), k=sum(k_rej)), by=.(dgp,method,N)][, rate := 100*k/n]
pooled <- byN[, .(rate = 100*sum(k)/sum(n),
                  worst = max(abs(rate-5)),
                  mad = mean(abs(rate-5)),
                  argmax_N = N[which.max(abs(rate-5))]), by=.(dgp,method)]
fwrite(pooled[order(dgp,mad)], "REV_31_pooled.csv")
print(pooled[order(dgp,mad)], nrows=40)

# ---- MAD by T (plain mean over design points)
byT <- H0[, .(rate=100*k_rej/n, dev=abs(100*k_rej/n - 5)), by=.(dgp,method,model,N,ntime)]
tstr <- byT[, .(mad = mean(dev), cells=.N), by=.(method, Tstrat = ifelse(ntime==2,"T2", ifelse(ntime==3,"T3","T4plus")))]
tw <- dcast(tstr, method ~ Tstrat, value.var="mad")
fwrite(tw, "REV_31_byT.csv"); print(tw)

# ---- cells won (closest to nominal per design point)
dp <- H0[, .(dev = abs(100*k_rej/n - 5)), by=.(dgp,model,N,ntime,method)]
won <- dp[, .(winner = method[which.min(dev)], nties = sum(dev==min(dev))), by=.(dgp,model,N,ntime)]
cat("\n--- cells won\n"); print(won[, .N, by=winner][order(-N)])
cat("design points:", nrow(won), " ties at minimum:", won[nties>1,.N], "\n")
fwrite(won, "REV_31_won.csv")

# ---- intersection shift
ish <- R[label=="manual", .(n=sum(n), k=sum(k_rej)), by=.(dgp,method,set)][, rate:=100*k/n]
iw <- dcast(ish, dgp+method ~ set, value.var=c("rate","n"))
iw[, shift := rate_intersection - rate_own]
fwrite(iw, "REV_33_intersection.csv"); print(iw[order(dgp,shift)], nrows=30)
cat("\nmatched replicates in intersection (all dgp):", R[label=="manual" & set=="intersection", sum(n)]/8, "\n")

# ---- empty intersection design points
ei <- unique(R[label=="manual", .(dgp,model,N,ntime,n_inter)])
cat("empty-intersection design points:", ei[n_inter==0,.N], " by ntime:\n"); print(ei[n_inter==0, .N, by=ntime])

# ---- usability / estimability
u <- unique(R[label=="manual", .(dgp,model,N,ntime,n_attempt)])
uu <- R[label=="manual" & set=="own", .(dgp,model,N,ntime,method,n_usable,n_present,n_attempt)]
# fill in method-cells with zero usable: expand grid
grid <- CJ(dgp=unique(u$dgp), model=unique(u$model), N=unique(u$N), ntime=unique(u$ntime),
           method=unique(uu$method), unique=TRUE)
grid <- merge(grid, u, by=c("dgp","model","N","ntime"), all.x=TRUE)
full <- merge(grid, uu[,.(dgp,model,N,ntime,method,n_usable)], by=c("dgp","model","N","ntime","method"), all.x=TRUE)
full[is.na(n_usable), n_usable := 0]
full[, pu := 100*n_usable/n_attempt]
est <- full[, .(pu = mean(pu)), by=.(method, ntime, slope = ifelse(model<=2,"RI","RS"))]
fwrite(dcast(est, method+slope ~ ntime, value.var="pu"), "REV_38_estimability.csv")
print(dcast(est, method+slope ~ ntime, value.var="pu"), nrows=30)
cat("\nGLMM Wald usable count at T=2, model>=3:", full[method=="glmer" & ntime==2 & model>=3, sum(n_usable)],
    "of", full[method=="glmer" & ntime==2 & model>=3, sum(n_attempt)], "\n")
cat("shrink:", full[method=="shrink" & ntime==2 & model>=3, sum(n_usable)], "\n")
cat("LRT usable T=2 RS:", full[method=="glmer_lrt" & ntime==2 & model>=3, mean(pu)], "\n")

# ---- unconditional validity at N<=20 (all 270 design points, unfitted at zero)
uv <- merge(full[N<=20], H0[,.(dgp,model,N,ntime,method,n,k_rej)], by=c("dgp","model","N","ntime","method"), all.x=TRUE)
uv[, corr := ifelse(is.na(n) | n==0, 0, 1 - k_rej/n)]
uv[, valid := (pu/100)*corr*100]
res <- uv[, .(usable=mean(pu), validity=mean(valid), cells=sum(n_usable>0)), by=method][order(-validity)]
fwrite(res, "REV_39_validity.csv"); print(res)
cat("design points at N<=20 per method:", uv[, .N, by=method][1]$N, "\n")

# ---- whole grid usable
cat("\nwhole-grid usable %:\n"); print(full[, .(pu=mean(pu)), by=method][order(-pu)])
