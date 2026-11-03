
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table)})
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
H0 <- R[label=="manual" & set=="own"]
H1 <- R[label=="mmlb"  & set=="own"]
att <- unique(H0[,.(dgp,model,N,ntime,n_attempt)])
cat("att rows:", nrow(att), " design points:", uniqueN(H0[,.(dgp,model,N,ntime)]), "\n")
if (nrow(att) != 1440) {
  cat("n_attempt not unique per design point -- taking max\n")
  att <- H0[, .(n_attempt = max(n_attempt)), by=.(dgp,model,N,ntime)]
}
grid <- CJ(dgp=unique(H0$dgp), model=1:5, N=sort(unique(H0$N)), ntime=sort(unique(H0$ntime)),
           method=sort(unique(H0$method)))
grid <- merge(grid, att, by=c("dgp","model","N","ntime"), all.x=TRUE)
grid <- merge(grid, H0[,.(dgp,model,N,ntime,method,n_usable,n,k_rej,k_cov,n_cov)],
              by=c("dgp","model","N","ntime","method"), all.x=TRUE)
grid[is.na(n_usable), n_usable := 0L]
grid[, pu := n_usable/n_attempt]
cat("grid rows:", nrow(grid), "\n")

MARG5 <- c("mmm","gee","gee_fg","gee_md","gee_kc")
diff <- grid[method %in% MARG5, .(pu_marg = mean(pu)), by=.(dgp,model,N,ntime)]
diff[, stratum := cut(pu_marg, breaks=c(-Inf,0.5,0.8,0.95,0.999,Inf),
                      labels=c("<50","50-80","80-95","95-99.9","100"), right=TRUE)]
cat("\n=== design points per stratum ===\n"); print(diff[, .N, by=stratum][order(stratum)])
gd <- merge(grid, diff[,.(dgp,model,N,ntime,stratum)], by=c("dgp","model","N","ntime"))
gd[, dev := ifelse(is.na(n)|n==0, NA_real_, abs(100*k_rej/n - 5))]
st <- gd[, .(mad=mean(dev, na.rm=TRUE)), by=.(method,stratum)]
st[, rank := frank(mad), by=stratum]
cat("\n=== MAD by stratum ===\n"); print(dcast(st, method ~ stratum, value.var="mad"))
cat("\n=== rank by stratum ===\n"); print(dcast(st, method ~ stratum, value.var="rank"))
fwrite(st, "REV_strata.csv")
w <- dcast(gd[method %in% c("mmm","gee_fg")], dgp+model+N+ntime+stratum ~ method, value.var="dev")
w <- w[is.finite(mmm) & is.finite(gee_fg)][, gap := mmm - gee_fg]
cat("\n=== signed gap MMM-FG, paired Wilcoxon ===\n")
print(w[, .(gap=mean(gap), n=.N, p=wilcox.test(mmm, gee_fg, paired=TRUE)$p.value), by=stratum][order(stratum)])
ws <- merge(w, diff[,.(dgp,model,N,ntime,pu_marg)], by=c("dgp","model","N","ntime"))
ct <- cor.test(ws$pu_marg, ws$gap, method="spearman")
cat(sprintf("Spearman(pu,gap) rho=%.4f p=%.4g n=%d\n", ct$estimate, ct$p.value, nrow(ws)))

cat("\n=== t REFERENCE ===\n")
tr <- H0[method!="glmer_lrt" & !is.na(k_rej_t), .(n=sum(n), kn=sum(k_rej), kt=sum(k_rej_t)), by=.(dgp,method,N)]
tr[, `:=`(norm=100*kn/n, tref=100*kt/n)]
pool <- tr[, .(norm=100*sum(kn)/sum(n), tref=100*sum(kt)/sum(n),
               wc_norm=max(abs(norm-5)), wc_t=max(abs(tref-5)),
               argmax_t=N[which.max(abs(tref-5))]), by=.(dgp,method)][, shift := tref-norm]
print(pool[order(dgp,wc_t)], nrows=30)
cat("\nshift at N<=25 (pooled over dgp):\n")
print(tr[N<=25, .(shift=100*sum(kt)/sum(n)-100*sum(kn)/sum(n)), by=method][order(shift)])
cat("\nKC under t, N=10 and 145:\n")
print(tr[method=="gee_kc" & N %in% c(10,145), .(norm=100*sum(kn)/sum(n), tref=100*sum(kt)/sum(n)), by=N])
fwrite(tr, "REV_tref_byN.csv")

cat("\n=== COMBINED RANKING (level+coverage+usability) ===\n")
mk <- function(nsub, csub, gsub) {
  lev <- nsub[, .(lev=abs(100*sum(k_rej)/sum(n)-5)), by=.(dgp,method)]
  cov <- csub[n_cov>0, .(cov=abs(100*sum(k_cov)/sum(n_cov)-95)), by=.(dgp,method)]
  usb <- gsub[, .(usable=mean(pu)), by=.(dgp,method)]
  x <- merge(merge(lev,cov,by=c("dgp","method")), usb, by=c("dgp","method"))
  x[, `:=`(r_lev=frank(lev), r_cov=frank(cov), r_use=frank(-usable)), by=dgp]
  x[, mean_rank := (r_lev+r_cov+r_use)/3]
  x[, .(mean_rank=round(mean(mean_rank),2)), by=method][order(mean_rank)]
}
cat("whole grid:\n"); print(mk(H0, H1, grid))
cat("N<=25:\n"); print(mk(H0[N<=25], H1[N<=25], grid[N<=25]))

cat("\n=== COVERAGE (nominal 95) ===\n")
cv <- H1[n_cov>0, .(cov=round(100*sum(k_cov)/sum(n_cov),2)), by=.(dgp,method)]
print(dcast(cv, method ~ dgp, value.var="cov"))

cat("\n=== FG / MD deviation above N=45 ===\n")
fgn <- H0[method %in% c("gee_fg","gee_md"), .(n=sum(n), k=sum(k_rej)), by=.(dgp,method,N)][, dev := round(abs(100*k/n-5),3)]
print(dcast(fgn[N>=45], N ~ dgp+method, value.var="dev"), nrows=20)
cat("\nMMM deviation above N=45:\n")
mmn <- H0[method=="mmm", .(n=sum(n), k=sum(k_rej)), by=.(dgp,N)][, dev := round(abs(100*k/n-5),3)]
print(dcast(mmn[N>=45], N ~ dgp, value.var="dev"))
cat("\n=== N=10 gee detail (sec 3.5) ===\n")
print(H0[method=="gee" & N==10, .(n=sum(n), k=sum(k_rej)), by=dgp][, `:=`(rate=100*k/n, mcse=100*sqrt((k/n)*(1-k/n)/n))])
cat("\nN=145 gee:\n")
print(H0[method=="gee" & N==145, .(n=sum(n), k=sum(k_rej)), by=dgp][, `:=`(rate=100*k/n, z=(k/n-0.05)/sqrt(0.05*0.95/n))])
cat("\ngee minima by dgp over N:\n")
gmin <- H0[method=="gee", .(n=sum(n), k=sum(k_rej)), by=.(dgp,N)][, rate := 100*k/n]
print(gmin[, .SD[which.min(rate)], by=dgp])
