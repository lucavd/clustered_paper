
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table)})
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
H0 <- R[label=="manual" & set=="own"]
H1 <- R[label=="mmlb"  & set=="own"]

## ---------- 1. multiplicity
cat("=== MULTIPLICITY ===\n")
cells <- H0[n>0]
cat("method-cells with >=1 usable replicate:", nrow(cells), "\n")
cat("of 1440*8 =", 1440*8, " -> missing:", 1440*8-nrow(cells), "\n")
miss <- CJ(dgp=unique(H0$dgp), model=1:5, N=unique(H0$N), ntime=unique(H0$ntime), method=unique(H0$method))
miss <- merge(miss, cells[,.(dgp,model,N,ntime,method,n)], all.x=TRUE)
cat("missing by method:\n"); print(miss[is.na(n), .N, by=method])
cat("missing by ntime:\n"); print(miss[is.na(n), .N, by=ntime])
cells[, p_bin := mapply(function(k,nn) binom.test(k, nn, 0.05)$p.value, k_rej, n)]
cat("raw <0.05:", round(100*mean(cells$p_bin<0.05),2), "%\n")
cells[, bh_pool := p.adjust(p_bin, "BH")]
cells[, bonf_pool := p.adjust(p_bin, "bonferroni")]
cat("BH pooled:", round(100*mean(cells$bh_pool<0.05),2), "%\n")
cat("Bonferroni pooled:", round(100*mean(cells$bonf_pool<0.05),2), "%\n")
cells[, bh_meth := p.adjust(p_bin,"BH"), by=method]
cells[, bonf_meth := p.adjust(p_bin,"bonferroni"), by=method]
cat("Bonferroni within-method:", round(100*mean(cells$bonf_meth<0.05),2), "%\n")
tab <- cells[, .(bh_within=round(100*mean(bh_meth<0.05),1), bh_pooled=round(100*mean(bh_pool<0.05),1),
                 ncell=.N), by=method][order(bh_within)]
print(tab)
fwrite(cells, "REV_multiplicity_cells.csv")

## ---------- 2. stratification by estimability
cat("\n=== STRATIFICATION ===\n")
MARG5 <- c("mmm","gee","gee_fg","gee_md","gee_kc")
att <- unique(H0[,.(dgp,model,N,ntime,n_attempt)])
grid <- CJ(dgp=unique(H0$dgp), model=1:5, N=unique(H0$N), ntime=unique(H0$ntime), method=unique(H0$method))
grid <- merge(grid, att, by=c("dgp","model","N","ntime"), all.x=TRUE)
grid <- merge(grid, H0[,.(dgp,model,N,ntime,method,n_usable,n,k_rej)], all.x=TRUE)
grid[is.na(n_usable), n_usable := 0]
grid[, pu := n_usable/n_attempt]
diff <- grid[method %in% MARG5, .(pu_marg = mean(pu)), by=.(dgp,model,N,ntime)]
diff[, stratum := cut(pu_marg, breaks=c(-Inf,0.5,0.8,0.95,0.999,Inf),
                      labels=c("<50","50-80","80-95","95-99.9","100"), right=TRUE)]
cat("design points per stratum:\n"); print(diff[, .N, by=stratum][order(stratum)])
gd <- merge(grid, diff[,.(dgp,model,N,ntime,stratum)], by=c("dgp","model","N","ntime"))
gd[, dev := ifelse(is.na(n)|n==0, NA_real_, abs(100*k_rej/n - 5))]
st <- gd[, .(mad=mean(dev, na.rm=TRUE)), by=.(method,stratum)]
st[, rank := frank(mad), by=stratum]
print(dcast(st, method ~ stratum, value.var="mad"))
print(dcast(st, method ~ stratum, value.var="rank"))
fwrite(st, "REV_strata.csv")
# signed MMM - FG gap per stratum + Wilcoxon paired
w <- dcast(gd[method %in% c("mmm","gee_fg")], dgp+model+N+ntime+stratum ~ method, value.var="dev")
w <- w[is.finite(mmm) & is.finite(gee_fg)]
w[, gap := mmm - gee_fg]
cat("\nsigned gap MMM-FG by stratum, with paired Wilcoxon:\n")
print(w[, .(gap=mean(gap), n=.N,
            p=wilcox.test(mmm, gee_fg, paired=TRUE)$p.value), by=stratum][order(stratum)])
cat("Spearman(pu_marg, gap):\n")
ws <- merge(w, diff[,.(dgp,model,N,ntime,pu_marg)], by=c("dgp","model","N","ntime"))
ct <- cor.test(ws$pu_marg, ws$gap, method="spearman")
cat(sprintf("rho=%.4f p=%.3g n=%d\n", ct$estimate, ct$p.value, nrow(ws)))

## ---------- 3. t reference, recomputed from k_rej_t
cat("\n=== t REFERENCE ===\n")
tr <- H0[method!="glmer_lrt" & !is.na(k_rej_t), .(n=sum(n), kn=sum(k_rej), kt=sum(k_rej_t)), by=.(dgp,method,N)]
tr[, `:=`(norm=100*kn/n, tref=100*kt/n)]
pool <- tr[, .(norm=100*sum(kn)/sum(n), tref=100*sum(kt)/sum(n),
               wc_norm=max(abs(norm-5)), wc_t=max(abs(tref-5)),
               argmax_t_N=N[which.max(abs(tref-5))]), by=.(dgp,method)]
pool[, shift := tref-norm]
print(pool[order(dgp,wc_t)], nrows=30)
cat("\nshift at N<=25:\n")
print(tr[N<=25, .(shift=100*sum(kt)/sum(n) - 100*sum(kn)/sum(n)), by=method][order(shift)])
fwrite(tr, "REV_tref_byN.csv")
cat("\nKC under t at N=10 and N=145 (replicate-weighted across dgp):\n")
print(tr[method=="gee_kc" & N %in% c(10,145), .(tref=100*sum(kt)/sum(n), norm=100*sum(kn)/sum(n)), by=N])

## ---------- 4. combined ranking: level + coverage + usability
cat("\n=== COMBINED RANKING ===\n")
lev <- H0[, .(lev = abs(100*sum(k_rej)/sum(n) - 5)), by=.(dgp,method)]
cov <- H1[n_cov>0, .(cov = abs(100*sum(k_cov)/sum(n_cov) - 95)), by=.(dgp,method)]
usb <- grid[, .(usable = mean(pu)), by=.(dgp,method)]
cmb <- merge(merge(lev, cov, by=c("dgp","method")), usb, by=c("dgp","method"))
cmb[, `:=`(r_lev=frank(lev), r_cov=frank(cov), r_use=frank(-usable)), by=dgp]
cmb[, mean_rank := (r_lev+r_cov+r_use)/3]
print(cmb[, .(mean_rank=mean(mean_rank)), by=method][order(mean_rank)])
cat("\n--- restricted to N<=25\n")
lev2 <- H0[N<=25, .(lev = abs(100*sum(k_rej)/sum(n) - 5)), by=.(dgp,method)]
cov2 <- H1[N<=25 & n_cov>0, .(cov = abs(100*sum(k_cov)/sum(n_cov) - 95)), by=.(dgp,method)]
usb2 <- grid[N<=25, .(usable = mean(pu)), by=.(dgp,method)]
c2 <- merge(merge(lev2, cov2, by=c("dgp","method")), usb2, by=c("dgp","method"))
c2[, `:=`(r_lev=frank(lev), r_cov=frank(cov), r_use=frank(-usable)), by=dgp]
c2[, mean_rank := (r_lev+r_cov+r_use)/3]
print(c2[, .(mean_rank=mean(mean_rank)), by=method][order(mean_rank)])
fwrite(rbind(cmb[,scope:="all"], c2[,scope:="Nle25"]), "REV_combined_ranking.csv")

## ---------- 5. coverage table
cat("\n=== COVERAGE (nominal 95, generating parameter) ===\n")
cv <- H1[n_cov>0, .(cov=100*sum(k_cov)/sum(n_cov), n=sum(n_cov)), by=.(dgp,method)]
print(dcast(cv, method ~ dgp, value.var="cov")[order(-GenBinaryY)])

## ---------- 6. FG above 45 clusters: "below 0.1 pp and stays there"
cat("\n=== FG and MD deviation by N above 45 ===\n")
fgn <- H0[method %in% c("gee_fg","gee_md"), .(n=sum(n), k=sum(k_rej)), by=.(dgp,method,N)][, dev := abs(100*k/n-5)]
print(dcast(fgn[N>=45], N ~ dgp+method, value.var="dev"), nrows=20)
cat("\nMMM deviation by N above 45:\n")
mmn <- H0[method=="mmm", .(n=sum(n), k=sum(k_rej)), by=.(dgp,N)][, dev := abs(100*k/n-5)]
print(dcast(mmn[N>=45], N ~ dgp, value.var="dev"))
