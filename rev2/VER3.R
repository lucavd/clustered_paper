
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
H <- fread(path.expand("~/MMM/rev2/REV_H1_targets.csv"))
G <- H[method=="glmer"]
cat("=== A. T-3.13-04 (27,2%) contro M-14 (23,8%): bias ASSOLUTO medio per cella ===\n")
G[, `:=`(mae_marg=sae_marg/n, mae_cond=sae_cond/n)]
cat(sprintf("media per cella:  marg %.5f  cond %.5f  riduzione %.4f  (celle=%d)\n",
  mean(G$mae_marg), mean(G$mae_cond), 1-mean(G$mae_cond)/mean(G$mae_marg), nrow(G)))
cat(sprintf("pesata per repliche: marg %.5f cond %.5f riduzione %.4f\n",
  sum(G$sae_marg)/sum(G$n), sum(G$sae_cond)/sum(G$n), 1-sum(G$sae_cond)/sum(G$sae_marg)))
cat(sprintf("|bias| del bias medio per cella: marg %.5f cond %.5f riduzione %.4f\n",
  mean(abs(G$bias_marg)), mean(abs(G$bias_cond)), 1-mean(abs(G$bias_cond))/mean(abs(G$bias_marg))))
cat("\n=== B. MCSE della coverage: 93,16 contro 93,19 e' distinguibile? ===\n")
for (dg in c("GenBinaryY","GLMM-native","GLMM-AR1")) {
  x <- G[dgp==dg]; n <- sum(x$n); p <- sum(x$k_cond)/n
  cat(sprintf("%-12s n=%d  cov_cond=%.4f%%  MCSE=%.4f pp\n", dg, n, 100*p, 100*sqrt(p*(1-p)/n)))
}
cat("t_cond per AR1 in REV_H1_targets (per model x ntime):\n")
print(dcast(unique(G[dgp=="GLMM-AR1",.(model,ntime,t_cond)]), model~ntime, value.var="t_cond"))
A <- fread(path.expand("~/MMM/rev2/R2_ar1_cond_true.csv"))
print(dcast(A, model~ntime, value.var="t_cond_true"))
cat("\n=== C. T-3.14-04: RMSE dello shrink, media per cella ===\n")
H[, rmse_marg := sqrt(sse_marg/n)]
S <- H[, .(rmse_mean=mean(rmse_marg), ncell=.N, rmse_pooled=sqrt(sum(sse_marg)/sum(n))), by=.(method,dgp)]
print(dcast(S[method %in% c("shrink","glmer","mmm","gee_fg","gee_kc","gee_md")], method~dgp,
            value.var=c("rmse_mean","ncell")))
cat("\nspread relativo (max-min)/min fra i quattro marginali, media per cella:\n")
q <- S[method %in% c("mmm","gee_fg","gee_kc","gee_md"), .(spread_rel=100*(max(rmse_mean)/min(rmse_mean)-1),
      spread_abs=100*(max(rmse_mean)-min(rmse_mean)), argmin=method[which.min(rmse_mean)]), by=dgp]
print(q)
cat("\n=== D. R2-01: identita, differenza fra le due letture ===\n")
P <- fread(path.expand("~/MMM/rev2/PAIRED_mmm_fg_cells.csv"))
pp <- P[, .(npair=sum(npair), b=sum(b), cc=sum(cc), km=sum(km), kf=sum(kf)), by=dgp]
pp[, `:=`(gap_pair=100*(km-kf)/npair, gap_bc=100*(b-cc)/npair)]
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
o <- dcast(R[label=="manual" & set=="own" & method %in% c("mmm","gee_fg"),
             .(rate=100*sum(k_rej)/sum(n)), by=.(dgp,method)], dgp~method, value.var="rate")
m <- merge(pp[,.(dgp,npair,gap_pair,gap_bc)], o, by="dgp")
m[, diff_own := mmm-gee_fg][, discrepanza := diff_own-gap_pair]
print(m)
cat("\nmax |gap_pair - gap_bc| (identita b-c):", max(abs(pp$gap_pair-pp$gap_bc)), "\n")
cat("\n=== E. D2-31r: worst case pooled per N (non per cella) sotto i cinque gdl ===\n")
DF <- fread(path.expand("~/MMM/rev2/DFREF_cells.csv"))
byN <- DF[, .(n=sum(n), r_norm=100*sum(k_norm)/sum(n), r_bw=100*sum(k_bw)/sum(n),
              r_k2=100*sum(k_k2)/sum(n), r_win=100*sum(k_win)/sum(n), r_bwin=100*sum(k_bwin)/sum(n)),
          by=.(dgp,method,N)]
WCN <- byN[, .(w_norm=max(abs(r_norm-5)), w_bw=max(abs(r_bw-5)), w_k2=max(abs(r_k2-5)),
               w_win=max(abs(r_win-5)), w_bwin=max(abs(r_bwin-5))), by=.(dgp,method)]
print(WCN[method %in% c("mmm","gee_fg")][order(dgp,method)])
for (dg in unique(WCN$dgp)) for (cc in c("w_norm","w_bw","w_k2","w_win","w_bwin"))
  cat(sprintf("%-12s %-7s primato -> %-7s %.3f\n", dg, cc,
      WCN[dgp==dg][which.min(get(cc))]$method, min(WCN[dgp==dg][[cc]])))
cat("\n=== F. §6.1: 864 punti di disegno a pendenza casuale ===\n")
cat("punti di disegno con model in 3,4,5:",
    nrow(unique(R[model %in% c(3,4,5), .(dgp,model,N,ntime)])), " totali:",
    nrow(unique(R[, .(dgp,model,N,ntime)])), "\n")
cat("\n=== G. nobs per cluster: la regola 70/30 ===\n")
cat("attesi a T=4: 0.7*4+0.3*mean(1:3) =", 0.7*4+0.3*mean(1:3), " -> N=10 danno", 10*(0.7*4+0.3*2), "\n")
cat("df_win e df_bwin nel file a N=10,T=4:", unique(DF[N==10 & ntime==4 & model==1]$df_win),
    unique(DF[N==10 & ntime==4 & model==1]$df_bwin), "\n")
cat("\nDONE3\n")
