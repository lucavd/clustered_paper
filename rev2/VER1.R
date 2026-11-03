
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
H <- fread(path.expand("~/MMM/rev2/REV_H1_targets.csv"))
cat("== REV_H1_targets: t_cond per dgp x ntime (glmer) ==\n")
print(dcast(unique(H[method=="glmer",.(dgp,model,ntime,t_cond)])[model==1], dgp~ntime, value.var="t_cond"))
G <- H[method=="glmer"]
cat("\n== §3.13 BIAS del glmer (dal file dei bersagli) ==\n")
cat(sprintf("pesato per repliche: cond %.4f  marg %.4f   (n=%d, celle=%d)\n",
  sum(G$bias_cond*G$n)/sum(G$n), sum(G$bias_marg*G$n)/sum(G$n), sum(G$n), nrow(G)))
cat(sprintf("media per cella:     cond %.4f  marg %.4f\n", mean(G$bias_cond), mean(G$bias_marg)))
cat(sprintf("mediana per cella:   cond %.4f  marg %.4f\n", median(G$bias_cond), median(G$bias_marg)))
cat("\n== coverage del glmer per dgp: k_marg/k_cond/k_plim su n ==\n")
print(G[, .(n=sum(n), vs_marg=round(100*sum(k_marg)/sum(n),2), vs_cond=round(100*sum(k_cond)/sum(n),2),
            vs_plim=round(100*sum(k_plim)/sum(n),2), vs_margdgp=round(100*sum(k_margdgp)/sum(n),2)), by=dgp])
cat("\n== §3.14 RMSE: media per cella (usa sse = somma errori quadratici?) ==\n")
print(names(H))
H[, rmse_marg := sqrt(sse_marg/n)]
S <- H[method %in% c("mmm","gee_fg","gee_md","gee_kc","glmer"),
       .(rmse=mean(rmse_marg), ncell=.N), by=.(method,dgp)]
print(dcast(S, method~dgp, value.var="rmse"))
sp <- S[method %in% c("mmm","gee_fg","gee_md","gee_kc"), .(spread=100*(max(rmse)-min(rmse))), by=dgp]
cat("\nspread percentuale fra i quattro metodi marginali (media per cella della RMSE):\n"); print(sp)
cat("\nfattore glmer / minimo dei quattro:\n")
print(S[, .(fattore=.SD[method=="glmer",rmse]/min(.SD[method!="glmer",rmse])), by=dgp])

cat("\n\n== M-33 / D2-31r: gradi di liberta ==\n")
DF <- fread(path.expand("~/MMM/rev2/DFREF_cells.csv"))
cat("colonne:", paste(names(DF),collapse=","), "\n")
print(unique(DF[N==10 & ntime==4 & model==1, .(N,ntime,model,df_bw,df_k2,df_win,df_bwin)]))
print(unique(DF[N %in% c(10,45) & ntime %in% c(2,4,10) & model==1, .(N,ntime,df_bw,df_k2,df_win,df_bwin)])[order(N,ntime)])
WC <- DF[, .(worst_norm=max(abs(100*k_norm/n-5)), worst_bw=max(abs(100*k_bw/n-5)),
             worst_k2=max(abs(100*k_k2/n-5)), worst_win=max(abs(100*k_win/n-5)),
             worst_bwin=max(abs(100*k_bwin/n-5))), by=.(method,dgp)]
cat("\nworst case per metodo x meccanismo sotto i cinque riferimenti:\n")
print(WC[method %in% c("mmm","gee_fg","gee_md")][order(dgp,method)])
cat("\nprimato per riferimento:\n")
for (dg in unique(WC$dgp)) for (cc in c("worst_norm","worst_bw","worst_k2","worst_win","worst_bwin")) {
  x <- WC[dgp==dg]; cat(sprintf("%-12s %-11s -> %-7s %.3f\n", dg, cc, x[which.min(get(cc))]$method, min(x[[cc]])))
}
cat("\n== C. R2-01: identita algebrica del divario appaiato ==\n")
P <- fread(path.expand("~/MMM/rev2/PAIRED_mmm_fg_cells.csv"))
cat("colonne:", paste(names(P),collapse=","), "\n")
pp <- P[, .(npair=sum(npair), b=sum(b), cc=sum(cc), km=sum(km), kf=sum(kf)), by=dgp]
pp[, `:=`(gap_paired=100*(km-kf)/npair, gap_bc=100*(b-cc)/npair,
          rate_mmm_pair=100*km/npair, rate_fg_pair=100*kf/npair)]
print(pp)
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
H0 <- R[label=="manual" & set=="own" & method %in% c("mmm","gee_fg")]
own <- H0[, .(n=sum(n), k=sum(k_rej), rate=100*sum(k_rej)/sum(n)), by=.(dgp,method)]
print(dcast(own, dgp~method, value.var=c("n","rate")))
o <- dcast(own, dgp~method, value.var="rate")
cat("\ndifferenza dei tassi marginali sugli insiemi PROPRI (mmm - fg):\n")
print(o[, .(dgp, diff_own=mmm-gee_fg)])
cat("\nconfronto con il divario appaiato:\n")
print(merge(pp[,.(dgp,gap_paired)], o[,.(dgp,diff_own=mmm-gee_fg)], by="dgp"))
