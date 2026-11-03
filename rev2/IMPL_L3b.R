.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); OUT <- path.expand("~/MMM/rev2")
Cg <- fread(file.path(OUT,"IMPL_L3a_glmer_cells.csv"))
S  <- fread(file.path(OUT,"IMPL_L3a_cells.csv"))
H  <- fread(file.path(OUT,"REV_H1_targets.csv"))[method=="glmer"]

cat("=== A. |bias| a livello di MEDIA DI CELLA (|mean_est - t|) ===\n")
cat(sprintf("celle=%d  |bias|_marg=%.6f  |bias|_cond=%.6f  |bias|_emp=%.6f  riduzione=%.4f%%\n",
  nrow(Cg), mean(abs(Cg$b_marg)), mean(abs(Cg$b_cond)), mean(abs(Cg$b_emp)),
  100*(1-mean(abs(Cg$b_cond))/mean(abs(Cg$b_marg)))))
print(Cg[, .(celle=.N, ab_marg=mean(abs(b_marg)), ab_cond=mean(abs(b_cond)),
             rid=100*(1-mean(abs(b_cond))/mean(abs(b_marg)))), by=dgp][order(dgp)])
cat(sprintf("pesata per repliche: |bias|_marg=%.6f |bias|_cond=%.6f rid=%.4f%%\n",
  sum(abs(Cg$b_marg)*Cg$n)/sum(Cg$n), sum(abs(Cg$b_cond)*Cg$n)/sum(Cg$n),
  100*(1-sum(abs(Cg$b_cond)*Cg$n)/sum(abs(Cg$b_marg)*Cg$n))))

cat("\n=== B. T-3.13-06 con denominatori POOLED per meccanismo ===\n")
P <- S[method %in% c("glmer","gee_fg"), .(cov_marg=100*sum(k_marg)/sum(n),
        cov_cond=100*sum(k_cond)/sum(n)), by=.(dgp,method)]
W <- dcast(P, dgp~method, value.var=c("cov_marg","cov_cond"))
W[, `:=`(gain=cov_cond_glmer-cov_marg_glmer, gap_fg=cov_marg_gee_fg-cov_marg_glmer,
         gap_95=95-cov_marg_glmer)]
print(W[, .(dgp, cov_marg_glmer, cov_cond_glmer, cov_marg_gee_fg, gain, gap_fg, gap_95)])
cat(sprintf("medie sui tre meccanismi: gain=%.4f  gap_fg=%.4f  gap_95=%.4f\n",
  mean(W$gain), mean(W$gap_fg), mean(W$gap_95)))
cat(sprintf("quota gain/gap_fg = %.4f%% ; gain/gap_95 = %.4f%%\n",
  100*mean(W$gain)/mean(W$gap_fg), 100*mean(W$gain)/mean(W$gap_95)))

cat("\n=== C. T-3.13-08 pooled contro il marginale specifico per meccanismo ===\n")
cat(sprintf("REV_H1_targets: pooled k_marg=%.4f  k_margdgp=%.4f  k_cond=%.4f (n=%d)\n",
  100*sum(H$k_marg)/sum(H$n), 100*sum(H$k_margdgp)/sum(H$n), 100*sum(H$k_cond)/sum(H$n), sum(H$n)))
cat(sprintf("mio pooled dai chunk: marg=%.4f  cond=%.4f (n=%d)\n",
  100*sum(S[method=="glmer",k_marg])/sum(S[method=="glmer",n]),
  100*sum(S[method=="glmer",k_cond])/sum(S[method=="glmer",n]), sum(S[method=="glmer",n])))

cat("\n=== D. fattore GLMM/comparatori, RMSE media per cella ===\n")
R <- S[, .(rp=mean(sqrt(sq_marg/n))), by=.(dgp,method)]
for (dd in sort(unique(R$dgp))) {
  g <- R[dgp==dd & method=="glmer", rp]
  v <- R[dgp==dd & method %in% c("mmm","gee_fg","gee_md","gee_kc")]
  cat(sprintf("%-12s glmer=%.6f ; fattori: %s ; range %.4f-%.4f\n", dd, g,
    paste(sprintf("%s=%.4f", v$method, g/v$rp), collapse=" "), min(g/v$rp), max(g/v$rp)))
}
cat("\n=== E. T-COH-03: 91.91 confronto con tutti i pooled ===\n")
allp <- S[, .(cov_marg=100*sum(k_marg)/sum(n), cov_cond=100*sum(k_cond)/sum(n),
              cov_emp=100*sum(k_emp)/sum(n)), by=.(method,dgp)]
print(allp[abs(cov_marg-91.91)<0.3 | abs(cov_cond-91.91)<0.3 | abs(cov_emp-91.91)<0.3])
cat(sprintf("glmer pooled tre meccanismi: marg=%.4f cond=%.4f emp_cella=%.4f\n",
  100*sum(S[method=="glmer",k_marg])/sum(S[method=="glmer",n]),
  100*sum(S[method=="glmer",k_cond])/sum(S[method=="glmer",n]),
  100*sum(S[method=="glmer",k_emp])/sum(S[method=="glmer",n])))
cat(sprintf("REV_H1_targets glmer pooled plim = %.4f\n", 100*sum(H$k_plim)/sum(H$n)))
cat("\n=== F. MCSE copertura ===\n")
G <- S[method=="glmer"]
for (nm in c("k_marg","k_cond")) {
  p <- sum(G[[nm]])/sum(G$n); cat(sprintf("%s pooled p=%.6f mcse=%.6f pp\n", nm, 100*p,
    100*sqrt(p*(1-p)/sum(G$n))))
}
Ga <- G[dgp=="GLMM-AR1"]; p <- sum(Ga$k_cond)/sum(Ga$n)
cat(sprintf("AR1 cond p=%.5f mcse=%.6f pp ; |93.16313-93.19300|=%.5f pp = %.3f MCSE\n",
  100*p, 100*sqrt(p*(1-p)/sum(Ga$n)), abs(93.16313-93.19300),
  abs(93.16313-93.19300)/(100*sqrt(p*(1-p)/sum(Ga$n)))))
cat("FATTO_B\n")
