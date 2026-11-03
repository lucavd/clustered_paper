
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
H <- fread(path.expand("~/MMM/rev2/REV_H1_targets.csv")); G <- H[method=="glmer"]
cat("=== T-3.13-04 (27,2%) vs M-14 (23,8%): quale quantita da 0,0799/0,1048 ? ===\n")
cat(sprintf("A |media bias| per cella          : marg %.5f cond %.5f -> %.4f\n",
  mean(abs(G$bias_marg)), mean(abs(G$bias_cond)), 1-mean(abs(G$bias_cond))/mean(abs(G$bias_marg))))
cat(sprintf("B mediana |bias| per cella        : marg %.5f cond %.5f -> %.4f\n",
  median(abs(G$bias_marg)), median(abs(G$bias_cond)), 1-median(abs(G$bias_cond))/median(abs(G$bias_marg))))
cat(sprintf("C |media bias| pesata per repliche: marg %.5f cond %.5f -> %.4f\n",
  abs(sum(G$bias_marg*G$n)/sum(G$n)), abs(sum(G$bias_cond*G$n)/sum(G$n)),
  1-abs(sum(G$bias_cond*G$n)/sum(G$n))/abs(sum(G$bias_marg*G$n)/sum(G$n))))
cat(sprintf("D |media bias| per cella, M1-M2   : marg %.5f cond %.5f -> %.4f\n",
  mean(abs(G[model<=2]$bias_marg)), mean(abs(G[model<=2]$bias_cond)),
  1-mean(abs(G[model<=2]$bias_cond))/mean(abs(G[model<=2]$bias_marg))))
cat(sprintf("E |media bias| per cella, T>=4    : marg %.5f cond %.5f -> %.4f\n",
  mean(abs(G[ntime>=4]$bias_marg)), mean(abs(G[ntime>=4]$bias_cond)),
  1-mean(abs(G[ntime>=4]$bias_cond))/mean(abs(G[ntime>=4]$bias_marg))))
cat(sprintf("F |media bias| per cella, escl M4  : marg %.5f cond %.5f -> %.4f\n",
  mean(abs(G[model!=4]$bias_marg)), mean(abs(G[model!=4]$bias_cond)),
  1-mean(abs(G[model!=4]$bias_cond))/mean(abs(G[model!=4]$bias_marg))))
cat("\n=== T-3.13-08: la coverage pooled contro il marginale ===\n")
cat(sprintf("pooled vs_marg = %.4f%% ; vs_margdgp = %.4f%% ; vs_cond = %.4f%% (n=%d)\n",
  100*sum(G$k_marg)/sum(G$n), 100*sum(G$k_margdgp)/sum(G$n), 100*sum(G$k_cond)/sum(G$n), sum(G$n)))
cat("\n=== T-3.13-05 MAPE ===\n")
cat(sprintf("MAPE marg %.2f%%  cond %.2f%%  riduzione %.4f (spe/n medio per cella)\n",
  100*mean(G$spe_marg/G$n)/100, 100*mean(G$spe_cond/G$n)/100, 1-mean(G$spe_cond/G$n)/mean(G$spe_marg/G$n)))
cat(sprintf("MAPE pesata: marg %.2f  cond %.2f  riduzione %.4f\n",
  sum(G$spe_marg)/sum(G$n), sum(G$spe_cond)/sum(G$n), 1-sum(G$spe_cond)/sum(G$spe_marg)))
cat(sprintf("escludendo M4: riduzione %.4f\n",
  1-sum(G[model!=4]$spe_cond)/sum(G[model!=4]$spe_marg)))
cat("\n=== M-22 / T-3.9-09 / R2-13: livello del LRT, quale sottoinsieme ===\n")
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
L <- R[method=="glmer_lrt" & label=="manual"]
f <- function(d, tag, setv="own") {
  x <- d[set==setv]
  o <- x[, .(n=sum(n), rate=100*sum(k_rej)/sum(n)), by=dgp][order(dgp)]
  cat(sprintf("%-28s [%s] %s\n", tag, setv,
    paste(sprintf("%s %.2f (n=%d)", o$dgp, o$rate, o$n), collapse=" | ")))
}
f(L, "tutta la griglia")
f(L[N<=20], "N<=20")
f(L[N<=25], "N<=25")
f(L, "tutta la griglia", "intersection")
f(L[N<=20], "N<=20", "intersection")
cat("\n=== D2-09: vantaggio del LRT a N=10, due denominatori ===\n")
P10 <- R[N==10 & label=="manual" & set=="own" & method %in% c("glmer_lrt","gee_fg","gee_kc","gee_md")]
w <- P10[, .(n=sum(n), n_att=sum(n_attempt), k=sum(k_rej)), by=.(dgp,method)]
print(dcast(w, dgp~method, value.var=c("n","n_att")))
cat("\nDONE4\n")
