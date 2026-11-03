.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
setwd(path.expand("~/MMM"))
DG <- fread("rev2/T3_degeneracy.csv")
P <- 207; K <- 25
cat("=== degenerazione osservata nel replay col seme 20260807 vs ipergeometrica ===\n")
for (n in c(10,15,20,25,50)) {
  p0 <- dhyper(0, K, P-K, n)
  obs <- sum(DG[N==n]$lev_inv==1); m <- 500
  cat(sprintf("N=%2d  p0=%.6f  attesi=%.1f  osservati=%d (%.2f%%)  tetto atteso=%.4f%%  osservato=%.2f%%  z=%+.2f  P(<=oss)=%.4f\n",
      n, p0, m*p0, obs, 100*obs/m, 100*(1-p0), 100*(1-obs/m),
      (obs-m*p0)/sqrt(m*p0*(1-p0)), pbinom(obs, m, p0)))
}
cat("\n=== il 78.2% pubblicato e' compatibile col disegno? ===\n")
p0 <- dhyper(0,K,P-K,10)
cat(sprintf("tetto atteso a N=10 = %.4f%% (media su repliche di 500 estrazioni)\n", 100*(1-p0)))
cat(sprintf("78.2%% = 391/500 estrazioni identificabili; degeneri = 109\n"))
cat(sprintf("P(degeneri <= 109 | 500 estrazioni, p0=%.6f) = %.5f\n", p0, pbinom(109, 500, p0)))
cat(sprintf("P(identificabili >= 391) = %.5f  -> il 78.2%% e' una fluttuazione a %.2f sd, non un'impossibilita'\n",
            pbinom(109, 500, p0), (109-500*p0)/sqrt(500*p0*(1-p0))))
cat(sprintf("intervallo di predizione 95%% del tetto empirico su 500 estrazioni: %.2f%% - %.2f%%\n",
            100*(1-qbinom(0.975,500,p0)/500), 100*(1-qbinom(0.025,500,p0)/500)))
cat(sprintf("\n71.6%% (=358/500, file di arbitrato REV_subsample_all): degeneri=142, P(>=142)=%.5f\n",
            1-pbinom(141,500,p0)))
cat("\n=== quanti pazienti hanno righe invasive=1 e come sono distribuiti ===\n")
cat("K (pazienti con almeno una riga invasive=1) = 25; con TUTTE le righe invasive=1 = 19; misti = 6\n")
cat(sprintf("se il criterio fosse 'nessun paziente interamente invasivo' (K=19): p0(N=10)=%.6f, tetto=%.4f%%\n",
            dhyper(0,19,207-19,10), 100*(1-dhyper(0,19,207-19,10))))
cat("\nDONE_T8\n")
