.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
# FINALCHK.R -- rilettura di ogni valore consegnato dai file che lo producono.
suppressMessages(library(data.table))
OUT <- path.expand("~/MMM/rev2")
sink(file.path(OUT, "FINALCHK_out.txt"), split = TRUE)

cat("== A. MMMROB: livello per N, filtri separati (file IMPL_MMMROB_byN.csv) ==\n")
byN <- fread(file.path(OUT, "IMPL_MMMROB_byN.csv"))
print(as.data.frame(byN[N==10, .(dgp, N, n_mod, lev_mod = round(lev_mod,2),
      n_rob, lev_rob = round(lev_rob,2))]), row.names = FALSE)
cat("worst case per meccanismo (max su N):\n")
print(as.data.frame(byN[, .(worst_mod = round(max(abs(lev_mod-5)),2),
                            worst_rob = round(max(abs(lev_rob-5)),2)), by = dgp]), row.names = FALSE)
cat("intervalli: mod", sprintf("%.2f-%.2f", min(byN[, max(abs(lev_mod-5)), by=dgp]$V1),
    max(byN[, max(abs(lev_mod-5)), by=dgp]$V1)), "| rob",
    sprintf("%.2f-%.2f", min(byN[, max(abs(lev_rob-5)), by=dgp]$V1),
    max(byN[, max(abs(lev_rob-5)), by=dgp]$V1)), "\n")

cat("\n== B. worst case MMM pubblicato vs ricalcolato (IMPL_MMMROB_worstcase.csv) ==\n")
print(as.data.frame(fread(file.path(OUT, "IMPL_MMMROB_worstcase.csv"))), row.names = FALSE)

cat("\n== C. rapporto SE_rob/SE_mod a N=10 (IMPL_MMMROB_raw.csv) ==\n")
Mn <- fread(file.path(OUT, "IMPL_MMMROB_raw.csv"))
U <- Mn[is.finite(est) & abs(est)<10 & is.finite(se_mod) & se_mod>0 & se_mod<10 &
        is.finite(se_rob) & se_rob>0 & se_rob<10]
cat(sprintf("N=10: rapporto medio %.4f (mediano %.4f), n=%d -> robusta %.1f%% piu' piccola\n",
    mean(U[N==10]$se_rob/U[N==10]$se_mod), median(U[N==10]$se_rob/U[N==10]$se_mod),
    nrow(U[N==10]), 100*(1-mean(U[N==10]$se_rob/U[N==10]$se_mod))))
cat(sprintf("fit totali nel file: %d ; utilizzabili con entrambe: %d\n", nrow(Mn), nrow(U)))

cat("\n== D. DGRID: livello e MCSE, quattro metodi richiesti (IMPL_DGRID_publishable.csv) ==\n")
A <- fread(file.path(OUT, "IMPL_DGRID_publishable.csv"))
for (mm in c("mmm","gee_fg","gee_md","gee_kc")) {
  z <- A[method==mm][order(Dvar)]
  cat(sprintf("%-9s : %s\n", mm, paste(sprintf("%.2f (%.3f)", z$level, z$mcse), collapse = " / ")))
}
cat("\n|livello-5| per metodo, ordine ICC 0.071/0.233/0.549:\n")
for (mm in sort(unique(A$method))) {
  z <- A[method==mm][order(Dvar)]
  cat(sprintf("%-9s : %s\n", mm, paste(sprintf("%.2f", z$dev), collapse = " / ")))
}
W <- dcast(A, method ~ Dvar, value.var = "rk")
cat("\nSpearman fra ordinamenti (ranghi su |livello-5|):\n")
cat(sprintf("  0.071 vs 0.233 = %.4f | 0.071 vs 0.549 = %.4f | 0.233 vs 0.549 = %.4f\n",
    cor(W[[2]],W[[3]],method="spearman"), cor(W[[2]],W[[4]],method="spearman"),
    cor(W[[3]],W[[4]],method="spearman")))
cat("n_use per D-metodo, intervallo:", range(A$n_use), "\n")

cat("\n== E. test di deterioramento (IMPL_DGRID_tests.csv) ==\n")
print(as.data.frame(fread(file.path(OUT, "IMPL_DGRID_tests.csv"))[, .(method,
      lev_low = round(lev_low,2), lev_high = round(lev_high,2), delta = round(delta,2),
      p_2camp = signif(p_2camp,3), p_trend = signif(p_trend,3))]), row.names = FALSE)

cat("\n== F. correlazione di rango fra MECCANISMI, griglia intera (REV_cells.csv) ==\n")
R <- fread(file.path(OUT, "REV_cells.csv"))
H <- R[label=="manual" & set=="own"]
z <- H[, .(n = sum(n), k = sum(k_rej)), by = .(dgp, method)][, dev := abs(100*k/n - 5)]
z[, rk := frank(dev, ties.method = "first"), by = dgp]
Q <- dcast(z, method ~ dgp, value.var = "rk")
print(as.data.frame(Q), row.names = FALSE)
nmv <- setdiff(names(Q), "method"); sp <- c()
for (a in seq_along(nmv)) for (b in seq_along(nmv)) if (a<b) {
  v <- cor(Q[[nmv[a]]], Q[[nmv[b]]], method="spearman"); sp <- c(sp, v)
  cat(sprintf("  %s vs %s = %.4f\n", nmv[a], nmv[b], v)) }
cat(sprintf("intervallo fra meccanismi: %.4f - %.4f\n", min(sp), max(sp)))

cat("\n== G. M5: diagonale della matrice D dei due generatori a effetti casuali ==\n")
gn <- readLines(path.expand("~/MMM/.rundir/manual/generation_model5.R"))
ga <- readLines(path.expand("~/MMM/.rundir/ar1/generation_models_ar1.R"))
cat("GLMM-native M5:\n"); cat(paste(grep("D <- matrix", gn, value=TRUE), collapse="\n"), "\n")
i <- grep("D <- matrix", gn)[1]; cat(paste(gn[i:(i+3)], collapse="\n"), "\n")
j <- grep("generate_data_model5", ga)[1]; cat("GLMM-AR1 M5:\n"); cat(paste(ga[j:(j+5)], collapse="\n"), "\n")

cat("\n== H. mediana tentativi per cella-metodo, file vecchio e nuovo ==\n")
for (fn in c("DGRID_raw.csv","IMPL_DGRID_raw.csv")) {
  X <- fread(file.path(OUT, fn))
  X[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  X[, usable := ifelse(method=="glmer_lrt", converged & is.finite(p),
      converged & is.finite(estimate) & is.finite(se) & is.finite(p) & abs(estimate)<10 & se>0 & se<10)]
  a <- X[, .(rows=.N, use=sum(usable)), by=.(Dvar,Ncl,Tt,method)]
  cat(sprintf("%-20s righe totali %d | righe/cella-metodo: mediana %.1f media %.1f | utilizzabili: mediana %.1f media %.1f\n",
      fn, nrow(X), median(a$rows), mean(a$rows), median(a$use), mean(a$use)))
}
cat("\n== I. calibrazione M-30 sotto i due generatori pseudocasuali ==\n")
for (tag in c("IMPL_MMMROB_n60_mmlb.csv","IMPL_MMMROB_n60_manual.csv")) {
  V <- fread(file.path(OUT, tag))
  Vu <- V[is.finite(est) & abs(est)<10 & is.finite(se_mod) & se_mod>0 & se_mod<10]
  cat(sprintf("%-28s (L'Ecuyer) n=%d | SE media %.5f | SD %.5f | rapporto %.4f\n",
      tag, nrow(Vu), mean(Vu$se_mod), sd(Vu$est), mean(Vu$se_mod)/sd(Vu$est)))
}
for (tag in c("VER_mmmrob_n60_mmlb.csv","VER_mmmrob_n60_manual.csv")) {
  if (!file.exists(file.path(OUT, tag))) next
  V <- fread(file.path(OUT, tag))
  Vu <- V[is.finite(est) & abs(est)<10 & is.finite(se_mod) & se_mod>0 & se_mod<10]
  cat(sprintf("%-28s (Mersenne) n=%d | SE media %.5f | SD %.5f | rapporto %.4f\n",
      tag, nrow(Vu), mean(Vu$se_mod), sd(Vu$est), mean(Vu$se_mod)/sd(Vu$est)))
}
cat("\nDONE_FINALCHK\n")
sink()
