LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages(library(data.table))
D <- as.data.table(readRDS(path.expand("~/MMM/rev2/R2_arm_raw.rds")))
cat("righe:", nrow(D), " celle:", uniqueN(D[, .(N,ntime,mode)]), " metodi:", paste(sort(unique(D$method)), collapse=","), "\n")
cat("repliche per cella e modo:\n")
print(dcast(D[, .(reps=uniqueN(rep)), by=.(N,ntime,mode)], N+ntime~mode, value.var="reps"))
S <- D[usable==TRUE, .(n=.N, k=sum(rej)), by=.(method,N,ntime,mode)]
S[, `:=`(rate=100*k/n, mcse=100*sqrt((k/n)*(1-k/n)/n))]
W <- dcast(S, method+N+ntime~mode, value.var=c("n","k","rate","mcse"))
W[, `:=`(diff_pp=rate_hybrid-rate_asrun, mcse_diff=sqrt(mcse_asrun^2+mcse_hybrid^2))]
W[, z := diff_pp/mcse_diff]
cat("\n===== LIVELLO GenBinaryY M3: H0 come eseguito (sigma=c(1,0)) vs H0 con la dipendenza di H1 (sigma=c(1,0.01)) =====\n")
print(as.data.frame(W[order(method,ntime,N), .(method,N,ntime,n_asrun,k_asrun,rate_asrun=round(rate_asrun,2),
  n_hybrid,k_hybrid,rate_hybrid=round(rate_hybrid,2), diff_pp=round(diff_pp,2),
  mcse_diff=round(mcse_diff,2), z=round(z,2))]), row.names=FALSE)
P <- S[, .(n=sum(n), k=sum(k)), by=.(method,mode)]
P[, `:=`(rate=100*k/n, mcse=100*sqrt((k/n)*(1-k/n)/n))]
PW <- dcast(P, method~mode, value.var=c("n","k","rate","mcse"))
PW[, `:=`(diff_pp=round(rate_hybrid-rate_asrun,2), mcse_diff=round(sqrt(mcse_asrun^2+mcse_hybrid^2),2))]
PW[, z := round(diff_pp/mcse_diff,2)]
cat("\n===== pooled sulle 9 celle =====\n")
print(as.data.frame(PW[order(method), .(method,n_asrun,rate_asrun=round(rate_asrun,2),
  n_hybrid,rate_hybrid=round(rate_hybrid,2),diff_pp,mcse_diff,z)]), row.names=FALSE)
cat(sprintf("\n|diff| massima per cella %.2f pp ; MCSE mediana della differenza %.2f pp ; celle |z|>1.96: %d su %d\n",
  max(abs(W$diff_pp),na.rm=TRUE), median(W$mcse_diff,na.rm=TRUE), sum(abs(W$z)>1.96,na.rm=TRUE), sum(is.finite(W$z))))
cat(sprintf("pooled: |diff| massima %.2f pp ; MCSE tipica %.2f ; metodi con |z|>1.96: %s\n",
  max(abs(PW$diff_pp)), median(PW$mcse_diff), paste(PW[abs(z)>1.96]$method, collapse=", ")))
cat("\n===== utilizzabilita nelle due versioni (il difetto tocca anche il denominatore?) =====\n")
U <- D[, .(tentativi=uniqueN(rep), utilizzabili=sum(usable)), by=.(method,mode)]
print(as.data.frame(dcast(U, method~mode, value.var=c("tentativi","utilizzabili"))), row.names=FALSE)
fwrite(W, path.expand("~/MMM/rev2/R2_arm_level_bycell.csv"))
fwrite(PW, path.expand("~/MMM/rev2/R2_arm_level_pooled.csv"))
cat("\nDONE_ARMAGG\n")
