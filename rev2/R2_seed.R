# R2_seed.R -- GIRO 2, contestazione 6: correlazione fra stime a N diverso con lo stesso
# (model, sim_id), e fra H0 e H1 nella stessa cella. Conseguenze su MCSE e bootstrap.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method %in% c("mmm","gee","gee_fg","glmer")]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10,
    .(dgp,method,model,N,ntime,label,sim_id,estimate,p)]
}
D <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)), fill=TRUE)
cat("righe:", nrow(D), "\n")

cat("\n############ A. stesso (model, sim_id), N DIVERSI: le stime sono correlate?\n")
pairs <- list(c(10,15), c(10,20), c(10,145), c(45,55), c(135,145), c(20,45))
out <- list(); k <- 0
for (mth in c("mmm","gee_fg","glmer")) for (pp in pairs) {
  w <- dcast(D[method==mth & label=="manual" & ntime==6 & N %in% pp,
               .(dgp,model,sim_id,N,estimate)],
             dgp+model+sim_id~N, value.var="estimate")
  setnames(w, as.character(pp), c("a","b"))
  w <- w[is.finite(a) & is.finite(b)]
  if (nrow(w) < 100) next
  k <- k+1
  out[[k]] <- data.table(method=mth, N1=pp[1], N2=pp[2], n=nrow(w),
    r_pearson=cor(w$a,w$b), r_spearman=cor(w$a,w$b,method="spearman"),
    p_val=cor.test(w$a,w$b)$p.value)
}
A <- rbindlist(out)
print(as.data.frame(A[, .(method,N1,N2,n,r_pearson=round(r_pearson,4),
  r_spearman=round(r_spearman,4), p=signif(p_val,3))]), row.names=FALSE)
cat(sprintf("\n|r| massimo fra N diversi: %.4f ; mediano %.4f\n", max(abs(A$r_pearson)), median(abs(A$r_pearson))))

cat("\n############ B. stessa cella, H0 vs H1: le stime sono correlate?\n")
out2 <- list(); k <- 0
for (mth in c("mmm","gee_fg","glmer")) for (nn in c(10,45,145)) for (T in c(2,6,10)) {
  w <- dcast(D[method==mth & N==nn & ntime==T, .(dgp,model,sim_id,label,estimate)],
             dgp+model+sim_id~label, value.var="estimate")
  w <- w[is.finite(manual) & is.finite(mmlb)]
  if (nrow(w) < 100) next
  k <- k+1
  out2[[k]] <- data.table(method=mth, N=nn, ntime=T, n=nrow(w),
    r=cor(w$manual,w$mmlb), p=cor.test(w$manual,w$mmlb)$p.value)
}
B <- rbindlist(out2)
print(as.data.frame(B[, .(method,N,ntime,n,r=round(r,4),p=signif(p,3))]), row.names=FALSE)
cat(sprintf("\n|r| H0-H1 massimo %.4f ; mediano %.4f ; celle con p<0.05: %d su %d\n",
  max(abs(B$r)), median(abs(B$r)), sum(B$p<0.05), nrow(B)))

cat("\n############ C. le decisioni di rifiuto sono correlate fra N diversi? (tocca il bootstrap)\n")
out3 <- list(); k <- 0
for (mth in c("mmm","gee_fg")) for (pp in list(c(10,15),c(10,20),c(135,145))) {
  w <- dcast(D[method==mth & label=="manual" & ntime==6 & N %in% pp,
               .(dgp,model,sim_id,N,rej=as.integer(p<0.05))],
             dgp+model+sim_id~N, value.var="rej")
  setnames(w, as.character(pp), c("a","b"))
  w <- w[is.finite(a) & is.finite(b)]
  if (nrow(w)<100) next
  k <- k+1
  out3[[k]] <- data.table(method=mth, N1=pp[1], N2=pp[2], n=nrow(w),
    phi=cor(w$a,w$b), p=suppressWarnings(cor.test(w$a,w$b)$p.value))
}
C <- rbindlist(out3)
print(as.data.frame(C[, .(method,N1,N2,n,phi=round(phi,4),p=signif(p,3))]), row.names=FALSE)

cat("\n############ D. conseguenza sul confronto fra MECCANISMI (disegni condivisi)\n")
out4 <- list(); k <- 0
for (mth in c("mmm","gee_fg","glmer")) for (nn in c(10,45)) for (T in c(6,10)) {
  w <- dcast(D[method==mth & label=="manual" & N==nn & ntime==T, .(model,sim_id,dgp,estimate)],
             model+sim_id~dgp, value.var="estimate")
  nmz <- intersect(c("GenBinaryY","GLMM-native","GLMM-AR1"), names(w))
  if (length(nmz)<3) next
  w <- w[is.finite(get("GenBinaryY")) & is.finite(get("GLMM-native")) & is.finite(get("GLMM-AR1"))]
  if (nrow(w)<100) next
  k <- k+1
  out4[[k]] <- data.table(method=mth, N=nn, ntime=T, n=nrow(w),
    r_GB_native=cor(w$GenBinaryY, w$`GLMM-native`),
    r_native_AR1=cor(w$`GLMM-native`, w$`GLMM-AR1`),
    r_GB_AR1=cor(w$GenBinaryY, w$`GLMM-AR1`))
}
Dd <- rbindlist(out4)
print(as.data.frame(Dd[, .(method,N,ntime,n,r_GB_native=round(r_GB_native,4),
  r_native_AR1=round(r_native_AR1,4), r_GB_AR1=round(r_GB_AR1,4))]), row.names=FALSE)
cat(sprintf("\ncorrelazione native-AR1 fra le stime: mediana %.4f, massimo %.4f\n",
  median(Dd$r_native_AR1), max(Dd$r_native_AR1)))
cat("=> se non nulla, la differenza di 0,22 pp fra i due meccanismi e' appaiata e la sua\n")
cat("   MCSE non pesata sovrastima l'incertezza (il confronto e' piu' preciso di quanto dichiarato),\n")
cat("   mentre l'indipendenza fra meccanismi assunta in §2.2/§4.4 non regge.\n")
fwrite(A, path.expand("~/MMM/rev2/R2_seed_acrossN.csv"))
fwrite(Dd, path.expand("~/MMM/rev2/R2_seed_acrossDGP.csv"))
cat("\nDONE_SEED\n")
