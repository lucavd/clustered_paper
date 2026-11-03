# R2_kc.R -- GIRO 2, contestazione 4: distribuzione COMPLETA del rapporto KC/non corretto,
# incluse le righe che il filtro scarta (se>=10, |stima|>=10).
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method %in% c("gee","gee_kc","gee_md","gee_fg")]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  # NESSUN filtro su se<10 o |stima|<10: prendo tutte le righe con SE finita e positiva
  y <- x[converged & is.finite(se) & se>0, .(dgp,model,N,ntime,label,sim_id,method,se,estimate)]
  if (!nrow(y)) return(NULL)
  w <- tryCatch(dcast(y, dgp+model+N+ntime+label+sim_id~method, value.var=c("se","estimate")),
                error=function(e) NULL)
  if (is.null(w)) return(NULL)
  for (nm in c("se_gee","se_gee_kc","se_gee_md","se_gee_fg","estimate_gee"))
    if (!nm %in% names(w)) w[[nm]] <- NA_real_
  setnames(w, c("se_gee","se_gee_kc","se_gee_md","se_gee_fg"), c("gee","gee_kc","gee_md","gee_fg"))
  w <- w[is.finite(gee) & gee>0 & is.finite(gee_kc)]
  if (!nrow(w)) return(NULL)
  w[, .(dgp,model,N,ntime,label,sim_id,gee,gee_kc,gee_md,gee_fg,estimate_gee)]
}
W <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)), fill=TRUE)
W[, `:=`(r_kc=gee_kc/gee, r_md=gee_md/gee, r_fg=gee_fg/gee)]
cat("righe con SE gee e gee_kc entrambe finite e positive, NESSUN filtro:", nrow(W), "\n")
cat(sprintf("di cui con se_gee>=10: %d ; se_kc>=10: %d ; |stima|>=10: %d\n",
  sum(W$gee>=10), sum(W$gee_kc>=10), sum(abs(W$estimate_gee)>=10, na.rm=TRUE)))
cat("\n===== DISTRIBUZIONE del rapporto KC / non corretto, su TUTTE le righe =====\n")
q <- quantile(W$r_kc, c(0.01,0.25,0.5,0.75,0.95,0.99,0.999,0.9999,1), na.rm=TRUE)
print(round(q, 6))
cat(sprintf("mediana %.4f | 99o pct %.4f | 99,99o pct %.4g | massimo %.6g\n",
  median(W$r_kc,na.rm=TRUE), q[["99%"]], q[["99.99%"]], max(W$r_kc,na.rm=TRUE)))
cat("\nconteggi sopra soglia e celle coinvolte:\n")
for (thr in c(10,1e2,1e4,1e6,1e9,1e12,1e14,1e15,1e18,1e21)) {
  s <- W[r_kc>thr]
  cat(sprintf("  > %-7.0e : righe %7d | celle %4d | punti di disegno %4d\n", thr, nrow(s),
    uniqueN(s[, paste(dgp,model,N,ntime,label)]), uniqueN(s[, paste(dgp,model,N,ntime)])))
}
cat("\n===== esiste QUALCHE cella con un rapporto di ordine 10^21? =====\n")
cat(sprintf("massimo assoluto su %d righe: %.6g (log10 = %.2f)\n", nrow(W), max(W$r_kc,na.rm=TRUE),
  log10(max(W$r_kc,na.rm=TRUE))))
top <- W[order(-r_kc)][1:8, .(dgp,model,N,ntime,label,sim_id, se_gee=round(gee,6),
   se_kc=round(gee_kc,6), r_kc=signif(r_kc,4), est=round(estimate_gee,3))]
print(as.data.frame(top), row.names=FALSE)
cat("\nmassimo per cella: le dieci celle peggiori\n")
cel <- W[, .(n=.N, max_r=max(r_kc,na.rm=TRUE), med_r=round(median(r_kc,na.rm=TRUE),4)),
         by=.(dgp,model,N,ntime,label)][order(-max_r)][1:10]
print(as.data.frame(cel[, .(dgp,model,N,ntime,label,n,max_r=signif(max_r,4),med_r)]), row.names=FALSE)
cat("\n===== confronto con MD e FG sugli stessi fit =====\n")
cat(sprintf("massimo: KC %.4g | MD %.4g | FG %.4g\n", max(W$r_kc,na.rm=TRUE),
  max(W$r_md,na.rm=TRUE), max(W$r_fg,na.rm=TRUE)))
cat(sprintf("99,99o pct: KC %.4g | MD %.4g | FG %.4g\n", quantile(W$r_kc,.9999,na.rm=TRUE),
  quantile(W$r_md,.9999,na.rm=TRUE), quantile(W$r_fg,.9999,na.rm=TRUE)))
cat(sprintf("righe sopra 1e12: KC %d | MD %d | FG %d\n", sum(W$r_kc>1e12,na.rm=TRUE),
  sum(W$r_md>1e12,na.rm=TRUE), sum(W$r_fg>1e12,na.rm=TRUE)))
cat("\n===== il rapporto si presenta anche sulle righe che il filtro SCARTA? =====\n")
W[, dropped := !(gee<10 & gee_kc<10 & abs(estimate_gee)<10)]
print(as.data.frame(W[, .(righe=.N, mediana=round(median(r_kc,na.rm=TRUE),4),
   p99=round(quantile(r_kc,.99,na.rm=TRUE),4), massimo=signif(max(r_kc,na.rm=TRUE),4)), by=dropped]),
   row.names=FALSE)
fwrite(cel, path.expand("~/MMM/rev2/R2_kc_worstcells.csv"))
cat("\nDONE_KC\n")
