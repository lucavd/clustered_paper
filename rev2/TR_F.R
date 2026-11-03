
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)}); setDTthreads(1)
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
DG <- c("GenBinaryY","GLMM-native","GLMM-AR1")
fx <- function(x) fifelse(x=="GLMMnative","GLMM-native", fifelse(x=="GLMMAR1","GLMM-AR1", x))
cut0 <- as.POSIXct("2026-08-08 10:18:26")

cat("###### A. ramo 1: media sui PUNTI DI DISEGNO di |tasso-5| (convenzione §2.4) ######\n")
R <- fread(file.path(OUT,"REV_cells.csv")); R[, dgp := fx(dgp)]
H0 <- R[label=="manual" & set=="own"]
H0[, rate := 100*k_rej/n]
for (tt in c(2,3)) {
  x <- H0[ntime==tt, .(mad_punti = mean(abs(rate-5)), n_punti=.N,
                       usable_pct = 100*sum(n_usable)/sum(n_attempt),
                       pooled_dev = abs(100*sum(k_rej)/sum(n)-5)), by=method]
  cat("\n  T =", tt, "  (media su punti di disegno, pooled sui 3 meccanismi)\n")
  print(x[order(mad_punti)][, .(method, n_punti, mad_punti=round(mad_punti,4),
        pooled_dev=round(pooled_dev,4), usable_pct=round(usable_pct,4))])
}
cat("\n###### B. mtime per indice di chunk ######\n")
fi <- as.data.table(file.info(list.files(RES, full.names=TRUE, pattern="[.]rds$")), keep.rownames="path")
fi[, bn := basename(path)][, cidx := as.integer(sub(".*_c([0-9]+)[.]rds$","\\1",bn))]
fi[, grp := sub("_c[0-9]+[.]rds$","",bn)][, post := mtime > cut0]
print(fi[, .(n=.N, post=sum(post)), by=cidx][order(cidx)])
cat("\ngruppi per numero di chunk post-cut:\n")
print(fi[, .(npost=sum(post)), by=grp][, .N, by=npost][order(npost)])

cat("\n###### C. tasso di discordanza mmm/FG su chunk PRE-cut contro POST-cut ######\n")
fs <- list.files(RES, pattern="[.]rds$"); key <- sub("_c[0-9]+[.]rds$","",fs)
H0g <- unique(key[grepl("_manual$", key)])
mt <- fi[, .(bn, post)]; setkey(mt, bn)
one <- function(g) {
  f <- sort(fs[key==g])
  L <- lapply(seq_along(f), function(i) {
    x <- try(as.data.table(readRDS(file.path(RES,f[i]))), silent=TRUE)
    if (inherits(x,"try-error")) return(NULL)
    x <- x[term=="time:trt" & method %in% c("mmm","gee_fg")]
    if (!nrow(x)) return(NULL)
    x[, `:=`(chunk=i, post=mt[f[i], post])]
    x[, .(method, sim_id, chunk, post, converged, estimate, se, p)]
  })
  d <- rbindlist(L, fill=TRUE); if (!nrow(d)) return(NULL)
  d[, converged := as.logical(converged)]; d[is.na(converged), converged := FALSE]
  d[, usable := converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
                abs(estimate)<10 & se>0 & se<10]
  d[, rid := paste0(chunk,"_",sim_id)]
  u <- d[usable==TRUE]
  W <- dcast(u[, .(rid, post, method, rej=p<0.05)], rid + post ~ method,
             value.var="rej", fun.aggregate=function(z) z[1])
  if (!all(c("mmm","gee_fg") %in% names(W))) return(NULL)
  W <- W[!is.na(mmm) & !is.na(gee_fg)]
  W[, .(npair=.N, b=sum(mmm & !gee_fg), cc=sum(!mmm & gee_fg),
        km=sum(mmm), kf=sum(gee_fg)), by=post][, grp := g]
}
Z <- rbindlist(mclapply(H0g, function(g) tryCatch(one(g), error=function(e) NULL), mc.cores=45), fill=TRUE)
Z[, dgp := fx(sub("_m[0-9]+_.*$","",grp))]
S <- Z[, .(npair=sum(npair), b=sum(b), cc=sum(cc), km=sum(km), kf=sum(kf)), by=.(dgp,post)]
S[, `:=`(disc_pct = 100*(b+cc)/npair, gap_pp = 100*(b-cc)/npair,
         r_mmm=100*km/npair, r_fg=100*kf/npair)]
print(S[order(dgp,post)][, .(dgp, post, npair, b, cc, disc_pct=round(disc_pct,4),
      gap_pp=round(gap_pp,4), r_mmm=round(r_mmm,4), r_fg=round(r_fg,4))])
cat("\nfile pubblicato: disc_pct e gap per meccanismo\n")
F <- data.table(dgp=DG, npair=c(430139,434919,433219), b=c(9825,8415,9056), cc=c(7811,6548,6371))
F[, `:=`(disc_pct=100*(b+cc)/npair, gap_pp=100*(b-cc)/npair)]
print(F[, .(dgp, npair, b, cc, disc_pct=round(disc_pct,4), gap_pp=round(gap_pp,4))])
fwrite(S, file.path(OUT,"TR_prepost.csv"))
cat("\nDONE_F\n")
