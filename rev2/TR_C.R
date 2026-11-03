
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(8)
OUT <- path.expand("~/MMM/rev2")
CELL <- fread(file.path(OUT,"TR_cells.csv"))
DG <- c("GenBinaryY","GLMM-native","GLMM-AR1")
CELL[, dgp := fifelse(dgp=="GLMMnative","GLMM-native", fifelse(dgp=="GLMMAR1","GLMM-AR1", dgp))]

cat("###### 3. GRADI DI LIBERTA ######\n")
cat("--- gdl a N=10, T=4 per scenario (dai miei calcoli) ---\n")
print(unique(CELL[N==10 & ntime==4, .(model, nobs_exp, df_bw, df_k2, df_win, df_bwin, df_dec, df_decm)])[order(model)])
cat("\n--- confronto con DFREF_cells.csv ---\n")
DF <- fread(file.path(OUT,"DFREF_cells.csv"))
DF[, dgp := fifelse(dgp=="GLMMnative","GLMM-native", fifelse(dgp=="GLMMAR1","GLMM-AR1", dgp))]
NP <- c(4,5,4,4,5)
DF[, nobs_exp := N*(0.7*ntime + 0.3*ntime/2)]
DF[, `:=`(f_bw=N-NP[model], f_win=nobs_exp-NP[model], f_bwin=nobs_exp-N-NP[model]+1)]
cat("df_bw == N-p :", all(DF$df_bw==DF$f_bw), "\n")
cat("df_win == nobs_exp-p :", all(DF$df_win==DF$f_win), " (max scarto ", max(abs(DF$df_win-DF$f_win)), ")\n", sep="")
cat("df_bwin == nobs_exp-N-p+1 :", all(DF$df_bwin==DF$f_bwin), " (max scarto ", max(abs(DF$df_bwin-DF$f_bwin)), ")\n", sep="")
cat("righe DFREF:", nrow(DF), "  righe TR_cells:", nrow(CELL), "\n")
cat("a N=10,T=4,model=1 nel file: nobs_exp/df_bw/df_win/df_bwin =",
    unique(DF[N==10&ntime==4&model==1]$nobs_exp), unique(DF[N==10&ntime==4&model==1]$df_bw),
    unique(DF[N==10&ntime==4&model==1]$df_win), unique(DF[N==10&ntime==4&model==1]$df_bwin), "\n")
cat("DEC-DF (nobs_exp - N - 2) a N=10,T=4:", 34-10-2, "\n")

cat("\n--- concordanza dei conteggi con DFREF su chiave comune ---\n")
K <- merge(CELL[, .(dgp,method,N,ntime,model,n,k_norm,k_bw,k_k2,k_win,k_bwin)],
           DF[, .(dgp,method,N,ntime,model,nD=n,kn=k_norm,kb=k_bw,kk=k_k2,kw=k_win,kbw=k_bwin)],
           by=c("dgp","method","N","ntime","model"))
cat("righe unite:", nrow(K), " n uguali:", sum(K$n==K$nD), " k_norm uguali:", sum(K$k_norm==K$kn),
    " k_bw uguali:", sum(K$k_bw==K$kb, na.rm=TRUE), " k_win uguali:", sum(K$k_win==K$kw, na.rm=TRUE), "\n")

WC <- function(refcol, byN=TRUE, weighted=TRUE) {
  if (byN) {
    if (weighted) {
      x <- CELL[, .(n=sum(n), k=sum(get(refcol))), by=.(dgp,method,N)][, rate := 100*k/n]
    } else {
      x <- CELL[, .(rate=mean(100*get(refcol)/n)), by=.(dgp,method,N)]
    }
  } else {
    x <- CELL[, .(rate=100*get(refcol)/n), by=.(dgp,method,N,ntime,model)]
  }
  x[, .(worst=max(abs(rate-5), na.rm=TRUE)), by=.(dgp,method)]
}
refs <- c(norm="k_norm", t_bw="k_bw", t_k2="k_k2", t_win="k_win", t_bwin="k_bwin",
          t_dec="k_dec", t_decm="k_decm")
cat("\n###### worst case AGGREGATO PER N, PESATO PER REPLICHE (16 valori, max) ######\n")
for (nm in names(refs)) {
  w <- WC(refs[[nm]], byN=TRUE, weighted=TRUE)
  cat("\n[", nm, "] primato per meccanismo:\n", sep="")
  for (g in DG) { x <- w[dgp==g]; cat(sprintf("  %-12s best=%-8s %.4f | mmm %.4f | gee_fg %.4f | gee_md %.4f\n",
      g, x[which.min(worst)]$method, min(x$worst), x[method=="mmm"]$worst,
      x[method=="gee_fg"]$worst, x[method=="gee_md"]$worst)) }
}
cat("\n###### worst case AGGREGATO PER N, MEDIA NON PESATA DELLE CELLE ######\n")
for (nm in names(refs)) {
  w <- WC(refs[[nm]], byN=TRUE, weighted=FALSE)
  cat("\n[", nm, "] :\n", sep="")
  for (g in DG) { x <- w[dgp==g]; cat(sprintf("  %-12s best=%-8s %.4f | mmm %.4f | gee_fg %.4f | gee_md %.4f\n",
      g, x[which.min(worst)]$method, min(x$worst), x[method=="mmm"]$worst,
      x[method=="gee_fg"]$worst, x[method=="gee_md"]$worst)) }
}
cat("\n###### worst case PER CELLA (max sui 1440 punti) ######\n")
for (nm in names(refs)) {
  w <- WC(refs[[nm]], byN=FALSE)
  cat("\n[", nm, "] :\n", sep="")
  for (g in DG) { x <- w[dgp==g]; cat(sprintf("  %-12s best=%-8s %.4f | mmm %.4f | gee_fg %.4f | gee_md %.4f\n",
      g, x[which.min(worst)]$method, min(x$worst), x[method=="mmm"]$worst,
      x[method=="gee_fg"]$worst, x[method=="gee_md"]$worst)) }
}
cat("\n###### tabella completa 8 metodi, pesata per N, sotto i riferimenti chiave ######\n")
for (nm in c("norm","t_bw","t_k2","t_dec")) {
  w <- WC(refs[[nm]], byN=TRUE, weighted=TRUE)
  cat("\n[", nm, "]\n", sep=""); print(dcast(w, method ~ dgp, value.var="worst")[, lapply(.SD, function(z) if(is.numeric(z)) round(z,4) else z)])
}
cat("\n###### tassi pooled e |livello-5| sotto t per N>=45 (ramo 4) ######\n")
for (nm in c("norm","t_bw","t_dec")) {
  x <- CELL[N>=45, .(n=sum(n), k=sum(get(refs[[nm]]))), by=.(method)][, rate := 100*k/n][, dev := abs(rate-5)]
  cat("\n[", nm, "] N>=45 pooled su tutti i meccanismi:\n", sep=""); print(x[order(dev)][, .(method,n,rate=round(rate,4),dev=round(dev,4))])
}
cat("\nDONE_C\n")
