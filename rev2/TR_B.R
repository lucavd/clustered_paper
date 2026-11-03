
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(8)
OUT <- path.expand("~/MMM/rev2")
CELL <- fread(file.path(OUT,"TR_cells.csv")); PD <- fread(file.path(OUT,"TR_paired_defs.csv"))
OWN <- fread(file.path(OUT,"TR_own.csv"))
DG <- c("GenBinaryY","GLMM-native","GLMM-AR1")
fix <- function(x) fifelse(x=="GLMMnative","GLMM-native", fifelse(x=="GLMMAR1","GLMM-AR1", x))
for (D in list(CELL,PD,OWN)) D[, dgp := fix(dgp)]
ord <- function(D, col) sapply(DG, function(g) D[dgp==g][[col]][1])

cat("###### 1. DIVARIO APPAIATO ######\n")
agg <- PD[, .(npair=sum(npair), b=sum(b), cc=sum(cc), km=sum(km), kf=sum(kf)), by=.(def,dgp)]
agg[, `:=`(gap_bc = 100*(b-cc)/npair, gap_k = 100*(km-kf)/npair,
           se_paired = 100*sqrt(b+cc-(b-cc)^2/npair)/npair,
           chisq = (abs(b-cc)-1)^2/(b+cc))]
agg[, p := pchisq(chisq,1,lower.tail=FALSE)]
agg[, z_paired := gap_bc/se_paired]
agg[, chisq_nocc := (b-cc)^2/(b+cc)]
setkey(agg, def, dgp)
for (dd in unique(agg$def)) {
  cat("\n-- def:", dd, "\n")
  print(agg[def==dd][match(DG,dgp), .(dgp,npair,b,cc,km,kf,gap_bc=round(gap_bc,4),
       se_paired=round(se_paired,4), z=round(z_paired,3), chisq=round(chisq,2),
       chisq_nocc=round(chisq_nocc,2), p=signif(p,3))])
}

cat("\n###### identita b-cc == km-kf per CELLA (def own_norm e inter_norm) ######\n")
for (dd in c("own_norm","inter_norm")) {
  X <- PD[def==dd]
  cat(dd, ": celle=", nrow(X), " max|(b-cc)-(km-kf)| = ", max(abs((X$b-X$cc)-(X$km-X$kf))),
      "  n. celle con scarto !=0: ", sum((X$b-X$cc)!=(X$km-X$kf)), "\n", sep="")
}

cat("\n###### SE NON APPAIATA e confronto ######\n")
ow <- OWN[, .(n_mmm=sum(n_mmm), k_mmm=sum(k_mmm), n_fg=sum(n_fg), k_fg=sum(k_fg)), by=dgp]
ow[, `:=`(r_mmm=100*k_mmm/n_mmm, r_fg=100*k_fg/n_fg)]
ow[, diff_own := r_mmm - r_fg]
ow[, se_unp := 100*sqrt((k_mmm/n_mmm)*(1-k_mmm/n_mmm)/n_mmm + (k_fg/n_fg)*(1-k_fg/n_fg)/n_fg)]
ow[, z_unp := diff_own/se_unp]
print(ow[match(DG,dgp), .(dgp,n_mmm,n_fg,r_mmm=round(r_mmm,4),r_fg=round(r_fg,4),
     diff_own=round(diff_own,4), se_unp=round(se_unp,4), z_unp=round(z_unp,3))])
M <- merge(agg[def=="own_norm", .(dgp,npair,gap_bc,se_paired,z_paired)], ow[,.(dgp,diff_own,se_unp,z_unp)], by="dgp")
M[, `:=`(rapporto_se = se_unp/se_paired, scarto_gap_vs_31 = gap_bc - diff_own)]
cat("\nappaiato vs non appaiato:\n")
print(M[match(DG,dgp), .(dgp, gap_appaiato=round(gap_bc,4), diff_31=round(diff_own,4),
   scarto_pp=round(scarto_gap_vs_31,4), se_app=round(se_paired,4), se_nonapp=round(se_unp,4),
   rapporto=round(rapporto_se,3), z_app=round(z_paired,2), z_nonapp=round(z_unp,2))])

cat("\n###### segno del divario sotto/sopra 25 cluster (def own_norm) ######\n")
for (thr in list(c("N<25","N < 25"), c("N<=25","N <= 25"))) NULL
X <- PD[def=="own_norm"]
sub <- X[N<25, .(npair=sum(npair), b=sum(b), cc=sum(cc)), by=dgp][, gap := 100*(b-cc)/npair]
sup <- X[N>=25, .(npair=sum(npair), b=sum(b), cc=sum(cc)), by=dgp][, gap := 100*(b-cc)/npair]
s25 <- X[N<=25, .(npair=sum(npair), b=sum(b), cc=sum(cc)), by=dgp][, gap := 100*(b-cc)/npair]
cat("N < 25 (solo N=10,15,20):\n"); print(sub[match(DG,dgp)])
cat("N <= 25:\n"); print(s25[match(DG,dgp)])
cat("N >= 25:\n"); print(sup[match(DG,dgp)])
cat("\ndivario per N (def own_norm):\n")
byN <- X[, .(npair=sum(npair), b=sum(b), cc=sum(cc)), by=.(dgp,N)][, gap := 100*(b-cc)/npair]
print(dcast(byN, N ~ dgp, value.var="gap")[, lapply(.SD, function(z) if(is.numeric(z)) round(z,3) else z)])
cat("\ndivario per T (def own_norm):\n")
byT <- X[, .(npair=sum(npair), b=sum(b), cc=sum(cc)), by=.(dgp,ntime)][, gap := 100*(b-cc)/npair]
print(dcast(byT, ntime ~ dgp, value.var="gap")[, lapply(.SD, function(z) if(is.numeric(z)) round(z,3) else z)])
cat("\nDONE_B1\n")
