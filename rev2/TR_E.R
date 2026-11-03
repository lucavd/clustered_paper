
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table)}); setDTthreads(8)
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
DG <- c("GenBinaryY","GLMM-native","GLMM-AR1")
fx <- function(x) fifelse(x=="GLMMnative","GLMM-native", fifelse(x=="GLMMAR1","GLMM-AR1", x))

cat("###### FORENSICA sui mtime dei chunk ######\n")
fi <- as.data.table(file.info(list.files(RES, full.names=TRUE, pattern="[.]rds$")), keep.rownames="path")
fi[, bn := basename(path)]
fi[, pref := sub("_m[0-9]+_.*$","",bn)]
fi[, arm := fifelse(grepl("_manual_", bn), "manual", "mmlb")]
cut0 <- as.POSIXct("2026-08-08 10:18:26")
fi[, post := mtime > cut0]
print(fi[, .(n=.N, post_cut=sum(post), frac=round(mean(post),4)), by=.(pref,arm)][order(pref,arm)])
cat("\ndistribuzione delle date (giorno):\n")
print(fi[, .N, by=.(giorno=as.Date(mtime))][order(giorno)])

cat("\n###### b+c sotto ogni definizione, contro il file (GenBinaryY b+c=17636) ######\n")
PD <- fread(file.path(OUT,"TR_paired_defs.csv")); PD[, dgp := fx(dgp)]
ag <- PD[, .(npair=sum(npair), b=sum(b), cc=sum(cc)), by=.(def,dgp)][, bc := b+cc]
FILEV <- data.table(dgp=DG, npair_f=c(430139,434919,433219), b_f=c(9825,8415,9056), c_f=c(7811,6548,6371))
FILEV[, bc_f := b_f+c_f]
M <- merge(ag, FILEV, by="dgp")
M[, `:=`(d_npair=npair-npair_f, d_b=b-b_f, d_c=cc-c_f, d_bc=bc-bc_f)]
print(dcast(M, def ~ dgp, value.var="d_bc"))
cat("\nscarto su npair:\n"); print(dcast(M, def ~ dgp, value.var="d_npair"))

cat("\n###### D2-34: i cinque rami ######\n")
R <- fread(file.path(OUT,"REV_cells.csv")); R[, dgp := fx(dgp)]
H0 <- R[label=="manual" & set=="own"]
cat("-- ramo 1: T<=3, livello e usabilita per metodo (pooled sui 3 meccanismi) --\n")
for (tt in c(2,3)) {
  x <- H0[ntime==tt, .(n=sum(n), k=sum(k_rej), att=sum(n_attempt), us=sum(n_usable)), by=method]
  x[, `:=`(rate=100*k/n, dev=abs(100*k/n-5), usable_pct=100*us/att)]
  cat("\n  T =", tt, "\n"); print(x[order(dev)][, .(method, n, rate=round(rate,4), dev=round(dev,4),
      usable_pct=round(usable_pct,4))])
}
cat("\n-- ramo 3: N<25 & T>=4: usabilita, validita non condizionata, livello --\n")
x <- H0[N<25 & ntime>=4, .(n=sum(n), k=sum(k_rej), att=sum(n_attempt), us=sum(n_usable)), by=method]
x[, `:=`(rate=100*k/n, usable_pct=100*us/att)]
x[, unc_val := usable_pct/100*(1-k/n)*100]
print(x[order(-usable_pct)][, .(method, n, att, usable_pct=round(usable_pct,4),
    unc_val=round(unc_val,4), livello=round(rate,4), volte_nominale=round(rate/5,4))])
cat("\n  (variante N<=25):\n")
x2 <- H0[N<=25 & ntime>=4, .(n=sum(n), k=sum(k_rej), att=sum(n_attempt), us=sum(n_usable)), by=method]
x2[, `:=`(rate=100*k/n, usable_pct=100*us/att)]; x2[, unc_val := usable_pct/100*(1-k/n)*100]
print(x2[order(-usable_pct)][, .(method, usable_pct=round(usable_pct,4), unc_val=round(unc_val,4),
    livello=round(rate,4), volte_nominale=round(rate/5,4))])
cat("\n-- ramo 2: pendenza casuale (model 3,4,5) & T=2: usabilita glmer contro glmer_lrt e mmm --\n")
x <- H0[model %in% c(3,4,5) & ntime==2, .(att=sum(n_attempt), us=sum(n_usable), n=sum(n), k=sum(k_rej)), by=method]
x[, `:=`(usable_pct=100*us/att, rate=100*k/n)]
print(x[order(-usable_pct)][, .(method, att, us, usable_pct=round(usable_pct,4), n, rate=round(rate,4))])
cat("\n  solo intercetta casuale (model 1,2) & T=2:\n")
x <- H0[model %in% c(1,2) & ntime==2, .(att=sum(n_attempt), us=sum(n_usable), n=sum(n), k=sum(k_rej)), by=method]
x[, `:=`(usable_pct=100*us/att, rate=100*k/n)]
print(x[order(-usable_pct)][, .(method, att, us, usable_pct=round(usable_pct,4), rate=round(rate,4))])
cat("\n-- ramo 5: rango medio combinato su livello, copertura, usabilita (tutta la griglia) --\n")
H1 <- R[label=="mmlb" & set=="own"]
lev <- H0[, .(n=sum(n), k=sum(k_rej), att=sum(n_attempt), us=sum(n_usable)), by=method]
lev[, `:=`(dev=abs(100*k/n-5), usable_pct=100*us/att)]
cov <- H1[, .(ncov=sum(n_cov), kcov=sum(k_cov)), by=method][, covdev := abs(100*kcov/ncov - 95)]
Z <- merge(lev[, .(method, dev, usable_pct)], cov[, .(method, covdev, cov=100*kcov/ncov)], by="method", all.x=TRUE)
Z[, `:=`(r_lev=rank(dev), r_cov=rank(covdev, na.last="keep"), r_us=rank(-usable_pct))]
Z[, rango_medio := rowMeans(cbind(r_lev,r_cov,r_us), na.rm=TRUE)]
print(Z[order(rango_medio)][, .(method, dev=round(dev,4), cov=round(cov,4), covdev=round(covdev,4),
    usable_pct=round(usable_pct,4), r_lev, r_cov, r_us, rango_medio=round(rango_medio,4))])
cat("\n  ristretto a N<=25:\n")
lev2 <- H0[N<=25, .(n=sum(n), k=sum(k_rej), att=sum(n_attempt), us=sum(n_usable)), by=method]
lev2[, `:=`(dev=abs(100*k/n-5), usable_pct=100*us/att)]
cov2 <- H1[N<=25, .(ncov=sum(n_cov), kcov=sum(k_cov)), by=method][, covdev := abs(100*kcov/ncov-95)]
Z2 <- merge(lev2[, .(method,dev,usable_pct)], cov2[, .(method,covdev)], by="method", all.x=TRUE)
Z2[, `:=`(r_lev=rank(dev), r_cov=rank(covdev, na.last="keep"), r_us=rank(-usable_pct))]
Z2[, rango_medio := rowMeans(cbind(r_lev,r_cov,r_us), na.rm=TRUE)]
print(Z2[order(rango_medio)][, .(method, rango_medio=round(rango_medio,4))])
cat("\nDONE_E\n"
)