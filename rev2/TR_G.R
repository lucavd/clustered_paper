
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(8)
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
DG <- c("GenBinaryY","GLMM-native","GLMM-AR1")
fx <- function(x) fifelse(x=="GLMMnative","GLMM-native", fifelse(x=="GLMMAR1","GLMM-AR1", x))
cut0 <- as.POSIXct("2026-08-08 10:18:26")

cat("###### A. era dei chunk per fattore di disegno ######\n")
fi <- as.data.table(file.info(list.files(RES, full.names=TRUE, pattern="[.]rds$")), keep.rownames="path")
fi[, bn := basename(path)][, grp := sub("_c[0-9]+[.]rds$","",bn)]
fi[, post := mtime > cut0]
G <- fi[, .(npost=sum(post)), by=grp][, era := fifelse(npost==10,"post","pre")]
G[, `:=`(pref=sub("_m[0-9]+_.*$","",grp),
         model=as.integer(sub(".*_m([0-9]+)_.*","\\1",grp)),
         N=as.integer(sub(".*_N([0-9]+)_.*","\\1",grp)),
         T=as.integer(sub(".*_T([0-9]+)_.*","\\1",grp)),
         arm=fifelse(grepl("_manual$",grp),"manual","mmlb"))]
cat("era x T:\n"); print(dcast(G[arm=="manual"], T ~ era, value.var="grp", fun.aggregate=length))
cat("era x N:\n"); print(dcast(G[arm=="manual"], N ~ era, value.var="grp", fun.aggregate=length))
cat("era x model:\n"); print(dcast(G[arm=="manual"], model ~ era, value.var="grp", fun.aggregate=length))
cat("era x meccanismo:\n"); print(dcast(G[arm=="manual"], pref ~ era, value.var="grp", fun.aggregate=length))

cat("\n###### B. decomposizione: file = [1152 gruppi pre, invariati] + [1728 gruppi, versione vecchia] ######\n")
S <- fread(file.path(OUT,"TR_prepost.csv"))
F <- data.table(dgp=DG, npair_f=c(430139,434919,433219), b_f=c(9825,8415,9056), c_f=c(7811,6548,6371))
pre  <- S[post==FALSE, .(dgp, npair_pre=npair, b_pre=b, c_pre=cc)]
post <- S[post==TRUE,  .(dgp, npair_post=npair, b_post=b, c_post=cc)]
X <- merge(merge(F, pre, by="dgp"), post, by="dgp")
X[, `:=`(npair_imp = npair_f - npair_pre, b_imp = b_f - b_pre, c_imp = c_f - c_pre)]
X[, `:=`(disc_imp = 100*(b_imp+c_imp)/npair_imp, gap_imp = 100*(b_imp-c_imp)/npair_imp,
         disc_post = 100*(b_post+c_post)/npair_post, gap_post = 100*(b_post-c_post)/npair_post,
         disc_pre  = 100*(b_pre+c_pre)/npair_pre,   gap_pre  = 100*(b_pre-c_pre)/npair_pre)]
print(X[match(DG,dgp), .(dgp, npair_imp, npair_post, scarto_npair=npair_imp-npair_post,
   b_imp, b_post, c_imp, c_post,
   disc_imp=round(disc_imp,4), disc_post=round(disc_post,4), disc_pre=round(disc_pre,4),
   gap_imp=round(gap_imp,4), gap_post=round(gap_post,4))])

cat("\n###### C. ramo 1: usabilita su TUTTI i 240 punti (celle mai stimate a zero) ######\n")
R <- fread(file.path(OUT,"REV_cells.csv")); R[, dgp := fx(dgp)]
H0 <- R[label=="manual" & set=="own"]
pts <- unique(H0[, .(dgp,model,N,ntime,n_attempt)])
meth <- unique(H0$method)
FULL <- CJ(key=1:nrow(pts), method=meth)[, .(idx=key, method)]
grid <- cbind(pts[FULL$idx], method=FULL$method)
Z <- merge(grid, H0[, .(dgp,model,N,ntime,method,n_usable,n,k_rej)],
           by=c("dgp","model","N","ntime","method"), all.x=TRUE)
Z[is.na(n_usable), n_usable := 0L]
for (tt in c(2,3)) {
  x <- Z[ntime==tt, .(n_punti=.N, att=sum(n_attempt), us=sum(n_usable),
                      usable_pct=100*sum(n_usable)/sum(n_attempt)), by=method]
  cat("\n  T =", tt, " usabilita su tutti i punti:\n"); print(x[order(-usable_pct)][, .(method,n_punti,att,us,usable_pct=round(usable_pct,4))])
}
cat("\n  N<25 & T>=4, usabilita su tutti i punti:\n")
x <- Z[N<25 & ntime>=4, .(n_punti=.N, usable_pct=100*sum(n_usable)/sum(n_attempt)), by=method]
print(x[order(-usable_pct)][, .(method,n_punti,usable_pct=round(usable_pct,4))])
cat("\n  N<=20 (i 270 punti di §3.9), usabilita su tutti i punti:\n")
x <- Z[N<=20, .(n_punti=.N, usable_pct=100*sum(n_usable)/sum(n_attempt)), by=method]
print(x[order(-usable_pct)][, .(method,n_punti,usable_pct=round(usable_pct,4))])
cat("\nDONE_G\n")
