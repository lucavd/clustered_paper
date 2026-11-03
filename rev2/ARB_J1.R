
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
z <- qnorm(0.975)
MARG <- c(0.200,0.200,0.200,0.050,0.200)
CFIX <- c(0.232,0.232,0.232,0.058,0.232)

## ---- le tre convenzioni di bersaglio condizionale in circolazione ----
arF <- fread(path.expand("~/MMM/tables/F_truth_ar1.csv"))[, .(tF=mean(cond_used)), by=.(model,ntime)]
arR <- fread(path.expand("~/MMM/rev2/R2_ar1_cond_true.csv"))[, .(model,ntime,tR=t_cond_true)]
gn  <- fread(path.expand("~/MMM/tables/F_truth_glmmnative.csv"))
PL  <- fread(path.expand("~/MMM/tables/O_glmer_plim.csv"))[, .(dgp,model,ntime,plim)]

grid <- CJ(dgp=c("GenBinaryY","GLMM-AR1","GLMM-native"), model=1:5, ntime=c(2,3,4,6,8,10))
grid[, tm := MARG[model]]
grid[, c1 := CFIX[model]]
grid <- merge(grid, arF, by=c("model","ntime"), all.x=TRUE)
grid <- merge(grid, arR, by=c("model","ntime"), all.x=TRUE)
grid[, c2 := ifelse(dgp=="GLMM-AR1" & is.finite(tF), tF, CFIX[model])]
grid[, c3 := ifelse(dgp=="GLMM-AR1" & is.finite(tR), tR, CFIX[model])]
## marginale per meccanismo (GLMM-native varia con T)
gnl <- melt(gn, id.vars="ntime", measure.vars=intersect(c("m_main","m_m4"), names(gn)))
grid[, tmd := tm]
for (tt in unique(gn$ntime)) {
  mm <- gn[ntime==tt]
  grid[dgp=="GLMM-native" & ntime==tt & model!=4, tmd := mm$m_main[1]]
  if ("m_m4" %in% names(gn)) grid[dgp=="GLMM-native" & ntime==tt & model==4, tmd := mm$m_m4[1]]
}
grid <- merge(grid, PL, by=c("dgp","model","ntime"), all.x=TRUE)
fwrite(grid, file.path(OUT,"ARB_grid_targets.csv"))
cat("griglia bersagli:", nrow(grid), "\n")
cat("AR1: F_truth vs R2_ar1_cond_true, scarto per (model,ntime):\n")
print(grid[dgp=="GLMM-AR1", .(model,ntime,tF=round(tF,5),tR=round(tR,5),diff=round(tR-tF,5))][order(model,ntime)])

MSET <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
fs <- list.files(RES, pattern="_mmlb_c[0-9]+[.]rds$")
cat("chunk H1 (mmlb):", length(fs), "\n")

one <- function(f) {
  x <- try(as.data.table(readRDS(file.path(RES,f))), silent=TRUE)
  if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method %in% MSET]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  u <- x[converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
         abs(estimate) < 10 & se > 0 & se < 10]
  if (!nrow(u)) return(NULL)
  u <- merge(u, grid, by=c("dgp","model","ntime"), all.x=TRUE)
  u[, tp := ifelse(method=="glmer" & is.finite(plim), plim, tm)]
  u[, .(n=.N, s_est=sum(estimate), s_est2=sum(estimate^2), s_se=sum(se),
        k_marg=sum(abs(estimate-tm)<=z*se), k_mdgp=sum(abs(estimate-tmd)<=z*se),
        k_c1=sum(abs(estimate-c1)<=z*se), k_c2=sum(abs(estimate-c2)<=z*se),
        k_c3=sum(abs(estimate-c3)<=z*se), k_plim=sum(abs(estimate-tp)<=z*se),
        sb_marg=sum(estimate-tm), sb_c1=sum(estimate-c1), sb_c2=sum(estimate-c2), sb_c3=sum(estimate-c3),
        sse_marg=sum((estimate-tm)^2), sse_c1=sum((estimate-c1)^2),
        sse_c2=sum((estimate-c2)^2), sse_c3=sum((estimate-c3)^2),
        sae_marg=sum(abs(estimate-tm)), sae_c1=sum(abs(estimate-c1)),
        sae_c2=sum(abs(estimate-c2)), sae_c3=sum(abs(estimate-c3)),
        spe_marg=sum(abs(estimate-tm)/abs(tm)), spe_c1=sum(abs(estimate-c1)/abs(c1)),
        spe_c2=sum(abs(estimate-c2)/abs(c2)), spe_c3=sum(abs(estimate-c3)/abs(c3))),
    by=.(method,dgp,model,N,ntime)]
}
L <- mclapply(fs, function(f) try(one(f), silent=TRUE), mc.cores=44)
cat("chunk falliti:", sum(!sapply(L, is.data.table)), "\n")
P <- rbindlist(L[sapply(L, is.data.table)], fill=TRUE)
key <- c("method","dgp","model","N","ntime")
C <- P[, lapply(.SD, sum), by=key, .SDcols=setdiff(names(P), key)]
fwrite(C, file.path(OUT,"ARB_H1_cells.csv"))
cat("celle:", nrow(C), " per metodo:\n"); print(C[, .N, by=method])

G <- C[method=="glmer"]
cat("\n=== GLMER: celle", nrow(G), " repliche", sum(G$n), "===\n")
f6 <- function(x) sprintf("%.6f", x)
for (cv in c("marg","c1","c2","c3")) {
  sb <- G[[paste0("sb_",cv)]]; sae <- G[[paste0("sae_",cv)]]
  cat(sprintf("%-5s bias pesato %s | media per cella %s | |bias| medio per cella %s\n",
      cv, f6(sum(sb)/sum(G$n)), f6(mean(sb/G$n)), f6(mean(abs(sb/G$n)))))
}
cat("\nriduzione |bias| per cella (marg -> cond), per convenzione:\n")
am <- mean(abs(G$sb_marg/G$n))
for (cv in c("c1","c2","c3")) {
  ac <- mean(abs(G[[paste0("sb_",cv)]]/G$n))
  cat(sprintf("  %s: %.5f -> %.5f  riduzione %.4f%%\n", cv, am, ac, 100*(am-ac)/am))
}
cat("\nriduzione MAPE (pesata), con e senza M4:\n")
for (cv in c("c1","c2","c3")) {
  a <- 100*sum(G$spe_marg)/sum(G$n); b <- 100*sum(G[[paste0("spe_",cv)]])/sum(G$n)
  G4 <- G[model!=4]
  a4 <- 100*sum(G4$spe_marg)/sum(G4$n); b4 <- 100*sum(G4[[paste0("spe_",cv)]])/sum(G4$n)
  cat(sprintf("  %s: MAPE %.4f -> %.4f  riduzione %.4f%%   | senza M4: %.4f -> %.4f  riduzione %.4f%%\n",
      cv, a, b, 100*(a-b)/a, a4, b4, 100*(a4-b4)/a4))
}
cat("\n=== COPERTURA glmer per meccanismo (%) ===\n")
cv <- G[, .(n=sum(n), marg=100*sum(k_marg)/sum(n), mdgp=100*sum(k_mdgp)/sum(n),
            c1=100*sum(k_c1)/sum(n), c2=100*sum(k_c2)/sum(n), c3=100*sum(k_c3)/sum(n),
            plim=100*sum(k_plim)/sum(n)), by=dgp]
print(cv, digits=6)
cat("\npooled glmer:\n")
print(G[, .(n=sum(n), marg=100*sum(k_marg)/sum(n), mdgp=100*sum(k_mdgp)/sum(n),
            c1=100*sum(k_c1)/sum(n), c2=100*sum(k_c2)/sum(n), c3=100*sum(k_c3)/sum(n),
            plim=100*sum(k_plim)/sum(n))], digits=6)
cat("MCSE copertura (pooled, c3):", sprintf("%.4f", 100*sqrt(0.9319*(1-0.9319)/sum(G$n))), "pp\n")
cat("MCSE copertura AR1 (c3):", sprintf("%.4f", 100*sqrt(0.9316*(1-0.9316)/G[dgp=="GLMM-AR1", sum(n)])), "pp\n")

cat("\n=== RMSE media per cella (bersaglio marginale) ===\n")
C[, rm := sqrt(sse_marg/n)]
rr <- C[, .(rmse=mean(rm), celle=.N, n=sum(n)), by=.(method,dgp)]
print(dcast(rr, method~dgp, value.var="rmse"), digits=6)
print(dcast(rr, method~dgp, value.var="celle"))
cat("\nspread fra i quattro marginali (mmm,gee_fg,gee_md,gee_kc):\n")
for (dd in unique(rr$dgp)) {
  v <- rr[dgp==dd & method %in% c("mmm","gee_fg","gee_md","gee_kc"), rmse]
  cat(sprintf("  %-12s diff*100 %.4f | rapporto (max/min-1)*100 %.4f\n", dd, 100*(max(v)-min(v)), 100*(max(v)/min(v)-1)))
}
cat("\nshrink RMSE e celle:\n"); print(rr[method=="shrink"], digits=6)
cat("\n=== §3.15 glmer media a N=145,T=10 ===\n")
print(C[method=="glmer" & N==145 & ntime==10, .(dgp,model,n,mean_est=s_est/n)], digits=6)
print(C[method=="glmer" & N==145 & ntime==10 & model!=4, .(n=sum(n), mean_est=sum(s_est)/sum(n)), by=dgp], digits=6)
cat("\nDONE_JOB1\n")
