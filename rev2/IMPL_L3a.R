.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
setDTthreads(1)
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
z <- qnorm(0.975)
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)
COND <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232)
AR1 <- fread(file.path(OUT,"R2_ar1_cond_true.csv"))
cat("== R2_ar1_cond_true.csv colonne ==\n"); print(names(AR1)); print(dim(AR1))

fs <- list.files(RES, pattern="_mmlb_c[0-9]+[.]rds$", full.names=TRUE)
key <- sub("_c[0-9]+[.]rds$","",basename(fs))
cat(sprintf("file _mmlb_ = %d ; celle = %d ; celle con 10 chunk = %d\n",
            length(fs), length(unique(key)), sum(table(key)==10)))
cells <- names(which(table(key)==10))

work <- function(cl) {
  d <- rbindlist(lapply(fs[key==cl], readRDS), fill=TRUE)
  d <- d[term=="time:trt" & method %in% WALD & converged &
         is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10 & is.finite(p)]
  if (!nrow(d)) return(NULL)
  d[, model:=as.integer(model)]
  tm <- MARG[as.character(d$model)]; tc <- COND[as.character(d$model)]
  if (d$dgp[1]=="GLMM-AR1") {
    m <- AR1[match(paste(d$model,d$ntime), paste(AR1$model,AR1$ntime)), t_cond_true]
    tc <- m
  }
  d[, `:=`(t_marg=tm, t_cond=tc)]
  # limite empirico PER CELLA = media delle stime del metodo in quella cella
  d[, t_emp := mean(estimate), by=.(method)]
  s <- d[, .(n=.N, mean_est=mean(estimate), sd_est=sd(estimate), mean_se=mean(se),
             t_marg=t_marg[1], t_cond=t_cond[1], t_emp=t_emp[1],
             k_marg=sum(abs(estimate-t_marg)<=z*se),
             k_cond=sum(abs(estimate-t_cond)<=z*se),
             k_emp =sum(abs(estimate-t_emp )<=z*se),
             sbias_marg=sum(estimate-t_marg), sbias_cond=sum(estimate-t_cond),
             sabs_marg=sum(abs(estimate-t_marg)), sabs_cond=sum(abs(estimate-t_cond)),
             sq_marg=sum((estimate-t_marg)^2), sq_cond=sum((estimate-t_cond)^2),
             sq_emp=sum((estimate-t_emp)^2),
             spe_marg=sum(abs(estimate-t_marg)/abs(t_marg)),
             spe_cond=sum(abs(estimate-t_cond)/abs(t_cond))),
         by=.(dgp,model,N,ntime,method)]
  raw <- d[method=="glmer", .(dgp,model,N,ntime,estimate,se,t_marg,t_cond)]
  list(s=s, raw=raw)
}
res <- mclapply(cells, work, mc.cores=40)
bad <- vapply(res, is.null, TRUE); cat("celle vuote:", sum(bad), "\n")
S <- rbindlist(lapply(res[!bad], `[[`, "s"))
G <- rbindlist(lapply(res[!bad], `[[`, "raw"))
fwrite(S, file.path(OUT,"IMPL_L3a_cells.csv"))
cat(sprintf("\nS: %d righe ; G (glmer repliche): %d righe ; celle glmer: %d\n",
            nrow(S), nrow(G), nrow(S[method=="glmer"])))
print(S[method=="glmer", .(celle=.N, n=sum(n)), by=dgp][order(dgp)])
print(S[, .(celle=.N, n=sum(n)), by=method][order(method)])

## ---------- targets extra sul glmer ----------
# limite empirico "plim" convenzione agg_all5: media stime a N massimo per (dgp,model,ntime)
pl <- G[, .SD[N==max(N)], by=.(dgp,model,ntime)][, .(t_plim=mean(estimate)), by=.(dgp,model,ntime)]
G <- merge(G, pl, by=c("dgp","model","ntime"))
G[, t_empcell := mean(estimate), by=.(dgp,model,N,ntime)]
G[, `:=`(c_marg=abs(estimate-t_marg)<=z*se, c_cond=abs(estimate-t_cond)<=z*se,
         c_plim=abs(estimate-t_plim)<=z*se, c_emp=abs(estimate-t_empcell)<=z*se)]

pooled <- function(x) 100*mean(x)
cat("\n===== 1. GLMM: bias e copertura contro i tre bersagli =====\n")
cat(sprintf("n repliche glmer = %d su %d celle\n", nrow(G), nrow(S[method=="glmer"])))
# (i) pesata per repliche
cat("\n-- (i) POOLED sulle repliche --\n")
cat(sprintf("bias marginale   = %.6f\nbias condizionale= %.6f\nbias vs emp.cell = %.6f\n",
  mean(G$estimate-G$t_marg), mean(G$estimate-G$t_cond), mean(G$estimate-G$t_empcell)))
cat(sprintf("|bias| marginale = %.6f\n|bias| condiz.   = %.6f\n",
  mean(abs(G$estimate-G$t_marg)), mean(abs(G$estimate-G$t_cond))))
cat(sprintf("copertura pooled tre meccanismi: marg=%.4f cond=%.4f plim=%.4f empcell=%.4f\n",
  pooled(G$c_marg), pooled(G$c_cond), pooled(G$c_plim), pooled(G$c_emp)))
print(G[, .(n=.N, cov_marg=pooled(c_marg), cov_cond=pooled(c_cond),
            cov_plim=pooled(c_plim), cov_empcell=pooled(c_emp),
            mcse_marg=100*sqrt(mean(c_marg)*(1-mean(c_marg))/.N),
            mcse_cond=100*sqrt(mean(c_cond)*(1-mean(c_cond))/.N),
            bias_marg=mean(estimate-t_marg), bias_cond=mean(estimate-t_cond)), by=dgp][order(dgp)])
# (ii) media per cella
Cg <- G[, .(n=.N, cov_marg=pooled(c_marg), cov_cond=pooled(c_cond), cov_plim=pooled(c_plim),
            cov_emp=pooled(c_emp), b_marg=mean(estimate-t_marg), b_cond=mean(estimate-t_cond),
            b_emp=mean(estimate-t_empcell),
            ab_marg=mean(abs(estimate-t_marg)), ab_cond=mean(abs(estimate-t_cond)),
            pe_marg=mean(abs(estimate-t_marg)/abs(t_marg)),
            pe_cond=mean(abs(estimate-t_cond)/abs(t_cond))),
        by=.(dgp,model,N,ntime)]
cat("\n-- (ii) MEDIA PER CELLA (non pesata) --\n")
cat(sprintf("celle = %d\nbias marginale   = %.6f\nbias condizionale= %.6f\nbias vs emp.cell = %.6f\n",
  nrow(Cg), mean(Cg$b_marg), mean(Cg$b_cond), mean(Cg$b_emp)))
cat(sprintf("|bias| marginale = %.6f\n|bias| condiz.   = %.6f\nriduzione = %.4f%%\n",
  mean(Cg$ab_marg), mean(Cg$ab_cond), 100*(1-mean(Cg$ab_cond)/mean(Cg$ab_marg))))
cat(sprintf("copertura media per cella: marg=%.4f cond=%.4f plim=%.4f empcell=%.4f\n",
  mean(Cg$cov_marg), mean(Cg$cov_cond), mean(Cg$cov_plim), mean(Cg$cov_emp)))
print(Cg[, .(celle=.N, cov_marg=mean(cov_marg), cov_cond=mean(cov_cond), cov_plim=mean(cov_plim),
             cov_emp=mean(cov_emp), b_marg=mean(b_marg), b_cond=mean(b_cond),
             ab_marg=mean(ab_marg), ab_cond=mean(ab_cond)), by=dgp][order(dgp)])
cat("\n-- quota della differenza di bias segnato attribuibile alla convenzione --\n")
bp_m <- mean(G$estimate-G$t_marg); bp_c <- mean(G$estimate-G$t_cond)
bc_m <- mean(Cg$b_marg); bc_c <- mean(Cg$b_cond)
cat(sprintf("pooled : (%.6f-%.6f)/%.6f = %.4f%%\n", bp_m,bp_c,bp_m, 100*(bp_m-bp_c)/bp_m))
cat(sprintf("percella:(%.6f-%.6f)/%.6f = %.4f%%\n", bc_m,bc_c,bc_m, 100*(bc_m-bc_c)/bc_m))
cat(sprintf("MAPE: pooled marg=%.4f%% cond=%.4f%% quota=%.4f%% | percella marg=%.4f%% cond=%.4f%% quota=%.4f%%\n",
  100*mean(abs(G$estimate-G$t_marg)/abs(G$t_marg)), 100*mean(abs(G$estimate-G$t_cond)/abs(G$t_cond)),
  100*(1-mean(abs(G$estimate-G$t_cond)/abs(G$t_cond))/mean(abs(G$estimate-G$t_marg)/abs(G$t_marg))),
  100*mean(Cg$pe_marg), 100*mean(Cg$pe_cond), 100*(1-mean(Cg$pe_cond)/mean(Cg$pe_marg))))

cat("\n===== 2. T-3.13-06 denominatori =====\n")
FG <- S[method=="gee_fg", .(dgp,model,N,ntime,cov_fg=100*k_marg/n)]
GL <- S[method=="glmer", .(dgp,model,N,ntime,cov_gl=100*k_marg/n)]
J <- merge(FG,GL,by=c("dgp","model","N","ntime"))
gain <- Cg[, .(dgp,model,N,ntime,g=cov_cond-cov_marg)]
J <- merge(J, gain, by=c("dgp","model","N","ntime"))
cat(sprintf("(a) guadagno marg->cond, media per cella = %.4f pp\n", mean(J$g)))
cat(sprintf("(b) divario FG - GLMM (marginale), media per cella = %.4f pp\n", mean(J$cov_fg-J$cov_gl)))
cat(sprintf("(c) scarto GLMM dal nominale 95 (marginale), media per cella = %.4f pp\n", mean(95-J$cov_gl)))
print(J[, .(celle=.N, gain=mean(g), gap_fg=mean(cov_fg-cov_gl), gap_95=mean(95-cov_gl),
            q_su_fg=100*mean(g)/mean(cov_fg-cov_gl), q_su_95=100*mean(g)/mean(95-cov_gl)), by=dgp][order(dgp)])
cat(sprintf("quota su (b) = %.4f%% ; quota su (c) = %.4f%%\n",
  100*mean(J$g)/mean(J$cov_fg-J$cov_gl), 100*mean(J$g)/mean(95-J$cov_gl)))

cat("\n===== 3. T-3.13-08 pooled + marginale per meccanismo =====\n")
H <- fread(file.path(OUT,"REV_H1_targets.csv"))[method=="glmer"]
cat("da REV_H1_targets.csv (riferimento incrociato):\n")
print(H[, .(n=sum(n), cells=.N, marg=100*sum(k_marg)/sum(n), margdgp=100*sum(k_margdgp)/sum(n),
            cond=100*sum(k_cond)/sum(n), plim=100*sum(k_plim)/sum(n)), by=dgp][order(dgp)])
cat(sprintf("pooled marg (mio calcolo) = %.4f ; pooled cond = %.4f\n", pooled(G$c_marg), pooled(G$c_cond)))
# recupero del bersaglio marginale specifico per meccanismo su GLMM-native via griglia
Hn <- H[dgp=="GLMM-native", .(k=sum(k_margdgp), n=sum(n)), by=model]
Gn <- G[dgp=="GLMM-native"]
cat("\nrecupero t_marg specifico GLMM-native per scenario (griglia 1e-5):\n")
for (m in sort(unique(Gn$model))) {
  gg <- Gn[model==m]; tgt <- Hn[model==m, k]
  grid <- seq(if (m==4) 0.030 else 0.150, if (m==4) 0.080 else 0.250, by=1e-5)
  kk <- vapply(grid, function(t) sum(abs(gg$estimate-t)<=z*gg$se), 1L)
  ok <- grid[kk==tgt]
  cat(sprintf("  model %d: k_margdgp=%d ; t compatibili in [%.5f, %.5f] (n=%d valori)\n",
      m, tgt, if(length(ok)) min(ok) else NA_real_, if(length(ok)) max(ok) else NA_real_, length(ok)))
}

cat("\n===== 4. §3.14 RMSE =====\n")
S[, `:=`(rmse_cell=sqrt(sq_marg/n))]
R1 <- S[, .(celle=.N, n=sum(n), rmse_percella=mean(rmse_cell),
            rmse_pooled=sqrt(sum(sq_marg)/sum(n))), by=.(dgp,method)]
print(dcast(R1, method~dgp, value.var="rmse_percella"))
print(dcast(R1, method~dgp, value.var="rmse_pooled"))
print(dcast(R1, method~dgp, value.var="celle"))
CORR <- c("mmm","gee_fg","gee_md","gee_kc")
for (nm in c("rmse_percella","rmse_pooled")) {
  cat("\nspread MMM vs tre GEE corretti,", nm, ":\n")
  for (dd in sort(unique(R1$dgp))) {
    v <- R1[dgp==dd & method %in% CORR][[nm]]; g <- R1[dgp==dd & method=="glmer"][[nm]]
    cat(sprintf("  %-12s min=%.6f max=%.6f spread=%.4f%% ; glmer=%.6f fattore=%.4f\n",
                dd, min(v), max(v), 100*(max(v)-min(v))/min(v), g, g/min(v)))
  }
}
cat("\nShrinkage RMSE:\n")
print(R1[method=="shrink"])
cat("\nRMSE contro il bersaglio condizionale (per cella), glmer e shrink:\n")
print(S[method %in% c("glmer","shrink"), .(rmse_cond_percella=mean(sqrt(sq_cond/n)),
      rmse_cond_pooled=sqrt(sum(sq_cond)/sum(n))), by=.(dgp,method)][order(method,dgp)])

cat("\n===== 5. T-3.15-02 media stime glmer N=145 T=10 GenBinaryY =====\n")
print(S[method=="glmer" & dgp=="GenBinaryY" & N==145 & ntime==10,
        .(model,n,mean_est=round(mean_est,6),sd_est=round(sd_est,6))][order(model)])
cat(sprintf("media sui cinque scenari = %.6f\n",
  mean(S[method=="glmer" & dgp=="GenBinaryY" & N==145 & ntime==10, mean_est])))
cat(sprintf("N massimi presenti: %s\n", paste(sort(unique(S$N)), collapse=",")))

cat("\n===== 6. tabella §3.15 completa (bersaglio marginale per OGNI metodo) =====\n")
T315 <- S[, .(cov_pooled=100*sum(k_marg)/sum(n), cov_percella=mean(100*k_marg/n),
              celle=.N, n=sum(n)), by=.(method,dgp)]
print(dcast(T315, method~dgp, value.var="cov_pooled"))
cat("\nla stessa tabella se il glmer fosse scorato sul condizionale (gli altri sul marginale):\n")
T315b <- copy(T315); cc <- S[method=="glmer", .(cov_pooled=100*sum(k_cond)/sum(n)), by=dgp]
print(cc)
cat("\ncopertura shrink contro marginale e contro condizionale:\n")
print(S[method=="shrink", .(cov_marg_pooled=100*sum(k_marg)/sum(n), cov_cond_pooled=100*sum(k_cond)/sum(n),
                            cov_marg_percella=mean(100*k_marg/n)), by=dgp][order(dgp)])
fwrite(R1, file.path(OUT,"IMPL_L3a_rmse.csv"))
fwrite(T315, file.path(OUT,"IMPL_L3a_cov315.csv"))
fwrite(Cg, file.path(OUT,"IMPL_L3a_glmer_cells.csv"))
cat("\nFATTO\n")
