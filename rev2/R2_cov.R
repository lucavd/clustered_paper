# R2_cov.R -- GIRO 2, contestazione 2: coverage di glmer e shrink contro TRE bersagli.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel); library(fastGHQuad)})
RES <- path.expand("~/MMM/results_full"); zq <- qnorm(0.975)
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)
COND <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232)

# verita condizionale VERA per scenario e per T. Per GLMM-native beta_c e' fisso: la
# verita condizionale E' 0.232 (0.058 su M4) per costruzione, non varia con T.
# Per AR1 e' calibrata per replica: la recupero con la struttura del generatore.
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_ar1_data.R"), envir=globalenv())
sys.source(file.path(RD,"ar1/generation_models_ar1.R"), envir=globalenv())
D2 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
D3 <- matrix(c(1,.01,.01, .01,.001,0, .01,0,.01),3,byrow=TRUE)
TL <- c(2,3,4,6,8,10)
ar1c <- rbindlist(lapply(1:5, function(m) rbindlist(lapply(TL, function(T) {
  ex <- if (m==1) "none" else if (m==5) "full" else "age"
  re <- if (m %in% c(1,2)) ~1 else if (m==5) ~1+time+dose else ~1+time
  Dm <- if (m %in% c(1,2)) 1.0 else if (m==5) D3 else D2
  mf <- if (m==2) ~time*trt+age else if (m==5) ~time*trt+dose else ~time*trt
  vals <- numeric(0)
  for (rep in 1:6) {
    set.seed(20260805L + 1000L*m + rep)
    d <- tryCatch(.ar1_design(75, T, ex), error=function(e) NULL); if (is.null(d)) next
    b <- tryCatch(as.numeric(calibrate_beta_c_ar1(mf, re, .ar1_target(m,"mmlb"), Dm, AR1_SIGMA2, d)),
                  error=function(e) NULL)
    if (is.null(b)) next
    nm <- colnames(model.matrix(mf, model.frame(mf, d)))
    j <- match("time:trt", nm); if (!is.na(j)) vals <- c(vals, b[j])
  }
  data.table(model=m, ntime=T, t_cond_true=if (length(vals)) mean(vals) else NA_real_,
             n_calib=length(vals))
}))))
fwrite(ar1c, path.expand("~/MMM/rev2/R2_ar1_cond_true.csv"))
cat("verita condizionale AR1 ricalcolata (media su 6 disegni, N=75):\n")
print(as.data.frame(dcast(ar1c, model~ntime, value.var="t_cond_true")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,5) else x)]), row.names=FALSE)

P <- fread(path.expand("~/MMM/tables/O_glmer_plim.csv"))[, .(dgp,model,ntime,plim)]
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & label=="mmlb" & method %in% c("glmer","shrink","mmm","gee_fg")]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10,
    .(dgp,method,model,N,ntime,estimate,se)]
}
D <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)), fill=TRUE)
D <- merge(D, P, by=c("dgp","model","ntime"), all.x=TRUE)
D <- merge(D, ar1c[, .(model,ntime,t_cond_ar1=t_cond_true)], by=c("model","ntime"), all.x=TRUE)
D[, t_marg := MARG[as.character(model)]]
# condizionale vero: fisso per GLMM-native e GenBinaryY (0.232/0.058), calibrato per AR1
D[, t_cond := ifelse(dgp=="GLMM-AR1" & is.finite(t_cond_ar1), t_cond_ar1, COND[as.character(model)])]
D[, `:=`(c_marg = abs(estimate-t_marg)<=zq*se,
         c_cond = abs(estimate-t_cond)<=zq*se,
         c_plim = abs(estimate-plim)<=zq*se)]
cat("\n===== COVERAGE contro i TRE bersagli, con n =====\n")
tb <- D[, .(n=.N, vs_marg=round(100*mean(c_marg),2), vs_cond=round(100*mean(c_cond),2),
            vs_plim=round(100*mean(c_plim),2)), by=.(method,dgp)]
print(as.data.frame(tb[order(method,dgp)]), row.names=FALSE)
cat("\n===== solo glmer e shrink, affiancati ai valori pubblicati in §3.15 =====\n")
pub <- data.table(method=c("glmer","glmer","glmer","shrink","shrink","shrink"),
                  dgp=c("GenBinaryY","GLMM-native","GLMM-AR1","GenBinaryY","GLMM-native","GLMM-AR1"),
                  pubblicato=c(88.32,93.03,92.28,90.06,92.94,92.59))
m <- merge(tb[method %in% c("glmer","shrink")], pub, by=c("method","dgp"))
m[, `:=`(d_marg=round(vs_marg-pubblicato,2), d_cond=round(vs_cond-pubblicato,2),
         d_plim=round(vs_plim-pubblicato,2))]
print(as.data.frame(m[order(method,dgp)]), row.names=FALSE)
cat("\n===== bias segnato contro i tre bersagli (per §3.13) =====\n")
bb <- D[, .(n=.N, b_marg=round(mean(estimate-t_marg),4), b_cond=round(mean(estimate-t_cond),4),
            b_plim=round(mean(estimate-plim),4)), by=.(method,dgp)]
print(as.data.frame(bb[method %in% c("glmer","shrink")][order(method,dgp)]), row.names=FALSE)
cat(sprintf("\nglmer pooled: marg %.4f | cond %.4f | plim %.4f (n=%d)\n",
  mean(D[method=="glmer"]$estimate-D[method=="glmer"]$t_marg),
  mean(D[method=="glmer"]$estimate-D[method=="glmer"]$t_cond),
  mean(D[method=="glmer"]$estimate-D[method=="glmer"]$plim, na.rm=TRUE), nrow(D[method=="glmer"])))
cat(sprintf("shrink pooled: marg %.4f | cond %.4f | plim %.4f (n=%d)\n",
  mean(D[method=="shrink"]$estimate-D[method=="shrink"]$t_marg),
  mean(D[method=="shrink"]$estimate-D[method=="shrink"]$t_cond),
  mean(D[method=="shrink"]$estimate-D[method=="shrink"]$plim, na.rm=TRUE), nrow(D[method=="shrink"])))
cat("\n===== shrink: e' il glmer riscalato? rapporto stima/SE sulle repliche appaiate =====\n")
one2 <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & label=="mmlb" & method %in% c("glmer","shrink")]
  if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  w <- dcast(x[, .(dgp,model,N,ntime,sim_id,method,estimate,se)],
             dgp+model+N+ntime+sim_id~method, value.var=c("estimate","se"))
  w[is.finite(estimate_glmer) & is.finite(estimate_shrink) & is.finite(se_glmer) & is.finite(se_shrink)]
}
W <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one2, mc.cores=45)), fill=TRUE)
W[, `:=`(r_est=estimate_shrink/estimate_glmer, r_se=se_shrink/se_glmer)]
cat(sprintf("coppie: %d ; mediana stima_shrink/stima_glmer %.6f ; mediana SE_shrink/SE_glmer %.6f\n",
  nrow(W), median(W$r_est, na.rm=TRUE), median(W$r_se, na.rm=TRUE)))
cat(sprintf("|r_est - r_se| : mediana %.2e, 99o pct %.2e, max %.2e -> la statistica di Wald e' invariata\n",
  median(abs(W$r_est-W$r_se), na.rm=TRUE), quantile(abs(W$r_est-W$r_se), .99, na.rm=TRUE),
  max(abs(W$r_est-W$r_se), na.rm=TRUE)))
cat(sprintf("mediana del fattore di riscalatura per T: "))
print(W[, .(r=round(median(r_est,na.rm=TRUE),4)), by=ntime][order(ntime)])
cat("\n===== shrink: il suo bersaglio corretto e' MARGINALE (e' il glmer riscalato verso il marginale) =====\n")
cat("agg_all5.R non distingue: t_gen = t_marg per shrink come per tutti. Coverage di shrink contro\n")
cat("il marginale (0.200) e' quindi il bersaglio GIUSTO; il difetto M-05 tocca solo glmer.\n")
fwrite(tb, path.expand("~/MMM/rev2/R2_cov_targets.csv"))
cat("\nDONE_COV\n")
