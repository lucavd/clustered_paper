
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
setDTthreads(1)
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink"); LRT <- "glmer_lrt"
ALL8 <- c(WALD, LRT)
NP  <- c("1"=4,"2"=5,"3"=4,"4"=4,"5"=5)
NPW <- c("1"=2,"2"=3,"3"=2,"4"=2,"5"=3)
zq  <- qnorm(0.975)

fs  <- list.files(RES, pattern="[.]rds$")
key <- sub("_c[0-9]+[.]rds$", "", fs)
tab <- table(key); grps <- names(tab)
H0g <- grps[grepl("_manual$", grps)]
cat("gruppi totali:", length(grps), " gruppi H0:", length(H0g),
    " con 10 chunk:", sum(tab[H0g]==10), " file H0:", sum(key %in% H0g), "\n")
cat("prefissi:", paste(sort(unique(sub("_m[0-9]+_.*$","",H0g))), collapse=" | "), "\n")

one <- function(g) {
  f <- sort(file.path(RES, fs[key==g]))
  L <- vector("list", length(f))
  for (i in seq_along(f)) {
    x <- try(as.data.table(readRDS(f[i])), silent=TRUE)
    if (inherits(x,"try-error")) return(NULL)
    x <- x[term=="time:trt"]; if (!nrow(x)) next
    x[, chunk := i]
    kp <- intersect(c("method","N","ntime","label","model","converged","sim_id",
                      "estimate","se","p","dgp","chunk"), names(x))
    L[[i]] <- x[, ..kp]
  }
  d <- rbindlist(L, fill=TRUE); if (!nrow(d)) return(NULL)
  d[, model := as.integer(model)]
  d[, rid := paste0(chunk,"_",sim_id)]
  d[, converged := as.logical(converged)]; d[is.na(converged), converged := FALSE]
  d[, fin := is.finite(estimate) & is.finite(se) & is.finite(p)]
  d[, usable  := ifelse(method==LRT, converged & is.finite(p),
                        converged & fin & abs(estimate) < 10 & se > 0 & se < 10)]
  d[, us_conv := converged & is.finite(p)]
  d[, us_fin  := ifelse(method==LRT, converged & is.finite(p), converged & fin)]
  d[, us_nose := ifelse(method==LRT, converged & is.finite(p),
                        converged & fin & abs(estimate) < 10 & se > 0)]
  d[, us_nobig:= ifelse(method==LRT, converged & is.finite(p),
                        converged & fin & se > 0 & se < 10)]
  m <- d$model[1]; nn <- d$N[1]; tt <- d$ntime[1]; dg <- d$dgp[1]
  p_tot <- NP[[as.character(m)]]; p_win <- NPW[[as.character(m)]]
  nobs_exp <- nn*(0.7*tt + 0.3*tt/2)
  df_bw <- nn - p_tot; df_k2 <- nn - 2; df_win <- nobs_exp - p_tot
  df_bwin <- nobs_exp - nn - p_tot + 1; df_dec <- nobs_exp - nn - 2
  df_decm <- nobs_exp - nn - p_win
  n_attempt <- uniqueN(d$rid)
  u <- d[usable==TRUE]
  cnt <- u[, .(nm=uniqueN(method)), by=rid]; inter <- cnt[nm==length(ALL8), rid]
  d[, in_inter := rid %in% inter]
  cellrow <- u[, {
    t_ <- abs(estimate/se); isl <- method[1]==LRT
    .(n=.N, k_norm=sum(p<0.05),
      k_bw  = if (isl) sum(p<0.05) else if (df_bw  >0) sum(2*pt(-t_,df_bw )<0.05) else NA_integer_,
      k_k2  = if (isl) sum(p<0.05) else if (df_k2  >0) sum(2*pt(-t_,df_k2 )<0.05) else NA_integer_,
      k_win = if (isl) sum(p<0.05) else if (df_win >0) sum(2*pt(-t_,df_win)<0.05) else NA_integer_,
      k_bwin= if (isl) sum(p<0.05) else if (df_bwin>0) sum(2*pt(-t_,df_bwin)<0.05) else NA_integer_,
      k_dec = if (isl) sum(p<0.05) else if (df_dec >0) sum(2*pt(-t_,df_dec)<0.05) else NA_integer_,
      k_decm= if (isl) sum(p<0.05) else if (df_decm>0) sum(2*pt(-t_,df_decm)<0.05) else NA_integer_,
      k_z196= if (isl) sum(p<0.05) else sum(t_ > zq))
  }, by=method]
  cellrow[, `:=`(dgp=dg, model=m, N=nn, ntime=tt, n_attempt=n_attempt,
                 nobs_exp=nobs_exp, df_bw=df_bw, df_k2=df_k2, df_win=df_win,
                 df_bwin=df_bwin, df_dec=df_dec, df_decm=df_decm)]
  # ---- 2x2 appaiate mmm vs gee_fg sotto piu definizioni
  P <- d[method %in% c("mmm","gee_fg")]
  P[, t_ := abs(estimate/se)]
  P[, rej_norm := p < 0.05]
  P[, rej_z196 := t_ > zq]
  P[, rej_tbw  := if (df_bw  > 0) 2*pt(-t_, df_bw ) < 0.05 else NA]
  P[, rej_tdec := if (df_dec > 0) 2*pt(-t_, df_dec) < 0.05 else NA]
  mk2 <- function(flag, rej, keyv=c("rid")) {
    S <- P[get(flag)==TRUE, c(keyv,"method",rej), with=FALSE]
    setnames(S, rej, "R")
    W <- dcast(S, as.formula(paste(paste(keyv,collapse="+"), "~ method")),
               value.var="R", fun.aggregate=function(z) z[1])
    if (!all(c("mmm","gee_fg") %in% names(W))) return(data.table(npair=0L,b=0L,cc=0L,km=0L,kf=0L))
    W <- W[!is.na(mmm) & !is.na(gee_fg)]
    data.table(npair=nrow(W), b=sum(W$mmm & !W$gee_fg), cc=sum(!W$mmm & W$gee_fg),
               km=sum(W$mmm), kf=sum(W$gee_fg))
  }
  defs <- list(
    own_norm   = mk2("usable",  "rej_norm"),
    inter_norm = mk2("in_inter","rej_norm"),
    own_z196   = mk2("usable",  "rej_z196"),
    own_tbw    = mk2("usable",  "rej_tbw"),
    own_tdec   = mk2("usable",  "rej_tdec"),
    conv_norm  = mk2("us_conv", "rej_norm"),
    fin_norm   = mk2("us_fin",  "rej_norm"),
    nose_norm  = mk2("us_nose", "rej_norm"),
    nobig_norm = mk2("us_nobig","rej_norm"),
    own_norm_simid = mk2("usable","rej_norm", keyv="sim_id")
  )
  PD <- rbindlist(lapply(names(defs), function(nm) cbind(def=nm, defs[[nm]])))
  PD[, `:=`(dgp=dg, model=m, N=nn, ntime=tt, n_attempt=n_attempt)]
  # k per metodo sull'insieme proprio, per il confronto con §3.1
  own_k <- u[method %in% c("mmm","gee_fg"), .(n=.N, k=sum(p<0.05)), by=method]
  OK <- data.table(dgp=dg, model=m, N=nn, ntime=tt,
                   n_mmm=own_k[method=="mmm", n][1], k_mmm=own_k[method=="mmm", k][1],
                   n_fg =own_k[method=="gee_fg", n][1], k_fg=own_k[method=="gee_fg", k][1])
  list(cell=cellrow, paired=PD, own=OK,
       simid=data.table(dgp=dg, model=m, N=nn, ntime=tt, n_uniq_simid=uniqueN(d$sim_id),
                        n_rid=n_attempt, n_inter=length(inter),
                        max_dup=max(d[method=="mmm", .N, by=sim_id]$N)))
}
t0 <- Sys.time()
Z <- mclapply(H0g, function(g) try(one(g), silent=TRUE), mc.cores=45)
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),2), "min\n")
isbad <- vapply(Z, function(z) inherits(z,"try-error") || is.null(z), TRUE)
cat("falliti:", sum(isbad), "\n"); if (any(isbad)) print(Z[isbad][[1]])
Z <- Z[!isbad]
CELL <- rbindlist(lapply(Z,`[[`,"cell"), fill=TRUE)
PD   <- rbindlist(lapply(Z,`[[`,"paired"), fill=TRUE)
OWN  <- rbindlist(lapply(Z,`[[`,"own"), fill=TRUE)
SID  <- rbindlist(lapply(Z,`[[`,"simid"), fill=TRUE)
fwrite(CELL, file.path(OUT,"TR_cells.csv")); fwrite(PD, file.path(OUT,"TR_paired_defs.csv"))
fwrite(OWN, file.path(OUT,"TR_own.csv")); fwrite(SID, file.path(OUT,"TR_simid.csv"))
cat("celle:", nrow(unique(CELL[,.(dgp,model,N,ntime)])), " righe CELL:", nrow(CELL),
    " righe PD:", nrow(PD), "\n")
cat("\n=== sim_id: duplicati entro cella ===\n")
cat("celle con n_uniq_simid == n_rid:", sum(SID$n_uniq_simid==SID$n_rid), "/", nrow(SID), "\n")
cat("max replicazioni dello stesso sim_id entro cella (metodo mmm):", max(SID$max_dup), "\n")
cat("\nDONE_A\n")
