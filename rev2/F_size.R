
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
rd <- function(x){ z <- as.data.table(readRDS(x))
  z[term=="time:trt", .(estimate,se,p,method,N,ntime,model,label,converged,dgp)] }
fs <- list.files(RES, pattern="[.]rds$", full.names=TRUE)
k <- sub("_c[0-9]+[.]rds$","",basename(fs)); fs <- fs[k %in% names(which(table(k)==10))]
D <- rbindlist(mclapply(fs, rd, mc.cores=40)); D[, model := as.integer(model)]
keep <- function(d) rbind(
  d[method %in% WALD & converged & is.finite(estimate) & is.finite(se) &
    abs(estimate)<10 & se>0 & se<10 & is.finite(p)],
  d[method=="glmer_lrt" & converged & is.finite(p)])
H0 <- keep(D[label=="manual"]); H1 <- keep(D[label=="mmlb"]); rm(D); gc()

CUT <- H0[, .(q05 = quantile(p, 0.05, na.rm=TRUE, type=7), nH0=.N,
              n_p0 = sum(p==0), ties_at_q = sum(p==quantile(p,0.05,na.rm=TRUE,type=7))),
          by=.(dgp,model,N,ntime,method)]
cat("### cells (of 1440 per method-set) whose q05 is EXACTLY 0 -> strict p<q05 can never reject\n")
print(CUT[, .(ncell=.N, cells_q05_eq0=sum(q05==0), cells_ties=sum(ties_at_q>1),
              rows_p_eq_0=sum(n_p0)), by=method][order(method)])

S0 <- merge(H0, CUT, by=c("dgp","model","N","ntime","method"))
cat("\n### SIZE VALIDATION, strict (<) vs inclusive (<=), row-weighted, per dgp x method\n")
SV <- S0[, .(n=.N, size_lt=round(100*mean(p<q05),4), size_le=round(100*mean(p<=q05),4)),
         by=.(dgp,method)][order(dgp,method)]
print(SV)
cat("\nrange of size_lt over the 24:", round(min(SV$size_lt),4),"-",round(max(SV$size_lt),4),
    "\nrange of size_le over the 24:", round(min(SV$size_le),4),"-",round(max(SV$size_le),4),"\n")
cat("\n### 20 (dgp,method) EXCLUDING glmer and shrink, strict convention\n")
print(SV[!method %in% c("glmer","shrink"), .(n=.N, min=min(size_lt), max=max(size_lt))])

S1 <- merge(H1, CUT, by=c("dgp","model","N","ntime","method"))
PW <- S1[, .(nH1=.N, adj_lt = 100*mean(p<q05), adj_le = 100*mean(p<=q05),
             raw = 100*mean(p<0.05)), by=.(dgp,method)]
cat("\n### POWER: DEC-CRIT strict vs the tie-inclusive repair\n")
print(PW[order(dgp,-adj_lt), .(dgp,method,nH1, raw=round(raw,3),
        adj_lt=round(adj_lt,3), adj_le=round(adj_le,3), d=round(adj_le-adj_lt,3))])
cat("\n### spread and glmer-vs-FG gap under each convention\n")
print(PW[, .(min_lt=round(min(adj_lt),3), max_lt=round(max(adj_lt),3),
             min_le=round(min(adj_le),3), max_le=round(max(adj_le),3)), by=dgp])
print(PW[method %in% c("glmer","gee_fg"), .(gap_lt=round(adj_lt[method=="glmer"]-adj_lt[method=="gee_fg"],4),
         gap_le=round(adj_le[method=="glmer"]-adj_le[method=="gee_fg"],4)), by=dgp])
fwrite(SV, file.path(OUT,"R2_size_validation.csv"))
fwrite(PW, file.path(OUT,"R2_power_tieconvention.csv"))
fwrite(CUT, file.path(OUT,"R2_cut_percell.csv"))
cat("\nDONE_F\n")
