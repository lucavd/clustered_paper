
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
rd <- function(x){ z <- as.data.table(readRDS(x))
  z[term=="time:trt", .(estimate,se,p,method,N,ntime,model,label,converged,dgp)] }
fs <- list.files(RES, pattern="[.]rds$", full.names=TRUE)
k <- sub("_c[0-9]+[.]rds$","",basename(fs)); fs <- fs[k %in% names(which(table(k)==10))]
cat("chunk files used:", length(fs), "\n")
D <- rbindlist(mclapply(fs, rd, mc.cores=40)); D[, model := as.integer(model)]
keep <- function(d) rbind(
  d[method %in% WALD & converged & is.finite(estimate) & is.finite(se) &
    abs(estimate)<10 & se>0 & se<10 & is.finite(p)],
  d[method=="glmer_lrt" & converged & is.finite(p)])
H0 <- keep(D[label=="manual"]); H1 <- keep(D[label=="mmlb"]); rm(D); gc(verbose=FALSE)
H0[, z := ifelse(method=="glmer_lrt", NA_real_, abs(estimate/se))]
H1[, z := ifelse(method=="glmer_lrt", NA_real_, abs(estimate/se))]

## ---- VARIANT A critical values, per cell -----------------------------------
## Wald methods: 0.95 quantile (type 7) of |estimate/se| on that cell's H0 draws.
## glmer_lrt has no Wald SE, so A is undefined on the |z| scale for it; the
## 0.05 quantile of its p is used instead and that is flagged in the output.
CUT <- H0[, {
  if (all(is.na(z))) .(q95 = NA_real_, q05p = quantile(p, 0.05, type=7), nH0=.N,
                       ties = sum(p == quantile(p, 0.05, type=7)), scale="p")
  else .(q95 = quantile(z, 0.95, type=7), q05p = NA_real_, nH0=.N,
         ties = sum(z == quantile(z, 0.95, type=7)), scale="z")
}, by=.(dgp,model,N,ntime,method)]
cat("\ncells per method, and cells whose critical value is tied/degenerate:\n")
print(CUT[, .(ncell=.N, scale=scale[1], cells_ties_gt1 = sum(ties>1),
              q95_min=round(min(q95),4), q95_max=round(max(q95),4)), by=method][order(method)])

## ---- VALIDATION: apply the cell's own critical value back to its H0 draws ---
S0 <- merge(H0, CUT, by=c("dgp","model","N","ntime","method"))
SV <- S0[, .(n=.N, size = 100*mean(if (scale[1]=="z") z > q95 else p < q05p)),
         by=.(dgp,method,scale)]
SV[, scale := NULL]
cat("\n### VARIANT A validation, all 24 (dgp, method) — empirical size (%)\n")
print(SV[order(dgp,method)], nrows=30)
cat("range over the 24:", sprintf("%.4f", min(SV$size)), "-", sprintf("%.4f", max(SV$size)), "\n")
cat("all 24 within 5.00-5.03:", all(SV$size >= 5.00 & SV$size <= 5.03), "\n")
cat("number outside 5.00-5.03:", sum(SV$size < 5.00 | SV$size > 5.03), "\n")

## ---- POWER under variant A --------------------------------------------------
S1 <- merge(H1, CUT, by=c("dgp","model","N","ntime","method"))
PW <- S1[, .(nH1=.N, raw = 100*mean(p < 0.05),
             adjA = 100*mean(if (scale[1]=="z") z > q95 else p < q05p)),
         by=.(dgp,method)]
cat("\n### TABLE 3.11 under variant A\n")
print(PW[order(dgp,-adjA), .(dgp,method,nH1, raw=round(raw,3), adjA=round(adjA,3))], nrows=30)
cat("\n### spread of size-adjusted power (variant A), per mechanism\n")
print(PW[, .(min=round(min(adjA),3), max=round(max(adjA),3),
             spread=round(max(adjA)-min(adjA),3),
             argmin=method[which.min(adjA)], argmax=method[which.max(adjA)]), by=dgp])
cat("\n### glmer vs gee_fg (variant A), per mechanism\n")
print(PW[method %in% c("glmer","gee_fg"), .(
  glmer = round(adjA[method=="glmer"],4), gee_fg = round(adjA[method=="gee_fg"],4),
  gap = round(adjA[method=="glmer"]-adjA[method=="gee_fg"],4)), by=dgp])

## per-N detail, for the record
PWN <- S1[, .(nH1=.N, raw=100*mean(p<0.05),
              adjA=100*mean(if (scale[1]=="z") z > q95 else p < q05p)), by=.(dgp,method,N)]
fwrite(PW,  file.path(OUT,"R2_311_variantA_pooled.csv"))
fwrite(SV,  file.path(OUT,"R2_311_variantA_validation.csv"))
fwrite(PWN, file.path(OUT,"R2_311_variantA_byN.csv"))
fwrite(CUT, file.path(OUT,"R2_311_variantA_cut.csv"))
cat("\nDONE_G\n")
