
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
S <- fread(path.expand("~/MMM/rev2/REV_subsample_all.csv"))
cat("colonne:", paste(names(S),collapse=","), " righe:", nrow(S), "\n")
cat("metodi:", paste(unique(S$method), collapse=","), "\n")
cat("N:", paste(sort(unique(S$N)), collapse=","), " unit:", paste(unique(S$unit),collapse=","), "\n")
cat("repliche per (method,N):\n"); print(S[, .N, by=.(method,N)][order(method,N)][1:12])
## rapporto SE/SD per N, come lo riporta §2.5/§3.16: media(se)/sd(estimate)
U <- S[!separated & !se_missing & is.finite(estimate) & is.finite(se)]
R <- U[, .(n=.N, sd_emp=sd(estimate), se_mean=mean(se)), by=.(method,N,unit)]
R[, ratio := se_mean/sd_emp]
cat("\n=== rapporto SE/SD per metodo e N (unit=patient) ===\n")
print(R[unit=="patient"][order(method,N)], digits=5)
cat("\n=== GLMM: la serie contestata da R3-R01/R3-R02 (0,5205-0,6351?) ===\n")
g <- R[method=="GLMM" & unit=="patient"][order(N)]
print(g, digits=6)
cat("serie GLMM arrotondata:", paste(sprintf("%.4f", g$ratio), collapse=", "), "\n")
cat("range:", sprintf("%.4f - %.4f", min(g$ratio), max(g$ratio)), "\n")
cat("rapporto ultimo/primo - 1 = ", sprintf("%.1f%%", 100*(g$ratio[nrow(g)]/g$ratio[1]-1)), "\n")
cat("\n=== MMM (per R3-R04, 0,8758-1,0948) ===\n")
m <- R[method=="MMM" & unit=="patient"][order(N)]
print(m, digits=6)
if (nrow(m)) cat("serie MMM:", paste(sprintf("%.4f", m$ratio), collapse=", "), "\n")
cat("\n=== usable fits per metodo e N su 500 (per R3-D06/R3-S06) ===\n")
tot <- S[, .(tent=.N), by=.(method,N)]
ok  <- S[!separated & !se_missing & is.finite(se), .(usable=.N), by=.(method,N)]
mg <- merge(tot, ok, by=c("method","N"), all.x=TRUE)
mg[is.na(usable), usable:=0][, pct := round(100*usable/tent,2)]
print(mg[order(method,N)], digits=5)
cat("\nDONE_JOB5\n")
