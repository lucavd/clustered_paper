.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
setwd(path.expand("~/MMM"))
A <- fread("rev2/T3_subsample_rerun.csv"); REF <- fread("rev2/T3_ref.csv")
cat("=== T3_ref.csv (riferimenti a coorte piena) ===\n"); print(REF, digits=6)
A[, sep := (is.finite(estimate) & abs(estimate) > 10) | (is.finite(se) & se > 10)]
A[, has_se := is.finite(se) & se > 0]
A[, usable := is.finite(estimate) & has_se & !sep & se < 10]
U <- A[usable == TRUE]
RR <- REF[, .(spec, method, ref = estimate)]
r_mmm <- REF[method=="MMM", .(spec, ref_mmm=estimate)]
r_gee <- REF[method=="GEE", .(spec, ref_gee=estimate)]
B <- U[, .(n=.N, mean_est=mean(estimate)), by=.(spec,method,N)]
B <- merge(B, RR, by=c("spec","method"), all.x=TRUE)
B <- merge(B, r_mmm, by="spec", all.x=TRUE); B <- merge(B, r_gee, by="spec", all.x=TRUE)
B[method %in% c("GEE-FG","GEE-KC","GEE-MD"), ref := ref_gee]
B[, `:=`(bias = mean_est - ref, bias_vs_mmm = mean_est - ref_mmm)]
fwrite(B, "rev2/T7_bias.csv")
cat("\n=== bias (contro il riferimento del proprio metodo) ===\n")
print(dcast(B, spec+method ~ N, value.var="bias"), digits=4)
cat("\n=== bias contro il riferimento marginale MMM ===\n")
print(dcast(B, spec+method ~ N, value.var="bias_vs_mmm"), digits=4)
for (sp in c("pieno","ridotto")) for (mm in c("MMM","GEE","GLMM","GEE-FG","GEE-KC","GEE-MD")) {
  g <- B[spec==sp & method==mm][order(N)]
  cat(sprintf("%-8s %-7s bias: %s | n=%s\n", sp, mm,
      paste(sprintf("%+.4f", g$bias), collapse=" / "), paste(g$n, collapse=",")))
}
cat("\n=== chi arriva a ~0.16 a N=10 e ~0.015 a N=50? ===\n")
print(B[N %in% c(10,50), .(spec,method,N,bias=round(bias,4))][order(spec,method,N)])
cat("\nDONE_T7\n")
