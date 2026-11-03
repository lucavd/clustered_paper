.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table))
setwd(path.expand("~/MMM"))
A <- fread("rev2/T3_subsample_rerun.csv")
DG <- fread("rev2/T3_degeneracy.csv")
REF <- fread("rev2/T3_ref.csv")
DEN <- 500
cat("=== righe per spec/method ===\n"); print(dcast(A[, .N, by=.(spec,method)], method ~ spec, value.var="N"))
cat("\n=== stati per spec/method (conteggi) ===\n")
print(dcast(A[, .N, by=.(spec,method,status)], spec+method ~ status, value.var="N", fill=0))

cat("\n=== TASSO DI TIMEOUT (% su 500 estrazioni per N; pooled su N) ===\n")
print(dcast(A[, .(pct=round(100*mean(status=="timeout"),3)), by=.(spec,method)], method ~ spec, value.var="pct"))
cat("\ntimeout per N (spec=ridotto):\n")
print(dcast(A[spec=="ridotto", .(pct=round(100*mean(status=="timeout"),3)), by=.(method,N)], method ~ N, value.var="pct"))
cat("\ntimeout per N (spec=pieno):\n")
print(dcast(A[spec=="pieno", .(pct=round(100*mean(status=="timeout"),3)), by=.(method,N)], method ~ N, value.var="pct"))

## flag di usabilita' allineati a §2.4
A[, sep := is.finite(estimate) & abs(estimate) > 10 | (is.finite(se) & se > 10)]
A[, has_se := is.finite(se) & se > 0]
A[, returned := status %in% c("ok","bad_var","no_term") & (is.finite(estimate) | has_se)]
A[, returned := is.finite(estimate)]      # "fit restituito" = stima puntuale disponibile
A[, usable := returned & has_se & !sep & is.finite(estimate) & abs(estimate) <= 10 & se < 10]

tab <- A[, .(returned=sum(returned), no_se=sum(returned & !has_se),
             sep_with_se=sum(returned & has_se & sep), usable=sum(usable),
             conv_only=sum(returned)), by=.(spec,method,N)]
tab[, `:=`(pct_not_returned=100*(DEN-returned)/DEN, pct_no_se_500=100*no_se/DEN,
           pct_sep_500=100*sep_with_se/DEN, pct_usable_500=100*usable/DEN,
           pct_no_se_returned=100*no_se/returned, pct_usable_returned=100*usable/returned,
           pct_convonly_500=100*conv_only/DEN)]
tab[, somma := pct_not_returned+pct_no_se_500+pct_sep_500+pct_usable_500]
fwrite(tab, "rev2/T5_tab316.csv")
for (sp in unique(tab$spec)) {
  cat("\n=== TABELLA 3.16, spec =", sp, ", N=10 ===\n")
  print(tab[spec==sp & N==10, .(method, not_returned=round(pct_not_returned,2),
     no_se_500=round(pct_no_se_500,2), sep_500=round(pct_sep_500,2),
     usable_500=round(pct_usable_500,2), somma=round(somma,4),
     no_se_of_returned=round(pct_no_se_returned,2), usable_of_returned=round(pct_usable_returned,2),
     conv_only_500=round(pct_convonly_500,2))][order(-usable_500)], digits=5)
}
cat("\n=== usable su 500 per tutti gli N ===\n")
print(dcast(tab, spec+method ~ N, value.var="pct_usable_500"), digits=5)
cat("\n=== SE presente su 500 (R3-S06) ===\n")
print(dcast(A[, .(pct=round(100*sum(has_se)/DEN,2)), by=.(spec,method,N)], spec+method ~ N, value.var="pct"))
cat("\n=== modo di fallimento a N=10: degenerazione a monte vs SE assente ===\n")
dg10 <- DG[N==10]
cat("estrazioni a N=10 con invasive a un solo livello:", sum(dg10$lev_inv==1), "=",
    round(100*mean(dg10$lev_inv==1),2), "%\n")
cat("estrazioni con disease < 4 livelli:", sum(dg10$lev_dis<4), "=", round(100*mean(dg10$lev_dis<4),2), "%\n")
M <- merge(A[N==10][, !c("lev_inv","lev_dis"), with=FALSE], dg10[, .(rep, lev_inv, lev_dis)], by="rep")
cat("\nusable a N=10 condizionato alla degenerazione di invasive (spec=pieno):\n")
print(M[spec=="pieno", .(n=.N, usable=sum(usable), pct=round(100*mean(usable),2)), by=.(method, lev_inv)][order(method,lev_inv)])
cat("\nusable a N=10 (spec=ridotto), per lev_inv (il modello ridotto non usa invasive):\n")
print(M[spec=="ridotto", .(n=.N, usable=sum(usable), pct=round(100*mean(usable),2)), by=.(method, lev_inv)][order(method,lev_inv)])

## ---- serie SE/SD ----
U <- A[usable == TRUE]
R <- U[, .(n=.N, mean_se=mean(se), sd_est=sd(estimate)), by=.(spec,method,N)][, ratio := mean_se/sd_est]
fwrite(R, "rev2/T5_seratio.csv")
cat("\n=== rapporto mean(se)/sd(estimate) ===\n")
for (sp in unique(R$spec)) for (mm in c("GLMM","MMM","GEE","GEE-FG","GEE-KC","GEE-MD")) {
  g <- R[spec==sp & method==mm][order(N)]
  if (!nrow(g)) next
  cat(sprintf("%-12s %-7s: %s | range %.4f-%.4f | %+.1f%% | monotona:%s | n=%s\n", sp, mm,
      paste(sprintf("%.4f", g$ratio), collapse=" / "), min(g$ratio), max(g$ratio),
      100*(g$ratio[nrow(g)]/g$ratio[1]-1), all(diff(g$ratio)>0),
      paste(g$n, collapse=",")))
}
## ---- bias ----
cat("\n=== riferimenti a coorte piena ===\n"); print(REF, digits=6)
B <- U[, .(n=.N, mean_est=mean(estimate), sd_est=sd(estimate)), by=.(spec,method,N)]
B <- merge(B, REF[, .(spec, method, ref=est)], by=c("spec","method"), all.x=TRUE)
refm <- REF[method=="MMM", .(spec, ref_mmm=est)]
B <- merge(B, refm, by="spec", all.x=TRUE)
B[is.na(ref), ref := ref_mmm]     # correzioni: riferimento = GEE marginale
gref <- REF[method=="GEE", .(spec, ref_gee=est)]
B <- merge(B, gref, by="spec", all.x=TRUE)
B[method %in% c("GEE-FG","GEE-KC","GEE-MD"), ref := ref_gee]
B[, bias := mean_est - ref]
fwrite(B, "rev2/T5_bias.csv")
cat("\n=== bias per spec/metodo/N ===\n")
print(dcast(B, spec+method ~ N, value.var="bias"), digits=4)
cat("\nserie MMM (spec=pieno):", paste(sprintf("%+.4f", B[spec=="pieno" & method=="MMM"][order(N)]$bias), collapse=" / "), "\n")
cat("serie MMM (spec=ridotto):", paste(sprintf("%+.4f", B[spec=="ridotto" & method=="MMM"][order(N)]$bias), collapse=" / "), "\n")
cat("\n=== MMM: ha restituito coefficienti? ===\n")
print(A[method=="MMM", .(n=.N, ok=sum(status=="ok"), finite_est=sum(is.finite(estimate)),
        finite_se=sum(is.finite(se)), timeout=sum(status=="timeout"),
        error=sum(status=="error")), by=.(spec,N)][order(spec,N)])
cat("\nDONE_T5\n")
