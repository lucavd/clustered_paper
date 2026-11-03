.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table)})
setwd(path.expand("~/MMM"))
S <- fread("rev2/REV_subsample_all.csv")
cat("=== REV_subsample_all.csv ===\n")
cat("righe:", nrow(S), " metodi:", paste(sort(unique(S$method)),collapse=","), "\n")
cat("unit:", paste(unique(S$unit),collapse=","), " N:", paste(sort(unique(S$N)),collapse=","), "\n")
S <- S[unit=="patient"]
cat("righe unit=patient:", nrow(S), "\n")

## quali repliche esistono affatto (= modello identificabile e almeno un fit)
pres <- S[, .(nmeth=uniqueN(method)), by=.(N,rep)]
cat("\nrepliche presenti per N (su 500):\n"); print(pres[, .(reps_presenti=.N, pct=round(100*.N/500,2)), by=N][order(N)])

## confronto con la degenerazione replicata dal seme (T3_degeneracy.csv)
DG <- fread("rev2/T3_degeneracy.csv")
M <- merge(DG, pres[, .(N, rep, presente=TRUE)], by=c("N","rep"), all.x=TRUE)
M[is.na(presente), presente:=FALSE]
cat("\n=== il flusso RNG del replay coincide con quello degli autori? ===\n")
print(M[, .(n=.N, degen=sum(lev_inv==1), presenti=sum(presente),
            accordo=sum((lev_inv==2)==presente), pct_accordo=round(100*mean((lev_inv==2)==presente),2)),
        by=N][order(N)])
cat("\ntabella di contingenza a N=10 (lev_inv x presente):\n"); print(table(M[N==10]$lev_inv, M[N==10]$presente))

## ---- tabella 3.16 sul MODELLO PIENO (dati esistenti) ----
den <- 500
tab <- S[, {
  ret <- .N
  us  <- sum(!separated & !se_missing & is.finite(se) & is.finite(estimate) & se > 0)
  conv_only <- sum(is.finite(estimate))
  .(returned=ret, pct_returned=100*ret/den,
    no_se=sum(se_missing | !is.finite(se)), pct_no_se_of_returned=100*sum(se_missing|!is.finite(se))/ret,
    pct_no_se_of_500=100*sum(se_missing|!is.finite(se))/den,
    sep=sum(separated), usable=us, pct_usable_500=100*us/den, pct_usable_returned=100*us/ret,
    pct_convonly_500=100*conv_only/den)
}, by=.(method,N)]
fwrite(tab, "rev2/T4_tab316_pieno.csv")
cat("\n=== §3.16 modello PIENO, N=10 ===\n")
print(tab[N==10][order(-pct_usable_500)], digits=5)
cat("\nsomma delle tre quote a N=10 (non restituito + restituito senza SE + usable):\n")
print(tab[N==10, .(method, non_ret=round(100-pct_returned,2), no_se=round(pct_no_se_of_500,2),
                   usable=round(pct_usable_500,2),
                   somma=round((100-pct_returned)+pct_no_se_of_500+pct_usable_500,4))])
cat("\n=== percentuale di sottocampioni (su 500) con SE presente, N=10 (R3-S06) ===\n")
print(S[N==10, .(se_presente_su500=round(100*sum(!se_missing & is.finite(se))/500,2),
                 restituiti=.N), by=method][order(-se_presente_su500)])
cat("\n=== tutti gli N ===\n"); print(tab[order(method,N)], digits=5)

## ---- serie SE/SD (modello pieno) ----
U <- S[!separated & !se_missing & is.finite(se) & is.finite(estimate)]
R <- U[, .(n=.N, mean_se=mean(se), sd_est=sd(estimate)), by=.(method,N)][, ratio := mean_se/sd_est]
fwrite(R, "rev2/T4_seratio_pieno.csv")
cat("\n=== rapporto mean(se)/sd(estimate), modello PIENO ===\n")
print(R[order(method,N)], digits=6)
for (mm in c("GLMM","MMM","GEE")) {
  g <- R[method==mm][order(N)]
  cat(sprintf("%s: %s | range %.4f-%.4f | +%.1f%% da N=10 a N=50 | monotona: %s\n", mm,
      paste(sprintf("%.4f", g$ratio), collapse=" / "), min(g$ratio), max(g$ratio),
      100*(g$ratio[nrow(g)]/g$ratio[1]-1), all(diff(g$ratio) > 0)))
}
## ---- bias contro riferimento a coorte piena ----
REF <- fread("rev2/T2_fullcohort.csv")
ref_mmm <- REF[tag=="A:pieno" & method=="MMM"]$est
ref_gee <- REF[tag=="A:pieno" & method=="GEE"]$est
ref_glmm<- REF[tag=="A:pieno" & method=="GLMM"]$est
cat(sprintf("\nriferimenti coorte piena (modello pieno): MMM=%.6f GEE=%.6f GLMM=%.6f\n", ref_mmm, ref_gee, ref_glmm))
B <- U[, .(n=.N, mean_est=mean(estimate)), by=.(method,N)]
B[, ref := ifelse(method=="GLMM", ref_glmm, ifelse(method=="MMM", ref_mmm, ref_gee))]
B[, bias := mean_est - ref]
B[, bias_vs_mmm := mean_est - ref_mmm]
fwrite(B, "rev2/T4_bias_pieno.csv")
cat("\n=== bias per metodo e N (modello PIENO) ===\n"); print(B[order(method,N)], digits=5)
cat("\nserie MMM:", paste(sprintf("%+.4f", B[method=="MMM"][order(N)]$bias), collapse=" / "), "\n")
cat("\nqualche metodo parte da ~0.16 a N=10?\n")
print(B[N==10, .(method, bias=round(bias,4), bias_vs_ref_mmm=round(bias_vs_mmm,4))])

## ---- D2-33: cella di simulazione corrispondente (T=10, N=10) ----
CE <- fread("rev2/REV_cells.csv")
CE <- CE[set=="own"]
K <- CE[N==10 & ntime==10]
cat("\n=== REV_cells.csv, N=10 ntime=10, set=own ===\n")
cat("righe:", nrow(K), " label:", paste(unique(K$label),collapse=","),
    " dgp:", paste(unique(K$dgp),collapse=","), " model:", paste(unique(K$model),collapse=","), "\n")
agg <- K[, .(n_attempt=sum(n_attempt), n_usable=sum(n_usable), n_no_se=sum(n_no_se),
             n_present=sum(n_present), n_conv=sum(n_conv), cells=.N), by=.(method,label)]
agg[, `:=`(pct_usable=100*n_usable/n_attempt, pct_no_se_of_present=100*n_no_se/n_present,
           pct_no_se_of_attempt=100*n_no_se/n_attempt)]
fwrite(agg, "rev2/T4_simcell.csv")
print(agg[order(label,method)], digits=5)
cat("\n=== solo H0 (label=manual), pooled dgp/model ===\n")
print(agg[label=="manual", .(method, n_attempt, pct_usable=round(pct_usable,2),
      pct_no_se_of_present=round(pct_no_se_of_present,4), pct_no_se_of_attempt=round(pct_no_se_of_attempt,4))][order(-pct_usable)])
cat("\nrange usabilita' H0:", sprintf("%.2f - %.2f", min(agg[label=="manual"]$pct_usable), max(agg[label=="manual"]$pct_usable)), "\n")
cat("\nDONE_T4\n")
