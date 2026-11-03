
suppressMessages(library(data.table))
A <- rbind(fread("tables/V4_S7d_dp_rates.csv"), fread("tables/V4_S7e_dp_rates.csv"))
A <- A[strategy %in% c("gee_ar1_fg","mmm_cc","study_mmm") & ref=="t"][strategy=="study_mmm", strategy := "mmm_cc"]
A <- unique(A, by=c("dgp","scen","N","ntime","strategy"))
cat("rows", nrow(A), " points per strategy:"); print(A[, .N, strategy])
W <- dcast(A, dgp+scen+N+ntime ~ strategy, value.var=c("rate","mc","assess","pu","upow","hi","lo"))
cat("points", nrow(W), "\n")
both <- W[assess_gee_ar1_fg & assess_mmm_cc]
cat("both assessable", nrow(both), "\n")
summ <- function(D, lab) {
  out <- rbindlist(lapply(c("gee_ar1_fg","mmm_cc"), function(s) {
    r <- D[[paste0("rate_",s)]]; a <- D[[paste0("assess_",s)]]
    data.table(set=lab, strategy=s, n_pts=nrow(D), n_assess=sum(a),
      mad=mean(abs(r[a]-5)), signed=mean(r[a]-5), worst=max(abs(r[a]-5)),
      worst_above=max(r[a]-5), n_hi=sum(D[[paste0("hi_",s)]]), n_lo=sum(D[[paste0("lo_",s)]]),
      pu=mean(D[[paste0("pu_",s)]]), upow=mean(D[[paste0("upow_",s)]]))
  }))
  out
}
R <- rbind(summ(both,"all, both assessable"),
  summ(both[ntime<=3],"T 2-3"), summ(both[ntime>=4],"T >= 4"),
  summ(both[N==10],"N 10"), summ(both[N>=15 & N<=25],"N 15-25"), summ(both[N>25],"N > 25"),
  summ(both[dgp=="Copula-AR1"],"copula"), summ(both[dgp!="Copula-AR1"],"three mechanisms"))
R[, `:=`(mad=round(mad,2), signed=round(signed,2), worst=round(worst,2), worst_above=round(worst_above,2), pu=round(pu,2), upow=round(upow,2))]
fwrite(R, "tables/V5_S7f_h2h.csv"); print(R)
# pointwise: which is closer to 5, among points where both assessable
both[, d := abs(rate_gee_ar1_fg-5) - abs(rate_mmm_cc-5)]
both[, se := sqrt(mc_gee_ar1_fg^2 + mc_mmm_cc^2)]
cat("MMM closer (any):", sum(both$d>0), " GEE closer:", sum(both$d<0), "\n")
cat("MMM clearly closer (>2se):", sum(both$d > 2*both$se), " GEE clearly closer:", sum(both$d < -2*both$se), "\n")
cat("pu diff (GEE-MMM) mean:", mean(W$pu_gee_ar1_fg - W$pu_mmm_cc), " T2-3:", W[ntime<=3, mean(pu_gee_ar1_fg - pu_mmm_cc)], " T>=4:", W[ntime>=4, mean(pu_gee_ar1_fg - pu_mmm_cc)], "\n")
cat("upow diff (GEE-MMM) T2-3:", W[ntime<=3, mean(upow_gee_ar1_fg - upow_mmm_cc)], " T>=4:", W[ntime>=4, mean(upow_gee_ar1_fg - upow_mmm_cc)], "\n")
