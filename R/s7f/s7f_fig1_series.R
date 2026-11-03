
suppressMessages(library(data.table))
A <- rbind(fread("tables/V5_S7f_cp_S7d_dp_rates.csv"), fread("tables/V5_S7f_cp_S7e_dp_rates.csv"), fill=TRUE)
C <- fread("tables/V4_S7b_rule_dp.csv")
k <- c("dgp","scen","N","ntime")
s1 <- A[ref=="t" & strategy %in% c("gee_ar1_fg","study_mmm","mmm_cc","gee_exch_fg","gee_exch_md"), c(k,"strategy","rate","mc","pu","upow","nd","lo"), with=FALSE]
s1[strategy=="study_mmm", strategy:="mmm_cc"]
s2 <- rbind(C[ref=="n" & method=="glmer_lrt", c(k,"rate","mc","pu","upow","nd","lo"), with=FALSE][, strategy:="glmer_lrt"],
            A[ref=="n" & strategy=="glmer_lrt", c(k,"strategy","rate","mc","pu","upow","nd","lo"), with=FALSE], use.names=TRUE)
D <- unique(rbind(s1, s2, use.names=TRUE), by=c(k,"strategy"))
print(D[, .N, strategy]); stopifnot(D[, .N, strategy][, all(N==2376L)], D[is.na(rate), all(strategy=="glmer_lrt")], D[is.na(rate), .N]==282L)
NA_LRT <- D[is.na(rate), .N]; cat("LRT non-assessable", NA_LRT, "\n")
D[, ok := !is.na(rate)]
D[, Nb := fifelse(N==10, "10", fifelse(N<=25, "15-25", ">25"))]
D[, Tb := fifelse(ntime==2, "2", "3+")]   # T = 2 vs T >= 3: the calibration of the two strategies ranks differently only at T = 2
# MCSE approximations (independence across points; conservative for paired differences)
D[, se_pu := sqrt(pu*(100-pu)/nd)]; D[, se_up := sqrt(upow*(100-upow)/nd)]
agg <- function(x) x[, .(n_pts=.N, n_ok=sum(ok), mad=mean(abs(rate[ok]-5)), mad_se=sqrt(sum(mc[ok]^2))/sum(ok),
   rate=mean(rate[ok]), n_lo=sum(lo[ok]), pu=mean(pu), pu_se=sqrt(sum(se_pu^2))/.N,
   upow=mean(upow, na.rm=TRUE), upow_se=sqrt(sum(se_up^2, na.rm=TRUE))/.N), by=strategy]
G <- rbind(agg(D)[, `:=`(Nb="all", Tb="all")],
  D[, agg(.SD), by=.(Nb)][, Tb:="all"], D[, agg(.SD), by=.(Tb)][, Nb:="all"], D[, agg(.SD), by=.(Nb,Tb)], D[N>=15, agg(.SD), by=.(Tb)][, Nb:="15+"], use.names=TRUE)
fwrite(G, "tables/V5_S7f_h2h_grid.csv")
print(G[strategy %in% c("gee_ar1_fg","mmm_cc")][order(Nb,Tb,strategy)], digits=3)
# panel series: by exact T and by N
PT <- D[, .(n_ok=sum(ok), mad=mean(abs(rate[ok]-5)), mad_se=sqrt(sum(mc[ok]^2))/sum(ok), pu=mean(pu), pu_se=sqrt(sum(se_pu^2))/.N,
            upow=mean(upow,na.rm=TRUE), n=.N), by=.(strategy, ntime)]
PN <- D[, .(n_ok=sum(ok), mad=mean(abs(rate[ok]-5)), mad_se=sqrt(sum(mc[ok]^2))/sum(ok), rate=mean(rate[ok]), n=.N), by=.(strategy, N, Tb)]
fwrite(PT, "tables/V5_S7f_fig1_byT.csv"); fwrite(PN, "tables/V5_S7f_fig1_byN.csv")
print(dcast(PT, ntime~strategy, value.var="mad"), digits=3); print(dcast(PT, ntime~strategy, value.var="pu"), digits=3)
print(dcast(PN[strategy %in% c("gee_ar1_fg","mmm_cc")], N~strategy+Tb, value.var="mad"), digits=3)
# paired differences GEE - MMM
W <- dcast(D[strategy %in% c("gee_ar1_fg","mmm_cc")], dgp+scen+N+ntime ~ strategy, value.var=c("pu","upow","se_pu","se_up"))
cat("pu diff", W[, mean(pu_gee_ar1_fg-pu_mmm_cc)], "se", W[, sqrt(sum(se_pu_gee_ar1_fg^2+se_pu_mmm_cc^2))/.N], "\n")
cat("upow diff", W[, mean(upow_gee_ar1_fg-upow_mmm_cc)], "se", W[, sqrt(sum(se_up_gee_ar1_fg^2+se_up_mmm_cc^2))/.N], "\n")

DF <- data.table(measure=c("pu_gee_minus_mmm","upow_mmm_minus_gee"),
  est=c(W[, mean(pu_gee_ar1_fg-pu_mmm_cc)], W[, mean(upow_mmm_cc-upow_gee_ar1_fg)]),
  se=c(W[, sqrt(sum(se_pu_gee_ar1_fg^2+se_pu_mmm_cc^2))/.N], W[, sqrt(sum(se_up_gee_ar1_fg^2+se_up_mmm_cc^2))/.N]))
fwrite(DF, "tables/V5_S7f_h2h_diff.csv"); print(DF)
# consistency with Table 4 (same points, same flags)
T4 <- fread("tables/V5_S7f_table4.csv")
chk <- merge(G[Nb=="all" & Tb=="all" & strategy %in% c("gee_ar1_fg","mmm_cc"), .(strategy, n_lo, pu, upow)], T4[, .(strategy, lo_all, pu4=pu, upow4=upow)], by="strategy")
print(chk); stopifnot(chk[, all(n_lo==lo_all) && all(abs(pu-pu4)<1e-9) && all(abs(upow-upow4)<1e-9)])
cat("GF-h2h: consistent with Table 4\n")

# panels e-f: usable results and power by T, the two strategies and the sign-flip test (same sources as Table 4)
FL <- fread("tables/V5_S7f_cp_S7d2_flip_rules_dp.csv")[version == "with D-57 rules"]
f0 <- FL[label == "manual", .(dgp, scen, N, ntime, pu, nd)]
f1 <- FL[label == "mmlb", .(dgp, scen, N, ntime, upow = pu * rate / 100, nd1 = nd)]
fl <- merge(f0, f1, by = k)[, `:=`(strategy = "flip", se_pu = sqrt(pu*(100-pu)/nd), se_up = sqrt(upow*(100-upow)/nd1))]
stopifnot(nrow(fl) == 2376L)
E <- rbind(D[strategy %in% c("gee_ar1_fg","mmm_cc"), .(strategy, ntime, pu, upow, se_pu, se_up)], fl[, .(strategy, ntime, pu, upow, se_pu, se_up)])
EF <- E[, .(pu=mean(pu), pu_se=sqrt(sum(se_pu^2))/.N, upow=mean(upow), upow_se=sqrt(sum(se_up^2))/.N, n=.N), by=.(strategy, ntime)]
fwrite(EF, "tables/V5_S7f_fig1_ef.csv")
cat("flip overall pu", fl[, mean(pu)], " Table 4:", T4[strategy=="flip", pu], "\n"); stopifnot(abs(fl[, mean(pu)] - T4[strategy=="flip", pu]) < 1e-9)
cat("max MCSE e/f:", EF[, max(c(pu_se, upow_se))], "\n"); print(EF, digits=3)

DF2 <- rbind(fread("tables/V5_S7f_h2h_diff.csv"), data.table(measure="max_mcse_ef", est=EF[, max(c(pu_se, upow_se))], se=NA_real_))
fwrite(DF2, "tables/V5_S7f_h2h_diff.csv")

# panel a (D-61): design map, Nb x T (2, 3, 4+): which strategy is closer to 5% (difference beyond 2 MCSE), usable results
D[, Tm := fifelse(ntime == 2, "2", fifelse(ntime == 3, "3", "4+"))]
M <- D[strategy %in% c("gee_ar1_fg", "mmm_cc"), .(mad = mean(abs(rate - 5)), mad_se = sqrt(sum(mc^2)) / .N, pu = mean(pu), n = .N), by = .(strategy, Nb, Tm)]
M <- dcast(M, Nb + Tm + n ~ strategy, value.var = c("mad", "mad_se", "pu"))
M[, diff := mad_gee_ar1_fg - mad_mmm_cc][, diff_se := sqrt(mad_se_gee_ar1_fg^2 + mad_se_mmm_cc^2)]
M[, closer := fifelse(diff < -2 * diff_se, "GEE", fifelse(diff > 2 * diff_se, "MMM", "either"))]
fwrite(M, "tables/V5_S7f_fig1_map.csv"); print(M[, .(Nb, Tm, n, g = round(mad_gee_ar1_fg, 2), m = round(mad_mmm_cc, 2), z = round(diff / diff_se, 1), closer, pg = round(pu_gee_ar1_fg, 1), pm = round(pu_mmm_cc, 1))])
