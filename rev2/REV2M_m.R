# REV2M_m.R -- AR(1) conditional truth under the ACTUAL vs the SPEC re.formula,
# and the pooled coverage/bias figures §3.13 quotes.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(fastGHQuad); library(MASS); library(data.table)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_ar1_data.R"), envir=globalenv())
sys.source(file.path(RD,"ar1/generation_models_ar1.R"), envir=globalenv())

cat("########## AR(1): generator re.formula vs the one F_truth_ar1.csv used\n")
cat("generator (generation_models_ar1.R):\n")
print(grep("\\.ar1_generate\\(", readLines(file.path(RD,"ar1/generation_models_ar1.R")), value=TRUE))
cat("\ntruth_ar1_conditional.R SPEC:\n")
print(grep("re = ~", readLines(path.expand("~/MMM/R/postrun/truth_ar1_conditional.R")), value=TRUE))

D3 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
D5 <- matrix(c(1,.01,.01, .01,.001,0, .01,0,.01),3,byrow=TRUE)
TL <- c(2,3,4,6,8,10)
out <- list(); k <- 0
for (m in c(1,3,5)) for (T in TL) {
  set.seed(20260805L + 1000L*m + 1)
  ex <- if (m==1) "none" else if (m==3) "age" else "full"
  d <- .ar1_design(75, T, ex)
  tgt <- .ar1_target(m, "mmlb")
  mf <- if (m==5) ~ time*trt + dose else ~ time*trt
  spec <- list("1"=list(re=~1, D=1.0, re2=~1, D2=1.0),
               "3"=list(re=~1+time, D=D3, re2=~1+age, D2=D3),
               "5"=list(re=~1+time+dose, D=D5, re2=~1+age+sex, D2=D5))[[as.character(m)]]
  bA <- tryCatch(as.numeric(calibrate_beta_c_ar1(mf, spec$re, tgt, spec$D, AR1_SIGMA2, d)), error=function(e) NA)
  # the SPEC used in F_truth_ar1.csv: mean model WITHOUT dose for m=5, re on age/sex
  mf2 <- ~ time*trt
  tgt2 <- if (m==5) c(-2,0.1,0.2,0.2) else tgt
  bS <- tryCatch(as.numeric(calibrate_beta_c_ar1(mf2, spec$re2, tgt2, spec$D2, AR1_SIGMA2, d)), error=function(e) NA)
  nmA <- colnames(model.matrix(mf, model.frame(mf, d)))
  nmS <- colnames(model.matrix(mf2, model.frame(mf2, d)))
  k <- k+1
  out[[k]] <- data.frame(model=m, ntime=T,
    cond_actual = if (length(bA)>1) bA[match("time:trt", nmA)] else NA,
    cond_Ftruth = if (length(bS)>1) bS[match("time:trt", nmS)] else NA,
    hardcoded = if (m==4) 0.058 else 0.232)
}
A <- bind_rows(out)
A$pct_err_actual <- round(100*(A$hardcoded-A$cond_actual)/A$cond_actual, 2)
A$pct_err_Ftruth <- round(100*(A$hardcoded-A$cond_Ftruth)/A$cond_Ftruth, 2)
print(as.data.frame(A %>% mutate(across(c(cond_actual,cond_Ftruth), ~round(.x,5)))), row.names=FALSE)
cat(sprintf("\nrange with the ACTUAL structure: %.4f - %.4f ; with the F_truth structure: %.4f - %.4f\n",
  min(A$cond_actual,na.rm=TRUE), max(A$cond_actual,na.rm=TRUE),
  min(A$cond_Ftruth,na.rm=TRUE), max(A$cond_Ftruth,na.rm=TRUE)))
Fp <- fread(path.expand("~/MMM/tables/F_truth_ar1.csv"))
cat(sprintf("F_truth_ar1.csv as stored: cond_used range %.4f - %.4f (M4: %.4f - %.4f)\n",
  min(Fp[model!=4]$cond_used), max(Fp[model!=4]$cond_used),
  min(Fp[model==4]$cond_used), max(Fp[model==4]$cond_used)))
write.csv(A, path.expand("~/MMM/rev2/REV2M_ar1_truth.csv"), row.names=FALSE)

cat("\n########## §3.13: pooled coverage and bias of the mixed model, both targets\n")
o <- lapply(c("GenBinaryY","GLMMnative","GLMMAR1"), function(p)
  readRDS(path.expand(paste0("~/MMM/rev2/REV2M_cov_", p, ".rds"))))
names(o) <- c("GenBinaryY","GLMM-native","GLMM-AR1")
for (nm in names(o)) {
  x <- o[[nm]]$pooled[method=="glmer"]
  cat(sprintf("%-12s glmer: cov vs conditional 0.232 = %.2f ; vs marginal 0.200 = %.2f ; vs T-varying marginal = %.2f (n=%d)\n",
    nm, x$cov_flat, x$cov_flat_marg, x$cov_Tmarg, x$n))
}
tot <- rbindlist(lapply(names(o), function(nm) o[[nm]]$pooled[method=="glmer"][, dgp := nm]))
cat(sprintf("\npooled over the three mechanisms: conditional %.2f | marginal %.2f | T-varying marginal %.2f\n",
  sum(tot$cov_flat*tot$n)/sum(tot$n), sum(tot$cov_flat_marg*tot$n)/sum(tot$n),
  sum(tot$cov_Tmarg*tot$n)/sum(tot$n)))
rm2 <- rbindlist(lapply(names(o), function(nm) o[[nm]]$rmse[method=="glmer"][, dgp := nm]))
cat(sprintf("glmer bias: pooled vs conditional %.4f | vs T-varying marginal %.4f\n",
  sum(rm2$bias_flat*rm2$n)/sum(rm2$n), sum(rm2$bias_Tmarg*rm2$n)/sum(rm2$n)))
cat("\nper mechanism (bias, RMSE):\n")
print(as.data.frame(rm2[, .(dgp, n, bias_cond=round(bias_flat,4), bias_Tmarg=round(bias_Tmarg,4),
   rmse_cond=round(rmse_flat,4), rmse_Tmarg=round(rmse_Tmarg,4))]), row.names=FALSE)
cat("\n########## §3.14 RMSE: shrinkage vs the corrected GEEs, both targets\n")
for (nm in names(o)) {
  x <- o[[nm]]$rmse
  cat(sprintf("%-12s shrink %.4f/%.4f | gee_fg %.4f/%.4f | gee_md %.4f/%.4f | mmm %.4f/%.4f (flat/Tmarg)\n", nm,
    x[method=="shrink"]$rmse_flat, x[method=="shrink"]$rmse_Tmarg,
    x[method=="gee_fg"]$rmse_flat, x[method=="gee_fg"]$rmse_Tmarg,
    x[method=="gee_md"]$rmse_flat, x[method=="gee_md"]$rmse_Tmarg,
    x[method=="mmm"]$rmse_flat, x[method=="mmm"]$rmse_Tmarg))
  cat(sprintf("             glmer/mmm RMSE ratio: %.3f (flat) %.3f (Tmarg)\n",
    x[method=="glmer"]$rmse_flat/x[method=="mmm"]$rmse_flat,
    x[method=="glmer"]$rmse_Tmarg/x[method=="mmm"]$rmse_Tmarg))
  fgs <- x[method %in% c("gee_fg","gee_md","gee_kc")]
  cat(sprintf("             mmm vs best corrected GEE, relative gap: %.2f%% (flat) %.2f%% (Tmarg)\n",
    100*(x[method=="mmm"]$rmse_flat/min(fgs$rmse_flat)-1),
    100*(x[method=="mmm"]$rmse_Tmarg/min(fgs$rmse_Tmarg)-1)))
}
cat("\nDONE_M\n")
