
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(geepack); library(lme4); library(MMLB)})
D <- as.data.table(readRDS(path.expand("~/MMM/ARTICOLO 2026/LURID_bis/lurid_long.rds")))
D <- D[!is.na(lus_severe)]
## region = posteriore vs anteriore, derivata dal suffisso di section (A/P)
D[, sec := as.character(section)]
D[, region := fifelse(grepl("P$", sec), "posterior", fifelse(grepl("A$", sec), "anterior", NA_character_))]
cat("region derivata da section:\n"); print(D[, .N, by=.(sec, region)][order(sec)])
D2 <- D[!is.na(region)]
D2[, region := factor(region, levels=c("anterior","posterior"))]
cat("\nrighe con region:", nrow(D2), " pazienti:", uniqueN(D2$id), "\n")
setorder(D2, id)
pick <- function(m, pat="^regionposterior") {
  cf <- try(summary(m)$coefficients, silent=TRUE); if (inherits(cf,"try-error")) return(c(NA,NA))
  j <- grep(pat, rownames(cf)); if (!length(j)) return(c(NA,NA)); c(cf[j[1],1], cf[j[1],2])
}
f_full <- lus_severe ~ region + disease + invasive + age_weeks
f_red  <- lus_severe ~ region
f_age  <- lus_severe ~ region + age_weeks
res <- list()
for (lab in c("full","reduced","region+age")) {
  f <- switch(lab, full=f_full, reduced=f_red, `region+age`=f_age)
  g  <- try(geeglm(f, id=D2$id, data=D2, family=binomial, corstr="exchangeable"), silent=TRUE)
  gl <- try(glmer(update(f, . ~ . + (1|id)), data=D2, family=binomial), silent=TRUE)
  mm <- try(mm(f, lv.formula = ~1, id=D2$id, data=as.data.frame(D2), verbose=FALSE), silent=TRUE)
  vg <- if (!inherits(g,"try-error")) pick(g) else c(NA,NA)
  vl <- if (!inherits(gl,"try-error")) pick(gl) else c(NA,NA)
  vm <- if (!inherits(mm,"try-error")) {
     cf <- try(mm$beta, silent=TRUE)
     se <- try(sqrt(diag(mm$mod.cov)), silent=TRUE)
     sr <- try(sqrt(diag(mm$rob.cov)), silent=TRUE)
     j <- grep("region", names(cf))
     if (length(j)) c(cf[j[1]], if(!inherits(se,"try-error")) se[j[1]] else NA,
                      if(!inherits(sr,"try-error")) sr[j[1]] else NA) else c(NA,NA,NA)
   } else c(NA,NA,NA)
  cat(sprintf("\n=== %s ===\n", lab))
  cat(sprintf("  GEE  est %.4f  SE %.4f\n", vg[1], vg[2]))
  cat(sprintf("  GLMM est %.4f  SE %.4f\n", vl[1], vl[2]))
  cat(sprintf("  MMM  est %.4f  SE(model) %.4f  SE(rob) %.4f\n", vm[1], vm[2], vm[3]))
  res[[lab]] <- list(gee=vg, glmm=vl, mmm=vm)
}
cat("\n=== scarto pieno -> ridotto ===\n")
for (meth in c("gee","glmm","mmm")) {
  a <- res[["full"]][[meth]][1]; b <- res[["reduced"]][[meth]][1]
  se <- res[["full"]][[meth]][2]
  cat(sprintf("  %-5s %.4f -> %.4f  scarto %+.4f  (%.1f%% della stima, %.2f SE)\n",
      meth, a, b, b-a, 100*(b-a)/a, (b-a)/se))
}
cat("\nDONE_JOB8\n")
