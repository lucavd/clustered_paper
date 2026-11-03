
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(geepack); library(lme4)})
D <- as.data.table(readRDS(path.expand("~/MMM/ARTICOLO 2026/LURID_bis/lurid_long.rds")))
cat("livelli di section:", paste(levels(factor(D$section)), collapse="|"), "\n")
cat("livelli di side:", paste(levels(factor(D$side)), collapse="|"), "\n")
D <- D[!is.na(lus_severe)]
cat("righe con esito:", nrow(D), " pazienti:", uniqueN(D$id), "\n")
D[, section := factor(section)]
D[, zona := section]
setorder(D, id)
grabc <- function(m) {
  cf <- try(summary(m)$coefficients, silent=TRUE); if (inherits(cf,"try-error")) return(NULL)
  rn <- rownames(cf); j <- grep("^(section|zona)", rn)
  if (!length(j)) return(NULL)
  data.table(term=rn[j], est=cf[j,1], se=cf[j,2])
}
f_full <- lus_severe ~ zona + disease + invasive + age_weeks
f_red  <- lus_severe ~ zona
for (lab in c("pieno","ridotto")) {
  f <- if (lab=="pieno") f_full else f_red
  cat(sprintf("\n=========== modello %s: %s ===========\n", lab, deparse(f)))
  g  <- try(geeglm(f, id=D$id, data=D, family=binomial, corstr="exchangeable"), silent=TRUE)
  gl <- try(glmer(update(f, . ~ . + (1|id)), data=D, family=binomial), silent=TRUE)
  if (!inherits(g,"try-error"))  { cat("GEE:\n");  print(grabc(g), digits=5) }  else cat("GEE FALLITO\n")
  if (!inherits(gl,"try-error")) { cat("GLMM:\n"); print(grabc(gl), digits=5) } else cat("GLMM FALLITO\n")
}
cat("\nDONE_JOB7\n")
