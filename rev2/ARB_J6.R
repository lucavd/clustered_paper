
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table)})
ok <- suppressWarnings(suppressMessages(require(mmmb, quietly=TRUE))) ||
      suppressWarnings(suppressMessages(require(MMLB, quietly=TRUE)))
cat("pacchetto MMM disponibile:", ok, "\n")
cat("candidati installati:", paste(intersect(c("MMLB","mmmb","mm","geepack","lme4"), rownames(installed.packages())), collapse=", "), "\n")
D <- readRDS(path.expand("~/MMM/ARTICOLO 2026/LURID_bis/lurid_long.rds"))
D <- as.data.table(D)
cat("\ncoorte: righe", nrow(D), " colonne:", paste(names(D), collapse=","), "\n")
cat("pazienti unici:", uniqueN(D$id), "\n")
for (v in intersect(c("region","disease","invasive","age_weeks","lus_severe","pt","hemi"), names(D)))
  cat(sprintf("  %-12s classe %s, livelli/NA: %s / %d\n", v, class(D[[v]])[1],
      paste(head(unique(D[[v]]),4), collapse="|"), sum(is.na(D[[v]]))))
## ---- V2-10: il modello ridotto cambia l'estimando? GEE e GLMM sono verificabili qui ----
suppressMessages({library(geepack); library(lme4)})
D2 <- D[!is.na(lus_severe) & !is.na(region)]
cat("\nrighe usate:", nrow(D2), "\n")
setorder(D2, id)
f_full <- lus_severe ~ region + disease + invasive + age_weeks
f_red  <- lus_severe ~ region
getc <- function(m, nm="regionposterior") {
  cf <- try(summary(m)$coefficients, silent=TRUE); if (inherits(cf,"try-error")) return(c(NA,NA))
  rn <- rownames(cf); j <- grep("^region", rn)
  if (!length(j)) return(c(NA,NA)); c(cf[j[1],1], cf[j[1],2])
}
for (lab in c("full","reduced")) {
  f <- if (lab=="full") f_full else f_red
  g <- try(geeglm(f, id=D2$id, data=D2, family=binomial, corstr="exchangeable"), silent=TRUE)
  gl <- try(glmer(update(f, . ~ . + (1|id)), data=D2, family=binomial), silent=TRUE)
  cat(sprintf("\n--- %s ---\n", lab))
  if (!inherits(g,"try-error")) { v<-getc(g); cat(sprintf("  GEE  est %.4f  SE %.4f\n", v[1], v[2])) } else cat("  GEE  FALLITO\n")
  if (!inherits(gl,"try-error")) { v<-getc(gl); cat(sprintf("  GLMM est %.4f  SE %.4f\n", v[1], v[2])) } else cat("  GLMM FALLITO\n")
}
cat("\nDONE_JOB6\n")
