.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(geepack); library(data.table)})
setwd(path.expand("~/MMM"))
TERM <- "regionposterior"
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region=factor(case_when(grepl("A$",as.character(section))~"anterior",
                                 grepl("P$",as.character(section))~"posterior",
                                 TRUE~"lateral"), levels=c("anterior","lateral","posterior")),
         id=droplevels(factor(id)))
draws <- readRDS("rev2/T3_draws.rds"); idx <- which(sapply(draws,function(z) z$N==10))
A <- fread("rev2/T9_authors_path.csv")[method=="gee_fg"]
A[, usable := is.finite(estimate) & is.finite(se) & se>0 & se<10 &
     !((is.finite(estimate)&abs(estimate)>10)|(is.finite(se)&se>10))]
B <- fread("rev2/T3_subsample_rerun.csv")[spec=="pieno" & N==10 & method=="GEE-FG"]
B[, usable := is.finite(estimate) & is.finite(se) & se>0 & se<10 &
     !((is.finite(estimate)&abs(estimate)>10)|(is.finite(se)&se>10))]
mm <- merge(A[,.(rep,u_aut=usable,se_aut=se)], B[,.(rep,u_mio=usable,se_mio=se,status)], by="rep")
seven <- mm[u_aut==FALSE & u_mio==TRUE]$rep
cat("repliche discordanti (autori non usable, mio usable):", paste(seven, collapse=", "), "\n")
print(mm[rep %in% seven])
f <- lus_severe ~ region + disease + invasive + age_weeks
cat("\n=== per ciascuna: il base geeglm converge con e senza scale.fix? ===\n")
for (rp in seven) {
  z <- draws[[idx[which(sapply(idx, function(i) draws[[i]]$rep)==rp)]]]
  if (is.null(z)) { z <- draws[[ idx[ which(sapply(draws[idx], function(w) w$rep)==rp) ] ]] }
  dat <- d[d$id %in% z$sel, ] %>% arrange(id) %>%
    mutate(id=droplevels(factor(id)), across(c(disease,invasive), droplevels))
  D <- as.data.frame(dat)
  g_fix <- try(geeglm(f, id=id, data=D, family=binomial, corstr="exchangeable", scale.fix=TRUE), silent=TRUE)
  g_free<- try(geeglm(f, id=id, data=D, family=binomial, corstr="exchangeable"), silent=TRUE)
  cat(sprintf("rep %3d: scale.fix=TRUE %s | scale.fix libero %s | livelli invasive %d, disease %d\n",
      rp, if (inherits(g_fix,"try-error")) "ERRORE" else "ok",
      if (inherits(g_free,"try-error")) "ERRORE" else "ok",
      length(unique(as.character(dat$invasive))), length(unique(as.character(dat$disease)))))
  if (inherits(g_fix,"try-error")) cat("   messaggio:", substr(as.character(g_fix),1,110), "\n")
}
cat("\nDONE_T10\n")
