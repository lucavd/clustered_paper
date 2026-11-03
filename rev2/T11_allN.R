.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(geepack); library(data.table); library(parallel)})
setwd(path.expand("~/MMM"))
source("R/analysis_smallsample.R")
TERM <- "regionposterior"
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region=factor(case_when(grepl("A$",as.character(section))~"anterior",
                                 grepl("P$",as.character(section))~"posterior",
                                 TRUE~"lateral"), levels=c("anterior","lateral","posterior")),
         id=droplevels(factor(id)))
draws <- readRDS("rev2/T3_draws.rds")
f <- lus_severe ~ region + disease + invasive + age_weeks
res <- mclapply(seq_along(draws), function(i) {
  z <- draws[[i]]
  dat <- d[d$id %in% z$sel, ] %>% arrange(id) %>%
    mutate(id=droplevels(factor(id)), across(c(disease,invasive), droplevels))
  ss <- try(analyze_gee_smallsample(as.data.frame(dat), f, N=z$N, ntime=10,
            label="lurid", term_of_interest=TERM), silent=TRUE)
  if (inherits(ss,"try-error") || is.null(ss))
    return(data.table(N=z$N, rep=z$rep, method=c("gee_md","gee_kc","gee_fg"),
                      estimate=NA_real_, se=NA_real_, returned=FALSE))
  data.table(N=z$N, rep=z$rep, method=ss$method, estimate=ss$estimate, se=ss$se, returned=TRUE)
}, mc.cores=40)
A <- rbindlist(res[sapply(res, is.data.table)])
fwrite(A, "rev2/T11_authors_allN.csv")
A[, sep := (is.finite(estimate)&abs(estimate)>10)|(is.finite(se)&se>10)]
A[, usable := is.finite(estimate)&is.finite(se)&se>0&se<10&!sep]
aut <- A[, .(usable_aut=sum(usable), pct_aut=round(100*sum(usable)/500,2)), by=.(method,N)]
aut[, method := c(gee_fg="GEE-FG", gee_kc="GEE-KC", gee_md="GEE-MD")[method]]
B <- fread("rev2/T3_subsample_rerun.csv")[spec=="pieno" & method %in% c("GEE-FG","GEE-KC","GEE-MD")]
B[, sep := (is.finite(estimate)&abs(estimate)>10)|(is.finite(se)&se>10)]
B[, usable := is.finite(estimate)&is.finite(se)&se>0&se<10&!sep]
cor <- B[, .(usable_cor=sum(usable), pct_cor=round(100*sum(usable)/500,2)), by=.(method,N)]
M <- merge(aut, cor, by=c("method","N"))
M[, `:=`(diff_pp=round(pct_cor-pct_aut,2), n_spurie=usable_cor-usable_aut)]
fwrite(M, "rev2/T11_confronto.csv")
cat("=== fit utilizzabili su 500: via autori (scale.fix=TRUE) vs via corretta ===\n")
print(M[order(method,N)])
cat("\n=== separazioni spurie per N (somma sui tre metodi) ===\n")
print(M[, .(spurie=sum(n_spurie)), by=N][order(N)])
cat("\nDONE_T11\n")
