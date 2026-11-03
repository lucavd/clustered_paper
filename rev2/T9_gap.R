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
idx <- which(sapply(draws, function(z) z$N==10))
f <- lus_severe ~ region + disease + invasive + age_weeks
## via ESATTA degli autori: analyze_gee_smallsample (base con scale.fix=TRUE, timeout 5s)
res <- mclapply(idx, function(i) {
  z <- draws[[i]]
  dat <- d[d$id %in% z$sel, ] %>% arrange(id) %>%
    mutate(id=droplevels(factor(id)), across(c(disease,invasive), droplevels))
  ss <- try(analyze_gee_smallsample(as.data.frame(dat), f, N=10, ntime=10,
            label="lurid", term_of_interest=TERM), silent=TRUE)
  if (inherits(ss,"try-error") || is.null(ss))
    return(data.table(rep=z$rep, method=c("gee_md","gee_kc","gee_fg"),
                      estimate=NA_real_, se=NA_real_, returned=FALSE))
  data.table(rep=z$rep, method=ss$method, estimate=ss$estimate, se=ss$se, returned=TRUE)
}, mc.cores=40)
A <- rbindlist(res[sapply(res, is.data.table)])
fwrite(A, "rev2/T9_authors_path.csv")
A[, sep := (is.finite(estimate) & abs(estimate)>10) | (is.finite(se) & se>10)]
A[, usable := is.finite(estimate) & is.finite(se) & se>0 & se<10 & !sep]
cat("=== via degli autori (analyze_gee_smallsample, base scale.fix=TRUE), N=10, modello pieno ===\n")
print(A[, .(righe=.N, restituiti=sum(returned), no_se=sum(returned & !is.finite(se)),
            usable=sum(usable), pct_usable_500=round(100*sum(usable)/500,2),
            pct_no_se_of_ret=round(100*sum(returned & !is.finite(se))/sum(returned),2)), by=method])
cat("\n=== confronto con la mia riesecuzione (T3_rerun.R, geesmv chiamata diretta) ===\n")
B <- fread("rev2/T3_subsample_rerun.csv")[spec=="pieno" & N==10 & method %in% c("GEE-MD","GEE-KC","GEE-FG")]
B[, sep := (is.finite(estimate) & abs(estimate)>10) | (is.finite(se) & se>10)]
B[, usable := is.finite(estimate) & is.finite(se) & se>0 & se<10 & !sep]
print(B[, .(usable=sum(usable), pct_usable_500=round(100*sum(usable)/500,2)), by=method])
cat("\npubblicato in §3.16: FG 32.8, KC 32.6, MD 14.6 (=164, 163, 73 su 500)\n")
cat("\n=== differenza appaiata, replica per replica (FG) ===\n")
mm <- merge(A[method=="gee_fg", .(rep, u_aut=usable, se_aut=se)],
            B[method=="GEE-FG", .(rep, u_mio=usable, se_mio=se)], by="rep")
print(table(autori=mm$u_aut, mio=mm$u_mio))
cat("SE identiche dove entrambe presenti:",
    sum(abs(mm$se_aut-mm$se_mio) < 1e-8, na.rm=TRUE), "su",
    sum(is.finite(mm$se_aut) & is.finite(mm$se_mio)), "\n")
cat("\nDONE_T9\n")
