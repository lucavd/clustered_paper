.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(lme4)})
setwd(path.expand("~/MMM"))
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN=seq_along)
d0$tp <- (d0$row_in_id-1L) %/% 10L + 1L
d <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region=factor(case_when(grepl("A$",as.character(section))~"anterior",
                                 grepl("P$",as.character(section))~"posterior",
                                 TRUE~"lateral"), levels=c("anterior","lateral","posterior")),
         pt=factor(paste(id,tp,sep="_")), id=droplevels(factor(id)))
nt <- d %>% group_by(id) %>% summarise(k=n_distinct(pt), .groups="drop")
multi <- as.character(nt$id[nt$k>1]); single <- as.character(nt$id[nt$k==1])
G <- function(f, dat) glmer(f, family=binomial, data=droplevels(dat),
      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
ch <- function(dat) { m1 <- G(lus_severe~region+(1|id), dat); m2 <- G(lus_severe~region+(1|id)+(1|pt), dat)
  as.numeric(2*(logLik(m2)-logLik(m1))) }
tot <- ch(d); c22 <- ch(d[d$id %in% multi,]); c185 <- ch(d[d$id %in% single,])
cat(sprintf("chi2 totale=%.6f  22 multi=%.6f (%.4f%%)  185 single=%.8f  somma parti=%.6f  non-additivita'=%.6f\n",
            tot, c22, 100*c22/tot, c185, c22+c185, tot-(c22+c185)))
## i coefficienti region e la varianza id sono stimati su tutti i 207: fissiamoli
mfull <- G(lus_severe~region+(1|id), d)
b <- fixef(mfull)
X <- model.matrix(~region, data=d[d$id %in% multi,])
off <- as.vector(X %*% b)
dm <- d[d$id %in% multi,]; dm$off <- off
m1o <- glmer(lus_severe ~ 0 + offset(off) + (1|id), family=binomial, data=droplevels(dm),
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m2o <- glmer(lus_severe ~ 0 + offset(off) + (1|id) + (1|pt), family=binomial, data=droplevels(dm),
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
cat(sprintf("chi2 sui 22 con coefficienti region fissati alla coorte piena = %.6f (%.2f%% del totale)\n",
            as.numeric(2*(logLik(m2o)-logLik(m1o))), 100*as.numeric(2*(logLik(m2o)-logLik(m1o)))/tot))
## quanta varianza pt resta se i 22 sono rimossi
m2s <- G(lus_severe~region+(1|id)+(1|pt), d[d$id %in% single,])
cat("var pt nei soli 185 single:", format(unname(unlist(VarCorr(m2s))), digits=6), "\n")
cat("prevalenza per timepoint nei 22 multi:\n")
print(d %>% filter(id %in% multi) %>% group_by(id) %>%
        summarise(k=n_distinct(pt), range_prev=diff(range(tapply(lus_severe, tp, mean))), .groups="drop") %>%
        summarise(n=n(), mean_range=mean(range_prev), max_range=max(range_prev)))
cat("\nDONE_T6\n")
