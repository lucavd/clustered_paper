
# PoW-C: the manuscript attributes the LRT's 98.8% usability at two occasions to
# "requiring no covariance matrix for the random effects". But do_analysis.R passes
# the LRT a RANDOM-INTERCEPT-ONLY formula (binY ~ mean_rhs + (1|id)) in every
# scenario, including M3-M5 where the GLMM-Wald arm fits (1 + time | id).
# So the 98.8% may be the specification, not the test. This refits the LRT with
# the SAME random structure the Wald arm uses, and reports usability and level.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(lme4); library(dplyr); library(MASS); library(fastGHQuad); library(parallel)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_glmm_data.R"), envir=globalenv())
for (f in c("generation_model3.R","generation_model5.R","generate_data.R"))
  sys.source(file.path(RD,"manual",f), envir=globalenv())

lrt_fit <- function(d, full_f, red_f) {
  d <- as.data.frame(d); d$id <- factor(d$id)
  ctl <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
  fu <- tryCatch(suppressMessages(glmer(full_f, data=d, family=binomial, control=ctl)), error=function(e) NULL)
  if (is.null(fu)) return(c(NA,NA))
  re <- tryCatch(suppressMessages(glmer(red_f, data=d, family=binomial, control=ctl)), error=function(e) NULL)
  if (is.null(re)) return(c(NA,NA))
  a <- tryCatch(anova(re, fu), error=function(e) NULL)
  if (is.null(a)) return(c(NA,NA))
  c(a$Chisq[2], a$`Pr(>Chisq)`[2])
}
CELLS <- expand.grid(N=c(10,20,45), ntime=c(2,3,4), model=3, stringsAsFactors=FALSE)
NREP <- 250
one <- function(k) {
  cc <- CELLS[k,]; acc <- list()
  for (i in seq_len(NREP)) {
    set.seed(66000L + 811L*k + i)
    d <- tryCatch(generate_data_model3(cc$N, cc$ntime, i, "manual"), error=function(e) NULL)
    if (is.null(d)) next
    # (a) as run: random intercept only in both nested fits
    a1 <- lrt_fit(d, binY ~ time*trt + (1|id), binY ~ time + trt + (1|id))
    # (b) as declared: random intercept + random slope on time in both fits
    a2 <- lrt_fit(d, binY ~ time*trt + (1 + time|id), binY ~ time + trt + (1 + time|id))
    acc[[length(acc)+1]] <- data.frame(N=cc$N, ntime=cc$ntime, rep=i,
      chisq_ri=a1[1], p_ri=a1[2], chisq_rs=a2[1], p_rs=a2[2])
  }
  bind_rows(acc)
}
res <- mclapply(seq_len(nrow(CELLS)), function(k) try(one(k), silent=TRUE), mc.cores=9)
R <- bind_rows(Filter(is.data.frame, res))
saveRDS(R, path.expand("~/MMM/rev2/PoW_C_raw.rds"))
cat("=== LRT usability and level: random intercept (as run) vs random slope (as declared) ===\n")
tab <- R %>% group_by(N, ntime) %>%
  summarise(attempts=n(),
            usable_ri=round(100*mean(is.finite(p_ri)),1),
            usable_rs=round(100*mean(is.finite(p_rs)),1),
            level_ri=round(100*mean(p_ri<0.05, na.rm=TRUE),2),
            level_rs=round(100*mean(p_rs<0.05, na.rm=TRUE),2), .groups="drop")
print(as.data.frame(tab), row.names=FALSE)
cat("\n=== pooled ===\n")
print(as.data.frame(R %>% summarise(n=n(),
  usable_ri=round(100*mean(is.finite(p_ri)),2), usable_rs=round(100*mean(is.finite(p_rs)),2),
  level_ri=round(100*mean(p_ri<0.05,na.rm=TRUE),2), level_rs=round(100*mean(p_rs<0.05,na.rm=TRUE),2))),
  row.names=FALSE)
cat("\n=== at T=2 only ===\n")
print(as.data.frame(R %>% filter(ntime==2) %>% summarise(n=n(),
  usable_ri=round(100*mean(is.finite(p_ri)),2), usable_rs=round(100*mean(is.finite(p_rs)),2),
  level_ri=round(100*mean(p_ri<0.05,na.rm=TRUE),2), level_rs=round(100*mean(p_rs<0.05,na.rm=TRUE),2))),
  row.names=FALSE)
write.csv(tab, path.expand("~/MMM/rev2/PoW_C_lrt_spec.csv"), row.names=FALSE)
