# REV2M_b.R -- B: LRT re-run with the DECLARED random-effects structure.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(lme4); library(parallel); library(MASS); library(fastGHQuad)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_glmm_data.R"), envir=globalenv())
for (f in c("manual/generation_model3.R","manual/generate_data.R","manual/make_all_data.R"))
  sys.source(file.path(RD,f), envir=globalenv())

CTRL <- lme4::glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))

lrt_one <- function(d, re) {
  ff <- as.formula(paste("binY ~ time*trt +", re))
  fr <- as.formula(paste("binY ~ time + trt +", re))
  full <- tryCatch(suppressMessages(glmer(ff, data=d, family=binomial, control=CTRL)), error=function(e) NULL)
  if (is.null(full)) return(c(ok=0, p=NA, chisq=NA, est=NA, se=NA))
  red <- tryCatch(suppressMessages(glmer(fr, data=d, family=binomial, control=CTRL)), error=function(e) NULL)
  if (is.null(red)) return(c(ok=0, p=NA, chisq=NA, est=NA, se=NA))
  a <- tryCatch(anova(red, full), error=function(e) NULL)
  if (is.null(a)) return(c(ok=0, p=NA, chisq=NA, est=NA, se=NA))
  p <- a$`Pr(>Chisq)`[2]
  if (!is.finite(p)) return(c(ok=0, p=NA, chisq=NA, est=NA, se=NA))
  s <- summary(full)$coefficients
  c(ok=1, p=p, chisq=a$Chisq[2], est=s["time:trt",1], se=s["time:trt",2])
}

NL <- c(10,20,45); TL <- c(2,3,4); NREP <- 200; MODEL <- 3
grid <- expand.grid(N=NL, ntime=TL, label=c("manual","mmlb"), stringsAsFactors=FALSE)

run_cell <- function(k) {
  g <- grid[k,]
  rows <- lapply(1:NREP, function(i) {
    set.seed(20260805L + 1000L*MODEL + i)
    d <- tryCatch(make_all_data(g$N, g$ntime, i, as.character(MODEL), g$label), error=function(e) NULL)
    if (is.null(d)) return(NULL)
    d <- as.data.frame(d); d$id <- factor(d$id)
    a <- lrt_one(d, "(1|id)")            # as run
    b <- lrt_one(d, "(1 + time|id)")     # as declared in Methods for M3
    data.frame(rep=i, ok_ri=a["ok"], p_ri=a["p"], ok_rs=b["ok"], p_rs=b["p"],
               est_rs=b["est"], se_rs=b["se"])
  })
  r <- bind_rows(rows)
  data.frame(N=g$N, ntime=g$ntime, label=g$label, n_attempt=nrow(r),
             usable_ri = 100*mean(r$ok_ri==1),
             usable_rs = 100*mean(r$ok_rs==1),
             rej_ri = 100*mean(r$p_ri[r$ok_ri==1] < 0.05),
             rej_rs = 100*mean(r$p_rs[r$ok_rs==1] < 0.05),
             n_ri = sum(r$ok_ri==1), n_rs = sum(r$ok_rs==1),
             mean_est_rs = mean(r$est_rs, na.rm=TRUE))
}
res <- mclapply(seq_len(nrow(grid)), function(k) try(run_cell(k), silent=TRUE), mc.cores=18)
ok <- sapply(res, is.data.frame)
cat("cells failed:", sum(!ok), "\n")
R <- bind_rows(res[ok])
cat("\n===== LRT: (1|id) as run  vs  (1+time|id) as declared -- GLMM-native, M3, 200 reps/cell =====\n")
print(as.data.frame(R %>% mutate(across(where(is.numeric), ~round(.x,2)))), row.names=FALSE)
cat("\n===== pooled by T, H0 only (level) =====\n")
print(as.data.frame(R %>% filter(label=="manual") %>% group_by(ntime) %>%
  summarise(usable_asrun=round(mean(usable_ri),1), usable_declared=round(mean(usable_rs),1),
            level_asrun=round(100*sum(rej_ri*n_ri/100)/sum(n_ri),2),
            level_declared=round(100*sum(rej_rs*n_rs/100)/sum(n_rs),2),
            n_asrun=sum(n_ri), n_declared=sum(n_rs), .groups="drop")), row.names=FALSE)
cat("\n===== overall pooled, H0 =====\n")
h0 <- R %>% filter(label=="manual")
cat(sprintf("usable as-run (1|id): %.1f%%   as-declared (1+time|id): %.1f%%\n",
            mean(h0$usable_ri), mean(h0$usable_rs)))
cat(sprintf("level  as-run: %.2f%% (n=%d)   as-declared: %.2f%% (n=%d)\n",
            100*sum(h0$rej_ri*h0$n_ri/100)/sum(h0$n_ri), sum(h0$n_ri),
            100*sum(h0$rej_rs*h0$n_rs/100)/sum(h0$n_rs), sum(h0$n_rs)))
write.csv(R, path.expand("~/MMM/rev2/REV2M_lrt_spec.csv"), row.names=FALSE)
cat("\nDONE_B\n")
