
# PoW-B: the limitations section calls the exchangeable-only working correlation
# "the principal gap". This measures it: on GLMM-AR1 data (which has a genuine
# AR(1) latent process, rho=0.6) refit the GEE and the three corrections under
# independence, exchangeable and ar1 working correlations, and compare level,
# coverage and the corrections' worst case.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(lme4); library(geepack); library(geesmv); library(dplyr)
                  library(MASS); library(fastGHQuad); library(parallel)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_ar1_data.R"), envir=globalenv())
sys.source(file.path(RD,"ar1","generation_models_ar1.R"), envir=globalenv())
sys.source(file.path(RD,"manual","generate_data.R"), envir=globalenv())

CORS <- c("independence","exchangeable","ar1")
CELLS <- expand.grid(N=c(10,20,45,105), ntime=c(4,6,10), stringsAsFactors=FALSE)
NREP <- 300
fit_all <- function(d, corstr) {
  d <- as.data.frame(d); d <- d[order(d$id, d$time),]; d$idf <- factor(d$id)
  m <- tryCatch(geeglm(binY ~ time*trt, id=idf, data=d, family=binomial,
                       corstr=corstr, scale.fix=TRUE), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)$coefficients
  if (!"time:trt" %in% rownames(s)) return(NULL)
  est <- s["time:trt",1]; se <- s["time:trt",2]
  out <- data.frame(corstr=corstr, method="gee", estimate=est, se=se)
  f <- binY ~ time*trt
  for (nm in c("md","kc","fg")) {
    fn <- get(paste0("GEE.var.", nm), asNamespace("geesmv"))
    v <- tryCatch({
      p <- parallel::mcparallel(fn(f, id="idf", family=binomial, data=d, corstr.cor=corstr))
      z <- parallel::mccollect(p, wait=FALSE, timeout=8)
      if (is.null(z)) { tools::pskill(p$pid, tools::SIGKILL); parallel::mccollect(p); NULL } else z[[1]]
    }, error=function(e) NULL)
    if (!is.null(v) && !inherits(v,"try-error")) {
      idx <- which(names(coef(m))=="time:trt")
      se2 <- sqrt(v$cov.beta[idx])
      out <- rbind(out, data.frame(corstr=corstr, method=paste0("gee_",nm), estimate=est, se=se2))
    }
  }
  out
}
one <- function(k) {
  cc <- CELLS[k,]; acc <- list()
  for (i in seq_len(NREP)) {
    set.seed(88000L + 733L*k + i)
    d0 <- tryCatch(generate_data_model1(cc$N, cc$ntime, i, "manual"), error=function(e) NULL)
    d1 <- tryCatch(generate_data_model1(cc$N, cc$ntime, i, "mmlb"),   error=function(e) NULL)
    for (cs in CORS) {
      if (!is.null(d0)) { a <- fit_all(d0, cs); if (!is.null(a)) { a$arm <- "H0"; acc[[length(acc)+1]] <- a } }
      if (!is.null(d1)) { b <- fit_all(d1, cs); if (!is.null(b)) { b$arm <- "H1"; acc[[length(acc)+1]] <- b } }
    }
  }
  if (!length(acc)) return(NULL)
  r <- bind_rows(acc); r$N <- cc$N; r$ntime <- cc$ntime; r
}
res <- mclapply(seq_len(nrow(CELLS)), function(k) try(one(k), silent=TRUE), mc.cores=12)
R <- bind_rows(Filter(is.data.frame, res))
saveRDS(R, path.expand("~/MMM/rev2/PoW_B_raw.rds"))
U <- R %>% filter(is.finite(estimate), is.finite(se), abs(estimate)<10, se>0, se<10) %>%
  mutate(p = 2*pnorm(-abs(estimate/se)))
cat("=== TYPE I ERROR by working correlation (nominal 5) ===\n")
t1 <- U %>% filter(arm=="H0") %>% group_by(method, corstr) %>%
  summarise(n=n(), rate=round(100*mean(p<0.05),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=c(rate,n))
print(as.data.frame(t1), row.names=FALSE)
cat("\n=== worst case across N (H0) ===\n")
wc <- U %>% filter(arm=="H0") %>% group_by(method, corstr, N) %>%
  summarise(rate=100*mean(p<0.05), .groups="drop") %>%
  group_by(method, corstr) %>% summarise(wc=round(max(abs(rate-5)),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=wc)
print(as.data.frame(wc), row.names=FALSE)
cat("\n=== level by N and corstr ===\n")
byN <- U %>% filter(arm=="H0") %>% group_by(method, corstr, N) %>%
  summarise(rate=round(100*mean(p<0.05),2), n=n(), .groups="drop")
print(as.data.frame(byN %>% tidyr::pivot_wider(names_from=corstr, values_from=c(rate,n))), row.names=FALSE)
cat("\n=== coverage (H1, target 0.2, nominal 95) ===\n")
cv <- U %>% filter(arm=="H1") %>% group_by(method, corstr) %>%
  summarise(n=n(), cov=round(100*mean(abs(estimate-0.2) <= qnorm(0.975)*se),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=c(cov,n))
print(as.data.frame(cv), row.names=FALSE)
write.csv(byN, path.expand("~/MMM/rev2/PoW_B_workingcorr.csv"), row.names=FALSE)
