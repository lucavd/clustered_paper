
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(geepack); library(geesmv); library(dplyr); library(MASS)
                  library(fastGHQuad); library(parallel)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_ar1_data.R"), envir=globalenv())
sys.source(file.path(RD,"ar1","generation_models_ar1.R"), envir=globalenv())
sys.source(file.path(RD,"manual","generate_data.R"), envir=globalenv())
CORS <- c("independence","exchangeable","AR-M")
guard <- function(expr, secs=8) {
  p <- parallel::mcparallel(expr)
  z <- parallel::mccollect(p, wait=FALSE, timeout=secs)
  if (is.null(z)) { tools::pskill(p$pid, tools::SIGKILL); parallel::mccollect(p); return(NULL) }
  z[[1]]
}
fit_all <- function(d, cs) {
  d <- as.data.frame(d); d <- d[order(d$id,d$time),]; d$idf <- factor(d$id)
  gcs <- if (cs=="AR-M") "ar1" else cs
  m <- tryCatch(geeglm(binY ~ time*trt, id=idf, data=d, family=binomial, corstr=gcs, scale.fix=TRUE),
                error=function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)$coefficients
  if (!"time:trt" %in% rownames(s)) return(NULL)
  idx <- which(names(coef(m))=="time:trt")
  est <- s["time:trt",1]
  out <- data.frame(corstr=cs, method="gee", estimate=est, se=s["time:trt",2])
  for (nm in c("md","kc","fg")) {
    fn <- get(paste0("GEE.var.",nm), asNamespace("geesmv"))
    v <- guard(tryCatch(fn(binY ~ time*trt, id="idf", family=binomial, data=d, corstr=cs),
                        error=function(e) NULL))
    if (!is.null(v) && is.list(v) && !is.null(v$cov.beta) && length(v$cov.beta)>=idx) {
      sev <- suppressWarnings(sqrt(v$cov.beta[idx]))
      out <- rbind(out, data.frame(corstr=cs, method=paste0("gee_",nm), estimate=est, se=sev))
    }
  }
  out
}
CELLS <- expand.grid(N=c(10,20,45,105), ntime=c(4,6,10), stringsAsFactors=FALSE)
NREP <- 250
one <- function(k) {
  cc <- CELLS[k,]; acc <- list()
  for (i in seq_len(NREP)) {
    set.seed(88000L + 733L*k + i)
    d0 <- tryCatch(generate_data_model1(cc$N, cc$ntime, i, "manual"), error=function(e) NULL)
    d1 <- tryCatch(generate_data_model1(cc$N, cc$ntime, i, "mmlb"),   error=function(e) NULL)
    for (cs in CORS) {
      if (!is.null(d0)) { a <- fit_all(d0,cs); if(!is.null(a)){a$arm<-"H0"; acc[[length(acc)+1]]<-a} }
      if (!is.null(d1)) { b <- fit_all(d1,cs); if(!is.null(b)){b$arm<-"H1"; acc[[length(acc)+1]]<-b} }
    }
  }
  if(!length(acc)) return(NULL)
  r <- bind_rows(acc); r$N <- cc$N; r$ntime <- cc$ntime; r
}
res <- mclapply(seq_len(nrow(CELLS)), function(k) try(one(k), silent=TRUE), mc.cores=12)
R <- bind_rows(Filter(is.data.frame,res))
saveRDS(R, path.expand("~/MMM/rev2/PoW_B2_raw.rds"))
U <- R %>% filter(is.finite(estimate), is.finite(se), abs(estimate)<10, se>0, se<10) %>%
  mutate(p=2*pnorm(-abs(estimate/se)))
cat("### rows per method x corstr\n"); print(as.data.frame(U %>% count(method,corstr)), row.names=FALSE)
cat("\n### TYPE I ERROR (nominal 5)\n")
print(as.data.frame(U %>% filter(arm=="H0") %>% group_by(method,corstr) %>%
  summarise(n=n(), rate=round(100*mean(p<0.05),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=c(rate,n))), row.names=FALSE)
cat("\n### WORST CASE across N\n")
print(as.data.frame(U %>% filter(arm=="H0") %>% group_by(method,corstr,N) %>%
  summarise(rate=100*mean(p<0.05), .groups="drop") %>% group_by(method,corstr) %>%
  summarise(wc=round(max(abs(rate-5)),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=wc)), row.names=FALSE)
cat("\n### LEVEL by N\n")
print(as.data.frame(U %>% filter(arm=="H0") %>% group_by(method,corstr,N) %>%
  summarise(rate=round(100*mean(p<0.05),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=rate)), row.names=FALSE)
cat("\n### COVERAGE (H1 target 0.2)\n")
print(as.data.frame(U %>% filter(arm=="H1") %>% group_by(method,corstr) %>%
  summarise(n=n(), cov=round(100*mean(abs(estimate-0.2)<=qnorm(0.975)*se),2), .groups="drop") %>%
  tidyr::pivot_wider(names_from=corstr, values_from=c(cov,n))), row.names=FALSE)
cat("\n### max |level difference| across corstr, per method\n")
mx <- U %>% filter(arm=="H0") %>% group_by(method,corstr) %>% summarise(rate=100*mean(p<0.05), .groups="drop") %>%
  group_by(method) %>% summarise(spread=round(max(rate)-min(rate),3))
print(as.data.frame(mx), row.names=FALSE)
