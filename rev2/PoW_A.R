
# PoW-A: on GenBinaryY, scenarios M3-M5 give the H0 arm a DIFFERENT dependence
# structure from the H1 arm (sigma and gamma differ between the mmlb and manual
# branches). The manuscript states the arms differ in the tested parameter only.
# This measures what that does: refit the null arm with the H1 dependence and
# compare the type I error.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(MMLB); library(lme4); library(geepack); library(geesmv)
                  library(dplyr); library(tibble); library(MASS); library(fastGHQuad)
                  library(parallel)})
SRC <- path.expand("~/MMM/R")
for (f in c("analysis_model1.R","analysis_model2.R","analysis_model3.R","analysis_model4.R",
            "analysis_model5.R","analysis_smallsample.R","do_analysis.R")) sys.source(file.path(SRC,f), envir=globalenv())

# --- generator with the dependence held FIXED across the two arms -------------
gen_m3_matched <- function(N, ntime, sim_id, arm) {
  beta_vec <- if (arm=="H1") c(-2,0.1,0.2,0.2) else c(-2,0.1,0.2,0)
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime-1), N, replace=TRUE))
  id <- rep(seq(N), nclust); trt <- rep(rbinom(N,1,0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x)-1))
  age <- scale(rep(rnorm(N,50,10), nclust))
  data <- data.frame(id,time,trt,age) |> arrange(id,time)
  GenBinaryY(mean.formula=~time*trt, lv.formula=~time, t.formula=~1,
             beta=beta_vec, sigma=c(1,0.01), gamma=1,   # <-- H1 dependence in BOTH arms
             id=id, data=data, q=20, Yname="binY")
}
# --- generator as actually run: H0 has sigma=c(1,0) ---------------------------
gen_m3_asrun <- function(N, ntime, sim_id, arm) {
  beta_vec <- if (arm=="H1") c(-2,0.1,0.2,0.2) else c(-2,0.1,0.2,0)
  sig <- if (arm=="H1") c(1,0.01) else c(1,0)
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime-1), N, replace=TRUE))
  id <- rep(seq(N), nclust); trt <- rep(rbinom(N,1,0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x)-1))
  age <- scale(rep(rnorm(N,50,10), nclust))
  data <- data.frame(id,time,trt,age) |> arrange(id,time)
  GenBinaryY(mean.formula=~time*trt, lv.formula=~time, t.formula=~1,
             beta=beta_vec, sigma=sig, gamma=1, id=id, data=data, q=20, Yname="binY")
}

CELLS <- expand.grid(N=c(10,20,45,105), ntime=c(2,4,10), spec=c("asrun","matched"),
                     stringsAsFactors=FALSE)
NREP <- 400
one <- function(k) {
  cc <- CELLS[k,]; gen <- if (cc$spec=="asrun") gen_m3_asrun else gen_m3_matched
  out <- vector("list", NREP)
  for (i in seq_len(NREP)) {
    set.seed(77000L + 991L*k + i)
    d <- tryCatch(gen(cc$N, cc$ntime, i, "H0"), error=function(e) NULL)
    if (is.null(d)) next
    rows <- lapply(c("mmm","gee","glmer","gee_ss","glmer_lrt"), function(tc)
      tryCatch(do_analysis(d, cc$N, cc$ntime, "manual", 3, tc, i), error=function(e) NULL))
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) out[[i]] <- bind_rows(rows)
  }
  r <- bind_rows(Filter(Negate(is.null), out))
  if (!nrow(r)) return(NULL)
  r$N <- cc$N; r$ntime <- cc$ntime; r$spec <- cc$spec
  r %>% filter(term=="time:trt")
}
res <- mclapply(seq_len(nrow(CELLS)), function(k) try(one(k), silent=TRUE), mc.cores=12)
R <- bind_rows(Filter(function(x) is.data.frame(x), res))
saveRDS(R, path.expand("~/MMM/rev2/PoW_A_raw.rds"))
usable <- function(d) bind_rows(
  d %>% filter(method!="glmer_lrt", converged, is.finite(estimate), is.finite(se), is.finite(p),
               abs(estimate)<10, se>0, se<10),
  d %>% filter(method=="glmer_lrt", converged, is.finite(p)))
U <- usable(R)
tab <- U %>% group_by(spec, method, N, ntime) %>%
  summarise(n=n(), rate=100*mean(p<0.05), .groups="drop")
w <- tab %>% tidyr::pivot_wider(names_from=spec, values_from=c(rate,n))
w$delta <- w$rate_matched - w$rate_asrun
print(as.data.frame(w %>% arrange(method,ntime,N)), row.names=FALSE)
cat("\n=== pooled by method ===\n")
p <- U %>% group_by(spec, method) %>% summarise(n=n(), rate=100*mean(p<0.05), .groups="drop") %>%
  tidyr::pivot_wider(names_from=spec, values_from=c(rate,n)) %>% mutate(delta=rate_matched-rate_asrun)
print(as.data.frame(p), row.names=FALSE)
write.csv(w, path.expand("~/MMM/rev2/PoW_A_h0h1_dependence.csv"), row.names=FALSE)
