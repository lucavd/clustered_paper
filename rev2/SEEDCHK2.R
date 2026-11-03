.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(MMLB); library(dplyr); library(MASS); library(fastGHQuad)})
RD <- path.expand("~/MMM/.rundir")
sp <- function(f) sys.source(file.path(RD, f), envir=globalenv())
for (f in c("generation_model1.R","generation_model3.R","generate_data.R","make_all_data.R")) sp(f)

ch <- as.data.table(readRDS(path.expand("~/MMM/results_full/GenBinaryY_m1_N10_T4_manual_c01.rds")))
cm <- ch[term=="time:trt" & method=="mmm"][order(sim_id)]
target <- cm[sim_id %in% 1:4]$estimate
cat("chunk (sim_id 1..4):", paste(sprintf("%.6f", target), collapse="  "), "\n\n")

fitest <- function(d) {
  mo <- tryCatch(suppressWarnings(mm(binY ~ time*trt, lv.formula=~1, t.formula=~1,
                                     data=as.data.frame(d), id=id, q=20)), error=function(e) NULL)
  if (is.null(mo)) return(NA_real_)
  summary(mo)$mean.table["time:trt",1]
}
try_kind <- function(kind, seedfun, tag) {
  out <- numeric(0); nobs <- integer(0)
  for (i in 1:4) {
    suppressWarnings(RNGkind(kind))
    set.seed(seedfun(i))
    d <- tryCatch(make_all_data(10, 4, i, 1L, "manual"), error=function(e) NULL)
    if (is.null(d)) { out <- c(out, NA_real_); nobs <- c(nobs, NA_integer_); next }
    out <- c(out, fitest(d)); nobs <- c(nobs, nrow(d))
  }
  cat(sprintf("%-42s nobs=%s  est=%s\n", tag, paste(nobs, collapse=","),
              paste(sprintf("%.6f", out), collapse="  ")))
  invisible(out)
}
S <- function(i) 20260805L + 1000L*1L + i
try_kind("Mersenne-Twister", S, "Mersenne-Twister (default --vanilla)")
try_kind("L'Ecuyer-CMRG",    S, "L'Ecuyer-CMRG (crew/mirai)")
try_kind("Wichmann-Hill",    S, "Wichmann-Hill")
try_kind("Marsaglia-Multicarry", S, "Marsaglia-Multicarry")
try_kind("Super-Duper",      S, "Super-Duper")
try_kind("Knuth-TAOCP-2002", S, "Knuth-TAOCP-2002")

cat("\n-- L'Ecuyer con semi alternativi --\n")
for (nm in names(list(a=1))) NULL
alt <- list("20260805+i"=function(i) 20260805L+i, "1000+i"=function(i) 1000L+i,
            "i"=function(i) i, "20260805+2000+i"=function(i) 20260805L+2000L+i)
for (nm in names(alt)) try_kind("L'Ecuyer-CMRG", alt[[nm]], paste("L'Ecuyer,", nm))
cat("\nRNGkind corrente:", paste(RNGkind(), collapse=" / "), "\n")
cat("DONE_SEEDCHK2\n")
