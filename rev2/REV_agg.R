
# REV_agg.R -- independent re-aggregation for review nr.2. Reads only the raw
# per-replicate chunks in ~/MMM/results_full. Nothing is taken from tables/.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
OUT <- path.expand("~/MMM/rev2"); dir.create(OUT, showWarnings = FALSE)

fs <- list.files(RES, pattern="[.]rds$")
key <- sub("_c[0-9]+[.]rds$", "", fs)
tab <- table(key)
cat("files:", length(fs), " groups:", length(tab), " groups with 10 chunks:", sum(tab==10), "\n")
cat("groups with != 10 chunks:", sum(tab!=10), "\n")
if (any(tab!=10)) print(head(tab[tab!=10], 20))

WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
LRT  <- "glmer_lrt"
ALL8 <- c(WALD, LRT)
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)
NP   <- c("1"=4,"2"=5,"3"=4,"4"=4,"5"=5)
zq <- qnorm(0.975)

grps <- names(tab)

one <- function(g) {
  f <- file.path(RES, fs[key==g]); f <- sort(f)
  L <- vector("list", length(f))
  for (i in seq_along(f)) {
    x <- try(as.data.table(readRDS(f[i])), silent=TRUE)
    if (inherits(x,"try-error")) return(NULL)
    x <- x[term=="time:trt"]
    if (!nrow(x)) next
    x[, chunk := i]
    keep <- intersect(c("method","N","ntime","label","model","converged","sim_id",
                        "estimate","se","p","dgp","chunk"), names(x))
    L[[i]] <- x[, ..keep]
  }
  d <- rbindlist(L, fill=TRUE)
  if (!nrow(d)) return(NULL)
  d[, model := as.integer(model)]
  d[, rid := paste0(chunk,"_",sim_id)]
  d[, converged := as.logical(converged)]
  d[is.na(converged), converged := FALSE]
  # ---- usability, exactly as documented in the methods section
  d[, fin := is.finite(estimate) & is.finite(se) & is.finite(p)]
  d[, usable := ifelse(method==LRT,
                       converged & is.finite(p),
                       converged & fin & abs(estimate) < 10 & se > 0 & se < 10)]
  # diagnostics on WHY a Wald row is unusable
  d[, no_se   := method!=LRT & converged & (!is.finite(se) | is.na(se))]
  d[, big_est := method!=LRT & converged & is.finite(estimate) & abs(estimate) >= 10]
  m <- d$model[1]; nn <- d$N[1]; tt <- d$ntime[1]; lb <- d$label[1]; dg <- d$dgp[1]
  dfree <- nn - NP[[as.character(m)]]
  truth <- MARG[[as.character(m)]]
  # attempted replicates present in the cell
  n_attempt <- uniqueN(d$rid)
  u <- d[usable==TRUE]
  # intersection: replicates usable for all 8 methods
  cnt <- u[, .(nm = uniqueN(method)), by=rid]
  inter <- cnt[nm==length(ALL8), rid]
  mk <- function(sub, setlab) {
    if (!nrow(sub)) return(NULL)
    sub[, .(set=setlab, n=.N,
            k_rej = sum(p < 0.05),
            k_rej_t = if (dfree > 0) sum(method!=LRT & 2*pt(-abs(estimate/se), df=dfree) < 0.05) +
                                     sum(method==LRT & p < 0.05) else NA_integer_,
            mean_est = mean(estimate), sd_est = sd(estimate), mean_se = mean(se),
            k_cov = sum(is.finite(se) & abs(estimate-truth) <= zq*se),
            n_cov = sum(is.finite(se))), by=method]
  }
  out <- rbindlist(list(mk(u,"own"), mk(u[rid %in% inter],"intersection")), fill=TRUE)
  diag <- d[, .(n_present=.N, n_conv=sum(converged), n_usable=sum(usable),
                n_no_se=sum(no_se), n_big=sum(big_est)), by=method]
  out <- merge(out, diag, by="method", all.x=TRUE)
  out[, `:=`(dgp=dg, model=m, N=nn, ntime=tt, label=lb,
             n_attempt=n_attempt, n_inter=length(inter), truth=truth, df=dfree)]
  out[]
}

res <- mclapply(grps, function(g) try(one(g), silent=TRUE), mc.cores = 45)
bad <- which(!sapply(res, is.data.table))
cat("failed groups:", length(bad), "\n")
if (length(bad)) print(grps[bad][1:min(5,length(bad))])
R <- rbindlist(res[sapply(res, is.data.table)], fill=TRUE)
fwrite(R, file.path(OUT, "REV_cells.csv"))
cat("rows:", nrow(R), " design points:", uniqueN(R[,.(dgp,model,N,ntime)]), "\n")
cat("total usable replicate-rows (own,H0):", R[set=="own" & label=="manual", sum(n)], "\n")
