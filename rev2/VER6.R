
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(MMLB); library(lme4); library(geepack); library(geesmv); library(parallel)})
d <- readRDS(path.expand("~/MMM/ARTICOLO 2026/LURID_bis/lurid_long.rds")) |>
  filter(!is.na(lus_severe)) |>
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior",
                                   grepl("P$", as.character(section)) ~ "posterior",
                                   TRUE ~ "lateral"), levels=c("anterior","lateral","posterior")),
         id = droplevels(factor(id)))
TERM <- "regionposterior"
tmo <- function(expr, s=20) { p <- mcparallel(expr); r <- mccollect(p, wait=FALSE, timeout=s)
  if (is.null(r)) { kill(p$pid, 9); mccollect(p); return(NULL) }; r[[1]] }
fitset <- function(dat, rhs) {
  dat <- dat |> arrange(id) |> mutate(id=droplevels(factor(id))) |> as.data.frame()
  f <- reformulate(rhs, response="lus_severe"); out <- list()
  m <- tmo(try(mm(f, t.formula=~1, data=dat, id=id, verbose=FALSE), silent=TRUE))
  out$MMM <- data.frame(method="MMM", se=if(!is.null(m) && !inherits(m,"try-error") && TERM %in% names(m$beta))
    sqrt(m$mod.cov[TERM,TERM]) else NA_real_)
  g <- tmo(try(geeglm(f, family=binomial, id=id, corstr="exchangeable", data=dat), silent=TRUE))
  ok <- !is.null(g) && !inherits(g,"try-error")
  out$GEE <- data.frame(method="GEE", se=if(ok) summary(g)$coefficients[TERM,2] else NA_real_)
  for (nm in c("md","kc","fg")) {
    fn <- switch(nm, md=GEE.var.md, kc=GEE.var.kc, fg=GEE.var.fg)
    v <- tmo(try(fn(f, id="id", family=binomial, data=dat, corstr="exchangeable"), silent=TRUE), 20)
    se <- NA_real_
    if (!is.null(v) && !inherits(v,"try-error")) {
      nmz <- names(v$cov.beta); if (is.null(nmz)) nmz <- colnames(model.matrix(f, dat))
      j <- match(TERM, nmz); if (!is.na(j)) se <- sqrt(v$cov.beta[j])
    }
    out[[toupper(nm)]] <- data.frame(method=paste0("GEE-",toupper(nm)), se=se)
  }
  q <- tmo(try(suppressWarnings(glmer(update(f, . ~ . + (1|id)), family=binomial, data=dat,
        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))), silent=TRUE), 30)
  out$GLMM <- data.frame(method="GLMM", se=if(!is.null(q) && !inherits(q,"try-error") &&
    TERM %in% rownames(summary(q)$coefficients)) summary(q)$coefficients[TERM,2] else NA_real_)
  bind_rows(out)
}
set.seed(20260807); ids <- levels(d$id); REPS <- 150; N <- 10
run <- function(rhs, tag) {
  set.seed(20260807)
  L <- lapply(seq_len(REPS), function(r) {
    sel <- sample(ids, N)
    fitset(d |> filter(id %in% sel), rhs) |> mutate(rep=r)
  })
  S <- bind_rows(L)
  cat(sprintf("\n=== %s (N=%d, %d sottocampioni, seme 20260807) ===\n", tag, N, REPS))
  print(as.data.frame(S |> group_by(method) |>
    summarise(pct_SE_presente = round(100*mean(is.finite(se)),1),
              pct_SE_assente  = round(100*mean(!is.finite(se)),1), .groups="drop")), row.names=FALSE)
}
run(c("region","disease","invasive","age_weeks"), "modello medio COMPLETO")
run(c("region"), "modello medio RIDOTTO (solo region) -- l'esito (a) di §6.2")
run(c("region","age_weeks"), "modello medio RIDOTTO (region + age_weeks)")
cat("\nDONE6\n")
