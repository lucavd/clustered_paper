.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(lme4); library(geepack); library(MMLB)
                  library(data.table); library(parallel)})
setwd(path.expand("~/MMM"))
TERM <- "regionposterior"
NCORE <- 40

## hard timeout: uccide il figlio se supera secs (protegge anche i loop in C)
hto <- function(expr, secs) {
  p <- parallel::mcparallel(expr)
  r <- parallel::mccollect(p, wait = FALSE, timeout = secs)
  if (is.null(r)) {
    tools::pskill(p$pid, tools::SIGKILL); parallel::mccollect(p, wait = FALSE)
    return(structure("timeout", class = "hard_timeout"))
  }
  r[[1]]
}
TO <- c(mmm = 120, gee = 60, glmm = 120, smv = 5)

d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along)
d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
d <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior",
                                   grepl("P$", as.character(section)) ~ "posterior",
                                   TRUE ~ "lateral"),
                         levels = c("anterior","lateral","posterior")),
         pt = factor(paste(id, tp, sep="_")), id = droplevels(factor(id)))

row1 <- function(meth, est, se, status)
  data.table(method=meth, estimate=est, se=se, status=status)

fit_set <- function(dat, f) {
  dat <- dat %>% arrange(id) %>%
    mutate(id = droplevels(factor(id)), across(where(is.factor), droplevels))
  D <- as.data.frame(dat)
  out <- list()

  ## --- MMM ---
  m <- hto(try(mm(f, t.formula = ~1, data = D, id = id, verbose = FALSE), silent=TRUE), TO["mmm"])
  if (inherits(m, "hard_timeout")) out$MMM <- row1("MMM", NA, NA, "timeout")
  else if (inherits(m, "try-error")) out$MMM <- row1("MMM", NA, NA, "error")
  else if (!TERM %in% names(m$beta)) out$MMM <- row1("MMM", NA, NA, "no_term")
  else { V <- if (!is.null(m$rob.cov)) m$rob.cov else m$mod.cov
         se <- try(sqrt(V[TERM,TERM]), silent=TRUE)
         out$MMM <- row1("MMM", unname(m$beta[TERM]),
                         if (inherits(se,"try-error")) NA_real_ else as.numeric(se), "ok") }

  ## --- GEE ---
  g <- hto(try(geeglm(f, family=binomial, id=id, corstr="exchangeable", data=D), silent=TRUE), TO["gee"])
  if (inherits(g, "hard_timeout")) { out$GEE <- row1("GEE", NA, NA, "timeout"); g <- NULL }
  else if (inherits(g, "try-error")) { out$GEE <- row1("GEE", NA, NA, "error"); g <- NULL }
  else { s <- try(summary(g)$coefficients, silent=TRUE)
    if (inherits(s,"try-error") || !TERM %in% rownames(s)) { out$GEE <- row1("GEE",NA,NA,"no_term"); g <- NULL }
    else out$GEE <- row1("GEE", s[TERM,1], s[TERM,2], "ok") }

  ## --- tre correzioni geesmv, sulla stima puntuale del GEE ---
  est_g <- if (!is.null(g)) unname(coef(g)[TERM]) else NA_real_
  for (nm in c("md","kc","fg")) {
    lab <- paste0("GEE-", toupper(nm))
    if (is.null(g)) { out[[lab]] <- row1(lab, NA, NA, "no_base_gee"); next }
    v <- hto(try(do.call(getExportedValue("geesmv", paste0("GEE.var.", nm)),
              list(f, id="id", family=binomial, data=D, corstr="exchangeable")), silent=TRUE), TO["smv"])
    if (inherits(v,"hard_timeout")) { out[[lab]] <- row1(lab, est_g, NA, "timeout"); next }
    if (inherits(v,"try-error")) { out[[lab]] <- row1(lab, est_g, NA, "error"); next }
    cn <- try(names(coef(g)), silent=TRUE); k <- match(TERM, cn)
    vv <- try(as.numeric(v$cov.beta)[k], silent=TRUE)
    if (inherits(vv,"try-error") || is.na(vv) || !is.finite(vv) || vv <= 0)
      out[[lab]] <- row1(lab, est_g, NA, "bad_var")
    else out[[lab]] <- row1(lab, est_g, sqrt(vv), "ok")
  }

  ## --- GLMM ---
  q <- hto(try(suppressWarnings(glmer(update(f, . ~ . + (1|id)), family=binomial, data=dat,
        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))), silent=TRUE), TO["glmm"])
  if (inherits(q,"hard_timeout")) out$GLMM <- row1("GLMM", NA, NA, "timeout")
  else if (inherits(q,"try-error")) out$GLMM <- row1("GLMM", NA, NA, "error")
  else { s <- try(summary(q)$coefficients, silent=TRUE)
    if (inherits(s,"try-error") || !TERM %in% rownames(s)) out$GLMM <- row1("GLMM",NA,NA,"no_term")
    else out$GLMM <- row1("GLMM", s[TERM,1], s[TERM,2], "ok") }

  rbindlist(out)
}

FS <- list(pieno   = lus_severe ~ region + disease + invasive + age_weeks,
           ridotto = lus_severe ~ region,
           ridotto_age = lus_severe ~ region + age_weeks)

## ---- riferimento a coorte piena ----
cat("=== riferimento coorte piena ===\n")
REF <- rbindlist(lapply(names(FS), function(fn) fit_set(d, FS[[fn]])[, spec := fn]))
fwrite(REF, "rev2/T3_ref.csv"); print(REF, digits=6)

## ---- 2500 estrazioni, seme 20260807 ----
draws <- readRDS("rev2/T3_draws.rds")
cat("\nestrazioni:", length(draws), " (seme 20260807)\n")
t0 <- Sys.time()
run_spec <- function(fn) {
  idx <- if (fn == "ridotto_age") which(sapply(draws, function(z) z$N == 10)) else seq_along(draws)
  res <- mclapply(idx, function(i) {
    z <- draws[[i]]
    s <- d[d$id %in% z$sel, ]
    r <- try(fit_set(s, FS[[fn]]), silent=TRUE)
    if (inherits(r,"try-error"))
      r <- rbindlist(lapply(c("MMM","GEE","GEE-MD","GEE-KC","GEE-FG","GLMM"),
                            function(mm) row1(mm, NA, NA, "outer_error")))
    r[, `:=`(spec=fn, N=z$N, rep=z$rep, nrows=nrow(s), nev=sum(s$lus_severe),
             lev_inv=length(unique(as.character(s$invasive))),
             lev_dis=length(unique(as.character(s$disease))))]
    r
  }, mc.cores = NCORE, mc.preschedule = TRUE)
  bad <- sapply(res, function(x) !is.data.table(x))
  cat(sprintf("spec=%s: figli falliti=%d\n", fn, sum(bad)))
  rbindlist(res[!bad])
}
ALL <- rbindlist(lapply(names(FS), run_spec))
cat(sprintf("tempo totale: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
fwrite(ALL, "rev2/T3_subsample_rerun.csv")
cat("righe scritte:", nrow(ALL), "\n")
print(ALL[, .N, by=.(spec,method)][order(spec,method)])
cat("\n=== tasso di timeout per metodo e spec (su 2500) ===\n")
print(dcast(ALL[, .(pct_timeout=round(100*mean(status=="timeout"),3),
                    n=.N), by=.(spec,method)], method ~ spec, value.var="pct_timeout"))
cat("\nDONE_T3\n")
