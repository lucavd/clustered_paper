
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(lme4); library(MMLB); library(geepack)})
d <- readRDS(path.expand("~/MMM/ARTICOLO 2026/LURID_bis/lurid_long.rds")) |>
  filter(!is.na(lus_severe)) |>
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior",
                                   grepl("P$", as.character(section)) ~ "posterior",
                                   TRUE ~ "lateral"),
                         levels=c("anterior","lateral","posterior")),
         id = droplevels(factor(id)))
cat("righe:", nrow(d), " pazienti:", nlevels(d$id), "\n")
cat("zone per paziente: "); print(summary(as.numeric(table(d$id))))
cat("invasive: "); print(table(d$invasive)); cat("pazienti invasive: ")
print(table(unique(d[,c("id","invasive")])$invasive))
TERM <- "regionposterior"
fitset <- function(dat, rhs, tag) {
  dat <- dat |> arrange(id) |> as.data.frame()
  f <- reformulate(rhs, response="lus_severe")
  out <- list()
  m <- try(mm(f, t.formula=~1, data=dat, id=id, verbose=FALSE), silent=TRUE)
  if (!inherits(m,"try-error") && TERM %in% names(m$beta)) {
    Vm <- m$mod.cov; Vr <- if (!is.null(m$rob.cov)) m$rob.cov else NULL
    out[[length(out)+1]] <- data.frame(spec=tag, method="MMM", est=unname(m$beta[TERM]),
      se_model=sqrt(Vm[TERM,TERM]), se_rob=if(is.null(Vr)) NA else sqrt(Vr[TERM,TERM]))
  }
  g <- try(geeglm(f, family=binomial, id=id, corstr="exchangeable", data=dat), silent=TRUE)
  if (!inherits(g,"try-error")) { s <- summary(g)$coefficients
    if (TERM %in% rownames(s)) out[[length(out)+1]] <-
      data.frame(spec=tag, method="GEE", est=s[TERM,1], se_model=s[TERM,2], se_rob=NA) }
  gl <- try(glmer(update(f, . ~ . + (1|id)), family=binomial, data=dat,
                  control=glmerControl(optimizer="bobyqa")), silent=TRUE)
  if (!inherits(gl,"try-error")) { s <- summary(gl)$coefficients
    if (TERM %in% rownames(s)) out[[length(out)+1]] <-
      data.frame(spec=tag, method="GLMM", est=s[TERM,1], se_model=s[TERM,2], se_rob=NA) }
  bind_rows(out)
}
full <- c("region","disease","invasive","age_weeks")
res <- bind_rows(
  fitset(d, full, "completo: region+disease+invasive+age_weeks"),
  fitset(d, c("region"), "ridotto: region"),
  fitset(d, c("region","age_weeks"), "ridotto: region+age_weeks"))
cat("\n===== §6.2: contrasto posteriore-anteriore, coorte piena, tre modelli medi =====\n")
print(as.data.frame(res |> mutate(across(where(is.numeric), ~round(.x,4)))), row.names=FALSE)
cat("\ndifferenze rispetto al modello completo, per metodo:\n")
w <- res |> select(spec,method,est,se_model) |> tidyr::pivot_wider(names_from=spec, values_from=c(est,se_model))
print(as.data.frame(w |> mutate(across(where(is.numeric), ~round(.x,4)))), row.names=FALSE)
cat("\n===== R3-S08: tetto ipergeometrico =====\n")
npat <- nlevels(d$id); ninv <- sum(unique(d[,c("id","invasive")])$invasive %in% c(1,TRUE,"yes","Invasive","invasive"))
cat("pazienti:", npat, "\n")
tab <- table(unique(as.data.frame(d)[,c("id","invasive")])$invasive); print(tab)
for (N in c(10,15,20,25,50)) {
  k <- as.numeric(tab[length(tab)])
  cat(sprintf("N=%d  P(nessun invasivo)=%.4f  tetto=%.1f%%\n", N,
    dhyper(0, k, npat-k, N), 100*(1-dhyper(0, k, npat-k, N))))
}
cat("\nDONE2\n")
