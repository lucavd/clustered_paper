# R2_misc.R -- GIRO 2: controlli 6 (semi), 7 (varianza MMLB), 8 (generalizzazione stretta),
#              9 (conteggi per cella pubblicati), 11 (p e gdl del riferimento t).
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(data.table); library(dplyr); library(MMLB); library(MASS)
                  library(parallel); library(fastGHQuad)})

cat("############ 7. quale varianza usa mm() per la mean.table?\n")
src <- readLines(file.path(find.package("MMLB"), "R", "MMLB.rdb"), warn=FALSE, n=1)
cat("codice R compilato in .rdb: ispeziono per deparse.\n")
f <- MMLB::mm
b <- paste(deparse(f), collapse="\n")
cat("argomenti di mm():\n"); print(names(formals(f)))
for (pat in c("sandwich","robust","vcov","hessian","Hessian","solve","crossprod","score",
              "information","observed","expected","mod.cov","varMat","cov")) {
  n <- length(gregexpr(pat, b, fixed=TRUE)[[1]][gregexpr(pat, b, fixed=TRUE)[[1]]>0])
  if (n) cat(sprintf("  '%s' compare %d volte in mm()\n", pat, n))
}
g <- paste(deparse(MMLB::MMLongit), collapse="\n")
cat("\nin MMLongit():\n")
for (pat in c("sandwich","robust","solve(","crossprod","hess","Hess","numDeriv","optim",
              "mod.cov","rob.cov","information","score")) {
  n <- length(gregexpr(pat, g, fixed=TRUE)[[1]][gregexpr(pat, g, fixed=TRUE)[[1]]>0])
  if (n) cat(sprintf("  '%s' compare %d volte\n", pat, n))
}
cat("\nrighe di MMLongit() che costruiscono la covarianza:\n")
gl <- deparse(MMLB::MMLongit)
print(grep("cov|vcov|solve|hess|Hess|sandwich|robust", gl, value=TRUE)[1:25])
cat("\nverifica empirica: la SE di mm() e' model-based o sandwich? Confronto con la SD empirica\n")
cat("(se model-based-only, la SE sottostima quando il modello e' misspecificato).\n")
cat("Vedi la riga mmm in R2_cov: SE/SD empirica per cella.\n")

cat("\n############ 8. GLMM-AR1 con rho=0, sigma2_ar=0 riproduce gen_glmm_data 'exactly'?\n")
RD <- path.expand("~/MMM/.rundir")
e1 <- new.env(); e2 <- new.env()
sys.source(file.path(RD,"gen_glmm_data.R"), envir=e1)
sys.source(file.path(RD,"gen_ar1_data.R"), envir=e2)
gg <- get("gen_glmm_data", envir=e1); ga <- get("gen_ar1_data", envir=e2)
N <- 45; T <- 6
set.seed(999); nc <- ifelse(runif(N)<0.7, T, sample(seq_len(T-1), N, replace=TRUE))
id <- rep(seq(N), nc); tt <- unlist(lapply(nc, function(x) seq(x)-1))
trt <- rep(rbinom(N,1,0.5), nc)
dat <- data.frame(id=id, time=tt, trt=trt)
bc <- c(-2.32, 0.116, 0.232, 0.232)
set.seed(4242); A <- gg(~time*trt, ~1, beta_c=bc, D=1.0, id=id, data=dat, Yname="binY")
set.seed(4242); B <- ga(~time*trt, ~1, beta_c=bc, D=1.0, rho=0, sigma2_ar=0,
                        time="time", id=id, data=dat, Yname="binY")
cat(sprintf("stesso seme, stesso beta_c, rho=0, sigma2_ar=0:\n  outcome identico? %s (concordanza %.4f su %d osservazioni)\n",
  identical(A$binY, B$binY), mean(A$binY==B$binY), nrow(A)))
cat(sprintf("  prevalenza A %.4f vs B %.4f\n", mean(A$binY), mean(B$binY)))
cat(sprintf("  attributi dgp: A '%s' vs B '%s'\n", attr(A,"dgp"), attr(B,"dgp")))
cat("  perche': gen_ar1_data estrae b_subj poi (se sigma2_ar>0) il processo AR1; con sigma2_ar=0\n")
cat("  il ciclo e' saltato, quindi la sequenza di numeri casuali dovrebbe coincidere.\n")
# ripeto su 20 semi per escludere un caso fortunato
ident <- sapply(1:20, function(s) {
  set.seed(s); a <- gg(~time*trt, ~1, beta_c=bc, D=1.0, id=id, data=dat, Yname="binY")
  set.seed(s); b <- ga(~time*trt, ~1, beta_c=bc, D=1.0, rho=0, sigma2_ar=0, time="time",
                       id=id, data=dat, Yname="binY")
  identical(a$binY, b$binY) })
cat(sprintf("  identico su %d semi su 20\n", sum(ident)))
cat("\nma i BRACCI dello studio usano beta_c diversi: native fisso 0.232, AR1 calibrato.\n")
sys.source(file.path(RD,"ar1/generation_models_ar1.R"), envir=globalenv())
sys.source(file.path(RD,"gen_ar1_data.R"), envir=globalenv())
set.seed(20260805L+1000L+1); d1 <- .ar1_design(75, 6, "none")
bcal0 <- as.numeric(calibrate_beta_c_ar1(~time*trt, ~1, .ar1_target(1,"mmlb"), 1.0, 0.0, d1))
bcal5 <- as.numeric(calibrate_beta_c_ar1(~time*trt, ~1, .ar1_target(1,"mmlb"), 1.0, AR1_SIGMA2, d1))
cat(sprintf("  beta_c calibrato con sigma2_ar=0:   %s\n", paste(round(bcal0,5), collapse=", ")))
cat(sprintf("  beta_c calibrato con sigma2_ar=0.5: %s\n", paste(round(bcal5,5), collapse=", ")))
cat(sprintf("  beta_c fisso di GLMM-native:        %s\n", paste(c(-2.32,0.116,0.232,0.232), collapse=", ")))
cat(sprintf("  differenza sull'interazione: calibrato(0) %.5f vs fisso 0.232 -> %.5f\n",
  bcal0[4], bcal0[4]-0.232))

cat("\n############ 11. p per i gdl del riferimento t\n")
for (m in 1:5) {
  mf <- switch(as.character(m), "1"=~time*trt, "2"=~time*trt+age, "3"=~time*trt,
               "4"=~time*trt, "5"=~time*trt+dose)
  dd <- data.frame(time=0:3, trt=c(0,1,0,1), age=rnorm(4), dose=rnorm(4))
  nm <- colnames(model.matrix(mf, dd))
  cat(sprintf("  M%d: %d parametri (%s) ; NP dichiarato in REV_agg.R = %d\n", m, length(nm),
    paste(nm, collapse=", "), c(4,5,4,4,5)[m]))
}
cat("\ngdl usati: N - p, con N = numero di CLUSTER (non di osservazioni). Verifica su REV_cells:\n")
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv"))
print(as.data.frame(unique(R[, .(model, N, df)])[N %in% c(10,45) & model %in% c(1,2,5)][order(model,N)]),
      row.names=FALSE)
cat("\nnota: per un effetto ENTRO cluster (time:trt) la regola between-within della letteratura\n")
cat("assegna i gdl al residuo entro-cluster, cioe' (numero di osservazioni) - (numero di cluster)\n")
cat("- (numero di covariate entro-cluster), NON N - p. Calcolo i due valori a confronto:\n")
att <- unique(R[set=="own", .(dgp,model,N,ntime,label,n_attempt)])
for (m in c(1,5)) for (nn in c(10,45)) for (T in c(4,10)) {
  # numero medio di osservazioni per replica: N * (0.7*T + 0.3*mean(1..T-1))
  nobs <- nn*(0.7*T + 0.3*mean(1:(T-1)))
  p <- c(4,5,4,4,5)[m]
  bw_within <- round(nobs - nn - 2)   # 2 = time e time:trt, covariate entro-cluster
  cat(sprintf("  M%d N=%d T=%d: N-p = %d ; between-within entro-cluster ~ nobs-N-2 = %d (nobs medio %.0f)\n",
    m, nn, T, nn-p, bw_within, nobs))
}

cat("\n############ 9. esistono conteggi PER CELLA pubblicati?\n")
tb <- list.files(path.expand("~/MMM/tables"), pattern="[.]csv$")
cat("tabelle in ~/MMM/tables:", length(tb), "\n")
for (f in tb) {
  h <- tryCatch(names(fread(file.path(path.expand("~/MMM/tables"), f), nrows=1)), error=function(e) character(0))
  has_cell <- all(c("model","ntime") %in% h) || all(c("model","N","ntime") %in% h)
  has_nonus <- any(grepl("n_attempt|n_usable|nonconv|converg|n_fail|n_no_se|n_big", h))
  if (has_cell || has_nonus)
    cat(sprintf("  %-32s colonne: %s\n", f, paste(h, collapse=",")))
}
cat("\nO_level_all.csv (la fonte dichiarata) ha n per cella ma non i NON utilizzabili:\n")
print(names(fread(path.expand("~/MMM/tables/O_level_all.csv"), nrows=1)))
cat("\ngeesmv_timeouts.csv:\n")
print(names(fread(path.expand("~/MMM/tables/geesmv_timeouts.csv"), nrows=1)))
print(head(as.data.frame(fread(path.expand("~/MMM/tables/geesmv_timeouts.csv"))), 4))
cat("\nDONE_MISC\n")
