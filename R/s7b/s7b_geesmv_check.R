# =============================================================================
# s7b_geesmv_check.R -- S7b: independent check of the geesmv small-sample corrections on geepack::respiratory.
# Shows the `id` column trap of geesmv (GEE.var.* uses data$id whenever that column exists, ignoring the `id`
# argument) and compares, with the cluster column named `id`, the three corrections with clubSandwich (CR2 = KC,
# CR3 = MD) and with a hand-written Fay-Graubard sandwich on the geepack fit. Log: reports/morris/S7b_logs/s7b_geesmv_check.log
# =============================================================================

source(path.expand("~/MMM/R/s7b/s7b_common.R"))
suppressMessages({library(geepack); library(clubSandwich)})
data(respiratory, package = "geepack"); d <- as.data.frame(respiratory)
d$pid <- factor(paste(d$center, d$id, sep = "_")); d$treat <- factor(as.character(d$treat), levels = c("P", "A")); d$sex <- factor(d$sex); d$center <- factor(d$center)
d <- d[order(d$pid, d$visit), ]; d$pidn <- as.integer(d$pid)
cat("--- with the original respiratory `id` column present (within-centre id), geesmv keyed on it ---\n")
f <- outcome ~ treat + baseline + visit + age + sex + center
g <- geeglm(f, family = binomial, id = pidn, corstr = "exchangeable", data = d)
cat("geepack alpha:", g$geese$alpha, " scale:", g$geese$gamma, "\n")
for (ty in c("CR0","CR1","CR2","CR3")) cat(ty, sqrt(diag(as.matrix(vcovCR(g, cluster = d$pidn, type = ty))))[2], "\n")
# manual sandwich for exchangeable GEE with the geepack alpha: V_i = A^1/2 R A^1/2, B = sum D' V^-1 D, M = sum D' V^-1 e e' V^-1 D
X <- model.matrix(f, d); mu <- fitted(g); e <- d$outcome - mu; ids <- split(seq_len(nrow(d)), d$pidn)
alpha <- g$geese$alpha; B <- 0; Ms <- list(); Hs <- list(); Dl <- list(); Vl <- list(); el <- list()
for (k in seq_along(ids)) { i <- ids[[k]]; a <- mu[i] * (1 - mu[i]); Dk <- X[i, , drop = FALSE] * a; Rk <- matrix(alpha, length(i), length(i)); diag(Rk) <- 1
  Vk <- sqrt(a) %o% sqrt(a) * Rk; Vi <- solve(Vk); B <- B + t(Dk) %*% Vi %*% Dk; Dl[[k]] <- Dk; Vl[[k]] <- Vi; el[[k]] <- e[i] }
Bi <- solve(B); M0 <- 0; Mmd <- 0; Mkc <- 0; Mfg <- 0; b <- 0.75
for (k in seq_along(ids)) { Dk <- Dl[[k]]; Vi <- Vl[[k]]; ek <- el[[k]]; Hk <- Dk %*% Bi %*% t(Dk) %*% Vi; I <- diag(nrow(Dk))
  M0 <- M0 + t(Dk) %*% Vi %*% ek %*% t(ek) %*% Vi %*% Dk
  r1 <- solve(I - Hk) %*% ek; Mmd <- Mmd + t(Dk) %*% Vi %*% r1 %*% t(r1) %*% Vi %*% Dk
  Qk <- t(Dk) %*% Vi %*% Dk %*% Bi; qk <- pmin(b, diag(Qk)); Ak <- diag(1 / sqrt(1 - qk)); Mfg <- Mfg + Ak %*% t(Dk) %*% Vi %*% ek %*% t(ek) %*% Vi %*% Dk %*% Ak }
cat("manual CR0 SE:", sqrt(diag(Bi %*% M0 %*% Bi))[2], " MD:", sqrt(diag(Bi %*% Mmd %*% Bi))[2], " FG:", sqrt(diag(Bi %*% Mfg %*% Bi))[2], "\n")
v <- geesmv::GEE.var.fg(f, id = "pidn", family = binomial, data = d, corstr = "exchangeable"); cat("geesmv FG:", sqrt(v$cov.beta)[2], "\n")
v <- geesmv::GEE.var.md(f, id = "pidn", family = binomial, data = d, corstr = "exchangeable"); cat("geesmv MD:", sqrt(v$cov.beta)[2], "\n")
# geesmv with the design as in the simulation (numeric covariates only)
d$treatn <- as.integer(d$treat == "A"); d$sexn <- as.integer(d$sex == "M"); d$centern <- as.integer(d$center == "2")
f2 <- outcome ~ treatn + baseline + visit + age + sexn + centern
v <- geesmv::GEE.var.fg(f2, id = "pidn", family = binomial, data = d, corstr = "exchangeable"); cat("geesmv FG numeric design:", sqrt(v$cov.beta)[2], "\n")
f3 <- outcome ~ treatn + baseline
v <- geesmv::GEE.var.fg(f3, id = "pidn", family = binomial, data = d, corstr = "exchangeable"); cat("geesmv FG reduced:", sqrt(v$cov.beta)[2], "\n")
g3 <- geeglm(f3, family = binomial, id = pidn, corstr = "exchangeable", data = d); cat("geepack reduced:", summary(g3)$coefficients[2, 2], "\n")

cat("--- cluster column renamed to `id` ---\n")
d$id <- d$pidn
for (nm in c("md", "kc", "fg")) { v <- do.call(getExportedValue("geesmv", paste0("GEE.var.", nm)), list(f, id = "id", family = binomial, data = d, corstr = "exchangeable")); cat("geesmv", nm, "treatA SE:", sqrt(v$cov.beta)[2], "\n") }
