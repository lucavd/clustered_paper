# R/analysis_model4.R

analyze_model4_mmm <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  mod <- mm(
    binY ~ time*trt,
    lv.formula = ~ time,
    t.formula  = ~ trt,
    data = data,
    id   = id,
    q    = 20
  )
  s <- summary(mod)$mean.table
  df_mean <- data.frame(
    term     = rownames(s),
    estimate = s[,1],
    se       = s[,2],
    stat     = s[,3],
    p        = s[,4]
  )
  a <- summary(mod)$assoc.table
  df_assoc <- data.frame(
    term     = rownames(a),
    estimate = a[,1],
    se       = a[,2],
    stat     = a[,3],
    p        = a[,4]
  )
  df_final <- rbind(df_mean, df_assoc)
  df_final$method <- "mmm"
  df_final$N      <- N
  df_final$ntime  <- ntime
  df_final$label  <- label
  df_final
}

analyze_model4_gee <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  data$id <- factor(data$id)
  mod <- geeglm(
    binY ~ time*trt,
    id = id,
    data = data,
    family = binomial,
    corstr = "exchangeable",
    scale.fix = TRUE
  )
  s <- summary(mod)$coefficients
  data.frame(
    term     = rownames(s),
    estimate = s[,1],
    se       = s[,2],
    stat     = s[,3],
    p        = s[,4],
    method   = "gee",
    N        = N,
    ntime    = ntime,
    label    = label
  )
}

analyze_model4_glmer <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  data$id <- factor(data$id)
  mod <- glmer(
    binY ~ time*trt + (1 + time | id),
    data = data,
    family = binomial
  )
  s <- summary(mod)$coefficients
  D <- VarCorr(mod)$id[[1]]
  data.frame(
    term     = rownames(s),
    estimate = s[,1],
    se       = s[,2],
    stat     = s[,3],
    p        = s[,4],
    method   = "glmer",
    D        = D,
    N        = N,
    ntime    = ntime,
    label    = label
  )
}
