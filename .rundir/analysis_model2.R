# R/analysis_model2.R

analyze_model2_mmm <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
    mod <- mm(
      binY ~ time*trt + age,
      lv.formula = ~ 1,
      t.formula  = ~ 1,
      data = data,
      id   = id,
      q    = 20
    )
    s <- summary(mod)$mean.table
    data.frame(
      term     = rownames(s),
      estimate = s[,1],
      se       = s[,2],
      stat     = s[,3],
      p        = s[,4],
      method   = "mmm",
      N        = N,
      ntime    = ntime,
      label    = label
    )
}

analyze_model2_gee <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  data$id <- factor(data$id)
    mod <- geeglm(
      binY ~ time*trt + age,
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

analyze_model2_glmer <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  data$id <- factor(data$id)
    mod <- glmer(
      binY ~ time*trt + age + (1 | id),
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
