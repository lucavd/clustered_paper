#!/usr/bin/env Rscript
# S7f registry builder (D-38): every number of the v5 manuscript and supplement (copy of R/s7/s7_numbers.R, paths changed).
# Input : manuscript/v5/numbers_spec.csv (one row per number: id, file, filter as JSON, field, op, scale,
#         abs, fmt, kind, mcse_field), produced from the tokens of the manuscript source.
# Output: tables/V5_S7f_numbers.csv = spec + value (raw), mcse, text (as printed) and md5 of the source file.
# The independent check is tools/s7_verify_registry.py (G3). Read-only on every source file.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(jsonlite)}); setDTthreads(1L)
setwd("/home/user/MMM")
S <- as.data.table(read.csv("manuscript/v5/numbers_spec.csv", colClasses = "character", na.strings = character(0), check.names = FALSE))
cache <- list()
getf <- function(f) { if (is.null(cache[[f]])) cache[[f]] <<- fread(f, colClasses = "character", na.strings = c("NA", "")); cache[[f]] }
isnum <- function(x) suppressWarnings(!is.na(as.numeric(x)))
match_eq <- function(col, v) {
  if (v == "") return(is.na(col) | col == "")
  if (isnum(v)) { cn <- suppressWarnings(as.numeric(col)); return(!is.na(cn) & cn == as.numeric(v)) }
  !is.na(col) & col == v
}
select_rows <- function(D, filt) {
  keep <- rep(TRUE, nrow(D))
  for (k in names(filt)) {
    v <- as.character(filt[[k]])
    if (grepl("__lt$", k)) { col <- D[[sub("__lt$", "", k)]]; keep <- keep & !is.na(col) & as.numeric(col) < as.numeric(v) }
    else if (grepl("__gt$", k)) { col <- D[[sub("__gt$", "", k)]]; keep <- keep & !is.na(col) & as.numeric(col) > as.numeric(v) }
    else if (grepl("__in$", k)) { col <- as.numeric(D[[sub("__in$", "", k)]]); keep <- keep & !is.na(col) & col %in% as.numeric(strsplit(v, ",")[[1]]) }
    else { stopifnot(k %in% names(D)); keep <- keep & match_eq(D[[k]], v) }
  }
  D[keep]
}
sig2 <- function(v, mcse) {
  if (v == 0) return("0")
  e <- floor(log10(abs(v))); d <- max(0, 1 - e)
  if (is.finite(mcse) && mcse > 0) d <- min(d, max(0, -floor(log10(mcse))))
  s <- sprintf("%.*f", d, v)
  if (abs(as.numeric(s)) >= 10^(e + 1) && d > 0) s <- sprintf("%.*f", d - 1, v)
  s
}
fmt_int <- function(v) formatC(round(v), format = "d", big.mark = ",")
fmt_tab <- function(v, mcse) { d <- if (is.finite(mcse) && mcse >= 0.5) 0 else 1; sprintf("%.*f", d, v) }
out <- vector("list", nrow(S))
for (i in seq_len(nrow(S))) {
  s <- S[i]; D <- getf(s$file); filt <- fromJSON(s$filter); R <- select_rows(D, filt)
  op <- s$op; mc <- NA_real_
  if (op == "count") val <- nrow(R)
  else {
    x <- R[[s$field]]; stopifnot(!is.null(x))
    if (op == "value") { if (nrow(R) != 1) stop(sprintf("%s: %d rows", s$id, nrow(R))); val <- x
      if (nzchar(s$mcse_field)) mc <- as.numeric(R[[s$mcse_field]]) }
    else { xn <- as.numeric(x); stopifnot(length(xn) > 0, all(!is.na(xn)))
      val <- switch(op, min = min(xn), max = max(xn), sum = sum(xn), stop("op")) }
  }
  if (s$kind == "text") { txt <- as.character(val); raw <- txt }
  else {
    v <- as.numeric(val) * as.numeric(s$scale); if (s$abs == "TRUE") v <- abs(v)
    mc <- mc * as.numeric(s$scale)
    txt <- if (s$kind == "int") fmt_int(v) else if (s$kind == "plain") sprintf("%.15g", v) else if (s$fmt == "tab") fmt_tab(v, mc) else sig2(v, mc)
    raw <- sprintf("%.17g", v)
  }
  out[[i]] <- data.table(id = s$id, value = raw, mcse = if (is.finite(mc)) sprintf("%.17g", mc) else "", text = txt,
                         md5 = unname(tools::md5sum(s$file)))
}
REG <- cbind(S, rbindlist(out)[, -"id"])
fwrite(REG, "tables/V5_S7f_numbers.csv")
cat(sprintf("registry: %d numbers from %d files -> tables/V5_S7f_numbers.csv\n", nrow(REG), uniqueN(REG$file)))
