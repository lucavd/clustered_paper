# =============================================================================
# s3_build_v4_code.R -- builds R/v4/rundir/, the code executed by S3, from the
# code that produced results_full/ (.rundir/). Every edit asserts the original
# line before replacing it. Edits:
#   D-01  GenBinaryY M3-M5: H0 sigma/gamma set to the H1 values (same edits as
#         the S2 pilot, R/s2/pilot_gen/).
#   D-02  GLMM-AR1 M5: D[3,3] 0.01 -> 0.001 (as GLMM-native M5).
#   D-08  do_analysis: glmer_lrt uses the random structure of the Wald glmer;
#         lrt_re_structure() is exported so rows can be tagged.
#   D-13/D-14  scenario M1c generators for the three families (new files).
# =============================================================================
REPO <- path.expand("~/MMM"); SRC <- file.path(REPO, ".rundir"); V4 <- file.path(REPO, "R/v4/rundir")
NEWSRC <- file.path(REPO, "R/v4/src")
unlink(V4, recursive = TRUE); dir.create(V4, recursive = TRUE)
for (f in list.files(SRC, recursive = TRUE)) {
  dir.create(dirname(file.path(V4, f)), showWarnings = FALSE, recursive = TRUE)
  stopifnot(file.copy(file.path(SRC, f), file.path(V4, f)))
}
edit <- function(f, edits, header) {
  p <- file.path(V4, f); L <- readLines(p, encoding = "UTF-8")
  for (e in edits) {
    k <- as.integer(e[[1]]); stopifnot(startsWith(L[k], e[[2]]))  # same rule as the S2 pilot
    L[k] <- e[[3]]
  }
  writeLines(c(header, L), p, useBytes = TRUE)
}
H01 <- "# V4 (S3, D-01): copy of .rundir/%s with the H0 sigma/gamma set to the H1 values."
edit("generation_model3.R", list(list(34, "    manual = c(1, 0))", "    manual = c(1, 0.01))")), sprintf(H01, "generation_model3.R"))
edit("generation_model4.R", list(list(35, "    manual = c(1, 0) )", "    manual = c(1, 0.01) ) # D-01: H0 sigma = H1 sigma"),
                                 list(40, "    manual = c(1, 0) )", "    manual = c(1, 0.01) )")), sprintf(H01, "generation_model4.R"))
edit("generation_model5.R", list(list(42, "    manual = c(1, 0, 0))", "    manual = c(1, 0.01, 0.01))"),
                                 list(47, "     manual = c(1, 0, 0))", "     manual = c(1, 0.1, 0.05))")), sprintf(H01, "generation_model5.R"))
# D-02
p <- file.path(V4, "ar1/generation_models_ar1.R"); L <- readLines(p, encoding = "UTF-8")
old <- "                0.01, 0.0,   0.01), nrow = 3, byrow = TRUE)"
k <- which(L == old); stopifnot(length(k) == 1L)
L[k] <- "                0.01, 0.0,   0.001), nrow = 3, byrow = TRUE)  # V4 (S3, D-02): was 0.01, now as GLMM-native M5"
writeLines(c("# V4 (S3, D-02): copy of .rundir/ar1/generation_models_ar1.R with D[3,3] of M5 = 0.001.", L), p, useBytes = TRUE)
# D-08
p <- file.path(V4, "do_analysis.R"); L <- readLines(p, encoding = "UTF-8")
old <- "  f_glmer <- as.formula(paste(\"binY ~\", mean_rhs, \"+ (1|id)\"))"
k <- which(L == old); stopifnot(length(k) == 1L)
L[k] <- "  f_glmer <- as.formula(paste(\"binY ~\", mean_rhs, \"+\", lrt_re_structure(model)))  # V4 (S3, D-08)"
L <- c("# V4 (S3, D-08): copy of .rundir/do_analysis.R; glmer_lrt uses the random structure of the Wald glmer.",
       "lrt_re_structure <- function(model) switch(as.character(model),",
       "  \"1\" = \"(1|id)\", \"2\" = \"(1|id)\", \"3\" = \"(1 + time|id)\", \"4\" = \"(1 + time|id)\",",
       "  \"5\" = \"(1 + time + dose|id)\")", "", L)
writeLines(L, p, useBytes = TRUE)
# D-13 / D-14: M1c
stopifnot(file.copy(file.path(NEWSRC, "m1c_params.R"), file.path(V4, "m1c_params.R")),
          file.copy(file.path(NEWSRC, "generation_model1c.R"), file.path(V4, "generation_model1c.R")),
          file.copy(file.path(NEWSRC, "manual_generation_model1c.R"), file.path(V4, "manual/generation_model1c.R")),
          file.copy(file.path(NEWSRC, "ar1_generation_model1c.R"), file.path(V4, "ar1/generation_model1c.R")))
# check: GenBinaryY v4 files equal the S2 pilot copies except the header line
for (f in sprintf("generation_model%d.R", 3:5)) {
  a <- readLines(file.path(V4, f))[-1]; b <- readLines(file.path(REPO, "R/s2/pilot_gen", f))[-1]
  cat(sprintf("%s identical to S2 pilot (body): %s\n", f, identical(a, b)))
}
d <- tools::md5sum(file.path(V4, list.files(V4, recursive = TRUE)))
o <- tools::md5sum(file.path(SRC, list.files(SRC, recursive = TRUE)))
names(d) <- sub(paste0(V4, "/"), "", names(d)); names(o) <- sub(paste0(SRC, "/"), "", names(o))
chg <- names(d)[!names(d) %in% names(o) | d[names(d)] != o[names(d)]]
cat("files changed or added vs .rundir:\n"); print(chg)
writeLines(sprintf("%s  %s", d, names(d)), file.path(REPO, "R/v4/rundir_md5.txt"))
