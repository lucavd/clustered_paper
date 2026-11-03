# setup_lib.R -- flipscores (PPM binary) and fewclust 0.1.0 (~/MMM-package, commit recorded) into ~/R/s7c_lib
lib <- path.expand("~/R/s7c_lib"); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, "~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
        HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)),
        install.packages.compile.from.source = "never")
if (!requireNamespace("flipscores", quietly = TRUE)) install.packages("flipscores", lib = lib)
install.packages(path.expand("~/MMM-package/fewclust"), repos = NULL, type = "source", lib = lib)
cat("flipscores", as.character(packageVersion("flipscores", lib.loc = lib)), "fewclust", as.character(packageVersion("fewclust", lib.loc = lib)), "\n")
