.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
sm <- function(a, b) if (is.numeric(a)) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(a))) else (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
dif <- function(f) { x <- fread(file.path("tables/variant_nc", f)); y <- fread(file.path("tables/variant_cv", f)); d <- Reduce(`|`, lapply(names(x), function(cl) !sm(x[[cl]], y[[cl]]))); list(x = x, y = y, d = d) }
z <- dif("V4_perf_all.csv.gz"); X <- z$x[z$d & !grepl("mmm", z$x$method)]; cat("perf_all outside mmm by set:\n"); print(X[, .N, by = .(set, usab)])
z <- dif("V4_S6_table_cells.csv"); print(names(z$x)); X <- z$x[z$d]; print(X[, .N, by = .(pre = sub("\\|.*", "", get(names(X)[1])))])
z <- dif("V4_S5_values.csv"); print(names(z$x)); X <- z$x[z$d & !grepl("mmm|MMM|marginal", do.call(paste, z$x[, 1:3]))]; print(head(X[, 1:4], 8)); print(X[, .N, by = .(substr(get(names(X)[1]), 1, 4))])
x <- fread("tables/variant_nc/V4_S6_m1c_unusable_reasons.csv"); y <- fread("tables/variant_cv/V4_S6_m1c_unusable_reasons.csv"); print(fsetdiff(y, x)); print(fsetdiff(x, y))
z <- dif("V4_S5_claims.csv"); print(z$x[z$d, .(name, class)][, cls_nc := class][, cls_cv := z$y$class[z$d]][cls_nc != cls_cv])
