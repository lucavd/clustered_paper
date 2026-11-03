# s7d2_width_groups.R -- width ratio interval of the test / GEE AR(1) + Fay-Graubard + t under H1 by N group and T,
# on the same basis as prediction PG7 (exp of the mean over points of the per-point mean log ratio; tables/V4_S7d2_dp.csv).
suppressMessages(library(data.table)); D <- fread(path.expand('~/MMM/tables/V4_S7d2_dp.csv'))[label == 'mmlb' & is.finite(mlr)]
D[, Ng := cut(N, c(0, 10, 25, 65, 145), labels = c('10', '15-25', '35-65', '75-145'))]
W <- rbind(D[, .(group = 'all', points = .N, ratio = exp(mean(mlr)))], D[, .(points = .N, ratio = exp(mean(mlr))), by = .(group = paste0('N ', Ng))][order(group)],
           D[, .(points = .N, ratio = exp(mean(mlr))), by = .(group = paste0('T = ', ntime))][order(group)], D[small == TRUE, .(group = 'small designs', points = .N, ratio = exp(mean(mlr)))])
print(W, digits = 3); fwrite(W, path.expand('~/MMM/tables/V4_S7d2_width_groups.csv'))
