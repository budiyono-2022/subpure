suppressMessages(library(subpure)); library(jsonlite)
S <- readRDS('sim_results.rds'); cond <- S$cond
rows <- list(); figd <- list()
for (cc in seq_along(S$res)) {
  p <- cond[cc,]; r <- as.data.frame(S$res[[cc]])
  pop <- pure_population(simulate_pure(n = 10, lambda_s = p$ls, lambda_w = p$lw, lambda_g = .6)$loadings)
  oh <- pop$global[['omega_h']]; vv <- mean(pop$facets[, 'var']); v1 <- pop$facets[1, 'var']; vcl <- mean(pop$facets[, 'var_cl'])
  ok <- r$ok == 1; pr <- ok & r$proper == 1; mc <- pr & r$mc == 1
  rows[[cc]] <- data.frame(n = p$n, lw = p$lw, ls = p$ls, reps = nrow(r),
    proper = mean(pr), se_ok = mean(mc),
    pop_oh = oh, est_oh = mean(r$omega_h[pr]), cov_oh = mean(r$omega_h_lo[mc] <= oh & r$omega_h_hi[mc] >= oh),
    blind_oh = mean(r$omega_h_blind[ok]),
    pop_var = vv, est_var = median(r$var[pr]), cov_var = mean(r$var_lo[mc] <= v1 & r$var_hi[mc] >= v1),
    reject = mean(r$reject[mc]), pop_var_cl = vcl, hab = mean(r$hab_added[ok]), fd_s = mean(r$fd_s[pr]))
  figd[[cc]] <- data.frame(cc = cc, n = p$n, lw = p$lw, ls = p$ls, oh = r$omega_h[pr], ob = r$omega_h_blind[pr], pop = oh)
}
T <- do.call(rbind, rows); print(format(round(T, 3)), row.names = FALSE)
saveRDS(list(T = T, figd = do.call(rbind, figd)), 'simsum.rds')
writeLines(toJSON(T, digits = 6, dataframe = 'rows'), 'simT.json')
