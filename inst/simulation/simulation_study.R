suppressMessages(library(subpure)); library(parallel)
cond <- expand.grid(n = c(500, 2000), lw = c(.25, .5), ls = c(.15, .25, .40))
reps <- 50
haberman <- function(X, facets) { # observed-score classical PRMSE ratio (Haberman, 2008)
  alpha <- function(M) { k <- ncol(M); k/(k-1)*(1-sum(apply(M,2,var))/var(rowSums(M))) }
  Tt <- rowSums(X)
  sapply(facets, function(v){ s <- rowSums(X[v]); a <- alpha(X[v]); vs <- var(s); vt <- var(Tt)
    cov_tsT <- a*vs + cov(s, Tt-s); px <- cov_tsT^2/(a*vs*vt); a/px })
}
nm <- c('omega_h','omega_h_lo','omega_h_hi','omega_h_blind','var','var_lo','var_hi','var1',
        'reject','fd_s','proper','mc','hab_added','ok')
one_rep <- function(i, cc) {
  set.seed(1e5*cc + i); p <- cond[cc,]
  sim <- simulate_pure(n = p$n, lambda_s = p$ls, lambda_w = p$lw, lambda_g = .6)
  y <- setNames(rep(NA_real_, length(nm)), nm); y['ok'] <- 0
  r <- tryCatch(suppressWarnings(pure_value(sim$data, sim$facets, sim$method, R = 100)), error = function(e) NULL)
  if (is.null(r)) return(y)
  f <- r$facets; g <- r$global
  y[] <- c(g$omega_h, g$omega_h_lo, g$omega_h_hi, g$omega_h_blind, mean(f$var), f$var_lo[1], f$var_hi[1],
           f$var[1], mean(f$p_added < .05), mean(f$fd_s), g$proper, g$intervals == 'mc',
           mean(haberman(sim$data, sim$facets) > 1), 1)
  y
}
res <- if (file.exists('sim_results.rds')) readRDS('sim_results.rds')$res else list()
for (cc in seq_len(nrow(cond))) {
  if (length(res) >= cc && !is.null(res[[cc]])) next
  res[[cc]] <- do.call(rbind, mclapply(1:reps, one_rep, cc = cc, mc.cores = 2))
  saveRDS(list(cond = cond, res = res), 'sim_results.rds')
  cat(cc, date(), '\n')
}
