## Internal: indices from standardized loadings -------------------------------

# lG: general loadings (n); lS: n x K facet-specific loadings (zero outside facet);
# lW: wording loadings (n, zero for non-method items); th: residual variances (n);
# fidx: list of item positions per facet.
pure_indices <- function(lG, lS, lW, th, fidx) {
  n <- length(lG)
  one <- rep(1, n)
  cv <- function(a, b) {
    sum(a * lG) * sum(b * lG) + sum(crossprod(a, lS) * crossprod(b, lS)) +
      sum(a * lW) * sum(b * lW) + sum(a * b * th)
  }
  varX <- cv(one, one)
  glob <- c(omega_t = 1 - sum(th) / varX,
            omega_h = sum(lG)^2 / varX,
            omega_w = sum(lW)^2 / varX,
            ecv_g = sum(lG^2) / (sum(lG^2) + sum(lS^2) + sum(lW^2)))
  fac <- t(vapply(seq_along(fidx), function(k) {
    e <- numeric(n); e[fidx[[k]]] <- 1
    tG <- sum(e * lG); tS <- sum(e * lS[, k]); tW <- sum(e * lW)
    vXk <- cv(e, e); cXkX <- cv(e, one)
    covT <- function(a, withW) sum(a * lG) * tG + sum(a * lS[, k]) * tS +
      if (withW) sum(a * lW) * tW else 0
    res <- numeric(0)
    for (withW in c(FALSE, TRUE)) {
      vT <- tG^2 + tS^2 + if (withW) tW^2 else 0
      cs <- covT(e, withW); ct <- covT(one, withW)
      p_s <- cs^2 / (vXk * vT); p_x <- ct^2 / (varX * vT)
      S <- matrix(c(vXk, cXkX, cXkX, varX), 2)
      cc <- c(cs, ct)
      p_a <- drop(crossprod(cc, solve(S, cc))) / vT
      res <- c(res, p_s, p_x, p_a, p_s / p_x)
    }
    c(res, omega_hs = tS^2 / vXk, method_share = tW^2 / vXk)
  }, numeric(10)))
  colnames(fac) <- c("prmse_sub", "prmse_tot", "prmse_aug", "var",
                     "prmse_sub_cl", "prmse_tot_cl", "prmse_aug_cl", "var_cl",
                     "omega_hs", "method_share")
  rownames(fac) <- names(fidx)
  list(global = glob, facets = fac)
}

## Internal: loadings from a fitted lavaan model -------------------------------

extract_loadings <- function(fit, items, facets, method, ordered, par = NULL) {
  pe <- lavaan::parameterTable(fit)
  pe <- pe[pe$op == "=~", ]
  est <- if (is.null(par)) pe$est else {
    v <- pe$est; free <- pe$free > 0; v[free] <- par[pe$free[free]]; v
  }
  n <- length(items); K <- length(facets)
  lG <- lW <- stats::setNames(numeric(n), items)
  lS <- matrix(0, n, K, dimnames = list(items, names(facets)))
  for (r in seq_len(nrow(pe))) {
    lv <- pe$lhs[r]; it <- pe$rhs[r]
    if (lv == "G") lG[it] <- est[r] else if (lv == "W") lW[it] <- est[r] else
      lS[it, sub("^S_", "", lv)] <- est[r]
  }
  com <- lG^2 + rowSums(lS^2) + lW^2
  if (ordered) {
    th <- 1 - com
  } else {
    pt <- lavaan::parameterTable(fit)
    rv <- pt[pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% items, ]
    v <- if (is.null(par)) rv$est else {
      x <- rv$est; f <- rv$free > 0; x[f] <- par[rv$free[f]]; x
    }
    th <- stats::setNames(v, rv$lhs)[items]
    s <- sqrt(com + th)
    lG <- lG / s; lS <- lS / s; lW <- lW / s; th <- th / s^2
  }
  list(lG = lG, lS = lS, lW = lW, th = th)
}

fit_one <- function(data, facets, method, ordered, structure, ...) {
  items <- unlist(facets, use.names = FALSE)
  args <- list(model = pure_syntax(facets, method, structure), data = data,
               std.lv = TRUE, orthogonal = TRUE, ...)
  if (ordered) { args$ordered <- items; args$estimator <- "WLSMV" }
  fit <- suppressWarnings(do.call("cfa", args, envir = asNamespace("lavaan")))
  if (!lavaan::lavInspect(fit, "converged"))
    stop("The ", structure, " model did not converge.", call. = FALSE)
  fit
}

determinacy <- function(L) {
  Lam <- cbind(G = L$lG, L$lS, W = L$lW)
  Sig <- tcrossprod(Lam) + diag(L$th, length(L$th))
  fd <- tryCatch(sqrt(pmax(diag(crossprod(Lam, solve(Sig, Lam))), 0)), error = function(e) rep(NA_real_, ncol(Lam)))
  stats::setNames(fd, c("G", colnames(L$lS), "W"))
}

mc_draws <- function(fit, items, facets, method, ordered, R) {
  fidx <- lapply(facets, function(v) match(v, items))
  V <- tryCatch(suppressWarnings(lavaan::lavInspect(fit, "vcov")), error = function(e) NULL)
  ok <- !is.null(V) && all(is.finite(V)) && min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > 0
  if (!ok) return(NULL)
  b <- lavaan::coef(fit)
  P <- MASS::mvrnorm(R, b, V)
  lapply(seq_len(R), function(r) {
    L <- extract_loadings(fit, items, facets, method, ordered, par = P[r, ])
    pure_indices(L$lG, L$lS, L$lW, pmax(L$th, 1e-6), fidx)
  })
}

boot_draws <- function(data, fit, items, facets, method, ordered, B, ...) {
  fidx <- lapply(facets, function(v) match(v, items))
  out <- vector("list", B)
  for (b in seq_len(B)) {
    db <- data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
    fb <- tryCatch(fit_one(db, facets, method, ordered, "bifactor_method", start = fit, ...),
                   error = function(e) NULL)
    if (is.null(fb)) next
    L <- extract_loadings(fb, items, facets, method, ordered)
    out[[b]] <- pure_indices(L$lG, L$lS, L$lW, pmax(L$th, 1e-6), fidx)
  }
  out[!vapply(out, is.null, TRUE)]
}

#' Method-purified added value of facet subscores
#'
#' Fits a bifactor model with a wording factor (and, for comparison, the same
#' bifactor model without it) and computes, for each facet, the proportional
#' reduction in mean squared error (PRMSE) obtained when the facet's
#' method-purified true score is predicted from the facet subscore, the total
#' score, or their best linear combination (augmented score). The purified
#' target contains general and facet-specific variance but no wording variance.
#'
#' Composites are unit-weighted sums on the latent-response metric when
#' \code{ordered = TRUE}, and on the standardized item metric otherwise.
#' Uncertainty is propagated by Monte Carlo draws from the asymptotic sampling
#' distribution of the model parameters. The value-added ratio is
#' \code{prmse_sub / prmse_tot}; values above 1 indicate that the subscore
#' predicts its own true score better than the total score does.
#'
#' @param data Data frame with reverse-scored items.
#' @param facets Named list of item-name vectors.
#' @param method Character vector of wording-factor items.
#' @param group Optional grouping variable; models are fitted separately in each
#'   group (configural) and between-group contrasts are reported.
#' @param ordered Treat items as ordinal (WLSMV)? Default \code{TRUE}.
#' @param R Number of Monte Carlo draws, or bootstrap resamples when
#'   \code{ci = "boot"} (default 1000; 0 skips intervals).
#' @param ci Interval method: \code{"mc"} (Monte Carlo from the asymptotic
#'   covariance matrix, default) or \code{"boot"} (nonparametric bootstrap with
#'   refitting). When a facet-specific factor is very weak the information
#'   matrix can be singular (empirical underidentification); Monte Carlo
#'   intervals are then unavailable, a warning is issued, and
#'   \code{ci = "boot"} can be used instead.
#' @param level Confidence level.
#' @param seed Optional seed for the Monte Carlo draws.
#' @param ... Passed to \code{lavaan::cfa()}.
#' @return An object of class \code{pure_value} with elements \code{facets}
#'   (per-group, per-facet indices with intervals, the factor determinacy
#'   \code{fd_s} of each facet-specific factor and the share \code{neg_share}
#'   of its items with specific loadings below -0.10), \code{global} (including
#'   determinacies of the general and wording factors, whether the solution is
#'   proper, and which interval method was used),
#'   \code{misattribution} (method-blind minus purified facet omega),
#'   \code{contrasts} (if \code{group} is given) and the fitted models.
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_pure(n = 600)
#' pv <- pure_value(sim$data, sim$facets, sim$method, R = 200)
#' pv
#' }
#' @export
pure_value <- function(data, facets, method, group = NULL, ordered = TRUE,
                       R = 1000, ci = c("mc", "boot"), level = 0.95, seed = NULL, ...) {
  ci <- match.arg(ci)
  check_facets(facets, method)
  if (length(method) < 2) stop("`method` must list at least two items.", call. = FALSE)
  items <- unlist(facets, use.names = FALSE)
  if (!is.null(seed)) set.seed(seed)
  groups <- if (is.null(group)) list(All = data) else {
    if (!group %in% names(data)) stop("`group` not found in `data`.", call. = FALSE)
    g <- data[[group]]
    split(data, droplevels(as.factor(g)))
  }
  a <- (1 - level) / 2
  fidx <- lapply(facets, function(v) match(v, items))
  per <- lapply(names(groups), function(gn) {
    dg <- groups[[gn]]
    fit_m <- fit_one(dg, facets, method, ordered, "bifactor_method", ...)
    fit_b <- fit_one(dg, facets, method, ordered, "bifactor", ...)
    Lm <- extract_loadings(fit_m, items, facets, method, ordered)
    Lb <- extract_loadings(fit_b, items, facets, method, ordered)
    pm <- pure_indices(Lm$lG, Lm$lS, Lm$lW, Lm$th, fidx)
    proper <- all(Lm$th > 0)
    fdm <- if (proper) determinacy(Lm) else stats::setNames(rep(NA_real_, length(facets) + 2),
                                                    c("G", names(facets), "W"))
    negsh <- vapply(names(facets), function(k) {
      v <- Lm$lS[facets[[k]], k]; mean(v < -0.1) }, 0)
    if (!proper) warning("Group '", gn, "': improper (Heywood) solution in the bifactor-plus-wording ",
      "model; facet indices are not interpretable.", call. = FALSE)
    pb <- pure_indices(Lb$lG, Lb$lS, Lb$lW, Lb$th, fidx)
    draws <- NULL; how <- "none"
    if (R > 0) {
      if (ci == "mc") {
        draws <- mc_draws(fit_m, items, facets, method, ordered, R)
        if (is.null(draws)) warning("Group '", gn, "': the information matrix of the bifactor-plus-wording ",
          "model could not be inverted (likely a near-empty facet-specific factor). Point estimates ",
          "are reported without intervals; consider ci = \"boot\".", call. = FALSE) else how <- "mc"
      } else {
        draws <- boot_draws(dg, fit_m, items, facets, method, ordered, R, ...)
        if (length(draws) < 2) draws <- NULL else how <- "boot"
      }
    }
    list(group = gn, n = nrow(dg), fit_method = fit_m, fit_blind = fit_b,
         point = pm, blind = pb, draws = draws, how = how, proper = proper,
         fd = fdm, negsh = negsh)
  })
  names(per) <- names(groups)
  q <- function(x) stats::quantile(x, c(a, 1 - a), names = FALSE, na.rm = TRUE)
  fac <- do.call(rbind, lapply(per, function(p) {
    f <- p$point$facets
    out <- data.frame(group = p$group, n = p$n, facet = rownames(f),
                      f[, c("prmse_sub", "prmse_tot", "prmse_aug", "var",
                            "var_cl", "omega_hs", "method_share")],
                      omega_hs_blind = p$blind$facets[, "omega_hs"],
                      fd_s = unname(p$fd[rownames(f)]), neg_share = unname(p$negsh),
                      row.names = NULL, check.names = FALSE)
    if (!is.null(p$draws)) {
      arr <- simplify2array(lapply(p$draws, `[[`, "facets"))
      ci <- function(col) t(apply(arr[, col, , drop = FALSE], 1, q))
      out$var_lo <- ci("var")[, 1]; out$var_hi <- ci("var")[, 2]
      out$p_added <- apply(arr[, "var", , drop = FALSE], 1, function(x) mean(x <= 1))
      out$omega_hs_lo <- ci("omega_hs")[, 1]; out$omega_hs_hi <- ci("omega_hs")[, 2]
    } else if (R > 0) {
      out$var_lo <- out$var_hi <- out$p_added <- out$omega_hs_lo <- out$omega_hs_hi <- NA_real_
    }
    out
  }))
  glob <- do.call(rbind, lapply(per, function(p) {
    out <- data.frame(group = p$group, n = p$n, t(p$point$global),
                      omega_h_blind = p$blind$global[["omega_h"]], fd_g = p$fd[["G"]],
                      fd_w = p$fd[["W"]], proper = p$proper, intervals = p$how,
                      check.names = FALSE)
    if (!is.null(p$draws)) {
      gm <- vapply(p$draws, function(d) d$global, numeric(4))
      out$omega_h_lo <- q(gm["omega_h", ])[1]; out$omega_h_hi <- q(gm["omega_h", ])[2]
      out$omega_w_lo <- q(gm["omega_w", ])[1]; out$omega_w_hi <- q(gm["omega_w", ])[2]
    } else if (R > 0) {
      out$omega_h_lo <- out$omega_h_hi <- out$omega_w_lo <- out$omega_w_hi <- NA_real_
    }
    out
  }))
  rownames(glob) <- NULL
  mis <- fac[, c("group", "facet")]
  mis$omega_hs_blind <- fac$omega_hs_blind
  mis$omega_hs_pure <- fac$omega_hs
  mis$misattributed <- fac$omega_hs_blind - fac$omega_hs
  mis$omega_h_blind_minus_pure <- glob$omega_h_blind[match(mis$group, glob$group)] -
    glob$omega_h[match(mis$group, glob$group)]
  con <- NULL
  ok_draws <- vapply(per, function(p) !is.null(p$draws), TRUE)
  if (length(per) > 1 && R > 0 && sum(ok_draws) > 1) {
    per_all <- per; per <- per[ok_draws]
    pr <- utils::combn(names(per), 2, simplify = FALSE)
    con <- do.call(rbind, lapply(pr, function(z) {
      d1 <- per[[z[1]]]$draws; d2 <- per[[z[2]]]$draws
      m <- min(length(d1), length(d2))
      res <- lapply(c("var", "omega_hs"), function(ix) {
        a1 <- vapply(d1[seq_len(m)], function(d) d$facets[, ix], numeric(length(facets)))
        a2 <- vapply(d2[seq_len(m)], function(d) d$facets[, ix], numeric(length(facets)))
        dd <- matrix(a1 - a2, nrow = length(facets))
        pt <- per[[z[1]]]$point$facets[, ix] - per[[z[2]]]$point$facets[, ix]
        data.frame(contrast = paste(z, collapse = " - "), facet = names(facets),
                   index = ix, diff = pt,
                   lo = apply(dd, 1, q)[1, ], hi = apply(dd, 1, q)[2, ],
                   p = apply(dd, 1, function(x) 2 * min(mean(x <= 0), mean(x >= 0))),
                   row.names = NULL)
      })
      do.call(rbind, res)
    }))
    hg <- lapply(per, function(p) vapply(p$draws, function(d) d$global[["omega_h"]], 0))
    for (zz in pr) {
      m <- min(length(hg[[zz[1]]]), length(hg[[zz[2]]]))
      dd <- hg[[zz[1]]][seq_len(m)] - hg[[zz[2]]][seq_len(m)]
      con <- rbind(con, data.frame(contrast = paste(zz, collapse = " - "), facet = "(total)",
        index = "omega_h",
        diff = per[[zz[1]]]$point$global[["omega_h"]] - per[[zz[2]]]$point$global[["omega_h"]],
        lo = q(dd)[1], hi = q(dd)[2], p = 2 * min(mean(dd <= 0), mean(dd >= 0))))
    }
    rownames(con) <- NULL
    per <- per_all
  }
  out <- list(facets = fac, global = glob, misattribution = mis, contrasts = con,
              models = lapply(per, function(p) list(method = p$fit_method, blind = p$fit_blind)),
              level = level, R = R, ci = ci, ordered = ordered)
  class(out) <- "pure_value"
  out
}

#' @export
print.pure_value <- function(x, digits = 3, ...) {
  r <- function(d) { num <- vapply(d, is.numeric, TRUE); d[num] <- lapply(d[num], round, digits); d }
  cat("Method-purified subscore value (", if (x$ordered) "ordinal, WLSMV" else "continuous",
      if (x$R > 0) paste0("; ", x$R, if (x$ci == "boot") " bootstrap resamples" else " Monte Carlo draws") else "", ")\n\n", sep = "")
  cat("General factor and wording shares of the total score:\n")
  g <- x$global
  keep <- intersect(c("group", "n", "omega_t", "omega_h", "omega_h_lo", "omega_h_hi",
                      "omega_w", "omega_h_blind", "fd_g", "fd_w", "proper"), names(g))
  print(r(g[, keep]), row.names = FALSE)
  if (any(!g$proper)) cat("\nWARNING: improper solution in group(s):", g$group[!g$proper], "\n")
  cat("\nFacet added value (purified target; var > 1 = added value;\n",
      "fd_s < .80 or neg_share > 0 suggest a poorly defined facet-specific factor):\n", sep = "")
  keep <- intersect(c("group", "facet", "prmse_sub", "prmse_tot", "prmse_aug", "var",
                      "var_lo", "var_hi", "p_added", "omega_hs", "omega_hs_blind", "fd_s",
                      "neg_share"),
                    names(x$facets))
  print(r(x$facets[, keep]), row.names = FALSE)
  if (!is.null(x$contrasts)) {
    cat("\nBetween-group contrasts:\n")
    print(r(x$contrasts), row.names = FALSE)
  }
  invisible(x)
}
