#' Simulate ordinal responses from a bifactor model with a wording factor
#'
#' Generates Likert-type responses from a model with a general factor, one
#' specific factor per facet and a wording factor that loads on the negatively
#' keyed items (after reverse scoring). All factors are orthogonal and
#' standard normal; responses are obtained by cutting latent responses at the
#' supplied thresholds.
#'
#' @param n Sample size.
#' @param n_facets Number of facets.
#' @param items_per_facet Items per facet.
#' @param n_neg Number of negatively keyed (method) items per facet.
#' @param lambda_g General-factor loading (scalar or vector of length items).
#' @param lambda_s Facet-specific loading (scalar or vector).
#' @param lambda_w Wording loading for method items (scalar or vector).
#' @param thresholds Thresholds on the latent-response scale; the default
#'   gives five categories with a mild ceiling.
#' @return A list with \code{data} (data frame of integer responses),
#'   \code{facets}, \code{method} and the population \code{loadings}.
#' @examples
#' sim <- simulate_pure(n = 200)
#' str(sim$data[, 1:4])
#' @export
simulate_pure <- function(n = 1000, n_facets = 4, items_per_facet = 6, n_neg = 3,
                          lambda_g = 0.6, lambda_s = 0.3, lambda_w = 0.4,
                          thresholds = c(-2, -1.2, -0.3, 0.6)) {
  if (n_neg >= items_per_facet) stop("`n_neg` must be smaller than `items_per_facet`.", call. = FALSE)
  J <- n_facets * items_per_facet
  rep_to <- function(x) if (length(x) == 1) rep(x, J) else x
  lG <- rep_to(lambda_g); lS <- rep_to(lambda_s)
  fac <- rep(seq_len(n_facets), each = items_per_facet)
  isneg <- rep(c(rep(FALSE, items_per_facet - n_neg), rep(TRUE, n_neg)), n_facets)
  lW <- ifelse(isneg, rep_to(lambda_w), 0)
  com <- lG^2 + lS^2 + lW^2
  if (any(com >= 1)) stop("Loadings imply communalities >= 1.", call. = FALSE)
  G <- stats::rnorm(n); S <- matrix(stats::rnorm(n * n_facets), n); W <- stats::rnorm(n)
  E <- matrix(stats::rnorm(n * J), n) %*% diag(sqrt(1 - com))
  Y <- outer(G, lG) + S[, fac] * rep(lS, each = n) + outer(W, lW) + E
  X <- matrix(findInterval(Y, sort(thresholds)) + 1L, n, J)
  nm <- paste0(LETTERS[fac], stats::ave(fac, fac, FUN = seq_along))
  colnames(X) <- nm
  facets <- split(nm, LETTERS[fac])
  list(data = as.data.frame(X), facets = facets, method = nm[isneg],
       loadings = data.frame(item = nm, facet = LETTERS[fac], general = lG,
                             specific = lS, wording = lW))
}

#' Population values of the purified indices
#'
#' Computes the population values of the indices returned by
#' \code{\link{pure_value}} for a given set of standardized loadings, for
#' example those stored in \code{simulate_pure()$loadings}.
#'
#' @param loadings Data frame with columns \code{facet}, \code{general},
#'   \code{specific} and \code{wording}.
#' @return A list with \code{global} and \code{facets} indices.
#' @export
pure_population <- function(loadings) {
  f <- loadings$facet
  K <- unique(f)
  lS <- vapply(K, function(k) ifelse(f == k, loadings$specific, 0), numeric(length(f)))
  lS <- matrix(lS, ncol = length(K))
  th <- 1 - loadings$general^2 - rowSums(lS^2) - loadings$wording^2
  fidx <- lapply(K, function(k) which(f == k)); names(fidx) <- K
  pure_indices(loadings$general, lS, loadings$wording, th, fidx)
}
