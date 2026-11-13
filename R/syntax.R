#' Build lavaan syntax for the five competing structures
#'
#' Generates model syntax for a one-factor model (\code{"one"}), correlated
#' facets (\code{"facets"}), a bifactor model (\code{"bifactor"}), a one-factor
#' model with a wording factor (\code{"one_method"}) and a bifactor model with a
#' wording factor (\code{"bifactor_method"}).
#'
#' @param facets Named list; each element is a character vector of item names
#'   belonging to one facet. Every item must belong to exactly one facet.
#' @param method Character vector of items that load on the wording (method)
#'   factor, typically the negatively keyed items after reverse scoring.
#' @param structure One of the five structure names above.
#' @return A character string with lavaan model syntax. Factor names are
#'   \code{G} (general), \code{S_<facet>} (facet-specific), \code{W} (wording)
#'   and, for the correlated-facet model, the facet names themselves.
#' @examples
#' f <- list(A = paste0("a", 1:3), B = paste0("b", 1:3))
#' cat(pure_syntax(f, method = c("a3", "b3"), structure = "bifactor_method"))
#' @export
pure_syntax <- function(facets, method = NULL,
                        structure = c("bifactor_method", "bifactor", "one_method",
                                      "facets", "one")) {
  structure <- match.arg(structure)
  check_facets(facets, method)
  items <- unlist(facets, use.names = FALSE)
  line <- function(f, v) paste0(f, " =~ ", paste(v, collapse = " + "))
  g <- line("G", items)
  s <- vapply(names(facets), function(k) line(paste0("S_", k), facets[[k]]), "")
  w <- if (length(method)) line("W", method) else character(0)
  if (structure %in% c("one_method", "bifactor_method") && !length(method))
    stop("`method` items are required for structure '", structure, "'.", call. = FALSE)
  out <- switch(structure,
    one = g,
    facets = vapply(names(facets), function(k) line(k, facets[[k]]), ""),
    bifactor = c(g, s),
    one_method = c(g, w),
    bifactor_method = c(g, s, w))
  paste(out, collapse = "\n")
}

check_facets <- function(facets, method) {
  if (!is.list(facets) || is.null(names(facets)) || any(names(facets) == ""))
    stop("`facets` must be a named list of item-name vectors.", call. = FALSE)
  if (length(facets) < 2) stop("At least two facets are required.", call. = FALSE)
  items <- unlist(facets, use.names = FALSE)
  if (anyDuplicated(items)) stop("Each item must belong to exactly one facet.", call. = FALSE)
  if (any(lengths(facets) < 2)) stop("Each facet needs at least two items.", call. = FALSE)
  if (length(method) && !all(method %in% items))
    stop("All `method` items must also appear in `facets`.", call. = FALSE)
  if (length(method) == 1) stop("The wording factor needs at least two items.", call. = FALSE)
  invisible(TRUE)
}

#' Fit and compare the competing structures
#'
#' Fits the five structures produced by \code{\link{pure_syntax}} with
#' \code{lavaan::cfa()}, using orthogonal factors for the bifactor and
#' method-augmented models and correlated factors for the facet model.
#'
#' @param data A data frame containing the items (already reverse-scored so
#'   that higher values mean more of the attribute).
#' @param facets,method See \code{\link{pure_syntax}}.
#' @param ordered Logical; if \code{TRUE} (default) items are treated as ordinal
#'   and estimated with WLSMV.
#' @param group Optional name of a grouping variable for multiple-group models.
#' @param group.equal Passed to lavaan; default \code{"loadings"}.
#' @param structures Which structures to fit.
#' @param ... Further arguments passed to \code{lavaan::cfa()}.
#' @return An object of class \code{pure_fits}: a list of fitted lavaan models
#'   plus a \code{table} element with scaled fit indices.
#' @export
pure_fit <- function(data, facets, method = NULL, ordered = TRUE, group = NULL,
                     group.equal = "loadings",
                     structures = c("one", "facets", "bifactor", "one_method",
                                    "bifactor_method"), ...) {
  check_facets(facets, method)
  items <- unlist(facets, use.names = FALSE)
  miss <- setdiff(c(items, group), names(data))
  if (length(miss)) stop("Variables not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!length(method)) structures <- setdiff(structures, c("one_method", "bifactor_method"))
  fits <- list()
  for (st in structures) {
    args <- list(model = pure_syntax(facets, method, st), data = data, std.lv = TRUE,
                 orthogonal = st != "facets", ...)
    if (ordered) { args$ordered <- items; args$estimator <- "WLSMV" }
    if (!is.null(group)) { args$group <- group; args$group.equal <- group.equal }
    fits[[st]] <- do.call("cfa", args, envir = asNamespace("lavaan"))
  }
  idx <- if (ordered) c("chisq.scaled", "df.scaled", "cfi.scaled", "tli.scaled",
                        "rmsea.scaled", "srmr") else
    c("chisq", "df", "cfi", "tli", "rmsea", "srmr")
  tab <- t(vapply(fits, function(f) {
    if (!lavaan::lavInspect(f, "converged")) return(rep(NA_real_, length(idx)))
    unname(lavaan::fitMeasures(f, idx))
  }, numeric(length(idx))))
  colnames(tab) <- c("chisq", "df", "CFI", "TLI", "RMSEA", "SRMR")
  out <- list(fits = fits, table = as.data.frame(tab), facets = facets, method = method,
              ordered = ordered, group = group)
  class(out) <- "pure_fits"
  out
}

#' @export
print.pure_fits <- function(x, digits = 3, ...) {
  cat("Competing structures (", if (x$ordered) "WLSMV, scaled indices" else "ML",
      ")\n", sep = "")
  tab <- x$table
  tab$chisq <- round(tab$chisq, 1)
  tab[, 3:6] <- round(tab[, 3:6], digits)
  print(tab)
  invisible(x)
}
