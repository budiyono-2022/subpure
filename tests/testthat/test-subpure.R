test_that("syntax builder covers all structures", {
  f <- list(A = paste0("a", 1:3), B = paste0("b", 1:3))
  s <- pure_syntax(f, method = c("a3", "b3"), structure = "bifactor_method")
  expect_match(s, "G =~ a1 \\+ a2")
  expect_match(s, "S_A =~ a1 \\+ a2 \\+ a3")
  expect_match(s, "W =~ a3 \\+ b3")
  expect_match(pure_syntax(f, structure = "facets"), "^A =~")
  expect_error(pure_syntax(f, structure = "one_method"), "method")
})

test_that("input checks catch malformed facets", {
  expect_error(pure_syntax(list(c("a", "b"), c("c", "d"))), "named list")
  expect_error(pure_syntax(list(A = c("a", "b"), B = c("b", "c"))), "exactly one")
  expect_error(pure_syntax(list(A = c("a", "b"), B = c("c", "d")), method = "z",
                           structure = "bifactor_method"), "also appear|at least two")
})

test_that("population indices behave as expected", {
  sim <- simulate_pure(n = 10, lambda_s = 0, lambda_w = 0)
  p <- pure_population(sim$loadings)
  expect_equal(unname(p$global["omega_w"]), 0)
  expect_true(all(p$facets[, "omega_hs"] == 0))
  expect_true(all(p$facets[, "var"] < 1))
  sim2 <- simulate_pure(n = 10, lambda_s = 0.5, lambda_w = 0)
  expect_true(all(pure_population(sim2$loadings)$facets[, "var"] > 1))
  # wording inflates nothing in the purified target but lowers omega_h
  sim3 <- simulate_pure(n = 10, lambda_s = 0, lambda_w = 0.5)
  expect_lt(pure_population(sim3$loadings)$global[["omega_h"]], p$global[["omega_h"]])
})

test_that("augmented PRMSE is never below subscore or total PRMSE", {
  sim <- simulate_pure(n = 10, lambda_s = 0.3, lambda_w = 0.4)
  f <- pure_population(sim$loadings)$facets
  expect_true(all(f[, "prmse_aug"] >= pmax(f[, "prmse_sub"], f[, "prmse_tot"]) - 1e-12))
})

test_that("simulated data have the documented shape", {
  set.seed(3)
  sim <- simulate_pure(n = 50, n_facets = 3, items_per_facet = 4, n_neg = 2)
  expect_equal(dim(sim$data), c(50, 12))
  expect_length(sim$method, 6)
  expect_true(all(unlist(sim$data) %in% 1:5))
})

test_that("pure_value recovers population values", {
  skip_on_cran()
  set.seed(11)
  sim <- simulate_pure(n = 2000, lambda_s = 0.35, lambda_w = 0.4)
  pv <- pure_value(sim$data, sim$facets, sim$method, R = 100)
  pop <- pure_population(sim$loadings)
  expect_equal(pv$global$omega_h, pop$global[["omega_h"]], tolerance = 0.05)
  expect_equal(mean(pv$facets$var), mean(pop$facets[, "var"]), tolerance = 0.05)
  expect_gt(pv$global$omega_h_blind, pv$global$omega_h)
  expect_true(all(c("var_lo", "var_hi", "p_added", "fd_s", "neg_share") %in% names(pv$facets)))
  expect_true(pv$global$proper)
  expect_identical(pv$global$intervals, "mc")
  expect_gt(pv$global$fd_g, 0.9)
})

test_that("group contrasts are returned", {
  skip_on_cran()
  set.seed(5)
  sim <- simulate_pure(n = 800, lambda_s = 0.3, lambda_w = 0.4)
  d <- sim$data
  d$grp <- rep(c("a", "b"), each = 400)
  pv <- pure_value(d, sim$facets, sim$method, group = "grp", R = 50)
  expect_equal(sort(unique(pv$facets$group)), c("a", "b"))
  expect_true(any(pv$contrasts$index == "omega_h"))
})

test_that("bundled data are reverse-scored and complete", {
  expect_equal(nrow(rmbangka), 4197)
  items <- names(rmbangka)[7:42]
  r <- stats::cor(rmbangka[items])
  expect_true(all(r[upper.tri(r)] > -0.2))
  expect_false(anyNA(rmbangka[items]))
})

test_that("intervals can be skipped and determinacy is bounded", {
  skip_on_cran()
  set.seed(8)
  sim <- simulate_pure(n = 1000, lambda_s = 0.4, lambda_w = 0.4)
  pv <- pure_value(sim$data, sim$facets, sim$method, R = 0)
  expect_identical(pv$global$intervals, "none")
  expect_true(all(pv$facets$fd_s > 0 & pv$facets$fd_s <= 1))
  expect_null(pv$contrasts)
})

test_that("invariance sequence runs and reports nested comparisons", {
  skip_on_cran()
  skip_if_not_installed("semTools")
  set.seed(3)
  sim <- simulate_pure(n = 1000, lambda_s = 0.3, lambda_w = 0.4)
  d <- sim$data
  d$grp <- rep(c("a", "b"), each = 500)
  inv <- pure_invariance(d, sim$facets, sim$method, group = "grp")
  expect_s3_class(inv, "pure_invariance")
  expect_equal(inv$table$step, c("configural", "thresholds", "loadings"))
  expect_true(all(diff(inv$table$df) > 0))
  expect_true(all(!is.na(inv$table$p[-1])))
  expect_true(all(inv$table$supported[-1]))
  expect_error(pure_invariance(d, sim$facets, sim$method, group = "nope"), "group")
})
