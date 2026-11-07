# draw_uniforms ---------------------------------------------------------------

test_that("draw_uniforms returns a named matrix of the right shape in (0, 1)", {
  rng <- draw_uniforms(7, c("a.x", "a.y", "b.z"), seed = 1)
  expect_true(is.matrix(rng))
  expect_equal(dim(rng), c(7L, 3L))
  expect_identical(colnames(rng), c("a.x", "a.y", "b.z"))
  expect_true(all(rng > 0 & rng < 1))
})

test_that("draw_uniforms is reproducible for a given seed", {
  a <- draw_uniforms(10, c("p", "q"), seed = 99)
  b <- draw_uniforms(10, c("p", "q"), seed = 99)
  c <- draw_uniforms(10, c("p", "q"), seed = 100)
  expect_identical(a, b)
  expect_false(identical(a, c))
})

test_that("draw_uniforms does not disturb the global RNG state", {
  set.seed(1)
  before <- runif(3)
  set.seed(1)
  invisible(draw_uniforms(50, c("a", "b"), seed = 12345))
  after <- runif(3)
  expect_identical(before, after)
})

test_that("draw_uniforms restores a missing .Random.seed cleanly", {
  # Make sure .Random.seed exists first, remove it, then check the function
  # leaves it absent (rather than creating one as a side effect).
  set.seed(1)
  rm(".Random.seed", envir = globalenv())
  invisible(draw_uniforms(5, "a", seed = 1))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("draw_uniforms validates its inputs", {
  expect_error(draw_uniforms(0, "a", seed = 1), "n_draws")
  expect_error(draw_uniforms(5, c("a", "a"), seed = 1), "unique")
  expect_error(draw_uniforms(5, character(0), seed = 1), "stream_names")
  expect_error(draw_uniforms(5, "a", seed = "x"), "seed")
})

# mvnorm_draws ----------------------------------------------------------------

named_vc <- function(se, nms, rho = 0) {
  cor_mat <- matrix(rho, length(nms), length(nms))
  diag(cor_mat) <- 1
  m <- outer(se, se) * cor_mat
  dimnames(m) <- list(nms, nms)
  m
}

test_that("mvnorm_draws reproduces means and covariance", {
  nms <- c("shape", "scale")
  mu <- c(shape = log(1.5), scale = log(3))
  vc <- named_vc(c(0.1, 0.2), nms, rho = -0.4)
  rng <- draw_uniforms(20000, nms, seed = 3)
  d <- mvnorm_draws(rng, mu, vcov = vc)

  expect_equal(dim(d), c(20000L, 2L))
  expect_identical(colnames(d), nms)
  expect_true(max(abs(colMeans(d) - mu)) < 0.01)
  expect_true(max(abs(stats::cov(d) - vc)) < 0.005)
  expect_true(abs(stats::cor(d)[1, 2] - (-0.4)) < 0.02)
})

test_that("mvnorm_draws accepts a pre-computed Cholesky factor and matches vcov path", {
  nms <- c("a", "b", "c")
  mu <- c(a = 1, b = 2, c = 3)
  vc <- named_vc(c(1, 1, 1), nms, rho = 0.3)
  R <- chol(vc, pivot = TRUE)
  R <- R[, order(attr(R, "pivot"))]
  rng <- draw_uniforms(100, nms, seed = 4)
  expect_identical(
    mvnorm_draws(rng, mu, vcov = vc),
    mvnorm_draws(rng, mu, R = R)
  )
})

test_that("mvnorm_draws dimname checks actually fire", {
  mu <- c(rate = log(0.1))
  rng <- draw_uniforms(5, "rate", seed = 1)
  # no dimnames at all
  expect_error(mvnorm_draws(rng, mu, vcov = matrix(0.01)), "names")
  # wrong names
  expect_error(
    mvnorm_draws(rng, mu, vcov = named_vc(0.1, "shape")),
    "names"
  )
  # unnamed params
  expect_error(mvnorm_draws(rng, 0.1, vcov = named_vc(0.1, "rate")), "named")
})

test_that("mvnorm_draws requires exactly one of vcov or R", {
  mu <- c(rate = 0)
  rng <- draw_uniforms(5, "rate", seed = 1)
  vc <- named_vc(0.1, "rate")
  expect_error(mvnorm_draws(rng, mu), "exactly one")
  expect_error(mvnorm_draws(rng, mu, vcov = vc, R = chol(vc)), "exactly one")
  expect_error(mvnorm_draws(rng[, 0, drop = FALSE], mu, vcov = vc), "rng")
})

# dirichlet_draws -------------------------------------------------------------

test_that("dirichlet_draws rows sum to 1 and are centred on alpha / sum(alpha)", {
  alpha <- c(failed_surg = 35, responder = 55, die = 10)
  rng <- draw_uniforms(20000, names(alpha), seed = 5)
  d <- dirichlet_draws(rng, alpha)
  expect_equal(dim(d), c(20000L, 3L))
  expect_identical(colnames(d), names(alpha))
  expect_true(all(abs(rowSums(d) - 1) < 1e-12))
  expect_true(max(abs(colMeans(d) - alpha / sum(alpha))) < 0.01)
  expect_true(all(d >= 0 & d <= 1))
})

test_that("dirichlet_draws validates its inputs", {
  rng <- draw_uniforms(5, c("a", "b"), seed = 1)
  expect_error(dirichlet_draws(rng, c(1, 0)), "positive")
  expect_error(dirichlet_draws(rng, c(1, 1, 1)), "one column per")
})

# cum_haz_to_tps --------------------------------------------------------------

test_that("cum_haz_to_tps matches the closed form for a single exponential", {
  rate <- 0.1
  H <- list(die = rate * 0:5)
  p <- cum_haz_to_tps(H)
  expect_equal(dim(p), c(5L, 1L))
  expect_identical(colnames(p), "die")
  expect_equal(unname(p[, "die"]), rep(1 - exp(-rate), 5))
})

test_that("cum_haz_to_tps shares exits in proportion to hazard increments", {
  H <- list(a = 0.1 * 0:3, b = 0.3 * 0:3)
  p <- cum_haz_to_tps(H)
  total <- 1 - exp(-0.4)
  expect_equal(unname(rowSums(p)), rep(total, 3))
  expect_equal(unname(p[, "a"] / p[, "b"]), rep(1 / 3, 3))
})

test_that("cum_haz_to_tps never exceeds 1 even for very large hazards", {
  H <- list(a = 50 * 0:4, b = 50 * 0:4)
  p <- cum_haz_to_tps(H)
  expect_true(all(rowSums(p) <= 1))
  expect_true(all(p >= 0))
})

test_that("cum_haz_to_tps returns zeros for cycles with no hazard", {
  H <- list(a = c(0, 0, 0.5, 0.5), b = c(0, 0, 0, 0.2))
  p <- cum_haz_to_tps(H)
  expect_equal(unname(p[1, ]), c(0, 0))
  expect_true(p[2, "a"] > 0)
  expect_equal(unname(p[2, "b"]), 0)
})

test_that("cum_haz_to_tps validates its inputs", {
  expect_error(cum_haz_to_tps(list(0:3)), "named")
  expect_error(cum_haz_to_tps(list(a = 0:3, b = 0:2)), "same length")
  expect_error(cum_haz_to_tps(list(a = c(0, 1, 0.5))), "non-decreasing")
  expect_error(cum_haz_to_tps(list(a = c(0, Inf))), "finite")
})

# flexsurv interface ----------------------------------------------------------

test_that("flexsurv_dist_info is case-insensitive and accepts 'exponential'", {
  skip_if_not_installed("flexsurv")
  expect_identical(treatseqr:::flexsurv_dist_info("WeibullPH")$name, "weibullPH")
  expect_identical(treatseqr:::flexsurv_dist_info("exponential")$name, "exp")
  expect_identical(treatseqr:::flexsurv_dist_info("weibull")$pars, c("shape", "scale"))
  expect_error(treatseqr:::flexsurv_dist_info("notadist"), "Unsupported distribution")
})

test_that("flexsurv_inv_transform reproduces a fitted model's natural-scale estimates", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  fit <- flexsurv::flexsurvreg(
    survival::Surv(recyrs, censrec) ~ 1,
    data = flexsurv::bc,
    dist = "weibull"
  )
  est <- flexsurv_inv_transform(t(coef(fit)), "weibull")
  expect_equal(as.numeric(est), unname(fit$res[, "est"]))
  expect_identical(colnames(est), c("shape", "scale"))
})

test_that("flexsurv_inv_transform applies the right transform per parameter", {
  skip_if_not_installed("flexsurv")
  # gompertz: shape is identity, rate is log
  d <- matrix(c(0.1, log(0.2)), 1, dimnames = list(NULL, c("shape", "rate")))
  out <- flexsurv_inv_transform(d, "gompertz")
  expect_equal(unname(out[1, ]), c(0.1, 0.2))
  # wrong column order / names is an error
  expect_error(
    flexsurv_inv_transform(d[, 2:1, drop = FALSE], "gompertz"),
    "colnames"
  )
})

test_that("flexsurv_cum_haz routes to the right function", {
  skip_if_not_installed("flexsurv")
  H <- flexsurv_cum_haz("weibull")
  expect_identical(H, flexsurv::Hweibull)
  expect_identical(flexsurv_cum_haz("Exponential"), flexsurv::Hexp)
  expect_equal(flexsurv_cum_haz("exp")(0:3, rate = 0.1), 0.1 * 0:3)
})

test_that("cum_haz_to_tps + flexsurv_cum_haz agree with the exponential closed form", {
  skip_if_not_installed("flexsurv")
  tuota <- 30.44 / 365.25
  cl <- 7 / 365.25
  times <- (0:5) * cl / tuota
  p <- cum_haz_to_tps(list(die = flexsurv_cum_haz("exp")(times, rate = 0.1)))
  expect_equal(unname(p[, "die"]), rep(1 - exp(-0.1 * cl / tuota), 5))
})
