named_vc <- function(se, nms, rho = 0) {
  cor_mat <- matrix(rho, length(nms), length(nms))
  diag(cor_mat) <- 1
  m <- outer(se, se) * cor_mat
  dimnames(m) <- list(nms, nms)
  m
}
exp_spec <- function(rate, se = 0.1) {
  list(dist = "exp", coefs = c(rate = log(rate)), vcov = named_vc(se, "rate"))
}
wei_spec <- function(shape, scale, se = c(0.1, 0.1), rho = -0.3) {
  list(
    dist = "weibull",
    coefs = c(shape = log(shape), scale = log(scale)),
    vcov = named_vc(se, c("shape", "scale"), rho)
  )
}

th <- 20L
m_spec <- specify_m(tunnel_lengths = c(3L, 2L, th, th), pre_tunnel_states = 1L)
state_names <- c("on_trt", "surgery", "failed_surg", "responder", "palliative", "die")
tuota <- 30.44 / 365.25
cl <- 7 / 365.25
n_cycles <- c(on_trt = 20L, surgery = 3L, failed_surg = 2L, responder = 20L, palliative = 20L)

# surgery ends in a holding state powered by Weibull -> will warn unless a node
# is injected. failed_surg is exponential so its holding state is fine.
surv_all <- list(
  on_trt = list(surgery = exp_spec(0.05), die = exp_spec(0.01)),
  surgery = list(
    failed_surg = wei_spec(1.5, 3),
    responder = wei_spec(1.2, 5),
    die = exp_spec(0.05)
  ),
  failed_surg = list(palliative = exp_spec(0.27), die = exp_spec(0.02)),
  responder = list(palliative = wei_spec(1.1, 7), die = exp_spec(0.01)),
  palliative = list(die = exp_spec(0.25))
)

make_spec <- function(surv = surv_all) {
  specify_psa(surv, m_spec, state_names, th, tuota, cl)
}

# node at the end of surgery, as the vignette does it
make_node <- function(n, alpha = c(failed_surg = 35, responder = 55, die = 10), seed = 54321) {
  d <- dirichlet_draws(
    draw_uniforms(n, paste("surgery.node", names(alpha), sep = "."), seed = seed),
    alpha
  )
  draws <- lapply(setNames(names(alpha), names(alpha)), function(k) {
    rbind(matrix(NA_real_, 2, n), d[, k])
  })
  point <- lapply(setNames(names(alpha), names(alpha)), function(k) {
    c(NA, NA, unname(alpha[k] / sum(alpha)))
  })
  list(tp_draws = list(surgery = draws), tp_point = list(surgery = point))
}

# ---------------------------------------------------------------------------
# survival-only path
# ---------------------------------------------------------------------------

test_that("draw_psa draws natural-scale parameters of the right shape", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  n <- 50
  draws <- suppressWarnings(draw_psa(spec, n, seed = 1))

  expect_s3_class(draws, "treatseqr_psa_draws")
  expect_identical(draws$n_draws, 50L)
  expect_identical(colnames(draws$rng), spec$stream_names)
  expect_equal(dim(draws$rng), c(50L, 13L))

  p <- draws$params$surgery$failed_surg
  expect_identical(p$dist, "weibull")
  expect_equal(dim(p$params), c(2L, 50L))
  expect_identical(rownames(p$params), c("shape", "scale"))
  expect_equal(dim(draws$params$on_trt$die$params), c(1L, 50L))
  expect_identical(rownames(draws$params$on_trt$die$params), "rate")

  # natural scale: strictly positive
  all_params <- unlist(lapply(draws$params, function(s) lapply(s, `[[`, "params")))
  expect_true(all(all_params > 0))

  # nothing injected
  expect_true(all(!unlist(draws$injected_rows)))
  expect_null(draws$tp_draws)
})

test_that("draw_psa is reproducible and the rng path equals the seed path", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  a <- suppressWarnings(draw_psa(spec, 20, seed = 7))
  b <- suppressWarnings(draw_psa(spec, 20, seed = 7))
  expect_identical(a$params, b$params)

  rng <- draw_uniforms(20, spec$stream_names, seed = 7)
  c <- suppressWarnings(draw_psa(spec, 20, rng = rng))
  expect_identical(a$params, c$params)

  # extra streams in rng are fine; missing ones are not
  rng_extra <- cbind(rng, other.stream = runif(20))
  d <- suppressWarnings(draw_psa(spec, 20, rng = rng_extra))
  expect_identical(a$params, d$params)
  expect_error(
    suppressWarnings(draw_psa(spec, 20, rng = rng[, -1])),
    "Missing: 'on_trt.surgery.rate'"
  )
  expect_error(suppressWarnings(draw_psa(spec, 19, rng = rng)), "n_draws rows")
})

test_that("draw_psa draws are centred on the point estimates (log scale)", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  draws <- suppressWarnings(draw_psa(spec, 20000, seed = 3))
  p <- draws$params$surgery$failed_surg$params
  expect_true(abs(median(p["shape", ]) - 1.5) < 0.02)
  expect_true(abs(median(p["scale", ]) - 3) < 0.03)
  expect_true(abs(sd(log(p["shape", ])) - 0.1) < 0.005)
  expect_true(abs(cor(log(p["shape", ]), log(p["scale", ])) + 0.3) < 0.02)
})

test_that("seed/rng arguments are checked", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  expect_error(suppressWarnings(draw_psa(spec, 5)), "exactly one")
  expect_error(
    suppressWarnings(draw_psa(spec, 5, seed = 1, rng = matrix(0.5, 5, 13))),
    "exactly one"
  )
  expect_error(draw_psa(spec, 0, seed = 1), "n_draws")
  expect_error(draw_psa(list(), 5, seed = 1), "specify_psa")
})

test_that("tp_source_psa output feeds the deterministic pipeline", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  draws <- suppressWarnings(draw_psa(spec, 5, seed = 1))
  tp <- tp_source_psa(draws, iter = 3)

  expect_identical(names(tp), head(state_names, -1))
  expect_identical(names(tp$surgery), c("failed_surg", "responder", "die"))
  expect_identical(lengths(tp$on_trt), c(surgery = 20L, die = 20L))
  expect_identical(lengths(tp$surgery), c(failed_surg = 3L, responder = 3L, die = 3L))

  valid <- validate_tp_source(tp, tunnel_lengths = m_spec$tunnel_lengths)
  m <- generate_m_list(m_spec, tp, state_names = state_names)
  trace <- extrapolate_treatseqr(m, m_spec)
  expect_equal(unname(colSums(trace)), rep(1, th + 1))
})

test_that("tp_source_psa matches a hand calculation for one iteration", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  draws <- suppressWarnings(draw_psa(spec, 5, seed = 1))
  tp <- tp_source_psa(draws, iter = 2)

  times <- (0:20) * cl / tuota
  rate_s <- draws$params$on_trt$surgery$params["rate", 2]
  rate_d <- draws$params$on_trt$die$params["rate", 2]
  expected <- cum_haz_to_tps(list(
    surgery = flexsurv::Hexp(times, rate_s),
    die = flexsurv::Hexp(times, rate_d)
  ))
  expect_equal(tp$on_trt$surgery, unname(expected[, "surgery"]))
  expect_equal(tp$on_trt$die, unname(expected[, "die"]))
})

test_that("tp_source_psa(iter = NULL) gives the base case at the point estimates", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  draws <- suppressWarnings(draw_psa(spec, 5, seed = 1))
  bc <- tp_source_psa(draws)

  times <- (0:3) * cl / tuota
  expected <- cum_haz_to_tps(list(
    failed_surg = flexsurv::Hweibull(times, shape = 1.5, scale = 3),
    responder = flexsurv::Hweibull(times, shape = 1.2, scale = 5),
    die = flexsurv::Hexp(times, 0.05)
  ))
  expect_equal(bc$surgery$failed_surg, unname(expected[, "failed_surg"]))
  expect_equal(bc$surgery$responder, unname(expected[, "responder"]))
  expect_equal(bc$surgery$die, unname(expected[, "die"]))
  # the base case is not any particular iteration
  expect_false(identical(bc, tp_source_psa(draws, 1)))
})

test_that("tp_source_psa validates iter", {
  skip_if_not_installed("flexsurv")
  draws <- suppressWarnings(draw_psa(make_spec(), 5, seed = 1))
  expect_error(tp_source_psa(draws, 0), "1:5")
  expect_error(tp_source_psa(draws, 6), "1:5")
  expect_error(tp_source_psa(draws, 1.5), "1:5")
  expect_error(tp_source_psa(list(), 1), "draw_psa")
})

# ---------------------------------------------------------------------------
# holding-state warning
# ---------------------------------------------------------------------------

test_that("draw_psa warns when a fixed tunnel ends in a non-exponential holding state", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  expect_warning(
    draw_psa(spec, 5, seed = 1),
    "State 'surgery' is a fixed-duration tunnel \\(3 cycles\\).*failed_surg \\[weibull\\], responder \\[weibull\\]"
  )
  # failed_surg (exp only) and the full-length tunnels do not warn
  w <- testthat::capture_warnings(draw_psa(spec, 5, seed = 1))
  expect_length(w, 1L)
  expect_false(any(grepl("failed_surg' is a fixed", w)))
  expect_false(any(grepl("responder' is a fixed", w)))
})

test_that("no warning when the final cycle is injected or exponential", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  node <- make_node(5)
  expect_no_warning(suppressMessages(
    draw_psa(spec, 5, seed = 1, tp_draws = node$tp_draws, tp_point = node$tp_point)
  ))

  surv_exp <- surv_all
  surv_exp$surgery$failed_surg <- exp_spec(0.1)
  surv_exp$surgery$responder <- exp_spec(0.1)
  expect_no_warning(draw_psa(make_spec(surv_exp), 5, seed = 1))
})

# ---------------------------------------------------------------------------
# injection
# ---------------------------------------------------------------------------

test_that("injected node overrides the final row and zero-fills absent destinations", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  n <- 200
  node <- make_node(n)
  expect_message(
    draws <- draw_psa(spec, n, seed = 1, tp_draws = node$tp_draws, tp_point = node$tp_point),
    "state 'surgery' cycle\\(s\\) 3 injected"
  )
  expect_identical(draws$injected_rows$surgery, c(FALSE, FALSE, TRUE))
  expect_identical(draws$injected_rows$on_trt, rep(FALSE, 20))

  # every iteration: p_stay = 0 at cycle 3, survival rows unchanged at 1-2
  for (i in c(1, 77, n)) {
    tp <- tp_source_psa(draws, i)
    p_stay <- 1 - Reduce(`+`, tp$surgery)
    expect_equal(p_stay[3], 0)
    expect_true(all(p_stay[1:2] > 0.5))
    expect_equal(
      unname(c(tp$surgery$failed_surg[3], tp$surgery$responder[3], tp$surgery$die[3])),
      unname(c(
        node$tp_draws$surgery$failed_surg[3, i],
        node$tp_draws$surgery$responder[3, i],
        node$tp_draws$surgery$die[3, i]
      ))
    )
  }

  # base case uses tp_point
  bc <- tp_source_psa(draws)
  expect_equal(bc$surgery$failed_surg[3], 0.35)
  expect_equal(bc$surgery$responder[3], 0.55)
  expect_equal(bc$surgery$die[3], 0.10)
  # survival rows of the base case are at point estimates
  expect_equal(bc$surgery$die[1:2], tp_source_psa(suppressWarnings(draw_psa(spec, 2, seed = 9)))$surgery$die[1:2])
})

test_that("merge rule: NA destinations in an injected row become 0", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  n <- 10
  # inject only 'die' at surgery cycle 3 -> failed_surg/responder become 0 there
  tp_draws <- list(surgery = list(die = rbind(NA, NA, rep(0.4, n))))
  tp_point <- list(surgery = list(die = c(NA, NA, 0.4)))
  draws <- suppressMessages(draw_psa(spec, n, seed = 1, tp_draws = tp_draws, tp_point = tp_point))
  tp <- tp_source_psa(draws, 1)
  expect_equal(tp$surgery$die[3], 0.4)
  expect_equal(tp$surgery$failed_surg[3], 0)
  expect_equal(tp$surgery$responder[3], 0)
  expect_true(tp$surgery$failed_surg[1] > 0)
  # a destination only in the injection appears with zeros in survival rows
  tp_draws2 <- list(surgery = list(palliative = rbind(NA, NA, rep(0.4, n))))
  tp_point2 <- list(surgery = list(palliative = c(NA, NA, 0.4)))
  draws2 <- suppressMessages(draw_psa(spec, n, seed = 1, tp_draws = tp_draws2, tp_point = tp_point2))
  tp2 <- tp_source_psa(draws2, 1)
  expect_identical(names(tp2$surgery), c("failed_surg", "responder", "palliative", "die"))
  expect_equal(tp2$surgery$palliative, c(0, 0, 0.4))
})

test_that("a fully injected state needs no survival models; gaps are an error", {
  skip_if_not_installed("flexsurv")
  surv_no_fs <- surv_all[setdiff(names(surv_all), "failed_surg")]
  spec <- make_spec(surv_no_fs)
  n <- 4
  node <- make_node(n)
  fs_draws <- list(
    palliative = matrix(0.10, 2, n),
    die = matrix(c(0.15, 0.20), 2, n)
  )
  fs_point <- list(palliative = c(0.10, 0.10), die = c(0.15, 0.20))
  tp_draws <- c(node$tp_draws, list(failed_surg = fs_draws))
  tp_point <- c(node$tp_point, list(failed_surg = fs_point))

  draws <- suppressMessages(draw_psa(spec, n, seed = 1, tp_draws = tp_draws, tp_point = tp_point))
  expect_identical(draws$injected_rows$failed_surg, c(TRUE, TRUE))
  tp <- tp_source_psa(draws, 2)
  expect_equal(tp$failed_surg$die, c(0.15, 0.20))
  expect_equal(tp$failed_surg$palliative, c(0.10, 0.10))
  expect_no_error(generate_m_list(m_spec, tp, state_names = state_names))

  # leave cycle 2 of failed_surg uncovered
  fs_gap <- fs_draws
  fs_gap$palliative[2, ] <- NA
  fs_gap$die[2, ] <- NA
  tp_gap_pt <- fs_point
  tp_gap_pt$palliative[2] <- NA
  tp_gap_pt$die[2] <- NA
  expect_error(
    suppressMessages(draw_psa(
      spec, n, seed = 1,
      tp_draws = c(node$tp_draws, list(failed_surg = fs_gap)),
      tp_point = c(node$tp_point, list(failed_surg = tp_gap_pt))
    )),
    "State 'failed_surg' has no source of transition probabilities for cycle\\(s\\) 2"
  )
  # no survival models at all and no injection for a state
  expect_error(
    suppressMessages(draw_psa(spec, n, seed = 1)),
    "State 'failed_surg' has no source .* cycle\\(s\\) 1-2"
  )
})

test_that("surv = NULL with everything injected reproduces the deterministic vignette exactly", {
  n <- 3
  tp_source <- list(
    on_trt = list(
      surgery = seq(0.05, 0.02, length.out = th),
      palliative = rep(0.02, th),
      die = rep(0.01, th)
    ),
    surgery = list(
      failed_surg = c(0, 0, 0.35),
      responder = c(0, 0, 0.55),
      palliative = c(0.01, 0.01, 0),
      die = c(0.04, 0.04, 0.10)
    ),
    failed_surg = list(
      palliative = c(0.10, 0.10),
      die = c(0.15, 0.20)
    ),
    responder = list(
      palliative = seq(0.02, 0.05, length.out = th),
      die = seq(0.01, 0.05, length.out = th)
    ),
    palliative = list(die = rep(0.25, th))
  )
  tp_draws <- lapply(tp_source, function(s) lapply(s, function(v) matrix(v, length(v), n)))

  spec <- specify_psa(NULL, m_spec, state_names, th, tuota, cl)
  draws <- suppressMessages(draw_psa(spec, n, tp_draws = tp_draws, tp_point = tp_source))
  expect_identical(tp_source_psa(draws), tp_source)
  expect_identical(tp_source_psa(draws, 2), tp_source)

  # seed/rng are rejected when there is nothing to draw
  expect_error(
    suppressMessages(draw_psa(spec, n, seed = 1, tp_draws = tp_draws, tp_point = tp_source)),
    "no survival models"
  )

  # and the full deterministic pipeline agrees with the direct route
  m_direct <- generate_m_list(m_spec, tp_source, state_names = state_names)
  m_psa <- generate_m_list(m_spec, tp_source_psa(draws, 1), state_names = state_names)
  expect_identical(m_direct$m1, m_psa$m1)
  expect_identical(m_direct$m2, m_psa$m2)
})

test_that("injected inputs are validated", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  n <- 4
  node <- make_node(n)
  run <- function(tp_draws = node$tp_draws, tp_point = node$tp_point) {
    suppressMessages(draw_psa(spec, n, seed = 1, tp_draws = tp_draws, tp_point = tp_point))
  }

  expect_error(run(tp_point = NULL), "'tp_point'.*required")
  # structure
  expect_error(run(tp_draws = list(list(die = matrix(1, 3, n)))), "named list")
  expect_error(run(tp_draws = list(nope = node$tp_draws$surgery)), "'nope'")
  expect_error(run(tp_draws = list(die = node$tp_draws$surgery)), "'die'")
  bad <- node$tp_draws
  bad$surgery$surgery <- bad$surgery$die
  expect_error(run(tp_draws = bad), "'surgery'")
  # dims
  bad <- node$tp_draws
  bad$surgery$die <- bad$surgery$die[1:2, ]
  expect_error(run(tp_draws = bad), "3 x 4 matrix")
  bad <- node$tp_draws
  bad$surgery$die <- bad$surgery$die[, 1:3]
  expect_error(run(tp_draws = bad), "3 x 4 matrix")
  bad <- node$tp_draws
  bad$surgery$die <- as.vector(bad$surgery$die)
  expect_error(run(tp_draws = bad), "matrix")
  bad <- node$tp_point
  bad$surgery$die <- c(NA, 0.1)
  expect_error(run(tp_point = bad), "length 3")
  # values
  bad <- node$tp_draws
  bad$surgery$die[3, 1] <- 1.2
  expect_error(run(tp_draws = bad), "outside \\[0, 1\\]")
  bad <- node$tp_draws
  bad$surgery$die[3, 1] <- -0.1
  expect_error(run(tp_draws = bad), "outside \\[0, 1\\]")
  # row sums
  bad <- node$tp_draws
  bad$surgery$die[3, 2] <- 0.9
  expect_error(run(tp_draws = bad), "tp_draws\\$surgery: probabilities sum to > 1 in cycle\\(s\\) 3")
  bad <- node$tp_point
  bad$surgery$die[3] <- 0.9
  expect_error(run(tp_point = bad), "tp_point\\$surgery: probabilities sum to > 1")
  # NA pattern consistency
  bad <- node$tp_draws
  bad$surgery$die[3, 2] <- NA
  expect_error(run(tp_draws = bad), "mix NA and non-NA values across draws")
  bad <- node$tp_point
  bad$surgery$die[2] <- 0.1
  expect_error(run(tp_point = bad), "tp_point\\$surgery\\$die has values in cycle\\(s\\) 2")
  # states / destinations must agree between tp_draws and tp_point
  bad <- node$tp_point
  bad$surgery$palliative <- c(NA, NA, 0)
  expect_error(run(tp_point = bad), "same destinations")
  expect_error(
    run(tp_point = c(node$tp_point, list(palliative = list(die = rep(0.2, th))))),
    "same states"
  )
})

test_that("print methods summarise the draws", {
  skip_if_not_installed("flexsurv")
  spec <- make_spec()
  node <- make_node(5)
  draws <- suppressMessages(draw_psa(spec, 5, seed = 1, tp_draws = node$tp_draws, tp_point = node$tp_point))
  out <- capture.output(print(draws))
  expect_match(out[1], "treatseqr PSA draws")
  expect_true(any(grepl("5 iteration", out)))
  expect_true(any(grepl("surgery.*survival; injected at cycle\\(s\\) 3", out)))
  expect_true(any(grepl("on_trt.*survival$", out)))
})

test_that("format_cycles compacts runs", {
  expect_identical(treatseqr:::format_cycles(integer(0)), "")
  expect_identical(treatseqr:::format_cycles(3L), "3")
  expect_identical(treatseqr:::format_cycles(c(1L, 2L, 3L, 5L, 7L, 8L)), "1-3, 5, 7-8")
})
