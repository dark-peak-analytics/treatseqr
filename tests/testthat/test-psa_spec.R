named_vc <- function(se, nms, rho = 0) {
  cor_mat <- matrix(rho, length(nms), length(nms))
  diag(cor_mat) <- 1
  m <- outer(se, se) * cor_mat
  dimnames(m) <- list(nms, nms)
  m
}

th <- 20L
m_spec <- specify_m(tunnel_lengths = c(3L, 2L, th, th), pre_tunnel_states = 1L)
state_names <- c("on_trt", "surgery", "failed_surg", "responder", "palliative", "die")
tuota <- 30.44 / 365.25
cl <- 7 / 365.25

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

surv_ok <- list(
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

test_that("specify_psa builds a spec from dist/coefs/vcov lists", {
  skip_if_not_installed("flexsurv")
  spec <- specify_psa(surv_ok, m_spec, state_names, th, tuota, cl)

  expect_s3_class(spec, "treatseqr_psa_spec")
  expect_identical(spec$state_names, state_names)
  expect_identical(spec$death_name, "die")
  expect_identical(
    spec$n_cycles,
    c(on_trt = 20L, surgery = 3L, failed_surg = 2L, responder = 20L, palliative = 20L)
  )
  expect_identical(spec$th, 20L)
  expect_equal(spec$tuota, tuota)
  expect_equal(spec$cl, cl)

  # normalised entries keep dist/coefs/vcov
  expect_identical(spec$surv$surgery$failed_surg$dist, "weibull")
  expect_identical(names(spec$surv$surgery$failed_surg$coefs), c("shape", "scale"))
  expect_identical(spec$surv$on_trt$die, exp_spec(0.01))
})

test_that("stream names are state.destination.parameter, unique and in order", {
  skip_if_not_installed("flexsurv")
  spec <- specify_psa(surv_ok, m_spec, state_names, th, tuota, cl)
  expect_identical(
    spec$stream_names[1:5],
    c(
      "on_trt.surgery.rate", "on_trt.die.rate",
      "surgery.failed_surg.shape", "surgery.failed_surg.scale",
      "surgery.responder.shape"
    )
  )
  expect_false(anyDuplicated(spec$stream_names) > 0)
  # 2 + (2 + 2 + 1) + (1 + 1) + (2 + 1) + 1 = 13
  expect_length(spec$stream_names, 13L)
})

test_that("specify_psa accepts flexsurvreg fits and canonicalises them", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  fit_exp <- flexsurv::flexsurvreg(
    survival::Surv(recyrs, censrec) ~ 1,
    data = flexsurv::bc, dist = "exp"
  )
  fit_wei <- flexsurv::flexsurvreg(
    survival::Surv(recyrs, censrec) ~ 1,
    data = flexsurv::bc, dist = "weibull"
  )
  spec <- specify_psa(
    list(on_trt = list(surgery = fit_wei, die = fit_exp)),
    m_spec, state_names, th, tuota, cl
  )
  # "weibull.quiet" -> "weibull"
  expect_identical(spec$surv$on_trt$surgery$dist, "weibull")
  expect_identical(
    unname(spec$surv$on_trt$surgery$coefs),
    unname(fit_wei$res.t[, "est"])
  )
  expect_identical(names(spec$surv$on_trt$surgery$coefs), c("shape", "scale"))
  expect_equal(unname(spec$surv$on_trt$surgery$vcov), unname(vcov(fit_wei)))
  expect_identical(dimnames(spec$surv$on_trt$surgery$vcov), list(c("shape", "scale"), c("shape", "scale")))

  # single-parameter exponential: coef()/vcov() drop names, we restore them
  expect_identical(spec$surv$on_trt$die$dist, "exp")
  expect_identical(names(spec$surv$on_trt$die$coefs), "rate")
  expect_identical(dimnames(spec$surv$on_trt$die$vcov), list("rate", "rate"))
  expect_equal(unname(spec$surv$on_trt$die$coefs), unname(fit_exp$res.t[, "est"]))
})

test_that("specify_psa rejects flexsurvreg fits with covariates", {
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")
  fit_cov <- flexsurv::flexsurvreg(
    survival::Surv(recyrs, censrec) ~ group,
    data = flexsurv::bc, dist = "weibull"
  )
  expect_error(
    specify_psa(list(on_trt = list(die = fit_cov)), m_spec, state_names, th, tuota, cl),
    "covariates"
  )
})

test_that("specify_psa allows surv = NULL (everything injected)", {
  spec <- specify_psa(NULL, m_spec, state_names, th, tuota, cl)
  expect_length(spec$surv, 0L)
  expect_identical(spec$stream_names, character(0))
  expect_identical(unname(spec$n_cycles), c(20L, 3L, 2L, 20L, 20L))
})

test_that("specify_psa validates the model structure", {
  skip_if_not_installed("flexsurv")
  expect_error(
    specify_psa(surv_ok, list(), state_names, th, tuota, cl),
    "specify_m"
  )
  expect_error(
    specify_psa(surv_ok, m_spec, state_names[-1], th, tuota, cl),
    "length\\(state_names\\)"
  )
  expect_error(
    specify_psa(surv_ok, m_spec, state_names, 1L, tuota, cl),
    "'th'"
  )
  expect_error(
    specify_psa(surv_ok, m_spec, state_names, 20.5, tuota, cl),
    "'th'"
  )
  expect_error(
    specify_psa(surv_ok, m_spec, state_names, th, -1, cl),
    "tuota"
  )
})

test_that("specify_psa validates state and destination names", {
  skip_if_not_installed("flexsurv")
  bad_state <- surv_ok
  bad_state$not_a_state <- list(die = exp_spec(0.1))
  expect_error(
    specify_psa(bad_state, m_spec, state_names, th, tuota, cl),
    "not_a_state"
  )

  # death state cannot be a source
  bad_death <- surv_ok
  bad_death$die <- list(on_trt = exp_spec(0.1))
  expect_error(
    specify_psa(bad_death, m_spec, state_names, th, tuota, cl),
    "'die'"
  )

  bad_dest <- surv_ok
  bad_dest$on_trt$nowhere <- exp_spec(0.1)
  expect_error(
    specify_psa(bad_dest, m_spec, state_names, th, tuota, cl),
    "nowhere"
  )

  # self-transition is not a destination
  self_dest <- surv_ok
  self_dest$on_trt$on_trt <- exp_spec(0.1)
  expect_error(
    specify_psa(self_dest, m_spec, state_names, th, tuota, cl),
    "surv\\$on_trt.*'on_trt'"
  )

  # unnamed destination list
  unnamed <- surv_ok
  unnamed$palliative <- list(exp_spec(0.1))
  expect_error(
    specify_psa(unnamed, m_spec, state_names, th, tuota, cl),
    "named list of destinations"
  )
})

test_that("specify_psa validates each survival specification", {
  skip_if_not_installed("flexsurv")
  with_dest <- function(x) {
    s <- surv_ok
    s$palliative$die <- x
    s
  }
  run <- function(x) specify_psa(with_dest(x), m_spec, state_names, th, tuota, cl)

  expect_error(run(list(dist = "exp", coefs = c(rate = 0))), "'vcov'")
  expect_error(run(list(dist = "exp", coefs = c(rate = 0), vcov = named_vc(0.1, "rate"), extra = 1)), "unexpected")
  expect_error(run(list(dist = "notadist", coefs = c(rate = 0), vcov = named_vc(0.1, "rate"))), "Unsupported distribution")
  # wrong parameter names for the distribution
  expect_error(run(list(dist = "weibull", coefs = c(rate = 0), vcov = named_vc(0.1, "rate"))), "'shape', 'scale'")
  # wrong order
  expect_error(
    run(list(dist = "weibull", coefs = c(scale = 1, shape = 0), vcov = named_vc(c(0.1, 0.1), c("scale", "shape")))),
    "'shape', 'scale'"
  )
  # unnamed coefs
  expect_error(run(list(dist = "exp", coefs = 0, vcov = named_vc(0.1, "rate"))), "named exactly")
  # vcov without dimnames
  expect_error(run(list(dist = "exp", coefs = c(rate = 0), vcov = matrix(0.01))), "vcov")
  # vcov wrong size
  expect_error(run(list(dist = "exp", coefs = c(rate = 0), vcov = named_vc(c(0.1, 0.1), c("rate", "x")))), "vcov")
  # non-symmetric vcov
  vc <- named_vc(c(0.1, 0.1), c("shape", "scale"))
  vc[1, 2] <- 0.5
  expect_error(run(list(dist = "weibull", coefs = c(shape = 0, scale = 0), vcov = vc)), "symmetric")
  # non-finite
  expect_error(run(list(dist = "exp", coefs = c(rate = NA_real_), vcov = named_vc(0.1, "rate"))), "finite")
})

test_that("specify_psa error messages name the offending transition", {
  skip_if_not_installed("flexsurv")
  s <- surv_ok
  s$responder$palliative$vcov <- matrix(0.01, 2, 2)
  expect_error(
    specify_psa(s, m_spec, state_names, th, tuota, cl),
    "surv\\$responder\\$palliative"
  )
})

test_that("print method summarises the spec", {
  skip_if_not_installed("flexsurv")
  spec <- specify_psa(surv_ok[c("on_trt", "surgery")], m_spec, state_names, th, tuota, cl)
  out <- capture.output(print(spec))
  expect_match(out[1], "treatseqr PSA specification")
  expect_true(any(grepl("surgery.*failed_surg \\[weibull\\]", out)))
  expect_true(any(grepl("failed_surg.*inject via tp_draws", out)))
  # on_trt: 1 + 1; surgery: 2 + 2 + 1
  expect_true(any(grepl("7 survival parameter stream", out)))
})
