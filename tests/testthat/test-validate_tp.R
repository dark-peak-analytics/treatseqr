test_that("validate_tp_source works with valid input (vignette example)", {
  th_cyc <- 10

  tp_source <- list(
    on_1l = list(
      off_1l = seq(0.1, 0.05, length.out = th_cyc),
      on_2l = seq(0.05, 0.1, length.out = th_cyc),
      off_2l = rep(0, th_cyc),
      bsc = NULL,
      die = rep(0.01, th_cyc)
    ),
    off_1l = list(
      on_2l = seq(0.08, 0.02, length.out = th_cyc),
      off_2l = seq(0.01, 0.03, length.out = th_cyc),
      bsc = rep(0.02, th_cyc),
      die = seq(0.05, 0.1, length.out = th_cyc)
    ),
    on_2l = list(
      off_2l = seq(0.05, 0.02, length.out = th_cyc),
      bsc = rep(0.01, th_cyc),
      die = seq(0.1, 0.15, length.out = th_cyc)
    ),
    off_2l = list(
      bsc = seq(0.075, 0.05, length.out = th_cyc),
      die = seq(0.15, 0.2, length.out = th_cyc)
    ),
    bsc = list(
      die = seq(0.2, 0.3, length.out = th_cyc)
    )
  )

  res <- validate_tp_source(tp_source)

  # Check NULL replacement
  expect_equal(res$on_1l$bsc, rep(0, th_cyc))

  # Check structure integrity: all inner lists now uniform (canonical set minus self)
  expect_equal(names(res), names(tp_source))
  expect_equal(length(res), 5)
  expect_true(all(vapply(res, length, numeric(1)) == 5L))

  # Backwards destinations added for states that previously had fewer destinations
  expect_equal(names(res$off_1l), c("on_1l", "on_2l", "off_2l", "bsc", "die"))
  expect_equal(res$off_1l$on_1l, rep(0, th_cyc))

  expect_equal(names(res$on_2l), c("on_1l", "off_1l", "off_2l", "bsc", "die"))
  expect_equal(res$on_2l$on_1l, rep(0, th_cyc))
  expect_equal(res$on_2l$off_1l, rep(0, th_cyc))
})

test_that("validate_tp_source auto-expands reducing structure to uniform", {
  th_cyc <- 5
  tp_reducing <- list(
    a = list(
      b = rep(0.1, th_cyc),
      c = rep(0.05, th_cyc),
      d = rep(0.02, th_cyc)
    ),
    b = list(c = rep(0.10, th_cyc), d = rep(0.05, th_cyc)),
    c = list(d = rep(0.10, th_cyc))
  )
  res <- validate_tp_source(tp_reducing)

  # All inner lists now uniform length = 3 (canonical 4 names minus self)
  expect_true(all(vapply(res, length, numeric(1)) == 3L))

  # b and c gain backwards destinations as zero vectors
  expect_equal(res$b$a, rep(0, th_cyc))
  expect_equal(res$c$a, rep(0, th_cyc))
  expect_equal(res$c$b, rep(0, th_cyc))

  # Names follow canonical order
  expect_equal(names(res$b), c("a", "c", "d"))
  expect_equal(names(res$c), c("a", "b", "d"))
})

test_that("validate_tp_source errors when tp_source names missing from state_names", {
  tp <- list(
    a = list(b = 0.1, c = 0.05),
    b = list(c = 0.10)
  )
  expect_error(
    validate_tp_source(tp, state_names = c("x", "b", "c")),
    "'state_names' must match names\\(tp_source\\) in the same order"
  )
})

test_that("validate_tp_source errors on reordered state_names", {
  th_cyc <- 2
  tp <- list(
    pre = list(
      tun1 = rep(0.11, th_cyc),
      tun2 = rep(0.22, th_cyc),
      die = rep(0.03, th_cyc)
    ),
    tun1 = list(tun2 = rep(0.44, th_cyc), die = rep(0.04, th_cyc)),
    tun2 = list(tun1 = rep(0.55, th_cyc), die = rep(0.05, th_cyc))
  )
  # Same set of names, but tunnels swapped relative to names(tp): the matrix
  # layout is positional, so this must be rejected rather than silently
  # misplacing probabilities
  expect_error(
    validate_tp_source(tp, state_names = c("pre", "tun2", "tun1", "die")),
    "'state_names' must match names\\(tp_source\\) in the same order"
  )
})

test_that("validate_tp_source errors when death name in state_names mismatches", {
  tp <- list(
    a = list(b = 0.1, die = 0.05),
    b = list(die = 0.10)
  )
  # 'dead' is not the death name used in tp destinations ('die'); previously
  # this silently zeroed the death transition probabilities
  expect_error(
    validate_tp_source(tp, state_names = c("a", "b", "dead")),
    "'state_names' must match names\\(tp_source\\) in the same order"
  )
})

test_that("validate_tp_source catches cycle length mismatches", {
  bad_cyc <- list(
    a = list(b = rep(0.1, 10), c = rep(0.1, 10)),
    b = list(c = rep(0.1, 9))
  )
  expect_error(validate_tp_source(bad_cyc), "Length mismatch")
})

test_that("validate_tp_source catches invalid values", {
  bad_val <- list(
    a = list(b = 1.1, c = 0.1),
    b = list(c = -0.1)
  )
  expect_error(
    validate_tp_source(bad_val),
    "Invalid values.*Must be numeric and >= 0|Probabilities for state"
  )
})

test_that("validate_tp_source catches sum > 1", {
  bad_sum_valid_struct <- list(
    a = list(b = 0.6, c = 0.5),
    b = list(c = 0.1)
  )
  expect_error(
    validate_tp_source(bad_sum_valid_struct),
    "Probabilities for state 'a' sum to > 1"
  )
})

test_that("validate_tp_source should error on single state case", {
  single <- list(a = list(die = 0.5))
  expect_error(validate_tp_source(tp_source = single))
})

# ------------------------------------------------------------
# Variable-length tunnel tests
# ------------------------------------------------------------

test_that("validate_tp_source accepts variable-length tunnels", {
  # 1 pre-tunnel state, tunnel 1 length 3, tunnel 2 length 5
  # pre-tunnel vectors length 5 (= max, the model horizon)
  tp <- list(
    pre = list(
      tun1 = rep(0.10, 5),
      tun2 = rep(0.05, 5),
      die = rep(0.02, 5)
    ),
    tun1 = list(
      tun2 = rep(0.15, 3),
      die = rep(0.10, 3)
    ),
    tun2 = list(
      die = rep(0.20, 5)
    )
  )

  res <- validate_tp_source(tp, tunnel_lengths = c(3, 5))

  expect_equal(names(res), names(tp))
  # pre-tunnel vectors unchanged
  expect_length(res$pre$tun1, 5)
  # tunnel 1 vectors length 3
  expect_length(res$tun1$die, 3)
  # tunnel 2 vectors length 5
  expect_length(res$tun2$die, 5)
})

test_that("validate_tp_source fills NULL in variable-length tunnel with correct length", {
  tp <- list(
    pre = list(
      tun1 = rep(0.10, 6),
      tun2 = NULL, # should become rep(0, 6)
      die = rep(0.02, 6)
    ),
    tun1 = list(
      tun2 = NULL, # should become rep(0, 4)
      die = rep(0.10, 4)
    ),
    tun2 = list(
      die = rep(0.20, 7)
    )
  )

  res <- validate_tp_source(tp, tunnel_lengths = c(4, 7))

  expect_equal(res$pre$tun2, rep(0, 6))
  expect_equal(res$tun1$tun2, rep(0, 4))
})

test_that("validate_tp_source errors when tunnel vector has wrong length", {
  tp <- list(
    pre = list(
      tun1 = rep(0.10, 5),
      tun2 = rep(0.05, 5),
      die = rep(0.02, 5)
    ),
    tun1 = list(
      tun2 = rep(0.15, 3),
      die = rep(0.10, 3)
    ),
    tun2 = list(
      die = rep(0.20, 3) # wrong: tunnel 2 length should be 5, not 3
    )
  )

  expect_error(
    validate_tp_source(tp, tunnel_lengths = c(3, 5)),
    "Length mismatch"
  )
})

test_that("validate_tp_source errors when tunnel_lengths leaves no room for pre-tunnel", {
  tp <- list(
    a = list(b = rep(0.1, 5), die = rep(0.05, 5)),
    b = list(die = rep(0.2, 5))
  )

  # length(tunnel_lengths) == length(tp_source) — no pre-tunnel state
  expect_error(
    validate_tp_source(tp, tunnel_lengths = c(5, 5)),
    "need at least 1 pre-tunnel state"
  )
})

# ------------------------------------------------------------
# state_names argument
# ------------------------------------------------------------

test_that("validate_tp_source accepts and uses state_names argument", {
  th_cyc <- 5
  tp <- list(
    pre = list(
      tun1 = rep(0.1, th_cyc),
      tun2 = rep(0.05, th_cyc),
      die = rep(0.02, th_cyc)
    ),
    tun1 = list(tun2 = rep(0.1, th_cyc), die = rep(0.05, th_cyc)),
    tun2 = list(die = rep(0.1, th_cyc))
  )
  sn <- c("pre", "tun1", "tun2", "die")
  res <- validate_tp_source(tp, state_names = sn)

  # Canonical order is state_names order, excluding self
  expect_equal(names(res$tun1), c("pre", "tun2", "die"))
  expect_equal(names(res$tun2), c("pre", "tun1", "die"))

  # Backwards destinations filled with zeros
  expect_equal(res$tun1$pre, rep(0, th_cyc))
  expect_equal(res$tun2$pre, rep(0, th_cyc))
  expect_equal(res$tun2$tun1, rep(0, th_cyc))
})

test_that("validate_tp_source warns on explicitly supplied non-zero self-transition", {
  tp <- list(
    a = list(a = 0.1, b = 0.2, c = 0.05),
    b = list(c = 0.10)
  )
  expect_warning(
    validate_tp_source(tp),
    "Self-transition for state 'a'"
  )
})
