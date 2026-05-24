test_that("validate_tp_source works with valid input (vignette example)", {
  th_cyc <- 10

  tp_source <- list(
    on_1l = list(
      off_1l = seq(0.1, 0.05, length.out = th_cyc),
      on_2l = seq(0.05, 0.1, length.out = th_cyc),
      off_2l = rep(0, th_cyc),
      bsc = NULL, # Should be converted to 0
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

  # Check structure integrity
  expect_equal(names(res), names(tp_source))
  expect_equal(length(res), 5)
  expect_equal(length(res$on_1l), 5)
  expect_equal(length(res$bsc), 1)
})

test_that("validate_tp_source catches length errors", {
  # Last element is 1, but the sequence is wrong (3, 3, 1 instead of 3, 2, 1)
  bad_lengths <- list(
    a = list(b = 1, c = 1, d = 1),
    b = list(c = 1, d = 1, e = 1), # Error: length 3, expected 2
    c = list(d = 1)
  )
  expect_error(
    validate_tp_source(bad_lengths),
    "Lengths of tp_source elements must reduce by 1"
  )
})

test_that("validate_tp_source catches topological errors", {
  # 'b' has transition to 'z', which 'a' does not have
  bad_topo <- list(
    a = list(b = 1, c = 1),
    b = list(z = 1)
  )
  expect_error(
    validate_tp_source(bad_topo),
    "Topology error: State 'b' contains transitions to \\{z\\} which are not present in the previous state 'a'"
  )
})

test_that("validate_tp_source catches cycle length mismatches", {
  bad_cyc <- list(
    a = list(b = rep(1, 10), c = rep(1, 10)),
    b = list(c = rep(1, 9)) # Length 9 instead of 10
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
      die  = rep(0.02, 5)
    ),
    tun1 = list(
      tun2 = rep(0.15, 3),
      die  = rep(0.10, 3)
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
      tun2 = NULL,        # should become rep(0, 6)
      die  = rep(0.02, 6)
    ),
    tun1 = list(
      tun2 = NULL,        # should become rep(0, 4)
      die  = rep(0.10, 4)
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
      die  = rep(0.02, 5)
    ),
    tun1 = list(
      tun2 = rep(0.15, 3),  # correct: tunnel 1 length is 3
      die  = rep(0.10, 3)
    ),
    tun2 = list(
      die = rep(0.20, 3)    # wrong: tunnel 2 length should be 5, not 3
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
