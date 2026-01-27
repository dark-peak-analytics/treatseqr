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
    a = list(b = -0.1)
  )
  expect_error(
    validate_tp_source(bad_val),
    "Invalid values.*Must be numeric and >= 0"
  )
})

test_that("validate_tp_source catches sum > 1", {
  bad_sum <- list(
    a = list(b = 0.6, c = 0.5) # Sums to 1.1
  )
  # NOTE: validate_tp_source requires the last element to have length 1.
  # If we only provide one element, it satisfies length 1 but might fail other checks
  # if we intend it to be a multi-state example.
  # Let's make a valid 1-state list that fails the sum check.

  bad_sum_single <- list(
    a = list(b = 0.6, c = 0.5)
  )
  # But wait, if length is 1, n_states=1. element_lengths[1]=2.
  # expected_lengths would be seq(2, by=-1, length=1) = 2.
  # But we also assert last element has length 1.
  # So a single-element list MUST have length 1 to pass the length check.
  # So to test sum > 1, we need a list where the sum check fails but structure passes.

  bad_sum_valid_struct <- list(
    a = list(b = 0.6, c = 0.5),
    b = list(c = 0.1)
  )
  expect_error(
    validate_tp_source(bad_sum_valid_struct),
    "Probabilities for state 'a' sum to > 1"
  )
})

test_that("validate_tp_source handles single state case? (Edge case, might not be valid for sequence model)", {
  # A sequence model usually implies transition between states, but if length is 1...
  # The function requires reducing to 1. If length is 1, it satisfies the condition.
  single <- list(a = list(die = 0.5))
  res <- validate_tp_source(single)
  expect_equal(res$a$die, 0.5)
})
