# Test the m_size() function for correct output structure and values.
# - Checks that the return type is a list of length 2.
# - Verifies that the m1 vector contains the correct indices for the top row of
# the transition matrix.
# - Ensures that the matrix_size is correct for the given number of cycles and
# tunnel states.
# - Confirms that invalid input (n_cycles = 0) triggers an error.
test_that("m_size returns correct values", {
  cols <- m_size(n_cycles = 5, n_tunnels = 5)
  expect_type(cols, "list")
  expect_equal(length(cols), 2)
  expect_equal(
    cols$m1_ijx,
    structure(
      c(1, 1, 1, 1, 1, 1, 1, 1, 2, 7, 12, 17, 22, 28, 0, 0, 0, 0, 0, 0, 0),
      dim = c(7L, 3L),
      dimnames = list(NULL, c("i", "j", "x"))
    )
  )
  expect_equal(cols$matrix_size, 28)
  expect_error(m_size(n_cycles = 0))
})

test_that("m_size handles minimal valid input", {
  cols <- m_size(n_cycles = 1, n_tunnels = 1)
  expect_equal(
    cols$m1_ijx,
    structure(
      c(1, 1, 1, 1, 2, 4, 0, 0, 0),
      dim = c(3L, 3L),
      dimnames = list(
        NULL,
        c("i", "j", "x")
      )
    )
  )
  expect_equal(cols$matrix_size, 4)
})

test_that("m_size errors on non-integer input", {
  expect_error(
    m_size(n_cycles = 2.5, n_tunnels = 2),
    "n_cycles must be a single positive integer"
  )
  expect_error(
    m_size(n_cycles = 2, n_tunnels = 2.5),
    "n_tunnels must be a single positive integer"
  )
})

test_that("m_size errors on negative input", {
  expect_error(
    m_size(n_cycles = -1, n_tunnels = 2),
    "n_cycles must be a single positive integer"
  )
  expect_error(
    m_size(n_cycles = 2, n_tunnels = -2),
    "n_tunnels must be a single positive integer"
  )
})

test_that("m_size errors on vector input", {
  expect_error(
    m_size(n_cycles = c(2, 3), n_tunnels = 2),
    "n_cycles must be a single positive integer"
  )
  expect_error(
    m_size(n_cycles = 2, n_tunnels = c(2, 3)),
    "n_tunnels must be a single positive integer"
  )
})

test_that("m_size works for lots of pre-tunnel states", {
  cols <- m_size(n_cycles = 10, n_tunnels = 3, pre_tunnel_states = 10)
  expect_equal(nrow(cols$m1_ijx), 95)
  expect_false(
    any(cols$m1_ijx[cols$m1_ijx[, 1] == 10, "j"] < 10),
    label = paste0(
      "The last pre-tunnel state should not have any pre-tunnel",
      " transitions except to itself."
    )
  )
})

test_that("m_size works for large valid input", {
  cols <- m_size(n_cycles = 100, n_tunnels = 10)
  expect_equal(nrow(cols$m1_ijx), 12)
  expect_equal(
    cols$m1_ijx,
    structure(
      c(
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        102,
        202,
        302,
        402,
        502,
        602,
        702,
        802,
        902,
        1003,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0
      ),
      dim = c(12L, 3L),
      dimnames = list(NULL, c("i", "j", "x"))
    )
  )
  expect_equal(cols$matrix_size, 1003)
})

test_that("m_size errors if arguments are missing", {
  expect_error(m_size(n_cycles = 2), "argument \"n_tunnels\" is missing")
  expect_error(m_size(n_tunnels = 2), "argument \"n_cycles\" is missing")
})
