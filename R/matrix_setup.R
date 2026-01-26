#' Calculate Matrix Size and Top Row Indices for Sequence Model
#'
#' Computes the indices of the top row (`m1`) of the main transition matrix and
#' the overall matrix size, given the number of cycles and tunnel states.
#'
#' @param n_cycles Integer. Number of cycles in the model (must be > 0).
#' @param n_tunnels Integer. Number of tunnel states in the model (must be > 0).
#'
#' @return A list with:
#'   \describe{
#'     \item{m1}{Integer vector of column indices for the top row of the
#'     matrix.}
#'     \item{matrix_size}{Integer. The total number of columns (states) in the
#'     matrix.}
#'   }
#' @examples
#' m_size(5, 2)
#' @export
m_size <- function(n_cycles, n_tunnels) {
  assertthat::assert_that(
    is.numeric(n_cycles),
    length(n_cycles) == 1,
    n_cycles > 0,
    n_cycles == as.integer(n_cycles),
    msg = "n_cycles must be a single positive integer"
  )
  assertthat::assert_that(
    is.numeric(n_tunnels),
    length(n_tunnels) == 1,
    n_tunnels > 0,
    n_tunnels == as.integer(n_tunnels),
    msg = "n_tunnels must be a single positive integer"
  )
  # m1_columns gives the indices for the top row of the transition matrix:
  # 1 (stay in first state), then each tunnel state entry (every n_cycles
  #   steps).
  m1_columns <- c(1, seq(2, (n_tunnels * n_cycles) + 1, by = n_cycles))
  matrix_size <- max(m1_columns)

  list(
    m1 = m1_columns,
    matrix_size = matrix_size
  )
}
