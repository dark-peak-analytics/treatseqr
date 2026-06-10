#' Validate and Sanitize Transition Probability List for Sequence Modelling
#'
#' Validates that the transition probability list (`tp_source`) follows the
#' strict structural rules required for sparse matrix extrapolation. It checks
#' for topological structure, which can be "reducing" or "uniform", consistency
#' of state names, consistent time horizons, and valid probability values. It
#' also sanitizes the list by replacing `NULL` entries with zero vectors of the
#' correct length, respecting fixed-duration (rather than time-horizon length)
#' tunnels if they are included.
#'
#' Specifically, this function ensures:
#'   \itemize{
#'     \item All inner lists are expanded to a uniform destination set (the
#'     canonical set) in canonical order. Missing destinations are inserted as
#'     NULL, which are subsequently replaced with zero vectors respecting
#'     `tunnel_lengths`.
#'     \item If state_names is supplied, it defines the canonical set and order.
#'     Otherwise it is inferred from tp_source.
#'     \item All probability vectors must have consistent lengths (see
#'     \code{tunnel_lengths} below).
#'     \item Probabilities must be non-negative and sum to <= 1 (row-wise).
#'     \item Any `NULL` transition vectors are replaced with numeric vectors of
#'     zeros of the correct length.
#'     \item A warning is emitted if a state's own name appears as a destination
#'     with a non-zero value, since p_stay is normally computed as 1 − sum(other
#'     transitions).
#'   }
#'
#' @param tp_source A named list of lists containing transition probability
#'   vectors. The structure must follow specific rules as described above.
#' @param tunnel_lengths An optional integer vector giving the length of each
#'   tunnel state block, in order (one element per tunnel). When \code{NULL}
#'   (default), all probability vectors must have the same length and the
#'   function behaves as in the uniform-tunnel case. When supplied, the number
#'   of pre-tunnel states is inferred as
#'   \code{length(tp_source) - length(tunnel_lengths)}, pre-tunnel TP vectors
#'   must all share a common length, and each tunnel \code{k}'s vectors must
#'   have length \code{tunnel_lengths[k]}.
#' @param state_names An optional string vector which is used as the "ground
#' truth" for validating the names of the inner lists in `tp_source`, in order.
#' Both the names of the states and their ordering will use this if supplied. If
#' not supplied, then state names AND ORDERING will be inferred from the names
#' of \code{tp_source} itself.The last element must be the absorbing (death)
#' state. All names present in \code{tp_source} must appear in state_names.
#'
#' @return The sanitized `tp_source` list where all `NULL` entries have been
#'   replaced with vectors of zeros of the correct length.
#' @export
validate_tp_source <- function(
  tp_source,
  tunnel_lengths = NULL,
  state_names = NULL
) {
  assertthat::assert_that(
    is.list(tp_source),
    !is.null(names(tp_source)),
    msg = "tp_source must be a named list"
  )

  n_states <- length(tp_source)
  assertthat::assert_that(
    n_states > 1,
    msg = paste0(
      "'tp_source' must have length greater than 1 otherwise there are no",
      " tunnel states. The death state should not be included as there is",
      " nowhere to transition to from dead"
    )
  )

  # All inner lists must be named
  assertthat::assert_that(
    all(vapply(tp_source, function(x) !is.null(names(x)), logical(1L))),
    msg = "All inner lists in tp_source must be named"
  )

  # Determine canonical state order
  if (!is.null(state_names)) {
    assertthat::assert_that(
      is.character(state_names),
      length(state_names) >= n_states + 1L,
      msg = paste0(
        "'state_names' must be a character vector with at least ",
        "length(tp_source) + 1 elements (all states plus the death state)"
      )
    )
    assertthat::assert_that(
      all(names(tp_source) %in% state_names),
      msg = "All names in 'tp_source' must be present in 'state_names'"
    )
    canonical_names <- state_names
  } else {
    assertthat::assert_that(
      length(tp_source[[1L]]) >= 1L,
      msg = "First element of tp_source must have at least one destination (the death state)"
    )
    # Infer: first-level names + last destination name of first state (= death)
    death_name <- tail(names(tp_source[[1L]]), 1L)
    canonical_names <- c(names(tp_source), death_name)
  }

  # Cycle through and check names. Warn if uer gives non-zero p_stay.
  for (state in seq_along(tp_source)) {
    self_name <- names(tp_source)[state]
    inner <- tp_source[[state]]
    if (self_name %in% names(inner)) {
      vec <- inner[[self_name]]
      if (!is.null(vec) && any(vec != 0)) {
        warning(sprintf(
          paste0(
            "Self-transition for state '%s' manually supplied with non-zero ",
            "values. p_stay is computed as 1 - sum(other transitions). ",
            "Verify that row sums are still <= 1."
          ),
          self_name
        ))
      }
    }
  }

  # Expand tp list to the full set for every state. Put NULL for empty.
  # Result is list elements include all except self.
  orig_names <- names(tp_source)
  tp_source <- lapply(
    setNames(seq_along(tp_source), orig_names),
    function(state) {
      self_name <- orig_names[state]
      dest_names <- canonical_names[canonical_names != self_name]
      inner <- tp_source[[state]]
      expanded <- lapply(dest_names, function(nm) inner[[nm]])
      names(expanded) <- dest_names
      expanded
    }
  )

  # Determine expected vector lengths per state. respects tunnel_lengths if
  # supplied. Otherwise all are assumed to be time horizon length.
  if (is.null(tunnel_lengths)) {
    # Uniform case: all vectors must share a single common length. Infer it
    # from the first non-null, non-empty vector.
    all_vectors <- unlist(tp_source, recursive = FALSE)
    vec_lens <- lengths(all_vectors)
    n_cyc <- vec_lens[vec_lens > 0L][1L]
    assertthat::assert_that(
      !is.na(n_cyc),
      msg = "Could not determine number of cycles (all vectors are NULL or empty)"
    )
    # One expected length per state
    expected_vec_lengths <- rep(n_cyc, n_states)
  } else {
    # Variable-length case.
    n_tunnels <- length(tunnel_lengths)
    assertthat::assert_that(
      n_tunnels >= 1,
      n_tunnels <= n_states - 1,
      msg = sprintf(
        paste0(
          "length(tunnel_lengths) must be >= 1 and <= length(tp_source) - 1.",
          " Got %d tunnels for %d states (need at least 1 pre-tunnel state)."
        ),
        n_tunnels,
        n_states
      )
    )

    pre_tunnels <- n_states - n_tunnels

    # Infer th (pre-tunnel horizon) from first non-null, non-empty vector in
    # the pre-tunnel states.
    pre_tun_vectors <- unlist(
      tp_source[seq_len(pre_tunnels)],
      recursive = FALSE
    )
    vec_lens <- lengths(pre_tun_vectors)
    th <- vec_lens[vec_lens > 0L][1L]
    assertthat::assert_that(
      !is.na(th),
      msg = "Could not determine pre-tunnel horizon (all pre-tunnel vectors are NULL or empty)"
    )

    # Pre-tunnel states all share th; tunnel k has length tunnel_lengths[k]
    expected_vec_lengths <- c(rep(th, pre_tunnels), tunnel_lengths)
  }

  # Sanitize: validate lengths/values and replace NULLs
  sanitized_tp <- lapply(seq_along(tp_source), function(state_idx) {
    state_name <- names(tp_source)[state_idx]
    state_list <- tp_source[[state_idx]]
    exp_len <- expected_vec_lengths[state_idx]

    new_state_list <- lapply(names(state_list), function(trans_name) {
      vec <- state_list[[trans_name]]

      if (is.null(vec)) {
        return(numeric(exp_len))
      }

      assertthat::assert_that(
        length(vec) == exp_len,
        msg = sprintf(
          "Length mismatch in transition %s -> %s. Expected %d, found %d.",
          state_name,
          trans_name,
          exp_len,
          length(vec)
        )
      )

      assertthat::assert_that(
        is.numeric(vec),
        all(vec >= 0),
        msg = sprintf(
          "Invalid values in transition %s -> %s. Must be numeric and >= 0.",
          state_name,
          trans_name
        )
      )

      return(vec)
    })
    names(new_state_list) <- names(state_list)
    return(new_state_list)
  })
  names(sanitized_tp) <- names(tp_source)

  # Summation constraint: element-wise sum of all transitions out of a state
  # must be <= 1 for every cycle. Compute all row-sums at once, then report
  # the first violating state.
  total_probs <- lapply(sanitized_tp, function(state_list) {
    Reduce(`+`, state_list)
  })
  violations <- which(vapply(
    total_probs,
    function(tp) any(tp > 1 + 1e-10),
    logical(1L)
  ))
  if (length(violations) > 0L) {
    stop(sprintf(
      "Probabilities for state '%s' sum to > 1 in some cycles.",
      names(sanitized_tp)[violations[1L]]
    ))
  }

  sanitized_tp
}
