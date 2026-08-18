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
#'     \item If state_names is supplied, it is used purely as a cross-check:
#'     it must equal names(tp_source) in order, with the death state appended
#'     last. Otherwise the canonical set and order are inferred from tp_source.
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
#' @param state_names An optional character vector used to validate
#' \code{tp_source}. It must equal \code{names(tp_source)} in the same order,
#' with the absorbing (death) state appended as the last element. It cannot be
#' used to reorder states: the matrix layout produced by \code{specify_m()} is
#' positional. If not supplied, state names and ordering are inferred from
#' \code{tp_source} itself.
#'
#' @return The sanitized `tp_source` list where all `NULL` entries have been
#'   replaced with vectors of zeros of the correct length.
#' @examples
#' # 1 pre-tunnel state (pre) and 2 tunnels of 3 cycles each
#' tp <- list(
#'   pre  = list(tun1 = rep(0.2, 3), tun2 = rep(0.1, 3), die = rep(0.05, 3)),
#'   tun1 = list(tun2 = rep(0.3, 3), die = rep(0.1, 3)),
#'   tun2 = list(die = rep(0.2, 3))
#' )
#' state_names <- c("pre", "tun1", "tun2", "die")
#'
#' valid_tp <- validate_tp_source(
#'   tp_source      = tp,
#'   tunnel_lengths = c(3, 3),
#'   state_names    = state_names
#' )
#'
#' # every state is expanded to the full destination set, so the backwards
#' # slot tun1 -> pre is added as a zero vector
#' names(valid_tp$tun1)
#' valid_tp$tun1$pre
#' @importFrom utils tail
#' @importFrom stats setNames
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

  if (!is.null(tunnel_lengths)) {
    if (any(tunnel_lengths == 1)) {
      n_pre <- length(tp_source) - length(tunnel_lengths)
      which_one_tunnel <- which(tunnel_lengths == 1) + n_pre
      msg <- paste0(
        "Tunnel lengths of 1 are allowed, but the p_stay will be 0. Patients ",
        "will transition out of the tunnel after 1 cycle. You have assigned ",
        paste(sQuote(names(tp_source)[which_one_tunnel]), collapse = ",  "),
        " to have a tunnel length of 1."
      )
      message(msg)
    }
  }

  # Death name makes sense. it is the only thing in tp_source that doesn't
  # have a state (as it is absorbing)
  all_dests <- unique(unlist(lapply(tp_source, names)))
  death_name <- setdiff(all_dests, names(tp_source))
  assertthat::assert_that(
    length(death_name) == 1L,
    msg = "Death state does not appear only as a destination (absorbing)"
  )

  # Determine canonical state order
  if (!is.null(state_names)) {
    assertthat::assert_that(
      is.character(state_names),
      length(state_names) == n_states + 1L,
      msg = paste0(
        "'state_names' must have exactly length(tp_source) + 1 elements ",
        "(all transient states plus the death state)"
      )
    )
    assertthat::assert_that(
      identical(state_names, c(names(tp_source), death_name)),
      msg = paste0(
        "'state_names' must match names(tp_source) in the same order, with the ",
        "absorbing (death) state appended last. The matrix layout produced by ",
        "specify_m() is positional and cannot be reordered via state_names."
      )
    )

    canonical_names <- state_names
  } else {
    assertthat::assert_that(
      length(tp_source[[1L]]) >= 1L,
      msg = "First element of tp_source must have at least one destination (the death state)"
    )

    # Because user hasn't given the state names explicitly we must infer. We
    # have already deduced the name of the death transition so we now must
    # validate that it is in all elements of tp_source, and is the last element
    # in each.
    assertthat::assert_that(
      all(vapply(
        tp_source,
        function(x) {
          death_name == names(x)[length(names(x))]
        },
        logical(1L)
      )),
      msg = paste0(
        "The state ",
        sQuote(death_name),
        " must be the LAST element in each element of ",
        sQuote("tp_source")
      )
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
