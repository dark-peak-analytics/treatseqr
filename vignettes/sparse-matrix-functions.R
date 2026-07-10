## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(treatseqr)

## ----specify_m----------------------------------------------------------------
m_spec <- treatseqr::specify_m(
  n_cycles = 5,
  n_tunnels = 3,
  pre_tunnel_states = 2
)

## ----print_spec---------------------------------------------------------------
print(m_spec$m1_ijx)

## ----m2_example_first_20------------------------------------------------------
print(m_spec$m2_ijx[1:20, ])

## ----m2_example_top_5---------------------------------------------------------
print(m_spec$m2_ijx[1:5, ])

## ----m2_example_tunnel1_to_2--------------------------------------------------
print(m_spec$m2_ijx[6:10, ])

## ----m2_example_tunnel1_to_3--------------------------------------------------
print(m_spec$m2_ijx[11:15, ])

## ----m2_example_tunnel1_to_death----------------------------------------------
print(m_spec$m2_ijx[16:20, ])

## ----m2_example_tunnel2-------------------------------------------------------
print(m_spec$m2_ijx[21:35, ])

## ----m2_example_tunnel3-------------------------------------------------------
print(m_spec$m2_ijx[36:45, ])

## ----m2_example_death---------------------------------------------------------
print(m_spec$m2_ijx[46, ])

## ----transition_probabilities-------------------------------------------------
th_cyc <- 50

# some simple time varying transition probabilities. in reality these will come
# from hazard extrapolations e.g. from `flexsurv`
tp_source <- list(
  on_1l = list(
    off_1l = seq(0.1, 0.05, length.out = th_cyc), # Example: linearly decreasing
    on_2l = seq(0.05, 0.1, length.out = th_cyc), # Example: linearly increasing
    off_2l = rep(0, th_cyc),
    bsc = rep(0, th_cyc),
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

## ----view_data_tree-----------------------------------------------------------
lobstr::tree(
  rapply(tp_source, function(x) round(x, 2), how = "list")
)

## ----validate_transition_probabilities----------------------------------------
valid_tp <- validate_tp_source(tp_source = tp_source)
lobstr::tree(
  rapply(valid_tp, function(x) round(x, 2), how = "list")
)

## ----transition_probabilities_with_nulls--------------------------------------
th_cyc <- 50

# some simple time varying transition probabilities. in reality these will come
# from hazard extrapolations e.g. from `flexsurv`
tp_source <- list(
  on_1l = list(
    off_1l = seq(0.1, 0.05, length.out = th_cyc), # Example: linearly decreasing
    on_2l = seq(0.05, 0.1, length.out = th_cyc), # Example: linearly increasing
    off_2l = NULL,
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
lobstr::tree(
  rapply(tp_source, function(x) round(x, 2), how = "list")
)

## ----validate_transition_probabilities_with_nulls-----------------------------
# Fills in the blanks if they are NULLs
valid_tp <- validate_tp_source(tp_source = tp_source)
lobstr::tree(
  rapply(valid_tp, function(x) round(x, 2), how = "list")
)

## ----specify_m_from_valid_tp--------------------------------------------------
m_spec <- specify_m(
  n_cycles = 50,
  n_tunnels = length(valid_tp) - 1,
  pre_tunnel_states = 1
)

## ----derive_m1_data-----------------------------------------------------------
p_off_1l <- valid_tp$on_1l$off_1l
p_on_2l <- valid_tp$on_1l$on_2l
p_off_2l <- valid_tp$on_1l$off_2l
p_bsc <- valid_tp$on_1l$bsc
p_die <- valid_tp$on_1l$die
m_size <- m_spec$matrix_size
m1_coord <- m_spec$m1_ijx

## ----print_m1_ijx-------------------------------------------------------------
print(m_spec$m1_ijx)

## ----manual_m1_list-----------------------------------------------------------
m1_list <- lapply(seq_len(th_cyc), function(model_cycle) {
  m1_coord[, "x"] <- c(
    1 - sum(
      p_off_1l[model_cycle],
      p_on_2l[model_cycle],
      p_off_2l[model_cycle],
      p_bsc[model_cycle],
      p_die[model_cycle]
    ),
    p_off_1l[model_cycle],
    p_on_2l[model_cycle],
    p_off_2l[model_cycle],
    p_bsc[model_cycle],
    p_die[model_cycle]
  )
  m1_coord <- m1_coord[m1_coord[, "x"] > 0, , drop = FALSE]
  Matrix::sparseMatrix(
    i = m1_coord[, "i"],
    j = m1_coord[, "j"],
    x = m1_coord[, "x"],
    dims = c(1, m_size)
  )
})
print(m1_list[[1]])

## ----assert_m1_sums-----------------------------------------------------------
assertthat::assert_that(
  all(unlist(lapply(m1_list, function(x) round(sum(x), 1e-14))) == 1),
  msg = "m1 rows do not sum to 1"
)

## ----generate_m_list----------------------------------------------------------
# validate inputs and generate the full matrix m from one function
m <- generate_m_list(
  m_specification = m_spec,
  transition_prob_list = tp_source
)

## ----extrapolate_manual-------------------------------------------------------
# output population matrix:
pop <- matrix(0, nrow = m_spec$matrix_size, ncol = th_cyc + 1)
pop[1, 1] <- 1

# get the TPs
m1 <- m$m1
m2 <- m$m2

# extrapolate the model:
for (cyc in seq_len(th_cyc) + 1) {
  pop[, cyc] <-
    as.numeric(pop[1, cyc - 1] %*% m1[[cyc - 1]]) +
    as.numeric(pop[, cyc - 1] %*% m2)
}

## ----plot_os------------------------------------------------------------------
plot(1 - pop[nrow(pop), ],
  type = "l", xlab = "Cycle", ylab = "Overall Survival",
  main = "Overall Survival Curve"
)
