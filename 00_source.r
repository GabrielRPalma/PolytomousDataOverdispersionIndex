####################################################################################################
###
### File:    00_source.r
### Purpose: Load the required packages and define the main functions of the
###          project (simulation of longitudinal polytomous data, observed and
###          expected variances, dispersion index, simulation study and case
###          study helpers) together with the plot settings shared by all scripts.
### Authors: Gabriel Rodrigues Palma
### Date:    19/03/23
###
####################################################################################################

# Packages required -----
packages <- c("nnet",    # multinom(): multinomial log-linear models
              "mclogit", # mblogit(): multinomial logit models with random effects
              "dplyr",   # Data processing
              "tidyr",   # Data reshaping (pivot_longer / pivot_wider)
              "ggplot2") # Data visualisation

install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

# Main functions -----
# Function index (search the name to jump to its definition):
#
#   Simulating longitudinal polytomous data
#     linear_predictors_to_probabilities(linear_predictors)
#     create_prob_vector_grouped(n, j, theta)
#     create_prob_vector_individual(n, j, theta)
#     obtain_longitudinal_data_per_individual(id, j, k, theta, m, prob_vector_fun)
#     create_longitudinal_grouped_polytomous_data(n, j, k, theta, m)
#     create_longitudinal_polytomous_data(n, j, k, theta)
#
#   Observed and expected variances (per category and time)
#     get_observed_variance(data, time, categories)
#     get_expected_variance_from_probabilities(probabilities, m)
#     get_expected_variance_closed_form(data, time, categories, m)
#     make_multinomial_formula(j, response_prefix, predictor)
#     newdata_observed_covariate(data)
#     get_multinomial_variance(data, formula, newdata, m)
#     get_expected_variance_multinomial(data, time, formula, newdata_fun, m)
#     get_expected_variance_mixed(data, time, formula, random, m, conditional)
#
#   Dispersion index (three forms of the expected variance)
#     compute_dispersion_index(var_obs, var_exp, categories)
#     get_dispersion_index_closed_form(data, time, categories, m)
#     get_dispersion_index_multinomial(data, time, categories, m, formula, newdata_fun)
#     get_dispersion_index_mixed(data, time, categories, m, formula, random, conditional)
#
#   Simulation study
#     run_dispersion_index_replications(n, j, k, m, theta, theta_label, n_simulations)
#     get_simulation_results(n, j, k, m, n_simulations, theta_super)
#     get_index_comparison_results(n, j, k, m, theta, n_simulations, random)
#     get_auc(x, reference)
#     get_discrimination_summary(results, theta_reference, false_alarm_rate)
#
#   Case study helpers
#     get_probabilities(etas)
#
#   Plot settings
#     theme_new(base_size, base_family)
#     plot_simulation_results(final_results, n_value)
#     index_form_labels
#     label_comparison_factors(data)
#     plot_index_comparison(results)
#     plot_discrimination_summary(discrimination_summary)

## Simulating longitudinal polytomous data -----
linear_predictors_to_probabilities <- function(linear_predictors) {
  # This function converts the linear predictors of the non-baseline categories
  # into the probabilities of all the categories (baseline category first).
  # Input:
  #     linear_predictors: n x (j - 1) matrix with the linear predictors z2, ..., zj
  # Output:
  #     prob: n x j matrix with the probabilities p1, ..., pj
  exp_linear_predictors <- exp(linear_predictors)
  den <- 1 + rowSums(exp_linear_predictors)
  prob <- cbind(1, exp_linear_predictors) / den
  colnames(prob) <- paste0("p", seq_len(ncol(prob)))
  prob
}

create_prob_vector_grouped <- function(n, j, theta) {
  # This function creates the vector of category probabilities used to simulate
  # grouped polytomous data. The random effect u is shared by all the linear
  # predictors.
  # Input:
  #     n: number of covariate values to draw
  #     j: number of categories (3, 4 or 5)
  #     theta: standard deviation of the random effect
  # Output:
  #     result: list with prob (n x j matrix of probabilities) and x (covariate)
  stopifnot(j %in% 3:5)
  x <- rnorm(n)                       # Covariate
  u <- rnorm(1, mean = 0, sd = theta) # Random effect
  # NOTE: z2 was originally written 1 + 0.5**x + u (i.e. 0.5^x) and was
  # corrected to 0.5 * x on 20/08/26; output_data/final_results.csv and the
  # simulation plots produced before that date used 0.5^x.
  linear_predictors <- cbind(z2 = 1 + 0.5 * x + u,
                             z3 = 0.5 + x + u,
                             z4 = 1.5 - x + u,
                             z5 = 1 - 0.7 * x + u)
  prob <- linear_predictors_to_probabilities(linear_predictors[, seq_len(j - 1), drop = FALSE])
  list(prob = prob, x = x)
}

create_prob_vector_individual <- function(n, j, theta) {
  # This function creates the vector of category probabilities used to simulate
  # individual polytomous data. The random effect u enters only the linear
  # predictor of the second category.
  # Input:
  #     n: number of covariate values to draw
  #     j: number of categories (3, 4 or 5)
  #     theta: standard deviation of the random effect
  # Output:
  #     result: list with prob (n x j matrix of probabilities) and x (covariate)
  stopifnot(j %in% 3:5)
  x <- rnorm(n)                       # Covariate
  u <- rnorm(1, mean = 0, sd = theta) # Random effect
  linear_predictors <- cbind(z2 = 1.5 - 3 * x + u,
                             z3 = 3 - 5 * x,
                             z4 = 2 - 4 * x,
                             z5 = 4 - 7 * x)
  prob <- linear_predictors_to_probabilities(linear_predictors[, seq_len(j - 1), drop = FALSE])
  list(prob = prob, x = x)
}

obtain_longitudinal_data_per_individual <- function(id, j, k, theta, m = 1,
                                                    prob_vector_fun = create_prob_vector_individual) {
  # This function creates the longitudinal dataset of one individual (or group
  # of individuals) evaluated at k periods.
  # Input:
  #     id: identifier of the individual (or group)
  #     j: number of categories
  #     k: number of evaluated periods
  #     theta: standard deviation of the random effect
  #     m: number of individuals in the group (1 for individual data)
  #     prob_vector_fun: function used to create the category probabilities
  #                      (create_prob_vector_individual or create_prob_vector_grouped)
  # Output:
  #     longitudinal_data: data frame with the columns id, X, Y.1, ..., Y.j and Time
  longitudinal_data <- lapply(seq_len(k), function(period) {
    polytomous_data <- prob_vector_fun(n = 1, j = j, theta = theta)
    y <- t(apply(polytomous_data$prob, MARGIN = 1, FUN = rmultinom, n = 1, size = m))
    data.frame(id = id, X = polytomous_data$x, Y = y)
  })
  longitudinal_data <- do.call(rbind, longitudinal_data)
  longitudinal_data$Time <- seq_len(k)
  longitudinal_data
}

create_longitudinal_grouped_polytomous_data <- function(n, j, k, theta, m) {
  # This function creates a dataset of grouped longitudinal polytomous data
  # (n groups of m individuals evaluated at k periods).
  # Input:
  #     n: number of groups
  #     j: number of categories
  #     k: number of evaluated periods
  #     theta: standard deviation of the random effect
  #     m: number of individuals in each group
  # Output:
  #     longitudinal_data: data frame with the columns id, X, Y.1, ..., Y.j and Time
  longitudinal_data <- lapply(seq_len(n), function(id) {
    obtain_longitudinal_data_per_individual(id = id, j = j, k = k, theta = theta, m = m,
                                            prob_vector_fun = create_prob_vector_grouped)
  })
  do.call(rbind, longitudinal_data)
}

create_longitudinal_polytomous_data <- function(n, j, k, theta) {
  # This function creates a dataset of individual longitudinal polytomous data
  # (n individuals evaluated at k periods, m = 1).
  # Input:
  #     n: number of individuals
  #     j: number of categories
  #     k: number of evaluated periods
  #     theta: standard deviation of the random effect
  # Output:
  #     longitudinal_data: data frame with the columns id, X, Y.1, ..., Y.j and Time
  longitudinal_data <- lapply(seq_len(n), function(id) {
    obtain_longitudinal_data_per_individual(id = id, j = j, k = k, theta = theta, m = 1,
                                            prob_vector_fun = create_prob_vector_individual)
  })
  do.call(rbind, longitudinal_data)
}

## Observed and expected variances -----
get_observed_variance <- function(data, time, categories) {
  # This function computes the observed variance of each category at each time.
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     categories: names of the category columns
  # Output:
  #     var_obs: tibble with one row per time and one column per category
  data %>%
    pivot_longer(cols = all_of(categories), names_to = "category", values_to = "value") %>%
    group_by(.data[[time]], category) %>%
    summarise(var = var(value), .groups = "drop") %>%
    pivot_wider(names_from = category, values_from = var)
}

get_expected_variance_from_probabilities <- function(probabilities, m) {
  # This function computes the variance expected for grouped polytomous data
  # given the category probabilities: m * p * (1 - p).
  # Input:
  #     probabilities: vector or matrix of category probabilities (one column per category)
  #     m: number of individuals in each group
  # Output:
  #     expected_variance: expected variance with the same shape as probabilities
  probabilities * m * (1 - probabilities)
}

get_expected_variance_closed_form <- function(data, time, categories, m) {
  # This function computes the expected variance of each category at each time
  # with the closed form expression mean(y) * (m - mean(y)) / m, i.e. the
  # multinomial variance with the probabilities estimated by mean(y) / m.
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     categories: names of the category columns
  #     m: number of individuals in each group
  # Output:
  #     var_exp: tibble with one row per time and one column per category
  data %>%
    pivot_longer(cols = all_of(categories), names_to = "category", values_to = "value") %>%
    group_by(.data[[time]], category) %>%
    summarise(var = mean(value) * (.env$m - mean(value)) / .env$m, .groups = "drop") %>%
    pivot_wider(names_from = category, values_from = var)
}

make_multinomial_formula <- function(j, response_prefix = "Y", predictor = "X") {
  # This function builds the multinomial model formula for the simulated data
  # layout: cbind(Y.1, ..., Y.j) ~ X.
  # Input:
  #     j: number of response categories (Y.1, ..., Y.j)
  #     response_prefix: prefix of the response variables (default "Y")
  #     predictor: right-hand side of the formula (default "X")
  # Output:
  #     formula: formula object of the form cbind(Y.1, ..., Y.j) ~ X
  responses <- paste0(response_prefix, ".", seq_len(j))
  lhs <- paste0("cbind(", paste(responses, collapse = ", "), ")")
  as.formula(paste(lhs, "~", predictor))
}

newdata_observed_covariate <- function(data) {
  # This function builds the prediction grid for the simulated data layout:
  # the observed values of the covariate X.
  # Input:
  #     data: data frame with the column X
  # Output:
  #     newdata: data frame with the column X
  data.frame(X = data$X)
}

get_multinomial_variance <- function(data, formula, newdata, m) {
  # This function fits a multinomial model and computes the expected variance
  # of each category as the mean of m * p * (1 - p) over the predicted
  # probabilities p of the prediction grid.
  # Input:
  #     data: data frame in wide format (one column per category)
  #     formula: multinomial model formula, e.g. cbind(Y.1, Y.2, Y.3) ~ X
  #     newdata: data frame with the values of the predictors used to predict
  #     m: number of individuals in each group
  # Output:
  #     variance: named vector with the expected variance of each category
  fit <- multinom(formula, data = data, trace = FALSE)
  predicted_probabilities <- predict(fit, newdata = newdata, type = "probs")
  expected_variance <- get_expected_variance_from_probabilities(predicted_probabilities, m)
  apply(expected_variance, MARGIN = 2, FUN = mean)
}

get_expected_variance_multinomial <- function(data, time, formula, newdata_fun, m) {
  # This function computes the expected variance of each category at each
  # time from a multinomial model fitted separately at each time.
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     formula: multinomial model formula, e.g. cbind(Y.1, Y.2, Y.3) ~ X
  #     newdata_fun: function that receives the data of one time and returns
  #                  the prediction grid (e.g. newdata_observed_covariate)
  #     m: number of individuals in each group
  # Output:
  #     var_exp: tibble with one row per time and one column per category
  time_groups <- split(data, data[[time]])
  lapply(names(time_groups), function(group_name) {
    group_data <- time_groups[[group_name]]
    variance <- get_multinomial_variance(data = group_data, formula = formula,
                                         newdata = newdata_fun(group_data), m = m)
    c(setNames(list(group_name), time), as.list(variance))
  }) %>%
    bind_rows()
}

get_expected_variance_mixed <- function(data, time, formula, 
                                        random, m, conditional = TRUE) {
  # This function computes the expected variance of each category at each
  # time from a single random effects multinomial model (mclogit::mblogit)
  # fitted to the whole longitudinal data set. The probabilities of every
  # observation are predicted from the fitted model and the expected variance
  # m * p * (1 - p) is averaged over the observations of each time. Include
  # the time variable in the fixed effects formula to obtain probabilities
  # that change over time.
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     formula: fixed effects formula, e.g. cbind(Y.1, Y.2, Y.3) ~ X
  #     random: random effects formula, e.g. ~ 1 | id
  #     m: number of individuals in each group
  #     conditional: if TRUE (default) the predicted probabilities include the
  #                  estimated random effects of each subject; if FALSE they
  #                  are the population-level probabilities (random effects
  #                  set to zero)
  # Output:
  #     var_exp: tibble with one row per time and one column per category
  fit <- mblogit(formula, data = data, random = random,
                 control = mmclogit.control(trace = FALSE))
  predicted_probabilities <- predict(fit, type = "response", conditional = conditional)
  stopifnot(nrow(predicted_probabilities) == nrow(data))
  expected_variance <- get_expected_variance_from_probabilities(predicted_probabilities, m)
  expected_variance <- as.data.frame(expected_variance)
  expected_variance[[time]] <- data[[time]]
  expected_variance %>%
    group_by(.data[[time]]) %>%
    summarise(across(everything(), mean), .groups = "drop")
}

## Dispersion index -----
compute_dispersion_index <- function(var_obs, var_exp, categories) {
  # This function computes the proposed dispersion index: the ratio between the
  # observed and the expected variances averaged over time and then over the
  # categories.
  # Input:
  #     var_obs: observed variances (one row per time, one column per category)
  #     var_exp: expected variances (one row per time, one column per category)
  #     categories: names of the category columns
  # Output:
  #     index: value of the dispersion index
  stopifnot(nrow(var_obs) == nrow(var_exp))
  ratio <- as.matrix(var_obs[, categories]) / as.matrix(var_exp[, categories])
  mean(apply(ratio, MARGIN = 2, FUN = mean))
}

get_dispersion_index_closed_form <- function(data, time, categories, m) {
  # This function computes the proposed dispersion index with the expected
  # variance obtained from the closed form expression (see
  # get_expected_variance_closed_form).
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     categories: names of the category columns
  #     m: number of individuals in each group
  # Output:
  #     index: value of the dispersion index
  var_obs <- get_observed_variance(data = data, time = time, categories = categories)
  var_exp <- get_expected_variance_closed_form(data = data, time = time,
                                               categories = categories, m = m)
  compute_dispersion_index(var_obs = var_obs, var_exp = var_exp, categories = categories)
}

get_dispersion_index_multinomial <- function(data, time, categories, m, formula, newdata_fun) {
  # This function computes the proposed dispersion index with the expected
  # variance obtained from a multinomial model fitted at each time (see
  # get_expected_variance_multinomial).
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     categories: names of the category columns
  #     m: number of individuals in each group
  #     formula: multinomial model formula, e.g. cbind(Y.1, Y.2, Y.3) ~ X
  #     newdata_fun: function that receives the data of one time and returns
  #                  the prediction grid (e.g. newdata_observed_covariate)
  # Output:
  #     index: value of the dispersion index
  var_obs <- get_observed_variance(data = data, time = time, categories = categories)
  var_exp <- get_expected_variance_multinomial(data = data, time = time, formula = formula,
                                               newdata_fun = newdata_fun, m = m)
  compute_dispersion_index(var_obs = var_obs, var_exp = var_exp, categories = categories)
}

get_dispersion_index_mixed <- function(data, time, categories, m, formula, random,
                                       conditional = TRUE) {
  # This function computes the proposed dispersion index with the expected
  # variance obtained from a random effects multinomial model fitted to the
  # whole data set (see get_expected_variance_mixed).
  # Input:
  #     data: data frame in wide format (one column per category)
  #     time: name of the time column
  #     categories: names of the category columns
  #     m: number of individuals in each group
  #     formula: fixed effects formula, e.g. cbind(Y.1, Y.2, Y.3) ~ X
  #     random: random effects formula, e.g. ~ 1 | id
  #     conditional: if TRUE (default) the predicted probabilities include the
  #                  estimated random effects (see get_expected_variance_mixed)
  # Output:
  #     index: value of the dispersion index
  var_obs <- get_observed_variance(data = data, time = time, categories = categories)
  var_exp <- get_expected_variance_mixed(data = data, time = time, formula = formula,
                                         random = random, m = m, conditional = conditional)
  compute_dispersion_index(var_obs = var_obs, var_exp = var_exp, categories = categories)
}

## Simulation study -----
run_dispersion_index_replications <- function(n, j, k, m, theta, theta_label, n_simulations) {
  # This function simulates grouped longitudinal polytomous data n_simulations
  # times and computes the dispersion index (multinomial model expected
  # variance) of each replication.
  # Input:
  #     n: number of groups
  #     j: number of categories
  #     k: number of evaluated periods
  #     m: number of individuals in each group
  #     theta: standard deviation of the random effect
  #     theta_label: label of the scenario stored in the theta column
  #     n_simulations: number of replications
  # Output:
  #     results: data frame with the index of each replication and the
  #              scenario (n, j, k, m and theta as character columns)
  categories <- paste0("Y.", seq_len(j))
  model_formula <- make_multinomial_formula(j = j)
  index <- replicate(n_simulations, {
    data <- create_longitudinal_grouped_polytomous_data(n = n, j = j, k = k, theta = theta, m = m)
    get_dispersion_index_multinomial(data = data, time = "Time", categories = categories, m = m,
                                     formula = model_formula,
                                     newdata_fun = newdata_observed_covariate)
  })
  data.frame(index = index,
             n = as.character(n),
             j = as.character(j),
             k = as.character(k),
             m = as.character(m),
             theta = theta_label)
}

get_simulation_results <- function(n, j, k, m, n_simulations, theta_super) {
  # This function is the wrapper of the simulation study: it computes the
  # dispersion index for the scenario without overdispersion (theta = 0.01)
  # and for the scenario with overdispersion (theta = theta_super).
  # Input:
  #     n: number of groups
  #     j: number of categories
  #     k: number of evaluated periods
  #     m: number of individuals in each group
  #     n_simulations: number of replications of each scenario
  #     theta_super: standard deviation of the random effect in the
  #                  overdispersion scenario
  # Output:
  #     results: data frame with the overall results of the simulation study
  # NOTE: the label of the overdispersion scenario is kept as "10" (the value
  # used in the simulation study) regardless of theta_super, to preserve the
  # original behaviour.
  rbind(run_dispersion_index_replications(n = n, j = j, k = k, m = m, theta = 0.01,
                                          theta_label = "0.01", n_simulations = n_simulations),
        run_dispersion_index_replications(n = n, j = j, k = k, m = m, theta = theta_super,
                                          theta_label = "10", n_simulations = n_simulations))
}

get_index_comparison_results <- function(n, j, k, m, theta, n_simulations, random = ~ 1 | id) {
  # This function simulates grouped longitudinal polytomous data n_simulations
  # times for one scenario and computes, for each replication, the dispersion
  # index with three forms of the expected variance: the multinomial model
  # fitted at each time ("multinomial", get_dispersion_index_multinomial) and
  # the random effects multinomial model (get_dispersion_index_mixed) with the
  # probabilities conditional on the random effects ("mixed") and at the
  # population level ("mixed_population"). A model fit that fails gives NA.
  # Input:
  #     n: number of groups
  #     j: number of categories
  #     k: number of evaluated periods
  #     m: number of individuals in each group
  #     theta: standard deviation of the random effect
  #     n_simulations: number of replications
  #     random: random effects formula of the mixed model (default ~ 1 | id)
  # Output:
  #     results: data frame in long format with the columns n, j, k, m, theta,
  #              replication, index_form and index
  categories <- paste0("Y.", seq_len(j))
  model_formula <- make_multinomial_formula(j = j)
  compute_safely <- function(expression) {
    suppressWarnings(tryCatch(expression, error = function(e) NA_real_))
  }
  results <- lapply(seq_len(n_simulations), function(replication) {
    data <- create_longitudinal_grouped_polytomous_data(n = n, j = j, k = k, theta = theta, m = m)
    index <- c(
      multinomial = compute_safely(
        get_dispersion_index_multinomial(data = data, time = "Time", categories = categories,
                                         m = m, formula = model_formula,
                                         newdata_fun = newdata_observed_covariate)),
      mixed = compute_safely(
        get_dispersion_index_mixed(data = data, time = "Time", categories = categories, m = m,
                                   formula = model_formula, random = random,
                                   conditional = TRUE)),
      mixed_population = compute_safely(
        get_dispersion_index_mixed(data = data, time = "Time", categories = categories, m = m,
                                   formula = model_formula, random = random,
                                   conditional = FALSE))
    )
    data.frame(replication = replication, index_form = names(index), index = unname(index))
  })
  data.frame(n = n, j = j, k = k, m = m, theta = theta, bind_rows(results))
}

get_auc <- function(x, reference) {
  # This function computes the probability that a value of x exceeds a value
  # of reference (area under the ROC curve of the index as a discriminator of
  # the two scenarios; ties count one half): 0.5 = no discrimination,
  # 1 = perfect discrimination.
  # Input:
  #     x: index values of the scenario with dispersion
  #     reference: index values of the scenario without dispersion
  # Output:
  #     auc: value between 0 and 1
  mean(outer(x, reference, ">")) + 0.5 * mean(outer(x, reference, "=="))
}

get_discrimination_summary <- function(results, theta_reference = 0.01, false_alarm_rate = 0.05) {
  # This function summarises how well each form of the index discriminates the
  # scenarios with dispersion (theta > theta_reference) from the reference
  # scenario without dispersion (theta = theta_reference), separately for each
  # combination of n, j, k and m.
  # Input:
  #     results: output of get_index_comparison_results for several values of
  #              theta (bound by rows), including theta_reference
  #     theta_reference: theta of the scenario without dispersion
  #     false_alarm_rate: the detection threshold is the (1 - false_alarm_rate)
  #                       quantile of the index under theta_reference
  # Output:
  #     summary: data frame with one row per n, j, k, m, index_form and theta
  #              (theta_reference excluded) and the columns
  #              n_failed: replications whose model fit failed (NA index)
  #              mean_reference: mean index under theta_reference
  #              mean_index, sd_index: mean and sd of the index under theta
  #              auc: probability that the index under theta exceeds the index
  #                   under theta_reference (see get_auc)
  #              power: proportion of index values under theta above the
  #                     detection threshold (detection rate at false_alarm_rate)
  reference <- results %>%
    filter(theta == theta_reference, !is.na(index)) %>%
    group_by(n, j, k, m, index_form) %>%
    summarise(reference_index = list(index), .groups = "drop")
  results %>%
    filter(theta != theta_reference) %>%
    group_by(n, j, k, m, index_form, theta) %>%
    summarise(index_values = list(index), .groups = "drop") %>%
    inner_join(reference, by = c("n", "j", "k", "m", "index_form")) %>%
    rowwise() %>%
    mutate(n_failed = sum(is.na(index_values)),
           mean_reference = mean(reference_index),
           mean_index = mean(index_values, na.rm = TRUE),
           sd_index = sd(index_values, na.rm = TRUE),
           auc = get_auc(x = na.omit(index_values), reference = reference_index),
           power = mean(na.omit(index_values) >
                          quantile(reference_index, 1 - false_alarm_rate))) %>%
    ungroup() %>%
    select(-index_values, -reference_index)
}

## Case study helpers -----
get_probabilities <- function(etas) {
  # This function computes the category probabilities from the linear
  # predictors (etas) of the non-baseline categories of a multinomial logit
  # model. The probability of the baseline category is returned last.
  # Input:
  #     etas: vector with the linear predictors of the non-baseline categories
  # Output:
  #     probs: vector of probabilities (non-baseline categories, then baseline)
  den <- 1 + sum(exp(etas))
  unname(c(exp(etas) / den, 1 / den))
}

# Plot settings -----
theme_new <- function(base_size = 20, base_family = "Arial") {
  # ggplot2 theme shared by all the plots of the project.
  # Input:
  #     base_size: base font size
  #     base_family: base font family
  # Output:
  #     theme: ggplot2 theme object
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.text = element_text(size = 20, colour = "grey30"),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "grey20"),
          plot.title.position = "plot",
          legend.position = "bottom")
}

plot_simulation_results <- function(final_results, n_value) {
  # This function draws the boxplots of the index values of the simulation
  # study for a given number of groups n, by m, theta, j and k.
  # Input:
  #     final_results: data frame with the simulation results (columns index
  #                    and the factors n, j, k, m and theta)
  #     n_value: number of groups (n) to plot
  # Output:
  #     simulation_plot: ggplot object
  final_results %>%
    filter(n == n_value) %>%
    ggplot(mapping = aes(x = m, y = index, colour = theta)) +
    geom_boxplot() +
    facet_wrap(k ~ j) +
    theme_new() +
    ylab("Proposed index values") +
    scale_color_manual(values = c("#A3C4D9", "#043259")) +
    labs(colour = expression(theta))
}

index_form_labels <- c(multinomial = "Multinomial per time",
                       mixed = "Random effects (conditional)",
                       mixed_population = "Random effects (population)")

label_comparison_factors <- function(data) {
  # This function converts the factors of the comparison study into the
  # labelled factors used by its plots: theta as a factor, m, j and k labelled
  # and kept in numeric order, and the form of the expected variance labelled
  # as in index_form_labels.
  # Input:
  #     data: output of get_index_comparison_results or of
  #           get_discrimination_summary
  # Output:
  #     data: the same data with theta, m, j, k and index_form as factors
  data %>%
    mutate(theta = factor(theta),
           m = factor(paste("m =", m), levels = paste("m =", sort(unique(m)))),
           j = factor(paste("j =", j), levels = paste("j =", sort(unique(j)))),
           k = factor(paste("k =", k), levels = paste("k =", sort(unique(k)))),
           index_form = factor(index_form, levels = names(index_form_labels),
                               labels = index_form_labels))
}

plot_index_comparison <- function(results) {
  # This function draws the boxplots of the index values by theta for each
  # form of the expected variance, with one panel per group size m (rows) and
  # per combination of the number of categories j and of periods k (columns).
  # The y scale is free between rows because the index increases with m.
  # The size of the axis text is reduced with respect to theme_new() to suit
  # the number of panels.
  # Input:
  #     results: output of get_index_comparison_results (several thetas, m,
  #              j and k)
  # Output:
  #     comparison_plot: ggplot object
  results %>%
    filter(!is.na(index)) %>%
    label_comparison_factors() %>%
    ggplot(mapping = aes(x = theta, y = index, colour = index_form)) +
    geom_boxplot() +
    facet_grid(m ~ j + k, scales = "free_y") +
    theme_new() +
    xlab(expression(theta)) +
    ylab("Dispersion index") +
    scale_color_manual(name = "Expected variance",
                       values = c("#A3C4D9", "#043259", "#A34567")) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(axis.text = element_text(size = 12, colour = "grey30"),
          legend.title.position = "top")
}

plot_discrimination_summary <- function(discrimination_summary) {
  # This function draws the AUC and the power of each form of the expected
  # variance as a function of theta, with one panel per group size m (rows)
  # and per combination of the number of categories j and of periods k
  # (columns). The two metrics are distinguished by the line type and the
  # shape of the points. The size of the axis text is reduced with respect to
  # theme_new() to suit the number of panels.
  # Input:
  #     discrimination_summary: output of get_discrimination_summary
  # Output:
  #     discrimination_plot: ggplot object
  discrimination_summary %>%
    pivot_longer(cols = c(auc, power), names_to = "metric", values_to = "value") %>%
    label_comparison_factors() %>%
    mutate(metric = factor(metric, levels = c("auc", "power"),
                           labels = c("AUC", "Power"))) %>%
    ggplot(mapping = aes(x = theta, y = value, colour = index_form,
                         linetype = metric, shape = metric,
                         group = interaction(index_form, metric))) +
    geom_line(position = position_dodge(width = 0.3)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
    facet_grid(m ~ j + k) +
    theme_new() +
    xlab(expression(theta)) +
    ylab("Discrimination from no dispersion") +
    ylim(0, 1) +
    scale_color_manual(name = "Expected variance",
                       values = c("#A3C4D9", "#043259", "#A34567")) +
    scale_linetype_manual(name = "Metric", values = c(AUC = "solid", Power = "dashed")) +
    scale_shape_manual(name = "Metric", values = c(AUC = 16, Power = 17)) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
           linetype = guide_legend(nrow = 2, order = 2),
           shape = guide_legend(nrow = 2, order = 2)) +
    theme(axis.text = element_text(size = 12, colour = "grey30"),
          legend.title.position = "top")
}
