####################################################################################################
###
### File:    02_comparing_multinomial_and_mixed_indices.r
### Purpose: Simulation study comparing the dispersion index computed with the
###          multinomial model fitted at each time and with the random effects
###          multinomial model, focusing on the discrimination between the
###          scenarios with dispersion (theta = 0.1, 1, 10) and the scenario
###          without dispersion (theta = 0.01), for several numbers of
###          categories j, of evaluated periods k and group sizes m.
### Authors: Gabriel Rodrigues Palma
### Date:    20/08/26
###
####################################################################################################
# Loading packages and functions -----
source("00_source.r")

# Simulation settings -----
# The number of groups n is fixed and the compared factors are the standard
# deviation of the random effect (theta = 0.01 is the reference scenario
# without dispersion), the number of categories j, the number of evaluated
# periods k and the group size m.
n <- 100
j_values <- c(3, 5)
k_values <- c(3, 6)
m_values <- c(5, 15, 30)
theta_values <- c(0.01, 0.1, 1, 10)
n_simulations <- 100
# NOTE: in the simulated data a new random effect is drawn at every period, so
# the random intercept per individual (random = ~ 1 | id, the structure used
# in the case studies) does not reproduce the data generating process exactly.

# Simulation -----
# Each replication computes the index with the multinomial model per time and
# with the random effects model (conditional and population probabilities)
scenarios <- expand.grid(theta = theta_values, m = m_values, j = j_values, k = k_values)
comparison_results <- Map(function(theta, m, j, k) {
  cat("j = ", j, " k = ", k, " m = ", m, " theta = ", theta, "\n")
  get_index_comparison_results(n = n, j = j, k = k, m = m, theta = theta,
                               n_simulations = n_simulations, random = ~ 1 | id)
}, scenarios$theta, scenarios$m, scenarios$j, scenarios$k) %>%
  bind_rows()
write.csv(comparison_results, "output_data/index_comparison_results.csv", row.names = FALSE)

# Discrimination between dispersion and no dispersion -----
comparison_results <- read.csv("output_data/index_comparison_results.csv")
discrimination_summary <- get_discrimination_summary(comparison_results,
                                                     theta_reference = 0.01,
                                                     false_alarm_rate = 0.05)
as.data.frame(discrimination_summary)
write.csv(discrimination_summary, "output_data/index_comparison_summary.csv", row.names = FALSE)

# Visualising the results -----
comparison_plot <- plot_index_comparison(comparison_results)
comparison_plot
ggsave("Plots/index_comparison.png", plot = comparison_plot, dpi = 300, height = 9, width = 12)

discrimination_plot <- plot_discrimination_summary(discrimination_summary)
discrimination_plot
ggsave("Plots/index_comparison_discrimination.png", plot = discrimination_plot,
       dpi = 300, height = 9, width = 12)
