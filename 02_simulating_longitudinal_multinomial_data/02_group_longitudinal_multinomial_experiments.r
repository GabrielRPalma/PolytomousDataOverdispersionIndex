####################################################################################################
###
### File:    02_group_longitudinal_multinomial_experiments.r
### Purpose: Simulation study of the proposed dispersion index for grouped
###          longitudinal polytomous data. The expected variance is obtained
###          from a multinomial model fitted at each evaluated period.
### Authors: Gabriel Rodrigues Palma
### Date:    07/09/23
###
####################################################################################################
# Loading packages and functions -----
source("00_source.r")

# Example: simulating grouped data and computing the index -----
example_data <- create_longitudinal_grouped_polytomous_data(n = 4, j = 3, k = 3,
                                                            theta = 10, m = 20)
example_index <- get_dispersion_index_multinomial(
  data = example_data,
  time = "Time",
  categories = paste0("Y.", 1:3),
  m = 20,
  formula = make_multinomial_formula(j = 3),
  newdata_fun = newdata_observed_covariate
)
example_index

# Simulation study -----
# Scenarios: j varies fastest, then k, m and n (this keeps the original row
# order of final_results)
simulation_scenarios <- expand.grid(j = c(3, 4, 5),
                                    k = c(3, 4),
                                    m = c(5, 10, 15),
                                    n = c(100, 200, 500))

simulation_results <- Map(function(n, j, k, m) {
  get_simulation_results(n = n, j = j, k = k, m = m,
                         n_simulations = 1000, theta_super = 10)
}, simulation_scenarios$n, simulation_scenarios$j, simulation_scenarios$k, simulation_scenarios$m)
names(simulation_results) <- with(simulation_scenarios,
                                  paste0("n", n, "_m", m, "_j", j, "_k", k))
# Individual scenarios can be inspected as, e.g., simulation_results$n100_m5_j3_k3

final_results <- bind_rows(simulation_results)
write.csv(final_results, "output_data/final_results.csv")

# Visualising the results -----
final_results <- read.csv("output_data/final_results.csv")[, -1] %>%
  mutate(n = factor(n),
         m = factor(m),
         theta = factor(theta),
         j = factor(paste("j =", j)),
         k = factor(paste("k =", k)))

for (n_value in c(100, 200, 500)) {
  simulation_plot <- plot_simulation_results(final_results, n_value = n_value)
  ggsave(paste0("Plots/simulation_results_n", n_value, ".png"), plot = simulation_plot,
         dpi = 100, height = 6, width = 9)
}
