####################################################################################################
###
### File:    02_exploring_group_case.r
### Purpose: Exploring the behaviour of the proposed dispersion index as a
###          function of the group size m for grouped longitudinal
###          polytomous data.
### Authors: Gabriel Rodrigues Palma
### Date:    01/05/24
###
####################################################################################################
# Loading packages and functions -----
source("00_source.r")

# Quick check of the index range for a large group size -----
group_size <- 100
quick_check <- get_simulation_results(n = 100, j = 3, k = 10, m = group_size,
                                      n_simulations = 10, theta_super = 100)
quick_check %>%
  group_by(theta) %>%
  summarise(min = min(index),
            max = max(index),
            m_min = min(index) / group_size,
            m_max = max(index) / group_size,
            .groups = "drop")

# Behaviour of the index for m = 1, ..., 50 -----
m_behaviour_data <- lapply(1:10, function(group_size) {
  cat("m = ", group_size, "\n")
  get_simulation_results(n = 100, j = 5, k = 10, m = group_size,
                         n_simulations = 100, theta_super = 10) %>%
    group_by(theta) %>%
    summarise(m_min = min(index),
              m_mean = mean(index),
              m_max = max(index),
              .groups = "drop") %>%
    mutate(m = group_size)
}) %>%
  bind_rows()

m_behaviour_plot <- m_behaviour_data %>%
  ggplot(aes(x = m, y = m_mean, colour = theta)) +
  geom_line() +
  geom_point(aes(x = m, y = m_mean)) +
  geom_ribbon(aes(x = m, ymin = m_min, ymax = m_max),
              fill = "#232323", alpha = 0.4) +
  theme_new() +
  xlab(expression(italic(m))) +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "red") +
  ylab("Corrected longitudinal \n multinomial dispersion index") +
  scale_color_manual(name = "Clear overdispersion",
                     values = c("0.01" = "#232323", "10" = "#A34567"),
                     label = c("0.01" = "No", "10" = "Yes"))
m_behaviour_plot
ggsave("Plots/m_behaviour.png", plot = m_behaviour_plot, dpi = 300, height = 6, width = 8)
