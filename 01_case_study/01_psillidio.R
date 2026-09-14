####################################################################################################
###
### File:    01_psillidio.R
### Purpose: Case study for longitudinal multinomial data analysis (psyllid movement
###          between plants, input_data/psyllid_movement_full.csv). The expected
###          variance of the dispersion index is obtained from a multinomial
###          model fitted at each evaluation time.
### Authors: Gabriel Rodrigues Palma and ...
### Date:    20/08/26
###
####################################################################################################
# Loading packages and functions -----
source("00_source.r")

# Reading data -----
psyllid_data <- read.csv2("input_data/psyllid_movement_full.csv") %>%
  mutate(across(c(id, hour, density, plant), factor))
str(psyllid_data)

# Wide format: one column per plant
psyllid_data_wide <- psyllid_data %>%
  pivot_wider(names_from = plant, values_from = y, values_fill = 0)
str(psyllid_data_wide)
head(psyllid_data_wide)

# Descriptive analysis -----
# Mean and variance per treatment (density), time (hour) and plant
summary_by_density_hour <- psyllid_data %>%
  group_by(density, hour, plant) %>%
  summarise(media = mean(y),
            variancia = var(y),
            .groups = "drop")
summary_by_density_hour

# Mean proportion of insects per plant (90 insects per cage)
psyllid_data %>%
  group_by(hour, density, plant) %>%
  summarise(prop = mean(y / 90), .groups = "drop") %>%
  ggplot(aes(x = hour, y = prop, color = plant, group = plant)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.2) +
  facet_wrap(. ~ density, ncol = 3) +
  labs(x = "Tempo (horas)",
       y = "Proporção média de insetos",
       color = "Tratamento") +
  theme_new(base_size = 14) +
  theme(legend.position = "bottom")

# Dispersion index: three forms of the expected variance -----
m <- 90 # Number of females per cage
categories <- c("A", "B", "C", "D", "none")

# 1. Closed form: mean(y) * (m - mean(y)) / m at each hour
get_dispersion_index_closed_form(data = psyllid_data_wide,
                                 time = "hour",
                                 categories = categories,
                                 m = m)

# 2. Multinomial model fitted at each hour (nnet::multinom), predicted at the
#    levels of the density treatment
get_dispersion_index_multinomial(
  data = psyllid_data_wide,
  time = "hour",
  categories = categories,
  m = m,
  formula = cbind(A, B, C, D, none) ~ density,
  newdata_fun = function(data) expand.grid(density = levels(data$density))
)
# 3. Random effects multinomial model fitted to the whole data set
#    (mclogit::mblogit, random intercept per cage), probabilities of each
#    observation conditional on the estimated random effects
get_dispersion_index_mixed(
  data = psyllid_data_wide,
  time = "hour",
  categories = categories,
  m = m,
  formula = cbind(A, B, C, D, none) ~ density,
  random = ~ 1 | id
)
