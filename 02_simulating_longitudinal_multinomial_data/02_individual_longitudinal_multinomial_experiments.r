####################################################################################################
###
### File:    02_individual_longitudinal_multinomial_experiments.r
### Purpose: Simple example for simulating individual (m = 1) polytomous
###          longitudinal data and comparing multinomial models with and
###          without random effects.
### Authors: Gabriel Rodrigues Palma
### Date:    07/09/23
###
####################################################################################################
# Loading packages and functions -----
source("00_source.r")

# Simulating individual data -----
individual_longitudinal_data <- create_longitudinal_polytomous_data(n = 4, j = 3,
                                                                    k = 3, theta = 0)

# Fitting multinomial models (nnet) -----
y <- individual_longitudinal_data %>%
  dplyr::select(Y.1, Y.2, Y.3) %>%
  as.matrix()
null_model_fit <- multinom(y ~ 1, data = individual_longitudinal_data, trace = FALSE)
x_model_fit <- multinom(y ~ X + Time, data = individual_longitudinal_data, trace = FALSE)
correct_model_fit <- multinom(y ~ X + 1 | Time, data = individual_longitudinal_data, trace = FALSE)
anova(x_model_fit, correct_model_fit)

# Fitting multinomial models with random effects (mclogit) -----
individual_longitudinal_data <- create_longitudinal_polytomous_data(n = 4, j = 3,
                                                                    k = 10, theta = 1)
y <- individual_longitudinal_data %>%
  dplyr::select(Y.1, Y.2, Y.3) %>%
  as.matrix()
m1 <- mblogit(y ~ 1, data = individual_longitudinal_data)
m2 <- mblogit(y ~ X, data = individual_longitudinal_data)
m3 <- mblogit(y ~ X, random = ~ 1 | id, data = individual_longitudinal_data)
anova(m1, m2, test = "Chisq")
anova(m1, m3, test = "Chisq")
anova(m2, m3, test = "Chisq")
