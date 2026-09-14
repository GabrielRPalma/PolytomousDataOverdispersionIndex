####################################################################################################
###
### File:    01_case_study.R
### Purpose: Case study for longitudinal multinomial data analysis (pig behaviour,
###          input_data/Dados2.csv). The expected variance of the dispersion
###          index is obtained with the closed form expression.
### Authors: Maria Leticia and Gabriel Rodrigues Palma and ...
### Date:    01/11/23
###
####################################################################################################
# Loading packages and functions -----
source("00_source.r")

# Reading data -----
pig_data <- read.csv2("input_data/Dados2.csv") %>%
  mutate(across(c(baia, trat, dia, hora), factor))
str(pig_data)

# Long format: one row per pen (baia), treatment, day, hour and behaviour category
pig_data_long <- pig_data %>%
  pivot_longer(cols = c(repousar, comer, explorar),
               names_to = "categoria", values_to = "valor")
str(pig_data_long)

# Descriptive analysis -----
# Range, mean and variance per pen (baia), treatment and category
summary_by_pen <- pig_data_long %>%
  group_by(baia, trat, categoria) %>%
  summarise(amp = max(valor) - min(valor),
            me = mean(valor),
            va = var(valor),
            .groups = "drop")

# Range, mean and variance per day (dia), treatment and category
summary_by_day <- pig_data_long %>%
  group_by(dia, trat, categoria) %>%
  summarise(amp = max(valor) - min(valor),
            me = mean(valor),
            va = var(valor),
            .groups = "drop")
if (interactive()) View(summary_by_day)

# Multinomial models with random effects (mclogit) -----
mod1 <- mblogit(cbind(repousar, comer, explorar) ~ 1,
                data = pig_data, random = ~ 1 | baia)
mod1$info.coef
predict(mod1, type = "response")
summary(mod1)

mod2 <- mblogit(cbind(repousar, comer, explorar) ~ trat,
                data = pig_data, random = ~ 1 | baia)
predict(mod2, type = "link")

mod3 <- mblogit(cbind(repousar, comer, explorar) ~ trat,
                data = pig_data, random = ~ 1 | dia)
predict(mod3, type = "link")

# Expected variances from the fitted probabilities -----
m <- 16 # Number of times the response variables were recorded per group
etas <- predict(mod3, type = "link")
probs <- t(apply(as.matrix(etas), MARGIN = 1, FUN = get_probabilities))
get_expected_variance_from_probabilities(probabilities = probs, m = m)

# Dispersion index: three forms of the expected variance -----
categories <- c("repousar", "comer", "explorar")

# 1. Closed form: mean(y) * (m - mean(y)) / m at each day
get_dispersion_index_closed_form(data = pig_data,
                                 time = "dia",
                                 categories = categories,
                                 m = m)

# 2. Multinomial model fitted at each day (nnet::multinom), predicted at the
#    levels of the treatment
get_dispersion_index_multinomial(
  data = pig_data,
  time = "dia",
  categories = categories,
  m = m,
  formula = cbind(repousar, comer, explorar) ~ trat,
  newdata_fun = function(data) expand.grid(trat = levels(data$trat))
)

# 3. Random effects multinomial model fitted to the whole data set
#    (mclogit::mblogit, random intercept per pen, as mod1), probabilities of
#    each observation conditional on the estimated random effects
get_dispersion_index_mixed(
  data = pig_data,
  time = "dia",
  categories = categories,
  m = m,
  formula = cbind(repousar, comer, explorar) ~ 1,
  random = ~ 1 | baia
)