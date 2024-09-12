####################################################################################################
###
### File:    01_case_study.R
### Purpose: Case study for longitudinal multinomial data analysis
### Authors: Maria Leticia and Gabriel Rodrigues Palma and ...
### Date:    01/11/23
###
####################################################################################################
# source packages ----
source('00_source.r')
library(ggplot2)
library(magrittr)
library(mclogit)
library(tidyverse)

compute_disersion_index <- function(){
  #This function will compute the dispersion index among longitudinal multivariate data
  #Input:
  #     observed_varience: A vector with means of the observed varience per individual
  #                       computed per category. 
  #     expected_variance: The variance expected by the multinomial model with random effects
  #Output:
  #     index: Value of the index
}
# Reading files ----
Dados <- read.csv2("input_data/Dados2.csv") %>% 
  mutate(baia = factor(baia),
         trat = factor(trat),
         dia = factor(dia),
         hora = factor(hora))
str(Dados)

# Passa do formato largo (wide) pro formato longo (long)
dados_long <- Dados %>% 
  gather(key = "categoria", value = "valor",
         -baia, -trat, -dia, -hora)
str(dados_long)


# ANÁLISE DESCRITIVA ------------------------------------------------------
#Media - tempo
resumo <- dados_long %>% 
  group_by(baia,trat, categoria) %>% 
  summarise(amp= max(valor)-min(valor),
            me = mean(valor), 
            va = var(valor)
  ) %>% 
  ungroup()

#Media - Baia
resumo <- dados_long %>% 
  group_by(dia, trat, categoria) %>% 
  summarise(amp= max(valor)-min(valor),
            me = mean(valor), 
            va = var(valor)
  ) %>% 
  ungroup()


View(resumo)


# MODELOS -----------------------------------------------------------------
mod1 <- mblogit(cbind(repousar, comer, explorar) ~ 1, 
                data = Dados, random = ~1|baia)
mod1$info.coef
predict(mod1, type = 'response')
summary(mod1)

mod2 <- mblogit(cbind(repousar, comer, explorar) ~ trat, 
                data = Dados, random = ~1|baia)
predict(mod2, type = 'link')


mod3 <- mblogit(cbind(repousar, comer, explorar) ~ trat, 
                data = Dados, random = ~1|dia)
predict(mod3, type = 'link')

1-exp(eta)/1 + sum()

eta1 <- -1.512194 
eta2 <- -0.7599371

get_probabilities <- function(x){
  # This function computes the probabilities for a given value of etas for the case study
  eta1 <- x[1]
  eta2 <- x[2]
  probs1 <- exp(eta1) / (1 + sum(exp(eta1), exp(eta2)))
  probs2 <- exp(eta2) / (1 + sum(exp(eta1), exp(eta2)))
  probs3 <- 1/(1 + sum(exp(eta1), exp(eta2)))
  return(c(probs1, probs2, probs3))
  
}

etas <- predict(mod3, type = 'link')
probs <- t(apply(as.matrix(etas), MARGIN = 1, FUN = get_probabilities))
probs_rep1
m <- 16 # Computed based on the number of times that the response varaibles was computed

var1 <- probs_rep1[1] * (1-probs_rep1[1]) * m
var2 <- probs_rep1[2] * (1-probs_rep1[2]) * m
var3 <- probs_rep1[3] * (1-probs_rep1[3]) * m

get_variances <- function(probs, m){
  # This function computes the expected variance based on a grouped polytomous data
  #Input: 
  #     probs: a vector of probabilities for each category, J
  #     m: number of individuals per group
  #Output:
  #     variances: a vector of variances per category, J
  variances <- apply(probs, MARGIN = 1, FUN = function(pi){
                    var1 <- pi[1] * (1-pi[1]) * m
                    var2 <- pi[2] * (1-pi[2]) * m
                    var3 <- pi[3] * (1-pi[3]) * m
                    return(c(var1, var2, var3))
                    
                  })
  return(t(variances))
  
}
get_variances(probs = probs, m = 16)

#######################################################################################################
######################################## Main code for the dispersion #################################
#######################################################################################################
m <- 16
var_obs <- Dados %>% 
           pivot_longer(cols = 5:7) %>%
           group_by(dia, name) %>%
           summarise(var = var(value))%>%
             pivot_wider(names_from = name, values_from = var)

var_exp <- Dados %>% 
  pivot_longer(cols = 5:7) %>%
  group_by(dia, name) %>%
  summarise(var = mean(value)*(m - mean(value))/m) %>%
    pivot_wider(names_from = name, values_from = var)
mean(apply(as.matrix(var_obs[2:4])/as.matrix(var_exp[2:4]), MARGIN = 2, mean))

get_dispersion_index <- function(m, data, j){
  # This function computes the proposed dispersion index
  
  var_obs <- data %>% 
    pivot_longer(cols = 3:(3+j-1)) %>%
    group_by(Time, name) %>%
    summarise(var = var(value))%>%
    pivot_wider(names_from = name, values_from = var)
  
  var_exp <- data %>% 
    pivot_longer(cols = 3:(3+j-1)) %>%
    group_by(Time, name) %>%
    summarise(var = mean(value)*(m - mean(value))/m) %>%
    pivot_wider(names_from = name, values_from = var)
  mean(apply(as.matrix(var_obs[2:4])/as.matrix(var_exp[2:4]), MARGIN = 2, mean))
  
}
get_dispersion_index(m = 16, data = Dados, j = 3)
# Simulate grouped data with and without dispersion
