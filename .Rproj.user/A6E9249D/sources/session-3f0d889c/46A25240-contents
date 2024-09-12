####################################################################################################
###
### File:    05_individual_longitudinal_multinomial_experiments.r
### Purpose: Simple example for simulating data for a 
###         scenario with individual
###         polytomous longitudinal data
### Authors: Gabriel Rodrigues Palma
### Date:    07/09/23
###
####################################################################################################
# Loading packages -----
source('00_source.r')
library(mclogit)
# Creating dataset -----
create_prob_vector_cont_data_rf <- function(n, j, theta){
  # This function receives the number of categories 
  # and creates a vector of probability for each category
  if(j == 3){
    x <<- rnorm(n) # Covariate
    z2 <- 1.5-3*x + rnorm(1, mean = 0, sd = theta)
    z3 <- 3-5*x 
    den <- 1+exp(z2)+exp(z3)
    p1 <- 1/den; p2<-exp(z2)/den; p3<-exp(z3)/den
    prob <- cbind(p1,p2,p3)
    
  }
  else if(j == 4){
    
    x <<- rnorm(n) # Covariate
    z2 <- 1.5-3*x + rnorm(1, mean = 0, sd = theta)
    z3 <- 3-5*x
    z4 <- 2 - 4*x
    den <- 1+exp(z2)+exp(z3)+exp(z4)
    p1 <- 1/den; p2<-exp(z2)/den; p3<-exp(z3)/den; p4<-exp(z4)/den
    prob <- cbind(p1,p2,p3, p4)
    
  }
  else{
    x <<- rnorm(n) # Covariate
    z2 <- 1.5-3*x + rnorm(1, mean = 0, sd = theta)
    z3 <- 3-5*x
    z4 <- 2 - 4*x
    z5 <- 4 - 7*x
    den <- 1+exp(z2)+exp(z3)+exp(z4)+exp(z5)
    p1 <- 1/den; p2<-exp(z2)/den; p3<-exp(z3)/den; p4<-exp(z4)/den;p5 <- exp(z5)/den
    prob <- cbind(p1,p2,p3, p4, p5)
  }
  result <- list()
  result$prob <- prob
  result$x <- x
  
  return(result)
}

obtain_longitudinal_data_per_individual <- function(id, j, k,
                                                    theta){
  # This function creates a longitudinal dataset
  #for a individual with k time steps
  # Input: 
  #       j = Number of categories
  #       k = Number of evaluated periods
  #       sigma = standard deviation of the random effect     
  
  periods <- as.matrix(1:k)
  longitudinal_polytomous_data_frame <- apply(periods, MARGIN = 1, FUN = function(x){
    
    cont_polytomous_data <- create_prob_vector_cont_data_rf(n = 1, 
                                                            j = j, 
                                                            theta = theta)
    x <- cont_polytomous_data$x
    prob <- cont_polytomous_data$prob
    
    y <- t(apply(prob, 1, rmultinom, n=1, size = 1))
    #y_index <-factor(apply(y, 1, function(x) which(x==1))) 
    
    longitudinal_polytomous_data_frame <- data.frame(id = id,
                                                     X = x,
                                                     Y = y)
    return(longitudinal_polytomous_data_frame)
    
  })
  
  longitudinal_polytomous_data_frame <- do.call(rbind, longitudinal_polytomous_data_frame)
  longitudinal_polytomous_data_frame$Time <- 1:k
  return(longitudinal_polytomous_data_frame)
}

cont_polytomous_data <- create_prob_vector_cont_data_rf(n = 50, 
                                                     j = 4, 
                                                     theta = 0.1)

create_longitudinal_polytomous_data <- function(n, j, 
                                                k, theta){
  # This function creates a dataset of longitudinal 
  #polytomous data.
  #Inputs:
  #       n = Number of individuals
  #       j = Number of categories
  #       k = Number of evaluated periods
  longitudinal_data <- obtain_longitudinal_data_per_individual(id = 1, 
                                                               j, k,
                                                               theta)
  for (id in 2:n){
    temp <- obtain_longitudinal_data_per_individual(id = id, 
                                                    j, k,
                                                   theta)
    longitudinal_data <- rbind(longitudinal_data, temp)
  }
  return(longitudinal_data)
  }

# Individual data
individual_longitudinal_data <- create_longitudinal_polytomous_data(n = 4, j = 3, 
                                                                    k = 3, theta = 0)

# Fittings null and correct models
y <- individual_longitudinal_data %>%
  dplyr::select(Y.1, Y.2, Y.3)
y <- as.matrix(y)
null_model_fit <- multinom(y ~ 1,data = individual_longitudinal_data,
                           trace="FALSE") 
x_model_fit <- multinom(y ~ X + Time,
                        data = individual_longitudinal_data,
                        trace="FALSE") 

correct_model_fit <- multinom(y ~ X + 1|Time,
                              data = individual_longitudinal_data,
                              trace="FALSE") 

anova(x_model_fit, correct_model_fit)

# Fitting the models using mclogit -----
individual_longitudinal_data <- create_longitudinal_polytomous_data(n = 4, j = 3, 
                                                                    k = 10, theta = 1)

# Fittings null and correct models
y <- individual_longitudinal_data %>%
  dplyr::select(Y.1, Y.2, Y.3)
y <- as.matrix(y)


m1 <- mblogit(formula = y~1,
              data = individual_longitudinal_data)
m2 <- mblogit(y~X,
              data = individual_longitudinal_data)
m3 <- mblogit(y~X, random = ~1|id,
              data = individual_longitudinal_data)

anova(m1, m2, test = 'Chisq')
anova(m1, m3, test = 'Chisq')
anova(m2, m3, test = 'Chisq')

