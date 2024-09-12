####################################################################################################
###
### File:    02_exploring_group_case.r
### Purpose: Simple example for simulating data for a 
###         scenario with individual
###         polytomous longitudinal data
### Authors: Gabriel Rodrigues Palma
### Date:    01/05/24
###
####################################################################################################
# Loading packages -----
source('00_source.r')

# Main functions used
create_prob_vector_cont_data_rf <- function(n, j, theta){
  # This function receives the number of categories 
  # and creates a vector of probability for each category
  if(j == 3){
    x <<- rnorm(n) # Covariate
    u <- rnorm(1, mean = 0, sd = theta)
    z2 <- 1+0.5**x + u
    z3 <- 0.5+x + u
    den <- 1+exp(z2)+exp(z3)
    p1 <- 1/den; p2<-exp(z2)/den; p3<-exp(z3)/den
    prob <- cbind(p1,p2,p3)
    
  }
  else if(j == 4){
    
    x <<- rnorm(n) # Covariate
    u <- rnorm(1, mean = 0, sd = theta)
    z2 <- 1+0.5**x + u
    z3 <- 0.5+x + u
    z4 <- 1.5 - x + u
    den <- 1+exp(z2)+exp(z3)+exp(z4)
    p1 <- 1/den; p2<-exp(z2)/den; p3<-exp(z3)/den; p4<-exp(z4)/den
    prob <- cbind(p1,p2,p3, p4)
    
  }
  else{
    x <<- rnorm(n) # Covariate
    u <- rnorm(1, mean = 0, sd = theta)
    z2 <- 1+0.5**x + u
    z3 <- 0.5+x + u
    z4 <- 1.5 - x + u
    z5 <- 1 - 0.7*x + u
    den <- 1+exp(z2)+exp(z3)+exp(z4)+exp(z5)
    p1 <- 1/den; p2<-exp(z2)/den; p3<-exp(z3)/den; p4<-exp(z4)/den;p5 <- exp(z5)/den
    prob <- cbind(p1,p2,p3, p4, p5)
  }
  result <- list()
  result$prob <- prob
  result$x <- x
  
  return(result)
}

obtain_longitudinal_data_per_individal_grouped <- function(id, j, k,
                                                           theta, m){
  # This function creates a longitudinal dataset
  #for a individual with k time steps
  # Input: 
  #       j = Number of categories
  #       k = Number of evaluated periods
  #       theta = standard deviation of the random effect     
  #       m = group of individuals
  
  periods <- as.matrix(1:k)
  longitudinal_polytomous_data_frame <- apply(periods, MARGIN = 1, FUN = function(x){
    
    cont_polytomous_data <- create_prob_vector_cont_data_rf(n = 1, 
                                                            j = j, 
                                                            theta = theta)
    x <- cont_polytomous_data$x
    prob <- cont_polytomous_data$prob
    
    y <- t(apply(prob, 1, rmultinom, n=1, size = m))
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

cont_polytomous_data <- create_prob_vector_cont_data_rf(n = 1, 
                                                        j = 4, 
                                                        theta = 0.01)

create_longitudinal_grouped_polytomous_data <- function(n, j, 
                                                        k, theta, m){
  # This function creates a dataset of longitudinal 
  #polytomous data.
  #Inputs:
  #       n = Number of individuals
  #       j = Number of categories
  #       k = Number of evaluated periods
  #       m = group of individuals
  longitudinal_data <- obtain_longitudinal_data_per_individal_grouped(id = 1, 
                                                                      j = j, 
                                                                      k = k,
                                                                      theta = theta, 
                                                                      m = m)
  for (id in 2:n){
    temp <- obtain_longitudinal_data_per_individal_grouped(id = id, 
                                                           j = j, 
                                                           k = k,
                                                           theta = theta, 
                                                           m = m)
    longitudinal_data <- rbind(longitudinal_data, temp)
  }
  return(longitudinal_data)
}

# Individual data
individual_longitudinal_data <- create_longitudinal_grouped_polytomous_data(n = 4, j = 3, 
                                                                            k = 3, theta = 10, m = 20)

####################################################################################################
################################### Simulating overdispersion ######################################
####################################################################################################
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
#get_dispersion_index(m = 16, data = Dados, j = 3)
###########################################################################################################################
#################################### Preparing the simulation study #######################################################
###########################################################################################################################

get_simulation_results <- function(n, j, k, m, n_simulations, theta_super) {
  # This function is a wrapper that contains the main functions related to the simulation study
  #Inputs:
  #       n = Number of individuals
  #       j = Number of categories
  #       k = Number of evaluated periods
  #       m = group of individuals
  #       n_simulations: The number of times the simulation is repeated
  # Output:
  #       results: A dataframe containing the overall results of the simulation study
  
  # Obtaining results for theta = 0.001
  dispersion_index_theta_0001 <- numeric(n_simulations)
  theta <- 0.01
  for (rep in 1:n_simulations){
    data <- create_longitudinal_grouped_polytomous_data(n = n, j = j, k = k,
                                                        theta = theta, m = m)
    dispersion_index_theta_0001[rep] <- get_dispersion_index(m = m, data = data, j = j)
  }
  
  dispersion_index_theta_0001 <- data.frame(index = dispersion_index_theta_0001, 
                                            n = rep(as.character(n), n_simulations), 
                                            j = rep(as.character(j), n_simulations), 
                                            k = rep(as.character(k), n_simulations), 
                                            m = rep(as.character(m), n_simulations), 
                                            theta = rep('0.01', n_simulations))
  
  # Obtaining results for theta = 10
  dispersion_index_theta_10 <- numeric(n_simulations)
  theta <- theta_super
  for (rep in 1:n_simulations){
    data <- create_longitudinal_grouped_polytomous_data(n = n, j = j, k = k,
                                                        theta = theta, m = m)
    dispersion_index_theta_10[rep] <- get_dispersion_index(m = m, data = data, j = j)
  }
  
  dispersion_index_theta_10 <- data.frame(index = dispersion_index_theta_10, 
                                          n = rep(as.character(n), n_simulations), 
                                          j = rep(as.character(j), n_simulations), 
                                          k = rep(as.character(k), n_simulations), 
                                          m = rep(as.character(m), n_simulations), 
                                          theta = rep('10', n_simulations))
  results <- rbind(dispersion_index_theta_0001, 
                   dispersion_index_theta_10)
}

# Creating multiple values of the index
options(dplyr.summarise.inform = FALSE)
m <- 10000
a <- get_simulation_results(n = 100, j = 3, k = 10, 
                            m = m, n_simulations = 10, theta_super = 100)
a$m <- rep(m, nrow(a))
a %>%
  group_by(theta) %>%
  summarise(min = min(index), 
            max = max(index), 
            m_min = min(index)/unique(m),
            m_max = max(index)/unique(m))

# Checking multiple values of m
m_behaviour_data <- data.frame(theta = NaN, m_min = NaN,
                               m_mean = NaN,m_max = NaN, m = NaN)
for (m in 1:50){
  cat('m = ', m, '\n')
  a <- get_simulation_results(n = 100, j = 5, k = 10, 
                              m = m, n_simulations = 100, theta_super = 10)
  a$m <- rep(m, nrow(a))
  result <- a %>%
            group_by(theta) %>%
            summarise(m_min = min(index)/unique(m),
                      m_mean = mean(index)/unique(m),
                      m_max = max(index)/unique(m))
  result$m <- rep(m, nrow(result))
  m_behaviour_data <- rbind(m_behaviour_data, result)
}
m_behaviour_data <- m_behaviour_data %>% drop_na()

m_behaviour_data %>%
  ggplot(aes(x = m, y = m_mean, colour = theta)) +
  geom_line() +
  geom_point(aes(x = m, y = m_mean)) +
  geom_ribbon(aes(x = m, ymin = m_min, ymax = m_max), 
              fill="#232323", alpha = 0.4) +
  #facet_wrap(~as.factor(Costumer), scales = 'free') +
  theme_new() +
  xlab(expression(italic(m))) +
  #geom_hline(yintercept=0.6, linetype="dashed", colour = 'red') +
  geom_vline(xintercept=2, linetype="dashed", colour = 'red') +
  ylab("Corrected longitudinal \n multinomial dispersion index") +
  scale_color_manual(name = "Clear overdispersion",
                     values=c("0.01" = "#232323", "10" = "#A34567"), 
                     label=c("0.01" = "No", "10" = "Yes")) 
  
ggsave('Plots/m_behaviour.png', dpi = 300, 
       height = 6, width = 8)  


