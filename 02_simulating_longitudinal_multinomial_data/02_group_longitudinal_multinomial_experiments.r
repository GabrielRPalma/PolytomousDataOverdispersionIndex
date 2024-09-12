####################################################################################################
###
### File:    02_group_longitudinal_multinomial_experiments.r
### Purpose: Simple example for simulating data for a 
###         scenario with individual
###         polytomous longitudinal data
### Authors: Gabriel Rodrigues Palma
### Date:    07/09/23
###
####################################################################################################
# Loading packages -----
source('00_source.r')

# Creating dataset -----
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
options(dplyr.summarise.inform = FALSE)

################################################################################################
###################################### Simulation n = 100 ######################################
################################################################################################
n100_m5_j3_k3 <- get_simulation_results(n = 100, j = 3, k = 3, m = 5, n_simulations = 1000)
sink('output_data/n100_m5_j3_k3.csv')
write.csv(n100_m5_j3_k3)
sink()
n100_m5_j4_k3 <- get_simulation_results(n = 100, j = 4, k = 3, m = 5, n_simulations = 1000)
n100_m5_j5_k3 <- get_simulation_results(n = 100, j = 5, k = 3, m = 5, n_simulations = 1000)
n100_m5_j3_k4 <- get_simulation_results(n = 100, j = 3, k = 4, m = 5, n_simulations = 1000)
n100_m5_j4_k4 <- get_simulation_results(n = 100, j = 4, k = 4, m = 5, n_simulations = 1000)
n100_m5_j5_k4 <- get_simulation_results(n = 100, j = 5, k = 4, m = 5, n_simulations = 1000)
n100_m10_j3_k3 <- get_simulation_results(n = 100, j = 3, k = 3, m = 10, n_simulations = 1000)
n100_m10_j4_k3 <- get_simulation_results(n = 100, j = 4, k = 3, m = 10, n_simulations = 1000)
n100_m10_j5_k3 <- get_simulation_results(n = 100, j = 5, k = 3, m = 10, n_simulations = 1000)
n100_m10_j3_k4 <- get_simulation_results(n = 100, j = 3, k = 4, m = 10, n_simulations = 1000)
n100_m10_j4_k4 <- get_simulation_results(n = 100, j = 4, k = 4, m = 10, n_simulations = 1000)
n100_m10_j5_k4 <- get_simulation_results(n = 100, j = 5, k = 4, m = 10, n_simulations = 1000)
n100_m15_j3_k3 <- get_simulation_results(n = 100, j = 3, k = 3, m = 15, n_simulations = 1000)
n100_m15_j4_k3 <- get_simulation_results(n = 100, j = 4, k = 3, m = 15, n_simulations = 1000)
n100_m15_j5_k3 <- get_simulation_results(n = 100, j = 5, k = 3, m = 15, n_simulations = 1000)
n100_m15_j3_k4 <- get_simulation_results(n = 100, j = 3, k = 4, m = 15, n_simulations = 1000)
n100_m15_j4_k4 <- get_simulation_results(n = 100, j = 4, k = 4, m = 15, n_simulations = 1000)
n100_m15_j5_k4 <- get_simulation_results(n = 100, j = 5, k = 4, m = 15, n_simulations = 1000)

################################################################################################
###################################### Simulation n = 200 ######################################
################################################################################################
n200_m5_j3_k3 <- get_simulation_results(n = 200, j = 3, k = 3, m = 5, n_simulations = 1000)
n200_m5_j4_k3 <- get_simulation_results(n = 200, j = 4, k = 3, m = 5, n_simulations = 1000)
n200_m5_j5_k3 <- get_simulation_results(n = 200, j = 5, k = 3, m = 5, n_simulations = 1000)
n200_m5_j3_k4 <- get_simulation_results(n = 200, j = 3, k = 4, m = 5, n_simulations = 1000)
n200_m5_j4_k4 <- get_simulation_results(n = 200, j = 4, k = 4, m = 5, n_simulations = 1000)
n200_m5_j5_k4 <- get_simulation_results(n = 200, j = 5, k = 4, m = 5, n_simulations = 1000)
n200_m10_j3_k3 <- get_simulation_results(n = 200, j = 3, k = 3, m = 10, n_simulations = 1000)
n200_m10_j4_k3 <- get_simulation_results(n = 200, j = 4, k = 3, m = 10, n_simulations = 1000)
n200_m10_j5_k3 <- get_simulation_results(n = 200, j = 5, k = 3, m = 10, n_simulations = 1000)
n200_m10_j3_k4 <- get_simulation_results(n = 200, j = 3, k = 4, m = 10, n_simulations = 1000)
n200_m10_j4_k4 <- get_simulation_results(n = 200, j = 4, k = 4, m = 10, n_simulations = 1000)
n200_m10_j5_k4 <- get_simulation_results(n = 200, j = 5, k = 4, m = 10, n_simulations = 1000)
n200_m15_j3_k3 <- get_simulation_results(n = 200, j = 3, k = 3, m = 15, n_simulations = 1000)
n200_m15_j4_k3 <- get_simulation_results(n = 200, j = 4, k = 3, m = 15, n_simulations = 1000)
n200_m15_j5_k3 <- get_simulation_results(n = 200, j = 5, k = 3, m = 15, n_simulations = 1000)
n200_m15_j3_k4 <- get_simulation_results(n = 200, j = 3, k = 4, m = 15, n_simulations = 1000)
n200_m15_j4_k4 <- get_simulation_results(n = 200, j = 4, k = 4, m = 15, n_simulations = 1000)
n200_m15_j5_k4 <- get_simulation_results(n = 200, j = 5, k = 4, m = 15, n_simulations = 1000)
  
################################################################################################
###################################### Simulation n = 500 ######################################
################################################################################################
n500_m5_j3_k3 <- get_simulation_results(n = 500, j = 3, k = 3, m = 5, n_simulations = 1000)
n500_m5_j4_k3 <- get_simulation_results(n = 500, j = 4, k = 3, m = 5, n_simulations = 1000)
n500_m5_j5_k3 <- get_simulation_results(n = 500, j = 5, k = 3, m = 5, n_simulations = 1000)
n500_m5_j3_k4 <- get_simulation_results(n = 500, j = 3, k = 4, m = 5, n_simulations = 1000)
n500_m5_j4_k4 <- get_simulation_results(n = 500, j = 4, k = 4, m = 5, n_simulations = 1000)
n500_m5_j5_k4 <- get_simulation_results(n = 500, j = 5, k = 4, m = 5, n_simulations = 1000)
n500_m10_j3_k3 <- get_simulation_results(n = 500, j = 3, k = 3, m = 10, n_simulations = 1000)
n500_m10_j4_k3 <- get_simulation_results(n = 500, j = 4, k = 3, m = 10, n_simulations = 1000)
n500_m10_j5_k3 <- get_simulation_results(n = 500, j = 5, k = 3, m = 10, n_simulations = 1000)
n500_m10_j3_k4 <- get_simulation_results(n = 500, j = 3, k = 4, m = 10, n_simulations = 1000)
n500_m10_j4_k4 <- get_simulation_results(n = 500, j = 4, k = 4, m = 10, n_simulations = 1000)
n500_m10_j5_k4 <- get_simulation_results(n = 500, j = 5, k = 4, m = 10, n_simulations = 1000)
n500_m15_j3_k3 <- get_simulation_results(n = 500, j = 3, k = 3, m = 15, n_simulations = 1000)
n500_m15_j4_k3 <- get_simulation_results(n = 500, j = 4, k = 3, m = 15, n_simulations = 1000)
n500_m15_j5_k3 <- get_simulation_results(n = 500, j = 5, k = 3, m = 15, n_simulations = 1000)
n500_m15_j3_k4 <- get_simulation_results(n = 500, j = 3, k = 4, m = 15, n_simulations = 1000)
n500_m15_j4_k4 <- get_simulation_results(n = 500, j = 4, k = 4, m = 15, n_simulations = 1000)
n500_m15_j5_k4 <- get_simulation_results(n = 500, j = 5, k = 4, m = 15, n_simulations = 1000)

################################################################################################
###################################### Visualising the results #################################
################################################################################################
final_results <- rbind(n100_m5_j3_k3, n100_m5_j4_k3, n100_m5_j5_k3, n100_m5_j3_k4, n100_m5_j4_k4, n100_m5_j5_k4,
                       n100_m10_j3_k3, n100_m10_j4_k3, n100_m10_j5_k3, n100_m10_j3_k4, n100_m10_j4_k4, n100_m10_j5_k4,
                       n100_m15_j3_k3, n100_m15_j4_k3, n100_m15_j5_k3, n100_m15_j3_k4, n100_m15_j4_k4, n100_m15_j5_k4,
                       n200_m5_j3_k3, n200_m5_j4_k3, n200_m5_j5_k3, n200_m5_j3_k4, n200_m5_j4_k4, n200_m5_j5_k4,
                       n200_m10_j3_k3, n200_m10_j4_k3, n200_m10_j5_k3, n200_m10_j3_k4, n200_m10_j4_k4, n200_m10_j5_k4,
                       n200_m15_j3_k3, n200_m15_j4_k3, n200_m15_j5_k3, n200_m15_j3_k4, n200_m15_j4_k4, n200_m15_j5_k4,
                       n500_m5_j3_k3, n500_m5_j4_k3, n500_m5_j5_k3, n500_m5_j3_k4, n500_m5_j4_k4, n500_m5_j5_k4,
                       n500_m10_j3_k3, n500_m10_j4_k3, n500_m10_j5_k3, n500_m10_j3_k4, n500_m10_j4_k4, n500_m10_j5_k4,
                       n500_m15_j3_k3, n500_m15_j4_k3, n500_m15_j5_k3, n500_m15_j3_k4, n500_m15_j4_k4, n500_m15_j5_k4)
sink('output_data/final_results.csv')
write.csv(final_results)
sink()

final_results <- read.csv('output_data/final_results.csv')[, -1]
final_results$n <- factor(final_results$n)
final_results$j <- factor(final_results$j)
levels(final_results$j) <- c('j = 3', 'j = 4', 'j = 5')
final_results$k <- factor(final_results$k)
levels(final_results$k) <- c('k = 3', 'k = 4', 'k = 5')
final_results$m <- factor(final_results$m)
final_results$theta <- factor(final_results$theta)


final_results %>%
  filter(n == 100) %>%
  ggplot(mapping = aes(x = m, y = index, colour = theta)) +
  geom_boxplot() +
  facet_wrap(k~j) +
  theme_new() +
  ylab('Proposed index values') +
  scale_color_manual(values = c("#A3C4D9", "#043259")) +
  labs(colour = expression(theta))
ggsave('Plots/simulation_results_n100.png', dpi = 100, height = 6, width = 9)

final_results %>%
  filter(n == 200) %>%
  ggplot(mapping = aes(x = m, y = index, colour = theta)) +
  geom_boxplot() +
  facet_wrap(k~j) +
  theme_new() +
  ylab('Proposed index values') +
  scale_color_manual(values = c("#A3C4D9", "#043259"))+
  labs(colour = expression(theta))
ggsave('Plots/simulation_results_n200.png', dpi = 100, height = 6, width = 9)

final_results %>%
  filter(n == 500) %>%
  ggplot(mapping = aes(x = m, y = index, colour = theta)) +
  geom_boxplot() +
  facet_wrap(k~j) +
  theme_new() +
  ylab('Proposed index values') +
  scale_color_manual(values = c("#A3C4D9", "#043259"))+
  labs(colour = expression(theta))
ggsave('Plots/simulation_results_n500.png', dpi = 100, height = 6, width = 9)

