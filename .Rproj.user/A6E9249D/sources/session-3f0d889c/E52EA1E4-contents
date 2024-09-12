####################################################################################################
###
### File:    00_source.R
### Purpose: Load required packages and functions used to 
###          create interactive maps.
### Authors: Gabriel Rodrigues Palma
### Date:    19/03/23
###
####################################################################################################

# packages required ---------------------------

packages <- c('nnet', #Fit Multinomial Log-linear Models
              'hnp', #half-normal plot using envelope simulation
              'pmultinom', #Calculate cdf of multinomial dist
              'ggplot2', # Data visualisation
              'dplyr', # Data processing
              'tidyr', # Data processing
              'moments', # Functions to calculate: moments, Pearson's kurtosis, Geary's kurtosis and skewness; tests related to them
              'gridExtra', # Provides a number of user-level functions to work with ``grid'' graphics, notably to arrange multiple grid-based plots on a page, and draw tables.
              'here', 
              'scatterplot3d', 
              'VGAM',
              'haven', 
              'mclogit'
)

install.packages(setdiff(packages, rownames(installed.packages())), dependencies = T)
lapply(packages, library, character.only = TRUE)
if (!require("ggflags")) {
  devtools::install_github("rensa/ggflags")
  library(ggflags)
}

# Main functions ---------------------------


# plot settings ---------------------------
pallete = RColorBrewer::brewer.pal(9, "Set1")[ c(3, 1, 9, 6, 8, 5, 2) ]

theme_new <- function(base_size = 20, base_family = "Arial"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text = element_text(size = 20, colour = "grey30"),
      legend.key=element_rect(colour=NA, fill =NA),
      axis.line = element_line(colour = 'black'),
      axis.ticks =         element_line(colour = "grey20"),
      plot.title.position = 'plot',
      legend.position = "bottom"
    )
}
