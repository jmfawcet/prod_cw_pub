
# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(ez)
library(openxlsx)
library(brms)

# Useful Functions --------------------------------------------------------

short_angle = function(x, y)
{
  return(((x - y + 180) %% 360) - 180)
}

degree_to_rad = function(x)
{
  return(x*pi / 180)
}

rad_to_degrees = function(x)
{
  return(x*180 / pi)
}

kappa_to_sd = function(x)
{
  return(sqrt(1/x))
}

convert_log_kappa_to_degrees = function(x)
{
  return(rad_to_degrees(kappa_to_sd(exp(x))))
}

# Graphing constants ------------------------------------------------------

main_col='black'
APA_options = theme(panel.grid.major=element_line(colour = NA)
                    , panel.grid.minor=element_line(colour = NA)
                    , axis.line=element_line(colour=main_col)
                    , panel.background = element_rect(fill = 'transparent', colour=NA)
                    , plot.background = element_rect(fill = 'transparent', colour=NA)
                    , axis.text.x = element_text(colour = main_col, vjust=1, )
                    , axis.text.y = element_text(colour = main_col, hjust=1)
                    , axis.title.x = element_text(colour = main_col, hjust=0.5)
                    , axis.title.y = element_text(colour = main_col, angle=90, )
                    , axis.ticks.x =element_line(colour = main_col) 
                    , axis.ticks.y =element_line(colour = main_col) 
                    , legend.background = element_rect(colour=NA, fill=NA)
                    
)


