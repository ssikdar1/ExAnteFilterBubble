library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)

WORKING_DIR <- "/Users/guyaridor/Desktop/recommender_systems/rec_sys_conf_paper/"
rec_data <- read.csv(paste(WORKING_DIR, "rec_data.csv", sep=""))

# like above but get mean std, n , confidence intervals for graph
get_stats_diversity <- function(rec_data, var){
  df <- rec_data %>%
    group_by(regime, UQ(var)) %>% 
    mutate(diversity_mean = mean(diversity_score),
           diversity_sd = sd(diversity_score),
           diversity_n = n()) %>% 
    mutate(diversity_se =  diversity_sd/ sqrt(diversity_n),
           lower_ci_diversity = diversity_mean - qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se,
           upper_ci_diveristy = diversity_mean + qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se) %>% 
    slice(1)
  return(df)
}

# graph w/ themes
graph_stats_diversity <- function(df, var){
  g <- ggplot(df2, aes(x=UQ(var), y=diversity_mean)) +
    geom_line(aes(colour=regime)) +
    geom_errorbar(aes(ymin=lower_ci_diversity, ymax=upper_ci_diveristy), width=.02,
                  position=position_dodge(.9)) + 
    labs(x=var, y="diversity",
         title=paste("N=500,T=25, Diversity",sep=" ")) +
    theme_ipsum_rc()
  return(g)
}

get_stats_welfare <- function(rec_data, var){
  df <- rec_data %>%
    group_by(regime, UQ(var)) %>%
    mutate(welfare_mean = mean(welfare),
          welfare_sd = sd(welfare),
           welfare_n = n()) %>% 
    mutate(welfare_se =  welfare_sd/ sqrt(welfare_n),
           lower_ci_welfare = welfare_mean - qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se,
           upper_ci_welfare = welfare_mean + qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se) %>%
    slice(1)
  return(df)
}

# graph w/ themes
graph_stats_welfare <- function(d, var){
  g <- ggplot(d, aes(x=UQ(var), y=welfare_mean)) +
    geom_line(aes(colour=regime)) +
    geom_errorbar(aes(ymin=lower_ci_welfare, ymax=upper_ci_welfare), width=.02, position=position_dodge(.9)) + 
    labs(x=var, y="welfare",
         title=paste("N=500,T=25, Welfare",sep=" "))+
    theme_ipsum_rc()
  return(g)
}

#df2 <- get_stats_diversity(rec_data, 0.9)
#g <- graph_stats_diversity(df2)
#ggsave(filename=paste(WORKING_DIR, "figures/rho_09_diversity.jpeg", sep=""), plot=g)


  # diversity

g <- graph_stats_diversity(get_stats_diversity(rec_data, "rho"), "rho")
ggsave(filename=paste(WORKING_DIR, "figures/rho_diversity.jpeg", sep=""), plot=g)
  
g <- graph_stats_diversity(get_stats_diversity(rec_data, "beta"), "beta")
ggsave(filename=paste(WORKING_DIR, "figures/beta_diversity.jpeg", sep=""), plot=g)
  #welfare
g <- graph_stats_welfare(get_stats_welfare(rec_data, "rho"), "rho")
ggsave(filename=paste(WORKING_DIR, "figures/rho_welfare.jpeg", sep=""), plot=g)

g <- graph_stats_welfare(get_stats_welfare(rec_data, "beta"), "beta")
ggsave(filename=paste(WORKING_DIR, "figures/beta_welfare.jpeg", sep=""), plot=g)
  


# TODO circle graphs


