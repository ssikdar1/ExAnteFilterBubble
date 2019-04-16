library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)

WORKING_DIR <- "/Users/ssikdar/Downloads/"
rec_data <- read.csv(paste(WORKING_DIR, "rec_data.csv", sep=""))

get_mean_diversity <- function(df) {
  d <- df %>% group_by(regime, pop_idx, rho, beta) %>% 
    mutate(diversity = mean()) %>%
    slice(1) %>%
    group_by(regime, rho, beta) %>% 
    summarize(diversity = mean(diversity))
  return(d)
}

# like above but get mean std, n , confidence intervals for graph
get_stats_diversity <- function(rec_data, rho_rate){
  df <- rec_data %>% 
    filter(rho == rho_rate) %>%
    group_by(regime, pop_idx, rho, beta) %>% 
    mutate(diversity = mean(diversity_score)) %>%
    slice(1) %>%
    group_by(regime, rho, beta) %>% 
    mutate(diversity_mean = mean(diversity),
           diversity_sd = sd(diversity),
           diversity_n = n()) %>% 
    mutate(diversity_se =  diversity_sd/ sqrt(diversity_n),
           lower_ci_diversity = diversity_mean - qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se,
           upper_ci_diveristy = diversity_mean + qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se)
  return(df)
}


# graph w/ themes
graph_stats_diversity <- function(df,rho){
  g <- ggplot(df2, aes(x=beta, y=diversity_mean)) +
    geom_line(aes(colour=regime)) +
    geom_errorbar(aes(ymin=lower_ci_diversity, ymax=upper_ci_diveristy), width=.02,
                  position=position_dodge(.9)) + 
    labs(x="beta", y="diversity",
         title=paste("Rho=",rho,"N=500,T=25, Diversity",sep=" "),
         #subtitle="A plot that is only useful for demonstration purposes",
         caption="Columbia University") + 
    theme_ipsum_rc()
  return(g)
}

#TODO welfare
get_mean_welfare <- function(df) {
  d <- df %>% group_by(regime, pop_idx, rho, beta) %>% 
    mutate(welfare = mean(welfare)) %>%
    slice(1) %>%
    group_by(regime, rho, beta) %>% 
    summarize(welfare = mean(welfare))
  return(d)
}

get_stats_welfare <- function(rec_data, rho_rate){
  df <- rec_data %>% 
    filter(rho == rho_rate) %>%
    group_by(regime, pop_idx, rho, beta) %>% 
    mutate(welfare_calc = mean(welfare)) %>%
    slice(1) %>%
    group_by(regime, rho, beta) %>% 
    mutate(welfare_mean = mean(welfare_calc),
          welfare_sd = sd(welfare_calc),
           welfare_n = n()) %>% 
    mutate(welfare_se =  welfare_sd/ sqrt(welfare_n),
           lower_ci_welfare = welfare_mean - qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se,
           upper_ci_welfare = welfare_mean + qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se)
  return(df)
}

# graph w/ themes
graph_stats_welfare <- function(d,rho){
  g <- ggplot(d, aes(x=beta, y=welfare_mean)) +
    geom_line(aes(colour=regime)) +
    geom_errorbar(aes(ymin=lower_ci_welfare, ymax=upper_ci_welfare), width=.02, position=position_dodge(.9)) + 
    labs(x="beta", y="welfare",
         title=paste("Rho=",rho,"N=500,T=25, Welfare",sep=" "),
         #subtitle="A plot that is only useful for demonstration purposes",
         caption="Columbia University") + 
    theme_ipsum_rc()
  return(g)
}

#df2 <- get_stats_diversity(rec_data, 0.9)
#g <- graph_stats_diversity(df2)
#ggsave(filename=paste(WORKING_DIR, "figures/rho_09_diversity.jpeg", sep=""), plot=g)


for (rho in c(0.9, 0.5, 0.1)){
  # diversity
  rho_name <- paste("rho_", rho, sep="")
  print(rho_name)
  df2 <- get_stats_diversity(rec_data, rho)
  g <- graph_stats_diversity(df2,rho)
  ggsave(filename=paste(WORKING_DIR, "figures2/", rho_name, "_diversity.jpeg", sep=""), plot=g)
  
  #welfare
  g <- graph_stats_welfare(get_stats_welfare(rec_data,rho), rho)
  ggsave(filename=paste(WORKING_DIR, "figures2/", rho_name, "_welfare.jpeg", sep=""), plot=g)
  
}


# TODO circle graphs


