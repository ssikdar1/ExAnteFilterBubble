library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)


get_stats_diversity <- function(rec_data, var){
  df <- rec_data
  if (var == "rho") { 
    df <- df %>% group_by(regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(regime, beta) 
  } else if (var == "sigma") {
    df <- rec_data %>% group_by(regime, sigma) 
  }
  
  df <- df %>% 
    mutate(diversity_mean = mean(diversity_score),
           diversity_sd = sd(diversity_score),
           diversity_n = n()) %>%
    mutate(diversity_se =  diversity_sd/ sqrt(diversity_n),
           lower_ci_diversity = diversity_mean - qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se,
           upper_ci_diveristy = diversity_mean + qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se) %>% 
    slice(1)
  return(df)
  
}

get_stats_rec <- function(rec_data, var){
  rec_data$follow_recommendation <- as.numeric(levels(rec_data$follow_recommendation)[rec_data$follow_recommendation])
  
  df <- rec_data
  
  if (var == "rho") { # there must be a better way
    df <- df %>% group_by(rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(sigma) 
  }
  
  df <- df %>% 
    mutate(rec_mean = mean(follow_recommendation),
           rec_sd = sd(follow_recommendation),
           rec_n = n()) %>%
    mutate(rec_se =  rec_sd/ sqrt(rec_n),
           lower_ci_rec = rec_mean - qt(1 - (0.05 / 2), rec_n - 1) * rec_se,
           upper_ci_rec = rec_mean + qt(1 - (0.05 / 2), rec_n - 1) * rec_se) %>% 
    slice(1)
  return(df)
}


# graph w/ themes
graph_stats_diversity <- function(df, var){
  T_val <- df$T[1]
  N <- df$N[1]
  if (var == "rho") {
    g <- ggplot(df, aes(x=rho, y=diversity_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_diversity, ymax=upper_ci_diveristy), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="rho", y="diversity",
           title=paste("N =", N, "T =", T_val, "Diversity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "beta") {
    g <- ggplot(df, aes(x=beta, y=diversity_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_diversity, ymax=upper_ci_diveristy), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="beta", y="diversity",
           title=paste("N =", N, "T =", T_val, "Diversity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "sigma") {
    g <- ggplot(df, aes(x=sigma, y=diversity_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_diversity, ymax=upper_ci_diveristy), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="sigma", y="diversity",
           title=paste("N =", N, "T =", T_val, "Diversity",sep=" ")) + theme_ipsum_rc()
  }
  
}

graph_stats_rec <- function(df, var){
  T_val <- df$T[1]
  N <- df$N[1]
  if (var == "rho") {
    g <- ggplot(df, aes(x=rho, y=rec_mean)) +
      geom_line() +
      geom_errorbar(aes(ymin=lower_ci_rec, ymax=upper_ci_rec), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="rho", y="diversity",
           title=paste("N =", N, "T =", T_val, "Follow Rec",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "beta") {
    g <- ggplot(df, aes(x=beta, y=rec_mean)) +
      geom_line() +
      geom_errorbar(aes(ymin=lower_ci_rec, ymax=upper_ci_rec), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="beta", y="rec",
           title=paste("N =", N, "T =", T_val, "Follow Rec",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "sigma") {
    g <- ggplot(df, aes(x=sigma, y=rec_mean)) +
      geom_line() +
      geom_errorbar(aes(ymin=lower_ci_rec, ymax=upper_ci_rec), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="sigma", y="rec",
           title=paste("N =", N, "T =", T_val, "Follow Rec",sep=" ")) + theme_ipsum_rc()
    return(g)
  }
  
}

get_stats_welfare <- function(rec_data, var){
  df <- rec_data
  if (var == "rho") {
    df <- df %>% group_by(regime, rho) 
    
  } else if (var == "beta") {
    df <- df %>% group_by(regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(regime, sigma) 
  }
  
  df <- df %>%
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
  T_val <- d$T[1]
  N <- d$N[1]
  if (var == "rho") {
    g <- ggplot(d, aes(x=rho, y=welfare_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_welfare, ymax=upper_ci_welfare), width=.02, position=position_dodge(.9)) + 
      labs(x="rho", y="welfare",
           title=paste("N =", N, "T =", T_val, "Welfare",sep=" ")) + theme_ipsum_rc()
    # theme_ipsum_rc()
    return(g)
  } else if (var == "beta") {
    g <- ggplot(d, aes(x=beta, y=welfare_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_welfare, ymax=upper_ci_welfare), width=.02, position=position_dodge(.9)) + 
      labs(x="beta", y="welfare",
           title=paste("N =", N, "T =", T_val, "Welfare",sep=" ")) + theme_ipsum_rc()
    #      theme_ipsum_rc()
    return(g)
  }else if (var == "sigma") {
    g <- ggplot(d, aes(x=sigma, y=welfare_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_welfare, ymax=upper_ci_welfare), width=.02, position=position_dodge(.9)) + 
      labs(x="sigma", y="welfare",
           title=paste("N =", N, "T =", T_val, "Welfare",sep=" ")) + theme_ipsum_rc()
    #      theme_ipsum_rc()
    return(g) 
  }
}

#df2 <- get_stats_diversity(rec_data, 0.9)
#g <- graph_stats_diversity(df2)
#ggsave(filename=paste(WORKING_DIR, "figures/rho_09_diversity.jpeg", sep=""), plot=g)



## HOMOGENEITY DATA

get_stats_homogeneity <- function(stats_df, var){
  if (var == "rho") { # there must be a better way
    df <- stats_df %>%
      group_by(regime, rho) %>% 
      mutate(jaccard_mean = mean(jaccard),
             jaccard_sd = sd(jaccard),
             jaccard_n = n()) %>%
      mutate(jaccard_se =  jaccard_sd/ sqrt(jaccard_n),
             lower_ci_jaccard = jaccard_mean - qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
             upper_ci_jaccard = jaccard_mean + qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se) %>% 
      slice(1)
    return(df)
  } else if (var == "beta") {
    df <- stats_df %>%
      group_by(regime, beta) %>% 
      mutate(jaccard_mean = mean(jaccard),
             jaccard_sd = sd(jaccard),
             jaccard_n = n()) %>%
      mutate(jaccard_se =  jaccard_sd/ sqrt(jaccard_n),
             lower_ci_jaccard = jaccard_mean - qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
             upper_ci_jaccard = jaccard_mean + qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se) %>% 
      slice(1)
    return(df)
  } else if (var == "sigma") {
    df <- stats_df %>%
      group_by(regime, sigma) %>% 
      mutate(jaccard_mean = mean(jaccard),
             jaccard_sd = sd(jaccard),
             jaccard_n = n()) %>%
      mutate(jaccard_se =  jaccard_sd/ sqrt(jaccard_n),
             lower_ci_jaccard = jaccard_mean - qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
             upper_ci_jaccard = jaccard_mean + qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se) %>% 
      slice(1)
    return(df)
  }
}


graph_stats_homo <- function(df, var){
  T_val <- df$T[1]
  N <- df$N[1]
  if (var == "rho") {
    g <- ggplot(df, aes(x=rho, y=jaccard_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="rho", y="homogeneity",
           title=paste("N =", N, "T =", T_val, "Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "beta") {
    g <- ggplot(df, aes(x=beta, y=jaccard_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="beta", y="homogeneity",
           title=paste("N =", N, "T =", T_val, "Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "sigma") {
    g <- ggplot(df, aes(x=sigma, y=jaccard_mean)) +
      geom_line(aes(colour=regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="sigma", y="homogeneity",
           title=paste("N =", N, "T =", T_val, "Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  }
}

#WORKING_DIR <- "/Users/guyaridor/Desktop/recommender_systems/rec_sys_conf_paper/"
WORKING_DIR <- "/Users/ssikdar/Downloads/econ/bubbles/"

rec_data <- read.csv(paste(WORKING_DIR, "rec_data_t_15.csv", sep=""))
homogeneity <- read.csv(paste(WORKING_DIR, "homogeneity_data_t_15.csv", sep=""))


# CREATE GRAPHS
use_hrbrthemes <- TRUE
variables=list("rho", "beta", "sigma")
metrics=list("diversity", "welfare", "homogeneity")

for(metric in metrics){
  
  print(metric)
  for(variable in variables){
    print(variable)
    
    file_name <- paste(WORKING_DIR,"figures/",variable,"_", metric, "_t_15.jpeg", sep="")
    print(file_name)
    
    if (metric == "diversity"){
      g <- graph_stats_diversity(get_stats_diversity(rec_data, variable), variable)
      ggsave(filename=file_name, plot=g)
    }
    else if (metric == "welfare"){
      #welfare
      g <- graph_stats_welfare(get_stats_welfare(rec_data, variable), variable)
      ggsave(filename=file_name, plot=g)
      
      rec_file_name <- paste(WORKING_DIR, "figures/", variable, "_rec_t_15.jpeg", sep="")
      print(rec_file_name)
      g <- graph_stats_rec(get_stats_rec(filter(rec_data, regime == "partial"), variable), variable)
      ggsave(filename=rec_file_name, plot=g)
    } 
    else{
      # homo
      g <- graph_stats_homo(get_stats_homogeneity(homogeneity, variable), variable)
      ggsave(filename=file_name, plot=g)
    }
  }
}



