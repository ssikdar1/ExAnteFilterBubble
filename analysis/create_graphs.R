library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)
library(xtable)
library(gridExtra)


#' Calculate diversity mean w/ Confidence Intervals 
#'
#' @param rec_data Dataframe
#' @param var  one of rho | beta | sigma
#' @return Dataframe w/ mean and confidence intervals
#' @examples
#' get_stats_diversity(rec_data, "rho")
get_stats_diversity <- function(df, var){
  if (var == "rho") { 
    df <- df %>% group_by(formatted_regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(formatted_regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(formatted_regime, sigma) 
  } else if (var == "alpha") {
    df <- df %>% group_by(formatted_regime, alpha)
  }
  
  df <- df %>% 
    mutate(diversity_mean = mean(pop_diversity_avg),
           diversity_sd = sd(pop_diversity_avg),
           diversity_n = n()) %>%
    mutate(diversity_se =  diversity_sd/ sqrt(diversity_n),
           lower_ci_diversity = diversity_mean - qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se,
           upper_ci_diveristy = diversity_mean + qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se) %>% 
    slice(1)
  return(df)
  
}

#' Calculate follow recommendation mean w/ Confidence Intervals 
#'
#' @param rec_data Dataframe
#' @param var  one of rho | beta | sigma
#' @return Dataframe w/ mean and confidence intervals
#' @examples
#' get_stats_rec(rec_data, "rho")
get_stats_rec <- function(df, var){
  if (var == "rho") { # there must be a better way
    df <- df %>% group_by(formatted_regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(formatted_regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(formatted_regime, sigma) 
  } else if (var == "alpha") {
    df <- df %>% group_by(formatted_regime, alpha)
  }
  
  df <- df %>% 
    mutate(rec_mean = mean(pop_follow_rec_avg),
           rec_sd = sd(pop_follow_rec_avg),
           rec_n = n()) %>%
    mutate(rec_se =  rec_sd/ sqrt(rec_n),
           lower_ci_rec = rec_mean - qt(1 - (0.05 / 2), rec_n - 1) * rec_se,
           upper_ci_rec = rec_mean + qt(1 - (0.05 / 2), rec_n - 1) * rec_se) %>% 
    slice(1)
  return(df)
}


#' Calculate welfare mean w/ Confidence Intervals 
#'
#' @param rec_data Dataframe
#' @param var  one of rho | beta | sigma
#' @return Dataframe w/ mean and confidence intervals
#' @examples
#' get_stats_welfare(rec_data, "rho")
get_stats_welfare <- function(df, var){
  if (var == "rho") {
    df <- df %>% group_by(formatted_regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(formatted_regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(formatted_regime, sigma) 
  } else if (var == "alpha") {
    df <- df %>% group_by(formatted_regime, alpha)
  }
  
  df <- df %>%
    mutate(welfare_mean = mean(pop_welfare_avg),
           welfare_sd = sd(pop_welfare_avg),
           welfare_n = n()) %>% 
    mutate(welfare_se =  welfare_sd/ sqrt(welfare_n),
           lower_ci_welfare = welfare_mean - qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se,
           upper_ci_welfare = welfare_mean + qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se) %>%
    slice(1)
  return(df)
}


#' Calculate homogeneity w/ Confidence Intervals 
#'
#' @param stats_df Dataframe
#' @param var  one of rho | beta | sigma
#' @return Dataframe w/ mean and confidence intervals
#' @examples
#' get_stats_homogeneity(stats_df, "rho")
get_stats_homogeneity <- function(df, var){
  ## HOMOGENEITY DATA
  if (var == "rho") {
    df <- df %>% group_by(formatted_regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(formatted_regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(formatted_regime, sigma) 
  } else if (var == "alpha") {
    df <- df %>% group_by(formatted_regime, alpha)
  }
  df <- df %>%
    mutate(jaccard_mean = mean(pop_jaccard_avg),
           jaccard_sd = sd(pop_jaccard_avg),
           jaccard_n = n()) %>%
    mutate(jaccard_se =  jaccard_sd/ sqrt(jaccard_n),
           lower_ci_jaccard = jaccard_mean - qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
           upper_ci_jaccard = jaccard_mean + qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se) %>% 
    slice(1)
  return(df)
}

#' Graph diversity as function of `var`
#'
#' @param df Dataframe
#' @param var  one of rho | beta | sigma
#' @param use_hrbrthemes bool to toggle hrbrthemes
#' @return ggplot graph
graph_stats_diversity <- function(df, var, N, T_val, use_hrbrthemes){
  # https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
  # graph w/ themes
  
  g <- ggplot(df, aes_string(x=var, y="diversity_mean")) +
    geom_line(aes(colour=formatted_regime)) +
    geom_errorbar(aes(ymin=lower_ci_diversity, ymax=upper_ci_diveristy), width=.02,
                  position=position_dodge(.9)) + 
    labs(x=var, y="diversity",
         title=paste("N =", N, "T =", T_val, "Diversity",sep=" "))
  if(use_hrbrthemes){
    g <- g + theme_ipsum_rc()
  }
  return(g)
}

#' Graph rec mean as function of `var`
#'
#' @param df Dataframe
#' @param var  one of rho | beta | sigma
#' @param use_hrbrthemes bool to toggle hrbrthemes
#' @return ggplot graph
graph_stats_rec <- function(df, var, N, T_val, use_hrbrthemes){
  g <- ggplot(df, aes_string(x=var, y="rec_mean")) +
    geom_line() +
    geom_errorbar(aes(ymin=lower_ci_rec, ymax=upper_ci_rec), width=.02,
                  position=position_dodge(.9)) + 
    labs(x=var, y="diversity",
         title=paste("N =", N, "T =", T_val, "Follow Rec",sep=" "))
  if(use_hrbrthemes){
    g <- g + theme_ipsum_rc()
  }
  return(g)
}


#' Graph welfare as function of `var`
#'
#' @param df Dataframe
#' @param var  one of rho | beta | sigma
#' @param use_hrbrthemes bool to toggle hrbrthemes
#' @return ggplot graph
graph_stats_welfare <- function(d, var, N, T_val, use_hrbrthemes){
  g <- ggplot(d, aes_string(x=var, y="welfare_mean")) +
    geom_line(aes(colour=formatted_regime)) +
    geom_errorbar(aes(ymin=lower_ci_welfare, ymax=upper_ci_welfare), width=.02, position=position_dodge(.9)) + 
    labs(x=var, y="welfare",
         title=paste( "Welfare",sep=" ")) 
  if(use_hrbrthemes){
    g <- g + theme_ipsum_rc()
  }
  return(g)
  
}


#' Graph homogeneity as function of `var`
#'
#' @param df Dataframe
#' @param var  one of rho | beta | sigma
#' @param use_hrbrthemes bool to toggle hrbrthemes
#' @return ggplot graph
graph_stats_homo <- function(df, var, N, T_val){
  
  if (var == "rho") {
    g <- ggplot(df, aes(x=rho, y=jaccard_mean)) +
      geom_line(aes(colour=formatted_regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="rho", y="homogeneity",
           title=paste("N =", N, "T =", T_val, "Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "beta") {
    g <- ggplot(df, aes(x=beta, y=jaccard_mean)) +
      geom_line(aes(colour=formatted_regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="beta", y="homogeneity",
           title=paste("Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "sigma") {
    g <- ggplot(df, aes(x=sigma, y=jaccard_mean)) +
      geom_line(aes(colour=formatted_regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="sigma", y="homogeneity",
           title=paste("N =", N, "T =", T_val, "Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  } else if (var == "alpha") {
    g <- ggplot(df, aes(x=alpha, y=jaccard_mean)) +
      geom_line(aes(colour=formatted_regime)) +
      geom_errorbar(aes(ymin=lower_ci_jaccard, ymax=upper_ci_jaccard), width=.02,
                    position=position_dodge(.9)) + 
      labs(x="alpha", y="homogeneity",
           title=paste("N =", N, "T =", T_val, "Homogeneity",sep=" ")) + theme_ipsum_rc()
    return(g)
  }
}



scatter <- function(df, var_x, var_y, use_hrbrthemes, title){
  g <- ggplot(df, aes_string(var_x, var_y)) +
    geom_smooth() + #stat_summary_bin(fun.y = "mean", geom="point", bins = 10) +
    labs(x=var_x, y=var_y
         #subtitle="A plot that is only useful for demonstration purposes",
         #caption="Brought to you by the letter 'g'"
    )
  if(use_hrbrthemes){
    g <- g + theme_ipsum_rc(axis_title_size = 14,plot_title_size = 24)
  }
  return(g)
}


process_rec_homo_data <- function(N, t, use_hrbrthemes){
  
  rec_data <- read.csv(paste(WORKING_DIR, "data/rec_data_N_",N,"_t_",t,".csv", sep=""))
  rec_data <- rec_data %>% mutate(formatted_regime = ifelse(regime == "omni", "Omniscient", ifelse(regime == "no_rec", "No Rec", "Partial")))
  rec_data$follow_recommendation <- as.numeric(levels(rec_data$follow_recommendation)[rec_data$follow_recommendation])
  rec_data <- rec_data %>% group_by(pop_idx, formatted_regime, rho, beta, sigma, alpha) %>% summarize(pop_welfare_avg = mean(welfare), pop_diversity_avg = mean(diversity_score), pop_follow_rec_avg = mean(follow_recommendation))
  
  homogeneity <- read.csv(paste(WORKING_DIR, "data/homogeneity_data_N_",N,"_t_",t,".csv", sep=""))
  homogeneity <- homogeneity %>% mutate(formatted_regime = ifelse(regime == "omni", "Omniscient", ifelse(regime == "no_rec", "No Rec", "Partial")))
  homogeneity <- homogeneity %>% group_by(pop_idx, formatted_regime, rho, beta, sigma, alpha) %>% summarize(pop_jaccard_avg = mean(jaccard))
  
  # Calculate the marginal variables
  variables=list("rho", "beta", "sigma", "alpha")
  metrics=list("diversity", "welfare", "homogeneity")
  
  for (metric in metrics){
    
    print(metric)
    for(variable in variables){
      print(variable)
      
      file_name <- paste(WORKING_DIR,"figures/",variable,"_", metric,"_N_", N, "_T_", t, ".jpeg", sep="")
      print(file_name)
      
      if (metric == "diversity"){
        g <- graph_stats_diversity(get_stats_diversity(rec_data, variable), variable, N, t, use_hrbrthemes)
        ggsave(filename=file_name, plot=g)
      }
      else if (metric == "welfare"){
        #welfare
        g <- graph_stats_welfare(get_stats_welfare(rec_data, variable), variable, N, t, use_hrbrthemes)
        ggsave(filename=file_name, plot=g)
        
        rec_file_name <- paste(WORKING_DIR, "figures/", variable,"_", metric,"_N_", N, "_T_", t, "_rec.jpeg", sep="")
        print(rec_file_name)
        g <- graph_stats_rec(get_stats_rec(filter(rec_data, formatted_regime == "Partial"), variable), variable, N, t,use_hrbrthemes)
        ggsave(filename=rec_file_name, plot=g)
      } 
      else{
        # homo
        g <- graph_stats_homo(get_stats_homogeneity(homogeneity, variable), variable, N, t)
        ggsave(filename=file_name, plot=g)
      }
    } # END for(variable in variables)
    
    
  }
  rec_abbrv_to_name <- function(abbrv) {
    if (abbrv == "Omniscient") {
      return("Omniscient")
    } else if (abbrv == "No Rec") {
      return("No Recommendation")
    } else if (abbrv == "Partial") {
      return("Partial")
    }
  }
  
  rec_policies <- c("Omniscient", "No Rec", "Partial")
  for (rec_policy in rec_policies) {
    title <- paste("Diversity_vs_Welfare_N_",N, "_T_",T,"_", rec_abbrv_to_name(rec_policy), sep="")
    ggsave(
      filename=paste(WORKING_DIR, "figures/", title, ".jpeg", sep=""), 
      plot=scatter(filter(rec_data, formatted_regime == rec_policy), "pop_diversity_avg", "pop_welfare_avg", use_hrbrthemes, title) 
    )
  }
  
  df <- rec_data %>% group_by(formatted_regime) %>%
    mutate(welfare_mean = mean(pop_welfare_avg),
           welfare_sd = sd(pop_welfare_avg),
           welfare_n = n(),
           diversity_mean = mean(pop_diversity_avg),
           diversity_sd = sd(pop_diversity_avg),
           diversity_n = n()) %>% 
    mutate(welfare_se =  welfare_sd/ sqrt(welfare_n),
           welfare_mul = qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se,
           welfare_str = paste(round(welfare_mean, 2), "\\pm", round(welfare_mul, 4))) %>%
    mutate(diversity_se  =  diversity_sd/ sqrt(diversity_n),
           diversity_mul = qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se,
           diversity_str = paste(round(diversity_mean, 2), "\\pm", round(diversity_mul, 4))) %>%
    slice(1)
  
  df <- df[,c("formatted_regime", "welfare_str", "diversity_str")]
  
  homo <- homogeneity %>% group_by(formatted_regime) %>%
    mutate(jaccard_mean = mean(pop_jaccard_avg),
           jaccard_sd = sd(pop_jaccard_avg),
           jaccard_n = n()) %>%
    mutate(jaccard_se  =  jaccard_sd/ sqrt(jaccard_n),
           jaccard_mul = qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
           jaccard_str = paste(round(jaccard_mean, 3), "\\pm", round(jaccard_mul, 4))) %>%
    slice(1)
  homo <- homo[,c("formatted_regime", "jaccard_str")]
  
  df <- merge(df, homo, by="formatted_regime")
  xtable(df)
  
} # END process_rec__homo_data


### Consumption Diversity


get_stats_consumption_diversity_time_path <- function(df, var){
  if (var == "rho") {
    df <- df %>% group_by(formatted_regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(formatted_regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(formatted_regime, sigma) 
  } else if (var == "alpha") {
    df <- df %>% group_by(formatted_regime, alpha)
  }
  df <- df %>%
    mutate(local_move_mean = mean(local_move_05),
           local_move_sd = sd(local_move_05),
           local_move_n = n()) %>%
    mutate(local_move_se =  local_move_mean/ sqrt(local_move_sd),
           lower_ci_local_move = local_move_mean - qt(1 - (0.05 / 2), local_move_n - 1) * local_move_se,
           upper_ci_local_move = local_move_mean + qt(1 - (0.05 / 2), local_move_n - 1) * local_move_se) %>% 
    slice(1)
  return(df)
}

graph_stats_consumption_diversity_time_path <- function(df, var, N, T_val){
  
  g <- ggplot(df, aes_string(x=var, y="local_move_mean")) +
    geom_line() +
    geom_errorbar(aes(
      ymin=lower_ci_local_move, 
      ymax=upper_ci_local_move), width=.02,
      position=position_dodge(.9)) + 
    labs(x=var, y="local_move_mean",
         title=paste("N =", N, "T =", T_val, "local_move_mean",sep=" "))
  if(use_hrbrthemes){
    g <- g + theme_ipsum_rc()
  }
  return(g)
}

process_time_path <- function(N, T_val, use_hrbrthemes){
  
  # define the metric for consumption distribution
  time_data <- read.csv(paste(WORKING_DIR, "data/time_path_n_",N,"_t_", T_val, ".csv", sep=""))
  time_data <- time_data %>% mutate(formatted_regime = ifelse(regime == "omni", "Omniscient", ifelse(regime == "no_rec", "No Rec", "Partial")))
  # Look at consumption distribution held over all parameters
  time_data <- time_data %>% filter(t > 0) # local move isn't defined at the first time step so drop it to properly have smoothed plots
  
  # 
  # https://stackoverflow.com/questions/21066077/remove-fill-around-legend-key-in-ggplot
  
  no_rho_or_alpha <- filter(time_data, rho == 0.0 & alpha == 0.0)
  no_rho <- filter(time_data, rho == 0.0 & alpha > 0.0)
  no_alpha <- filter(time_data, alpha == 0.0 & rho > 0.0)
  g <- ggplot(no_rho_or_alpha, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=formatted_regime)) + theme_bw() +
    ggtitle(paste("Consumption Distance, No Correlation or Risk Aversion")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom", legend.title = element_blank())  + 
    guides(color=guide_legend(override.aes=list(fill=NA)))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"no_correlation_risk_aversion.jpeg", sep=""), 
    plot=g
  )
  
  g <- ggplot(no_rho, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=formatted_regime)) + theme_bw() +
    ggtitle(paste("Consumption Distance, No Correlation")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom", legend.title = element_blank())  + 
    guides(color=guide_legend(override.aes=list(fill=NA)))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"no_correlation.jpeg", sep=""), 
    plot=g
  )
  
  g <- ggplot(no_alpha, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=formatted_regime)) + theme_bw() +
    ggtitle(paste("Consumption Distance, No Risk Aversion")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom", legend.title = element_blank())  + 
    guides(color=guide_legend(override.aes=list(fill=NA)))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"no_risk_aversion.jpeg", sep=""), 
    plot=g
  )
  
  g <- ggplot(time_data, aes(x=t, y=instantaneous_welfare_average)) +
    geom_smooth(aes(colour=formatted_regime)) + theme_bw() +
    ggtitle(paste("Welfare")) + xlab("t") + ylab("Average Welfare")
  g <- g + theme(legend.position="bottom", legend.title = element_blank())  + 
    guides(color=guide_legend(override.aes=list(fill=NA)))
  ggsave(
    filename= paste(WORKING_DIR, "figures/welfare_N_", N,"T_", T_val,".jpeg", sep=""), 
    plot=g
  )
  
  no_rec <- filter(time_data, formatted_regime == "No Rec" & beta == 0.4)
  g <- ggplot(no_rec, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=as.factor(alpha))) + theme_bw() +
    ggtitle(paste("No Recommendation Consumption Distance")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom")  + 
    guides(color=guide_legend(override.aes=list(fill=NA), title="Risk Aversion"))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"_no_rec_ra.jpeg", sep=""), 
    plot=g
  )
  
  no_rec <- filter(time_data, formatted_regime == "No Rec")
  g <- ggplot(partial_rec, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=as.factor(beta))) + theme_bw() +
    ggtitle(paste("No Rec Recommendation Consumption Distance")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom")  + 
    guides(color=guide_legend(override.aes=list(fill=NA), title="Beta"))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"_partial_rec_ra.jpeg", sep=""), 
    plot=g
  )
  
  omni_rec <- filter(time_data, formatted_regime == "Omniscient"  & beta == 0.4)
  g <- ggplot(omni_rec, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=as.factor(alpha))) + theme_bw() +
    ggtitle(paste("Omniscient Recommendation Consumption Distance")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom")  + 
    guides(color=guide_legend(override.aes=list(fill=NA), title="Risk Aversion"))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"_omni_rec_ra.jpeg", sep=""), 
    plot=g
  )
  
  # vary alpha
  
  no_rec <- filter(time_data, formatted_regime == "No Rec" & rho > 0.0 & beta == 0.4)
  g <- ggplot(no_rec, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=as.factor(alpha))) + theme_bw() +
    ggtitle(paste("No Recommendation Consumption Distance")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom")  + 
    guides(color=guide_legend(override.aes=list(fill=NA), title="Risk Aversion"))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"_no_rec_ra.jpeg", sep=""), 
    plot=g
  )
  
  partial_rec <- filter(time_data, formatted_regime == "Partial" & rho > 0.0 & beta == 0.4)
  g <- ggplot(partial_rec, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=as.factor(alpha))) + theme_bw() +
    ggtitle(paste("Partial Recommendation Consumption Distance")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom")  + 
    guides(color=guide_legend(override.aes=list(fill=NA), title="Risk Aversion"))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"_partial_rec_ra.jpeg", sep=""), 
    plot=g
  )
  
  omni_rec <- filter(time_data, formatted_regime == "Omniscient" & rho > 0.0)
  g <- ggplot(omni_rec, aes(x=t, y=consumption_dist)) +
    geom_smooth(aes(colour=as.factor(beta))) + theme_bw() +
    ggtitle(paste("Omniscient Recommendation Consumption Distance")) + xlab("t") + ylab("Average Consumption Distance")
  g <- g + theme(legend.position="bottom")  + 
    guides(color=guide_legend(override.aes=list(fill=NA), title="Risk Aversion"))
  ggsave(
    filename= paste(WORKING_DIR, "figures/consumption_dist_N_", N,"T_", T_val,"_omni_rec_ra.jpeg", sep=""), 
    plot=g
  )
  ## Now look at the marginals of local_move over "rho", "beta", "sigma", "alpha"
  variables=list("rho", "beta", "sigma", "alpha")
  for(variable in variables){
    tmp <- get_stats_consumption_diversity_time_path(t, variable)
    g <- graph_stats_consumption_diversity_time_path(tmp, variable, N, t)
    file_name <- paste(WORKING_DIR, "figures/", variable, "_time_path_local_search_N_", N,"T_", t,".jpeg", sep="")
    ggsave(filename=file_name, plot=g)
  }
  
  ## TODO try different combinations of variables like rho and beta
  ## See how the fractional search changes
  
  partial <- filter(time_data, formatted_regime == "Partial")
  partial$follow_recommendation <- as.numeric(partial$follow_recommendation) - 1
  g <- ggplot(partial, aes(x=t, y=follow_recommendation)) +
    geom_smooth() +
    labs(x="t", y="Fraction Following Recommendation"#, title="Recommendation Effectiveness"
    )
  if(use_hrbrthemes){
    g <- g + theme_ipsum_rc(plot_title_size = 24, axis_title_size = 14)
  }
  ggsave(
    filename=paste(WORKING_DIR, "figures/rec_obedience_N_", N, "_T_", T, ".jpeg", sep=""), 
    plot=g
  )
}



### MAIN

USER <- Sys.getenv( "USER" )
WORKING_DIR <- paste("/Users/",USER, "/Desktop/ExAnteFilterBubble/", sep="")
setwd(WORKING_DIR)
use_hrbrthemes <- FALSE

N_s <- list(200)

for(N in N_s){
  #process_rec_homo_data(N, 20, use_hrbrthemes)
  process_time_path(N, 20, use_hrbrthemes)
}
