library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)

library(gridExtra)


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

# https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
# graph w/ themes
graph_stats_diversity <- function(df, var, use_hrbrthemes){
  T_val <- df$T[1]
  N <- df$N[1]
  
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

graph_stats_rec <- function(df, var, use_hrbrthemes){
  T_val <- df$T[1]
  N <- df$N[1]

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

# graph w/ themes
graph_stats_welfare <- function(d, var, use_hrbrthemes){
  T_val <- d$T[1]
  N <- d$N[1]
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



#df2 <- get_stats_diversity(rec_data, 0.9)
#g <- graph_stats_diversity(df2)
#ggsave(filename=paste(WORKING_DIR, "figures/rho_09_diversity.jpeg", sep=""), plot=g)



## HOMOGENEITY DATA

get_stats_homogeneity <- function(stats_df, var){
  df <- stats_df
  if (var == "rho") {
    df <- df %>% group_by(regime, rho) 
  } else if (var == "beta") {
    df <- df %>% group_by(regime, beta) 
  } else if (var == "sigma") {
    df <- df %>% group_by(regime, sigma) 
  }
  df <- df %>%
      mutate(jaccard_mean = mean(jaccard),
             jaccard_sd = sd(jaccard),
             jaccard_n = n()) %>%
      mutate(jaccard_se =  jaccard_sd/ sqrt(jaccard_n),
             lower_ci_jaccard = jaccard_mean - qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
             upper_ci_jaccard = jaccard_mean + qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se) %>% 
      slice(1)
  return(df)
}


graph_stats_homo <- function(df, var){
  T_val <- df$T[1]
  N <- df$N[1]
  
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
  }
}

#WORKING_DIR <- "/Users/guyaridor/Desktop/recommender_systems/rec_sys_conf_paper/"
WORKING_DIR <- "/Users/ssikdar/Downloads/econ/bubbles/"

rec_data <- read.csv(paste(WORKING_DIR, "rec_data_t_25.csv", sep=""))
rec_data <- rec_data %>% mutate(formatted_regime = ifelse(regime == "rec", "Omniscient", ifelse(regime == "no_rec", "No Rec", "Partial")))

homogeneity <- read.csv(paste(WORKING_DIR, "homogeneity_data_t_25.csv", sep=""))
homogeneity <- homogeneity%>% mutate(formatted_regime = ifelse(regime == "rec", "Omniscient", ifelse(regime == "no_rec", "No Rec", "Partial")))


# CREATE GRAPHS
use_hrbrthemes <- TRUE
variables=list("rho", "beta", "sigma")
metrics=list("diversity", "welfare", "homogeneity")

for(metric in metrics){
  
  print(metric)
  for(variable in variables){
    print(variable)
    
    file_name <- paste(WORKING_DIR,"figures/",variable,"_", metric, "_t_25.jpeg", sep="")
    print(file_name)
    
    if (metric == "diversity"){
      g <- graph_stats_diversity(get_stats_diversity(rec_data, variable), variable, use_hrbrthemes)
      ggsave(filename=file_name, plot=g)
    }
    else if (metric == "welfare"){
      #welfare
      g <- graph_stats_welfare(get_stats_welfare(rec_data, variable), variable,  use_hrbrthemes)
      ggsave(filename=file_name, plot=g)
      
      rec_file_name <- paste(WORKING_DIR, "figures/", variable, "_rec_t_25.jpeg", sep="")
      print(rec_file_name)
      g <- graph_stats_rec(get_stats_rec(filter(rec_data, regime == "partial"), variable), variable,  use_hrbrthemes)
      ggsave(filename=rec_file_name, plot=g)
    } 
    else{
      # homo
      g <- graph_stats_homo(get_stats_homogeneity(homogeneity, variable), variable)
      ggsave(filename=file_name, plot=g)
    }
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

rec_abbrv_to_name <- function(abbrv) {
  if (abbrv == "rec") {
    return("Omniscient")
  } else if (abbrv == "no_rec") {
    return("No Recommendation")
  } else if (abbrv == "partial") {
    return("Partial")
  }
}

rec_policies <- c("rec", "no_rec", "partial")
for (rec_policy in rec_policies) {
  title <- paste("Diversity vs Welfare -", rec_abbrv_to_name(rec_policy))
  ggsave(
    filename=paste(WORKING_DIR, "figures/", title, ".jpeg", sep=""), 
    plot=scatter(filter(rec_data, regime == rec_policy), "diversity_score", "welfare", use_hrbrthemes, title)
  )
}

df <- rec_data %>% group_by(regime) %>%
  mutate(welfare_mean = mean(welfare),
         welfare_sd = sd(welfare),
         welfare_n = n(),
         diversity_mean = mean(diversity_score),
         diversity_sd = sd(diversity_score),
         diversity_n = n()) %>% 
  mutate(welfare_se =  welfare_sd/ sqrt(welfare_n),
         welfare_mul = qt(1 - (0.05 / 2), welfare_n - 1) * welfare_se,
         welfare_str = paste(round(welfare_mean, 2), "\\pm", round(welfare_mul, 4))) %>%
  mutate(diversity_se  =  diversity_sd/ sqrt(diversity_n),
          diversity_mul = qt(1 - (0.05 / 2), diversity_n - 1) * diversity_se,
          diversity_str = paste(round(diversity_mean, 2), "\\pm", round(diversity_mul, 4))) %>%
  slice(1)

df <- df[,c("regime", "welfare_str", "diversity_str")]

homo <- homogeneity %>% group_by(regime) %>%
  mutate(jaccard_mean = mean(jaccard),
         jaccard_sd = sd(jaccard),
         jaccard_n = n()) %>%
  mutate(jaccard_se  =  jaccard_sd/ sqrt(jaccard_n),
         jaccard_mul = qt(1 - (0.05 / 2), jaccard_n - 1) * jaccard_se,
         jaccard_str = paste(round(jaccard_mean, 3), "\\pm", round(jaccard_mul, 4))) %>%
  slice(1)
homo <- homo[,c("regime", "jaccard_str")]

df <- merge(df, homo, by="regime")
library(xtable)
xtable(df)

time_dat <- read.csv(paste(WORKING_DIR, "time_path_25.csv", sep=""))
time_dat <- time_dat %>% mutate(formatted_regime = ifelse(regime == "rec", "Omniscient", ifelse(regime == "no_rec", "No Rec", "Partial")))
time_dat <- time_dat %>% mutate(local_move = as.numeric(consumption_dist < 10))

g <- ggplot(time_dat, aes(x=t, y=local_move)) +
    geom_smooth(aes(colour=formatted_regime)) +
    labs(x="t", y="Fraction of Local Search"
    )
if(use_hrbrthemes){
  g <- g + theme_ipsum_rc()
}
ggsave(
  filename=paste(WORKING_DIR, "figures/local_moves_25.jpeg", sep=""), 
  plot=g
)

partial <- filter(time_dat, regime == "partial")
partial$follow_recommendation <- as.numeric(partial$follow_recommendation) - 1
g <- ggplot(partial, aes(x=t, y=follow_recommendation)) +
  geom_smooth() +
  labs(x="t", y="Fraction Following Recommendation"#, title="Recommendation Effectiveness"
  )
if(use_hrbrthemes){
  g <- g + theme_ipsum_rc(plot_title_size = 24, axis_title_size = 14)
}
ggsave(
  filename=paste(WORKING_DIR, "figures/rec_obedience_25.jpeg", sep=""), 
  plot=g
)



g <- graph_stats_welfare(get_stats_welfare(rec_data, 'beta'), 'beta',  TRUE) + theme(legend.position="bottom" ,legend.title = element_blank())

h<- graph_stats_homo(get_stats_homogeneity(homogeneity, 'beta'), 'beta') + theme(legend.position="none" ,legend.title = element_blank())

# https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(g)

p3 <- grid.arrange(arrangeGrob(g + theme(legend.position="none"),
                               h + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))


ggsave(
  filename=paste(WORKING_DIR, "figures/welfare_homo_combo.jpeg", sep=""), 
  plot=p3
)


