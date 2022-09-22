library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(broom)
library(tidyr)
library(purrr)
library(RColorBrewer)
library(itsadug)
library(mgcv)
# Fit sigmoid curves along lines of equal Ks for simulations

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("4_EDA_gam.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)
probs$nbugs <- factor(probs$nbugs)
probs$spacing <- factor(probs$spacing)
linear_model <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+ s(ks_pct,mu_pct), data = ss_d)
summary(linear_model)
AIC(linear_model)
plot(linear_model)
vis.gam(linear_model, theta = 120, n.grid = 50, lwd = 0.4)
a<-data.frame(mu_pct=0.5,ks_pct=0.1,nbugs=4,spacing=5)
library(tidymv)
pg<-predict_gam(linear_model)

ss<-pg%>%filter(nbugs==4,spacing==2.5)
ss_d<-probs%>%filter(nbugs==4,spacing==2.5)
a<-inner_join(ss_d,ss)

