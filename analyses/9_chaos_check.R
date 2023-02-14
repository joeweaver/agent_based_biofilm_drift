# Compare results of 3x3_10 seeds with mu changed by a very small fraction
# Purpose is to check for chaos in system

library(readr) # handle csv file import
library(dplyr)  # work with tidy data
library(ggplot2) # figure generation
library(latex2exp)  # use LaTeX expressions for creating axis labels
library(here) # manage paths, keep @JennyBryan from incinerating our machine
library(broom) # helps fits work well with dataframes
library(tidyr) # used for pivots
library(RColorBrewer) # palette management in plots
library(itsadug) # GAM
library(mgcv) # GAM
library(tidymv) # GAME
library(ggpubr) # plotting
library(ggtext) # plotting
library(cowplot) # plotting
library(patchwork) # plotting
library(viridis) # plotting
library(ggpp) # plotting
library(tools) # md5 checksum



# precondition
# ./data should contain 'chaos_check_seeds.csv' and 'chaos_check_small_change_outcomes.csv'
#
# load common code
source('common.R')

# Set up and start logging ---------
start_logging("9_chaos_check.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','eda'))
create_dir_if_not_exist(here::here('output','eda','spread_analysis'))

# read the final results of the original runs
original_seeds<-read_csv(here::here("data","chaos_check_seeds.csv")) %>% group_by(seed) %>%
  filter(step == max(step)) %>%
  select(1, 6:14) %>% pivot_longer(
    cols = starts_with("het"),
    names_to = "colony",
    names_pattern = "het(.*)_rv",
    values_to = "rv_chaos") %>%
  mutate(colony = as.double(colony))

original_seeds<-original_seeds %>% group_by(seed) %>%
  mutate(rank_chaos = dense_rank(rv_chaos))

delta_outcomes<-read_csv(here::here("data","chaos_check_small_change_outcomes.csv")) %>%
  select(-ks, -mu, -yield, -expected, -cph)  %>%
  group_by(seed) %>%
  mutate(rank = dense_rank(rv))

combined <- inner_join(original_seeds,delta_outcomes)
which(combined$rank != combined$rank_chaos)

combined <- combined %>%
   mutate(delta_rv = rv-rv_chaos)

combined %>% group_by(rank) %>%
  summarise(av = mean(delta_rv),sd = sd(delta_rv),av_rv =mean(rv),sd_rv=sd(rv))
