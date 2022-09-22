library(readr) # handle csv file import
library(here)  # manage paths, keep @JennyBryan from incinerating our machine
library(tidyr)
library(dplyr)   # work with tidy data
library(gt) # table generation
library(logger) # record info
library(magrittr)
library(ggplot2) # figure generation

# Describe the distribution of biggest losers and thrive-survive-langquish
# classes within the seed simulations for the baseline

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("describe_losers.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# check for 1 and only 1 biggest loser per run
sim_results %<>% filter() %>% # baseline mu and ks
  group_by(nbugs,seed, spacing)  %>%   # for each simulation
  mutate(n_losers = sum(biggest_loser))

poor_losers <- sim_results %>% filter(n_losers == 121)  # 121 seeds per grouping
if( nrow(poor_losers) != 0){
  # TODO also check
  warning('Some runs did not have 1 and only 1 biggest loser. See poor_losers')

}

# TODO generate histogram of biggest loser positions, show it's evenly
# distributed
probs <- sim_results %>% filter(mu == 0.00028, ks == 3.50e-05) %>% # baseline mu and ks
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser))


create_dir_if_not_exist(here::here('output'))

# baseline mu and ks
baseline_sim_low_nutrients <- sim_results %>% filter(mu == 0.00028, ks == 3.50e-05)

# TODO need high nutrients as well

# create plot
p<- baseline_sim_low_nutrients %>%
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser)/121) %>%
  ggplot(aes(x=colony,y=times_lost, color=factor(spacing))) + geom_point() +
  facet_grid(cols = vars(nbugs))
ggsave(here::here('output','baseline_pct_losses.png'))

a<-baseline_sim_low_nutrients %>%
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser))  %>%
  group_by(nbugs,spacing) %>%
  summarize(sum(times_lost))
# Chisq test
chisq_tests <- baseline_sim_low_nutrients %>%
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser))  %>%
  group_by(nbugs,spacing) %>%
  summarize(cs = chisq.test(times_lost)$p.value) %>%
  mutate(sig = cs < 0.05/nrow(chisq_tests))


chisq_for_doc <- chisq_tests %>% pivot_wider(names_from = nbugs, values_from =cs) %>%
  gt() %>%
  fmt_number(
    columns = everything(),
    decimals = 2) %>%
  tab_spanner(
    label = "Initial Sites",
    columns = c(2:5))
gtsave(a,here::here('output','chisq_drift.html'))

# TODO histogram of thrive/not thrive
#
probs <- baseline_sim_low_nutrients %>% filter(nbugs==4) %>%
  filter(colony < 5) %>%# baseline mu and ks
  group_by(seed,nbugs, spacing)  %>%   # for each simulation
  count(category) %>%
  mutate(pct = n/nbugs)

probs %>% ggplot(aes(x=category,y=pct)) +
  geom_boxplot()

probs <- baseline_sim_low_nutrients %>% filter(nbugs==16) %>%

  group_by(seed,nbugs, spacing)  %>%   # for each simulation
  count(category) %>%
  mutate(pct = n/nbugs)

probs %>% ggplot(aes(x=category,y=n)) +
  geom_jitter(height=0)+ylim(0,9)
