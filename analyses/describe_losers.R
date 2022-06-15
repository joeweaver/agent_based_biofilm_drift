library(readr) # handle csv file import
library(here)  # manage paths, keep @JennyBryan from incinerating our machine
library(stringr) # extract values from strings. e.g. num bugs from filenames
library(tidyr)
library(dplyr)   # work with tidy data
library(gt) # table generation
library(logger) # record info
library(magrittr)
library(ggplot2) # figure generation

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first


# Set up logging ---------
# TODO send warnings to to their own file to make sure they're easy to browse
# and not lost in info clutter?
# Create the logs dir if it does not exist
if(!dir.exists(here::here('logs'))){
  dir.create(here::here('logs'))
}
logfile <- here::here("logs",
                      paste0(format(Sys.time(), "%Y-%m-%d-%H-%M-%S_"),
                             "describe_losers.txt"))
log_appender(appender_file(logfile))
log_info(strrep('#', 52)  )
log_info('Beginning run.')
log_info(strrep('#', 52)  )

# Keep track of session info ---------
log_info(strrep('#', 52)  )
log_info('Session information')
log_info(strrep('#', 52)  )
# hacky. multi-line comlex info is a pain to get nicely formatted into logger
sink(logfile, append=TRUE)
sessionInfo()
sink()

sweep_dir <- "sweep_colony_outcomes"

# read a simulation result csv and annotate the bug arrangements inferred
# from the filename
# sweep_colony_outcomes_2x3_4.csv would have a MxN 2x3 bug grid with
# a spacing of 4 bug diameters and a total number of 12 bugs
read_results <- function(filename){
  res <- read_csv(filename,col_types = "idddiclddd")
  file_meta <- str_match_all(
    filename,
    "sweep_colony_outcomes_(\\d*)x(\\d*)_(\\d*\\.*\\d*)\\.csv")
  m = as.integer(file_meta[[1]][2])
  n = as.integer(file_meta[[1]][3])
  spacing = as.double(file_meta[[1]][4])
  return(res %>% mutate(m = m,n = n,spacing = spacing,nbugs = m*n))
}

# gather all simulation results
sim_results <- list.files(here::here("data",sweep_dir), full.names = TRUE) %>%
  lapply(read_results) %>%
  bind_rows

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

create_dir <-function(location){
  if(!dir.exists(here::here(location))){
    dir.create(here::here(location))
    log_info('Created ',location)
  }
}

create_dir(here::here('output'))

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
