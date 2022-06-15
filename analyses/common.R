# Contains functions common to multiple analyses steps

library(logger) # logging
library(stringr) # string matching

base_mu = 0.00028
base_ks = 3.5e-5

# Create a directory if it does not exist, optionally log to log_info
create_dir_if_not_exist <- function(location, do_log = FALSE){
  if(!dir.exists(here::here('logs'))){
    dir.create(here::here('logs'))
    if(do_log){
      log_info('Created ',location)
    }
  }
}

# TODO send warnings to to their own file to make sure they're easy to browse
# and not lost in info clutter?

# Create a log file ending in log_suffix and prepended with a timestamp
# in the given logdir. Lists that a run is beginning and gives session info
start_logging <- function(log_suffix, logdir=here::here('logs')){
  create_dir_if_not_exist(logdir)
  logfile <- here::here("logs",
                        paste0(format(Sys.time(), "%Y-%m-%d-%H-%M-%S_"),
                               log_suffix))
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
  print(sessionInfo())
  sink()
}

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

# given a directory containing sweep results, read all of them and combine
# into one dataframe
get_all_results <- function(outcomes_directory){
  sweep_dir <- "sweep_colony_outcomes"
  # gather all simulation results
  sim_results <- list.files(outcomes_directory, full.names = TRUE) %>%
    lapply(read_results) %>%
    bind_rows
  return(sim_results)
}