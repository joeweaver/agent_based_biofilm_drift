create_dir_if_not_exist <- function(location, do_log = FALSE){
  if(!dir.exists(here::here('logs'))){
    dir.create(here::here('logs'))
    if(do_log){
      log_info('Created ',location)
    }
  }
}

library(logger) # logging

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