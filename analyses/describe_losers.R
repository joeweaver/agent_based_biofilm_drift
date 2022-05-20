library(readr) # handle csv file import
library(here)  # manage paths, keep @JennyBryan from incinerating our machine
library(stringr) # extract values from strings. e.g. num bugs from filenames

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

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
  print(file_meta)
  print(m)
  print(n)
  print(spacing)
  return(res %>% mutate(m = m,n = n,spacing = spacing,nbugs = m*n))
}

# gather all simulation results
sim_results <- list.files(here::here("data",sweep_dir), full.names = TRUE) %>%
  lapply(read_results) %>%
  bind_rows

# TODO ensure 1 and only 1 biggest loser per run

# TODO generate histogram of biggest loser positions, show it's evenly
# distributed

