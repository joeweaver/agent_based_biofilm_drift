# Simulation results are saved to osf.io.
# Before anything else can run, we need to grab them.
# During project development, we also need to run this to update
# the data dir with new results, maintain parity with how they're organized.

library(osfr)  # access data on osf
library(dplyr) # osf can be navigated using dplyr style stuff
library(here)  # manage paths, keep @JennyBryan from incinerating our machine

# Until we put this up on a preprint, the data repo is private, so
# we need to authenticate with an OSF PAT, stored in .Renvironment.
# This requirement can go away once it's public, but there's still benefits to
# using a PAT

osf_auth(Sys.getenv("osf_PAT"))

# this location should not change
project_guid <- "fch3z"

# Download simulation results into data dir
#
# by default we are having the download call error out instead of clobbering
# it's up to you to decide if it's best to
# a. clear the data/sweep_colony_outcomes directory  OR
# b. set the conflicts param to 'skip' or 'overwrite'
osf_retrieve_node(project_guid) %>%
  osf_ls_files() %>%
  filter(name == "sweep_colony_outcomes") %>%
  osf_download(here::here("data"))



