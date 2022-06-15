library(readr) # handle csv file import
library(dplyr) # work with tidy data
library(ggplot2) # figure generation
library(latex2exp) # generating better text labels for figures
library(here) # manage paths, keep @JennyBryan from incinerating our machine

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# TODO this is common code and should be extracted
# # Set up logging ---------
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

# TODO dry
# Create relevant output dirs ---------
# TODO logger should be param
create_dir <-function(location){
  if(!dir.exists(here::here(location))){
    dir.create(here::here(location))
    log_info('Created ',location)
  }
}

create_dir(here::here('output'))
create_dir(here::here('output', 'prob_maps'))

# TODO right now uses one hardcoded sim result, should be applied to all
# Read simulation results ---------
res <- read_csv(here::here("data","sweep_colony_outcomes","sweep_colony_outcomes_5x5_5.csv"))
filtered <- res %>% filter(biggest_loser)

both_modified_colony <- filtered %>%
  mutate(catnum = case_when(category == "Thriving" ~ 1,
                            TRUE ~ 0))
probs <- both_modified_colony %>%  group_by(ks,mu,yield) %>% summarise(prob_thrive = mean(catnum))

base_mu = 0.00028
base_ks = 3.5e-5
probs <- probs %>% mutate(mu_pct = (mu - base_mu)/base_mu) %>%
  mutate(ks_pct = (ks - base_ks)/base_ks)
p <- ggplot(data=probs, aes(x=mu_pct*100, y=ks_pct*100,fill=cut(prob_thrive,c(-0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  geom_tile()+
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  #ylab("Change in substrate affinity (ks)") +
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{max}$)")) +
  #xlab("Change in maximum specific growth rate (mu max)")+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave(here::here('output','prob_maps',"5x5_5_prob_map.png"),units="in",width=8,height=8,dpi=300)

