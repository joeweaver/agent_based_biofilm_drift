library(readr) # handle csv file import
library(dplyr) # work with tidy data
library(ggplot2) # figure generation
library(latex2exp) # generating better text labels for figures
library(here) # manage paths, keep @JennyBryan from incinerating our machine
library(logger) # logging

# Create a bunch of 2d probability maps

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("prob_maps.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output', 'prob_maps'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))


# process for plotting
# determine the number of times a biggest loser was categorized as thriving
prob_thrive <- sim_results %>% filter(biggest_loser) %>%
  mutate(catnum = case_when(category == "Thriving" ~ 1,
                            TRUE ~ 0)) %>%
  group_by(nbugs, spacing, ks,mu,yield) %>%
  summarise(prob_thrive = mean(catnum)) %>%
  mutate(mu_pct = (mu - base_mu)/base_mu) %>% # discuss mu and ks in terms of pct change
  mutate(ks_pct = (ks - base_ks)/base_ks)


p <- ggplot(data=prob_thrive, aes(x=mu_pct*100, y=ks_pct*100,fill=cut(prob_thrive,c(-0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
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
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + facet_grid(rows = vars(nbugs), cols=vars(spacing))
p
ggsave(here::here('output','prob_maps',"5x5_5_prob_map.png"),units="in",width=8,height=8,dpi=300)

