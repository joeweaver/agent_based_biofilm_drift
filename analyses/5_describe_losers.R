library(readr) # handle csv file import
library(here)  # manage paths, keep @JennyBryan from incinerating our machine
library(tidyr)
library(dplyr)   # work with tidy data
library(gt) # table generation
library(logger) # record info
library(magrittr)
library(ggplot2) # figure generation
library(ggh4x)

# Describe the distribution of biggest losers and thrive-survive-lanquish
# classes within the seed simulations for the baseline

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# TODO renumber this so it comes earlier in the sequence
# Set up and start logging ---------
start_logging("5_describe_losers.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# check for 1 and only 1 biggest loser per run
sim_results %<>% filter() %>% # baseline mu and ks
  group_by(nbugs,seed, spacing)  %>%   # for each simulation
  mutate(n_losers = sum(biggest_loser))

poor_losers <- sim_results %>% filter(n_losers != 121)  # 121 seeds per grouping
if( nrow(poor_losers) != 0){
  # TODO also check
  warning('Some runs did not have 1 and only 1 biggest loser. See poor_losers')
  unique(poor_losers$nbugs)
  fours<-poor_losers %>% filter(nbugs == 16)
  fives<-poor_losers %>% filter(nbugs == 25)
  unique(fours$seed)
  unique(fives$seed)
}
# TODO split this out into a data quality checker
fivesok<-sim_results %>% filter(n_losers == 121) %>% filter(nbugs==25) %>% filter(spacing == 5)
unique(fivesok$seed)
unique(fours$colony)
a<-fivesok %>% filter(seed == 1031)
foursok<-sim_results %>% filter(n_losers == 121) %>% filter(nbugs==16) %>% filter(spacing == 10)
unique(foursok$seed)
unique(foursok$colony)
# TODO generate histogram of biggest loser positions, show it's evenly
# distributed
probs <- sim_results %>% filter(mu == 0.00028, ks == 3.50e-05) %>% # baseline mu and ks
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser))


create_dir_if_not_exist(here::here('output'))

# baseline mu and ks
baseline_sim_low_nutrients <- sim_results %>% filter(mu == 0.00028, ks == 3.50e-05)

# TODO need high nutrients as well

short_nbugs_label <- function(string) {
  glue::glue("Initial Pop.: {string}")
}

# create plot
p <- baseline_sim_low_nutrients %>%
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser)) %>%
  mutate(expected = 1/nbugs*121) %>%
  ggplot(aes(x=colony,y=times_lost, color=factor(spacing))) +
  geom_hline(aes(yintercept=expected),alpha=0.5,linetype="dashed")+
  geom_point(alpha=0.6,size=0.9) +
  facet_grid(cols = vars(nbugs),scales="free_x", labeller = labeller( .cols = short_nbugs_label )) +
  xlab("Colony ID") +
  ylab("Occurrences as Lowest Abundance")+
  guides(color=guide_legend(title="Spacing (diameters):")) +
  ggh4x::facetted_pos_scales(x = list(
    nbugs == 4 ~ scale_x_continuous(breaks = c(1,2,3,4)),
    nbugs == 9 ~ scale_x_continuous(breaks = c(1,3,5,7,9)),
    nbugs == 16 ~ scale_x_continuous(breaks = c(4,8,12,16)),
    nbugs == 25 ~ scale_x_continuous(breaks = c(5,10,15,20,25))
  )) +
  theme(legend.position = c(0.625,0.95),#"top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.background=element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size=7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(color="black", fill="white", linetype="solid"),
        strip.text.x = element_text(size=6))

ggsave(here::here('output','baseline_pct_losses.tiff'),width=3.5,height=3.5, units="in",
       dpi=330)
ggsave(here::here('output','baseline_pct_losses.png'),width=3.5,height=3.5, units="in",
       dpi=330)
ggsave(here::here('output','baseline_pct_losses.pdf'),width=3.5,height=3.5, units="in",
       dpi=330)


# Chisq test
chisq_tests <- baseline_sim_low_nutrients %>%
  group_by(nbugs, spacing, colony)  %>%   # for each simulation
  summarize(times_lost = sum(biggest_loser)) %>%
  group_by(nbugs,spacing) %>%
  summarize(cs = chisq.test(times_lost,simulate.p.value=TRUE)$p.value) %>%
  mutate(sig = cs < 0.05/nrow(chisq_tests))

## convert to wider and save as CSV for a table in SI
chisq_tests %>% pivot_wider(names_from = spacing,
                            values_from = c(cs, sig)) %>%
  write_csv(here::here("output","chisq_biggest_loser.csv"))



# Get an idea of the percentages & abundances of survivorship classes
pctSurvivorClass<-baseline_sim_low_nutrients %>%
  group_by(seed,nbugs, spacing)  %>%   # for each simulation
  count(category) %>%
  mutate(pct=n/nbugs)%>%
  select(-n)%>%
  pivot_wider(names_from = category,
              values_from = pct) %>%
  ungroup(seed)%>%
  replace(is.na(.), 0) %>%
  summarize(pctThrive=mean(Thriving),
            pmThrive=sd(Thriving),
            pctSurvive=mean(Surviving),
            pmSurvive=sd(Surviving),
            pctLanguish=mean(Languishing),
            pmLanguish=sd(Languishing)) %>%
    mutate(absThrive=pctThrive*nbugs,
           abspmThrive=pmThrive*nbugs,
           absSurvive=pctSurvive*nbugs,
           abspmSurvive=pmSurvive*nbugs,
           absLanguish=pctLanguish*nbugs,
           abspmLanguishe=pmLanguish*nbugs)

pctSurvivorClass %>% select(-starts_with("abs")) %>%
  pivot_wider(names_from = spacing,
                            values_from = c(pctThrive,pmThrive,pctSurvive,pmSurvive,pctLanguish,pmLanguish)) %>%
  write_csv(here::here("output","precent_survivor_clases.csv"))





