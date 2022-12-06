# Draw a single probability map and illustrate sigmoid fit

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
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first
#
# need to have output/sigmoidfits.csv
# if not, you have to run 2_sigmoid_fitting.R
#
# also needs output/spread_fits.csv
# if not, you have to run 3_prob_map_and_spread_analysis.R

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("7_illustrate_sigmoid_fit.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','eda'))
create_dir_if_not_exist(here::here('output','eda','spread_analysis'))

# read the previously fit spreads
all_spreads<-read_csv(here::here("output","spread_fits.csv"))

# Create probability map --------------------------------------------------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))
probs <- thrive_probabilities(sim_results) %>% mutate(likely=prob_thrive/(1-prob_thrive+1e-6)+1e-6,
                                                      loglike=log(likely))

p_single_map<-ggplot(data=all_spreads %>% filter(nbugs==4) %>% filter(spacing==5), aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=probs,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(prob_thrive,c(-0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
   geom_segment(data=all_spreads%>% filter(nbugs==4) %>% filter(spacing==5) %>% filter(ks_pct==0.0) %>% filter(range_pct %in% c(0.95)),
                aes(x=lo*100,xend=hi*100,y=ks_pct*100,yend=ks_pct*100),linetype="dotted",alpha=0.5)+
   geom_point(data= all_spreads%>% filter(nbugs==4) %>% filter(spacing==5)  %>% filter(ks_pct==0.0), alpha=0.9)+
   geom_point(data=all_spreads%>% filter(nbugs==4) %>% filter(spacing==5) %>% filter(ks_pct==0.0) %>% filter(range_pct %in% c(0.68)),
              aes(x=lo*100),shape=3,alpha=0.75)+
   geom_point(data=all_spreads%>% filter(nbugs==4) %>% filter(spacing==5)  %>% filter(ks_pct==0.0)%>% filter(range_pct %in% c(0.68)),
              aes(x=hi*100),shape=3,alpha=0.75)+
  geom_point(data=all_spreads%>% filter(nbugs==4) %>% filter(spacing==5) %>% filter(ks_pct==0.0) %>% filter(range_pct %in% c(0.95)),
             aes(x=lo*100),shape=8,alpha=0.75)+
  geom_point(data=all_spreads%>% filter(nbugs==4) %>% filter(spacing==5)  %>% filter(ks_pct==0.0)%>% filter(range_pct %in% c(0.95)),
             aes(x=hi*100),shape=8,alpha=0.75)+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  guides(fill=guide_legend(title="Thrive\nProbability")) +
  theme(legend.position = "left",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size=7),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(color="black", fill="white", linetype="solid"))
#p_single_map


sf <- read_csv(here::here("output","sigmoid_fits.csv"))
df <- sf %>% filter(nbugs %in% c(4))%>%
  filter(spacing %in% c(5))
sig_fit <- df %>% filter(ks %in% c(base_ks))
mus=seq(-0.5,0.5,length.out=100)
preds<-logistic_fun(mus,sig_fit$mu_50,sig_fit$steep)
sig_preds<-data.frame(mu_pct=mus,prob_thrive=preds,ks=base_ks)

sig_run <-probs %>% filter(nbugs %in% c(4))%>%
  filter(spacing %in% c(5)) %>% filter(ks %in% c(base_ks))

p_fit<-ggplot(sig_run,aes(x=mu_pct*100,y=prob_thrive,group=ks,color=factor(ks)))+

  geom_line(data=sig_preds,linetype="solid",color="grey")+
  geom_point(size=1.2,color="black",pch=1)+

  geom_segment(x=-.55*100,xend=.119*100,y=0.5,yend=0.5,size=0.5,color="#1b9e77",linetype="dotted")+
  geom_segment(x=.119*100,xend=.119*100,y=0.5,yend=-0.05,size=0.5,color="#1b9e77",linetype="dotted")+
  geom_point(x=0.119*100,y=0.5,color="#1b9e77")+

  geom_segment(x=-0.003377495*100,xend=0.242*100,y=0.84,yend=0.84,size=0.5,color="#d95f02",linetype="dashed")+
  geom_segment(x=-0.003377495*100,xend=0.242*100,y=0.87,yend=0.87,size=0.25,color="#d95f02",linetype="solid",
               arrow=arrow(ends="both",type="closed",length=unit(1,"mm")))+
  geom_segment(x=-0.003377495*100,xend=-0.003377495*100,y=0.16,yend=0.84,size=0.5,color="#d95f02",linetype="dashed")+
  geom_segment(x=0.242*100,xend=0.242*100,y=0.84,yend=0.84,size=0.5,color="#d95f02",linetype="dashed")+
  geom_point(x=-0.003377495*100,y=0.16,color="#d95f02",pch=3)+
  geom_point(x=0.242*100,y=0.84,color="#d95f02",pch=3)+

  geom_segment(x=-0.15*100,xend=0.39*100,y=0.975,yend=0.975,size=0.5,color="#7570b3",linetype="dotdash")+
  geom_segment(x=-0.15*100,xend=0.39*100,y=1.01,yend=1.01,size=0.25,color="#7570b3",linetype="solid",
               arrow=arrow(ends="both",type="closed",length=unit(1,"mm")))+
  geom_segment(x=-0.15*100,xend=-0.15*100,y=0.025,yend=0.975,size=0.5,color="#7570b3",linetype="dotdash")+
  geom_segment(x=0.39*100,xend=0.39*100,y=0.975,yend=0.975,size=0.5,color="#7570b3",linetype="dotdash")+
  geom_point(x=-0.15*100,y=0.025,color="#7570b3",pch=8)+
  geom_point(x=0.39*100,y=0.975,color="#7570b3",pch=8)+
  annotate("text",size=3,label="mu[~~50]",x=20,y=0.5,color="#1b9e77",parse=TRUE)+
  annotate("text",size=3,label="spread[68]",x=10,y=0.9,color="#d95f02",parse=TRUE)+
  annotate("text",size=3,label="spread[95]",x=10,y=1.05,color="#7570b3",parse=TRUE)+
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +

  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  ylab("Probability of transitioning to thriving")+
  guides(color=guide_legend(title="Affinity (Ks)"),
         linetype=guide_legend(title="Affinity (Ks)")) +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line=element_line(),
        axis.text.x = element_text(angle=-45),
        axis.title = element_text(size = 7),
        axis.text = element_text(size=7),
        legend.text = element_text(size=4),
        legend.title = element_text(size=6))
#p_fit


p_illus<-p_single_map+p_fit+ plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 7)) # nb could use this more to reduce code in above

fname <- "sigmoid_illus.png"
floc <- here::here("output",fname)
ggsave(floc, p_illus, width=8,height=4,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "sigmoid_illus.pdf"
floc <- here::here("output",fname)
ggsave(floc, p_illus, width=8,height=4,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "sigmoid_illus.tiff"
floc <- here::here("output",fname)
ggsave(floc, p_illus, width=8,height=4,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))