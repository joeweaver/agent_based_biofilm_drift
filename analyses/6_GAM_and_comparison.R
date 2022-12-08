# Select a generalized additive model, attempt to predict simulation
# results, and compare with MLR

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
library(tools) #md5sum


# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first
#
# need to have output/sigmoidfits.csv
# if not, you have to run 1_sigmoid_fitting.R
#
# also needs MLR predictions output/mlr_predictions.csv
# if not, you have to run 5_multiple_logistic_regression.R

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("6_GAM_and_comparison.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','si'))
create_dir_if_not_exist(here::here('output','fits'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)
probs$nbugs <- as.double(probs$nbugs)#factor(probs$nbugs)
probs$spacing <- as.double(probs$spacing)# factor(probs$spacing)

# Screen models (% deviation and AIC) --------
# takes a while, uncomment as needed,note that fitted models are saved
# as RDS in output/fits
# TODO could break this out from the plotting code
# # no-interact, all smooth
# gam.no_i.all_smooth <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+s(spacing,k=3),
#                            data = probs,method = "REML")
# summary(gam.no_i.all_smooth )
# summary(gam.no_i.all_smooth )$p.table
# summary(gam.no_i.all_smooth )$s.table
# AIC(gam.no_i.all_smooth )
# k.check(gam.no_i.all_smooth)
# fname <- here::here("output","fits","gam.no_i.all_smooth")
# write_rds(gam.no_i.all_smooth,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.no_i.all_smooth"), ' MD5Sum: ',
#                md5sum(fname)))
#
# # # no-interact, all linear
# gam.no_i.all_linear <- gam(prob_thrive ~ mu_pct + ks_pct+nbugs+spacing, data = probs,method = "REML")
# summary(gam.no_i.all_linear )
# summary(gam.no_i.all_linear )$p.table
# summary(gam.no_i.all_linear )$s.table
# AIC(gam.no_i.all_linear )
# k.check(gam.no_i.all_linear)
# fname <- here::here("output","fits","gam.no_i.all_linear")
# write_rds(gam.no_i.all_linear,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.no_i.all_linear"), ' MD5Sum: ',
#                md5sum(fname)))
#
# # # no-interact, linear crowd
# gam.no_i.linear_crowd <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+nbugs+spacing, data = probs,method = "REML")
# summary(gam.no_i.linear_crowd )
# summary(gam.no_i.linear_crowd )$p.table
# summary(gam.no_i.linear_crowd )$s.table
# AIC(gam.no_i.linear_crowd )
# k.check(gam.no_i.linear_crowd)
# fname <- here::here("output","fits","gam.no_i.linear_crowd")
# write_rds(gam.no_i.linear_crowd,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.no_i.linear_crowd"), ' MD5Sum: ',
#                md5sum(fname)))
#
# # # all 2way interact, all smooth (except spacing in case of nbugsxspacing)
# gam.2way_i.all_smooth <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+s(spacing,k=3)+
# s(mu_pct,ks_pct) + s(mu_pct,nbugs) + s(mu_pct,spacing) +
#   s(ks_pct,nbugs) + s(ks_pct,spacing) +
#   s(nbugs,by=spacing,k=3), data = probs,method = "REML")
# gam.2way_i.all_smooth <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+
#                                s(mu_pct,ks_pct) + s(mu_pct,nbugs) +
#                                s(ks_pct,nbugs) , data = probs%>%filter(spacing==2.5),method = "REML")
# summary(gam.2way_i.all_smooth )
# summary(gam.2way_i.all_smooth )$p.table
# summary(gam.2way_i.all_smooth )$s.table
# AIC(gam.2way_i.all_smooth )
# k.check(gam.2way_i.all_smooth)
# fname <- here::here("output","fits","gam.2way_i.all_smooth")
# write_rds(gam.2way_i.all_smooth,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.2way_i.all_smooth"), ' MD5Sum: ',
#                md5sum(fname)))
#
# # # all 2way interact, linear crowd
# gam.2way_i.linear_crowd <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+nbugs+spacing+
#                                s(mu_pct,ks_pct) + s(mu_pct,by=nbugs) + s(mu_pct,by=spacing) +
#                                s(ks_pct,by=nbugs) + s(ks_pct,by=spacing) +
#                                nbugs*spacing, data = probs, method = "REML")
# summary(gam.2way_i.linear_crowd )
# summary(gam.2way_i.linear_crowd )$p.table
# summary(gam.2way_i.linear_crowd )$s.table
# AIC(gam.2way_i.linear_crowd )
# k.check(gam.2way_i.linear_crowd)
# fname <- here::here("output","fits","gam.2way_i.linear_crowd")
# write_rds(gam.2way_i.linear_crowd,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.2way_i.linear_crowd"), ' MD5Sum: ',
#                md5sum(fname)))
#
# # # w 3way interact, all smooth  (except spacing in case of nbugsxspacing)
# gam.3way_i.all_smooth  <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+s(spacing,k=3)+
#                                  s(mu_pct,ks_pct) + s(mu_pct,nbugs) + s(mu_pct,spacing) +
#                                  s(ks_pct,nbugs) + s(ks_pct,spacing) +
#                                  s(nbugs,by=spacing,k=3)+
#                                  te(mu_pct, ks_pct,spacing,d=c(1,2)) +
#                                  te(mu_pct, ks_pct,nbugs,d=c(1,2))
#                                  , data = probs, method = "REML")
# summary(gam.3way_i.all_smooth  )
# summary(gam.3way_i.all_smooth  )$p.table
# summary(gam.3way_i.all_smooth  )$s.table
# AIC(gam.3way_i.all_smooth  )
# k.check(gam.3way_i.all_smooth )
# fname <- here::here("output","fits","gam.3way_i.all_smooth")
# write_rds(gam.3way_i.all_smooth,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.3way_i.all_smooth"), ' MD5Sum: ',
#                md5sum(fname)))
#
# # # w 3way interact, linear crowd
# gam.3way_i.linear_crowd  <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+nbugs+spacing+
#                                 s(mu_pct,ks_pct,k=60) + s(mu_pct,by=nbugs) + s(mu_pct,by=spacing) +
#                                 s(ks_pct,by=nbugs) + s(ks_pct,by=spacing) +
#                                 nbugs*spacing+
#                                 + te(mu_pct, ks_pct,by=spacing) +
#                                 + te(mu_pct, ks_pct,by=nbugs)
#                               , data = probs, method = "REML")
# summary(gam.3way_i.linear_crowd  )
# summary(gam.3way_i.linear_crowd  )$p.table
# summary(gam.3way_i.linear_crowd  )$s.table
# AIC(gam.3way_i.linear_crowd  )
# k.check(gam.3way_i.linear_crowd )
# fname <- here::here("output","fits","gam.3way_i.linear_crowd")
# write_rds(gam.3way_i.linear_crowd,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.3way_i.linear_crowd"), ' MD5Sum: ',
#                md5sum(fname)))

# backwards step reduction, all smooth  (except spacing in case of nbugsxspacing)
# gam.bsr.3way_i.all_smooth  <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+s(spacing,k=3)+
#                                     s(mu_pct,ks_pct,k=60) + s(mu_pct,spacing) +
#                                     s(ks_pct,spacing) +
#                                     s(nbugs,by=spacing,k=3)+
#                                     te(mu_pct, ks_pct,spacing,d=c(1,2)) +
#                                     te(mu_pct, ks_pct,nbugs,d=c(1,2))
#                                   , data = probs, method = "REML")
# summary(gam.bsr.3way_i.all_smooth  )
# summary(gam.bsr.3way_i.all_smooth  )$p.table
# summary(gam.bsr.3way_i.all_smooth  )$s.table
# AIC(gam.bsr.3way_i.all_smooth  )
# k.check(gam.bsr.3way_i.all_smooth )
# fname <- here::here("output","fits","gam.bsr.3way_i.all_smooth")
# write_rds(gam.bsr.3way_i.all_smooth,fname)
# log_info(paste('Wrote', file.path("output","fits","gam.bsr.3way_i.all_smooth"), ' MD5Sum: ',
#                md5sum(fname)))

# read the model
gam.bsr.3way_i.all_smooth<-read_rds( here::here("output","fits","gam.bsr.3way_i.all_smooth"))

# Make predictions from selected model ------------------------------------
log_info("Using gam.bsr.3way_i.all_smooth as fitted model")

pg<-predict_gam(gam.bsr.3way_i.all_smooth ,values=list(mu_pct = unique(probs$mu_pct),
                                         ks_pct =  unique(probs$ks_pct),
                                 nbugs=unique(probs$nbugs),
                                 spacing=unique(probs$spacing)))

combined <-left_join(x=probs,y=pg) %>%
  mutate(pred_error = prob_thrive-fit)


# Plot histogram of errors ------------------------------------------------
spacing_label <- function(string) {
  glue::glue("<span style = 'color:#585858;'>Spacing:&nbsp; <span><span style = 'color:#000000;'> {string}<span> ")
}
nbugs_label <- function(string) {
  glue::glue("<span style = 'color:#585858;'>Ini. Pop.:&nbsp;<span><span style = 'color:#000000;'> {string}<span> ")
}

phist <- ggplot(data=combined,aes(x=pred_error)) +
  geom_histogram(binwidth=0.025,color="black",fill="skyblue")+
  ylab(TeX("Count")) +
  xlab(TeX("Error in Thrive Transition Chance (Prediction - Simulation)")) +
  #ylim(-0.3,0.7)+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1),
  #                   breaks = seq(-0.30,0.60,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.30,0.40,0.10))+
  theme(legend.position = c(0.71,0.85),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12,color="black"),
        axis.text = element_text(size=12,color="black"),
        axis.text.x= element_text(size=8,color="black",angle=360-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=10.5),
        strip.text.y = element_markdown(size=12)) +
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))

fname <- "gam_pred_error_hist.png"
floc <- here::here("output","si",fname)
ggsave(floc, phist, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "gam_pred_error_hist.pdf"
floc <- here::here("output","si",fname)
ggsave(floc, phist, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "gam_pred_error_hist.tiff"
floc <- here::here("output","si",fname)
ggsave(floc, phist, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))


# Create a large-scale standalone plot for SI ---------------
poplabs <- c("Ini. Pop.: 4", "Ini. Pop.: 9", "Ini. Pop.: 16")
names(poplabs) <- c("4","9","16")
p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
                facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "top",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 9,color="black"),
        axis.title.y = element_text(size = 9,color="black"),
        axis.text = element_text(size=6,color="black"),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=8),
        strip.text.y = element_markdown(size=10))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))


fname <- "gam_solo_predictions.png"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "gam_solo_predictions.tiff"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "gam_solo_predictions.pdf"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

# read in the mlr predictions and create composite
# TODO do we actually need to do the join and create composite?
mlr_res <- read_csv(here::here("output","mlr_predictions.csv"))
composite <- left_join(combined, mlr_res) %>%
  mutate(cpred = (fit+fit+preds)/3,
         cerror = prob_thrive - cpred)

# Create prediction plots for comparison figure ---------------
# TODO there's a lot boilerplate here that could be DRY'd
# p_comp_gam_pred<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   coord_fixed()+
#   ylab(TeX("Change in substrate affinity ($K_s$)")) +
#   xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
#   coord_fixed(xlim=c(-50,50))+
#   coord_fixed(ylim=c(-50,50))+
#   scale_y_continuous(labels = function(x) paste0(x, "%"),
#                      breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
#   scale_x_continuous(labels = function(x) paste0(x, "%"),
#                      breaks = c(-50,-30,-10,10,30,50)) +
#   guides(fill=guide_legend(nrow=1)) +
#
#   theme(legend.position = "top",
#         legend.direction = "horizontal",
#
#         legend.title = element_text(size = 5),
#         legend.text = element_text(size = 5),
#         axis.title.x = element_text(size = 7,color="black"),
#         axis.title.y = element_blank(),
#         axis.text = element_text(size=5,color="black"),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(angle=-45),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         strip.text.x = element_markdown(size=6),
#         strip.text.y = element_blank())+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))

#legend <- cowplot::get_legend(p_comp_gam_pred)


p<-ggplot(data=mlr_res, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=mlr_res,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(preds,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6,color="black"),
        axis.text = element_text(size=3.5,color="black"),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7,hjust = 0.5))+
  ggtitle("Multiple Linear Regression (MLR)")+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))

p_comp_gam_pred<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title.x = element_text(size = 5,color="black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=3.5,color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7,hjust = 0.5),
        legend.key.size=unit(4, 'mm'))+
  ggtitle("General Additive Model (GAM)")+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))
#p_comp_gam_pred

p_comp_sim_pred<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(prob_thrive,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "none",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=3.5,color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_markdown(size=4),
        strip.background = element_blank(),
        plot.title = element_text(size=7,hjust = 0.5))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))+
  ggtitle("Simulation 'Ground Truth'")


# plot3 is the collection of predictions from models and the simulation
plot3 <- p + p_comp_gam_pred +p_comp_sim_pred +
   plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A')
#
# ggsave(here::here("output","si","composite_predictions.png"),units="in",width=8,height=3.5,dpi=330)
#
# ggsave(here::here("output","si","composite_predictions.pdf"),units="in",width=8,height=4,dpi=330)

# Create error plots for comparison figure ---------------
p_err<-ggplot(data=mlr_res, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=mlr_res,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(resid,seq(-0.8,0.8,by=0.1)))) +
  scale_fill_manual(values = c("#424086FF",
                               "#3B528BFF", "#33638DFF", "#2C728EFF", "#26828EFF",
                               "#21908CFF", "#1F9F88FF", "#27AD81FF","#3EBC74FF",
                               "#5DC863FF", "#82D34DFF", "#AADC32FF", "#D5E21AFF",
                               "#FDE725FF"))+

  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(title="Error (Model-Simulation)")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6,color="black"),
        axis.text = element_text(size=3.5,color="black"),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_blank(),
        strip.background= element_blank(),
        plot.title = element_text(size=7,hjust = 0.5),
        legend.key.size=unit(4, 'mm'))+
  ggtitle("MLR Error")+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))
#p_err

# update theme for large format image
p_err_solo <- p_err+
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 8,color="black"),
        axis.title.x = element_text(size = 8,color="black"),
        axis.text = element_text(size=7,color="black"),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=7),
        strip.text.y = element_markdown(size=7),
        strip.background= element_blank(),
        plot.title = element_blank(),
        legend.key.size=unit(8, 'mm'))

fname <- "mlr_solo_error.png"
floc <- here::here("output","si",fname)
ggsave(floc, p_err_solo, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mlr_solo_error.pdf"
floc <- here::here("output","si",fname)
ggsave(floc, p_err_solo, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mlr_solo_error.tiff"
floc <- here::here("output","si",fname)
ggsave(floc, p_err_solo, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

# viridis colors for reference. Did some manual stuff to make sure that
# bin colors matched between plots. Breaks was causing weird issues and CBA
# "#440154FF" "#48186AFF" "#472D7BFF" "#424086FF" "#3B528BFF" "#33638DFF" "#2C728EFF" "#26828EFF" "#21908CFF" "#1F9F88FF" "#27AD81FF"
# "#3EBC74FF" "#5DC863FF" "#82D34DFF" "#AADC32FF" "#D5E21AFF" "#FDE725FF"
p_gam_err<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(pred_error,seq(-0.8,0.8,by=0.1)))) +
 scale_fill_manual(values = c("#33638DFF","#2C728EFF","#26828EFF","#21908CFF","#1F9F88FF","#27AD81FF","#3EBC74FF"))+
                    #labels = seq(-0.8,0.8,by=0.1))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(title="Error (Model-Simulation)")) +
  theme(legend.position = "none",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title.x = element_text(size = 5,color="black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=3.5,color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_blank(),
        strip.background= element_blank(),
        plot.title = element_text(size=7,hjust = 0.5))+
  ggtitle("GAM Error")+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))


# update theme for large format image
p_gam_err_solo <- p_gam_err+  theme(legend.position = "right",
                            legend.title = element_text(size = 8),
                            legend.text = element_text(size = 8),
                            axis.title.y = element_text(size = 8,color="black"),
                            axis.title.x = element_text(size = 8,color="black"),
                            axis.text = element_text(size=7,color="black"),
                            axis.text.x = element_text(angle=-45),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"),
                            strip.text.x = element_markdown(size=7),
                            strip.text.y = element_markdown(size=7),
                            strip.background= element_blank(),
                            plot.title = element_blank(),
                            legend.key.size=unit(8, 'mm'))

fname <- "gam_solo_error.png"
floc <- here::here("output","si",fname)
ggsave(floc, p_gam_err_solo, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "gam_solo_error.pdf"
floc <- here::here("output","si",fname)
ggsave(floc, p_gam_err_solo, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "gam_solo_error.tiff"
floc <- here::here("output","si",fname)
ggsave(floc, p_gam_err_solo, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

# Evaluate main effects only ----------------------------------------------

gam.no_i.all_smooth <- readRDS(here::here("output","fits","gam.no_i.all_smooth"))
pg.main<-predict_gam(gam.no_i.all_smooth ,values=list(mu_pct = unique(probs$mu_pct),
                                                       ks_pct =  unique(probs$ks_pct),
                                                       nbugs=unique(probs$nbugs),
                                                       spacing=unique(probs$spacing)))

combined.main <-left_join(x=probs,y=pg.main) %>%
  mutate(pred_error = prob_thrive-fit)

combined.main %>% mutate(sq_error = pred_error*pred_error)%>% ungroup() %>%
  summarise(RMSE = sqrt(mean(sq_error)))
summary(gam.no_i.all_smooth)

# Calculate per-crowding RMSEs ---------------

gam_rmse<-combined %>% mutate(sq_error = pred_error*pred_error) %>%
  group_by(nbugs,spacing) %>%
  summarise(RMSE = sqrt(mean(sq_error))) %>%
  mutate(Model="GAM")

mlr_rmse<-mlr_res %>% mutate(sq_error = resid*resid) %>%
  group_by(nbugs,spacing) %>%
  summarise(RMSE = sqrt(mean(sq_error))) %>%
  mutate(Model="MLR")

# also want to mention total RMSE in text
mlr_res %>% mutate(sq_error = resid*resid) %>% ungroup() %>%
  summarise(RMSE = sqrt(mean(sq_error)))

combined %>% mutate(sq_error = pred_error*pred_error)%>% ungroup() %>%
  summarise(RMSE = sqrt(mean(sq_error)))

rmse<-rbind(mlr_rmse,gam_rmse)%>%ungroup() %>%
  mutate_if(is.numeric,
            signif,
            digits = 3)

p_gam_mains_err<-ggplot(data=combined.main, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined.main,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($K_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title.x = element_text(size = 5,color="black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=3.5,color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7,hjust = 0.5),
        legend.key.size=unit(4, 'mm'))+
  ggtitle("General Additive Model (GAM)")+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))
#p_gam_mains_err

# Create table and plot-table for final panel ---------------
d <- tibble(x = c(0.95, 0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95), y = c(0.95, 0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95),
            nbugs = c(4,9,16,4,16,9,4,9,16),
            spacing= c(2.5,2.5,2.5,5,5,5,10,10,10),
            tb = list(rmse %>% filter(nbugs==4) %>% filter(spacing==2.5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==9) %>% filter(spacing==2.5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==16) %>% filter(spacing==2.5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==4) %>% filter(spacing==5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==9) %>% filter(spacing==5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==16) %>% filter(spacing==5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==4) %>% filter(spacing==10)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==9) %>% filter(spacing==10)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==9) %>% filter(spacing==10)%>%select(Model,RMSE)))

ptab<-ggplot(composite, aes(x=mu_pct*100, y=ks_pct*100)) +
  #geom_blank() +
  facet_grid(rows = vars(nbugs), cols=vars(spacing),
             labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))+
  geom_table(data=d,aes(x=x,y=y, label = tb),size=1.5,table.theme = ttheme_gtminimal)+
  xlab("")+ylab("")+
  theme(legend.position = "none",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 5,color="black"),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_markdown(size=5),
        strip.background= element_blank(),
        plot.title = element_text(size=7,hjust = 0.5))+
        ggtitle("RMSE comparison")

perrs<-p_err+p_gam_err+ptab+
  plot_layout(guides="collect")+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 8)) # nb could use this more to reduce code in above
#perrs
# ggsave(plot=perrs,here::here("output","model_errs.png"),units="in",width=8,height=4,dpi=330)

# Create final composition of all panels summarizing models and errors ---------------
p_mod_summ<-plot3/perrs+ plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 7),) # nb could use this more to reduce code in above

fname <- "model_summary.png"
floc <- here::here("output",fname)
ggsave(floc, p_mod_summ, width=8,height=5,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "model_summary.pdf"
floc <- here::here("output",fname)
ggsave(floc, p_mod_summ, width=8,height=5,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "model_summary.tiff"
floc <- here::here("output",fname)
ggsave(floc, p_mod_summ, width=8,height=5,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))
