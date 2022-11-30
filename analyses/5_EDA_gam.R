library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(broom)
library(tidyr)
library(purrr)
library(RColorBrewer)
library(itsadug)
library(mgcv)
library(tidymv)
library(RColorBrewer)
library(ggpubr)
library(ggtext)
library(cowplot)
library(viridis)
library(ggpp)
# Fit sigmoid curves along lines of equal Ks for simulations

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("4_EDA_gam.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)
probs$nbugs <- as.double(probs$nbugs)#factor(probs$nbugs)
probs$spacing <- as.double(probs$spacing)# factor(probs$spacing)

# quick scrren of models (% deviation and AIC) --------
# # no-interact, all smooth
# gam.no_i.all_smooth <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+s(spacing,k=3),
#                            data = probs,method = "REML")
# summary(gam.no_i.all_smooth )
# summary(gam.no_i.all_smooth )$p.table
# summary(gam.no_i.all_smooth )$s.table
# AIC(gam.no_i.all_smooth )
# k.check(gam.no_i.all_smooth)
#
# # no-interact, all linear
# gam.no_i.all_linear <- gam(prob_thrive ~ mu_pct + ks_pct+nbugs+spacing, data = probs,method = "REML")
# summary(gam.no_i.all_linear )
# summary(gam.no_i.all_linear )$p.table
# summary(gam.no_i.all_linear )$s.table
# AIC(gam.no_i.all_linear )
# k.check(gam.no_i.all_linear)
#
# # no-interact, linear crowd
# gam.no_i.linear_crowd <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+nbugs+spacing, data = probs,method = "REML")
# summary(gam.no_i.linear_crowd )
# summary(gam.no_i.linear_crowd )$p.table
# summary(gam.no_i.linear_crowd )$s.table
# AIC(gam.no_i.linear_crowd )
# k.check(gam.no_i.linear_crowd)
#
# # all 2way interact, all smooth (except spacing in case of nbugsxspacing)
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
#
# # all 2way interact, linear crowd
# gam.2way_i.linear_crowd <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+nbugs+spacing+
#                                s(mu_pct,ks_pct) + s(mu_pct,by=nbugs) + s(mu_pct,by=spacing) +
#                                s(ks_pct,by=nbugs) + s(ks_pct,by=spacing) +
#                                nbugs*spacing, data = probs, method = "REML")
# summary(gam.2way_i.linear_crowd )
# summary(gam.2way_i.linear_crowd )$p.table
# summary(gam.2way_i.linear_crowd )$s.table
# AIC(gam.2way_i.linear_crowd )
# k.check(gam.2way_i.linear_crowd)
#
# # w 3way interact, all smooth  (except spacing in case of nbugsxspacing)
# gam.3way_i.all_smooth  <- gam(prob_thrive ~ s(mu_pct) + s(ks_pct)+s(nbugs,k=3)+s(spacing,k=3)+
#                                  s(mu_pct,ks_pct,k=45) + s(mu_pct,nbugs) + s(mu_pct,spacing) +
#                                  s(ks_pct,nbugs) + s(ks_pct,spacing) +
#                                  s(nbugs,by=spacing,k=3)+
#                                  + te(mu_pct, ks_pct,spacing,d=c(1,2)) +
#                                  + te(mu_pct, ks_pct,nbugs,d=c(1,2))
#                                  , data = probs, method = "REML")
# summary(gam.3way_i.all_smooth  )
# summary(gam.3way_i.all_smooth  )$p.table
# summary(gam.3way_i.all_smooth  )$s.table
# AIC(gam.3way_i.all_smooth  )
# k.check(gam.3way_i.all_smooth )
#
#
#
# # w 3way interact, linear crowd
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

# remove low p w 3 way interact, all smooth  (except spacing in case of nbugsxspacing)
gam.bsr.3way_i.all_smooth  <- gam(prob_thrive ~ spacing +
                                    s(mu_pct,ks_pct,k=45) + s(mu_pct,spacing) +
                                    s(ks_pct,spacing) +
                                    + te(mu_pct, ks_pct,spacing,d=c(1,2)) +
                                    + te(mu_pct, ks_pct,nbugs,d=c(1,2))
                                  , data = probs, method = "REML")
summary(gam.bsr.3way_i.all_smooth  )
summary(gam.bsr.3way_i.all_smooth  )$p.table
summary(gam.bsr.3way_i.all_smooth  )$s.table
AIC(gam.bsr.3way_i.all_smooth  )
k.check(gam.bsr.3way_i.all_smooth )


#plot(linear_model)
#vis.gam(linear_model, theta = 40, n.grid = 50, lwd = 0.4)

pg<-predict_gam(gam.bsr.3way_i.all_smooth,values=list(mu_pct = unique(probs$mu_pct),
                                         ks_pct =  unique(probs$ks_pct),
                                 nbugs=unique(probs$nbugs),
                                 spacing=unique(probs$spacing)))

combined <-left_join(x=probs,y=pg) %>%
  mutate(pred_error = prob_thrive-fit)


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
ggsave(here::here("output","si","gam_pred_error_hist.png"),units="in",width=8,height=8,dpi=330)


spacing_label <- function(string) {
  glue::glue("<span style = 'color:#585858;'>Spacing:&nbsp; <span><span style = 'color:#000000;'> {string}<span> ")
}
nbugs_label <- function(string) {
  glue::glue("<span style = 'color:#585858;'>Ini. Pop.:&nbsp;<span><span style = 'color:#000000;'> {string}<span> ")
}


# spacing becomes more important at greater initial populations
poplabs <- c("Initial Population: 4", "Initial Population: 9", "Initial Population: 16","Initial Population: 25")
names(poplabs) <- c("4","9","16","25")
p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
                facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{max}$)"))+
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
p
ggsave(here::here("output","si","gam_solo_predictions.png"),units="in",width=8,height=8,dpi=330)


p_comp_gam_pred<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{~~max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  coord_fixed(ylim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-30,-10,10,30,50)) +
  guides(fill=guide_legend(nrow=1)) +

  theme(legend.position = "top",
        legend.direction = "horizontal",

        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title.x = element_text(size = 7,color="black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=5,color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=6),
        strip.text.y = element_blank())+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))
p_comp_gam_pred
legend <- cowplot::get_legend(p_comp_gam_pred)

# quick c&p from 5_multi
mlr_res <- read_csv(here::here("output","mlr_predictions.csv"))
composite <- left_join(combined, mlr_res) %>%
  mutate(cpred = (fit+fit+preds)/3,
         cerror = prob_thrive - cpred)


p<-ggplot(data=mlr_res, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=mlr_res,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(preds,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{max}$)"))+
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
        axis.text = element_text(size=4,color="black"),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7,hjust = 0.5))+
  ggtitle("Multiple Linear Regression (MLR)")+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))
#p

p_comp_gam_pred<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(fit,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","100-110%"))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
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
        axis.text = element_text(size=4,color="black"),
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
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
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
        axis.text = element_text(size=4,color="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_markdown(size=5),
        strip.text.y = element_markdown(size=5),
        strip.background = element_blank(),
        plot.title = element_text(size=7,hjust = 0.5))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))+
  ggtitle("Simulation 'Ground Truth'")
#p_comp_sim_pred

# p4<-plot_grid(p_comp_gam_pred,p_comp_gam_pred,p_comp_sim_pred,nrow=1,scale = c(1, 1, 1.1),axis="b",align="hv",
#               labels=c("A) MLR ","B) GAM"," C) Simulation 'Ground Truth'"),
#               label_size=8,
#               hjust=c(-0.5,-0.45,0.05))
# p5<-plot_grid(legend,p4,ncol=1,scale=c(0.5,1))
# p5
# ggsave(here::here("output","si","composite_predictions.png"),units="in",width=8,height=4.25,dpi=330)
#
# ggsave(here::here("output","si","composite_predictions.pdf"),units="in",width=8,height=4.25,dpi=330)

library(patchwork)

plot3 <- p + p_comp_gam_pred +p_comp_sim_pred +
  plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A')

ggsave(here::here("output","si","composite_predictions.png"),units="in",width=8,height=3.5,dpi=330)

ggsave(here::here("output","si","composite_predictions.pdf"),units="in",width=8,height=4,dpi=330)


p_err<-ggplot(data=mlr_res, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=mlr_res,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(resid,seq(-0.8,0.8,by=0.1)))) +
  scale_fill_manual(values = c("#424086FF",
                               "#3B528BFF", "#33638DFF", "#2C728EFF", "#26828EFF",
                               "#21908CFF", "#1F9F88FF", "#27AD81FF","#3EBC74FF",
                               "#5DC863FF", "#82D34DFF", "#AADC32FF", "#D5E21AFF",
                               "#FDE725FF"))+

  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{max}$)"))+
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
        axis.text = element_text(size=4,color="black"),
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
p_err

# "#440154FF" "#48186AFF" "#472D7BFF" "#424086FF" "#3B528BFF" "#33638DFF" "#2C728EFF" "#26828EFF" "#21908CFF" "#1F9F88FF" "#27AD81FF"
# "#3EBC74FF" "#5DC863FF" "#82D34DFF" "#AADC32FF" "#D5E21AFF" "#FDE725FF"
p_gam_err<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(pred_error,seq(-0.8,0.8,by=0.1)))) +
 scale_fill_manual(values = c("#33638DFF","#2C728EFF","#26828EFF","#21908CFF","#1F9F88FF","#27AD81FF","#3EBC74FF"))+
                    #labels = seq(-0.8,0.8,by=0.1))+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
  coord_fixed()+
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
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
        axis.title.x = element_text(size = 5,color="black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=4,color="black"),
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
p_gam_err

gam_rmse<-combined %>% mutate(sq_error = pred_error*pred_error) %>%
  group_by(nbugs,spacing) %>%
  summarise(RMSE = sqrt(mean(sq_error))) %>%
  mutate(Model="GAM")

mlr_rmse<-mlr_res %>% mutate(sq_error = resid*resid) %>%
  group_by(nbugs,spacing) %>%
  summarise(RMSE = sqrt(mean(sq_error))) %>%
  mutate(Model="MLR")

rmse<-rbind(mlr_rmse,gam_rmse)%>%ungroup() %>%
  mutate_if(is.numeric,
            signif,
            digits = 3)
rmse

d <- tibble(x = c(0.95, 0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95), y = c(0.95,0.95, 0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95),
            nbugs = c(4,9,16,25,4,16,9,4,9,16),
            spacing= c(2.5,2.5,2.5,5,5,5,5,10,10,10),
            tb = list(rmse %>% filter(nbugs==4) %>% filter(spacing==2.5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==9) %>% filter(spacing==2.5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==16) %>% filter(spacing==2.5)%>%select(Model,RMSE),
                      rmse %>% filter(nbugs==25) %>% filter(spacing==5)%>%select(Model,RMSE),
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
        legend.title = element_text(size = 8),
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
#ptab

perrs<-p_err+p_gam_err+ptab+
  plot_layout(guides="collect")+
  plot_annotation(tag_levels = 'A')
perrs
ggsave(plot=perrs,here::here("output","model_errs.png"),units="in",width=8,height=4,dpi=330)

p_mod_summ<-plot3/perrs+ plot_annotation(tag_levels = 'A')
ggsave(plot=p_mod_summ,here::here("output","model_summary.png"),units="in",width=8,height=6.5,dpi=330)
ggsave(plot=p_mod_summ,here::here("output","model_summary.tiff"),units="in",width=8,height=6.5,dpi=330)
ggsave(plot=p_mod_summ,here::here("output","model_summary.pdf"),units="in",width=8,height=6.5,dpi=330)


#
# mtcars %>% group_by(cyl) %>%
#   summarize(wt = mean(wt), mpg = mean(mpg)) %>%
#   ungroup() %>%
#   mutate(wt = sprintf("%.2f", wt),
#          mpg = sprintf("%.1f", mpg)) -> tb
#
#
# df <- tibble(x = 5.45, y = 34, tb = list(tb))
#
# ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
#   geom_point() +
#   geom_table(data = df,
#              aes(x = x, y = y, label = tb))
#
#
# mtcars %>% group_by(cyl) %>%
#   summarize(wt = mean(wt), mpg = mean(mpg)) %>%
#   ungroup() %>%
#   mutate(wt = sprintf("%.2f", wt),
#          mpg = sprintf("%.1f", mpg)) -> tb
#
# df <- tibble(x = 5.45, y = 34, cyl=cyl, tb = list(tb))
#
# tableA = tb %>% filter(cyl == 4)
# tableB = tb %>% filter(cyl == 6)
# tableC = tb %>% filter(cyl == 8)
# d <- tibble(x = c(0.95, 0.95,0.95), y = c(0.95, 0.95,0.95),
#             cyl = c(4,6,8), tb = list(tableA, tableB,tableC))
#
# ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
#   geom_point() +
#   geom_table(data = d,
#              aes(x = x, y = y, label = tb),table.theme = ttheme_gtbw, color="darkblue") +
#   facet_grid(rows=vars(cyl))
#
#
#
# ptab<-ggplot(composite, aes(x=mu_pct*100, y=ks_pct*100)) +
#   #geom_point() +
#   facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   geom_table(data=df,aes(x=x,y=y, label = tb))
# ptab
#
#
# mtcars %>%
#   group_by(cyl) %>%
#   summarize(wt = mean(wt), mpg = mean(mpg)) %>%
#   ungroup() %>%
#   mutate(wt = sprintf("%.2f", wt),
#          mpg = sprintf("%.1f", mpg)) -> tb
#
# df <- tibble(x=2, y = 2, tb = list(tb))
#
#
# gam_rmse<-combined %>% mutate(sq_error = pred_error*pred_error) %>%
#   group_by(nbugs,spacing) %>%
#   summarise(RMSE = sqrt(mean(sq_error))) %>%
#   mutate(Model="GAM")
#
# mlr_rmse<-mlr_res %>% mutate(sq_error = resid*resid) %>%
#   group_by(nbugs,spacing) %>%
#   summarise(RMSE = sqrt(mean(sq_error))) %>%
#   mutate(Model="MLR")
#
# rmse<-rbind(mlr_rmse,gam_rmse)%>%ungroup()
# rmse
#
# tableA = tb %>% filter(cyl == 4)
# tableB = tb %>% filter(cyl == 6)
# tableC = tb %>% filter(cyl == 8)
# d <- tibble(x = c(0.95, 0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95), y = c(0.95,0.95, 0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95),
#             nbugs = c(4,9,16,25,4,16,9,4,9,16),
#             spacing= c(2.5,2.5,2.5,5,5,5,5,10,10,10),
#             tb = list(rmse %>% filter(nbugs==4) %>% filter(spacing==2.5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==9) %>% filter(spacing==2.5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==16) %>% filter(spacing==2.5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==25) %>% filter(spacing==5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==4) %>% filter(spacing==5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==9) %>% filter(spacing==5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==16) %>% filter(spacing==5)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==4) %>% filter(spacing==10)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==9) %>% filter(spacing==10)%>%select(model,RMSE),
#                       rmse %>% filter(nbugs==9) %>% filter(spacing==10)%>%select(model,RMSE)))
#
# ptab<-ggplot(composite, aes(x=mu_pct*100, y=ks_pct*100)) +
#   #geom_blank() +
#   facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   geom_table(data=d,aes(x=x,y=y, label = tb))
# ptab
#
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(pred_error,c(-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25))))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values =
#                       c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8')
#          )+
#   coord_fixed()
# p
#
#
# p1<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=composite,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(composite$pred_error,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(composite$prob_thrive,breaks=10)))+
#   theme(legend.position="none")+
#   coord_fixed()
# p1
# p2<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=composite,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(composite$resid,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(composite$cpred,breaks=10)))+
#   theme(legend.position="none")+
#   coord_fixed()
# p2
# p3<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=composite,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(composite$cerror,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(composite$cpred,breaks=10)))+
#   theme(legend.position="none")+
#   coord_fixed()
# p1
#
#
#
# p1<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=composite,aes(x=mu_pct*100, y=ks_pct*100,
#                                  fill=cut(composite$prob_thrive/0.25,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(composite$prob_thrive/0.25,breaks=10)))+
# #  theme(legend.position="none")+
#   coord_fixed()
# p1
# p2<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=composite,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(composite$fit,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(composite$prob_thrive,breaks=10)))+
#   theme(legend.position="none")+
#   coord_fixed()
# p2
# p3<-ggplot(data=composite, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=composite,aes(x=mu_pct*100, y=ks_pct*100,fill=cut((composite$pred+composite$fit)/2,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(composite$prob_thrive,breaks=10)))+
#   theme(legend.position="none")+
#   coord_fixed()
# p1
# library(cowplot)
# p4<-plot_grid(p1,p2,p3,nrow=1)
# p4
#
# p2
# ggsave(plot=p4,"testcombined.png",width=8,height=9,units="in")
# data_plot <- ggplot(probs %>% filter(ks_pct == 0.0), aes(y = prob_thrive, x =mu_pct )) +
#   geom_point()+
#   theme_bw()+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing))
# data_plot
#
# fitted(gam1)
# #plot(gam1,all.terms = TRUE)
# #max(combined$se.fit)
# #min(combined$se.fit)
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(se.fit,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(combined$se.fit,breaks=10)))+
#   coord_fixed()
# #p
#
# combined <- combined %>%
#   mutate(l_95 = fit-se.fit*qnorm(0.975),
#          u_95 = fit-se.fit*qnorm(0.975),
#          se_rat = pred_error/se.fit,
#          in_95 = (prob_thrive>l_95),
#          in_2d = (pred_error< 2* 0.06271919),
#          err_cvish = pred_error/0.06271919)
#
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(se_rat,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(combined$se_rat,breaks=10)))+
#   coord_fixed()
# p
#
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=in_95))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   coord_fixed()
# p
#
# sqrt(summary(gam1)[["dispersion"]])
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=in_2d))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   coord_fixed()
# p
#
#
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(err_cvish,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(combined$err_cvish,breaks=10)))+
#   coord_fixed()
# p

# ### train on only n=4,16 and all spacings
# prb_sub <-probs%>%filter(nbugs %in% c(4,16))
# gam1 <- gam(prob_thrive ~ s(mu_pct,by=nbugs) +s(mu_pct,by=nbugs)+ s(ks_pct,by=nbugs)+ s(ks_pct,mu_pct,by=nbugs)+
#               s(mu_pct,by=spacing) +s(mu_pct,by=spacing)+ s(ks_pct,by=spacing)+ s(ks_pct,mu_pct,by=spacing) +nbugs+spacing, data = prb_sub,method = "REML")
#
#
#
# summary(gam1)
# AIC(gam1)
# #plot(linear_model)
# #vis.gam(linear_model, theta = 40, n.grid = 50, lwd = 0.4)
#
# pg<-predict_gam(gam1,values=list(mu_pct = unique(probs$mu_pct),
#                                  ks_pct =  unique(probs$ks_pct),
#                                  nbugs = unique(probs$nbugs)))
#
# combined <-left_join(x=probs,y=pg) %>%
#   mutate(pred_error = prob_thrive-fit)
#
# spacing_label <- function(string) {
#   glue::glue("<span style = 'color:#000000;'>{string}<span> <span style = 'color:#585858;'>diameter spacing<span>")
# }
# nbugs_label <- function(string) {
#   glue::glue("<span style = 'color:#000000;'>{string}<span> <span style = 'color:#585858;'>intial bacteria<span>")
# }
#
#
# p<-ggplot(data=combined, aes(x=mu_50*100, y=ks_pct*100)) +
#   geom_raster(data=combined,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(pred_error,breaks=10)))+
#   facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label ))+
#   scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
#                     labels = levels(cut(combined$pred_error,breaks=10)))+
#   coord_fixed()
# p



