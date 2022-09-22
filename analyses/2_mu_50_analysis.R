library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(broom)
library(tidyr)
library(purrr)
library(RColorBrewer)

# Checking linearity assumption for mu_50 and steepness

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first
#
# need to have output/sigmoidfits.csv
# if not, you have to run 1_sigmoid_fitting.R

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("2_mu_50_analysis.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','eda'))
create_dir_if_not_exist(here::here('output','eda','mu_50_fits'))

# read the values of the fitted sigmoids
sig_fits <- read_csv(here::here("output","sigmoid_fits.csv"))

# extracts p-values from lm
# https://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

linear_fit_summaries <- data.frame(b=c(),m=c(),m2=c(),
                                   b_p=c(),m_p=c(),m2_p=c(),
                                   ar2=c(),pval=c(),type=c(),
                                   N=c(),spacing=c())

# Checking linear fits ----------------------------------------------------
lmfits<-function(.data,N,s){
  b <- sig_fits %>% filter(nbugs %in% c(N)) %>% filter(spacing %in% c(s))
  b$ks_pct <- b$ks/base_ks-1
  lfit <- lm(mu_50~ks_pct, data=b)

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_lineplot.png")
  png(here::here("output","eda","mu_50_fits",fname), width=800,height=800,units="px")
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(x=b$ks_pct,y=b$mu_50,pch=16,col='blue',main=title)
  lines(x=b$ks_pct,predict(lfit,data=b))
  dev.off()

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_diag.png")
  png(here::here("output","eda","mu_50_fits",fname), width=800,height=800,units="px")
  par(mfrow=c(2,2))
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(lfit,main=title)
  dev.off()
  return(lfit)
}

fit2frame <-function(afit,type,N,spacing){
  if(type == 'lm'){
    return(data.frame(b=afit$coefficients[1],m=afit$coefficients[2],
                      m2=afit$coefficients[3],
                      b_p=summary(afit)$coefficients[1,4],
                      m_p=summary(afit)$coefficients[2,4],
                      m_2p=NA,
                      ar2=summary(afit)$adj.r.squared,pval=lmp(afit),
                      type=type,N=N,spacing=spacing))
  }
  if(type == 'lm_quad'){
    return(data.frame(b=afit$coefficients[1],m=afit$coefficients[2],
                      m2=afit$coefficients[3],
                      b_p=summary(afit)$coefficients[1,4],
                      m_p=summary(afit)$coefficients[2,4],
                      m_2p=summary(afit)$coefficients[3,4],
                      ar2=summary(afit)$adj.r.squared,pval=lmp(afit),
                      type=type,N=N,spacing=spacing))
  }
}

record_lmfit <-function(df,N,spacing){
  fit2frame(df %>%lmfits(N,spacing),type='lm',N,spacing)
}


linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,4,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,4,5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,4,10))

linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,9,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,9,5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,9,10))

linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,16,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,16,5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,16,10))

# TODO include when runs finished
#linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,25,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,25,5))
#linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,25,10))

# Checking quadratic term ----------------------------------------------------
lm_quad_fits<-function(.data,N,s){
  b <- sig_fits %>% filter(nbugs %in% c(N)) %>% filter(spacing %in% c(s))
  b$ks_pct <- b$ks/base_ks-1
  b$k2 <- b$ks_pct^2
  lfit <- lm(mu_50~ks_pct+k2, data=b)

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_lineplot_quad.png")
  png(here::here("output","eda","mu_50_fits",fname), width=800,height=800,units="px")
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(x=b$ks_pct,y=b$mu_50,pch=16,col='blue',main=title)
  pred <- predict(lfit)
  ix <- sort(b$ks_pct, index.return=TRUE)$ix
  lines(b$ks_pct[ix], pred[ix], col='red',lwd=2)
  dev.off()

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_diag_quad.png")
  png(here::here("output","eda","mu_50_fits",fname), width=800,height=800,units="px")
  par(mfrow=c(2,2))
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(lfit,main=title)
  dev.off()
  return(lfit)
}

record_quadfit <-function(df,N,spacing){
  fit2frame(df %>%lm_quad_fits(N,spacing),type='lm_quad',N,spacing)
}

linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,4,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,4,5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,4,10))

linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,9,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,9,5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,9,10))

linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,16,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,16,5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,16,10))

# TODO include when runs finished
#linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,25,2.5))
linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,25,5))
#linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,25,10))

# -------------------------------------------------------------------------
## checking to see if non linear issue is due to  using ks, mu pct
## this is interactive, so uncomment to run
# Calculate from observations a biggest loser's probabily of thriving
# probs <- thrive_probabilities(sim_results)
# a<-probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(2.5))%>% filter(ks %in% c(base_ks))
# plot(x=a$mu*1e5,y=a$prob_thrive)
# temp<- data.frame(ks=c(),xo=c())
# for (kval in unique(probs$ks)){
#   print(kval)
#   nls_pred<-nls( prob_thrive ~ 1/(1+exp(-k*(mu-xo))),
#                  data = probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(2.5))%>% filter(ks %in% c(kval)),
#                  start=list(xo=base_mu,k=1e5))
#   # sanity check the nls works
#   plot(predict(nls_pred))
#   readline(prompt="Press [enter] to continue")
#   temp <- rbind(temp, data.frame(ks=kval,xo=nls_pred$m$getPars()[1]))
# }
# # should still a parobola pattern in the residuals plot if it's not an artifact of pct, and we do
# lfit <- lm(xo~ks, data=temp)
# par(mfrow=c(2,2))
# plot(lfit)
# plot(predict(lfit))


to_plot <- sig_fits
to_plot$ks_pct <- to_plot$ks/base_ks-1
to_plot$k2 <- to_plot$ks_pct^2


p <- ggplot(to_plot %>% filter(spacing %in% c(2.5,5,10)), aes(x=ks_pct,y=mu_50,color=factor(nbugs),shape=factor(nbugs),linetype=factor(nbugs))) +
#  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="red",
#               linetype="solid",size=0.7,alpha=0.5) +
  geom_point(size=2.2)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab(TeX("Sigmoid midpoint ($\\mu_{50}$)")) +
  xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
  ylim(-0.3,0.7)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.30,0.60,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Dark2")+
  guides(color=guide_legend(title="Spacing"),
         shape=guide_legend(title="Spacing"),
         linetype=guide_legend(title="Spacing")) +
  theme(legend.position = c(0.1,0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18,color="black"),
        axis.text = element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(  color="black", fill="#FC4E07", size=1.5, linetype="solid"        ) )+
  facet_wrap(facets = vars(spacing),nrow=3,ncol=1)
p

ggsave("midpoints_1spacing.png",units="in",width=3.75,height=3*3.75,dpi=300)


p <- ggplot(to_plot %>% filter(spacing %in% c(2.5,5,10)), aes(x=ks_pct,y=mu_50,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
  #  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="red",
  #               linetype="solid",size=0.7,alpha=0.5) +
  geom_point(size=2.2)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab(TeX("Sigmoid midpoint ($\\mu_{50}$)")) +
  xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
  ylim(-0.3,0.7)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.30,0.60,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Dark2")+
  guides(color=guide_legend(title="Spacing"),
         shape=guide_legend(title="Spacing"),
         linetype=guide_legend(title="Spacing")) +
  theme(legend.position = c(0.1,0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18,color="black"),
        axis.text = element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(facets = vars(nbugs),nrow=4,ncol=1)
p

ggsave("midpoints_1nbugs.png",units="in",width=3.75,height=3*3.75,dpi=300)

p <- ggplot(rbind(a,b), aes(x=ks_pct,y=mu_50,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="black",
               linetype="solid",size=1.2,alpha=0.5) +
  geom_point(size=3.7)+
  geom_smooth(method='lm', formula= y~x, se=FALSE,size=1.4)+
  ylab(TeX("Sigmoid midpoint ($\\mu_{50}$)")) +
  xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
  ylim(-0.3,0.7)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.30,0.60,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Set2")+
  guides(color=guide_legend(title="Spacing"),
         shape=guide_legend(title="Spacing"),
         linetype=guide_legend(title="Spacing")) +
  theme(legend.position = c(0.1,0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18,color="black"),
        axis.text = element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave("midpoints_spacing_thikkk.png",units="in",width=8,height=8,dpi=300)

p <- ggplot(rbind(a,b), aes(x=ks_pct,y=steep,color=factor(spacing))) +
  geom_point()+
  geom_segment(aes(xend = ks_pct, yend = predicted_steep), color="red") +
  geom_smooth(method='lm', formula= y~x)
p

a<-fit_sigmoid((probs %>% filter(nbugs %in% c(4)) %>% filter(spacing %in% c(5))), 4, 5)
b<-fit_sigmoid((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(5))), 9, 5)
c<-fit_sigmoid((probs %>% filter(nbugs %in% c(16)) %>% filter(spacing %in% c(5))), 16, 5)
d<-fit_sigmoid((probs %>% filter(nbugs %in% c(25)) %>% filter(spacing %in% c(5))), 25, 5)

p <- ggplot(rbind(a,b,c,d), aes(x=ks_pct,y=mu_50,color=factor(nbugs),shape=factor(nbugs),linetype=factor(nbugs))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="black",
               linetype="solid",size=1.2,alpha=0.5) +
  geom_point(size=3.7)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, size=1.4)+
  ylab(TeX("Sigmoid midpoint ($\\mu_{50}$)")) +
  xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
  ylim(-0.3,0.9)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.30,0.90,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Set2")+
  scale_linetype_manual(values = c(1,4,6,5)) +
  scale_shape_manual(values = c(15,17,16,18)) +
  guides(color=guide_legend(title="Population Size"),
         shape=guide_legend(title="Population Size"),
         linetype=guide_legend(title="Population Size")) +
  theme(legend.position = c(0.17,0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18,color="black"),
        axis.text = element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave("midpoints_pop_thikkk.png",units="in",width=8,height=8,dpi=300)

p <- ggplot(rbind(a,b,c,d), aes(x=ks_pct,y=mu_50,color=factor(nbugs))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="red") +
    geom_point()+

  geom_smooth(method='lm', formula= y~x, se=FALSE)
p

p <- ggplot(rbind(a,b,c,d), aes(x=ks_pct,y=steep,color=factor(nbugs))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_steep), color="red") +
  geom_point()+

  geom_smooth(method='lm', formula= y~x, se=FALSE)
p

p <- ggplot(rbind(a,b,c,d), aes(x=ks_pct,y=steep,color=factor(nbugs))) +
  geom_line()

p



a<-fit_sigmoid((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(5))), 9, 5)
b<-fit_sigmoid((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(10))), 9, 10)
p <- ggplot(rbind(a,b), aes(x=ks_pct,y=steep,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_steep), color="black",
               linetype="solid",size=1.2,alpha=0.5) +
  geom_point(size=3.7)+
  geom_smooth(method='lm', formula= y~x, se=FALSE, spacing=1.4)+
  ylab(TeX("Steepness Factor ($S$)")) +
  xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
  #ylim(-0.3,0.7)+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1),
  #                   breaks = seq(-0.30,0.60,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Set2")+
  guides(color=guide_legend(title="Spacing"),
         shape=guide_legend(title="Spacing"),
         linetype=guide_legend(title="Spacing")) +
  theme(legend.position = c(0.71,0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18,color="black"),
        axis.text = element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave("steepness_spacing_thikkk.png",units="in",width=8,height=8,dpi=300)

a<-fit_sigmoid((probs %>% filter(nbugs %in% c(4)) %>% filter(spacing %in% c(5))), 4, 5)
b<-fit_sigmoid((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(5))), 9, 5)
c<-fit_sigmoid((probs %>% filter(nbugs %in% c(16)) %>% filter(spacing %in% c(5))), 16, 5)
d<-fit_sigmoid((probs %>% filter(nbugs %in% c(25)) %>% filter(spacing %in% c(5))), 25, 5)
p <- ggplot(rbind(a,b,c,d), aes(x=ks_pct,y=steep,color=factor(nbugs),shape=factor(nbugs),linetype=factor(nbugs))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_steep), color="black",
               linetype="solid",size=1.2,alpha=0.5) +
  geom_point(size=3.7)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
  ylab(TeX("Steepness Factor ($S$)")) +
  xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
  #ylim(-0.3,0.7)+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1),
  #                   breaks = seq(-0.30,0.60,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Set2")+
  guides(color=guide_legend(title="Population Size"),
         shape=guide_legend(title="Population Size"),
         linetype=guide_legend(title="Population Size")) +
  theme(legend.position = c(0.71,0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 18,color="black"),
        axis.text = element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave("steepness_pop_thikkk.png",units="in",width=8,height=8,dpi=300)
