library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(broom)
library(tidyr)
library(purrrwhich )
library(RColorBrewer)
library(ggpubr)
library(ggtext)

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
start_logging("3_spread_analysis.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','eda'))
create_dir_if_not_exist(here::here('output','eda','spread_analysis'))

# read the values of the fitted sigmoids
sig_fits <- read_csv(here::here("output","sigmoid_fits.csv"))

# add 1 to avoid 0 values which can blow up the optim call
delta_logistic<-function(mu_pct,delta,mu_50,steep){
  return(1+abs(delta-logistic_fun(mu_pct, mu_50, steep)))
}

sig_spread<-function(ub,lb,df){
  hi=optim(par=c(df$mu_50), delta_logistic, method="Brent",
           lower= c(-1.5), upper=c(1.5),
           delta=ub,mu_50=df$mu_50,steep=df$steep)
  if(hi$convergence != 0){
    message(glue::glue("HI did not converge {hi$convergence}"))
  }
  lo=optim(par=c(df$mu_50), delta_logistic, method="Brent",
          lower= c(-1.5), upper=c(1.5),
          delta=lb,mu_50=df$mu_50,steep=df$steep)
  if(hi$convergence != 0){
    message("LO did not converge")

  }
  return(data.frame(hi=hi$par[1],lo=lo$par[1]))
}

select_sig_fit<-function(.data,N,s,k_s){
  sig_fit <- .data %>% filter(nbugs %in% c(N)) %>%
    filter(spacing %in% c(s)) %>%
    filter(ks %in% c(k_s))
}

get_spreads<-function(.data,N,s,range_pct){
  lpct<-((1-range_pct)/2)
  upct<-1-lpct
  res<-data.frame()
  for(ks in unique(.data$ks)){
    sf<-select_sig_fit(.data,N,s,ks)
    spread<-sig_spread(upct,lpct,sf)
    #print(glue::glue("ks: {ks} {spread$hi} {spread$lo} {spread$hi-spread$lo} "))
    res<-rbind(res,cbind(sf,spread) %>% mutate(range_pct=range_pct, spread=hi-lo))


  }
  return(res)
}

all_spreads_2d<-rbind(get_spreads(sig_fits,4,2.5,0.95),
                    get_spreads(sig_fits,4,5,0.95),
                    get_spreads(sig_fits,4,10,0.95),

                    get_spreads(sig_fits,9,2.5,0.95),
                    get_spreads(sig_fits,9,5,0.95),
                    get_spreads(sig_fits,9,10,0.95),

                    get_spreads(sig_fits,16,2.5,0.95),
                    get_spreads(sig_fits,16,5,0.95),
                    get_spreads(sig_fits,16,10,0.95),

                    # TODO uncomment when runs complete
                    #get_spreads(sig_fits,25,2.5,0.90)
                    get_spreads(sig_fits,25,5,0.95))
                    #get_spreads(sig_fits,25,10,0.90)

all_spreads_1sd<-rbind(get_spreads(sig_fits,4,2.5,0.68),
                   get_spreads(sig_fits,4,5,0.68),
                   get_spreads(sig_fits,4,10,0.68),

                   get_spreads(sig_fits,9,2.5,0.68),
                   get_spreads(sig_fits,9,5,0.68),
                   get_spreads(sig_fits,9,10,0.68),

                   get_spreads(sig_fits,16,2.5,0.68),
                   get_spreads(sig_fits,16,5,0.68),
                   get_spreads(sig_fits,16,10,0.68),

                   # TODO uncomment when runs complete
                   #get_spreads(sig_fits,25,2.5,0.90)
                   get_spreads(sig_fits,25,5,0.68))

all_spreads <- rbind(all_spreads_1sd,all_spreads_2d) %>% mutate(ks_pct = ks/base_ks -1)

ggplot(all_spreads %>% filter(range_pct %in% c(0.95)),aes(x=ks_pct,y=spread))+geom_point()+
  geom_smooth(method="lm")+
  facet_grid(cols=vars(spacing),rows=vars(nbugs))+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")

ss<-all_spreads %>% filter(range_pct %in% c(0.95)) %>% filter(nbugs == 16) %>%filter(spacing == 10)
lmf <- lm(spread~ks_pct,ss)
lmf
summary(lmf)
plot(ss$ks_pct,ss$spread)
plot(lmf)

ggplot(all_spreads %>% filter(range_pct %in% c(0.95)) %>% filter(nbugs %in% c(9)),aes(x=spacing,y=mu_50))+geom_point()+
  geom_smooth(method="lm")+
  facet_grid(cols=vars(nbugs),rows=vars(ks_pct))+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")

ggplot(all_spreads %>% filter(range_pct %in% c(0.95)),aes(x=ks_pct,y=spread,color=factor(range_pct)))+geom_point()+
  geom_smooth(method="lm")+
  facet_grid(cols=vars(spacing),rows=vars(nbugs))+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")

write_csv(all_spreads,here::here("output","spread_fits.csv"))

spacing_label <- function(string) {
  glue::glue("<span style = 'color:#000000;'>{string}<span> <span style = 'color:#585858;'>diameter spacing<span>")
}
nbugs_label <- function(string) {
  glue::glue("<span style = 'color:#000000;'>{string}<span> <span style = 'color:#585858;'>intial bacteria<span>")
}


# spacing becomes more important at greater initial populations
poplabs <- c("Initial Population: 4", "Initial Population: 9", "Initial Population: 16","Initial Population: 25")
names(poplabs) <- c("4","9","16","25")
p <- ggplot(all_spreads%>%filter(range_pct==0.95), aes(x=ks_pct,y=spread,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
              #  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="black",
              #               linetype="solid",size=1.2,alpha=0.5) +
              geom_smooth(method='lm', formula= y~x+x*x, se=FALSE,size=0.4)+
              geom_point(size=1.1,alpha=0.9)+
              stat_poly_eq(formula=formula,
                           coef.digits = 3,
                           rr.digits = 2,
                           size=1.9,
                           eq.with.lhs = "italic(mu[~~50])~`=`~",
                           eq.x.rhs = "italic(K[s])",
                           label.y = "top", label.x = "left",
                           aes(label = paste(after_stat(eq.label),
                                             after_stat(rr.label), sep = "*\", \"*")))+
              #geom_smooth(method='lm', formula= y~x+x^4, se=FALSE,size=0.4)+
              ylab(TeX("Drift-Fitness co-dominance range ($spread_{95}$)")) +
              xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
              #ylim(-0.3,0.7)+
              coord_fixed()+
              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                 breaks = seq(-0.50,1.40,0.10))+
              scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                 breaks = seq(-0.50,0.50,0.10))+
              scale_color_brewer(palette = "Dark2")+
              guides(color=guide_legend(title="Spacing"),
                     shape=guide_legend(title="Spacing"),
                     linetype=guide_legend(title="Spacing")) +
              theme(legend.position ="top",#c(0.5,0.925),
                    legend.direction = "horizontal",
                    legend.title = element_text(size = 6),
                    legend.text = element_text(size = 6),
                    axis.title.x = element_text(size = 8,color="black"),
                    axis.title.y = element_text(size = 7,color="black"),
                    axis.text = element_text(size=6,color="black"),
                    axis.text.x = element_text(angle=-45),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    strip.text = element_text(size=8,color="black"),
                    strip.background = element_blank())+
              facet_grid(cols=vars(nbugs), labeller = labeller(nbugs = poplabs))

            ggsave(here::here('output','spread95_trend.tiff'),width=8,height=3.5, units="in",
                   dpi=330)
            ggsave(here::here('output','spread95_trend.png'),width=8,height=3.5, units="in",
                   dpi=330)
            ggsave(here::here('output','spread95_trend.pdf'),width=8,height=3.5, units="in",
                   dpi=330)

            # spacing becomes more important at greater initial populations
            poplabs <- c("Initial Population: 4", "Initial Population: 9", "Initial Population: 16","Initial Population: 25")
            names(poplabs) <- c("4","9","16","25")
            p <- ggplot(all_spreads%>%filter(range_pct==0.68), aes(x=ks_pct,y=spread,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
              #  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="black",
              #               linetype="solid",size=1.2,alpha=0.5) +
              geom_smooth(method='lm', formula= y~x+x*x, se=FALSE,size=0.4)+
              geom_point(size=1.1,alpha=0.9)+
              stat_poly_eq(formula=formula,
                           coef.digits = 3,
                           rr.digits = 2,
                           size=1.9,
                           eq.with.lhs = "italic(mu[~~50])~`=`~",
                           eq.x.rhs = "italic(K[s])",
                           label.y = "top", label.x = "left",
                           aes(label = paste(after_stat(eq.label),
                                             after_stat(rr.label), sep = "*\", \"*")))+
              #geom_smooth(method='lm', formula= y~x+x^4, se=FALSE,size=0.4)+
              ylab(TeX("Sigmoid midpoint ($\\μ_{50}$)")) +
              xlab(TeX("Change in substrate affnity ($\\k_{s}$)")) +
              #ylim(-0.3,0.7)+
              coord_fixed()+
              scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                 breaks = seq(-0.50,1.40,0.10))+
              scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                 breaks = seq(-0.50,0.50,0.10))+
              scale_color_brewer(palette = "Dark2")+
              guides(color=guide_legend(title="Spacing"),
                     shape=guide_legend(title="Spacing"),
                     linetype=guide_legend(title="Spacing")) +
              theme(legend.position = "top",#c(0.5,0.85),
                    legend.direction = "horizontal",
                    legend.title = element_text(size = 6),
                    legend.text = element_text(size = 6),
                    axis.title = element_text(size = 9,color="black"),
                    axis.text = element_text(size=6,color="black"),
                    axis.text.x = element_text(angle=-45),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    strip.text = element_text(size=8,color="black"),
                    strip.background = element_blank())+
              facet_grid(cols=vars(nbugs),rows=vars(range_pct), labeller = labeller(nbugs = poplabs))

            ggsave(here::here('output','spread68_trend.tiff'),width=8,height=3.5, units="in",
                   dpi=330)
            ggsave(here::here('output','spread68_trend.png'),width=8,height=3.5, units="in",
                   dpi=330)
            ggsave(here::here('output','spread68_trend.pdf'),width=8,height=3.5, units="in",
                   dpi=330)

# TODO break this out into its own
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))
probs <- thrive_probabilities(sim_results) %>% mutate(likely=prob_thrive/(1-prob_thrive+1e-6)+1e-6,
                                                   loglike=log(likely))

p<-ggplot(data=all_spreads, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=probs,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(prob_thrive,c(-0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
  scale_fill_manual(values = c("#d73027", "#f46d43", "#fdae61", "#fee090", '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#225ea8'),
                    labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  geom_segment(data=all_spreads %>% filter(range_pct %in% c(0.95)),
               aes(x=lo*100,xend=hi*100,y=ks_pct*100,yend=ks_pct*100),linetype="dotted",alpha=0.5)+
  geom_point(alpha=0.9)+
  geom_point(data=all_spreads %>% filter(range_pct %in% c(0.68)),
             aes(x=lo*100),shape=3,alpha=0.75)+
  geom_point(data=all_spreads %>% filter(range_pct %in% c(0.68)),
             aes(x=hi*100),shape=3,alpha=0.75)+
  facet_grid(rows = vars(nbugs), cols=vars(spacing), labeller = labeller(.rows = nbugs_label, .cols = spacing_label )) +
  ylab(TeX("Change in substrate affinity ($k_s$)")) +
  xlab(TeX("Change in maximum specific growth rate ($\\mu_{max}$)"))+
  coord_fixed(xlim=c(-50,50))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  guides(fill=guide_legend(title="Probability of becoming\na thriving colony\n")) +
  theme(legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 15),
        axis.text = element_text(size=8),
        axis.text.x = element_text(angle=-45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(color="black", fill="white", linetype="solid"),
        strip.text.x = element_markdown(size=10.5),
        strip.text.y = element_markdown(size=12))
ggsave(plot=p,here::here("output","prob_map.tiff"),width=8,height=9.5, units="in",
                     dpi=330)
ggsave(plot=p,here::here("output","prob_map.png"),width=8,height=9.5, units="in",
       dpi=330)
ggsave(plot=p,here::here("output","prob_map.pdf"),width=8,height=9.5, units="in",
       dpi=330)