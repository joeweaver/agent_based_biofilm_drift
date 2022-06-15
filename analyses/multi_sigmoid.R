library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(nls)
library(broom)
library(tidyr)
library(purrr)

# Fit sigmoid curves along lines of equal Ks for simulations

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("multi_sigmoid.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)



get_run <- function(df, N, spacing){
  return(probs %>% filter(N %in% c(N)) %>% filter(spacing %in% c(spacing)))
}

sigmoid <- function(df){
  retval <- df %>%
    nest(data = -ks_pct) %>%
    mutate(
      fit = map(data, ~ nls( prob_thrive ~ 1/(1+exp(-k*(mu_pct-xo))),
                             data = .x,
                             start=list(xo=0.2,k=4))),
      tidied = map(fit, tidy)
    ) %>%
    unnest(tidied)
  return(retval)
}

centers <- function(df){
  return(df %>%
           filter(term %in% c("xo")) %>%
           select(c("ks_pct", "estimate")) %>%
           rename(mu_50 = estimate))
}

slopes <- function(df){
  return(df %>%
           filter(term %in% c("k")) %>%
           select(c("ks_pct", "estimate")) %>%
           rename(steep = estimate))
}

library(dplyr)
fit_sigmoid <- function(df, N, spacing){
  a<-centers(sigmoid(df)) %>% mutate(nbugs=N, spacing = spacing)
  fit_a <- lm(a$mu_50 ~ ks_pct, data = a)
  a$predicted_mu_50 = predict(fit_a)
  a$residual_mu_50 = residuals(fit_a)

  b<-slopes(sigmoid(df))
  fit_b <- lm(b$steep ~ ks_pct, data = b)
  b$predicted_steep = predict(fit_b)
  b$residual_steep = residuals(fit_b)

  return(inner_join(a,b,by="ks_pct"))
}



a<-fit_sigmoid((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(5))), 9, 5)


library(RColorBrewer)
p <- ggplot(rbind(a), aes(x=ks_pct,y=mu_50,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="red",
               linetype="solid",size=0.7,alpha=0.5) +
  geom_point(size=2.2)+
  geom_smooth(method='lm', formula= y~x, se=FALSE)+
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
ggsave("midpoints_1spacing.png",units="in",width=8,height=8,dpi=300)

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
