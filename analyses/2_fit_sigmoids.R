library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(broom)
library(tidyr)
library(purrr)
library(RColorBrewer)

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

sigmoid <- function(df,xo_guess,k_guess){
  retval <- df %>%
    nest(data = -ks_pct) %>%
    mutate(
      fit = map(data, ~ nls( prob_thrive ~ 1/(1+exp(-k*(mu_pct-xo))),
                             data = .x,
                             start=list(xo=xo_guess,k=k_guess))),
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

# extracts p-values from lm
# https://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

fit_sigmoids <- function(df, N, spacing,xo_guess=0.2,k_guess=4){
  # assing nbugs and spacing with unique?
  a<-centers(sigmoid(df,xo_guess,k_guess)) %>% mutate(nbugs=N, spacing = spacing)
  fit_a <- lm(a$mu_50 ~ ks_pct, data = a)
  a$predicted_mu_50 = predict(fit_a)
  a$residual_mu_50 = residuals(fit_a)
  a$intercept <- fit_a$coefficients[1]
  a$slope_m50 <- fit_a$coefficients[2]
  a$a_r2_m50 <- summary(fit_a)$adj.r.squared
  a$p_value_m50r <- lmp(fit_a)

  b<-slopes(sigmoid(df,xo_guess,k_guess)) %>% mutate(nbugs=N, spacing = spacing)
  fit_b <- lm(b$steep ~ ks_pct, data = b)
  b$predicted_steep = predict(fit_b)
  b$residual_steep = residuals(fit_b)
  b$intercept_steep <- fit_b$coefficients[1]
  b$slope_steep <- fit_b$coefficients[2]
  b$a_r2_steep <- summary(fit_b)$adj.r.squared
  b$p_value_steep <- lmp(fit_b)

  return(inner_join(a,b,by=c("ks_pct","nbugs","spacing")))
}


bugs_4_2pt5<-(probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(2.5)))
a<-slopes(sigmoid(bugs_4_2pt5,0.2,4)) %>% mutate(nbugs=9, spacing =2.5)

a$log_k <- log(a$ks_pct+1)
a$k2 <- a$log_k^2
fit_a <- lm(a$steep ~ log_k + k2, data = a)
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(fit_a)
par(mfrow=c(1,1)) # Change back to 1 x 1

plot(predict(fit_a))
a<-fit_sigmoids((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(5))), 9, 5)
b<-fit_sigmoids((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(5))), 9, 5)
c<-fit_sigmoids((probs %>% filter(nbugs %in% c(9)) %>% filter(spacing %in% c(10))), 9, 10)


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
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +facet_grid(cols=vars(spacing))
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
