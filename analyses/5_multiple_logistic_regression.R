library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)


library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(here)
library(broom)
library(tidyr)

library(RColorBrewer)
library(ggpubr)
library(ggtext)

#
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



# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))


log_lik <- sim_results %>% filter(biggest_loser) %>%
  mutate(catnum = case_when(category == "Thriving" ~ 1,
                            TRUE ~ 0)) %>%
  group_by(ks,mu,yield,nbugs,spacing) %>%
  summarise(prob_thrive = mean(catnum))  %>%
  mutate(mu_pct = (mu - base_mu)/base_mu) %>%
  mutate(ks_pct = (ks - base_ks)/base_ks) %>%
  mutate(likely=prob_thrive/(1-prob_thrive+1e-6)+1e-6,
                        loglike=log(likely))

# fit a ff
fit_ff <- lm(loglike ~ mu_pct + ks_pct + nbugs + spacing +
            mu_pct*ks_pct + mu_pct*nbugs +mu_pct*spacing +
            ks_pct*nbugs +ks_pct*spacing +
            nbugs * spacing +
            mu_pct*ks_pct*nbugs  +
            mu_pct*ks_pct*spacing +
            mu_pct*ks_pct*nbugs*spacing, data=log_lik)
summary(fit_ff)

# fit main-only
fit_mains <- lm(loglike ~ mu_pct + ks_pct + nbugs + spacing
               , data=log_lik)
summary(fit_mains)

# results backwards step removal - highester order interactions with words p-value
fit_bsr <- lm(loglike ~ mu_pct + ks_pct + nbugs + spacing +
                mu_pct*spacing +
                ks_pct*spacing , data=log_lik)
summary(fit_bsr)

fit<-fit_bsr

preds<-exp(predict(fit))/(1+exp(predict(fit)))
log_lik$preds <- preds
log_lik$resid <- log_lik$prob_thrive - log_lik$preds

mlr_rmse<-log_lik %>% mutate(sq_error = resid*resid) %>%
  group_by(nbugs,spacing) %>%
  summarise(RMSE = sqrt(mean(sq_error))) %>%
  mutate(Model="MLR")

log_lik %>% mutate(sq_error = resid*resid) %>% ungroup() %>%
  summarise(RMSE = sqrt(mean(sq_error)))

spacing_label <- function(string) {
  glue::glue("<span style = 'color:#000000;'>{string}<span> <span style = 'color:#585858;'>diameter spacing<span>")
}
nbugs_label <- function(string) {
  glue::glue("<span style = 'color:#000000;'>{string}<span> <span style = 'color:#585858;'>intial bacteria<span>")
}

poplabs <- c("Initial Population: 4", "Initial Population: 9", "Initial Population: 16")
names(poplabs) <- c("4","9","16")
write_csv(log_lik,here::here("output","mlr_predictions.csv"))

phist <- ggplot(data=log_lik,aes(x=resid)) +
  geom_histogram(binwidth=0.025,color="black",fill="skyblue")+
  ylab(TeX("Count")) +
  xlab(TeX("Error in Thrive Transition Chance (Prediction - Simulation)")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.60,0.90,0.10))+
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
ggsave(here::here("output","si","mlr_pred_error_hist.png"),units="in",width=8,height=8,dpi=330)

#save the fitted model and results for comparison with GAM
write_rds(fit,here::here("output","mlr_fit.rds"))
write_csv(log_lik,here::here("output","mlr_predictions.csv"))

p<-ggplot(data=log_lik, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=log_lik,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(preds,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
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
ggsave(here::here("output","si","mlr_solo_predictions.png"),units="in",width=8,height=8,dpi=330)