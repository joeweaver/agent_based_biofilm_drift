# Select a multiple logistic regression model and attempt to predict simulation
# results.

library(readr) # handle csv file import
library(dplyr)  # work with tidy data
library(ggplot2) # figure generation
library(latex2exp)  # use LaTeX expressions for creating axis labels
library(here) # manage paths, keep @JennyBryan from incinerating our machine
library(broom) # helps fits work well with dataframes
library(tidyr) # used for pivots
library(RColorBrewer) # palette management in plots
library(ggpmisc) # add regression fits to plot
library(ggpubr) # help with plots
library(ggtext) # help with plots

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first
#
# need to have output/sigmoidfits.csv
# if not, you have to run 1_sigmoid_fitting.R

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("5_multiple_logistic_regression.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','si'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# conver probabilities to log likelihood

log_lik <- sim_results %>% filter(biggest_loser) %>%
  mutate(catnum = case_when(category == "Thriving" ~ 1,
                            TRUE ~ 0)) %>%
  group_by(ks,mu,yield,nbugs,spacing) %>%
  summarise(prob_thrive = mean(catnum))  %>%
  mutate(mu_pct = (mu - base_mu)/base_mu) %>%
  mutate(ks_pct = (ks - base_ks)/base_ks) %>%
  mutate(likely=prob_thrive/(1-prob_thrive+1e-6)+1e-6,
                        loglike=log(likely))


# Model selection and fitting ---------------------------------------------


# fit a ff linear model
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


# Make predictions from selected model ------------------------------------

preds<-exp(predict(fit))/(1+exp(predict(fit)))
log_lik$preds <- preds
log_lik$resid <- log_lik$prob_thrive - log_lik$preds

# calculate RMSE for entire model
log_lik %>% mutate(sq_error = resid*resid) %>% ungroup() %>%
  summarise(RMSE = sqrt(mean(sq_error)))


# Plot histogram of errors ------------------------------------------------

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

fname <- "mlr_pred_error_hist.png"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=4,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mlr_pred_error_hist.pdf"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=4,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mlr_pred_error_hist.tiff"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=4,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

# save the fitted model and results for comparison with GAM ---------------
write_rds(fit,here::here("output","mlr_fit.rds"))
log_info(paste('Wrote', file.path("output","mlr_fit.rds"), ' MD5Sum: ',
               md5sum(here::here("output","mlr_fit.rds"))))

write_csv(log_lik,here::here("output","mlr_predictions.csv"))
log_info(paste('Wrote', file.path("output","mlr_predictions.csv"), ' MD5Sum: ',
               md5sum(here::here("output","mlr_predictions.csv"))))

# Create a large-scale standalone plot for SI ---------------
p<-ggplot(data=log_lik, aes(x=mu_50*100, y=ks_pct*100)) +
  geom_raster(data=log_lik,aes(x=mu_pct*100, y=ks_pct*100,fill=cut(preds,c(-0.11,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1)))) +
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
        strip.text.y = element_markdown(size=10),
        strip.background= element_blank())+
  facet_grid(rows = vars(nbugs), cols=vars(spacing),labeller=labeller(.rows = nbugs_label, .cols = spacing_label ))

fname <- "mlr_solo_predictions.png"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mlr_solo_predictions.tiff"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mlr_solo_predictions.pdf"
floc <- here::here("output","si",fname)
ggsave(floc, p, width=8,height=8,units="in",dpi=330)
log_info(paste('Wrote', file.path("output","si",fname), ' MD5Sum: ',
               md5sum(floc)))
