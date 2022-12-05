# Analyze how mu_50 varies between conditions

library(readr) # handle csv file import
library(dplyr)  # work with tidy data
library(ggplot2) # figure generation
library(latex2exp)  # use LaTeX expressions for creating axis labels
library(here) # manage paths, keep @JennyBryan from incinerating our machine
library(broom) # helps fits work well with dataframes
library(tidyr) # used for pivots
library(purrr) # helps with some functional programming
library(RColorBrewer) # palette management in plots
library(ggpmisc) # add regression fits to plot


# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first
#
# need to have output/sigmoidfits.csv
# if not, you have to run 1_sigmoid_fitting.R

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("3_mu_50_analysis.log")

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

  floc <- here::here("output","eda","mu_50_fits",fname)
  log_info(paste('Wrote', file.path("output","eda","mu_50_fits",fname), ' MD5Sum: ',
                 md5sum(floc)))

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_diag.png")
  png(here::here("output","eda","mu_50_fits",fname), width=800,height=800,units="px")
  par(mfrow=c(2,2))
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(lfit,main=title)
  dev.off()

  floc <- here::here("output","eda","mu_50_fits",fname)
  log_info(paste('Wrote', file.path("output","eda","mu_50_fits",fname), ' MD5Sum: ',
                 md5sum(floc)))

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

  floc <- here::here("output","eda","mu_50_fits",fname)
  log_info(paste('Wrote', file.path("output","eda","mu_50_fits",fname), ' MD5Sum: ',
                 md5sum(floc)))

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_diag_quad.png")
  png(here::here("output","eda","mu_50_fits",fname), width=800,height=800,units="px")
  par(mfrow=c(2,2))
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(lfit,main=title)
  dev.off()

  floc <- here::here("output","eda","mu_50_fits",fname)
  log_info(paste('Wrote', file.path("output","eda","mu_50_fits",fname), ' MD5Sum: ',
                 md5sum(floc)))

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

sig_fits$ks_pct <- sig_fits$ks/base_ks-1
sig_fits$k2 <- sig_fits$ks_pct^2


# Create and save plot ----------------------------------------------------
formula <- y ~ poly(x, 1, raw = TRUE)
# spacing becomes more important at greater initial populations
poplabs <- c("Initial Population: 4", "Initial Population: 9", "Initial Population: 16")
names(poplabs) <- c("4","9","16")
p <- ggplot(sig_fits, aes(x=ks_pct,y=mu_50,color=factor(spacing),shape=factor(spacing),linetype=factor(spacing))) +
#  geom_segment(aes(xend = ks_pct, yend = predicted_mu_50), color="black",
#               linetype="solid",size=1.2,alpha=0.5) +
  geom_smooth(method='lm', formula= formula, se=FALSE,size=0.4)+
  stat_poly_eq(formula=formula,
               coef.digits = 3,
               rr.digits = 2,
               size=2.3,
               eq.with.lhs = "italic(mu[~~50])~`=`~",
               eq.x.rhs = "italic(K[s])",
               label.y = "bottom", label.x = "right",
               aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))+
  geom_point(size=1.1,alpha=0.9)+
  #geom_smooth(method='lm', formula= y~x+x^4, se=FALSE,size=0.4)+
  ylab(TeX("50% Thriving Odds ($\\μ_{50}$)")) +
  xlab(TeX("Change in substrate affnity ($\\K_{s}$)")) +
  #ylim(-0.3,0.7)+
  coord_fixed()+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.90,0.10))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.50,0.50,0.10))+
  scale_color_brewer(palette = "Dark2")+
  guides(color=guide_legend(title="Spacing"),
         shape=guide_legend(title="Spacing"),
         linetype=guide_legend(title="Spacing")) +
  theme(legend.position = c(0.5,0.95),
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
  facet_grid(cols=vars(nbugs), labeller = labeller(nbugs = poplabs))

fname <- "mu_50_trend.tiff"
floc <- here::here("output",fname)
ggsave(floc, p, width=8,height=3.5,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mu_50_trend.png"
floc <- here::here("output",fname)
ggsave(floc, p, width=8,height=3.5,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))

fname <- "mu_50_trend.pdf"
floc <- here::here("output",fname)
ggsave(floc, p, width=8,height=3.5,units="in",dpi=330)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))


# checking to see if non linear issue is due to  using ks, mu pct ------------
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


