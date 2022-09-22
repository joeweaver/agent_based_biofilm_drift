library(readr)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)

# Fit sigmoid curves to each run. Mu_pct vs prob thrive for each combination of
#  N bugs, spacing, and ks. Write plots of each fit to an eda directory
#  and save fit data to a csv.

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("1_sigmoid_fitting.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','si'))
create_dir_if_not_exist(here::here('output','eda'))
create_dir_if_not_exist(here::here('output','eda','sigmoid_fits'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)

# each run of n bugs with s spacing has a number of affinities (ks)
# this pulls out the data for a single n,s, affinity combination
select_sigmoid<-function(.data,n,s,affinity){
  .data %>%filter(nbugs %in% c(n)) %>%
    filter(spacing %in% c(s)) %>%
    filter(ks %in% c(affinity))
}

# perform a fit and calculate the SSE
# set up so that it can be used with optim
# We used to use nlm to fit, but some gradient stability issues gave errors
sigmoid_fit_error<-function(pars,prob){
  mu_g<-pars[1]
  k_g<-pars[2]

  mu=prob$mu_pct
  error = prob$prob_thrive - logistic_fun(mu,mu_g,k_g)
  return(sum(error*error))
}

# dataframe for recording fit data
sigmoid_fits <-data.frame(nbugs=c(),spacing=c(),ks=c(),mu_50=c(),steep=c(),
                          sse=c(),conv=c())

# fit all the sigmoids within a run for N bugs at s spacing (all ks_pct)
# returns a dataframe analagous to sigmoid_fits storing the fit data
# side effect, generates eda plots of all sigmoid fits
fit_sigmoids <- function(df,N,s){
  # storage dataframe
  rec <-data.frame(nbugs=c(),spacing=c(),ks=c(),mu_50=c(),steep=c(),
                            sse=c(),conv=c())
  # one sigmoid per ks
  for(ks in unique(df$ks)){
    spacing<-s
    sig_run <-select_sigmoid(df,N,spacing,ks)
    # generate initial guesses for mu_50 and steepness
    # mu_50 is closest mu_pct to 0.50 probability
    mu_guess<-df$mu_pct[which.min(abs(df$prob_thrive-0.5))]
    # k is trickier, right now it's the slope of the line connecting the mu_50
    # guess to the one right before it multiplied by 5 (empirically this works)
    dpt<-df$prob_thrive[which.min(abs(df$prob_thrive-0.5))]-df$prob_thrive[which.min(abs(df$prob_thrive-0.5))-1]
    dmu<-df$mu_pct[which.min(abs(df$prob_thrive-0.5))]-df$mu_pct[which.min(abs(df$prob_thrive-0.5))-1]
    k_guess=dpt/dmu*5
    # fit the simgoid using optim on sigmoid_fit_error
    m2=optim(par=c(mu_guess,k_guess), sigmoid_fit_error, method="L-BFGS-B",lower= c(-0.5,1.5), upper=c(1,100),prob=sig_run)
    # warn us if any of the fits don't converge
    if(m2$convergence != 0){
      message(glue::glue("fit did not converge n:{N} spacing:{spacing} ks:{ks} guess: {mu_guess} {k_guess}"))
    }
    # save the relevant fit information
    rec <- rbind(rec,
      data.frame(nbugs=N,spacing=spacing,ks=ks,mu_50=m2$par[1],
                               steep=m2$par[2],sse=m2$value,conv=m2$convergence))
    # make a nice little plot for EDA and actually visualizing the fits
    mus=seq(-0.5,0.5,length.out=100)
    fname<-glue::glue("sigmoid_fit_{sqrt(N)}x{sqrt(N)}_{spacing}_{ks}.png")
    png(here::here("output","eda",'sigmoid_fits',fname), width=800,height=800,units="px")
    title<-glue::glue("sigmoid_fit_{sqrt(N)}x{sqrt(N)}_{spacing}_{ks}")
    subt<-glue::glue("mu_50: {m2$par[1]} steepness:{m2$par[2]} SSE:{m2$value} conv:{m2$convergence}")
    plot(x=sig_run$mu_pct,y=sig_run$prob_thrive,type="p",pch=16,col="black",
         xlab="% mu change",ylab="Thrive probability",main=title,sub=subt)
    lines(mus,logistic_fun(mus,m2$par[1],m2$par[2]),col='blue')
    dev.off()
  }
  return(rec)
}

# TODO enable fitting on 5x5_2.5 and 5x5_10 when runs completed
# fit sigmoids for all runs
sigmoid_fits<-rbind(fit_sigmoids(probs,4,2.5),
      fit_sigmoids(probs,4,5),
      fit_sigmoids(probs,4,10),
      fit_sigmoids(probs,9,2.5),
      fit_sigmoids(probs,9,5),
      fit_sigmoids(probs,9,10),
      fit_sigmoids(probs,16,2.5),
      fit_sigmoids(probs,16,5),
      fit_sigmoids(probs,16,10),
      #fit_sigmoids(probs,25,2.5),
      fit_sigmoids(probs,25,5))
      #fit_sigmoids(probs,25,10))

# save the fits to a csv file
write_csv(sigmoid_fits,here::here("output","sigmoid_fits.csv"))

# creating some SI plots

per_run_sigmoid_plots <- function(sf,pthrive, N, s) {
  df <- sf %>% filter(nbugs %in% c(N))%>%
    filter(spacing %in% c(s))
  sig_preds <- data.frame(mu_pct=c(),prob_thrive=c(),ks=c())

  for(k_s in df$ks){
    sig_fit <- df %>% filter(ks %in% c(k_s))
    mus=seq(-0.5,0.5,length.out=100)
    preds<-logistic_fun(mus,sig_fit$mu_50,sig_fit$steep)
    sig_preds<-rbind(sig_preds,data.frame(mu_pct=mus,prob_thrive=preds,ks=k_s))
  }

  sig_run <-pthrive %>% filter(nbugs %in% c(N))%>%
    filter(spacing %in% c(s))

  p<-ggplot(sig_run,aes(x=mu_pct,y=prob_thrive,group=ks,color=factor(ks)))+
    geom_point(size=1.2)+geom_line(data=sig_preds,aes(linetype=factor(ks)))+
    xlab("% change to baseline mu_max")+
    ylab("Probability of transitioning to thriving")+
    guides(color=guide_legend(title="Affinity (Ks)"),
           linetype=guide_legend(title="Affinity (Ks)")) +
    theme(panel.background = element_blank(),
          axis.line=element_line(),
          axis.text = element_text(size=6),
          axis.title = element_text(size=6),
          legend.text = element_text(size=4),
          legend.title = element_text(size=6))

  ggsave(here::here("output","si",glue::glue("sigmoid_fits_{sqrt(N)}x{sqrt(N)}_{s}.tiff")),p,width=3.75,height=3.75,units="in",dpi=300)
}

per_run_sigmoid_plots(sigmoid_fits,probs,4,2.5)
per_run_sigmoid_plots(sigmoid_fits,probs,4,5)
per_run_sigmoid_plots(sigmoid_fits,probs,4,10)

per_run_sigmoid_plots(sigmoid_fits,probs,9,2.5)
per_run_sigmoid_plots(sigmoid_fits,probs,9,5)
per_run_sigmoid_plots(sigmoid_fits,probs,9,10)

per_run_sigmoid_plots(sigmoid_fits,probs,16,2.5)
per_run_sigmoid_plots(sigmoid_fits,probs,16,5)
per_run_sigmoid_plots(sigmoid_fits,probs,16,10)

# TODO include runs when complete
#per_run_sigmoid_plots(sigmoid_fits,probs,25,2.5)
per_run_sigmoid_plots(sigmoid_fits,probs,25,5)
#per_run_sigmoid_plots(sigmoid_fits,probs,25,10)