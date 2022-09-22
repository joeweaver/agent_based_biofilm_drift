library(readr)
library(dplyr)
library(here)

# precondition
# ./data/sweep_colony_outcomes should contain csvs describing simulation results
# if not, you may need to run 0_download_sim_results.R first

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("1_sigmoid_fitting.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','eda'))
create_dir_if_not_exist(here::here('output','eda','sigmoid_fits'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)

# each run of n bugs with s spacing has a number of affinities (ks_pct)
# this pulls out the data for a single n,s, affinity combination
select_sigmoid<-function(.data,n,s,affinity){
  .data %>%filter(nbugs %in% c(n)) %>%
    filter(spacing %in% c(s)) %>%
    filter(ks %in% c(affinity))
}

# the logistic/sigmoid function which we are fitting
logfun<-function(mu,xo,k){
  return(1/(1+exp(-k*(mu-xo))))
}

# perform a fit and calculate the SSE
# set up so that it can be used with optim
# We used to use nlm to fit these, but some gradient stability issues gave errors
sigmoid_fit_error<-function(pars,prob){
  mu_g<-pars[1]
  k_g<-pars[2]

  mu=prob$mu_pct
  error = prob$prob_thrive - logfun(mu,mu_g,k_g)
  return(sum(error*error))
}

# dataframe for recording fit data
sigmoid_fits <-data.frame(nbugs=c(),spacing=c(),ks=c(),mu_50=c(),steep=c(),
                          sse=c(),conv=c())

# fit all the sigmoids within a run for N bugs at s spacing (all ks_pct)
# returns a dataframe analagous to sigmoid_fits storing the fit data
# side effect, generates eda plots of all sigmoid fits
fit_sigmoids <- function(df,N,s){
  rec <-data.frame(nbugs=c(),spacing=c(),ks=c(),mu_50=c(),steep=c(),
                            sse=c(),conv=c())
  for(ks in unique(df$ks)){
    spacing<-s
    sig_run <-select_sigmoid(df,N,spacing,ks)
    mu_guess<-df$mu_pct[which.min(abs(df$prob_thrive-0.5))]
    dpt<-df$prob_thrive[which.min(abs(df$prob_thrive-0.5))]-df$prob_thrive[which.min(abs(df$prob_thrive-0.5))-1]
    dmu<-df$mu_pct[which.min(abs(df$prob_thrive-0.5))]-df$mu_pct[which.min(abs(df$prob_thrive-0.5))-1]
    k_guess=dpt/dmu*5
    m2=optim(par=c(mu_guess,k_guess), sigmoid_fit_error, method="L-BFGS-B",lower= c(-0.5,1.5), upper=c(1,100),prob=sig_run)
    if(m2$convergence != 0){
      message(glue::glue("fit did not converge n:{N} spacing:{spacing} ks:{ks} guess: {mu_guess} {k_guess}"))
    }

    rec <- rbind(rec,
      data.frame(nbugs=N,spacing=spacing,ks=ks,mu_50=m2$par[1],
                               steep=m2$par[2],sse=m2$value,conv=m2$convergence))
    mus=seq(-0.5,0.5,length.out=100)#c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
    fname<-glue::glue("sigmoid_fit_{sqrt(N)}x{sqrt(N)}_{spacing}_{ks}.png")
    png(here::here("output","eda",'sigmoid_fits',fname), width=800,height=800,units="px")
    title<-glue::glue("sigmoid_fit_{sqrt(N)}x{sqrt(N)}_{spacing}_{ks}")
    subt<-glue::glue("mu_50: {m2$par[1]} steepness:{m2$par[2]} SSE:{m2$value} conv:{m2$convergence}")
    plot(x=sig_run$mu_pct,y=sig_run$prob_thrive,type="p",pch=16,col="black",
         xlab="% mu change",ylab="Thrive probability",main=title,sub=subt)
    lines(mus,logfun(mus,m2$par[1],m2$par[2]),col='blue')
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



