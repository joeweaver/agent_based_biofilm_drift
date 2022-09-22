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

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("2_EDA_linear_vs_quad.log")

# Create relevant output dirs ---------
create_dir_if_not_exist(here::here('output'))
create_dir_if_not_exist(here::here('output','diagnostic'))

# Read simulation results ---------
sim_results <- get_all_results(here::here("data","sweep_colony_outcomes"))

# Calculate from observations a biggest loser's probabily of thriving --------
probs <- thrive_probabilities(sim_results)

sigmoid <- function(.data,xo_guess,k_guess){
  retval <- .data %>%
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

sigmoid_abs <- function(.data,xo_guess,k_guess){
  retval <- .data %>%
    nest(data = -ks) %>%
    mutate(
      fit = map(data, ~ nls( prob_thrive ~ 1/(1+exp(-k*(mu-xo))),
                             data = .x,
                             start=list(xo=xo_guess,k=k_guess))),
      tidied = map(fit, tidy)
    ) %>%
    unnest(tidied)
  return(retval)
}

centers <- function(.data){
  return(.data %>%
           filter(term %in% c("xo")) %>%
           select(c("ks_pct", "estimate")) %>%
           rename(mu_50 = estimate))
}

slopes <- function(.data){
  return(.data %>%
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

# TODO make plots of all sigmoid fits with r2
a<-probs%>%filter(nbugs %in% c(4)) %>% filter(spacing %in% c(2.5))

b<-a%>%filter(ks %in% c(base_ks))
plot(x=b$mu_pct,y=b$prob_thrive)
for(q in range(-10,10)){
  k_g=(1*10^(q))
  print(k_g)
  nls( prob_thrive ~ 1/(1+exp(-k*(mu-xo))),
       data = b,
       start=list(xo=0.3,k=1e-9))
  }
nls(prob_thrive ~ 1/(1+exp(-k*(mu-xo))),
     data = b,
     start=list(xo=0.05,k=25))

nls(prob_thrive ~ SSlogis(b$mu_pct,1,0.2,1), data=b)

linear_fit_summaries <- data.frame(b=c(),m=c(),m2=c(),
                                   b_p=c(),m_p=c(),m2_p=c(),
                                   ar2=c(),pval=c(),type=c(),
                                   N=c(),spacing=c())
logfun<-function(mu,xo,k){
  return(1/(1+exp(-k*(mu-xo))))
}
mus=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
logfun(mus,0.2,1)

k=28.2197589
mug=0.0796757
plot(x=b$mu_pct,y=b$prob_thrive)
mus=seq(-0.5,0.5,length.out=100)
points(mus,logfun(mus,mug,k),col='blue',pch=16)
points(mus,logfun(mus,mug,k),col='red',pch=17)
error = b$prob_thrive - logfun(mus,mug,k)
mean(error*error)

a1<-function(mu_g,k_g){

  mu=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
  error = b$prob_thrive - logfun(mu,mu_g,k_g)
  return(mean(error*error))
}
a1(0.075,43)
a1vec<-function(pars){
  mu_g<-pars[1]
  k_g<-pars[2]

  mu=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
  error = b$prob_thrive - logfun(mu,mu_g,k_g)
  return(sum(error*error))
}
a1vec(c(0.11908122,13.54))
logfun(b$mu_pct,0.0796757,68.2197589)
logfun(b$mu_pct,0.11908122,13.54)
b$prob_thrive

m1=optim(par=c(0.075,43),a1vec, method="L-BFGS-B",lower= c(-0.5,0.5), upper=c(1,100))
m1


# Checking linear fits ----------------------------------------------------
lmfits<-function(.data,N,s){
  b <- .data %>% filter(nbugs %in% c(N)) %>% filter(spacing %in% c(s))
  center_fits <-b %>% sigmoid(0.2,4) %>% centers() %>% mutate(nbugs=N, spacing =s)
  steep_fits <-b %>% sigmoid(0.2,4) %>% slopes() %>% mutate(nbugs=N, spacing =s)
  bugset_fits <- inner_join(center_fits,steep_fits,by=c("ks_pct","nbugs","spacing"))
  lfit <- lm(mu_50~ks_pct, data=bugset_fits)


  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_lineplot.png")
  png(here::here("output","diagnostic",fname), width=800,height=800,units="px")
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(x=bugset_fits$ks_pct,y=bugset_fits$mu_50,pch=16,col='blue',main=title)
  abline(lfit)
  dev.off()

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_diag.png")
  png(here::here("output","diagnostic",fname), width=800,height=800,units="px")
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

# TODO have to fit various sigmoids with diff init guesses for 2x2 2.5
#linear_fit_summaries<-rbind(linear_fit_summaries,record_lmfit(probs,4,2.5))
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
  b <- .data %>% filter(nbugs %in% c(N)) %>% filter(spacing %in% c(s))
  center_fits <-b %>% sigmoid(0.2,4) %>% centers() %>% mutate(nbugs=N, spacing =s)
  steep_fits <-b %>% sigmoid(0.2,4) %>% slopes() %>% mutate(nbugs=N, spacing =s)
  bugset_fits <- inner_join(center_fits,steep_fits,by=c("ks_pct","nbugs","spacing"))
  bugset_fits$k2 <- bugset_fits$ks_pct^2
  lfit <- lm(mu_50~ks_pct+k2, data=bugset_fits)

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_quadterm_lineplot.png")
  png(here::here("output","diagnostic",fname), width=800,height=800,units="px")
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(x=bugset_fits$ks_pct,y=bugset_fits$mu_50,pch=16,col='black',cex=1.5)
  pred <- predict(lfit)
  ix <- sort(bugset_fits$ks_pct, index.return=TRUE)$ix
  lines(bugset_fits$ks_pct[ix], pred[ix], col='red',lwd=2)
  dev.off()

  fname<-glue::glue("linfit_{sqrt(N)}x{sqrt(N)}_{s}_spacing_quadterm_diagnostic.png")
  png(here::here("output","diagnostic",fname), width=800,height=800,units="px")
  par(mfrow=c(2,2))
  title<-glue::glue("Diagnostic lm(mu_50~ks_pct+k2) {sqrt(N)}x{sqrt(N)} {s} spacing")
  plot(lfit,main=title)
  dev.off()
  return(lfit)
}

record_quadfit <-function(df,N,spacing){
  fit2frame(df %>%lm_quad_fits(N,spacing),type='lm_quad',N,spacing)
}

# TODO have to fit various sigmoids with diff init guesses for 2x2 2.5
#linear_fit_summaries<-rbind(linear_fit_summaries,record_quadfit(probs,4,2.5))
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
