# Plot the relative volumes of an example simulation over time, create an animated gif of plot
#
library(here) # manage paths, keep @JennyBryan from incinerating our machine
library(readr)  # handle csv file import
library(dplyr)  # work with tidy data
library(tidyr) # used for pivots
library(ggplot2) # figure generation
library(gganimate) # animation
library(transformr) # animation
library(gifski) # animation
library(tools) # md5sum

# load common code
source('common.R')

# Set up and start logging ---------
start_logging("8_show_relative_growth.log")

# log of relative volumes from a local run of 3x3_10 seed 1701
# accompanying visualization of cells growing was from VTK files generated
# and created in paraview (VTK was the reason for running locally)
# Associated still image in SI is the last frame.
# Gifs from paraview exports were assembled using ffmpeg
# ffmpeg was additionally used to combine movies side by side and extract the last frame
crv <- read_csv(here::here('data','cell_rel_volumes_seed_blank_axes.csv'))

# just want the heterotroph relative volumes for the relevant bugs
hrv <- crv %>% select(step,het1_rv,het2_rv,het3_rv,het4_rv,het5_rv,het6_rv,het7_rv,het8_rv,het9_rv) %>%
  pivot_longer(cols = starts_with("het"),
               names_to = "colony", names_pattern = "het(.*)_rv",
               values_to = "Relative Volume")

#add initial conditions
colonies = c("1","2","3","4","5","6","7","8","9")
ic<- data.frame(step=rep(0,9),`Relative Volume` = rep(1/9,9),colony = colonies,check.names = FALSE)
hrv <- bind_rows(ic,hrv)

# create the plot
theme_update(axis.title = element_text(size=20),
             axis.text  = element_text(size=18,color="black",angle=0),
             strip.text = element_text(size=13),
             title= element_text(size=16,hjust=0.5),
             plot.margin =  margin(0, 0, 0, 6, "pt"))
#8a4589
p <- ggplot(hrv %>% filter(step<2500), aes(x=step,y=`Relative Volume`,group=colony,color=colony)) +
  geom_line(size=1.2) +
  geom_hline(yintercept = 1/9,alpha=0.7,linetype="dashed") +
  scale_y_continuous(labels=scales::percent_format(accuracy = 2L)) +
  scale_color_manual(name="Colony",values = c("#2f6b93", "#3f8c38",  "#8a4589", "#FFA500", '#dad628','#8b471f','#d56d99','#7c7a74','#c91717'),
                    labels = c("1","2","3","4","5","6","7","8","9"))+
  xlab('Timestep') +
  ylab('Relative Biomass Volume') +
  #scale_color_discrete(name = "Colony ID") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=.5,colour="black"),
        legend.position="top",
        legend.text = element_text(size=18),
        legend.title= element_text(size=18))
#p

# animate the plot so that it is in sync with the video
panim <- p +   geom_point(size=5) + geom_text(aes(label=colony),color="black",check_overlap = TRUE)  + transition_reveal(step)
animate(panim, duration = 15, fps = 25, width = 800, height = 400, renderer = gifski_renderer(loop = F))
anim_save(here::here('output','cell_rel_volumes.gif'))

fname <- "cell_rel_volumes.gif"
floc <- here::here("output",fname)
log_info(paste('Wrote', file.path("output",fname), ' MD5Sum: ',
               md5sum(floc)))
