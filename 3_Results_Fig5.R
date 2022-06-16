
############################################################################################
   ## TRANSSIZ - bacterial communities ##
############################################################################################

# set working directory; save and load
setwd("D:/AWI_MPI/collaborations/TRANSSIZ/Rstats")
load("TRANSSIZ.Rdata")


###################################################################################
   ###  MeSH BY LATITUDE -- Fig. 5 ###
###################################################################################

ggplot(MeSH_DMS) +
  geom_point(aes(
    x=lat, y=MeSH_frac),
    color="palevioletred4", size=0.5)+
  scale_x_continuous(
    breaks=c(60,65,70,75,80))+
  scale_y_continuous(
    breaks=c(10,20,30,40,50))+
  labs(
    x="Latitude (°N)",
    y="MeSH (% of DMS+MeSH)") +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank())

###########################
# remove temporary data
rm(r, p)

# save
save.image("TRANSSIZ.Rdata")

