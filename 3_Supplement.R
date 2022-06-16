
############################################################################################
### TRANSSIZ - bacterial communities ###

## this script: Supplementary Figures
############################################################################################

## Chlorophyll

ggscatter(
  Chl, x="FerryBox_chla (ug/L)", y="HPLC_chla (ug/L)", 
  add = "reg.line", color="mediumseagreen",
  conf.int = T, size=3) +
  labs(
    x="Chla - Ferrybox [µg/L]",
    y="Chla - HPLC [µg/L]") +
  stat_cor(label.y=4.4) +
  #stat_regline_equation(label.y=4.2) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        panel.grid.major = element_line(size=1))

##########################

## Henry's Law

ggscatter(
  Henry, x="Henry's Law Cst [M/Atm]", 
  y="Calibration factor [nM/ppbv]", 
  add = "reg.line", color="gray16",
  conf.int = F) +
  labs(
    x="Henry's Law Constant [M/Atm]",
    y="Calibration factor [nM/ppbv]") +
  geom_point(aes(color=VOC), size=4.4) +
  scale_color_manual(values=c(
    "Acetonitrile"="darkorchid3",
    "Acetone"="gray22",
    "Methanol"="steelblue2",
    "Acetaldehyde"="gold3",
    "Isoprene"="lightseagreen")) +
  scale_x_reverse(labels = scientific) +
  scale_y_reverse(labels = scientific) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        panel.grid.major = element_line(size=1)) 

##########################

## Acetone

ggscatter(
  Acetone, x="mol / L water", 
  y="Headspace [ppb]", 
  add = "reg.line", color="gray22",
  conf.int = T) +
  labs(
    x="Acetone [mol/L in water]",
    y="Headspace [ppb]") +
  # stat_regline_equation(label.y = 4.2) +
  geom_point(color="gray22", size=4.4) +
  scale_x_continuous(labels = scientific) +
   theme_classic() +
  stat_cor(label.y=6) +
  theme(axis.ticks = element_blank(),
        panel.grid.major = element_line(size=1)) 

##############################

## Chl vs DMS >80°N

ggscatter(
  Chl_80N, x="Chla [ug/L]", y="DMS [nM]", 
  add = "reg.line", color="cadetblue4",
  conf.int = T, size=3) +
  labs(
    x="Chla [ug/L]",
    y="DMS [nM]") +
  stat_cor(label.y=34) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        panel.grid.major = element_line(size=1))

##############################

## Chl vs Isoprene >80°N

# without outlier station
ggscatter(data=subset(
  Chl_80N, Station!="0019-5"),
  x="Chla [ug/L]", y="Isoprene [pM]", 
  add = "reg.line", color="mediumturquoise",
  conf.int = T, size=3) +
  labs(
    x="Chla [ug/L]",
    y="Isoprene [pM]") +
  stat_cor(label.y=8) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        panel.grid.major = element_line(size=1))

# with outlier station
ggscatter(
  Chl_80N, 
  x="Chla [ug/L]", y="Isoprene [pM]", 
  add = "reg.line", color="mediumturquoise",
  conf.int = T, size=3) +
  labs(
    x="Chla [ug/L]",
    y="Isoprene [pM]") +
  stat_cor(label.y=8) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        panel.grid.major = element_line(size=1))

