
###################################################################################
   ###  LATITUDINAL COMMUNITY PATTERNS (Fig. 3a)
###################################################################################

# Hellinger transformation
hel = pseq.rel %>% transform_sample_counts(
  function(x) sqrt(x / sum(x)))

# Calculate NMDS
nmds <- ordinate(
  hel,"NMDS","bray") 

# export size 3x5
plot_ordination(
  hel, nmds) +
  geom_point(aes(
    color=lat, size=temperature)) + 
  scale_color_scico(
    palette="hawaii",
    begin=0.1, end=0.98,
    limits=c(58,83),
    breaks=c(60,70,80)) +
  scale_size(range = c(2,4)) +
  theme_bw() +
  theme(axis.ticks = element_blank())


############################################################################################
   ###  BACTERIA-VOC correlations (Fig. 3b)
############################################################################################

# Format ASVs; subset to most abundant
bac <- pseq.rel %>%
  filter_taxa(
    function(x) sum(x >= 0.1) >=5, T) %>% 
  transform_sample_counts(
    function(x) sqrt(x / sum(x))) 

# Extract taxdata
tax = as(tax_table(bac), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column("asv")

# reformat ASVs
asv = as(otu_table(bac), "matrix") %>% 
  as.data.frame() %>%
  rownames_to_column("asv") %>%
  reshape2::melt() %>%
  left_join(ENV, by=c("variable"="sample_title")) %>%
  group_by(asv, lat_short) %>%
  summarise(Abundance=mean(value)) %>%
  ungroup %>%
  spread(asv, Abundance) %>%
  drop_na()

# reformat VOCs
avg.voc <- VOCs %>% 
  mutate(lat_short=lat) %>%
  mutate(across(lat_short, round, 1)) %>%
  group_by(lat_short) %>%
  summarise(across(everything(), mean)) %>%
  ungroup  %>%
  filter(lat_short %in% asv$lat_short) %>%
  dplyr::select(c(
    "DMS","MeSH","MeSH_DMS","lat_short",
    "Acetonitrile","Acetaldehyde","CO",
    "Acetone","Isoprene","Chlorophyll"))

#######################

# Calculate pairwise correlations
asv %>% filter(
  lat_short %in% avg.voc$lat_short) %>%
corr.test(
  ., avg.voc, 
  use = "pairwise",
  method="spearman",
  adjust = "holm",
  alpha=.05, ci=T,
  minlength=5, normal=T) -> cor

# extract, subset most significant 
r <- cor$r 
r[abs(r) < -0.4 | abs(r) < 0.4 ] = NA 

# reformat; remove all-NA taxa
# remove taxa w/ latitude correlations
r <- as.data.frame(r) %>%
  rownames_to_column("asv") %>%
  filter(asv!="lat_short" & is.na(lat_short)) %>%
  filter_at(vars(2:11), any_vars(!is.na(.))) %>%
  left_join(tax) %>% 
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, round, 2) 

# reformat p-values
# remove taxa w/ latitude correlations
p <- as.data.frame(cor$p) %>%
  mutate(across(
    everything(), ~replace(., .>0.05, NA))) %>%
  rownames_to_column("asv") %>%
  filter(asv!="lat_short" & is.na(lat_short)) %>%
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, ~case_when(
    . < 0.05 & . > 0.01 ~ "*",
    . < 0.01 & . > 0.001 ~ "**",
    . < 0.001 ~ "***")) %>%
  dplyr::rename(pvalue=value)

# merge, reformat
# omit most unclassified taxa
# omit "double ASVs" per taxon (w/ same correlations)
cor <- right_join(r, p) %>%
  unite(sign, value, pvalue, sep="", remove=F) %>%
  mutate_at(c("sign"), ~gsub("NA", NA, .)) %>%
  mutate(across(value, as.numeric)) %>%
  filter(Genus !="Proteobacteria_uc" & !asv %in% c(
    "asv21","asv112","asv119","asv201","asv62",
    "asv50","asv147","asv88","asv216","asv117")) %>%
  mutate(asv = factor(asv, levels=c(
    "asv168","asv131","asv114","asv142","asv95",
    "asv72","asv118","asv151","asv81","asv125",
    "asv55","asv32","asv127","asv187","asv71",
    "asv178","asv124","asv189","asv45","asv80",
    "asv20","asv31","asv64","asv105","asv248"))) %>%
  mutate(id = paste(Genus,asv)) %>%
  mutate(id = str_replace(
    id, "_uc", "")) %>%
  mutate(id = str_replace(
    id, "_clade", "")) %>%
  mutate(id = str_replace(
    id, "_Clade_", " ")) %>%
  mutate(id = str_replace(
    id, "_cluster", "")) %>%
  mutate_at(c("Class"), ~gsub(
    "proteobacteria|microbiae|eroidia", "", .)) %>%
  drop_na(asv, pvalue) 

# order by Class
cor$id = factor(cor$id, levels=unique(
  cor$id[ order(cor$asv)]))

# plot; omit variables w-out significant correlations
# export size 3x4
ggplot(data=subset(cor, !variable %in% c(
  "lat_short","Acetaldehyde","CO")))  + 
  aes(x=id, y=variable, 
      fill=value, label=pvalue) +
  geom_tile(color="gray84") +
  geom_text(color="white", size=3.8) +
  scale_fill_scico(
    alpha = NULL,
    begin = 0.1,
    end = 0.9,
    direction = -1,
    na.value = "gray98",
    palette = "bamO",
    limits = c(-0.8, 0.8),
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    midpoint = 0) +
  facet_grid(
    Class ~ .,  
    scales="free", 
    space="free") +
  coord_flip() +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(
          angle=45, hjust=1, vjust=1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        #legend.position = "none",
        axis.title = element_blank())

# Figures A+B were exported as pdf
# pdfs imported to Inkscape, combined, and aesthetically modified


##################################################################################
   ###  MeSH + BACTERIA BY LATITUDE (Fig. 5)
###################################################################################

# Panel A -- MeSH cconcentrations
# export size 3x4
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
  theme_classic() + 
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank())

#######################

# Panel B -- latitudinal blocks

# DMS-MeSH fraction by block
avg <- MeSH_DMS %>%
  mutate(region=case_when(
    lat>60 & lat<70~"60-70",
    lat>70 & lat<80~"70-80",
    lat>80~"80-81",
    TRUE~"50-60")) %>%
  group_by(region) %>%
  drop_na(MeSH_frac) %>%
  summarize(MeSH_frac=mean(MeSH_frac)) %>%
  ungroup() 

# Major bacterial genera by block
psmelt(pseq.rel) %>% 
  filter(Genus %in% c(
    "Nitrincolaceae_uc",
    #"SAR11_Clade_Ia","Formosa",
    "OCS116_clade_uc",
    "Amylibacter","NS4",
    "Paraglaciecola")) %>%
  group_by(Genus, Sample, lat) %>%
  summarize(Abundance=sum(Abundance)) %>%
  mutate(region=case_when(
    lat>60 & lat<70~"60-70",
    lat>70 & lat<80~"70-80",
    lat>80~"80-81",
    TRUE~"50-60")) %>%
  group_by(Genus, region) %>%
  summarize(Abundance=mean(Abundance)) %>%
  ungroup() %>% 
  ggplot(aes(
    x=region, y=Abundance, 
    fill=Genus)) +
  geom_bar(width=0.8, stat="identity") +
  scale_fill_manual(
    values = scico(5, palette='vanimo')) +
  scale_x_discrete(expand = c(0.15,0.15)) +
  scale_y_continuous(expand = c(0.02,0.02)) +
  labs(y="Relative abundance",
       x="Latitudinal range") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        #axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust=0.5, vjust=0.5),
        axis.text = element_text(size=10),
        legend.position = "right")

# Figures were exported as pdf 
# pdfs imported to Inkscape, combined and aesthetically modified
