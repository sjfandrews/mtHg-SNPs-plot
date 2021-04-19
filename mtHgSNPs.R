# Script for plotting haplogroup defining mtSNPs on the mitochondrial genome 

library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(pals)
# library(ggrepel)

source("scripts/FamilyTreeDNA_mtSNPs.R")
source("scripts/mtDNAbdries.R")

# params 
y_heavy = 2   # y position for heavy strand
y_light = 1.8 # y position for light strand

# Import Haplogroup defining SNPs from HiMC and FamilyTreeDNA
snps <- read_csv('data/HiMC_mtSNPs.csv', skip = 1) %>% 
  mutate(pos = str_replace(mito_snp_id, "MT", "") %>% as.numeric()) %>% 
  rename(hg = associated_haplogroup) %>%
  select(-references) %>% 
  arrange(pos) %>% 
  left_join(mito_snps, by = c('pos', 'hg'))

# Data frame 
mtpos <- tibble(pos = seq(0,16567, by=1)) %>% 
  mutate(gene = addgenelabel(pos)) %>% 
  left_join(snps, by = "pos") 

heavy <- filter(bdries, strand == "+") %>% mutate(y = y_heavy)
light <- filter(bdries, strand == "-") %>% mutate(y = y_light)

# Format data for plotting (EUR only)
# code for angle from https://www.r-graph-gallery.com/296-add-labels-to-circular-barplot.html
dat.p = mtpos %>% 
  filter(hg %in% c("H", "HV", "I", "J", "JT", "K", "T", "U", "V", "W", "X", "R")) %>% 
  distinct(pos, .keep_all = T) %>% 
  mutate(val = y_heavy+0.2,
         angle = 90 - 360 * (pos-0.5) / 16567, 
         angle = ifelse(angle < -90, angle+180, angle), 
         angle = angle * -1, 
         hjust = ifelse( angle < -90, 1, 0)
  )

# Color pallets 
# n_features = distinct(bdries, feature) %>% nrow()
# gene_pals = pals::cols25(n = n_features)
gene_pals = c("#ff7f00", "#ff0000", "#33a02c", "#a6cee3", "#1F78C8", "#6A33C2", 'black')

n_hgs <- distinct(dat.p, hg) %>% nrow()
hg_pals = pals::cols25(n = n_hgs)

p1 <- ggplot() +
  coord_polar(direction = -1) +
  geom_text(aes(x = 0, y = 0, label = "Haplogroup \ndefining mtSNPs")) + 
  # Haplogroup mtSNPs
  geom_segment(data = dat.p, aes(x = pos, y = y_heavy, xend = pos, yend = val), size = 0.2, color = "black", alpha = 1) + 
  geom_point(data = dat.p, aes(x = pos, y = val, color = hg), size = 3) + 
  geom_text(data = dat.p, aes(x = pos, y = val+0.5, label = variant), angle = dat.p$angle, hjust = 0.5, size = 3) + 
  # Repel labels if I can get it working...
  # geom_text_repel(data = dat.p, aes(x = pos, y = val, label = variant), 
  #                 angle = dat.p$angle, 
  #                 direction = "x", 
  #                 hjust = dat.p$hjust,
  #                 ylim = c(2, NA)) + 
  scale_color_manual(values = hg_pals, "Haplogroup")  + 
  # Features - Heavy Strand
  new_scale_color() +
  geom_line(data = mtpos, aes(x = pos, y_heavy), color = 'black') + 
  geom_segment(data = heavy, aes(x = start, y = y, xend = end, yend = y, color = feature), size = 3) + 
  # Features - Light Strand
  geom_line(data = mtpos, aes(x = pos, y_light), color = 'black') + 
  geom_segment(data = light, aes(x = start, y = y, xend = end, yend = y, color = feature), size = 3) + 
  scale_color_manual(values = gene_pals, guide = 'none') +
  # Custom theme, remove everything
  theme_minimal() + 
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right", 
    legend.text = element_text(size=8), 
    text = element_text(size=8)
  ) 

ggsave("plots/mtDNA.png", plot = p1, w=6, h=6, dpi=300, units = "in")




































