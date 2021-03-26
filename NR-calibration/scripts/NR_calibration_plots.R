library(readxl)
library(fcs2R)
require(dplyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(ggpmisc)
library(ltxplot)
load_theme_ltx()

#===========================
# Import and set up dataset 
#===========================
setwd(getwd())

#Calculate mean FCS data
cell_density<-fcs_to_cell_density("../data/raw/fcs/","../data/clean/cell_density.txt")
fsc_ssc<-fcs_to_mean_fsc_ssc("../data/raw/fcs/","../data/clean/fcs_ssc.txt")
mean_PEA<-fcs_to_mean_PEA("../data/raw/fcs/","../data/clean/mean_PEA.txt")
mean_logPEA<-fcs_to_mean_logPEA("../data/raw/fcs/","../data/clean/mean_logPEA.txt") 

#Import cell density, FSC, and PE-A data from new gating method data
cell_density <- read.table("../data/clean/cell_density.txt")
fsc_ssc <- read.table("../data/clean/fcs_ssc.txt")
mean_PEA <- read.table("../data/clean/mean_PEA.txt")
mean_logPEA <- read.table("../data/clean/mean_logPEA.txt")%>% 
  rename(mean_logPEA = mean_PE_A)

#merge all raw data into one df 
raw_data<-left_join(cell_density,fsc_ssc, by = c("well_id","plate")) %>%
  left_join(.,mean_PEA, by = c("well_id","plate")) %>%
  left_join(.,mean_logPEA, by = c("well_id","plate")) %>%
  mutate(mean_FSC = as.numeric(as.character(mean_FSC)),
         mean_SSC = as.numeric(as.character(mean_SSC)),
         mean_PE_A = as.numeric(as.character(mean_PE_A)),
         mean_logPEA = as.numeric(as.character(mean_logPEA)),
         cell_density = as.numeric(as.character(cell_density)))
write.table(raw_data, "../data/clean/fcs_data.txt", quote = FALSE, sep = "\t")

#create table of cell density (cells/mL) per falcon tube
X96_well_plate_layout <- read_excel("../data/raw/96_well_plate_layout.xlsx")

X96_well_plate_layout$plate <- as.numeric(X96_well_plate_layout$plate)
mean_tube_cell_density_per_mL <- raw_data %>% 
  full_join(.,X96_well_plate_layout, by = c("well_id","plate")) %>%
  na.omit(.) %>% 
  filter(treatment != "15") %>%
  group_by(tube_id) %>%
  summarize(mean_tube_cell_density_per_mL = mean(cell_density)*1000)
write.table(mean_tube_cell_density_per_mL, "../data/clean/mean_tube_cell_density.txt", sep = "\t")
#copy and paste output into glass vial mass and record

#Calibrate mean_PEA to blank (0_uL Nile Red Treatment)
raw_data <- full_join(raw_data,X96_well_plate_layout, by = c("well_id","plate"))
write.table(raw_data, "../data/clean/fcs_data.txt", quote = FALSE, sep = "\t")

blank_data <- raw_data %>% dplyr::filter(treatment == "0") %>%
  group_by(plate,Strain) %>%
  #mutate(Blank = 10^mean(mean_logPEA)) %>%
  dplyr::summarize(Blank = mean(mean_PE_A))
clean_data <- raw_data %>% dplyr::filter(treatment != "0") %>% 
  full_join(., blank_data, by = c('plate','Strain')) %>%
  #mutate(mean_PE_A = 10^mean_logPEA) %>%
  mutate(mean_PEA_adjusted = mean_PE_A - Blank,
         ratio_PEA_to_blank = mean_PE_A/Blank) %>%
  select(., -c(Blank))

#import neutral lipid content info
Glass_vial_mass_and_record <- read_excel("../data/raw/neutral_lipids/Glass vial mass and record.xlsx", 
                                         col_types = c("numeric", "text", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric"))

#summarize neutral lipid info in glass vial mass and record and merge to clean data

neutral_lipid_raw <- Glass_vial_mass_and_record %>%
  group_by(tube_id) %>%
  summarize(ug_neutral_lipids_per_million_cells = mean(mass_of_neutral_lipids_per_cell_ng*1000),
            percent_neutral_lipids_of_dry_mass=mean(percent_neutral_lipids_of_dry_mass),
            mass_per_million_cells_ug = mean(mass_of_dried_sample_in_falcon_tube/number_of_cells_in_falcon_tube*10^9))
            

df <- left_join(clean_data, neutral_lipid_raw, by = c("tube_id")) %>% filter(treatment == "3")%>%
  na.omit(.) %>% filter(treatment != "15")
write.table(df, "../data/clean/clean_data.txt")

df$N_starvation <- factor(df$N_starvation, levels = c("Before", "After"))


#=============================================================
#NR fluorescence against mass of neutral lipids as ng per cell
#=============================================================
require(dplyr)
require(ggplot2)
require(scales)
require(ggpubr)
require(ggpmisc)
require(ltxplot)
load_theme_ltx()

p1 <- df %>% 
  ggscatter(data = . , y = "mean_PE_A", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x.npc = "left") +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     breaks = c(20000,60000,100000,140000)) +
  labs(x=expression( paste("Neutral lipid (", mu, "g ",10^-6," cells)")),
       y="NR fluorescence",
       tag="A",
       color = "N-starvation") +
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))


p2 <- df %>%
  ggscatter(data = . , y = "mean_PE_A", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x.npc = "left") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                #limits = c(10000, 10^6)) +
                limits = c(10000, 10^5.5))+
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size =2) +
  labs(x=expression( paste("Neutral lipid (", mu, "g ",10^-6," cells)")),
       y="NR fluorescence",
       tag = "B",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))
  #rremove("xlab")

p3 <- df %>% 
  ggscatter(data = . , y = "mean_PEA_adjusted", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x.npc = "left") +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001)) +
labs(x=expression( paste("Neutral lipid (", mu, "g ",10^-6," cells)")),
       y="Adjusted\nNR fluorescence",
       tag = "C",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

p4 <- df %>%
  ggscatter(data = . , y = "ratio_PEA_to_blank", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x.npc = "left") +
  geom_point(aes(shape=Strain, color = factor(N_starvation)), size = 2) +
  labs(x=expression( paste("Neutral lipid (", mu, "g ",10^-6," cells)")),
       y="NR fluorescence:blank\n(fold difference)",
       tag = "D",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

#==========================================================================
#Plotting NR fluorescence against neutral lipid content as percent dry mass
#==========================================================================
require(dplyr)
require(ggplot2)
require(scales)
require(ggpubr)
require(ggpmisc)
require(ltxplot)
load_theme_ltx()

p5 <- df %>% 
  ggscatter(data = . , y = "mean_PE_A", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     breaks = c(20000,60000,100000,140000)) +
  labs(x="Neutral lipid (% dry mass)",
       y="NR fluorescence",
       tag = "E",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

p6 <- df %>% 
  ggscatter(data = . , y = "mean_PE_A", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                #limits = c(10000, 10^6)) +
                limits = c(10000, 10^5.5)) +
  labs(x="Neutral lipid (% dry mass)",
       y="NR fluorescence",
       tag = "F",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

p7 <- df %>%
  ggscatter(data = . , y = "mean_PEA_adjusted", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x = 24)+
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001)) +
  labs(x="Neutral lipid (% dry mass)",
       y="Adjusted\nNR fluorescence",
       tag = "F",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))


p8 <- df %>%
  ggscatter(data = . , y = "ratio_PEA_to_blank", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 4, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  labs(x="Neutral lipid (% dry mass)",
       y="NR fluorescence:blank\n(fold difference)",
       tag = "F",
       color = "N-starvation")+
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))



P8 <- ggarrange(p1, p2, p3,p4,p5, p6, p7, p8,
                ncol = 2, nrow = 4,common.legend = TRUE)


#Figure 3
pdf("../analysis/mean_logPEA_all.pdf", height = 9, width = 6.5)
P8
dev.off()

#==============================
# Appendix
#==============================
#Figure S6: NR vs FSC  colored by neutral lipid content, size by dry mass ####
p1 <- ggplot(data = df, aes(y = mean_PE_A,
                            x = mean_FSC)) +
  geom_point(aes(colour = ug_neutral_lipids_per_million_cells, 
                 shape = Strain, 
                 size = mass_per_million_cells_ug)) +
  scale_color_gradient(high = "#E69F00", low = "#56B4E9", 
                       name = expression(atop("Neutral lipids", paste("(", mu, "g ",10^-6," cells)")))) +
  labs(x = "Cell size (FSC-A)",
       y = "NR fluorescence",
       size = expression(atop("dry mass",paste("(",mu, "g ",10^-6," cells)"))),
       tag = "A") +
  scale_x_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(0,500000)) +#,
                     #breaks = c(20000,60000,100000,140000)) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 10^-3),
                     limits = c(0,200000)) +
  guides(shape=FALSE, size = FALSE) + 
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = "top",
        legend.text=element_text(size=10),
        legend.margin=margin(t = 0,l= 0, unit='cm'),
        plot.margin = margin(r = 10,b=10, unit = 'pt'))

p2<-ggplot(data = df, aes(y = mean_PE_A,
                          x = mean_FSC))+
  geom_point(aes(colour = percent_neutral_lipids_of_dry_mass, 
                 shape = Strain,
                 size = mass_per_million_cells_ug)) +
  scale_color_gradient(high = "#D55E00", low = "#0072B2", 
                       name =  "Neutral lipids\n(% dry mass)") +
  labs(x = "Cell size (FSC-A)",
       y = "NR fluorescence",
       size = expression(atop("dry mass",paste("(",mu, "g ",10^-6," cells)"))),
       tag = "B") +
  scale_x_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(0,500000)) + #,
                     #breaks = c(20000,60000,100000,140000)) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 10^-3),
                     limits = c(0,200000))+
  guides(shape=FALSE, size = FALSE) + 
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = "top",
        legend.text=element_text(size=10),
        legend.margin=margin(t = 0,l= 0, unit ='pt'),
        plot.margin = margin(r = 10,b=10, unit = 'pt'))
  #ggtitle(" \n ")
p_empty<- ggplot(data = df, aes(y = mean_PE_A,
                                x = mean_FSC))+
  geom_point(aes(colour = percent_neutral_lipids_of_dry_mass, 
                 shape = Strain,
                 size = mass_per_million_cells_ug)) +
  scale_color_gradient(#high = "green", low = "blue", 
                       name =  "Neutral lipids\n(% dry mass)") +
  labs(x = "Cell size (FSC-A)",
       y = "NR fluorescence",
       size = expression(atop("dry mass",paste("(",mu, "g ",10^-6," cells)")))) +
  scale_x_continuous(labels = unit_format(unit = "K", scale = 0.001)) +#,
                     #breaks = c(20000,60000,100000,140000)) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 10^-3))+
  guides(colour = FALSE) + 
  theme_latex(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = "top",
        legend.text=element_text(size=10),
        legend.margin=margin(t = 0, l= 0, unit='cm'))

pdf("../analysis/Appendix FigS6.pdf", height = 6, width = 6.5)
lay <- rbind(c(1,2),
             c(3,3))
grid.arrange(p1, p2, p_empty, layout_matrix = lay)
#remove p_empty when taking screenshot
dev.off()

#Figure S7: PCA ####
PCA_data <- df %>%
  select(percent_neutral_lipids_of_dry_mass, ug_neutral_lipids_per_million_cells,mass_per_million_cells_ug, mean_FSC)

pca <- prcomp(PCA_data, scale = TRUE)
pca
summary(pca)
plot(pca)

PC1 <- pca[["x"]][,1]
PC2 <- pca[["x"]][,2]

df_temp <- df %>% cbind(., PC1,PC2)

pcplot1 <- ggscatter(data = df_temp , y = "PC2", x = "PC1", 
                     shape = "Strain",
                     conf.int = FALSE) +
  geom_point(aes(shape = Strain, size = mean_FSC, color = percent_neutral_lipids_of_dry_mass)) +
  scale_color_gradient(high = "green", low = "blue", 
                       name =  "Neutral lipid content\n(% dry mass)") +
  xlab("PC1") +
  ylab("PC2") +
  xlim(-2.5,4) +
  ylim(-2.5,2.5) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black")) +
  ggtitle(" ")

pcplot2 <- ggscatter(data = df_temp , y = "PC2", x = "PC1", 
                     shape = "Strain",
                     conf.int = FALSE) +
  geom_point(aes(shape = Strain, size = mass_per_million_cells_ug , color = ug_neutral_lipids_per_million_cells )) +
  scale_color_gradient(high = "red", low = "blue", 
                       name =  "Neutral lipid content\n(ng/cell)") +
  xlab("PC1") +
  ylab("PC2") +
  xlim(-2.5,4) +
  ylim(-2.5,2.5) +
  theme_classic(base_size = 10) +
  theme( axis.text = element_text( colour = "black")) +
  ggtitle(" ")

P2 <- ggarrange(pcplot1,pcplot2, labels = c("A","B"), ncol = 2, nrow = 1)
P2

pdf("../analysis/PCA_coloured.pdf", height = 4, width = 9.5)
ggarrange(pcplot1,pcplot2, labels = c("A","B"), ncol = 2, nrow = 1, legend = "right" )
dev.off()

pcplot3 <- ggscatter(data = df_temp , y = "PC2", x = "PC1",
                     shape = "Strain",
                     color= "N_starvation",
                     size = 4,
                     conf.int = FALSE) +
  geom_point(aes(shape = Strain, color = N_starvation)) +
  #scale_color_manual(labels=c("Before","After")) +
  labs(x = "PC1 (83.11%)",
       y = "PC2 (16.16%)",
       color = "N-starvation") +
  xlim(-2.5,4) +
  ylim(-2.5,2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")   +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = "right")
pdf("../analysis/Appendix FigS7.pdf", height = 5, width = 6.5)
pcplot3
dev.off()

