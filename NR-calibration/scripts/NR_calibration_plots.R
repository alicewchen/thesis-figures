
#===========================
# Import and set up dataset 
#===========================


setwd(getwd())

library(fcs2R)
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
require(dplyr)
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
library(readxl)
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
library(readxl)
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


#==========================================================================
#Plotting NR fluorescence against neutral lipid content as percent dry mass
#==========================================================================
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggpmisc)
library(ltxplot)
load_theme_ltx()
###### 
y_axes = c("mean_logPEA",
           "mean_logPEA_adjusted",
           "ratio_logPEA_to_blank")
y_labels = c("Mean log NR fluorescence",
             "Mean log NR fluorescence adjusted",
             "log NR fluorescence:Blank")
             
y_axes = c("mean_PE_A",
           "mean_PEA_adjusted",
           "ratio_PEA_to_blank")
y_labels = c("Mean NR fluorescence", 
             "Mean log NR fluorescence".
             "NR fluorescence:Blank")
x_axes = c("percent_neutral_lipids_of_dry_mass", 
           "ug_neutral_lipids_per_million_cells",
           "mass_per_million_cells_ug"  )
##### 
p1 <- df %>% 
  ggscatter(data = . , y = "mean_PE_A", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     breaks = c(20000,60000,100000,140000)) +
  xlab("Neutral lipid content\n(% dry mass)") +
  ylab("Mean NR fluorescence") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) #+

p2 <- df %>% 
  ggscatter(data = . , y = "mean_PE_A", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                #limits = c(10000, 10^6)) +
                limits = c(10000, 10^5.5)) +
  xlab("Neutral lipid content\n(% dry mass)") +
  ylab("Mean NR fluorescence") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) #+
  #rremove("xlab")

p3 <- df %>%
  ggscatter(data = . , y = "mean_PEA_adjusted", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001)) +
  xlab("Neutral lipid content\n(% dry mass)") +
  ylab("Mean adjusted\nNR fluorescence") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) #+
  #rremove("xlab")

p4 <- df %>%
  ggscatter(data = . , y = "ratio_PEA_to_blank", x = "percent_neutral_lipids_of_dry_mass", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x = 24) +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  xlab("Neutral lipid content\n(% dry mass)") +
  ylab("NR fluorescence:blank\n(fold difference)") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
  #rremove("xlab")

P4 <- ggarrange(p1, p2, p3,p4,labels = c("A","B","C","D"),
          ncol = 2, nrow = 2,common.legend = TRUE)
annotate_figure(P4, 
                bottom = text_grob("Neutral lipid content (% dry mass)",
                size = 18, family = 'lmroman'))

#=============================================================
#NR fluorescence against mass of neutral lipids as ng per cell
#=============================================================

p5 <- df %>% 
ggscatter(data = . , y = "mean_PE_A", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x.npc = "left") +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     breaks = c(20000,60000,100000,140000)) +
  xlab(expression( paste("Neutral lipids (", mu, "g ",10^-6," cells)"))) +
  ylab("Mean NR fluorescence") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())#+
  #rremove("xlab")

p6 <- df %>%
  ggscatter(data = . , y = "mean_PE_A", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x.npc = "left") +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size =2) +
  xlab(expression( paste("Neutral lipids (", mu, "g ",10^-6," cells)"))) +
  ylab("Mean NR fluorescence") +
  labs(color= "N Starvation") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                #limits = c(10000, 10^6)) +
                limits = c(10000, 10^5.5)) +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
  #rremove("xlab")

p7 <- df %>% 
  ggscatter(data = . , y = "mean_PEA_adjusted", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x.npc = "left") +
  geom_point(aes(shape = Strain, color = factor(N_starvation)), size = 2) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001)) +
  xlab(expression( paste("Neutral lipids (", mu, "g ",10^-6," cells)"))) +
  ylab("Mean adjusted\nNR fluorescence") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
  #rremove("xlab")

p8 <- df %>%
  ggscatter(data = . , y = "ratio_PEA_to_blank", x = "ug_neutral_lipids_per_million_cells", 
            shape = "Strain",
            size = 2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), family = 'lmroman'), 
           size = 5, label.x.npc = "left") +
  geom_point(aes(shape=Strain, color = factor(N_starvation)), size = 2) +
  xlab(expression( paste("Neutral lipids (", mu, "g ",10^-6," cells)"))) +
  ylab("NR fluorescence:blank\n(fold difference)") +
  labs(color = "N Starvation") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
  #rremove("xlab")

P4 <- ggarrange(p5, p6, p7, p8,labels = c("A","B","C","D"),
              ncol = 2, nrow = 2,common.legend = TRUE)
annotate_figure(P4, bottom = text_grob(expression( paste("Neutral lipids (", mu, "g ",10^-6," cells)")),size =18))

P8 <- ggarrange(p5, p6, p7, p8,p1, p2, p3,p4,
                labels = c("A","B","C","D","E","F","G","H"),
                ncol = 2, nrow = 4,common.legend = TRUE)

#pdf("../analysis/mean_logPEA_all.pdf", height = 11, width = 8.5)
pdf("../analysis/mean_logPEA_all.pdf", height = 11, width = 8.5)
P8
dev.off()

#=============================================
# NR vs FSC  colored by neutral lipid content 
#=============================================

#name = paste("Neutral lipids",paste(expression("(", mu, "g ",10^-6," cells)")))) +

p1 <- ggplot(data = df, aes(y = mean_FSC,
                                 x = mean_PE_A)) +
  geom_point(aes(colour = ug_neutral_lipids_per_million_cells, 
                           shape = Strain, 
                           size = mass_per_million_cells_ug)) +
  scale_color_gradient(high = "red", low = "blue", 
                       name = expression(atop("Neutral lipids", paste("(", mu, "g ",10^-6," cells)")))) +
  xlab("Mean NR fluorescence") +
  ylab("Cell size (FSC)") +
  labs(size = expression(atop("dry mass",paste("(",mu, "g ",10^-6," cells)")))) +
  scale_x_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     breaks = c(20000,60000,100000,140000)) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 10^-3)) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(colour = "black", size = 14)) +
  ggtitle(" \n ")

p2<-ggplot(data = df, aes(y = mean_FSC,
                        x = mean_PE_A))+
  geom_point(aes(colour = percent_neutral_lipids_of_dry_mass, 
                 shape = Strain,
                 size = mean_FSC)) +
  scale_color_gradient(high = "green", low = "blue", 
                       name =  "Neutral lipid content\n(% dry mass)") +
  xlab("Mean NR fluorescence") +
  ylab("Cell size (FSC)") +
  labs(size = "mean FSC") +
  scale_x_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     breaks = c(20000,60000,100000,140000)) +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 10^-3))+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black", size = 14))+
  ggtitle(" \n ")

P2<-ggarrange(p1, p2,labels = c("A","B"),
              ncol = 2, nrow = 1,common.legend = TRUE)
P2