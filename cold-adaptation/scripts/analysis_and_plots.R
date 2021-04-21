library(fcs2R)
library(ltxplot)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(MuMIn)
library(car)
library(lme4)
load_theme_ltx()

#=======================
# Data manipulation
#=======================
# Extract fluorescence data from gated .fcs files ####
mean_PEA <- fcs_to_mean_logPEA('../data/raw/fcs', '../data/clean/mean_logPEA.txt')

# Merge mean_PEA with treatment information using old clean dataframe ####

mean_PEA <- read.table("../data/clean/mean_logPEA.txt")
ids <- read.table("../data/raw/Data.csv", sep = ",", header = TRUE) %>% 
  filter(Stained.Unstained == 1) %>%
  select(Strain, Temperature, Light) %>%
  #rename(plate = plate_number) %>%
  distinct() %>%
  mutate(id = as.factor(row_number()))
df <- read.table("../data/raw/Data.csv", sep = ",", header = TRUE) %>% 
  select(well_id, plate_number, Strain, Stained.Unstained, Temperature, Light, Timepoint) %>%
  distinct() %>%
  rename(plate = plate_number) %>%
  left_join(., mean_PEA, by = c('well_id', 'plate')) %>% 
  left_join(., ids, by = c('Strain', 'Temperature', 'Light')) %>%
  drop_na() %>%
  mutate(
    mean_PE_A = 10^mean_PE_A,
    Strain = as.factor(Strain),
    Stained.Unstained = as.factor(Stained.Unstained),
    Light = as.factor(Light), 
    Temperature = as.factor(ordered(Temperature, c("4","25")))) %>%
  filter(Stained.Unstained == "1") %>%
  #Add latitude
  mutate(                            
    Latitude = case_when(
      Strain == "1690" ~ 42.4,
      Strain == "1691" ~ 42.4,
      Strain == "2343" ~ 28.1,
      Strain == "3064" ~ 45.3,
      Strain == "3082" ~ 45.4
    )) %>% 
  #Add geographical region
  mutate(
    N_S = as.factor(case_when(
      Latitude == "42.4" ~ "N",
      Latitude == "45.3" ~ "N",
      Latitude == "45.4" ~ "N",
      Latitude == "28.1" ~ "S"
      ))) %>%
  mutate(N_S = as.factor(ordered(N_S, c("S","N"))))

# Make a copy of dataframe without the outlier ####
clean <- df %>% filter(!(mean_PE_A >= 58625))

#==========================
# Linear model and ANOVA
#==========================
# check if random effects are needed ####
model <- lm(mean_PE_A ~ Latitude*Temperature*Light*Strain*Timepoint, data = clean)
aov<- Anova(model)
aov 
#Strain is a signficant term -> might need to set Strain as random effect


# Linear model including all two-way interactions (with outlier)####
sink(file = "../analysis/all_interactions_with_outlier.txt")
model <- lm(mean_PE_A ~ 
                Latitude*Light+ 
                Latitude*Timepoint+ 
                Latitude*Temperature + 
                Temperature*Timepoint +  
                Temperature*Light+ 
                Timepoint*Light, data = df)
aov<- Anova(model, test="F", type="III")
summary(model)
aov
sink(file = NULL)
pdf("../analysis/diagnostic_plots_with_outlier.pdf", height = 6, width = 6)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))
dev.off()

# Linear model including all two-way interactions (without outlier) ####
sink(file = "../analysis/all_interactions_no_outlier.txt")
model <- lm(mean_PE_A ~ 
              Latitude*Light+ 
              Latitude*Timepoint+ 
              Latitude*Temperature + 
              Temperature*Timepoint +  
              Temperature*Light+ 
              Timepoint*Light, data = clean)
aov<- Anova(model, test="F", type="III")
summary(model)
aov
sink(file = NULL)
pdf("../analysis/diagnostic_plots_no_outlier.pdf", height = 6, width = 6)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))
dev.off()

# Linear model including all two-way interactions (with outlier); exclude one strain ####
sink(file = "../analysis/no_3064.txt")
model <- df %>% 
  filter(Strain != "3064") %>%
  lm(mean_PE_A ~ Latitude + Temperature+ Light + 
       Latitude*Light+ 
       Latitude*Timepoint+ 
       Latitude*Temperature + 
       Temperature*Timepoint +  
       Temperature*Light+ 
       Timepoint*Light, data = .)
aov<- aov(model)
summary(model)
summary(aov)
sink(file = NULL)

sink(file = "../analysis/no_3082.txt")
model <- df %>% 
  filter(Strain != "3082") %>%
  lm(mean_PE_A ~ Latitude + Temperature+ Light + 
       Latitude*Light+ 
       Latitude*Timepoint+ 
       Latitude*Temperature + 
       Temperature*Timepoint +  
       Temperature*Light+ 
       Timepoint*Light, data = .)
aov<- aov(model)
summary(model)
summary(aov)
sink(file = NULL)

######
# Plots
######
p[[1]]<-df %>% group_by(Latitude, Temperature) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.character(Latitude), y = mean, colour = Temperature)) +
  geom_point()+ 
  geom_line(aes(group = Temperature))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.05)+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  scale_color_manual(labels = c(expression("4"*degree*"C"), expression("25"*degree*"C")), 
                     values = c("#56B4E9","#E69F00"))+
  labs(x=expression("Latitude ("*degree*"N)"),
       y="NR fluorescence",
       color = "Temperature",
       tag="A")+
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.title = element_blank())
  

p[[2]]<-df %>% group_by(Latitude, Light) %>%
  summarize(mean = mean(mean_PE_A),
             se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.character(Latitude), y = mean, color = Light)) +
  geom_point(position = position_dodge(width = 0.3))+ 
  geom_line(aes(group = Light),position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.05, position = position_dodge(width = 0.3))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  scale_color_discrete(labels = c("No light", "Light"))+
  labs(x=expression("Latitude ("*degree*"N)"),
       y="NR fluorescence",
       tag ="B",
       colour = "Light")+
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.title = element_blank())

p[[3]]<-df %>% group_by(Latitude, Timepoint) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.factor(Timepoint), y = mean, color = as.factor(Latitude))) +
  geom_point(position = position_dodge(width = 0.25))+ 
  geom_line(aes(group = Latitude),position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.05, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  xlab("Hours of treatment") +
  ylab("NR fluorescence") +
  labs(color = "Latitude") +
  theme_latex(base_size = 12)+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

p[[4]]<-df %>% group_by(Light, Timepoint) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.factor(Timepoint), y = mean, color = Light)) +
  geom_point(position = position_dodge(width = 0.25))+ 
  geom_line(aes(group = Light),position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.05, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  #scale_color_discrete(labels = c("No light", "Light"))+
  labs(x="Hours of treatment",
       y="NR fluorescence",
       #color ="Light",
       tag = "C")+
  theme_latex(base_size = 12)+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = 'none')

p[[5]]<-df %>% group_by(Light, Temperature) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.factor(Temperature), y = mean, color = Light)) +
  geom_point(position = position_dodge(width = 0.25))+ 
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.05, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  scale_color_discrete(labels = c("No light", "Light"))+
  labs(x=expression("Temperature ("*degree*"C)"),
       y="NR fluorescence",
       colour="Light",
       tag="D")+
  theme_latex(base_size = 12)+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = 'none')

p[[6]]<-df %>% group_by(Temperature, Timepoint) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.factor(Timepoint), y = mean, color = Temperature)) +
  geom_point(position = position_dodge(width = 0.25))+ 
  geom_line(aes(group = Temperature),position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.05, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  #scale_color_discrete(labels = c("No light", "Light"))+
  labs(x="Hours of treatment",
       y="NR fluorescence",
       color ="Temperature",
       tag = "E")+
  theme_latex(base_size = 12)+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

pdf("../analysis/interactions01.pdf", height = 6.5, width = 6.5)

P <- grid.arrange(p[[1]], p[[2]],p[[4]],p[[5]],
                  ncol = 2, nrow = 2)
dev.off()

pdf("../analysis/interactions02.pdf", height = 6, width = 6)
p[[3]]
dev.off()




