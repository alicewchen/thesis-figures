library(fcs2R)
library(ltxplot)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(MuMIn)
library(car)
load_theme_ltx()
#####
# Extract fluorescence data from gated .fcs files
#####

#mean_PEA <- fcs_to_mean_PEA('../data/raw/fcs', '../data/clean/mean_PEA.txt')
mean_PEA <- fcs_to_mean_logPEA('../data/raw/fcs', '../data/clean/mean_logPEA.txt')

#####
# Merge mean_PEA with treatment information using old clean dataframe
#####

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
      )))
    #
  
#####
# Dredge
#####

options(contrasts = c("contr.sum", "contr.poly")) #Set this option to run type 3 ANOVA
lmdiff <- lm(mean_PE_A~ Latitude*Light*Timepoint*Temperature, data = df)
lmdiff <- lmer(mean_PE_A ~ Latitude*Temperature+ Latitude*Light+ Light*Temperature+ (1|Strain) + (1|Timepoint), data = clean)
summary(lmdiff)
lmdiff_dredge<- 0
options(na.action = "na.fail") #Set this option to run dredge
lmdiff_dredge <- dredge(lmdiff, evaluate = TRUE, rank = "AICc")
lmdiff_subset <- get.models(lmdiff_dredge, subset = delta < 2)
lmdiff_avg <- model.avg(lmdiff_subset)
summary(lmdiff_avg)
options(na.action = "na.omit") #Set this option to restore to default option of R

model <- lm(mean_PE_A ~ Latitude*Temperature + Latitude*Light+ Latitude*Timepoint + Light*Temperature + Light*Timepoint , data = df)
aov<- Anova(model, type = 3)
aov
plot(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))

#####
# ANOVA
#####

library(car)
model <- lmer(mean_PE_A ~ Latitude*Temperature+ Latitude*Light+ Light*Temperature+ (1|Strain) + (1|Timepoint) , data = clean)
aov<- Anova(model, type = 3)
aov
plot(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))

#remove outlier
clean <- df %>% filter(!(mean_PE_A >= 58625))

#check if random effects are needed
model <- lm(mean_PE_A ~ Latitude*Temperature*Light*Strain*Timepoint, data = clean)

#all interactions
model <- lm(mean_PE_A ~ Latitude*Light+ Latitude*Timepoint+ Latitude*Temperature + Temperature*Light+ Temperature*Timepoint, data = clean)
plot(model)
aov<- Anova(model, type = 3)
aov

###Use this model (with outlier)
model <- lmer(mean_PE_A ~ 
                Latitude*Light+ 
                Latitude*Timepoint+ 
                Latitude*Temperature + 
                Temperature*Timepoint +  
                Temperature*Light+ 
                Timepoint*Light+
                (1|Strain), data = df)
aov<- Anova(model, test="F", type="III")
sink(file = "../analysis/all_interactions_with_outlier.txt")
summary(model)
aov
sink(file = NULL)
pdf("../analysis/diagnostic_plots_with_outlier.pdf", height = 6, width = 6)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))
dev.off()

#Same model as above without outlier
model <- lmer(mean_PE_A ~ 
                Latitude*Light+ 
                Latitude*Timepoint+ 
                Latitude*Temperature + 
                Temperature*Timepoint +  
                Temperature*Light+ 
                Timepoint*Light+
                (1|Strain), data = clean)
aov<- anova(model)
sink(file = "../analysis/all_interactions_no_outlier.txt")
summary(model)
aov
sink(file = NULL)
pdf("../analysis/diagnostic_plots_no_outlier.pdf", height = 6, width = 6)
plot(model)
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model))
dev.off()

#Check equal variance among groups ####
plot( model, resid(., type = "pearson") ~ fitted(.) | Latitude*Temperature,
      id = 0.05, adj = -0.3 )
plot( model, resid(., type = "pearson") ~ fitted(.) | Latitude*Light,
      id = 0.05, adj = -0.3 )
plot( model, resid(., type = "pearson") ~ fitted(.) | Light*Temperature,
      id = 0.05, adj = -0.3 )
plot( model, resid(., type = "pearson") ~ fitted(.) | Latitude,
      id = 0.05, adj = -0.3 )
#####


######
# Plots
######
p[[1]]<-df %>% group_by(Latitude, Temperature) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.character(Latitude), y = mean, colour = Temperature)) +
  geom_point()+ 
  geom_line(aes(group = Temperature))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.1)+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  xlab(expression("Latitude ("*degree*"N)")) +
  ylab("Mean\nNR fluorescence") +
  labs(color = "Temperature") +
  theme_latex(base_size = 12) +
  #scale_color_manual(values=c("#D55E00", "#0072B2")) +
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
  

p[[2]]<-df %>% group_by(Latitude, Light) %>%
  summarize(mean = mean(mean_PE_A),
             se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.character(Latitude), y = mean, color = Light)) +
  geom_point(position = position_dodge(width = 0.3))+ 
  geom_line(aes(group = Light),position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.1, position = position_dodge(width = 0.3))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  xlab(expression("Latitude ("*degree*"N)")) +
  ylab("Mean\nNR fluorescence") +
  labs(color = "Light") +
  theme_latex(base_size = 12) +
  #scale_colour_manual(values=c("#D55E00", "#0072B2"))+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

p[[3]]<-df %>% group_by(Latitude, Timepoint) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.factor(Timepoint), y = mean, color = as.factor(Latitude))) +
  geom_point(position = position_dodge(width = 0.25))+ 
  geom_line(aes(group = Latitude),position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.1, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  xlab("Hours of treatment") +
  ylab("Mean\nNR fluorescence") +
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
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.1, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  xlab("Hours of treatment") +
  ylab("Mean\nNR fluorescence") +
  labs(color = "Light") +
  theme_latex(base_size = 12)+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

p[[5]]<-df %>% group_by(Light, Temperature) %>%
  summarize(mean = mean(mean_PE_A),
            se = sd(mean_PE_A)/sqrt(n())) %>%
  ggplot(., aes(x = as.factor(Temperature), y = mean, color = Light)) +
  geom_point(position = position_dodge(width = 0.25))+ 
  #geom_line(aes(group = Light),position = position_dodge(width = 0.25))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.1, position = position_dodge(width = 0.25))+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001),
                     limits = c(12000, 30000)) +
  xlab(expression("Temperature ("*degree*"C)")) +
  ylab("Mean\nNR fluorescence") +
  labs(color = "Light") +
  theme_latex(base_size = 12)+
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_text(family = ),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())


P <- ggarrange(p[[1]], p[[2]],p[[4]],p[[5]], labels = c("A","B","C","D"),
                ncol = 2, nrow = 2)

pdf("../analysis/interactions01.pdf", height = 8.5, width = 8.5)
P
dev.off()

pdf("../analysis/interactions02.pdf", height = 6, width = 6)
p[[3]]
dev.off()




