library(readxl)
library(reshape)
library(dplyr)
library(stringi)
library(fcs2R)
library(ggplot2)
library(scales)
library(ggtext)
library(ltxplot)
load_theme_ltx()
############################
# Create well to Strain key#
############################

#Import files
w96_layout <- read_excel("../data/raw/w96_layout.xlsx")
w96_layout_key <- read_excel("../data/raw/w96_layout_key.xlsx")

#Rearrange data in plate format to dataframe
w96_layout_temp <- melt(as.data.frame(w96_layout), id=c("Plate","Row","N_Starvation","Stained"), na.rm = TRUE) 
w96_layout_temp$Key <- w96_layout_temp$value

#Assign Strains to each well
w96_key<-right_join(w96_layout_key,w96_layout_temp, by=c("Key"))

#Generate Well id and finalize Strain key
w96_key<- w96_key %>% 
  mutate(well_id = stri_join(Row,variable,sep="")) %>% 
  select(Strain, well_id, N_Starvation, Stained, Plate, Key)

#####################################################
# Extract FCS data: mean FL2 and cell concentration #
#####################################################

#Extract FCS data
fcs_to_cell_density("../data/raw/fcs/","../data/clean/cell_density.txt")
fcs_to_mean_logPEA("../data/raw/fcs/","../data/clean/mean_logPEA.txt")
fcs_to_mean_PEA("../data/raw/fcs/","../data/clean/mean_PEA.txt")
#fcs_to_median_PEA("../data/raw/fcs/","../data/clean/median_PEA.txt")
fcs_to_mean_fsc_ssc("../data/raw/fcs/","../data/clean/mean_FSC.txt")

#Import data from saved files
mean_PEA <- read.table("../data/clean/mean_logPEA.txt") %>% 
  mutate(well_id = as.character(well_id),
         mean_PE_A= 10^mean_PE_A) %>%
  dplyr::rename(Plate = plate)
cell_density <- read.table("../data/clean/cell_density.txt")
mean_FSC <- read.table("../data/clean/mean_FSC.txt")

###########################
#  Prep data for analysis #
###########################

#Assign strain info to FCS data 
#Remove Blanks and keep only stained data
mean_PEA <- right_join(mean_PEA, w96_key) %>%
  filter(Strain !="Blank",
         Stained == 1
  ) %>%
  mutate(Strain =as.factor(Strain)) %>%
  select(Strain, N_Starvation, mean_PE_A)

mean_PEA$N_Starvation <- factor(mean_PEA$N_Starvation, 
                                       levels = c("Before", "After"))

strain.to.region <- read.delim("../data/raw/strain to region.txt")
latitude_df<-strain.to.region %>% select (Strain, Latitude) %>%
  left_join(mean_PEA,.)
#====================================================================
# ANOVA: Comparison between Strains 
#     Look at Before and After N Starvation Treatments independently
#====================================================================
library(multcompView)
library(emmeans)
library(multcomp)
library(ggplot2)
#Compare strains before starvation

before.lm <- mean_PEA %>% filter(N_Starvation == "Before") %>%
  # filter(Strain !="GB117")%>%
  lm(log10(mean_PE_A)~Strain, data = .)
before.aov <- aov(before.lm)
summary(before.aov)
plot(before.aov)
hist(before.aov$residuals, breaks = 20)

Tukey<-TukeyHSD(x=before.aov, 'Strain', conf.level=0.95)
multcompLetters(extract_p(Tukey$Strain))
marginal<-lsmeans(before.aov,
                  ~ Strain)
tukey_output <- cld(marginal,
                    alpha=0.05,
                    Letters=letters,
                    adjust="sidak")
write.table(marginal,"../analysis/aov_tukey_before_model.txt", sep="\t", quote = FALSE)

#Compare strains after starvation
after.lm <- mean_PEA %>% filter(N_Starvation == "After") %>%
  # filter(Strain !="GB117")%>%
  lm(log10(mean_PE_A)~Strain, data = .)
after.aov <- aov(after.lm)
summary(after.aov)
plot(after.aov)
hist(after.aov$residuals, breaks = 20)

Tukey<-TukeyHSD(x=after.aov, 'Strain', conf.level=0.95)
multcompLetters(extract_p(Tukey$Strain))
marginal<-lsmeans(after.aov,
                  ~ Strain)
tukey_output <- cld(marginal,
                    alpha=0.05,
                    Letters=letters,
                    adjust="sidak")
write.table(tukey_output,"../analysis/aov_tukey_after_model.txt", sep="\t", quote = FALSE)

#Compare strains (whole data)
all.lm <- mean_PEA %>%
  # filter(Strain !="GB117")%>%
  lm(log10(mean_PE_A)~Strain*N_Starvation, data = .)
all.aov <- aov(all.lm)
summary(all.aov)
plot(all.aov)
hist(all.aov$residuals)

Tukey<-TukeyHSD(x=all.aov, 'Strain', conf.level=0.95)
multcompLetters(extract_p(Tukey$Strain))

marginal<-lsmeans(all.aov,
                   ~ Strain:N_Starvation)
tukey_output <- cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="sidak")
write.table(tukey_output,"../analysis/aov_tukey_whole_model.txt", sep="\t", quote = FALSE)

#====================================================================
#Plot mean and se as error bars
#====================================================================

df <- mean_PEA %>% group_by(Strain,N_Starvation) %>% 
  summarise( mean = mean(log10(mean_PE_A)), 
             se = sd(log10(mean_PE_A))/sqrt(n()), 
             max = mean+se, 
             min = mean - se) %>%
  mutate(mean = 10^mean,
         se = 10^se,
         max = 10^max,
         min = 10^min) %>%
  mutate(min = case_when( min <= 0 ~ 1,
                          TRUE ~ min)) %>%
  left_join(.,strain.to.region) %>% 
  na.omit(.) %>%  droplevels() %>%
  arrange(., N_Starvation, mean,Strain)

vjust<-rep(0.5)
df$group<-gsub(" ", "", tukey_output$.group)

color_by_region<-as.character(df$Color[1:26])  

p<-ggplot(data = df,aes( y = mean, x = factor(Strain, levels = df$Strain[1:26]))) + 
  geom_point(size = 1, aes(color = N_Starvation)) +
  geom_errorbar(aes(ymin = min, ymax = max, 
                    color = N_Starvation), width = 0, size = 0.5) +
  xlab("Strains") +
  ylab("Mean NR fluorescence") +
  labs(colour = "N Starvation") +
  #scale_color_manual(values=c("#9E4338", "#3B6C9D"), labels=c("Before","After"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #limits = c(10000, 10^6)) +
  # geom_text(hjust = -0.5,vjust=0.25,
  #           color   = "black") +
  theme_latex(base_size = 12) +
  theme(legend.position="top",
        axis.text.x = element_markdown(angle = 45, size =12, vjust = 1, hjust=1, 
                                       color= color_by_region),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

pdf("../analysis/NR fluorescence of all strains.pdf", height = 4, width = 6.5)
p
dev.off() 

#===================================================================
# Default dataframe: Effect of Latitude on mean PE A
#     Look at Before and After N Starvation Treatments independently
#====================================================================
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(lme4)
library(car)

#Fixed effects model
lm<-latitude_df %>% 
  lm(log10(mean_PE_A)~Latitude*N_Starvation,data =.)
summary(lm)
plot(lm)
hist(resid(lm))
latitude.aov<- Anova(lm, type = 2)
latitude.aov
#plot(latitude.aov)

#Mixed effects models
lm<-latitude_df %>% #filter(Strain!="GB117")%>% 
  lmer(log10(mean_PE_A)~Latitude*N_Starvation+(1|Strain),data =.)
summary(lm)
plot(lm)
hist(resid(lm), breaks = 15)
latitude.aov<- Anova(lm,type = 2)
latitude.aov

p1 <- latitude_df %>% filter(N_Starvation == "Before") %>%
  # filter(Strain!="GB117")%>%
  ggscatter(data = ., y = "mean_PE_A", x = "Latitude",
            size = 2,
            color = "#9E4338",
            alpha = 0.2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size = 5) +
  ylab("mean NR fluorescence") +
  ggtitle("Before N Starvation") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(colour = "black"))
p1

p2 <- latitude_df %>% filter(N_Starvation == "After") %>%
  # filter(Strain!="GB117")%>%
  ggscatter(data = ., y = "mean_PE_A", x = "Latitude",
            size = 2,
            color = "#3B6C9D",
            alpha = 0.2,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size = 5) +
  ylab("mean NR fluorescence") +
  ggtitle("After N Starvation") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(colour = "black"))
p2

P2 <- ggarrange(p1, p2,labels = c("A","B"),
                ncol = 2, nrow = 1)
P2
#-----------------
#Plot N starvations together
#-----------------
p1 <- latitude_df %>% 
  ggplot(., aes( y = mean_PE_A, x = Latitude)) +
  geom_point(size = 2, alpha = 0.2,
             aes(colour = N_Starvation)) +
  #scale_color_manual(values = c("#9E4338","#3B6C9D"), labels=c("Before","After")) +
  geom_smooth(method = "lm", se=F, aes(color = N_Starvation)) +
  ylab("NR fluorescence") +
  xlab(expression("Latitude ("*degree*"N)"))+
  labs(color = "N Starvation") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  #annotate("text", x = 28, y = 10^6, hjust = 0,
  #         label = "Latitude: p = 0.4699\nN Starvation: p << 0.0001\nLatitude\u00D7N Starvation: p = 0.4009",
  #         family = "lmroman") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
p1

pdf("../analysis/NR fluorescence across latitudes.pdf", height = 4, width = 3.5)
p1
dev.off() 

