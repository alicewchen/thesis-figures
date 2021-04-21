library(readxl)
library(reshape)
library(dplyr)
library(stringi)
library(fcs2R)
library(ggplot2)
library(scales)
library(ggtext)
library(ltxplot)
library(ggpubr)
library(lme4)
library(car)
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
# Figure 4 and S: Variation of NR fluorescence among strains
#====================================================================
model <- readRDS("../../NR-calibration/data/clean/model.rds")
predicted_NL<- latitude_df %>% mutate(predicted_NL = predict(model, newdata = .))
df <- predicted_NL %>% group_by(Strain,N_Starvation) %>% 
  summarise( mean = mean(predicted_NL), 
             se = sd(predicted_NL)/sqrt(n()), 
             max = mean+se, 
             min = mean - se) %>%
  mutate(mean = case_when(mean <= 0 ~0,
                          TRUE ~ mean),
         min = case_when( min <= 0 ~ 0,
                          TRUE ~ min),
         max = case_when(max<=0 ~ 0,
                         TRUE ~ max)
         ) %>%
  left_join(.,strain.to.region) %>% 
  na.omit(.) %>%  droplevels() %>%
  arrange(., N_Starvation, mean,Strain)

vjust<-rep(0.5)
df$group<-gsub(" ", "", tukey_output$.group)

color_by_region<-as.character(df$Color[1:26])  

p<-ggplot(data = df,aes( y = mean, x = factor(Strain, levels = df$Strain[1:26]))) + 
  geom_point(size = 1, aes(color = N_Starvation)) +
  geom_errorbar(aes(ymin = min, ymax = max, 
                    color = N_Starvation), width = 0, size = 0.4) +
  labs(x = "Strains",
       y=expression( paste("Neutral lipid (", mu, "g ",10^-6," cells)")),
       color = "N-starvation") +
  theme_latex(base_size = 12) +
  theme(legend.position="top",
        axis.text.x = element_markdown(angle = 45, size =12, vjust = 1, hjust=1, 
                                       color= color_by_region),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())
pdf("../analysis/NL of all strains.pdf", height = 4, width = 6.5)
p
dev.off() 
# Supplemental figure ####
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

p<-ggplot(data = df,aes( y = mean, x = factor(Strain, levels = df$Strain[1:26]))) + 
  geom_point(size = 1, aes(color = N_Starvation)) +
  geom_errorbar(aes(ymin = min, ymax = max, 
                    color = N_Starvation), width = 0, size = 0.4) +
  labs(x = "Strains",
       y="NR fluorescence",
       color = "N-starvation") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
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
# Effect of Latitude on mean NR fluorescence
#====================================================================
# Fixed effects model ####
sink("../analysis/lm_model_output.txt", type = "output")
lm<-latitude_df %>% 
  lm(log10(mean_PE_A)~Latitude*N_Starvation,data =.)
summary(lm)
latitude.aov<- Anova(lm, type = 2)
latitude.aov
sink(NULL)

pdf("../analysis/LatxStarvation_lm_plots.pdf", height = 7.5, width = 11)
par(mfrow=c(2,3))
plot(lm)
hist(resid(lm))
dev.off() 

# Mixed effects models #####
sink("../analysis/lmer_model_output.txt", type = "output")
lm<-latitude_df %>% #filter(Strain!="GB117")%>% 
  lmer(log10(mean_PE_A)~Latitude*N_Starvation+(1|Strain),data =.)
summary(lm)
latitude.aov<- Anova(lm,type = 2)
latitude.aov
sink(NULL)
pdf("../analysis/LatxStarvation_lmer_plots.pdf", height = 7.5, width = 11)
par(mfrow=c(1,1))
plot(lm)
hist(resid(lm), breaks = 15)
dev.off() 

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

# Resampling ####
library(purrr)
library(broom)
set.seed(1)
#Function for random stratified resampling (balanced)
strain_key<-strain.to.region %>% #filter(Strain!="GB117") %>% 
  dplyr::select(Strain, Latitude)
resampling <- function (df, nruns){
  output<-list()
  for (i in 1:nruns){
    set.seed(i)
    sample <- strain_key %>% group_by(Latitude) %>% sample_n(., 1)
    sample <- df %>% right_join(.,sample, by=c("Strain", "Latitude"))
    output <- append(output,list(sample))  
  }
  return(output)
}

#Function for running lm on each sample
run_lm <- function(x) {
  mod <- lm(log10(mean_PE_A) ~ Latitude+N_Starvation, x)
  r_sq<-summary(mod)$r.squared
  tibble <- as.data.frame(t(broom::tidy(mod)))[2,]
  tibble<-cbind(tibble,r_sq)
  colnames(tibble)<- c(as.data.frame(broom::tidy(mod))[,1], "R squared")
  return(tibble)
} 

run_lm_interaction <- function(x) {
  mod <- lm(log10(mean_PE_A) ~ Latitude*N_Starvation, x)
  tibble <- as.data.frame(t(broom::tidy(mod)))[2,]
  colnames(tibble)<- as.data.frame(broom::tidy(mod))[,1]
  return(tibble)
} 

to_numeric <- function(x){
  y<- lapply(x, as.character)
  y <- lapply(y, as.numeric)
  y<- as.data.frame(y)
  return(y)
}
#test ####
sample <- strain_key %>% group_by(Latitude) %>% sample_n(., 1)
sample <- latitude_df %>% #filter(Strain!="GB117") %>% 
  right_join(.,sample, by=c("Strain", "Latitude"))

mod <- lm(log10(mean_PE_A) ~ Latitude+N_Starvation, sample)
print(t(broom::tidy(mod)))
summary(mod)
tibble <- as.data.frame(t(broom::tidy(mod)))[2,]
colnames(tibble) <-c(as.data.frame(broom::tidy(mod))[,1])
print(tibble)
#end of test --------------

resample <- latitude_df %>% #filter(Strain!="GB117") %>%
  resampling(.,10000)
output <- to_numeric(map_df(resample[1:10000], run_lm))
output2 <- to_numeric(map_df(resample[1:10000], run_lm_interaction))
#utput3 <- map_df(resample[1:100], run_lm)
#names(output)<- col_names

#create summary of mean, std, CI for simulation statistics
resample_summary <- function (df){
  # for (i in 1:ncol(df)){
  #   df[,i] <- as.numeric(df[,i])
  # }
  col_names <- names(df)
  mean<-sapply(df, mean, na.rm=TRUE)
  std<-sapply(df, sd, na.rm=TRUE)
  lower_CI<-sapply(df,function (x) {quantile(x,.025)})
  upper_CI<-sapply(df,function (x) {quantile(x,.975)})
  p_value <- 2*pnorm(-abs((0-mean)/std))
  more_than <- sapply(df, function(x){sum(x>0)})
  less_than <- sapply(df, function(x){sum(x<0)})
  output_summary <- as.data.frame(cbind(mean, std, lower_CI, upper_CI, p_value, more_than, less_than))
  rownames(output_summary)<-col_names
  return(output_summary)
}

output_sum <- resample_summary(output)
output2_sum <- resample_summary(output2)

#Extract slope and intercept for Before/After of each resample
count = 0
merged_resamples <- data.frame(matrix(ncol = 9, nrow = 0))
for (i in resample[1:10000]) {
  count <- count + 1
  i<-cbind(i, count)
  merged_resamples <- rbind(merged_resamples, i)
}

# Plot histogram of slopes ####
p1<-ggplot(output2, aes(x=`Latitude.N_StarvationAfter`)) + 
  geom_histogram(color="black", fill="white", bins=15)+
  xlab("Slope estimate") +
  ylab("Frequency") +
  ggtitle("Latitude*N Starvation")+
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

#shapiro.test(output2$`Latitude:N_StarvationAfter`)
p2<-ggplot(output2, aes(x=Latitude)) + 
  # geom_histogram(color="black", fill="white", bins=15)+
  #geom_density(aes(x=Latitude,
  #                  y=..density.., alpha = 0.2)) +
  geom_histogram(color="black", fill="white", bins=15)+
  geom_vline(aes(xintercept=output2_sum[["Latitude","mean"]], linetype="dashed",
                 colour = "blue"))+
  xlab("Slope estimate") +
  ylab("Count") +
  ggtitle("Latitude")+
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'),
        legend.position = "none")
pdf("../analysis/histogram of latitude slopes.pdf", height = 4, width = 6.5)
p2
dev.off() 

p3<-ggplot(output2, aes(x=N_StarvationAfter)) + 
  # geom_histogram(color="black", fill="white", bins=15)+
  # geom_density(aes(x=N_StarvationAfter,
  #                  y=..density.., fill = "grey", alpha = 0.2)) +
  geom_histogram(color="black", fill="white", bins=15)+
  #scale_x_continuous(labels = unit_format(unit = "K", scale = 0.001)) +
  xlab("Slope estimate") +
  ylab("Frequency") +
  ggtitle("N Starvation")+
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

P3 <- ggarrange(p1,p2,p3,labels=c("A","B","C"), nrow = 2, ncol = 2)
P3

#===============================================================================
# Figure 5A: NR fluorescence across Latitudes grouped by N-starvation treatment
#===============================================================================
predict_df <- latitude_df %>% cbind(., predicted_NL) %>%
  mutate(compressed_predicted_NL = case_when( #compress negative NL values to 0
    predicted_NL< 0 ~ 0,
    predicted_NL >= 0 ~ predicted_NL
  ))

p1 <- latitude_df %>% 
  ggplot(., aes( y = mean_PE_A, x = Latitude)) +
  geom_point(size = 2, alpha = 0.2,
             aes(colour = N_Starvation)) +
  geom_smooth(method = "lm", se=F, aes(color = N_Starvation),size = 0.5) +
  labs(x = expression("Latitude ("*degree*"N)"),
       y = "NR fluorescence",
       color = "N-starvation",
       tag = "A") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))


#===============================================================================
# Figure 5B: NR fluorescence across latitudes (downsampled)
#===============================================================================
p2<- ggplot(merged_resamples, 
       aes(y = mean_PE_A, x = Latitude, color = N_Starvation,
           group = interaction(N_Starvation, count))) +
  geom_line(stat="smooth",method = "lm", formula = y ~x, alpha = 0.01, se = TRUE) +
  #scale_color_manual(values = c("#D55E00", "#0072B2"), labels=c("Before","After"), ) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x = expression("Latitude ("*degree*"N)"),
       y = "NR fluorescence",
       color = "N-starvation",
       tag = "B") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))


P <- ggarrange(p1,p2,nrow = 1, ncol = 2, common.legend = T)
pdf("../analysis/NR across latitudes.pdf", height = 4, width = 6.5)
P
dev.off()

#===============================================================================
# Figure S: NL across latitudes
#===============================================================================
p3 <- predict_df %>% 
  ggplot(., aes( y = compressed_predicted_NL, x = Latitude)) +
  geom_point(size = 2, alpha = 0.2,
             aes(colour = N_Starvation)) +
  geom_smooth(method = "lm", se=F, aes(color = N_Starvation),size = 0.5) +
  labs(x = expression("Latitude ("*degree*"N)"),
       y = expression( paste("Neutral lipid (", mu, "g ",10^-6," cells)")),
       color = "N-starvation",
       tag = "A") +
  theme_latex(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))

pdf("../analysis/NL across latitudes.pdf", height = 4.5, width = 4.5)
p3
dev.off()
