library(ggplot2)
library(dplyr)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(ltxplot)
load_theme_ltx()

Data <- read.delim("../data/clean/growth_curve_data.txt")
Data$Strain<-as.factor(Data$Strain)
Data$Plate<-as.factor(Data$Plate)
Data$Replicate<-as.factor(Data$Replicate)
Data$ABS650<-as.numeric(as.character(Data$ABS650))
Data$ABS750<-as.numeric(as.character(Data$ABS750))

#Calculate standard error
std <- function(x){
  sd(x)/sqrt(length(x))}

df_blank <- Data %>% filter(Strain=="Blank") %>%  group_by(Day,Plate,Replicate) %>%
  summarize(ABS650_blank=mean(ABS650), ABS750_blank=mean(ABS750))

df_temp<-full_join(Data,df_blank) %>% 
  mutate(Strain = case_when(
    Strain == "124" ~ "1373",
    TRUE ~ as.character(Strain)))

df_grouped_summary <- df_temp %>% 
  dplyr::filter(Strain!="Blank") %>%
  group_by(Strain,Day,Plate) %>%
  summarize(ABS650_std=std(ABS650),
            ABS750_std=std(ABS750),
            ABS650=mean(ABS650-ABS650_blank),
            ABS750=mean(ABS750-ABS750_blank)) %>%
  mutate(ABS650_std=replace(ABS650_std, is.na(ABS650_std),0),
         ABS750_std=replace(ABS750_std, is.na(ABS750_std),0),
         ABS650_min=ABS650-ABS650_std,
         ABS650_max=ABS650+ABS650_std,
         ABS750_min=ABS750-ABS750_std,
         ABS750_max=ABS750+ABS750_std,
         Strain = as.factor(Strain)) %>%
  arrange(Strain)
  

#grouped plot

p<-list()
count=0
for (strain in levels(df_grouped_summary$Strain) ){
  count<-count+1
  temp<-df_grouped_summary %>% filter(Strain==strain)
  p[[count]]<-ggplot(temp, aes(x = Day, y = ABS750, colour=Plate)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ABS750_min, ymax = ABS750_max)) +
    labs(title = strain)+
    scale_x_continuous(breaks=seq(0,7,1),limits=c(0,7))+ 
    scale_y_continuous(breaks=seq(-0.5,2,.5),limits=c(-0.5,2)) +
    theme_latex(base_size = 12)+
    theme(plot.title = element_text(size = 13), 
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank())
     #+
    # theme(legend.position = "none")
}
P1<- ggarrange(plotlist=p[1:12], labels = LETTERS[1:12],
          ncol = 3, nrow = 4, common.legend= TRUE)
P2 <- ggarrange(plotlist=p[13:count], labels = LETTERS[13:count],
                ncol = 3, nrow = 4, common.legend= TRUE)

pdf(file = '../analysis/growth_curve.pdf',  
    width = 6.5, 
    height = 9)
P1
P2
dev.off()
do.call(grid.arrange,p)

#Plot for one strain only
#p1<-df_grouped_summary %>% filter(Strain=="1691") %>%
#  ggplot(., aes(x = Day, y = ABS750, colour = Plate)) +
#  geom_point(size = 2) +
#  geom_errorbar(aes(ymin = ABS750_min, ymax = ABS750_max, width = 0)) +
#  labs(title = "CC-1691")+
#  scale_x_continuous(breaks=seq(0,8,1),limits=c(0,8))+ 
#  scale_y_continuous(breaks=seq(-0.5,2,.5),limits=c(-0.5,2)) +
#  theme_classic(base_size = 18) +
#  theme(legend.position = "right",
#        axis.text = element_text(colour = "black"))

