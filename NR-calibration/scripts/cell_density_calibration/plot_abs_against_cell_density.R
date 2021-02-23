####Plot 96 abs against cell density
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(ltxplot)
load_theme_ltx()

params_750_96 <- read.delim("../data/clean/abs750_96_model_parameters.txt")

p<-list()
count=0
for (strain in Strains ){
  count<-count+1
  temp<- data %>% filter(Strain==strain)
  p[[count]]<-ggscatter(data = temp , y = "ABS750_96", x = "cell_density",
              size = 2,
              add = "reg.line", 
              add.params = list(color = "black"), 
              conf.int = FALSE) +
    stat_cor(method = 'pearson', label.sep = "\n", 
             size = 4, label.x.npc = 0.4, label.y.npc = 0.3, family = 'lmroman')+
    labs(title = strain, y = "ABS750", x = "Cell density (cells/uL)")+
    scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))+ 
    scale_x_continuous(breaks=seq(0,11000,5000),limits=c(0,11000)) +
    theme_latex() +
    theme( title = element_text(size = 11), 
           panel.grid.major.x = element_blank(),
           panel.grid.major.y = element_blank())
}

P1<- ggarrange(plotlist=p[1:12], labels = LETTERS[1:12],
               ncol = 3, nrow = 4, common.legend= TRUE)
P2 <- ggarrange(plotlist=p[13:count], labels = LETTERS[13:count],
                ncol = 3, nrow = 4, common.legend= TRUE)

pdf(file = '../analysis/cell_density_calibration.pdf',  
    width = 6.5, 
    height = 9)
P1
P2
dev.off()

#plot CC-1691
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggpmisc)

p1 <- data %>% filter(Strain == "1691") %>%
  ggscatter(data = . , y = "ABS750_96", x = "cell_density",
            size = 4,
            add = "reg.line", 
            add.params = list(color = "black"), 
            conf.int = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size = 5, label.x = 24, family = 'lmroman') +
  labs(title = "CC-1691", y = "ABS750", x = "Cell density (cells/uL)") +
  # scale_y_continuous(breaks = seq(0,1,0.2),limits=c(0,1)) + 
  # scale_x_continuous(breaks = seq(0,11000,5000),limits = c(0,11000)) +
  theme_latex() + 
  theme(title = element_text(size = 13),  
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank())+
  theme(axis.text = element_text(colour = "black"))

pdf(file = '../analysis/CC-1691.pdf',  
    width = 4, 
    height = 4)
p1
dev.off()

#plot histogram of R-squared
params_750_96 <- params_750_96 %>% filter (strain != "Blank")
med <- median(params_750_96$r_squared)

p2 <- ggplot(data=params_750_96, aes(r_squared)) + 
  geom_histogram(bins = 10, color = "black", fill = "white") +
  geom_vline(xintercept = med, linetype = "dashed", color = "blue") +
  xlab(expression(R^2)) +
  ylab("Frequency") +
  annotate("text", x = 0.87, y = 5, 
           label = expression(paste("Median ", R^2, "= 0.9303")),
           size = 5, family = 'lmroman')+
  theme_latex() + 
  theme(title = element_text(size = 13),  
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank())+
  theme(axis.text = element_text(colour = "black"))

pdf(file = '../analysis/histogram_rsquared.pdf',  
    width = 4, 
    height = 4)
p2
dev.off()
