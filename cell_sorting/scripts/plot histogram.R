library(flowCore)
library(dplyr)
library(ggplot2)
library(ltxplot)
library(scales)
library(ggpubr)
load_theme_ltx()
ls <- c()
df <- data.frame(sample = as.character(),
                 val = as.numeric())
f.names <- dir("../data/raw/fcs/",pattern=".fcs")
for (i in seq(1,3)){
  f <- read.FCS(paste("../data/raw/fcs/", f.names[i], sep = ""))
  #ls[[f.names[i]]]<- list(f@exprs[,40])
  temp.df<- data.frame(val = f@exprs[,40], sample = f.names[i])
  df<- rbind(df, temp.df)
}
detach("package:flowCore")
levels(df$sample) <- c("CC-1952","GB13", "Recombinants")

summary.df <- df %>% group_by(sample) %>%
  summarize(mean = mean(val),
         variance = sd(val)^2)

p<-ggplot(df, aes(y = ..count.., x = val, color = sample)) +
  geom_density() +
  xlab(expression("log"[10]*"NR fluorescence:FSC-A*10K"))+
  ylab("Number of cells")+
  labs(fill = "Sample")+
  theme_latex(base_size = 12)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

pdf("../analysis/histogram.pdf", height = 4, width = 6.5)
p
dev.off()

#==============================================
# Castle-Wright Estimator for haploids
#==============================================
cw <- ((0.7794627-0.3690399)^2-0.012597381-0.007558088)/(4*0.006751404)
cw

#========================
# Figure S3
#========================
f <- read.FCS("../data/raw/fcs/parent 2_001_Center.fcs")
GB13_df <- data.frame(PE_A = as.numeric(f@exprs[,23])) %>% mutate(log_PEA = log10(PE_A))
p1<-ggplot(GB13_df, aes(y = ..count.., x = PE_A)) +
  geom_density() +
  labs(x="NR fluorescence",
       y = "Number of cells",
       tag = "A")+
  xlim(0,10^1.5)+
  theme_latex(base_size = 12)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))
p2<-ggplot(GB13_df, aes(y = ..count.., x = PE_A)) +
  geom_density() +
  labs(x="NR fluorescence",
       y = "Number of cells",
       tag = "B")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_latex(base_size = 12)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.tag = element_text(face = 'bold'))
P<- ggarrange(p1,p2, nrow = 1, ncol = 2)

pdf("../analysis/Figure S3.pdf", height = 3, width = 6.5)
P
dev.off()

