######Extract cell concentraion from fcs
library(fcs2R)
default_wd<-getwd()

Live_alga<- fcs_to_cell_density("../data/raw/fcs/", "../data/clean/Live_alga.txt")

#removes file extension ".fcs" from well id manually

######Assign plate layout to FACS reads and rename columns
FACS_plate<-c(3,4,5)
Plate_layout<-c("96-1","96-2","96-3")

library(dplyr)

Live_alga<-Live_alga %>% 
  mutate(Layout = case_when(plate==3~"96-1",
                            plate==4~ "96-2",
                            plate==5~"96-3")) 

w96layout<-read.table("../data/raw/96_Well_Plate_layout.txt", header = TRUE)

Live_alga<-merge(w96layout, Live_alga,c("Layout","well_id"))

##################################################################################
# Extract and compile 96 well abs information from multiple txt files in a folder#
##################################################################################

#Define work directory (where the files are at)
setwd("../data/raw/abs/")

#Define project directory (where you want to save output file)
dir<-"../../clean/"

out.file <- data.frame(matrix(ncol = 6, nrow = 0))
columns<-c("Plate_id","ABS650", "ABS750", "well_id", "Date", "plate_type")
colnames(out.file)<-columns
file.names <- dir(pattern ="Plate")
for(i in 1:length(file.names)){
  file<-file.names[i]
  temp_data<-read.delim(file, skip = 11,header=FALSE, col.names = c("Well","Read.1.650","Read.2.750"))
  Plate_id<-as.numeric(strsplit(strsplit(file, "_")[[1]][2]," ")[[1]][2])
  Date<- as.character(strsplit(strsplit(file, "_")[[1]][1]," ")[[1]][4])
  plate_type<-as.character(strsplit(strsplit(file, "_")[[1]][1]," ")[[1]][1])
  ABS650<-temp_data$Read.1.650
  ABS750<-temp_data$Read.2.750
  well_id<-as.character(temp_data$Well)
  df<-cbind(Plate_id, ABS650, ABS750, well_id, Date, plate_type)
  out.file<-rbind(out.file,df)
}

setwd(dir)
write.table(out.file, "96 well abs data.txt", quote=FALSE, sep="\t", row.names=FALSE)

w96abs_data<- out.file

w96abs_data<-w96abs_data %>% 
  mutate(Layout = case_when(Plate_id==1~"96-1",
                            Plate_id==2~"96-2",
                            Plate_id==3~"96-3"))

#####Standardize plates

w96_final<-left_join(w96layout, w96abs_data,c("Layout","well_id"))

w96_final$ABS650<-as.numeric(as.character(w96_final$ABS650))
w96_final$ABS750<-as.numeric(as.character(w96_final$ABS750))
temp<-w96_final %>% filter (Strain=="Blank") %>% group_by(Plate_id) %>%
  summarize(ABS650_blank=mean(ABS650), ABS750_blank=mean(ABS750))
w96_final<-full_join(w96_final, temp)
w96_final<-w96_final %>% mutate(ABS650_96=ABS650-ABS650_blank, ABS750_96=ABS750-ABS750_blank)

####Combine 96 well plate abs data and cell density data

w96layout<-read.table("../raw/96_Well_Plate_layout.txt", header=TRUE)

library(purrr)

cell_density_abs_data<-merge(Live_alga,w96_final, 
                             by=c("Strain",
                                  "Strain_replicate",
                                  "Layout", 
                                  "Concentration",
                                  "well_id")) %>% 
  droplevels(.) %>% 
  filter (Strain != "Blank")
  
write.table(cell_density_abs_data, "../clean/cell_density_abs_data.txt", quote=FALSE, sep="\t", row.names=FALSE)
#data<-select(cellsity_abs_data,Strain,cell_density,ABS650_48,ABS650_96,ABS750_48,ABS750_96)

####Plot 96 abs against cell density
data<-read.table("../clean/cell_density_abs_data.txt", header=TRUE)

attach(data)
params_650_96 <- data.frame(matrix(ncol = 4, nrow = 0))
params_750_96 <- data.frame(matrix(ncol = 4, nrow = 0))
Strains<-data %>% filter (Strain != "Blank")
Strains<-levels(Strains$Strain)
data$cell_density<-as.numeric(as.character(data$cell_density))
par(mfrow=c(5,6))
for (strain in Strains){
  
  #create cell density v abs750 plots for each strain, 
  #run lm model, 
  #save model parameters
  
  Temp<- data[ which(data$Strain== strain), ]
  #plot(Temp$cell_density, Temp$ABS750_96, type='p', main=strain, xlab="Number of cells per uL", ylab="ABS750_96")
  model<- lm(Temp$ABS750_96 ~ Temp$cell_density, data=Temp)
  output<-cbind(strain, as.numeric(model$coefficients[2]), 
                as.numeric(model$coefficients[1]), as.numeric(summary(model)$r.squared))
  params_750_96<-rbind(params_750_96,output)
  
  Temp<- data[ which(data$Strain== strain), ]
  #plot(Temp$cell_density, Temp$ABS650_96, type='p', main=strain, xlab="Number of cells per uL", ylab="ABS650_96")
  model<- lm(Temp$ABS650_96 ~ Temp$cell_density, data=Temp)
  output<-cbind(strain, as.numeric(model$coefficients[2]), 
                as.numeric(model$coefficients[1]), as.numeric(summary(model)$r.squared))
  params_650_96<-rbind(params_650_96,output)
  
}

columns<-c("strain","slope","intercept","r_squared")
colnames(params_750_96)<-columns
colnames(params_650_96)<-columns
write.table(params_750_96, "abs750_96_model_parameters.txt", sep="\t")
write.table(params_650_96, "abs650_96_model_parameters.txt", sep="\t")
