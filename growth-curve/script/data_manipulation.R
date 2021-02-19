library(dplyr)

##########################################################################
# Extract and compile abs information from multiple txt files in a folder#
##########################################################################

#Default directory
default_wd<-getwd()
#Define where the absorbance files are at)
setwd("../data/raw/absorbance/")

out.file <- data.frame(matrix(ncol = 6, nrow = 0))
columns<-c("Plate_id","ABS650", "ABS750", "Well_id", "File_create_date", "Plate_type")
colnames(out.file)<-columns
file.names <- dir(pattern ="Plate")
for(i in 1:length(file.names)){
  file<-file.names[i]
  temp_data<-read.delim(file, skip = 10)
  Plate_id<-as.numeric(strsplit(strsplit(file, "_")[[1]][2]," ")[[1]][2])
  Date<- as.character(strsplit(strsplit(file, "_")[[1]][1]," ")[[1]][4])
  plate_type<-as.character(strsplit(strsplit(file, "_")[[1]][1]," ")[[1]][1])
  ABS650<-temp_data$Read.1.650
  ABS750<-temp_data$Read.2.750
  Well_id<-as.character(temp_data$Well)
  df<-cbind(Plate_id, ABS650, ABS750, Well_id, Date, plate_type)
  out.file<-rbind(out.file,df)
}

setwd(default_wd)

write.table(out.file, "../data/clean/Growth curve data.txt", quote=FALSE, sep="\t", row.names=FALSE)

####################################################
# Extract unique plates and assign plate layout id #
####################################################
number_of_48_well_plates=3

out.file$Plate_id<-as.numeric(levels(out.file$Plate_id))[out.file$Plate_id]

unique_plates<- out.file %>% 
  group_by(Plate_id,Date,plate_type) %>% 
  select (-c(ABS650, ABS750, Well_id)) %>% 
  tally() %>%
  mutate(Plate = case_when(plate_type == "96" ~ "96-well-1",
                           plate_type == "48" && Plate_id<(number_of_48_well_plates*2+1) ~ "48-well-1",
                           TRUE ~ "48-well-2")) %>%
  
  
  #Change the numbers below for different number of days for growth curve and different number of 48 well plates
  
  
  mutate(Plate_layout = case_when(Plate=="96-well-1" ~ "96-1",
                            plate_type== "48" && Plate_id %in% seq(1,52,3) ~ "48-1",
                            plate_type== "48" && Plate_id %in% seq(2,52,3) ~ "48-2",
                            plate_type== "48" && Plate_id %in% seq(3,52,3) ~ "48-3")) %>%
  
  mutate(Day = case_when(Plate_id==1 && Plate=="96-well-1" ~ 0,
                         Plate_id==2 && Plate=="96-well-1" ~ 1,
                         Plate_id==3 && Plate=="96-well-1" ~ 2,
                         plate_type== "48" && Plate_id %in% seq(1,3) ~ 2,
                         plate_type== "48" && Plate_id %in% seq(4,6) ~ 4,
                         plate_type== "48" && Plate_id %in% seq(7,9) ~ 4,
                         plate_type== "48" && Plate_id %in% seq(10,12) ~ 5,
                         plate_type== "48" && Plate_id %in% seq(13,15) ~ 6))


#####################################################################
# Import plate layout: id strain for each well in each plate layout #
#####################################################################

#Import table format
#Plate layout id, well id, strain/BLANK

plate_layout_df <- read.delim("../data/raw/07242019.txt")

##################################################################################
# Assign strain to each used well in master data file & remove unecessary columns#
##################################################################################

df<-left_join(out.file,unique_plates, by=c("Plate_id", "Date","plate_type")) %>% 
  select(-c(Date, plate_type, n, Plate_id)) %>%
  mutate(Plate_layout = as.factor(Plate_layout)) %>%
  merge(.,plate_layout_df) %>%
  rename(Replicate=Plate_layout) %>% 
  select (-c(Well_id))

write.table(df, "../data/clean/growth_curve_data.txt", quote=FALSE, sep="\t", row.names=FALSE)




