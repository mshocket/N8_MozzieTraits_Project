#not sure what this should be for future use
#setwd()

#read in data
mdr_data <- read.csv("../data/MDR_Data.csv", stringsAsFactors = T)
pla_data <- read.csv("../data/pLA_Data.csv",stringsAsFactors = T)

names(mdr_data)
names(pla_data)

# mdr data has additional variables "includes.egg" "Increasing"

#make variable names consistent
colnames(pla_data)[which(names(pla_data) == "T")] <- "T.C"

mdr_vec <- as.vector(names(mdr_data))
pla_vec <- as.vector(names(pla_data))

#this should pop out any variables that are only in one dataset
union(setdiff(mdr_vec,pla_vec), setdiff(pla_vec,mdr_vec))
setdiff(mdr_vec,pla_vec) #this is just vars unique to mdr
setdiff(pla_vec,mdr_vec) #this is just vars unique to pla

#add additional vars to other dataset 
pla_data$includes.egg <- "NA"
pla_data$Increasing <- "NA"

mdr_and_pla <- rbind(mdr_data,pla_data)

#is there a better way of doing this final part?

