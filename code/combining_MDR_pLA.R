#not sure what this should be for future use
#setwd()

#read in data
mdr_data <- read.csv("../data/MDR_Data.csv", stringsAsFactors = T)
pla_data <- read.csv("../data/pLA_Data.csv",stringsAsFactors = T)

#convert all dev rates to days 
for(i in 1:nrow(mdr_data)){
  
  if(mdr_data$trait.name[i] == "1/MDR"){
    mdr_data$trait[i] <- 1/mdr_data$trait[i]
  }
}
summary(mdr_data)

#change trait name to reflect conversion
mdr_data$trait.name <- as.factor("MDR")
summary(mdr_data)

#see var names for both datasets
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
# pla_data$includes.egg <- "NA"
# pla_data$Increasing <- "NA"

#remove increasing and egg cols (egg stage sometimes included in dev time, increasing denotes portion before turnover)
mdr_data <- subset(mdr_data, select = -c(includes.egg,Increasing))
names(mdr_data)

mdr_and_pla <- rbind(mdr_data,pla_data)

summary(mdr_and_pla)
getwd()

#write.csv(mdr_and_pla,file='../data/MDR_and_pLA_Data.csv', row.names=FALSE)


