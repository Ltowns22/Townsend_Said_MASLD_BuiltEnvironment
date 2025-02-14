# Townsend & Said 2024. 
# R code for incorporating features of the built enviornment associated with subject zip codes. 
# to ensure patient privacy only a file with the list of the subject zipcodes is included here. 

library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(zipcodeR)
library(geosphere)

setwd("https://github.com/Ltowns22/Townsend_Said_MASLD_BuiltEnvironment")

# function to get the mode values for the zipcodes with only 5 digits 
MedianADI_function <- function(Zip5_list, ADI.df){
  Zip5Sum.df <- data.frame(matrix(ncol = 5, nrow =0, dimnames=list(NULL, c("Zip_4", "GISJOIN","ADI_2020state","ADI_2020national","FIPS_2020"))))
  ADI.df$Zip5<- substr(ADI.df$ZIP_4, 1,5)
  zipSubset.df <- as.data.table(filter(ADI.df, Zip5 %in% Zip5_list))
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  for(z in zipSubset.df$Zip5){
    if(z %in% Zip5_list){
      if(!z %in% Zip5Sum.df$ZIP_4){
        zipSmallSubset.df <- as.data.table(filter(zipSubset.df, Zip5 ==z))
        Sum.df <- zipSmallSubset.df[, list("ZIP_4" = Zip5,
                                           "GISJOIN" = getmode(GISJOIN),
                                           "FIPS" = getmode(FIPS),
                                           "ADI_NATRANK"= getmode(ADI_NATRANK),
                                           "ADI_STATERANK" = getmode(ADI_STATERANK),
                                           "TYPE" = getmode(TYPE))]
        Sum.df<-as.data.frame(Sum.df)
        Zip5Sum.df <- rbind(Zip5Sum.df, Sum.df)
        print(z)
      }
    }
  }
  print(unique(Zip5Sum.df))
  return(unique(Zip5Sum.df))
}

# Input data frames - combining the ADI data with the patient data 
NAFLD.df <- read_csv("./NAFLD_ZipCodes.csv")
NAFLD.df<-unique(NAFLD.df )
NAFLD.df$State[NAFLD.df$ZIP == 55974] <- "MN"
ZipALL<- NAFLD.df$ZIP_4
Zip5<- c("53562", "54656", "53538", "53527", "54449", "53190", "53963", "53548", "55974", "54498") # list of the zipcodes in the dataframe with only 5 digits 

WI_2020_ADI_9DigitZipCode_v3_2<-read.csv("WI_2020_ADI_9DigitZipCode_v3.2.csv")
MN_2020_ADI_9DigitZip_Code_v3_2 <-read.csv("MN_2020_ADI_9DigitZipCode_v3.2.csv")
IL_2020_ADI_9DigitZipCode_v3_2<-read.csv("IL_2020_ADI_9DigitZipCode_v3.2.csv")
CO_2020_ADI_9DigitZipCode_v3_2<-read.csv("CO_2020_ADI_9DigitZipCode_v3.2.csv")
FL_2020_ADI_9DigitZipCode_v3_2 <-read.csv("FL_2020_ADI_9DigitZipCode_v3.2.csv")

WI_ADI <- as.data.frame(filter(WI_2020_ADI_9DigitZipCode_v3_2, ZIP_4 %in% ZipALL))
WI5zip_ADI <- WI_2020_ADI_9DigitZipCode_v3_2
WI5zip_ADI$Zip_5<- substr(WI5zip_ADI$ZIP_4, 1,5)
WI5zip_ADI <- filter(WI5zip_ADI, Zip_5 %in% Zip5)
WI_5zip.df <- as.data.frame(MedianADI_function(Zip5, WI5zip_ADI))

MN_ADI <- as.data.frame(filter(MN_2020_ADI_9DigitZip_Code_v3_2, ZIP_4 %in% ZipALL))
MN5zip_ADI <- MN_2020_ADI_9DigitZip_Code_v3_2
MN5zip_ADI$Zip_5<- substr(MN5zip_ADI$ZIP_4, 1,5)
MN5zip_ADI <- filter(MN5zip_ADI, Zip_5 %in% Zip5)
MN_5zip.df <- as.data.frame(MedianADI_function(Zip5, MN5zip_ADI))

IL_ADI <- as.data.frame(filter(IL_2020_ADI_9DigitZipCode_v3_2, ZIP_4 %in% ZipALL))
CO_ADI <- as.data.frame(filter(CO_2020_ADI_9DigitZipCode_v3_2, ZIP_4 %in% ZipALL))
FL_ADI <- as.data.frame(filter(FL_2020_ADI_9DigitZipCode_v3_2, ZIP_4 %in% ZipALL))

ADI_allStates <- rbind(WI_ADI, WI_5zip.df)
ADI_allStates <- rbind(ADI_allStates, MN_ADI)
ADI_allStates <- rbind(ADI_allStates, MN_5zip.df)
ADI_allStates <- rbind(ADI_allStates, IL_ADI)
ADI_allStates <- rbind(ADI_allStates, CO_ADI)
ADI_allStates <- rbind(ADI_allStates, FL_ADI)
ADI_allSates<- ADI_allStates %>% rename("ADI_2020state" = "ADI_STATERANK",
                                        "ADI_2020national" = "ADI_NATRANK",
                                        "FIPS_2020" = "FIPS",
                                        "ZIPcode_Type"= "TYPE")

NAFLD_ADI.df <- merge(NAFLD.df, ADI_allStates, by = "ZIP_4", all =T )

# Combining the NAFLD ADI dataframe with the EPA dataframe 
EPA.df <- as.data.frame(read.csv("EPA_SmartLocationDatabase_V3_Jan_2021_Final.csv"))
KeyStateCodes <- c("55", "17", "8","27", "12", "41") # WI, Il, CO, MN, FL, OR
EPA.keyStates <- subset.data.frame(EPA.df, STATEFP %in% KeyStateCodes)
EPA.keyStates$STATEFP<- str_pad(EPA.keyStates$STATEFP, 2, pad = "0") # adding 0 to make "8" into "08" so all are 2 digits with a leading 0 if needed 
EPA.keyStates$COUNTYFP<- str_pad(EPA.keyStates$COUNTYFP, 3, pad = "0")
EPA.keyStates$TRACTCE<- str_pad(EPA.keyStates$TRACTCE, 6, pad = "0")
EPA.keyStates$GISJOIN2 <-paste("G",EPA.keyStates$STATEFP,"0",EPA.keyStates$COUNTYFP,"0",EPA.keyStates$TRACTCE,EPA.keyStates$BLKGRPCE, sep="")
head(EPA.keyStates$GISJOIN2)

NAFLD_ADI.df$GISJOIN2 <-NAFLD_ADI.df$GISJOIN 
NAFLD_ADI.df <- NAFLD_ADI.df %>% mutate(GISJOIN2 = case_when(GISJOIN2 == "G55011100001041" ~ "G55011100001004", # 8 of the GISjoins in the NAFLD dataset are NOt in the EPA dataset - so replacing these with the closest geographic equivelent GISjoin 
                                                             GISJOIN2  == "G55002500120033" ~ "G55002500120021",
                                                             GISJOIN2  == "G55007809401041" ~ "G55007809401011",
                                                             GISJOIN2  == "G55002500006004" ~ "G55002500006002",
                                                             GISJOIN2  == "G55002500120031" ~ "G55002500120021",
                                                             GISJOIN2  == "G12007100502131" ~ "G12007100502071",
                                                             GISJOIN2  == "G55010500013061" ~ "G55010500013041",
                                                             GISJOIN2  == "G55005501001005" ~ "G55005501001004",
                                                             GISJOIN2  == "G55002500004062" ~ "G55002500004061",
                                                             T ~ GISJOIN2))
NAFLD_ADI.df <- NAFLD_ADI.df %>% mutate(ATLAS_EPA_GIS_MATCH = if_else(GISJOIN == GISJOIN2, "Match", "No"))

NAFLD_ADI.df <- merge(NAFLD_ADI.df, EPA.keyStates, by = "GISJOIN2", all =F )
write.csv(NAFLD_ADI.df, "./NAFLD_ZIP_BuiltEnviornment_Merged.csv")

# distances to hospitals 
# UW_health zip = 53792
# base code for calculating distance (in meters between two zips)
# list of transplant centers: https://optn.transplant.hrsa.gov/about/search-membership/?memberType=Transplant+Centers&state=-1&region=0 
# filtered for states immediatly surrounding states of pt in the data set, filtered for liver transplant (listed as active), and removed childrens hospitals. 
# = 44 key transplant centers that could be close to pt
df <- data.frame("ZIP_START" = c(95051, 94534, 60193, 94591, 94128, 94015, 94553, 10994, 95008), 
                 "ZIP_END" = c(98053, 94128, 60666, 73344, 94128, 73344, 94128, "07105", 94128), 
                 stringsAsFactors = FALSE)

df$distance_meters<-apply(df, 1, function(x){
  startindex<-which(x[["ZIP_START"]]==zip_code_db$zipcode)
  endindex<-which(x[["ZIP_END"]]==zip_code_db$zipcode)
  distGeo(p1=c(zip_code_db[startindex, "lng"], 
               zip_code_db[startindex, "lat"]), 
          p2=c(zip_code_db[endindex, "lng"], 
               zip_code_db[endindex, "lat"]))
})

df <- df %>% mutate (distance_miles = (distance_meters /1609.34))

## distances
NAFLD.df <- NAFLD_ADI.df 
LiverTxCenterZips.df <- as.data.frame(read_csv("./LiverTxCenterZips.csv"))

NAFLD.df$DistanceTo_UWhospital_meters<-apply(NAFLD.df, 1, function(x){
  startindex<-which(x[["Zip_short"]]==zip_code_db$zipcode)
  endindex<- ("53792"==zip_code_db$zipcode)
  distGeo(p1=c(zip_code_db[startindex, "lng"], 
               zip_code_db[startindex, "lat"]), 
          p2=c(zip_code_db[endindex, "lng"], 
               zip_code_db[endindex, "lat"]))
})

NAFLD.df <- NAFLD.df %>% mutate (DistanceTo_UWhospital_miles = (DistanceTo_UWhospital_meters /1609.34))

UniTxZips <- unique(LiverTxCenterZips.df$Zip)

Pt.Zip.df <- NAFLD.df %>% select(c("PtID", "Zip_short"))
full.df<- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Pt.Zip", "Hos.Zip", "Distance"))))
for(pZip in Pt.Zip.df$Zip_short){
  dist.df<- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Pt.Zip", "Hos.Zip", "Distance"))))
  for(hZip in UniTxZips){
    startindex<-which(pZip==zip_code_db$zipcode)
    endindex<- (hZip==zip_code_db$zipcode)
    distance = (distGeo(p1=c(zip_code_db[startindex, "lng"], 
                             zip_code_db[startindex, "lat"]), 
                        p2=c(zip_code_db[endindex, "lng"], 
                             zip_code_db[endindex, "lat"]))/1609.34)
    dist.df<-rbind(dist.df, c(pZip, hZip, distance))
  }
  dist.df2 <- dist.df
  colnames(dist.df2) <- c("Pt.Zip", "Hos.Zip", "Distance")
  smallest_dist<- dist.df2[which.min(dist.df2$Distance),]
  #print(smallest_dist)
  full.df<-rbind(full.df, smallest_dist)
}
print(full.df)
NearestTxCenter.df <- full.df
NearestTxCenter.df <- NearestTxCenter.df %>% rename(Distance_To_ClosestLiverTxCenter_miles = Distance, 
                                                    LiverTxCenterZip = Hos.Zip)

NAFLD.LTxCenter.df <- unique(merge(NAFLD.df, NearestTxCenter.df, by.x = "Zip_short", by.y = "Pt.Zip", all = T))
NAFLD.LTxCenter.df

write.csv(NAFLD.LTxCenter.df, "./NAFLD_distance_to_Hospitals.csv")


# Distances to Nearest Clinic 
NAFLD.df <- NAFLD_ADI.df

Mini.df <- NAFLD.df %>% select(c("PtID", "Zip_short", "ClosestGIClinic_Zip"))
Mini2.DF <- subset.data.frame(Mini.df, Mini.df$Zip_short != Mini.df$ClosestGIClinic_Zip)
Mini2.DF <- Mini.df[45:56,]

UniClinZips <- unique(NAFLD.df$ClosestGIClinic_Zip)

Mini.df$DistanceTo_GIclinic_meters<-apply(Mini.df, 1, function(x){
  startindex<-which(x[["Zip_short"]]==zip_code_db$zipcode)
  endindex<- which(x[["ClosestGIClinic_Zip"]]==zip_code_db$zipcode)
  distGeo(p1=c(zip_code_db[startindex, "lng"], 
               zip_code_db[startindex, "lat"]), 
          p2=c(zip_code_db[endindex, "lng"], 
               zip_code_db[endindex, "lat"]))
})
Mini.df<- Mini.df %>% mutate (DistanceTo_UWGIclinic_miles = (DistanceTo_GIclinic_meters /1609.34))

write.csv(Mini.df, "./NAFLD_distance_toUWGI_Clinics.csv")

UniClinZips <- unique(NAFLD.df$ClosestGIClinic_Zip)
UniClinZips2 <- c(53705, 53066, 53533, 53715, 54601,52001,53548, 53901, 53511, 54154,
                  53095, 54449, 61107, 33907, 53121, 54935,53098,80401,61032)
full.df<- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Pt.Zip", "GIclinic.Zip", "Distance"))))
for(pZip in Mini.df$Zip_short){
  dist.df<- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Pt.Zip", "GIclinic.Zip", "Distance"))))
  for(hZip in UniClinZips2 ){
    startindex<-which(pZip==zip_code_db$zipcode)
    #endindex<- (hZip==zip_code_db$zipcode)
    endindex<- (hZip==zip_code_db$zipcode)
    #endindex<- (c(54601, 52001, 53548, 53901, 61703, 53706, 54154, 53095, 54449, 61107)==zip_code_db$zipcode)
    distance = (distGeo(p1=c(zip_code_db[startindex, "lng"], 
                             zip_code_db[startindex, "lat"]), 
                        p2=c(zip_code_db[endindex, "lng"], 
                             zip_code_db[endindex, "lat"]))/1609.34)
    dist.df<-rbind(dist.df, c(pZip, hZip, distance))
  }
  dist.df2 <- dist.df
  colnames(dist.df2) <- c("Pt.Zip", "GIclinic.Zip", "Distance")
  smallest_dist<- dist.df2[which.min(dist.df2$Distance),]
  #print(smallest_dist)
  full.df<-rbind(full.df, smallest_dist)
}
print(full.df)
NearestClinic.df <- full.df

write.csv(NearestClinic.df, "./NAFLD_distance_to_GIClinics.csv")

