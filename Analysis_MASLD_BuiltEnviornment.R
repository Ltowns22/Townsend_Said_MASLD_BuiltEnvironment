# Townsend & Said 2024. 
# R code for evaluating the assocaitions between the built enviornment and NAFLD severity and progression 
# to ensure patient privacy subject zipcodes are removed from these files 

library(tidyr)
library(readr)
library(dplyr)
library(rstatix) 
library(ggplot2)
library(data.table)
library(stringr)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(Hmisc)
library(WGCNA)
library(reshape2)
library(devtools)
library(lme4)
library(mixOmics)
library(Matrix)
library(lme4)
library(lmerTest)

setwd("https://github.com/Ltowns22/Townsend_Said_MASLD_BuiltEnvironment")

# demographics and descriptive stats 
NAFLD.dem.df <- as.data.frame(read_csv("./NAFLD_ADI_all_Lab_FU_data.csv"))
colnames(NAFLD.dem.df)

NAFLD.dem.df$DistanceTo_UWhospital_miles[NAFLD.dem.df$DistanceTo_UWhospital_miles >500] <-NA # omiting the two patients who live more than 500 miles from UW 

NAFLD.d.df <-subset.data.frame(NAFLD.dem.df, select = c("PtID",
                                                        "ADI_2020state" , "ADI_2020national","NWI" , "DistanceTo_UWhospital_miles" ,"Distance_To_ClosestLiverTxCenter_miles", "Distance_To_GIhepClinic",
                                                        "Age_StudyEntry", "BMI_Entry" , "SystolicBloodPressure_E" ,"Fibroscan_kPa_E",  "FibroscanCAP_E", "ASCVD_E", "FIB-4_E" ,
                                                        "AST_E","ALT_E","PLt_E","Ferritin_E","Albumin_E",
                                                        "TotalChol_E", "TG_E","TotalChol_E" ,"HDL_E","LDL_E",
                                                        "HbA1c_E","Insulin_E", "FBS_E",  
                                                        "BMI_FU1", "SystolicBloodPressure_FU1","FibroscanKPA_FU1","FibroscanCAP_FU1", "ASCVD__FU1" , "FIB-4_FU1",
                                                        "AST_FU1" ,"ALT_FU1", "PLt_FU1" ,"Ferritin_FU1" ,"Albumin_FU1",
                                                        "TotalChol_FU1", "TG_FU1","HDL_FU1","LDL_FU1" ,                     
                                                        "HbA1c_FU1" , "Insulin_FU1","FBS_FU1",
                                                        "BMI_FU5", "SystolicBloodPressure_FU5","FibroscanKPa_FU5","FibroscanCAP_FU5", "ASCVD_FU5", "FIB-4_FU5",
                                                        "AST_FU5","ALT_FU5","PLt_FU5","Ferritin_FU5","Albumin_FU5",  
                                                        "TotalChol_FU5","TG_FU5", "HDL_FU5","LDL_FU5",
                                                        "HbA1c_FU5", "Insulin_FU5" , "FBS_FU5",
                                                        "Fib4_Bto1_Change", "Fib4_Bto5_Change","Fib4_1yr_Actual-Expect","Fib4_5yr_Actual-Expect",
                                                        "1_Change_ASCVD", "5_Change_ASCVD","ASCVD_1yr_Actual-Expect", "ASCVD_5yr_Actual-Expect", "1_Change_BMI", "5_Change_BMI" )) 
row.names(NAFLD.d.df) <- NAFLD.d.df$PtID
NAFLD.d.df<-NAFLD.d.df[,-1]

Descriptive.df<- NAFLD.d.df %>%  get_summary_stats(  type = "full") 
Descriptive.df$MeanSD <- paste(Descriptive.df$mean, " +/- ",Descriptive.df$sd )
Descriptive.df$range<- paste(Descriptive.df$min, " - ",Descriptive.df$max )

write.csv(Descriptive.df, "./NAFLD_ADI_DescriptiveStats.csv")

NAFLD.dist.df <-subset.data.frame(NAFLD.dem.df, select = c("PtID", "ADI_2020state" , "ADI_2020national","NWI" , "DistanceTo_UWhospital_miles" ,"Distance_To_ClosestLiverTxCenter_miles", "Distance_To_GIhepClinic"))
Descriptive.dist.df<- NAFLD.dist.df %>%  get_summary_stats(  type = "full") 
Descriptive.dist.df # getting IQR for distances

write.csv(Descriptive.df, "./NAFLD_ADI_DistancesDescriptiveStats.csv")


# Figure 1 Correlations with built environment 
NAFLD.BE.df <- NAFLD.dem.df

cor.test(NAFLD.BE.df$DistanceTo_UWhospital_miles, NAFLD.BE.df$`FIB-4_E`, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=DistanceTo_UWhospital_miles, y = `FIB-4_E`))+
  geom_point( size = 3.5)+ # 
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  annotate(geom="text", x=38, y=7, color="black",
           label="p-value = 0.0001 / spearman rho = 0.50")+
  ggtitle("Fib-4 vs. distance to UW")

cor.test(NAFLD.BE.df$Distance_To_ClosestLiverTxCenter_miles, NAFLD.BE.df$`FIB-4_E`, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=Distance_To_ClosestLiverTxCenter_miles, y = `FIB-4_E`))+
  geom_point(size = 3.5)+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  annotate(geom="text", x=38, y=7, color="black",
           label="p-value = 7.7e-05 / spearman rho = 0.51")+
  ggtitle("Fib-4 vs. distance to liver Transplant center")

cor.test(NAFLD.BE.df$Distance_To_GIhepClinic , NAFLD.BE.df$`FIB-4_E`, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=Distance_To_GIhepClinic, y = `FIB-4_E`))+
  geom_point( size = 3.5)+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  annotate(geom="text", x=8, y=7, color="black",
           label="p-value = 0.012 / spearman rho = 0.34")+
  ggtitle("Fib-4 vs. distance to GI clinic")

cor.test(NAFLD.BE.df$ADI_2020state , NAFLD.BE.df$`FIB-4_E`, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=ADI_2020state, y = `FIB-4_E`))+
  geom_point( size = 3, position=position_dodge2(width=0.1))+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  annotate(geom="text", x=2.5, y=7, color="black",
           label="p-value = 0.03 / spearman rho = 0.29")+
  ggtitle("Fib-4 vs. ADI State")

cor.test(NAFLD.BE.df$ADI_2020national , NAFLD.BE.df$`FIB-4_E`, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=ADI_2020national, y = `FIB-4_E`))+
  geom_point( size = 3.5)+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  xlim(0,100)+
  annotate(geom="text", x=30, y=7, color="black",
           label="p-value = 0.005 / spearman rho = 0.37")+
  ggtitle("Fib-4 vs. ADI National")

cor.test(NAFLD.BE.df$NWI , NAFLD.BE.df$`FIB-4_E`, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=NWI, y = `FIB-4_E`))+
  geom_point( size = 3.5)+
  #  geom_smooth(method='lm', color = "black")+
  theme_light()+
  xlim(0,20)+
  ggtitle("Fib-4 vs. NWI")


cor.test(NAFLD.BE.df$DistanceTo_UWhospital_miles , NAFLD.BE.df$ASCVD_E, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=DistanceTo_UWhospital_miles, y = `ASCVD_E`))+
  geom_point(aes(color = ADI_state_group), size = 3)+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  annotate(geom="text", x=38, y=30, color="black",
           label="p-value = 0.013  / spearman rho = 0.35")+
  ggtitle("ASCVD vs. distance to UW")

cor.test(NAFLD.BE.df$Distance_To_ClosestLiverTxCenter_miles , NAFLD.BE.df$ASCVD_E, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=Distance_To_ClosestLiverTxCenter_miles, y = `ASCVD_E`))+
  geom_point( size = 3.5)+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  annotate(geom="text", x=38, y=30, color="black",
           label="p-value = 0.013 / spearman rho = 0.35")+
  ggtitle("ASCVD vs. distance to Liver Tx center")

cor.test(NAFLD.BE.df$Distance_To_GIhepClinic , NAFLD.BE.df$ASCVD_E, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=Distance_To_GIhepClinic, y = `ASCVD_E`))+
  geom_point( size = 3.5)+
  #geom_smooth(method='lm', color = "black")+
  theme_light()+
  #annotate(geom="text", x=10, y=41, color="black",
  #label="NS p-value = 0.06 / spearman rho = 0.27")+
  ggtitle("ASCVD vs. distance to GI clinic")

cor.test(NAFLD.BE.df$ADI_2020state , NAFLD.BE.df$ASCVD_E, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=ADI_2020state, y = `ASCVD_E`))+
  geom_point( size = 3)+
  #geom_smooth(method='lm', color = "black")+
  theme_light()+
  #annotate(geom="text", x=38, y=7, color="black",
  #label="p-value = *** / spearman rho = ***")+
  ggtitle("ASCVD vs ADI State")

cor.test(NAFLD.BE.df$ADI_2020national , NAFLD.BE.df$ASCVD_E, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=ADI_2020national, y = `ASCVD_E`))+
  geom_point( size = 3.5)+
  geom_smooth(method='lm', color = "black")+
  theme_light()+
  xlim(0,100)+
  annotate(geom="text", x=30, y=36, color="black",
           label="p-value = 0.020 / spearman rho = 0.32")+
  ggtitle("ASCVD vs ADI National")

cor.test(NAFLD.BE.df$NWI, NAFLD.BE.df$ASCVD_E, method=c("spearman"))
ggplot(NAFLD.BE.df, aes(x=NWI, y = `ASCVD_E`))+
  geom_point( size = 3.5)+
  # geom_smooth(method='lm', color = "black")+
  theme_light()+
  xlim(0,20)+
  ggtitle("ASCVD vs NWI")


# mixed linear regressions 
me<- lmerTest::lmer(`FIB-4_E` ~  Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles +ADI_2020national + NWI
                    + (1|Age_StudyEntry)+ (1|CirrhosisYN)+ (1| EnrollSite) , data =NAFLD.BE.df )
summary(me)
anova(me)
me<- lmerTest::lmer(`ASCVD_E` ~  Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles +ADI_2020national+ NWI + 
                      + (1|Age_StudyEntry) +(1|Cardiac_E)+ (1| EnrollSite) , data =NAFLD.BE.df )
summary(me)
anova(me)

m2<-lmerTest::lmer(Death ~ ADI_2020national + Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles + NWI 
                   + (1|Age_StudyEntry)+(1| CirrhosisYN )+(1|Cardiac_E) , data =NAFLD.BE.df  )
summary(m2)
m2<-lmerTest::lmer(`Fib4_5yr_Actual-Expect` ~ ADI_2020national + Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles + NWI + (1|Age_StudyEntry)+(1| CirrhosisYN ) , data =NAFLD.BE.df  )
summary(m2)
m2<-lmerTest::lmer(`ASCVD_5yr_Actual-Expect` ~ ADI_2020national + Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles + NWI + (1|Age_StudyEntry)+(1| CirrhosisYN ) , data =NAFLD.BE.df  )
summary(m2)

NAFLD.BE.df<-NAFLD.BE.df %>% mutate(Fib4_5yr_progression = case_when(`Fib4_5yr_Actual-Expect` <= -0.25 ~ 0,
                                                                     `Fib4_5yr_Actual-Expect` > -0.25 &`Fib4_5yr_Actual-Expect` < 0.25 ~ 0.5,
                                                                     `Fib4_5yr_Actual-Expect` >=0.25 ~ 1))
#Outcome_Death_Liver_Cardio == "Liver Death" ~ "Progress"))
NAFLD.BE.df<-NAFLD.BE.df %>% mutate(ASCVD_5yr_progression = case_when(`ASCVD_5yr_Actual-Expect`  <= -2.5 ~ 0,
                                                                      `ASCVD_5yr_Actual-Expect`>-2.5 & `ASCVD_5yr_Actual-Expect`<2.5 ~ 0.5,
                                                                      `ASCVD_5yr_Actual-Expect` >=2.5 ~ 1))
# Outcome_Death_Liver_Cardio == "Cardiac Death" ~ "Progress"))

m2<-lmerTest::lmer(Fib4_5yr_progression ~ ADI_2020national + Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles + NWI + (1|Age_StudyEntry) , data =NAFLD.BE.df  )
summary(m2)
m2<-lmerTest::lmer(ASCVD_5yr_progression ~ ADI_2020national + Distance_To_GIhepClinic + Distance_To_ClosestLiverTxCenter_miles + NWI + (1|Age_StudyEntry)+(1| CirrhosisYN ) , data =NAFLD.BE.df  )
summary(m2)




# Figure 2 Over time charts 
NAFLD.T.df <- as.data.frame(read_csv("./NAFLD_OVERTIME.csv"))
NAFLD.T.df$Outcome[NAFLD.T.df$Outcome == "Death / lost to followup 1 yr"] <- "Death"
NAFLD.T.df<-NAFLD.T.df %>% mutate(AWI_state_group = case_when(ADI_2020state <= 3 ~ "State AWI 1-3",
                                                              ADI_2020state > 3 & ADI_2020state< 7 ~ "State AWI 3-6",
                                                              ADI_2020state >=7 ~ "State AWI 7-10"))


ggplot(NAFLD.T.df, aes(x = Year, y= BMI, shape = Outcome))+
  geom_point(aes(shape=Outcome), fill = "black", size = 4, alpha = 0.5)+
  geom_line(aes(group = PtID))+
  geom_hline(yintercept = 18.5, linetype = "dashed")+
  geom_hline(yintercept = 25, linetype = "dashed")+
  geom_hline(yintercept = 30, linetype = "dashed")+
  geom_hline(yintercept = 35, linetype = "dashed")+
  geom_hline(yintercept = 40, linetype = "dashed")+
  theme_light()+
  scale_shape_manual(values = c(22,21, 24,25))+
  ylim(18,60)+
  ggtitle("BMI over time")

ggplot(NAFLD.T.df, aes(x = Year, y= HbA1c, shape = Outcome))+
  geom_point(aes(shape=Outcome), fill = "black", size = 4, alpha = 0.5)+
  geom_line(aes(group = PtID))+
  geom_hline(yintercept = 5.7, linetype = "dashed")+
  geom_hline(yintercept = 6.5, linetype = "dashed")+
  theme_light()+
  scale_shape_manual(values = c(22,21, 24,25))+
  ylim(4,11)+
  ggtitle("HbA1C over time")

ggplot(NAFLD.T.df, aes(x = Year, y=`FIB-4`))+
  geom_point(aes(shape=Outcome), fill = "black", size = 4, alpha = 0.5)+
  geom_line(aes(group = PtID))+
  geom_hline(yintercept = 1.33, linetype = "dashed")+
  geom_hline(yintercept = 2.67, linetype = "dashed")+
  theme_light()+
  scale_shape_manual(values = c(22,21, 24,25))+
  ylim(0,12)+
  ggtitle("FIB-4 over time")

ggplot(NAFLD.T.df, aes(x = Year, y=ASCVD))+
  geom_point(aes(shape=Outcome), fill = "black", size = 4, alpha = 0.5)+
  geom_line(aes(group = PtID))+
  geom_hline(yintercept = 5, linetype = "dashed")+
  geom_hline(yintercept = 7.5, linetype = "dashed")+
  geom_hline(yintercept = 20, linetype = "dashed")+
  theme_light()+
  scale_shape_manual(values = c(22,21, 24,25))+
  #ylim(0,12)+
  ggtitle("ASCVD over time")


# Figure 3: Partial Least Squares Discriminant Analysis 
### PLS-DA - determining which features help identify patients who die vs are alive at 5 yrs, progress vs. improve in FIB at 5 year, and progress vs. improve in ASCVD at 5 years.  
Outcome.df<- as.data.frame(read.csv("./Fib4_ASCVD_Outcomes.csv"))
rownames(Outcome.df) <- Outcome.df$PtID
Outcome.df<-Outcome.df %>% mutate(Fib4_5yr_progression = case_when(Fib4_5yr_Actual.Expect <= -0.25 ~ "Improve",
                                                                   Fib4_5yr_Actual.Expect > -0.25 & Fib4_5yr_Actual.Expect < 0.25 ~ "As expected with age",
                                                                   Fib4_5yr_Actual.Expect >=0.25 ~ "Progress",
                                                                   Outcome_Death_Liver_Cardio == "Liver Death" ~ "Progress"))
Outcome.df<-Outcome.df %>% mutate(ASCVD_5yr_progression = case_when(ASCVD_5yr_Actual.Expect  <= -2.5 ~ "Improve",
                                                                    ASCVD_5yr_Actual.Expect>-2.5 & ASCVD_5yr_Actual.Expect <2.5 ~ "as expected with age",
                                                                    ASCVD_5yr_Actual.Expect >=2.5 ~ "Progress",
                                                                    Outcome_Death_Liver_Cardio == "Cardiac Death" ~ "Progress"))



Clin.df<- as.data.frame(read.csv("./Clinics_ADI_NWI.csv"))
rownames(Clin.df) <- Clin.df$PtID
Clin.df <- Clin.df[,-1]
Clin.df<- subset(Clin.df, select = c("ADI_2020national" , "NWI",
                                     "Distance_To_ClosestLiverTxCenter_miles", "Distance_To_GIhepClinic"  ))
Clin_ADI_List <- colnames(Clin.df)

Vitals.df<- as.data.frame(read.csv("./VItals_Dem_Entry.csv"))
age.df <-subset(Vitals.df, select = c("PtID", "Age_StudyEntry"))
rownames(Vitals.df) <- Vitals.df$PtID
Vitals.df <- Vitals.df[,-1]
Vitals_List <- colnames(Vitals.df)

CoMorb.df<- as.data.frame(read.csv("./LiverDisease_CoMorbid.csv"))
rownames(CoMorb.df) <- CoMorb.df$PtID
CoMorb.df <- CoMorb.df[,-1]
CoMorb_List <- colnames(CoMorb.df)

Labs.df<- as.data.frame(read.csv("./EntryLabs.csv"))
Labs.df <- merge(Labs.df, age.df, by = "PtID")
rownames(Labs.df) <- Labs.df$PtID
Labs.df <- Labs.df[,-1]
Labs_List <- colnames(Labs.df)


# ASCVD 5yr 
Outcome.FIN.df <- subset.data.frame(Outcome.df, ASCVD_5yr_progression %in% c( "Improve","as expected with age", "Progress") ) 

Clin.FIN.df <- subset.data.frame(Clin.df, row.names(Clin.df) %in% row.names(Outcome.FIN.df))
CoMorb.FIN.df<- subset.data.frame(CoMorb.df, row.names(CoMorb.df) %in% row.names(Outcome.FIN.df))
Labs.FIN.df <- subset.data.frame(Labs.df, row.names(Labs.df) %in% row.names(Outcome.FIN.df))

X <- list(CLinic_ADI = Clin.FIN.df,
          #Vitals= Vitals.FIN.df,
          MedHX = CoMorb.FIN.df,
          Labs = Labs.FIN.df)

Y <- as.factor(Outcome.FIN.df$ASCVD_5yr_progression) ## 


Outcome.diablo <- block.plsda(X, Y, ncomp = 3)
plotDiablo(Outcome.diablo, ncomp = 1)

coordinates <- plotVar(Outcome.diablo , plot = F) # to save the variables
plotVar(Outcome.diablo , #var.names = c(FALSE,F,F),
        cutoff = 0.6,
        legend=TRUE, #pch=c(16,16,16),
        title = "Diablo variable plot - cutoff 0.5")

plotIndiv(Outcome.diablo ,comp  = c(1,2),
          ind.names = F, 
          legend=TRUE, cex=6,
          #centroid = T,
          #col =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'NAFLD with DIABLO')

plotArrow(Outcome.diablo, comp  = c(1,2),
          ind.names = F,
          group = Outcome.FIN.df$ASCVD_5yr_progression, ## 
          legend=TRUE, #cex=c(5,5,5,5),
          pch.size =5.5, 
          arrow.size = 0.5,
          #col.per.group =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-4,4), ylim = c(-5,5),
          title = 'NAFLD with DIABLO') +
  scale_shape_manual(values = c(centroid = 19,CLinic_ADI = 45, MedHX = 42, Labs = 43))


plotLoadings(Outcome.diablo, comp = 1, contrib = 'max', method = 'median')
plotLoadings(Outcome.diablo, comp = 2, contrib = 'max', method = 'median')

cord1_Load_ADI<-plotLoadings(Outcome.diablo, comp = 1, block = "CLinic_ADI", contrib = 'max', method = 'median')
cord1_Load_MedHX<-plotLoadings(Outcome.diablo, comp = 1, block = "MedHX", contrib = 'max', method = 'median')
cord1_Load_Labs<-plotLoadings(Outcome.diablo, comp = 1, block = "Labs", contrib = 'max', method = 'median')
Cord_1_Loadings_dataframe <- rbind(cord1_Load_ADI$CLinic_ADI, cord1_Load_MedHX$MedHX,cord1_Load_Labs$Labs)
Cord_1_Loadings_dataframe <- tibble::rownames_to_column(Cord_1_Loadings_dataframe, "Variable")

cord2_Load_ADI<-plotLoadings(Outcome.diablo, comp = 2, block = "CLinic_ADI", contrib = 'max', method = 'median')
cord2_Load_MedHX<-plotLoadings(Outcome.diablo, comp = 2, block = "MedHX", contrib = 'max', method = 'median')
cord2_Load_Labs<-plotLoadings(Outcome.diablo, comp = 2, block = "Labs", contrib = 'max', method = 'median')
Cord_2_Loadings_dataframe <- rbind(cord2_Load_ADI$CLinic_ADI, cord2_Load_MedHX$MedHX,cord2_Load_Labs$Labs)
Cord_2_Loadings_dataframe <- tibble::rownames_to_column(Cord_2_Loadings_dataframe, "Variable")

Cord_both_loadings.df <- merge(Cord_1_Loadings_dataframe,Cord_2_Loadings_dataframe, by.x = "Variable",by.y = "Variable" , all = T )
Cord_both_loadings.df$TotalImprotance <- sqrt((Cord_both_loadings.df$importance.x^2)+ (Cord_both_loadings.df$importance.y^2))
Cord_both_loadings.df$GroupContrib.combined<- paste(Cord_both_loadings.df$GroupContrib.x , Cord_both_loadings.df$GroupContrib.y)
Cord_both_loadings.df$GroupContrib.combined
Cord_both_loadings.df <- subset.data.frame(Cord_both_loadings.df, GroupContrib.combined %in% c("Improve Improve" ,"as expected with age as expected with age","Progress Progress") )
Cord_both_loadings.df <- Cord_both_loadings.df %>% mutate(VariableGroup= case_when(Variable %in% Clin_ADI_List ~ 'Built Enviornment', 
                                                                                   Variable %in% CoMorb_List ~ 'CoMorbidities',
                                                                                   Variable %in% Labs_List ~ 'Labs at Entry'))
Cord_both_loadings.df<- subset.data.frame(Cord_both_loadings.df,TotalImprotance > 0.4)
ggplot(Cord_both_loadings.df, aes(x = TotalImprotance, y = reorder(Variable, abs(TotalImprotance)), color = VariableGroup, fill = VariableGroup))+
  # geom_point(size = 6)+
  geom_col(width = 0.2)+
  xlim(0, 1)+
  # scale_color_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  #scale_fill_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  facet_grid(rows=vars(factor(GroupContrib.combined, levels = c("Improve Improve","as expected with age as expected with age" ,"Progress Progress"))), scales = "free_y",space = "free_y")+
  theme_bw()+
  ggtitle("predictive factors for ASCVD improve/ progress at 5 yr ")


# Fib 4 at 5 yr 
Outcome.FIN.df <- subset.data.frame(Outcome.df, Fib4_5yr_progression %in% c( "Improve", "Progress","As expected with age") ) # ASCVD 1yr 

Clin.FIN.df <- subset.data.frame(Clin.df, row.names(Clin.df) %in% row.names(Outcome.FIN.df))
CoMorb.FIN.df<- subset.data.frame(CoMorb.df, row.names(CoMorb.df) %in% row.names(Outcome.FIN.df))
Labs.FIN.df <- subset.data.frame(Labs.df, row.names(Labs.df) %in% row.names(Outcome.FIN.df))

X <- list(CLinic_ADI = Clin.FIN.df,
          # Vitals= Vitals.FIN.df,
          MedHX = CoMorb.FIN.df,
          Labs = Labs.FIN.df)

Y <- as.factor(Outcome.FIN.df$Fib4_5yr_progression) ## 


Outcome.diablo <- block.plsda(X, Y, ncomp = 3)
plotDiablo(Outcome.diablo, ncomp = 1)

coordinates <- plotVar(Outcome.diablo , plot = F) # to save the variables
plotVar(Outcome.diablo , #var.names = c(FALSE,F,F),
        cutoff = 0.6,
        legend=TRUE, #pch=c(16,16,16),
        title = "Diablo variable plot - cutoff 0.5")

plotIndiv(Outcome.diablo ,comp  = c(1,2),
          ind.names = F, 
          legend=TRUE, cex=6,
          #centroid = T,
          #col =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'NAFLD with DIABLO')

plotArrow(Outcome.diablo, comp  = c(1,2),
          ind.names = F,
          group = Outcome.FIN.df$Fib4_5yr_progression, ## 
          legend=TRUE, #cex=c(5,5,5,5),
          pch.size =5.5, 
          arrow.size = 0.5,
          #col.per.group =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-4,4), ylim = c(-5,5),
          title = 'NAFLD with DIABLO') +
  scale_shape_manual(values = c(centroid = 19,CLinic_ADI = 45, MedHX = 42, Labs = 43))


plotLoadings(Outcome.diablo, comp = 1, contrib = 'max', method = 'median')
plotLoadings(Outcome.diablo, comp = 2, contrib = 'max', method = 'median')

cord1_Load_ADI<-plotLoadings(Outcome.diablo, comp = 1, block = "CLinic_ADI", contrib = 'max', method = 'median')
cord1_Load_MedHX<-plotLoadings(Outcome.diablo, comp = 1, block = "MedHX", contrib = 'max', method = 'median')
cord1_Load_Labs<-plotLoadings(Outcome.diablo, comp = 1, block = "Labs", contrib = 'max', method = 'median')
Cord_1_Loadings_dataframe <- rbind(cord1_Load_ADI$CLinic_ADI, cord1_Load_MedHX$MedHX,cord1_Load_Labs$Labs)
Cord_1_Loadings_dataframe <- tibble::rownames_to_column(Cord_1_Loadings_dataframe, "Variable")

cord2_Load_ADI<-plotLoadings(Outcome.diablo, comp = 2, block = "CLinic_ADI", contrib = 'max', method = 'median')
cord2_Load_MedHX<-plotLoadings(Outcome.diablo, comp = 2, block = "MedHX", contrib = 'max', method = 'median')
cord2_Load_Labs<-plotLoadings(Outcome.diablo, comp = 2, block = "Labs", contrib = 'max', method = 'median')
Cord_2_Loadings_dataframe <- rbind(cord2_Load_ADI$CLinic_ADI,cord2_Load_MedHX$MedHX,cord2_Load_Labs$Labs)
Cord_2_Loadings_dataframe <- tibble::rownames_to_column(Cord_2_Loadings_dataframe, "Variable")

Cord_both_loadings.df <- merge(Cord_1_Loadings_dataframe,Cord_2_Loadings_dataframe, by.x = "Variable",by.y = "Variable" , all = T )
Cord_both_loadings.df$TotalImprotance <- sqrt((Cord_both_loadings.df$importance.x^2)+ (Cord_both_loadings.df$importance.y^2))
Cord_both_loadings.df$GroupContrib.combined<- paste(Cord_both_loadings.df$GroupContrib.x , Cord_both_loadings.df$GroupContrib.y)
Cord_both_loadings.df$GroupContrib.combined
Cord_both_loadings.df <- subset.data.frame(Cord_both_loadings.df, GroupContrib.combined %in% c("Improve Improve" ,"As expected with age As expected with age","Progress Progress") )
Cord_both_loadings.df <- Cord_both_loadings.df %>% mutate(VariableGroup= case_when(Variable %in% Clin_ADI_List ~ 'Built Enviornment', 
                                                                                   Variable %in% CoMorb_List ~ 'CoMorbidities',
                                                                                   Variable %in% Labs_List ~ 'Labs at Entry'))

Cord_both_loadings.df<- subset.data.frame(Cord_both_loadings.df,TotalImprotance > 0.5)
ggplot(Cord_both_loadings.df, aes(x = TotalImprotance, y = reorder(Variable, abs(TotalImprotance)), color = VariableGroup, fill = VariableGroup))+
  #geom_point(aes(size = abs(TotalImprotance)))+
  geom_col(width = 0.2)+
  xlim(0, 1)+
  # scale_color_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  #scale_fill_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  facet_grid(rows=vars(factor(GroupContrib.combined, levels = c("Improve Improve","As expected with age As expected with age" ,"Progress Progress"))), scales = "free_y",space = "free_y")+
  theme_bw()+
  ggtitle("predictive factors for FIB4 improve/ progress at 5 yr ")


# All cause mortality 
Outcome.FIN.df <- subset.data.frame(Outcome.df, Outcome_Death %in% c( "Death", "Alive")) 

Clin.FIN.df <- subset.data.frame(Clin.df, row.names(Clin.df) %in% row.names(Outcome.FIN.df))
CoMorb.FIN.df<- subset.data.frame(CoMorb.df, row.names(CoMorb.df) %in% row.names(Outcome.FIN.df))
Labs.FIN.df <- subset.data.frame(Labs.df, row.names(Labs.df) %in% row.names(Outcome.FIN.df))

X <- list(CLinic_ADI = Clin.FIN.df,
          MedHX = CoMorb.FIN.df,
          Labs = Labs.FIN.df)

Y <- as.factor(Outcome.FIN.df$Outcome_Death) ## 


Outcome.diablo <- block.plsda(X, Y, ncomp = 3)
plotDiablo(Outcome.diablo, ncomp = 1)

coordinates <- plotVar(Outcome.diablo , plot = F) # to save the variables
plotVar(Outcome.diablo , #var.names = c(FALSE,F,F),
        cutoff = 0.6,
        legend=TRUE, #pch=c(16,16,16),
        title = "Diablo variable plot - cutoff 0.5")

plotIndiv(Outcome.diablo ,comp  = c(1,2),
          ind.names = F, 
          legend=TRUE, cex=6,
          centroid = F,
          #col =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'NAFLD with DIABLO')
plotArrow(Outcome.diablo, comp  = c(1,2),
          ind.names = F,
          group = Outcome.FIN.df$Outcome_Death, ## 
          legend=TRUE, #cex=c(5,5,5,5),
          pch.size =5.5, 
          arrow.size = 0.5,
          #col.per.group =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-4,4), ylim = c(-5,5),
          title = 'NAFLD with DIABLO') +
  scale_shape_manual(values = c(centroid = 19,CLinic_ADI = 45, MedHX = 42, Labs = 43))+
  
  plotLoadings(Outcome.diablo, comp = 1, contrib = 'max', method = 'median')
plotLoadings(Outcome.diablo, comp = 2, contrib = 'max', method = 'median')

cord1_Load_ADI<-plotLoadings(Outcome.diablo, comp = 1, block = "CLinic_ADI", contrib = 'max', method = 'median')
cord1_Load_MedHX<-plotLoadings(Outcome.diablo, comp = 1, block = "MedHX", contrib = 'max', method = 'median')
cord1_Load_Labs<-plotLoadings(Outcome.diablo, comp = 1, block = "Labs", contrib = 'max', method = 'median')
Cord_1_Loadings_dataframe <- rbind(cord1_Load_ADI$CLinic_ADI,cord1_Load_MedHX$MedHX,cord1_Load_Labs$Labs)
Cord_1_Loadings_dataframe <- tibble::rownames_to_column(Cord_1_Loadings_dataframe, "Variable")

cord2_Load_ADI<-plotLoadings(Outcome.diablo, comp = 2, block = "CLinic_ADI", contrib = 'max', method = 'median')
cord2_Load_MedHX<-plotLoadings(Outcome.diablo, comp = 2, block = "MedHX", contrib = 'max', method = 'median')
cord2_Load_Labs<-plotLoadings(Outcome.diablo, comp = 2, block = "Labs", contrib = 'max', method = 'median')
Cord_2_Loadings_dataframe <- rbind(cord2_Load_ADI$CLinic_ADI, cord2_Load_MedHX$MedHX,cord2_Load_Labs$Labs)
Cord_2_Loadings_dataframe <- tibble::rownames_to_column(Cord_2_Loadings_dataframe, "Variable")

Cord_both_loadings.df <- merge(Cord_1_Loadings_dataframe,Cord_2_Loadings_dataframe, by.x = "Variable",by.y = "Variable" , all = T )
Cord_both_loadings.df$TotalImprotance <- sqrt((Cord_both_loadings.df$importance.x^2)+ (Cord_both_loadings.df$importance.y^2))
Cord_both_loadings.df$GroupContrib.combined<- paste(Cord_both_loadings.df$GroupContrib.x , Cord_both_loadings.df$GroupContrib.y)
Cord_both_loadings.df$GroupContrib.combined
Cord_both_loadings.df <- subset.data.frame(Cord_both_loadings.df, GroupContrib.combined %in% c("Death Death" ,"Alive Alive") )
Cord_both_loadings.df <- Cord_both_loadings.df %>% mutate(VariableGroup= case_when(Variable %in% Clin_ADI_List ~ 'Built Enviornment', 
                                                                                   Variable %in% CoMorb_List ~ 'CoMorbidities',
                                                                                   Variable %in% Labs_List ~ 'Labs at Entry'))

Cord_both_loadings.df<- subset.data.frame(Cord_both_loadings.df,TotalImprotance > 0.38)
ggplot(Cord_both_loadings.df, aes(x = TotalImprotance, y = reorder(Variable, abs(TotalImprotance)), color = VariableGroup, fill = VariableGroup))+
  #geom_point(aes(size = abs(TotalImprotance)))+
  geom_col(width = 0.2)+
  xlim(0, 1)+
  # scale_color_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  #scale_fill_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  facet_grid(rows=vars(factor(GroupContrib.combined, levels = c("Alive Alive" ,"Death Death"))), scales = "free_y",space = "free_y")+
  #facet_wrap(~factor(GroupContrib.combined, levels = c("Alive Alive" ,"Death Death")), scales = "free_y")+
  theme_bw()+
  ggtitle("predictive factors for Death within 5 yrs")

