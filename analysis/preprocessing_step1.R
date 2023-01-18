#   This file is the main program for KOMP/data/2021 preprocessing procedure (17 domains)
#
#   This R code converts the meta data and long-format KOMP2 data
#   files to the format used in the KOMP analysis
#   pipelines.
#
#   Main functions:
#
#   A: long-format to wide-format
#
#   B: outlier removal
#
#  Modified by ___Hao He________  Date ___10/03/2021_______


##-------Define enviroment and source all util codes----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#install packages
install.packages(c("grid", "ggplot2", "readxl", "stringr", "plyr", "reshape2", "lubridate",
                   "gdata", "glmnet", "nlme", "scales", "plotly", "gridExtra", "data.table",
                   "kableExtra", "pheatmap", "viridis", "cowplot", "RColorBrewer", "tidyverse",
                   "ComplexHeatmap", "circlize"
                   ))

#define the path where the LIMS data deposit
## Define enviroment and source all util codes
batch = "2022-05-13" #to change for different version of data
projectDir = "/projects/compsci/vmp/USERS/heh/KOMP2022/" #project dir
dataDir =  paste0(projectDir, "data/", batch, "/")
workDir =  paste0(projectDir, "output/", batch, "_out/")
if (!dir.exists(workDir))
  { dir.create(workDir) } #create folder

#plot functions
plotAllFile = paste0(projectDir, "code/RCodePlotEverything_06072018.R")
plotAllFile2 = paste0(projectDir, "code/RCodePlotEverything2.R") # for ECT that has no time of test

# READ Functions file
source(file = paste0(projectDir, "code/RCode_Functions_v3.R"), local  = TRUE)
setwd(workDir)

#   Body Weight --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   RUN THIS!!!!
# 	BODY WEIGHT  ## Run this
# 	How to call the function
#		Read in KOMP BW data
# 	bw.komp <- read.csv("BW_reformatted_12022015_noNA.csv",header=T)
#		loadPheno
#		Do all the functions that go after loadPheno
#		At the very end call the addBW function
#		It takes three parameters:
#		1. data = OFA (after it has been loaded and cleaned)
#		2. bw.data = bw.komp (KOMP BW data)
#		3. prefix =  "OFA" (This is needed to prepend "OFA" to the closestbw and predbw columns)
#		It returns a data object with two new columms: OFA_closestbw and OFA_predbw
#		Any negative pred bw  are set to NA
#		Any mutant group with n < 3 has a pred bw set to NA

#PreProcess of BW data
options(stringsAsFactors = FALSE)
#BW <- read.csv("BW.csv",header=T)
BW <- read.delim(paste0(dataDir ,"BodyWeight.csv"),header=T, sep = "\t")
#subset all rows with body weight
BW.weight <- BW[grep(BW$Name,pattern = "Body Weight"),]
#subset all rows with weeks
BW.weight.weeks <- BW[grep(BW$Name,pattern = "Comments  week"),]
#strip time from the date time field
tmp.date.time <- BW.weight$Date.of.Birth
#dob <-as.Date(sapply(strsplit(tmp.date.time,' '),'[',1),'%m/%d/%Y') # vivek P.'s doesnt work for me
dob <- MakeSameDateFormat(tmp.date.time) # use Don's function for fixing date
BW.weight$dob <- dob
#add two new columns one being weighed_date and weighed_date_MouseID
#weighed_date_MouseID is mostly for sanity checking to make sure that the weighed date and dob are for the same mouse
BW.weight$weighed_date <- MakeSameDateFormat(BW.weight.weeks$DateComplete)
#bw.weight$weighed_date_MouseID <- BW.weight.weeks$Mouse.ID
BW.weight$Value <- as.numeric(BW.weight$Value) # BW is a character for some reason
BW.weight$BW <- BW.weight$Value
#bw.weight$logBW <- log(BW.weight$Value,10)
#calculate the age in days
#BW.weight$age <- as.numeric(as.Date(BW.weight$weighed_date) - as.Date(BW.weight$dob))
BW.weight$age <- as.numeric(BW.weight$weighed_date - BW.weight$dob)
#remove NA rows and write results - remove all with BW that is NA or if logBW is NA
na.ids <- which(is.na(BW.weight$age) == TRUE | is.na(BW.weight$BW) == TRUE )
BW.weight.noNAs <- BW.weight[-na.ids,]
BW.weight.noNAs <- BW.weight.noNAs[BW.weight.noNAs$BW != 0,]
BW.weight.noNAs <- BW.weight.noNAs[!(BW.weight.noNAs$BW > 80),]
write.csv(file="BWb.csv",BW.weight.noNAs,quote=TRUE,row.names = FALSE)
bw.komp <- BW.weight.noNAs # bw.komp is used in the function below
rm(BW, BW.weight, BW.weight.weeks, BW.weight.noNAs, na.ids, dob, tmp.date.time)

# HB --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#	Read Data
HB <-	loadPheno(paste0(dataDir,"HB.csv"))
HB <- fixTesterNames(HB, "HB")
HB <- convertVars(data = HB, dict = DictData, "HB")  #change classes of factors from character to the appropriate
# the Test.Date field is not used for HB
#	The Start.Time has test date and time - So, here we recreate the test date data from the test time data
#	2012-06-13 09:00:00
# replace incompletely reported event with NA
# some start.time is not properly reported.  The start.time have date and time, but a few only have time
# Mike M. is fixing these, but meanwhile, this code will remove these incomplete data
HB$Start.Time.New <- strptime(HB$Start.Time, "%Y-%m-%d %H:%M")  # Convert Start.Time to date and time
HB$Date.of.test.New <- strptime(HB$Start.Time.New, "%Y-%m-%d")
HB$Start.time.New <- strftime(HB$Start.Time.New, "%H:%M:%S")
HB$Start.time.New <- strptime(HB$Start.time.New, "%H:%M:%S")
HB <- ExtractDateParameters(HB, 'HB')
HB <- ExtractTimeParameters(HB, 'HB')

# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(HB$Date.of.test.New))))
# Remove all data before 2013
HB <- HB[as.Date(HB$Date.of.test.New) > ymd("2012-12-31"),]
# Equipment Model is all the same - Chamber ID is not recorded.
HB$Equipment.Model[HB$Equipment.Model == "Accuscan Versamax"] <- "Versamax"
HB$Equipment.Model[HB$Equipment.Model == "VMX 1.4b"] <- "Versamax"
# Calculate HB parameters %corner pokes and repetitive pokes

HB$Corner <-	HB$Parameters.Total.Hole.Pokes.Hole.1 +
  HB$Parameters.Total.Hole.Pokes.Hole.4 +
  HB$Parameters.Total.Hole.Pokes.Hole.13 +
  HB$Parameters.Total.Hole.Pokes.Hole.16

HB$Center <- 	HB$Parameters.Total.Hole.Pokes.Hole.6 +
  HB$Parameters.Total.Hole.Pokes.Hole.7 +
  HB$Parameters.Total.Hole.Pokes.Hole.10 +
  HB$Parameters.Total.Hole.Pokes.Hole.11

# calculate % corner hole and %center hole
HB$Percent.Corner <- HB$Corner/HB$Parameters.Total.Hole.Pokes*100
HB$Percent.Center <- HB$Center/HB$Parameters.Total.Hole.Pokes*100
#HB$Percent.Center2 <- HB$Center/HB$Parameters.Total.Hole.Pokes
# calculate what is the max percent one hole is being explored.  I mean if an animal like to just poke one hole, can we detect it.
# I am finding the whole that is poked the most and then dividing by total number of nose pokes.

HB1 <- 	 HB[c("Parameters.Total.Hole.Pokes.Hole.1",
              "Parameters.Total.Hole.Pokes.Hole.2",
              "Parameters.Total.Hole.Pokes.Hole.3",
              "Parameters.Total.Hole.Pokes.Hole.4" ,
              "Parameters.Total.Hole.Pokes.Hole.5",
              "Parameters.Total.Hole.Pokes.Hole.6" ,
              "Parameters.Total.Hole.Pokes.Hole.7" ,
              "Parameters.Total.Hole.Pokes.Hole.8" ,
              "Parameters.Total.Hole.Pokes.Hole.9" ,
              "Parameters.Total.Hole.Pokes.Hole.10",
              "Parameters.Total.Hole.Pokes.Hole.11" ,
              "Parameters.Total.Hole.Pokes.Hole.12",
              "Parameters.Total.Hole.Pokes.Hole.13" ,
              "Parameters.Total.Hole.Pokes.Hole.14",
              "Parameters.Total.Hole.Pokes.Hole.15" ,
              "Parameters.Total.Hole.Pokes.Hole.16")]
HB1$MaxVal <- apply(HB1, 1, max)
HB1$SumVal <- apply(HB1, 1, sum)
HB$PercentMaxHole <- HB1$MaxVal/HB1$SumVal * 100
rm(HB1)

# REMOVE OUTLIERS -- NOT USED FOR HB
PhenoList <- c(
  'Percent.Corner',
  'Percent.Center',
  'PercentMaxHole',
  'Corner',
  'Parameters.Total.Hole.Pokes',
  'Center'
)


# EXPORT THE DATA before outlier removal
write.csv(HB, file = "HBOut0.csv",row.names = FALSE)
# Remove outliers from data
# Outliers NOT removed for this data since it is count data and not normal.
HB <- remove_outliers_main(komp.data = HB ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
HB <- addBW(HB, bw.komp, "HB", avg_aot = 63)
# also add name of file to col names for later merging
# the [-c(1:3) will not paste on OFA or LD to the first three columns.  I will use these for merging.
colnames(HB)[-c(1:5)] = paste("HB_", colnames(HB)[-c(1:5)], sep = "")

save(HB, file="HBout.RData")
# EXPORT THE DATA
write.csv(HB, file = "HBOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoHB"
out1Dir <- "ControlDiagnotsticPlots"
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  'HB_Percent.Corner',
  'HB_Percent.Center',
  'HB_PercentMaxHole',
  'HB_Corner',
  'HB_Parameters.Total.Hole.Pokes',
  'HB_Center'
)

DataSet <- HB

TestDates <- "HB_Date.of.test.New"
TestTime <- "HB_Start.time.New"
DataName <- "HB"
Year <- "HB_Year"
Month <- "HB_Month"
DayOfWeek <- "HB_DayOfWeek"
Hour <- "HB_Hour"
timeCat <- "HB_timeCat"
DayOfYear <- "HB_DayOfYear"
DayOfMonth <- "HB_DayOfMonth"
StrainGeno <- "HB_StrainGeno"
predBW <- "HB_PredictedBW_outs"
#EScell <- "HB_EScell"

ExperimentalFactors <- c(
  "Sex",
  "HB_Experimenter",
  "HB_Room.origin",
  "HB_Equipment.Model",
  "HB_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
  #EScell
)

control <- HB[HB$HB_StrainGeno == "C57BL/6NJ ",]  ## control is now by itself, has an space at end
#Run PlotEverything scrips
#source(plotAllFile)

# OFA --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#	Read Data
OFA <-	loadPheno(paste0(dataDir,"OFA.csv"))
OFA <- fixTesterNames(OFA, "OFA")
OFA <- convertVars(data = OFA, dict = DictData, "OFA")  # change classes of factors from character to the appropriate
OFA$Date.of.test.New <- ymd_hms(OFA$Date.of.test)  # 4 failed to parse
# Start.time has date and time OR just date, extract just the Time
data1 <- strsplit(OFA$Start.time," ")
data1 <- sapply(data1, "[",2)
data1 <- as.data.frame(data1)
colnames(data1) <- "Start.time.New"

OFA <- cbind(OFA, data1)
OFA$Start.time.New <- strptime(OFA$Start.time.New, "%H:%M:%S")

# remove rearing data from box Arena.ID == 4 & Outputs.Number.of.Rears.Total == 0
# there was equipment malfunction in 2012 and 2013 and led to no rearing in box 4 only
# I am removing all rearing data with 0 or  no rearing.
OFA$Outputs.Number.of.Rears.Total[OFA$Outputs.Number.of.Rears.Total == "0"] <- NA

# Add test date and time parameters
OFA <- ExtractDateParameters(OFA, 'OFA')
OFA <- ExtractTimeParameters(OFA, 'OFA')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(OFA$Date.of.test.New))))
# Remove all data before 2013
OFA <- OFA[!is.na(OFA$Date.of.test.New),]
OFA <- OFA[OFA$Date.of.test.New > ymd("2012-12-31"),]

# Testers missing and marked as Unknown - fixed and changed to the person testing at that time
OFA$Experimenter.ID[OFA$Date.of.test.New == ymd(20160822)] <- "Pamelia Fraungruber"
OFA$Experimenter.ID[OFA$Date.of.test.New == ymd(20130114)] <- "Zachery Seavey"

# REMOVE OUTLIERS
PhenoList <- 	  c(
  "Outputs.Whole.Arena.Resting.Time",
  "Outputs.Periphery.Permanence.Time",
  "Outputs.Distance.Traveled.Total",
  "Outputs.Center.Average.Speed",
  "Outputs.Distance.Traveled.First.Five.Minutes",
  "Outputs.Number.of.Rears.Total",
  "Outputs.Fecal.Boli",
  "Outputs.PctTime.Corners.Slope",
  "Outputs.PctTime.Center.Slope",
  "Outputs.Distance.Traveled.Slope"
)
# EXPORT THE DATA before outlier removal
write.csv(OFA, file = "OFAOut0.csv",row.names = FALSE)
# Remove outliers from data
OFA <- remove_outliers_main(komp.data = OFA ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
OFA <- addBW(OFA, bw.komp, "OFA", avg_aot = 63)

# also add name of file to col names for later merging
# the [-c(1:3) will not paste on OFA or LD to the first three columns.  I will use these for merging.
colnames(OFA)[-c(1:5)] = paste("OFA_", colnames(OFA)[-c(1:5)], sep = "")

save(OFA,file="OFAout.RData")
# EXPORT THE DATA
write.csv(OFA, file = "OFAOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoOFA"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  "OFA_Outputs.Whole.Arena.Resting.Time",
  "OFA_Outputs.Periphery.Permanence.Time",
  "OFA_Outputs.Distance.Traveled.Total",
  "OFA_Outputs.Center.Average.Speed",
  "OFA_Outputs.Distance.Traveled.First.Five.Minutes",
  "OFA_Outputs.Number.of.Rears.Total",
  "OFA_Outputs.Fecal.Boli",
  "OFA_Outputs.PctTime.Corners.Slope",
  "OFA_Outputs.PctTime.Center.Slope",
  "OFA_Outputs.Distance.Traveled.Slope"
)

DataSet <- OFA

TestDates <- "OFA_Date.of.test.New"
TestTime <- "OFA_Start.time.New"
DataName <- "OFA"
Year <- "OFA_Year"
Month <- "OFA_Month"
DayOfWeek <- "OFA_DayOfWeek"
Hour <- "OFA_Hour"
timeCat <- "OFA_timeCat"
DayOfYear <- "OFA_DayOfYear"
DayOfMonth <- "OFA_DayOfMonth"
StrainGeno <- "OFA_StrainGeno"
predBW <- "OFA_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "OFA_Experimenter.ID",
  "OFA_Arena.ID",
  "OFA_Room.origin",
  "OFA_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- OFA[OFA$OFA_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# LD --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LD <-	loadPheno(paste0(dataDir ,"LD.csv"))
LD <- fixTesterNames(LD,"LD")
LD <- convertVars(data = LD, dict = DictData, "LD")  # change classes of factors from character to the appropriate
# Time and Date
# some of the dates have "/" most have "-" need to fix this
LD$Date.of.Test_temp <- gsub("-","/",LD$Date.of.Test) # replace a "-" with "/" allows parsing in next statement
data1 <- strsplit(LD$Date.of.Test_temp," ")
data1 <- sapply(data1, "[",1)
data1 <- as.data.frame(data1)
colnames(data1) <- "Date.of.Measurement.New"
LD <- cbind(LD, data1)
LD$Date.of.test.New <- ymd(LD$Date.of.Measurement.New)
LD$Start.time.New <- strptime(LD$Start.Time, "%H:%M")
LD <- ExtractDateParameters(LD, 'LD')
LD <- ExtractTimeParameters(LD, 'LD')

# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(LD$Date.of.test.New))))
# Remove all data before 2013 #
LD <- LD[LD$Date.of.test.New > ymd("2012-12-31"),]
# Experimenter.ID "Unknown" on date 11/1/16 - Change to "Pamelia Fraungruber"
LD$Dark.Chamber.Light.Level[LD$Dark.Chamber.Light.Level == ""] <- 4
# Dark.Chamber.Light.Level is "" (empty) in early 2013 then changes to "4" for the rest of the data. I'm changing this all to 4 based on James' feedback.
# REMOVE OUTLIERS
PhenoList <- c(
  'Collected.Values.Pct.Time.in.Dark',
  'Collected.Values.Right.Side.Time.Spent',
  'Collected.Values.Right.Side.Mobile.Time.Spent',
  'Collected.Values.Left.Side.Mobile.Time.Spent',
  'Collected.Values.Side.Changes',
  'Collected.Values.Fecal.Boli',
  'Collected.Values.Reaction.Time'
)
# EXPORT THE DATA before outlier removal
write.csv(LD, file = "LDOut0.csv",row.names = FALSE)
# Remove outliers from data
LD <- remove_outliers_main(komp.data = LD ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
LD <- addBW(LD, bw.komp, "LD", avg_aot = 63)

# also add name of file to col names for later merging
# the [-c(1:3) will not paste on OFA or LD to the first three columns.  I will use these for merging.
colnames(LD)[-c(1:5)] = paste("LD_", colnames(LD)[-c(1:5)], sep = "")

save(LD,file="LDout.RData")
# EXPORT THE DATA
write.csv(LD, file = "LDOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoLD"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  'LD_Collected.Values.Pct.Time.in.Dark',
  'LD_Collected.Values.Right.Side.Time.Spent',
  'LD_Collected.Values.Right.Side.Mobile.Time.Spent',
  'LD_Collected.Values.Left.Side.Mobile.Time.Spent',
  'LD_Collected.Values.Side.Changes',
  'LD_Collected.Values.Fecal.Boli',
  'LD_Collected.Values.Reaction.Time'
)

DataSet <- LD

TestDates <- "LD_Date.of.test.New"
TestTime <- "LD_Start.time.New"
DataName <- "LD"
Year <- "LD_Year"
Month <- "LD_Month"
DayOfWeek <- "LD_DayOfWeek"
Hour <- "LD_Hour"
timeCat <- "LD_timeCat"
DayOfYear <- "LD_DayOfYear"
DayOfMonth <- "LD_DayOfMonth"
StrainGeno <- "LD_StrainGeno"
predBW <- "LD_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "LD_Experimenter.ID",
  "LD_Arena.ID",
  "LD_Room.origin",
  "LD_Dark.Chamber.Light.Level",
  "LD_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- LD[LD$LD_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# TST --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TST <-	loadPheno(paste0(dataDir ,"TST.csv"))
TST <- fixTesterNames(TST,"TST")
TST <- convertVars(data = TST, dict = DictData, "TST")  # change classes of factors from character to the appropriate
# fix Arena.ID
TST$Arena[TST$Arena == 1] <- 3  # according to James -there are only two chambers 1 was labeled as 3 and 2 was labeled as 4
TST$Arena[TST$Arena == 2] <- 4	# this substitutes 1 and 2 with 3 and 4
TST$Date.of.Test_temp <- gsub("-","/",TST$Test.Date) # replace a "-" with "/" allows parsing in next statement
data1 <- strsplit(TST$Date.of.Test_temp," ")
data1 <- sapply(data1, "[",1)
data1 <- as.data.frame(data1)
colnames(data1) <- "Date.of.Measurement.New"
TST <- cbind(TST, data1)
TST$Date.of.test.New <- ymd(TST$Date.of.Measurement.New)
TST$Start.time.New <- strptime(TST$Time.of.Test, "%H:%M")   #mostly 0:00:00 don't use in plotting

# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(TST$Date.of.test.New))))
# Remove all data before 2013 #
TST <- TST[!is.na(TST$Date.of.test.New),]
TST <- TST[TST$Date.of.test.New > ymd("2012-12-31"),]
TST <- ExtractDateParameters(TST, 'TST')
TST <- ExtractTimeParameters(TST, 'TST')
rm(data1)

# REMOVE OUTLIERS
PhenoList <- c(
  'Results.Time.immobile'
)
# EXPORT THE DATA before outlier removal
write.csv(TST, file = "TSTOut0.csv",row.names = FALSE)

# Remove outliers from data
TST <- remove_outliers_main(komp.data = TST ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
TST <- addBW(TST, bw.komp, "TST", avg_aot = 77)

colnames(TST)[-c(1:5)] = paste("TST_", colnames(TST)[-c(1:5)], sep = "")

save(TST,file="TSTout.RData")
# EXPORT THE DATA
write.csv(TST, file = "TSTOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoTST"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  'TST_Results.Time.immobile'
)

DataSet <- TST

TestDates <- "TST_Date.of.test.New"
#TestTime <- "TST_Start.time.New"
TestTime <- TestDates # since most of the times are 0:00:00, I'm just replotting the date as time.
DataName <- "TST"
Year <- "TST_Year"
Month <- "TST_Month"
DayOfWeek <- "TST_DayOfWeek"
Hour <- "TST_Hour"
timeCat <- "TST_timeCat"
DayOfYear <- "TST_DayOfYear"
DayOfMonth <- "TST_DayOfMonth"
StrainGeno <- "TST_StrainGeno"
predBW <- "TST_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "TST_Experimenter",
  "TST_Arena",
  "TST_Room.origin",
  "TST_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
  #timeCat # most of the time data is 0:00:00 so not used
)

control <- TST[TST$TST_StrainGeno == "C57BL/6NJ ",]  ## control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# RR --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
RR <-	read.delim(paste0(dataDir ,"RR.csv"), header = T, sep = "\t")
#	RR has the following column labels that are repeated.  I will remove all these rows.  The actual calculation is from QC.avg.Duration and slope
# Results:Outcome
# Results:Trial
# Results:Duration
RR <- RR[RR$Name != "Results:Outcome",]
RR <- RR[RR$Name != "Results:Trial",]
RR <- RR[RR$Name != "Results:Duration",]
write.table(RR,file = paste0(dataDir,"RR_b.csv"), sep = "\t")  # write the CSV file with these three things removed

RR <- 	loadPheno(paste0(dataDir,"RR_b.csv")) # now use the loadPheno Command with this new dataset
RR <- fixTesterNames(RR,"RR")
RR <- convertVars(data = RR, dict = DictData, "RR")  # change classes of factors from character to the appropriate
RR$Date.of.test.New <- ymd_hms(RR$Test.Date)
RR$Date.of.test.New <- strptime(RR$Date.of.test.New, "%Y-%m-%d")
RR$Start.time.New <- strptime(RR$Time.of.Test, "%H:%M:%S")
RR <- ExtractDateParameters(RR, 'RR')
RR <- ExtractTimeParameters(RR, 'RR')

# some of the arenaID is blank - change this to one since this was very early in the pipeline when there was only one arena
levels(RR$Arena.ID)[levels(RR$Arena.ID) == ""] <- "1"

# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(RR$Date.of.test.New))))
# Remove all data before 2013
RR <- RR[!(is.na(RR$Date.of.test.New)),]
RR <- RR[RR$Date.of.test.New > ymd("2012-12-31"),]

# create a new Equipment column -
#	this combines Arena.ID (the machine that is used) and Chamber.ID
RR$Hardware <- paste(RR$Arena.ID, "_" , RR$Chamber.ID, sep = "")
# Remove Outliers
PhenoList <-  c(
  'QC.Totals.Average.Duration',
  'QC.Totals.Learning.Slope'
)
# EXPORT THE DATA before outlier removal
write.csv(RR, file = "RROut0.csv",row.names = FALSE)

RR <- remove_outliers_main(komp.data = RR ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
RR <- addBW(RR, bw.komp, "RR", avg_aot = 77)

colnames(RR)[-c(1:5)] = paste("RR_", colnames(RR)[-c(1:5)], sep = "")

save(RR, file="RRout.RData")
# EXPORT THE DATA
write.csv(RR, file = "RROut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoRR"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  'RR_QC.Totals.Average.Duration',
  'RR_QC.Totals.Learning.Slope'
)

DataSet <- RR

TestDates <- "RR_Date.of.test.New"
TestTime <- "RR_Start.time.New"
DataName <- "RR"
Year <- "RR_Year"
Month <- "RR_Month"
DayOfWeek <- "RR_DayOfWeek"
Hour <- "RR_Hour"
timeCat <- "RR_timeCat"
DayOfYear <- "RR_DayOfYear"
DayOfMonth <- "RR_DayOfMonth"
StrainGeno <- "RR_StrainGeno"
predBW <- "RR_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "RR_Experimenter",
  "RR_Arena.ID",
  "RR_Chamber.ID",
  "RR_Room.origin",
  "RR_Acceleration",
  "RR_Equipment.Manufacturer",
  "RR_EScell",
  "RR_Hardware"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- RR[RR$RR_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# PPI --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PPI <-	loadPheno(paste0(dataDir ,"PPI.csv"))
PPI <- fixTesterNames(PPI,"PPI")
PPI <- convertVars(data = PPI, dict = DictData, "PPI")  # change classes of factors from character to the appropriate
# Process Time and Date
#PPI$Date.of.test.New <- mdy(PPI$Date)
#formats <- guess_formats(PPI$Date, c("mdY", "BdY", "Bdy", "bdY", "bdy", "mdy", "dby"))
#PPI$Date.of.test.New <- strptime(PPI$Date, formats)
PPI$Date.of.test.New <- MakeSameDateFormat(PPI$Date)
PPI$Start.time.New <- strptime(PPI$Time.of.day, "%I:%M:%S %p")
PPI <- ExtractDateParameters(PPI, 'PPI')
PPI <- ExtractTimeParameters(PPI, 'PPI')

# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(PPI$Date.of.test.New))))
# Remove all data before 2013 #
PPI <- PPI[!(is.na(PPI$Date.of.test.New)),]
PPI <- PPI[PPI$Date.of.test.New > ymd("2012-12-31"),]
# One days in which PPI3, S, and global are very high - remove this data
PPI$Startle.Amplitude.Corrected.Percent.PP3_vmax[PPI$Date.of.test.New == ymd(20160906)] <- NA
PPI$Startle.Amplitude.Corrected.PPI.Global_vmax[PPI$Date.of.test.New == ymd(20160906)] <- NA
PPI$Startle.Amplitude.Corrected.Percent.PP3_vmax[PPI$Date.of.test.New == ymd(20160927)] <- NA
PPI$Startle.Amplitude.Corrected.PPI.Global_vmax[PPI$Date.of.test.New == ymd(20160927)] <- NA
PPI$AVG.Corrected_S_avg[PPI$Date.of.test.New == ymd(20160530)] <- NA

# One Chamber has the id ""Chamber identification" as a level
# These animals ahve no data, there are 20 animals  - sent to James
PPI <- PPI[PPI$Chamber.Used != "Chamber identification", ]
PPI <- PPI[PPI$Chamber.Used != "Chamber identification (8boxes)", ]

# One animal with calibration meter is not assigned.
# Since this animal was tested when calibration meter 1 was used, I labeled this as 1
PPI$Calibration.Meter[PPI$Calibration.Meter == ""] <- 1

# REMOVE OUTLIERS
PhenoList <- 	  c(
  'AVG.Corrected.Percent.PP1_avg',
  'AVG.Corrected.Percent.PP2_avg',
  'AVG.Corrected.Percent.PP3_avg',
  'AVG.Corrected.PPI.Global_avg',
  'AVG.Corrected.Slope.Max_avg',
  'Startle.Amplitude.Corrected.Percent.PP1_vmax',
  'Startle.Amplitude.Corrected.Percent.PP2_vmax',
  'Startle.Amplitude.Corrected.Percent.PP3_vmax',
  'Startle.Amplitude.Corrected.PPI.Global_vmax',
  'Startle.Amplitude.Corrected.Slope.Max_vmax',
  'AVG.Corrected_S_avg',
  'Startle.Amplitude.Corrected_S_vmax'
)
# EXPORT THE DATA before outlier removal
write.csv(PPI, file = "PPIOut0.csv",row.names = FALSE)
# Remove outliers from data
PPI <- remove_outliers_main(komp.data = PPI ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data to PPI data
PPI <- addBW(PPI, bw.komp, "PPI", avg_aot = 70)

colnames(PPI)[-c(1:5)] = paste("PPI_", colnames(PPI)[-c(1:5)], sep = "")

save(PPI, file="PPIout.RData")
# EXPORT THE DATA
write.csv(PPI, file = "PPIOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoPPI"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <-  c(
  'PPI_AVG.Corrected.Percent.PP1_avg',
  'PPI_AVG.Corrected.Percent.PP2_avg',
  'PPI_AVG.Corrected.Percent.PP3_avg',
  'PPI_AVG.Corrected.PPI.Global_avg',
  'PPI_AVG.Corrected.Slope.Max_avg',
  'PPI_AVG.Corrected_S_avg',
  'PPI_Startle.Amplitude.Corrected.Percent.PP1_vmax',
  'PPI_Startle.Amplitude.Corrected.Percent.PP2_vmax',
  'PPI_Startle.Amplitude.Corrected.Percent.PP3_vmax',
  'PPI_Startle.Amplitude.Corrected.PPI.Global_vmax',
  'PPI_Startle.Amplitude.Corrected.Slope.Max_vmax',
  'PPI_Startle.Amplitude.Corrected_S_vmax'
)

DataSet <- PPI

TestDates <- "PPI_Date.of.test.New"
TestTime <- "PPI_Start.time.New"
DataName <- "PPI"
Year <- "PPI_Year"
Month <- "PPI_Month"
DayOfWeek <- "PPI_DayOfWeek"
Hour <- "PPI_Hour"
timeCat <- "PPI_timeCat"
DayOfYear <- "PPI_DayOfYear"
DayOfMonth <- "PPI_DayOfMonth"
StrainGeno <- "PPI_StrainGeno"
predBW <- "PPI_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "PPI_Room.origin",
  "PPI_Number.of.Trials",
  "PPI_Chamber.Used",
  "PPI_Calibration.Meter",
  "PPI_Experimenter",
  "PPI_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- PPI[PPI$PPI_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# ECT --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ECT <-	loadPheno(paste0(dataDir ,"ECT.csv"))
ECT <- fixTesterNames(ECT,"ECT") # fix names before converting the variables
ECT <- convertVars(data = ECT, dict = DictData, "ECT")  # change classes of factors from character to the appropriate
# New DB output fro Mike has dates in one format (2016-07-01)
ECT$Date.of.test.New <- ymd_hms(ECT$ECT.date1)
ECT <- ExtractDateParameters(ECT, 'ECT')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(ECT$Date.of.test.New))))
# Remove all data before 2013 #
#ECT <- ECT[ECT$Strain.Name != "C57BL/6NJ" | !is.na(ECT$Date.of.test.New),]
ECT <- ECT[!(is.na(ECT$Date.of.test.New)),]
ECT <- ECT[ECT$Date.of.test.New > ymd("2012-12-31"),]
#Fix the levels of Equipment.model to remove space from end
levels(ECT$Equipment.model)[levels(ECT$Equipment.model)=="ECT unit 7801 (#25138) "] <- "ECT unit 7801 (#25138)"

# REMOVE OUTLIERS
PhenoList <- c(
  'ECT.mA.Threshold'
)
# EXPORT THE DATA before outlier removal
write.csv(ECT, file = "ECTOut0.csv",row.names = FALSE)
# Remove outliers from data
ECT <- remove_outliers_main(komp.data = ECT ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  weight data
ECT <- addBW(ECT, bw.komp, "ECT", avg_aot = 119)

colnames(ECT)[-c(1:5)] = paste("ECT_", colnames(ECT)[-c(1:5)], sep = "")

save(ECT, file="ECTout.RData")
# EXPORT THE DATA
write.csv(ECT, file = "ECTOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoECT"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  'ECT_ECT.mA.Threshold'
)

DataSet <- ECT

TestDates <- "ECT_Date.of.test.New"
#TestTime <- "ECT_Start.time.New"
DataName <- "ECT"
Year <- "ECT_Year"
Month <- "ECT_Month"
DayOfWeek <- "ECT_DayOfWeek"
#Hour <- "ECT_Hour"
#timeCat <- "ECT_timeCat"
DayOfYear <- "ECT_DayOfYear"
DayOfMonth <- "ECT_DayOfMonth"
StrainGeno <- "ECT_StrainGeno"
predBW <- "ECT_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  #"ECT_Experimenter",
  #"ECT_Equipmnent",
  "ECT_Equipment.model",
  "ECT_Room.origin",
  "ECT_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
  #timeCat
)

control <- ECT[ECT$ECT_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
source(plotAllFile2)

# GRIP --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GRIP <-	loadPheno(paste0(dataDir ,"GRIP.csv"))
GRIP <- fixTesterNames(GRIP, "GRIP")
GRIP <- convertVars(data = GRIP, dict = DictData, "GRIP")  # change classes of factors from character to the appropriate
# Some of the dates and time in data are in different formats.  I fixed this in Excel GRIP.csv files.
# Genotype '-Dcaf10<tm1.1(KOMP)Vlcg> -/-' '-Dcaf10<tm1.1(KOMP)Vlcg> -/-' two instances that need to be fixed. fixed in excel.
# Mouse.ID A12442 has time of test at 23:31, removed by hand in GRIP.csv file.
GRIP$Date.of.test.New <- ymd_hms(GRIP$Test.date)
GRIP$Start.time.New <- strptime(GRIP$Time.of.Day,"%H:%M")  # creates a new date of test col - use for plotting
GRIP <- ExtractDateParameters(GRIP, 'GRIP')
GRIP <- ExtractTimeParameters(GRIP, 'GRIP')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(GRIP$Date.of.test.New))))
# Remove all data before 2013 ##
#GRIP <- GRIP[GRIP$Strain.Name != "C57BL/6NJ" | !is.na(GRIP$Date.of.test.New),]
GRIP <- GRIP[!(is.na(GRIP$Date.of.test.New)),]
GRIP <- GRIP[GRIP$Date.of.test.New > ymd("2012-12-31"),]
# Fix equipment inconsistencies - see email exchange between Jenn and Me  add ROOM OF TEST
GRIP$equip.New[GRIP$Date.of.test.New < ymd(20130805)] <- "DFIS 2 Elissa" # all tests before 8/14/2014 was on Elissa's DFIS
GRIP$equip.New[GRIP$Date.of.test.New >= ymd(20130805) & GRIP$Date.of.test.New < ymd(20131114)] <- "DFIS 2 Physio"
GRIP$equip.New[GRIP$Date.of.test.New >= ymd(20131114) & GRIP$Date.of.test.New < ymd(20131206)] <- "DFE 2 Beep"
GRIP$equip.New[GRIP$Date.of.test.New > ymd(20131206)] <- "DFE 2"
GRIP$TestRoom[GRIP$Date.of.test.New < ymd(20140205)] <- "G3140B"
GRIP$TestRoom[GRIP$Date.of.test.New >= ymd(20140205)] <- "G3120C"

# REMOVE OUTLIERS
PhenoList <- 	c(
  'Results.Fore.and.hind.limb.grip.mean',
  'Results.fore.and.hind.limb.grip.mean.body.weight'
)
# EXPORT THE DATA before outlier removal
write.csv(GRIP, file = "GRIPOut0.csv",row.names = FALSE)
# Remove outliers from data
GRIP <- remove_outliers_main(komp.data = GRIP ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
GRIP <- addBW(GRIP, bw.komp, "GRIP", avg_aot = 56)

colnames(GRIP)[-c(1:5)] = paste("GRIP_", colnames(GRIP)[-c(1:5)], sep = "")

save(GRIP, file="GRIPout.RData")
# EXPORT THE DATA
write.csv(GRIP, file = "GRIPOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoGRIP"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- 	c(
  'GRIP_Results.Fore.and.hind.limb.grip.mean',
  'GRIP_Results.fore.and.hind.limb.grip.mean.body.weight'
)

DataSet <- GRIP

TestDates <- "GRIP_Date.of.test.New"
TestTime <- "GRIP_Start.time.New"
DataName <- "GRIP"
Year <- "GRIP_Year"
Month <- "GRIP_Month"
DayOfWeek <- "GRIP_DayOfWeek"
Hour <- "GRIP_Hour"
timeCat <- "GRIP_timeCat"
DayOfYear <- "GRIP_DayOfYear"
DayOfMonth <- "GRIP_DayOfMonth"
StrainGeno <- "GRIP_StrainGeno"
predBW <- "GRIP_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "GRIP_Experimenter",
  "GRIP_equip.New",
  "GRIP_Date.last.calibrated",
  "GRIP_Room.origin",
  "GRIP_TestRoom",
  "GRIP_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- GRIP[GRIP$GRIP_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# Hematology  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Hematology <-	loadPheno(paste0(dataDir ,"Hematology.csv"))
# several animals have duplicated data - emailed James and Rick - fixed them manually in csv file.
# C5229 had two different sets of data.  I removed this completely
Hematology$Analyst.ID[Hematology$Analyst.ID == "SC"] <- "Steve Ciciotte"
Hematology$Experimenter.ID[Hematology$Experimenter.ID == "OH"] <- "Olivia Hon"
Hematology <- convertVars(data = Hematology, dict = DictData, "Hematology")  # change classes of factors from character to the appropriate
#	Date and Time Conversion
Hematology$Date.of.test.New <- ymd_hms(Hematology$Date.and.Time.of.Blood.Collection)  # creates a new date of test col - use for plotting
# Extract time from Date and Time - two step use strftime to convert to character and then convert back to date time, but with same date
Hematology$Time.of.Blood.Collection.New <- strftime(Hematology$Date.of.test.New, format="%H:%M:%S")
Hematology$Start.time.New <- strptime(Hematology$Time.of.Blood.Collection.New, format="%H:%M")
Hematology$Date.of.test.New <- ymd(as.Date(Hematology$Date.of.test.New))
Hematology <- ExtractDateParameters(Hematology, 'Hematology')
Hematology <- ExtractTimeParameters(Hematology, 'Hematology')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(Hematology$Date.of.test.New))))
# Remove all data before 2013 ##
Hematology <- Hematology[!is.na(Hematology$Date.of.test.New),]
Hematology <- Hematology[Hematology$Date.of.test.New > ymd("2012-12-31"),]
# REMOVE OUTLIERS
PhenoList <- 	  c('Hematology.Results.Mean.corpuscular.hemoglobin..CHg.',
                  'Hematology.Results.Mean.Cell.Volume..MCV.',
                  'Hematology.Results.Red.Cell.Distr..Width..RDW.',
                  'Hematology.Results.Red.Cell.Hem..Con..Mean..CHCM.',
                  'Hematology.Results.Lymphocytes..LYM.',
                  'Hematology.Results.Red.Blood.Cells..RBC.',
                  'Hematology.Results.Hemoglobin.Concentration.Distr..Width..HDW.',
                  'Hematology.Results.Mean.Platelet.Volume..MPV.',
                  'Hematology.Results.Measured.Hemoglobin..mHGB.',
                  'Hematology.Results.Platelet.Count..PLT.',
                  'Hematology.Results.Percent.Retic',
                  'Hematology.Results.Eosinophils..EOS.',
                  'Hematology.Results.White.Blood.Cells..WBC.',
                  'Hematology.Results.Neutrophils..NEUT.',
                  'Hematology.Results.Monocytes..MONO.'
)
# EXPORT THE DATA before outlier removal
write.csv(Hematology, file = "HematologyOut0.csv",row.names = FALSE)
# Remove outliers from data
Hematology <- remove_outliers_main(komp.data = Hematology ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
Hematology <- addBW(Hematology, bw.komp, "Hematology", avg_aot = 126)

# this has to be done before modifying the columnnames below.
colnames(Hematology)[-c(1:5)] = paste("Hematology_", colnames(Hematology)[-c(1:5)], sep = "")

save(Hematology, file="Hematologyout.RData")
# EXPORT THE DATA
write.csv(Hematology, file = "HematologyOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoHematology"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- 		c('Hematology_Hematology.Results.Mean.corpuscular.hemoglobin..CHg.',
                 'Hematology_Hematology.Results.Mean.Cell.Volume..MCV.',
                 'Hematology_Hematology.Results.Red.Cell.Distr..Width..RDW.',
                 'Hematology_Hematology.Results.Red.Cell.Hem..Con..Mean..CHCM.',
                 'Hematology_Hematology.Results.Lymphocytes..LYM.',
                 'Hematology_Hematology.Results.Red.Blood.Cells..RBC.',
                 'Hematology_Hematology.Results.Hemoglobin.Concentration.Distr..Width..HDW.',
                 'Hematology_Hematology.Results.Mean.Platelet.Volume..MPV.',
                 'Hematology_Hematology.Results.Measured.Hemoglobin..mHGB.',
                 'Hematology_Hematology.Results.Platelet.Count..PLT.',
                 #'Hematology_Hematology.Results.Percent.Retic', too many NA's remove
                 'Hematology_Hematology.Results.Eosinophils..EOS.',
                 'Hematology_Hematology.Results.White.Blood.Cells..WBC.',
                 'Hematology_Hematology.Results.Neutrophils..NEUT.',
                 'Hematology_Hematology.Results.Monocytes..MONO.'
)

DataSet <- Hematology

TestDates <- "Hematology_Date.of.test.New"
TestTime <- "Hematology_Start.time.New"
DataName <- "Hematology"
Year <- "Hematology_Year"
Month <- "Hematology_Month"
DayOfWeek <- "Hematology_DayOfWeek"
Hour <- "Hematology_Hour"
timeCat <- "Hematology_timeCat"
DayOfYear <- "Hematology_DayOfYear"
DayOfMonth <- "Hematology_DayOfMonth"
StrainGeno <- "Hematology_StrainGeno"
predBW <- "Hematology_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "Hematology_Room.origin",
  "Hematology_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
  #timeCat	 # time is not accurate in this field.
)

# Hematology_Analyst.ID not used only one right now
# Hematology_Experimenter.ID not used - only one
control <- Hematology[Hematology$Hematology_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# GTT  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GTT <-	loadPheno(paste0(dataDir ,"GTT.csv"))
GTT <- fixTesterNames(GTT, "GTT")
# No Test and time of date in this assay - use DateReceived field
GTT$Date.of.test.New <- ymd(GTT$DateReceived)
GTT <- ExtractDateParameters(GTT, 'GTT')
GTT <- convertVars(data = GTT, dict = DictData, "GTT")  # change classes of factors from character to the appropriate
# AUC using trapezoid method [(y1 + y2)/2] * Time] subtract the T0 *120 from this for account for baseline.
GTT$AUC <- (((GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.0 + GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.15.minutes)/2)*15) +
  (((GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.15.minutes + GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.30.minutes)/2)*15) +
  (((GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.30.minutes + GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.60.minutes)/2)*30) +
  (((GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.60.minutes + GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.120.minutes)/2)*60) -
  GTT$GTT.Data.Fields.Plasma.glucose.level.at.Time.0 * 120
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(GTT$Date.of.test.New))))
# Remove all data before 2013 ##
GTT <- GTT[GTT$Date.of.test.New > ymd("2012-12-31"),]
PhenoList <- 	c(#		'GTT.Data.Fields.Amount.glucose.injected',
  #'GTT.Data.Fields.Body.weight..g.',
  'GTT.Data.Fields.Plasma.glucose.level.at.Time.0', # this is fasting glucose
  'GTT.Data.Fields.Plasma.glucose.level.at.Time.15.minutes', # this is initial response
  'GTT.Data.Fields.Plasma.glucose.level.at.Time.60.minutes',
  'GTT.Data.Fields.Plasma.glucose.level.at.Time.120.minutes',
  'AUC'
)
# EXPORT THE DATA before outlier removal
write.csv(GTT, file = "GTTOut0.csv",row.names = FALSE)
# Remove outliers from data
GTT <- remove_outliers_main(komp.data = GTT ,list = PhenoList, stdv = 2, iqr = TRUE)
# Add  BW data
# GTT <- addBW(GTT, bw.komp, "GTT")
# BW is not used as a fator for modeling here  -the BW was collected as part of the measure.
# In addition - since there is no date of test, the predicted BW will be wrong
# this has to be done before modifying the columnnames below.
colnames(GTT)[-c(1:5)] = paste("GTT_", colnames(GTT)[-c(1:5)], sep = "")

save(GTT, file="GTTout.RData")
# EXPORT THE DATA
write.csv(GTT, file = "GTTOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoGTT"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- 	c(#'GTT_GTT.Data.Fields.Amount.glucose.injected',
  #'GTT_GTT.Data.Fields.Body.weight..g.',
  'GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.0', # this is fasting glucose
  'GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.15.minutes', # this is initial response
  'GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.60.minutes',
  'GTT_GTT.Data.Fields.Plasma.glucose.level.at.Time.120.minutes',
  'GTT_AUC'
)

DataSet <- GTT

TestDates <- "GTT_Date.of.test.New"
TestTime <- "GTT_Test.Time"
DataName <- "GTT"
Year <- "GTT_Year"
Month <- "GTT_Month"
DayOfWeek <- "GTT_DayOfWeek"
#Hour <- "GTT_Hour"
#timeCat <- "GTT_timeCat"
DayOfYear <- "GTT_DayOfYear"
DayOfMonth <- "GTT_DayOfMonth"
StrainGeno <- "GTT_StrainGeno"
predBW <- "GTT_GTT.Data.Fields.Body.weight..g."

ExperimentalFactors <- c(
  "Sex",
  "GTT_Experimenter.ID", # only Olivia - so not used in model
  #"GTT_Equipment.Name", # Only one equipment used
  "GTT_Room.origin",
  "GTT_Mouse.Restrained",
  "GTT_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
)

# GTT_Experimenter.ID not used - only one
control <- GTT[GTT$GTT_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# HeartWeight  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HeartWeight <-	loadPheno(paste0(dataDir ,"HeartWeight.csv"))
HeartWeight <- fixTesterNames(HeartWeight, "HeartWeight")
HeartWeight <- convertVars(data = HeartWeight, dict = DictData, "HeartWeight")  # change classes of factors from character to the appropriate
#Date is in multiple formats with and without hours.  # this is not true with Mike's newest export, the dates look uniform
# no time information
#HeartWeight$Date.of.experiment.New <- gsub("\\ .*","",HeartWeight$Date.of.experiment)  # replace delete everything after space with nothing
#HeartWeight$Date.of.test.New <- strptime(HeartWeight$Date.of.experiment.New,"%m/%d/%y") ## convert to date
HeartWeight$Date.of.test.New <- strptime(HeartWeight$Date.of.experiment, "%Y-%m-%d")
HeartWeight <- ExtractDateParameters(HeartWeight, 'HeartWeight')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(HeartWeight$Date.of.test.New))))
HeartWeight <- HeartWeight[!is.na(HeartWeight$Date.of.test.New),]
# K56071 had a heart weight that was an outlier - removed manually in HeartWeight.csv file.
# heart weight prior to 2013-11-01 is highly variable - remove all data prior to 2013-11-01
HeartWeight <- HeartWeight[HeartWeight$Date.of.test.New > "2013-10-31",]  # remove 2013-11-01 data
# Remove all data before 2013 ##
HeartWeight <- HeartWeight[HeartWeight$Date.of.test.New > "2012-12-31",]

# Remove Outliers
PhenoList <- 	c(
  'Weights.Heart.Weight',
  'Weights.HW.BW'
)
# EXPORT THE DATA before outlier removal
write.csv(HeartWeight, file = "HeartWeightOut0.csv",row.names = FALSE)
HeartWeight <- remove_outliers_main(komp.data = HeartWeight ,list = PhenoList, stdv = 2, iqr = TRUE)
# Add  BW data
# HeartWeight <- addBW(HeartWeight, bw.komp, "HeartWeight")
# predicted BW is not used - the actual BW at time of data collection is used.
colnames(HeartWeight)[-c(1:5)] = paste("HeartWeight_", colnames(HeartWeight)[-c(1:5)], sep = "")

save(HeartWeight, file="HeartWeightout.RData")
# EXPORT THE DATA
write.csv(HeartWeight, file = "HeartWeightOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoHeartWeight"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- 	c(
  'HeartWeight_Weights.Heart.Weight',
  'HeartWeight_Weights.HW.BW'
)

DataSet <- HeartWeight

TestDates <- "HeartWeight_Date.of.test.New"
TestTime <- "HeartWeight_Date.of.test.New"
DataName <- "HeartWeight"
Year <- "HeartWeight_Year"
Month <- "HeartWeight_Month"
DayOfWeek <- "HeartWeight_DayOfWeek"
Hour <- "HeartWeight_Hour"
timeCat <- "HeartWeight_timeCat"
DayOfYear <- "HeartWeight_DayOfYear"
DayOfMonth <- "HeartWeight_DayOfMonth"
StrainGeno <- "HeartWeight_StrainGeno"
predBW <- "HeartWeight_Weights.Body.Weight" ## use the BW at time of dissection

ExperimentalFactors <- c(
  "Sex",
  "HeartWeight_Room.origin",
  "HeartWeight_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
  #timeCat
)
# HeartWeight_Experimenter not used only one
control <- HeartWeight[HeartWeight$HeartWeight_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# Insulin  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Insulin <-	loadPheno(paste0(dataDir ,"Insulin.csv"))
Insulin <- fixTesterNames(Insulin, "Insulin")
Insulin <- convertVars(data = Insulin, dict = DictData, "Insulin")  # change classes of factors from character to the appropriate
# no time info
# fix dates - remove time and then formate as date
# fixed in Insulin.csv file.  There were Date and Time.of.blood, Date.of.Measurement in 2018 and 2915
# Insulin$Date.and.Time.of.blood.New <- gsub("\\ .*","",Insulin$Date.and.Time.of.blood)  ## replace everything after space with nothing
Insulin$Date.of.test.New <- parse_date_time(Insulin$Date.and.Time.of.blood, orders="ymd_hms")
Insulin$Date.of.test.New <- strptime(Insulin$Date.of.test.New, "%Y-%m-%d")
#Insulin$Date.of.test.New <- parse_date_time(Insulin$Date.of.Measurement, orders="mdy") ## dont use this date
Insulin <- ExtractDateParameters(Insulin, 'Insulin')
# Remove all 0 containing values - these must be mistakes?
Insulin <- Insulin[Insulin$Insulin.Value.Insulin != 0,]
# insulin is exponentially distributed - not normal at all
# log transform values
Insulin$Insulin.Value.Insulin.LOG <- log(Insulin$Insulin.Value.Insulin)  # make things non zero
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(Insulin$Date.of.test.New))))
# Remove all data before 2013 ##
Insulin <- Insulin[Insulin$Date.of.test.New > "2012-12-31",]
# Remove Outliers
PhenoList <- c(
  'Insulin.Value.Insulin',
  'Insulin.Value.Insulin.LOG'
)
# EXPORT THE DATA before outlier removal
write.csv(Insulin, file = "InsulinOut0.csv",row.names = FALSE)
Insulin <- remove_outliers_main(komp.data = Insulin ,list = PhenoList, stdv = 2, iqr = TRUE)
# Add  BW data
Insulin <- addBW(Insulin, bw.komp, "Insulin", avg_aot = 126)

colnames(Insulin)[-c(1:5)] = paste("Insulin_", colnames(Insulin)[-c(1:5)], sep = "")

save(Insulin, file="Insulinout.RData")
# EXPORT THE DATA
write.csv(Insulin, file = "InsulinOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoInsulin"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <- c(
  'Insulin_Insulin.Value.Insulin',
  'Insulin_Insulin.Value.Insulin.LOG'
)

DataSet <- Insulin

TestDates <- "Insulin_Date.of.test.New"
TestTime <- "Insulin_Date.of.test.New"
DataName <- "Insulin"
Year <- "Insulin_Year"
Month <- "Insulin_Month"
DayOfWeek <- "Insulin_DayOfWeek"
#Hour <- "Insulin_Hour"
#timeCat <- "Insulin_timeCat"
DayOfYear <- "Insulin_DayOfYear"
DayOfMonth <- "Insulin_DayOfMonth"
StrainGeno <- "Insulin_StrainGeno"
predBW <- "Insulin_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "Insulin_Blood.Collection.Experimenter",
  "Insulin_Room.origin",
  "Insulin_Kit.Lot.Number",
  "Insulin_Sample.Status", # fresh or frozen
  "Insulin_Plasma.Dilution",
  "Insulin_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
)

control <- Insulin[Insulin$Insulin_StrainGeno == "C57BL/6NJ ",]  ## control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# ClinicalBloodChemistry  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ClinicalBloodChemistry <-	loadPheno(paste0(dataDir ,"ClinicalBloodChenistry.csv"))
ClinicalBloodChemistry <- fixTesterNames(ClinicalBloodChemistry, "ClinicalBloodChemistry")
ClinicalBloodChemistry <- convertVars(data = ClinicalBloodChemistry, dict = DictData, "ClinicalBloodChemistry")  # change classes of factors from character to the appropriate
ClinicalBloodChemistry$Date.and.time.of.Blood.Collection.New <- ymd_hms(ClinicalBloodChemistry$Date.and.time.of.Blood.Collection)
# extract the date only
ClinicalBloodChemistry$Date.of.test.New <- round_date(ClinicalBloodChemistry$Date.and.time.of.Blood.Collection.New, unit = "day")
# Date.of.Measurement has date and time OR just date, extract just the date
#time  - not accurate
data1 <- strsplit(as.character(ClinicalBloodChemistry$Date.and.time.of.Blood.Collection.New)," ")
data1 <- sapply(data1, "[",2)
data1<- strptime(data1, "%H:%M:%S")
data1 <- as.data.frame(data1)
colnames(data1) <- "Start.time.New"
ClinicalBloodChemistry <- cbind(ClinicalBloodChemistry, data1)
ClinicalBloodChemistry <- ExtractDateParameters(ClinicalBloodChemistry, 'ClinicalBloodChemistry')
ClinicalBloodChemistry <- ExtractTimeParameters(ClinicalBloodChemistry, 'ClinicalBloodChemistry')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(ClinicalBloodChemistry$Date.of.test.New))))
# Remove all data before 5/3/13 ## This is different than the rest of the phenotypes. The data is not consistant up to this point.
#ClinicalBloodChemistry <- ClinicalBloodChemistry[ClinicalBloodChemistry$Strain.Name != "C57BL/6NJ" | !is.na(ClinicalBloodChemistry$Date.of.test.New),]
ClinicalBloodChemistry <- ClinicalBloodChemistry[!(is.na(ClinicalBloodChemistry$Date.of.test.New)),]
ClinicalBloodChemistry <- ClinicalBloodChemistry[ClinicalBloodChemistry$Date.of.test.New > ymd("2013-05-03"),]
# Remove all data after 2016 ##
#ClinicalBloodChemistry <- ClinicalBloodChemistry[ClinicalBloodChemistry$Date.of.test.New < ymd("2017-01-01"),]
#	Remove certain dates of protein that has very high values
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20151008)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20150716)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20151119)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20151124)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160304)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160317)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160414)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160421)] <- NA
ClinicalBloodChemistry$Results.Total.Protein..TPROT.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160428)] <- NA
#	Remove certain dates of Glucose that has very high values
ClinicalBloodChemistry$Results.Glucose..GLU.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160304)] <- NA
ClinicalBloodChemistry$Results.Glucose..GLU.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160317)] <- NA
ClinicalBloodChemistry$Results.Glucose..GLU.[ClinicalBloodChemistry$Date.of.test.New == ymd(20160331)] <- NA

# REMOVE OUTLIERS
PhenoList <-  c(
  'Results.Alanine.aminotransferase...ALT.',
  'Results.Albumin..ALB.',
  'Results.Alkaline.Phosphatase',
  'Results.Aspartate.aminotransferase...AST.',
  'Results.AST.ALT.Ratio',
  'Results.Bilirubin..TBIL.',
  'Results.Calcium..Ca.',
  #'Results.Creatinine.', # mostly empty
  'Results.Carbon.Dioxide', #
  'Results.Chloride',#
  'Results.Free.Fatty.Acids..NEFA.',#
  'Results.Glucose..GLU.',
  'Results.HDL.Cholesterol..HDLD.',
  'Results.HDLD.TCHOL.Ratio',
  'Results.Iron..Fe.', #
  'Results.Phosphorous..Phos.',
  'Results.Potassium',#
  'Results.Sodium',#
  'Results.Total.Cholesterol..T.CHOL.',
  'Results.Total.Protein..TPROT.',
  'Results.Triglycerides..TG.',
  'Results.Urea..BUN.'
)
# EXPORT THE DATA before outlier removal
write.csv(ClinicalBloodChemistry, file = "ClinicalBloodChemistryOut0.csv",row.names = FALSE)
ClinicalBloodChemistry <- remove_outliers_main(komp.data = ClinicalBloodChemistry ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
ClinicalBloodChemistry <- addBW(ClinicalBloodChemistry, bw.komp, "ClinicalBloodChemistry", avg_aot = 126)

colnames(ClinicalBloodChemistry)[-c(1:5)] = paste("ClinicalBloodChemistry_", colnames(ClinicalBloodChemistry)[-c(1:5)], sep = "")

save(ClinicalBloodChemistry, file="ClinicalBloodChemistryout.RData")
# EXPORT THE DATA
write.csv(ClinicalBloodChemistry, file = "ClinicalBloodChemistryOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoClinicalBloodChemistry"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <-  c(
  'ClinicalBloodChemistry_Results.Alanine.aminotransferase...ALT.',
  'ClinicalBloodChemistry_Results.Albumin..ALB.',
  'ClinicalBloodChemistry_Results.Alkaline.Phosphatase',
  'ClinicalBloodChemistry_Results.Aspartate.aminotransferase...AST.',
  'ClinicalBloodChemistry_Results.AST.ALT.Ratio',
  'ClinicalBloodChemistry_Results.Bilirubin..TBIL.',
  'ClinicalBloodChemistry_Results.Calcium..Ca.',
  #'ClinicalBloodChemistry_Results.Creatinine.',
  'ClinicalBloodChemistry_Results.Glucose..GLU.',
  'ClinicalBloodChemistry_Results.HDL.Cholesterol..HDLD.',
  'ClinicalBloodChemistry_Results.HDLD.TCHOL.Ratio',
  'ClinicalBloodChemistry_Results.Phosphorous..Phos.',
  'ClinicalBloodChemistry_Results.Total.Cholesterol..T.CHOL.',
  'ClinicalBloodChemistry_Results.Total.Protein..TPROT.',
  'ClinicalBloodChemistry_Results.Triglycerides..TG.',
  'ClinicalBloodChemistry_Results.Urea..BUN.'
)

DataSet <- ClinicalBloodChemistry

TestDates <- "ClinicalBloodChemistry_Date.of.test.New"
TestTime <- "ClinicalBloodChemistry_Start.time.New"
DataName <- "ClinicalBloodChemistry"
Year <- "ClinicalBloodChemistry_Year"
Month <- "ClinicalBloodChemistry_Month"
DayOfWeek <- "ClinicalBloodChemistry_DayOfWeek"
Hour <- "ClinicalBloodChemistry_Hour"
timeCat <- "ClinicalBloodChemistry_timeCat"
DayOfYear <- "ClinicalBloodChemistry_DayOfYear"
DayOfMonth <- "ClinicalBloodChemistry_DayOfMonth"
StrainGeno <- "ClinicalBloodChemistry_StrainGeno"
predBW <- "ClinicalBloodChemistry_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "ClinicalBloodChemistry_Room.origin",
  "ClinicalBloodChemistry_Storage.temperature.from.blood.collection.till.measurement",
  "ClinicalBloodChemistry_Hemolysis.Status..enum.",
  "ClinicalBloodChemistry_Equipment.Model",
  "ClinicalBloodChemistry_EScell"
)

# ClinicalBloodChemistry_Analyst.ID not used  - only Steve Cicciotte
# ClinicalBloodChemistry_Experimenter.ID not used - only Olivia Hon

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- ClinicalBloodChemistry[ClinicalBloodChemistry$ClinicalBloodChemistry_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# BodyComp  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BodyComp <-	loadPheno(paste0(dataDir ,"BodyComp.csv"))
BodyComp <- fixTesterNames(BodyComp, "BodyComp")
BodyComp <- convertVars(data = BodyComp, dict = DictData, "BodyComp")  # change classes of factors from character to the appropriate
# extract the date from image file name
# PIXIMUS:Image	9353 08-20-12 07'59'38.img
data1 <- strsplit(as.character(BodyComp$PIXIMUS.Image)," ")
data1 <- sapply(data1, "[",2)
data1<- mdy(data1)
data1 <- as.data.frame(data1)
colnames(data1) <- "Date.of.test.New"
BodyComp <- cbind(BodyComp, data1)
BodyComp <- ExtractDateParameters(BodyComp, 'BodyComp')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(BodyComp$Date.of.test.New))))
# Remove all data before 2013 #
BodyComp <- BodyComp[!is.na(BodyComp$Date.of.test.New),]
BodyComp <- BodyComp[BodyComp$Date.of.test.New > ymd("2012-12-31"),]
# REMOVE OUTLIERS
PhenoList <-  c(
  'PIXIMUS.B.Area',
  'PIXIMUS.BMC',
  'PIXIMUS.BMC..Body.weight',
  'PIXIMUS.BMD',
  'PIXIMUS.Bone.Area..BMC.BMD.',
  'PIXIMUS.Fat.mass',
  'PIXIMUS.Fat..Body.weight',
  'PIXIMUS.Lean.mass',
  'PIXIMUS.Lean..Body.weight',
  'PIXIMUS.Percent.Fat',
  'PIXIMUS.RST',
  'PIXIMUS.Subject.Length',
  'PIXIMUS.Subject.Weight',
  'PIXIMUS.TTM'
)
# EXPORT THE DATA before outlier removal
write.csv(BodyComp, file = "BodyCompOut0.csv",row.names = FALSE)
BodyComp <- remove_outliers_main(komp.data = BodyComp ,list = PhenoList, stdv = 2, iqr = TRUE)
# Add  BW data
# BodyComp has a BW field (subject weight) I will use the fitted predicted value here.
BodyComp <- addBW(BodyComp, bw.komp, "BodyComp", avg_aot = 98)

colnames(BodyComp)[-c(1:5)] = paste("BodyComp_", colnames(BodyComp)[-c(1:5)], sep = "")

save(BodyComp, file="BodyCompout.RData")
# EXPORT THE DATA
write.csv(BodyComp, file = "BodyCompOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoBodyComp"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <-  c(
  'BodyComp_PIXIMUS.B.Area',
  'BodyComp_PIXIMUS.BMC',
  'BodyComp_PIXIMUS.BMC..Body.weight',
  'BodyComp_PIXIMUS.BMD',
  'BodyComp_PIXIMUS.Bone.Area..BMC.BMD.',
  'BodyComp_PIXIMUS.Fat.mass',
  'BodyComp_PIXIMUS.Fat..Body.weight',
  'BodyComp_PIXIMUS.Lean.mass',
  'BodyComp_PIXIMUS.Lean..Body.weight',
  #'BodyComp_PIXIMUS...Fat',
  'BodyComp_PIXIMUS.RST',
  'BodyComp_PIXIMUS.Subject.Length',
  'BodyComp_PIXIMUS.Subject.Weight',
  'BodyComp_PIXIMUS.TTM'
)

DataSet <- BodyComp

TestDates <- "BodyComp_Date.of.test.New"
TestTime <- 9 # no test time - this will force a time of 9
DataName <- "BodyComp"
Year <- "BodyComp_Year"
Month <- "BodyComp_Month"
DayOfWeek <- "BodyComp_DayOfWeek"
#Hour <- 9 # no test time - this will force a time of 9
#timeCat <- "Morning"
DayOfYear <- "BodyComp_DayOfYear"
DayOfMonth <- "BodyComp_DayOfMonth"
StrainGeno <- "BodyComp_StrainGeno"
predBW <- "BodyComp_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "BodyComp_Room.origin",
  "BodyComp_Experimenter.ID",
  "BodyComp_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek
  #timeCat
)

control <- BodyComp[BodyComp$BodyComp_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# SLEEP  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# when the new code is done  - remove any animals with confidence level less than 0.6
# although it would be interesting to make sure all animals sleep
SLEEP <-	loadPheno(paste0(dataDir ,"SLEEP.csv"))
SLEEP <- fixTesterNames(SLEEP,"SLEEP")
SLEEP <- convertVars(data = SLEEP, dict = DictData, "SLEEP")
levels(SLEEP$Software)[levels(SLEEP$Software)=="Signal Solutions Sensory Systems, Lexington KY"] <- "SignalSolutions"
levels(SLEEP$Software)[levels(SLEEP$Software)=="SleepStats, MatLab, MouseRec"] <- "SleepStats"
SLEEP$Start.Date <- MakeSameDateFormat(SLEEP$Start.Date)
SLEEP$Date.of.test.New <- ymd(SLEEP$Start.Date)
#time does not matter
#SLEEP$Start.Time <- MakeSameTimeFormat(SLEEP$Start.Time)
#SLEEP$Start.time.New <- SLEEP$Start.Time
SLEEP <- ExtractDateParameters(SLEEP, 'SLEEP')
#SLEEP <- ExtractTimeParameters(SLEEP, 'SLEEP')
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(SLEEP$Date.of.test.New))))
# Remove all data before 2013 ##
SLEEP <- SLEEP[!is.na(SLEEP$Date.of.test.New),]
SLEEP <- SLEEP[SLEEP$Date.of.test.New > ymd("2012-12-31"),]

## Remove outliers from data
PhenoList <-  c(
  'Results.Sleep.Light.Phase.Percent',
  'Results.Sleep.Daily.Percent',
  'Results.Breath.Rate.During.Sleep.Mean',
  'Results.Breath.Rate.During.Sleep.Standard.Deviation',
  'Results.Light.Onset.Wake.Median',
  'Results.Diurnal.Wake.Ratio.Median',
  'Results.Sleep.Bout.Lengths.Mean',
  'Results.Light.Sleep.Bout.Lengths.Mean',
  'Results.Peak.Wake.wrt.Dark.Onset.Median',
  'Results.Light.Sleep.Bout.Lengths.Standard.Deviation'
)
# EXPORT THE DATA before outlier removal
write.csv(SLEEP, file = "SLEEPOut0.csv",row.names = FALSE)
SLEEP <- remove_outliers_main(komp.data = SLEEP ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data -
SLEEP <- addBW(SLEEP, bw.komp, "SLEEP", avg_aot = 105)

colnames(SLEEP)[-c(1:5)] = paste("SLEEP_", colnames(SLEEP)[-c(1:5)], sep = "")

save(SLEEP, file="SLEEPout.RData")
# EXPORT THE DATA
write.csv(SLEEP, file = "SLEEPOut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoSLEEP"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))

PhenoList <-  c(
  'SLEEP_Results.Sleep.Light.Phase.Percent',
  'SLEEP_Results.Sleep.Daily.Percent',
  'SLEEP_Results.Breath.Rate.During.Sleep.Mean',
  'SLEEP_Results.Breath.Rate.During.Sleep.Standard.Deviation',
  'SLEEP_Results.Light.Onset.Wake.Median',
  'SLEEP_Results.Diurnal.Wake.Ratio.Median',
  'SLEEP_Results.Sleep.Bout.Lengths.Mean',
  'SLEEP_Results.Light.Sleep.Bout.Lengths.Mean',
  'SLEEP_Results.Peak.Wake.wrt.Dark.Onset.Median',
  'SLEEP_Results.Light.Sleep.Bout.Lengths.Standard.Deviation'
)

DataSet <- SLEEP

TestDates <- "SLEEP_Date.of.test.New"
TestTime <- "SLEEP_Start.Time"
DataName <- "SLEEP"
Year <- "SLEEP_Year"
Month <- "SLEEP_Month"
DayOfWeek <- "SLEEP_DayOfWeek"
#Hour <- "SLEEP_Hour"
#timeCat <- "SLEEP_timeCat"
DayOfYear <- "SLEEP_DayOfYear"
DayOfMonth <- "SLEEP_DayOfMonth"
StrainGeno <- "SLEEP_StrainGeno"
predBW <- "SLEEP_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "SLEEP_Room.origin",
  "SLEEP_Experimenter",
  "SLEEP_Software",
  "SLEEP_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  predBW
  #DayOfWeek
  #timeCat
)

control <- SLEEP[SLEEP$SLEEP_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)

# ABR --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#	SOME COLUMNS HAVE DUPLICATED DATA FOR ABR
# change the minumum cutoff to 3 in the RCode_functions.R
loadPheno <- function(WorkBookName) {
  setwd(dataDir)
  # PhenoName <- read.csv(WorkBookName, header = TRUE) ## load CSV sheet
  # new data is labeled .csv but is tab delimted
  # PhenoName <- read.table(WorkBookName, sep="\t", header = TRUE) ## was giving errors
  PhenoName <- read.delim(WorkBookName, sep="\t", header = TRUE)

  PhenoName <- as.data.frame(PhenoName)
  #colnames(PhenoName) = make.names(colnames(PhenoName))


  if (WorkBookName == "FLOW") {
    ### FLOW data has columns with + and - (positive and neg) these get substituted with '.' making two columns with same name.  This should be changed
    PhenoName$Name <- gsub("\\-","_neg_",PhenoName$Name)
    PhenoName$Name <- gsub("\\+","_pos_",PhenoName$Name)
    PhenoName$Name <- gsub("\\%","Percent",PhenoName$Name)
  }

  ### 	Remove Duplicates from data
  ###	This is keep the first value and remove the second value
  ### code from http://stackoverflow.com/questions/9944816/unique-on-a-dataframe-with-only-selected-columns
  PhenoName <- PhenoName[!duplicated(PhenoName[,c("OrganismKey", "Mouse.ID", "Date.of.Birth", "JR.", "Sex", "Strain.Name","Genotype","Room.origin", "GenotypeSymbol","Job.ID", "TestCode","Life.Status","Exit.Reason", "DateReceived", "DateComplete","Name")]),]

  #### can also use - I did not use it.
  #PhenoName <- PhenoName[!duplicated(PhenoName[,c(1:16)]),]

  ####	Long to Wide here
  ####	Plyr will calculate summary stats if there are repeated rows - be careful
  PhenoName <- dcast(PhenoName, OrganismKey + Mouse.ID + Date.of.Birth + JR.+ Sex + Strain.Name + Genotype + Room.origin + GenotypeSymbol + Job.ID + TestCode + Life.Status + Exit.Reason + DateReceived + DateComplete ~ Name, value.var = "Value")

  colnames(PhenoName) = make.names(colnames(PhenoName))

  PhenoName$StrainGeno <- paste(PhenoName$Strain.Name, PhenoName$GenotypeSymbol) ### creates new factor that is combination of strain name and genotype symbol - use this for plotting

  ### Remove all animals with Genotype.Symbol is +/+
  remove_list <- which(PhenoName$GenotypeSymbol == "+/+")
  print(paste("Removing ", length(remove_list)," animals that don't have +/+ Genotype.", sep = ""))
  PhenoName <- PhenoName[-remove_list,]
  rm(remove_list)

  ### Change all +/- genotypes to -/+
  change_list <- which(PhenoName$GenotypeSymbol == "+/-")
  print(paste("Changing ", length(change_list)," animals that have +/- to -/+ Genotype.", sep = ""))
  PhenoName$GenotypeSymbol[change_list] <- "-/+"
  rm(change_list)

  ### Remove all animals with Genotype.Symbol is missing (means that genotype failed or is not complete)
  remove_list = intersect(which(PhenoName$GenotypeSymbol == ""), which(PhenoName$Strain.Name != "C57BL/6NJ"))
  print(paste("Removing ", length(remove_list)," animals that don't have any Genotype data.", sep = ""))
  PhenoName <- PhenoName[-remove_list,]
  rm(remove_list)

  ### DOB gets converted to number, check this - this does not happens all systems
  #PhenoName$Date.of.Birth <- as.Date(PhenoName$Date.of.Birth , origin = "1899-12-30")


  #### ADD ES Cell info
  #### HeartWeight and CBC - add ES cell information for modeling
  ### there is a private mutation in B6NJ that causes low heart weight in controls and may effect other phenotypes.
  ### to overcome this use the center or ES cell line as a covariate.
  PhenoName$EScell <- ifelse(grepl("^C57BL/6NJ$", PhenoName$Strain.Name), "C57BL/6NJ",
                             ifelse(grepl("Vlcg", PhenoName$Strain.Name), "Vlcg",
                                    ifelse(grepl("Wtsi", PhenoName$Strain.Name), "Wtsi",
                                           ifelse(grepl("em1J", PhenoName$Strain.Name), "em1J",
                                                  ifelse(grepl("Hmgu", PhenoName$Strain.Name), "Hmgu",
                                                         ifelse(grepl("Mbp", PhenoName$Strain.Name), "Mbp",
                                                                ifelse(grepl("tm1b", PhenoName$Strain.Name), "Misc",
                                                                       ifelse(grepl("tm1", PhenoName$Strain.Name), "Misc",
                                                                              ifelse(grepl("COIN", PhenoName$Strain.Name), "Misc",
                                                                                     ifelse(grepl("Lutzy", PhenoName$Strain.Name), "Misc",
                                                                                            "em1J"))))))))))

  PhenoName$EScell <- as.factor(PhenoName$EScell)

  ### remove StrainGeno < 3
  StrainGeno.Freq <- as.data.frame(table(PhenoName$StrainGeno))
  colnames(StrainGeno.Freq) <- c("StrainGeno","Counts")
  StrainGeno2rm <- as.character(StrainGeno.Freq[which(StrainGeno.Freq$Counts < 3),1])

  print(paste("Removing ", length(StrainGeno2rm)," Strains with n < 3.", sep = ""))


  fname <- gsub(WorkBookName,pattern = "csv",replacement = "StrainGenoFreq.csv")
  write.csv(file=fname, StrainGeno.Freq, quote = FALSE, row.names = FALSE)

  fname <- gsub(WorkBookName,pattern = "csv",replacement = "StrainGenoRemoved.csv")
  StrainGeno.Freq.rm <- StrainGeno.Freq[which(StrainGeno.Freq$Counts < 3),]
  write.csv(file=fname, StrainGeno.Freq.rm, quote = FALSE, row.names = FALSE)
  PhenoName <- PhenoName[-c(which(PhenoName$StrainGeno %in% StrainGeno2rm)),]

  setwd(workDir)
  return(PhenoName)
}
ABR <-	loadPheno(paste0(dataDir ,"ABR.csv")) ## there are duplicates that I fixed manually ::::: One time fix - need to fix database
ABR <- fixTesterNames(ABR, "ABR")
ABR <- convertVars(data = ABR, dict = DictData, "ABR")  # change classes of factors from character to the appropriate
# this is for removing rows from original csv file - need to do this once and then write ABR_2.csv file
# rows to drop
# 					ABR <- ABR[ABR$Name != "Minimum dB:12kHz No response to max stimuli",]
# 					ABR <- ABR[ABR$Name != "Minimum dB:18kHz No response to max stimuli",]
# 					ABR <- ABR[ABR$Name != "Minimum dB:24kHz No response to max stimuli",]
# 					ABR <- ABR[ABR$Name != "Minimum dB:30kHz No response to max stimuli",]
# 					ABR <- ABR[ABR$Name != "Minimum dB:6kHz No response to max stimuli",]
# 					ABR <- ABR[ABR$Name != "Waveform File:PDF report",]
#
# 					write.csv(ABR, "ABR_2.csv")

# http://stackoverflow.com/questions/37736131/lubridate-mdy-function
# lubridate failed when modifyng 2013 and 13 - I posted on Stack Overflow to sort this out
# my_select <-   function(trained){
#   n_fmts <- nchar(gsub("[^%]", "", names(trained))) + grepl("%y", names(trained))*1.5
#   names(trained[ which.max(n_fmts) ])
# }
ABR$Date.of.test.New <- parse_date_time(ABR$Test.date, orders = c("mdy","mdy_hms"), select_formats = my_select)
# How many NA in the Date.of.test.New
print(paste0("The Number of NA in the Date.of.test.New is ",sum(is.na(ABR$Date.of.test.New))))
# Remove all data before 2013 ##
ABR <- ABR[ABR$Date.of.test.New > ymd("2012-12-31"),]
# ABR TIME in multiple formats
# "10:05 AM" | "14:25"
data1 <-  parse_date_time(ABR$Time.of.day, c("%I:%M %p","%H:%M", "hm", "hms")) # lubridate works mostly some have today's date and time, others have "0000-01-01 05:28:58" "2015-10-27 07:05:00"
data1 <- strftime(data1) #convert to string
data1 <- strsplit(data1, " ") #split
data1 <- sapply(data1, "[",2) # take second element from list
data1 <- strptime(data1, "%H:%M:%S") # convert back to time
ABR$Start.time.New <- data1 #add to data frame
rm(data1)
ABR <- ExtractDateParameters(ABR, 'ABR')
ABR <- ExtractTimeParameters(ABR, 'ABR')
# Fix Software Name
ABR$Software[ABR$Software == "HIS Smart "] <- "HIS Smart"
# still data duplications - not easy to upload this data!!!! sent email with data duplication to James to fix
# 2015-10-03 remove all data duplicates manually from ABR.csv
#  Test.date has some elements with date and time together
#  set of dates in numberic format
#  set of dates with time
#  FIXED ALL IN EXCEL!!!!
#  REMOVE OUTLIERS
PhenoList <- 	c(
  'Minimum.dB.6kHz.Minimum.Threshold',
  'Minimum.dB.12kHz.Minimum.Threshold',
  'Minimum.dB.18kHz.Minimum.Threshold',
  'Minimum.dB.24kHz.Minimum.Threshold',
  'Minimum.dB.30kHz.Minimum.Threshold'
)
# EXPORT THE DATA before outlier removal
write.csv(ABR, file = "ABROut0.csv",row.names = FALSE)
## Remove outliers from data
ABR <- remove_outliers_main(komp.data = ABR ,list = PhenoList, stdv = 2, iqr = TRUE)

# Add  BW data
ABR <- addBW(ABR, bw.komp, "ABR", avg_aot = 112)

colnames(ABR)[-c(1:5)] = paste("ABR_", colnames(ABR)[-c(1:5)], sep = "")

save(ABR, file="ABRout.RData")
# EXPORT THE DATA
write.csv(ABR, file = "ABROut.csv",row.names = FALSE)

# Run the RCodePlotsDiagnostics.r and PlotEverything scrips
outDir<- "phenoABR"
out1Dir <- "ControlDiagnotsticPlots"
setwd(file.path(workDir))
dir.create(file.path(workDir, outDir, out1Dir), showWarnings = FALSE, recursive = TRUE) # will give warning if the directory already exists
setwd(file.path(workDir, outDir, out1Dir))


PhenoList <- 	c(
  'ABR_Minimum.dB.6kHz.Minimum.Threshold',
  'ABR_Minimum.dB.12kHz.Minimum.Threshold',
  'ABR_Minimum.dB.18kHz.Minimum.Threshold',
  'ABR_Minimum.dB.24kHz.Minimum.Threshold',
  'ABR_Minimum.dB.30kHz.Minimum.Threshold'
)

DataSet <- ABR

TestDates <- "ABR_Date.of.test.New"
TestTime <- "ABR_Start.time.New"
DataName <- "ABR"
Year <- "ABR_Year"
Month <- "ABR_Month"
DayOfWeek <- "ABR_DayOfWeek"
Hour <- "ABR_Hour"
timeCat <- "ABR_timeCat"
DayOfYear <- "ABR_DayOfYear"
DayOfMonth <- "ABR_DayOfMonth"
StrainGeno <- "ABR_StrainGeno"
predBW <- "ABR_PredictedBW_outs"

ExperimentalFactors <- c(
  "Sex",
  "ABR_Experimenter",
  "ABR_Equipment.ID.",
  #"ABR_Last.calibration.date", # only one calibration date - so I removed it
  "ABR_Room.origin",
  "ABR_EScell"
)

ExperimentalFactors2 <- c(
  ExperimentalFactors,
  Year,
  Month,
  DayOfWeek,
  timeCat
)

control <- ABR[ABR$ABR_StrainGeno == "C57BL/6NJ ",]  # control is now by itself, has an space at end

#Run PlotEverything scrips
#source(plotAllFile)
