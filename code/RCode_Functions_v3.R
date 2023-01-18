#packages
library(ggplot2)
library(readxl)
library(stringr)
library(plyr)
library(reshape2)
library(lubridate)
library(gdata)

#' For the input KOMP_DataDictionary.xlsx.
#' DATA DICTIONARY
#' OTENTIAL USE OF convertVars
#' @param KOMP_DataDictionary.xlsx - a dataframe for the domain
#' @return - DictData
#'
projectDir = "/projects/compsci/vmp/USERS/heh/KOMP2022/" #project dir

DictData <- read_excel(paste0(projectDir, "data/KOMP_DataDictionary.xlsx"), sheet = "DataDict") # read data dictionary
DictData$Parameter.fixed <- gsub("-|>|:|\\(|\\)| |\\#|%|\\/",".",DictData$Parameter) ## since we have not done MakeNames command yet, the future column names are not properly formated.  This command will fixt this.

#' For the input dataframe dt and vector list (phenotypes in dt), detect outliers and remove outliers as NA.
#'
#' @param dt - a dataframe for the domain
#' @param list - a vector for the phenotypes in the domain
#' @param stdv - Length as multiple of IQR/SD.
#' @param na.rm - a logical value indicating whether NA values should be stripped before the computation proceeds. Default: TRUE
#' @param iqr - a logical value indicating whether IQR should be used for outlier detection. If false, use SD.
#' @return - a dataframe dt with outlier removal for the domain
#'
remove_outliers <- function(dt,list, stdv,  na.rm = TRUE, iqr = TRUE, ...) {
  col.ids <- which(colnames(dt) %in% list)
  for(i in 1:length(col.ids)) {
    x <- dt[,col.ids[i]]
    #print(x)
    if (is.na(sd(x,na.rm = na.rm))) {
      #print("skip")
      next
    } else {
      qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
      iqr.val <- IQR(x, na.rm = na.rm)
      H <- ifelse(iqr, stdv * iqr.val, stdv * sd(x, na.rm = na.rm))
    }
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    dt[,col.ids[i]] <- y
  }
  dt
}

#' For the input dataframe komp.data and list (phenotypes in dt), detect outliers and remove outliers as NA.
#'
#' @param komp.data - a dataframe for the domain
#' @param list - a vector for the phenotypes in the domain
#' @param stdv - Length as multiple of IQR/SD.
#' @param iqr - a logical value indicating whether IQR should be used for outlier detection. If false, use SD.
#' @return - a list object which contatin a gxl estimate in estimates_for_gxl and a ANOVA table in anova.str
#'
remove_outliers_main <- function(komp.data,list,stdv,iqr = TRUE) {
  StrainColumn <- grepl("Strain.Name",colnames(komp.data))
  strain.list <- unique(komp.data[,StrainColumn])
  strain.list <- na.omit(strain.list)

  DateColumn <- grepl("Date.of.test.New",colnames(komp.data))

  ###working on the mutant group###
  strain.list.mutant <- strain.list[strain.list != "C57BL/6NJ"]
  dt.mutant <- NULL
  for (i in 1:length(strain.list.mutant)) {
    cat("Working on Strain = ", as.character(strain.list.mutant[i]),"\n")
    dt.sub <- komp.data[which(komp.data[,StrainColumn] == strain.list.mutant[i]),]
    if(dim(dt.sub)[1] >= 5) {		### this is key - will only work on genotypes with 5 or more individuals
      dt.sub.no.out <- remove_outliers(dt.sub, list=list, stdv=stdv, iqr)
      dt.mutant <- rbind.data.frame(dt.mutant,dt.sub.no.out)
    }
    else{
      dt.mutant <- rbind.data.frame(dt.mutant,dt.sub)
    }
  }

  ###working on the control group with NA and not NA date test ###
  dt.control <- komp.data[which(komp.data[,StrainColumn] == "C57BL/6NJ"),]
  dt.control <- remove_outliers(dt.control, list=list, stdv=3, iqr)
  ###extract rows with date test having NA.
  dt.control.date.NA <- dt.control[is.na(dt.control[,DateColumn]),]
  ###extract rows with date test not having NA.
  dt.control <- dt.control[!is.na(dt.control[,DateColumn]),]

  ctrl.test.date <- unique(as.character(dt.control[,DateColumn]))
  ctrl.test.date <- ctrl.test.date[!is.na(ctrl.test.date)]
  dt.final.control <- NULL
  for (k in 1:length(ctrl.test.date)){
    cat("Working on date test = ", ctrl.test.date[[k]],"\n")
    dt.sub2 <- dt.control[as.character(dt.control[,DateColumn]) == ctrl.test.date[[k]],]
    if(dim(dt.sub2)[1] >= 5) {		### this is key - will only work on controls with 5 or more individuals
      dt.sub.no.out <- remove_outliers(dt.sub2,list=list, stdv=stdv, iqr)
      dt.final.control <- rbind.data.frame(dt.final.control,dt.sub.no.out)
      col.ids <- which(colnames(dt.sub.no.out) %in% list)
      for(m in 1:length(col.ids)) {
        temp <- which(is.na(dt.sub.no.out[,col.ids[m]]))
        n.end <- length(temp)
        cat(n.end, " outliers Removed for phenotype ", colnames(dt.sub.no.out)[col.ids[m]], "\n")
      }
    }
    else{
      cat("n < 5. No Outliers Removed\n")
      dt.final.control <- rbind.data.frame(dt.final.control,dt.sub2)
    }
  }

  ###combine mutant, control group without NA date test and with NA date test
  dt.final <- rbind.data.frame(dt.final.control,dt.control.date.NA, dt.mutant)
  colnames(dt.final) <- colnames(komp.data)
  dt.final[order(dt.final$OrganismKey),]
  dt.final
}

#' For the input WorkBookName,.
#'
#' @param WorkBookName - a csv file for the domain to load
#' @return - a dataframe objet
#'
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

  ### remove StrainGeno < 8
  StrainGeno.Freq <- as.data.frame(table(PhenoName$StrainGeno))
  colnames(StrainGeno.Freq) <- c("StrainGeno","Counts")
  StrainGeno2rm <- as.character(StrainGeno.Freq[which(StrainGeno.Freq$Counts < 8),1])

  print(paste("Removing ", length(StrainGeno2rm)," Strains with n < 8.", sep = ""))


  fname <- gsub(WorkBookName,pattern = "csv",replacement = "StrainGenoFreq.csv")
  write.csv(file=fname, StrainGeno.Freq, quote = FALSE, row.names = FALSE)

  fname <- gsub(WorkBookName,pattern = "csv",replacement = "StrainGenoRemoved.csv")
  StrainGeno.Freq.rm <- StrainGeno.Freq[which(StrainGeno.Freq$Counts < 8),]
  write.csv(file=fname, StrainGeno.Freq.rm, quote = FALSE, row.names = FALSE)
  PhenoName <- PhenoName[-c(which(PhenoName$StrainGeno %in% StrainGeno2rm)),]

  setwd(workDir)
  return(PhenoName)
}

#' For the input DataSet, ExtractDateParameters.
#'
#' @param DataSet - a dataframe for the domain
#' @param DataName - a character name for the domain
#' @return - a dataframe object
#'
ExtractDateParameters <- function (DataSet, DataName) {
  DateColumnNumber <- which(colnames(DataSet)=="Date.of.test.New")
  DateColumn <- DataSet[,DateColumnNumber]
  Year <- as.factor(year(DateColumn))
  Month <- month(DateColumn, label = TRUE)
  DayOfWeek <- wday(DateColumn, label = TRUE)
  DayOfYear <- yday(DateColumn)
  DayOfMonth <- mday(DateColumn)
  DataSet <- cbind(DataSet, Year, Month, DayOfWeek, DayOfYear, DayOfMonth)
  return(DataSet)
}

#' For the input DataSet, ExtractTimeParameters
#'
#' @param DataSet - a dataframe for the domain
#' @param DataName - a character name for the domain
#' @return - a dataframe object
#'
ExtractTimeParameters <-	function (DataSet, DataName) {
  TimeColumnNumber <- which(colnames(DataSet)== "Start.time.New")
  TimeColumn <- DataSet[,TimeColumnNumber]
  Hour <- hour(TimeColumn)
  timeCat = ifelse(Hour<9.00, "EarlyMorning", ifelse(Hour >=9.00 & Hour <12, "Morning", ifelse(Hour >=12 & Hour <15, "Afternoon", "LateAfternoon")))
  ## before 9 is EarlyMorning,
  ## 9 to noon is Morning,
  ## noon to 3 is afternoon
  ## 3 to midnight is late afternoon
  DataSet <- cbind(DataSet, Hour, timeCat)
  return(DataSet)
}

#' For the input DataSet, ExtractTimeParameters
#'
#' @param trained
#' @param drop
#' @return
#'
my_select <-   function(trained,drop){
  n_fmts <- nchar(gsub("[^%]", "", names(trained))) + grepl("%y", names(trained))*1.5
  names(trained[ which.max(n_fmts) ])
}

#' For the input date.list, MakeSameDateFormat
#'
#' @param date.list - a date.list
#' @return - a date.list formatted
#'
MakeSameDateFormat <-function(date.list){
  date.list <- parse_date_time(x = date.list, orders=c("mdy", "m/d/y", "Ymd", "Ymd H:M:S", "Y/m/d", "m/d/Y",  "mdy H:M","mdY"), select_formats = my_select)
  date.list
}

#' For the input time.list
#'
#' @param time.list - a time.list
#' @return - a time.list
#'
MakeSameTimeFormat <-function(time.list){
  time.list <- parse_date_time(x = time.list, orders=c("H:M:S", "H:M"))
  time.list
}

#' For the Split Dates and Times
#'
#' @param data.komp - a data.komp
#' @param dataColumn - a dataColumn
#' @return - a time.list
#'
extractDates <- function (data.komp, dataColumn) {
  data1 <- do.call(rbind, strsplit(data.komp$dataColumn," "))
  data1[1] <- mdy(data1[1])
  data1[2] <- strptime(data1[2], "%H:%M")
  colnames(data1) <- c(paste(data.komp, ".",dataColumn, ".Date.New", sep = ""), paste(data.komp, ".",dataColumn, ".Time.New", sep = ""))
  data.komp <- cbind(data.komp, data1)
  data.komp
}

#' For  TO CONVERT CHARACTER FACTORS TO NUMBERS, FACTORS OR DATE,
#' WRITTEN BY VIVEK PHILLIP
#' RELIES ON A DATA DICTIONARY, THAT HAS INFORMATION ABOUT WHAT THE FACTORS SHOULD BE
#' Dictionary_KOMP.csv' has 4 columns samples below
#' Phenotype	Parameter					Field		Sample
#' HB			Parameters:Total Hole Pokes	Numeric		37
#'
#' @param data - a data.komp
#' @param dict - a dataColumn
#' @param pheno - a pheno
#' @return - a data
#'
convertVars <- function(data,dict,pheno) {
  dict.sub <- dict[which(dict$Phenotype == pheno | dict$Phenotype == "ALL"),]  ### this will select OFA OR ALL, important for common factors
  data.colnames <- colnames(data)
  for(i in 1:length(data.colnames)) {
    idx <- which(dict.sub$Parameter.fixed == data.colnames[i])
    if (length(idx) > 0) {
      if(dict.sub[idx,3] == "Numeric") {
        data[,i] <- as.numeric(data[,i])
      }
      if(dict.sub[idx,3] == "Factor") {
        data[,i] <- as.factor(data[,i])
      }
      if(dict.sub[idx,3] == "Date_ymd_hms") {
        data[,i] <- ymd_hms(data[,i])
      }
      if(dict.sub[idx,3] == "Date_ymd") {
        data[,i] <- ymd(data[,i])
      }
      if(dict.sub[idx,3] == "Date_mdy") {
        data[,i] <- mdy(data[,i])
      }
      if(dict.sub[idx,3] == "Date_mdy_hm") {
        data[,i] <- mdy_hm(data[,i])
      }
      if(dict.sub[idx,3] == "Time_hms") {
        data[,i] <- hms(data[,i])
      }
    }
  }

  data
}

#' FIX TESTER NAMES FUNCTION
#'
#' @param dataToFix - a dataToFix
#' @param dataName - a dataName
#' @return - a data
#'
fixTesterNames <- function(dataToFix, dataName) {

  ###### Fix Tester names
  # OFA LD
  if (dataName == "OFA"| dataName == "LD" | dataName == "GTT" | dataName == "BodyComp") {
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "Zachery Seavey "] <- "Zachery Seavey"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "ZS"] <- "Zachery Seavey"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "REP"] <- "Rose Presby"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "RP"] <- "Rose Presby"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "Accuscan"] <- "Unknown"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "OH"] <- "Olivia Hon"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "JC"] <- "James Clark"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "WR"] <- "Wilson Roper"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "PF"] <- "Pamelia Fraungruber"
  }
  # HB PPI ABR GRIP HeartWeight ECT SLEEP
  if (dataName == "RR" | dataName == "HB"| dataName == "PPI" | dataName == "ABR" | dataName == "GRIP" | dataName == "HeartWeight" | dataName == "ECT" | dataName == "SHIRPA" | dataName == "SLEEP") {
    dataToFix$Experimenter[dataToFix$Experimenter == "ZS"] <- "Zachery Seavey"
    dataToFix$Experimenter[dataToFix$Experimenter == "GW"] <- "Gaylynn Wells"
    dataToFix$Experimenter[dataToFix$Experimenter == "SB"] <- "Samantha Burrill"
    dataToFix$Experimenter[dataToFix$Experimenter == "RP"] <- "Rose Presby"
    dataToFix$Experimenter[dataToFix$Experimenter == "KL"] <- "Kate Lachapelle"
    dataToFix$Experimenter[dataToFix$Experimenter == "CR"] <- "Christine Rosales"
    dataToFix$Experimenter[dataToFix$Experimenter == "Christine Rosales "] <- "Christine Rosales"
    dataToFix$Experimenter[dataToFix$Experimenter == "NC"] <- "Neil Cole"
    dataToFix$Experimenter[dataToFix$Experimenter == "NB"] <- "Nicholas Brown"
    dataToFix$Experimenter[dataToFix$Experimenter == "KL"] <- "Kate Lachapelle"
    dataToFix$Experimenter[dataToFix$Experimenter == "jr"] <- "Jennifer Ryan"
    dataToFix$Experimenter[dataToFix$Experimenter == "JR"] <- "Jennifer Ryan"
    dataToFix$Experimenter[dataToFix$Experimenter == "LH"] <- "Leslie Haynes"
    dataToFix$Experimenter[dataToFix$Experimenter == "LHl"] <- "Leslie Haynes"
    dataToFix$Experimenter[dataToFix$Experimenter == "OH"] <- "Olivia Hon"
    dataToFix$Experimenter[dataToFix$Experimenter == "SB/NC"] <- "Cole | Burrill"
    dataToFix$Experimenter[dataToFix$Experimenter == "TW"] <- "Troy Wilcox"
    dataToFix$Experimenter[dataToFix$Experimenter == "twilcox"] <- "Troy Wilcox"
    dataToFix$Experimenter[dataToFix$Experimenter == "lcanders"] <- "Laura Anderson"
    dataToFix$Experimenter[dataToFix$Experimenter == "LRB"] <- "Lindsay Bates"
    dataToFix$Experimenter[dataToFix$Experimenter == "LRB"] <- "Lindsay Bates"
  }
  # TST
  if (dataName == "TST") {
    dataToFix$Experimenter[dataToFix$Experimenter == "ZS"] <- "Zachery Seavey"
    dataToFix$Experimenter[dataToFix$Experimenter == "TW"] <- "Troy Wilcox"
  }

  # Insulin
  if (dataName == "Insulin") {
    dataToFix$Blood.Collection.Experimenter[dataToFix$Blood.Collection.Experimenter == "OH"] <- "Olivia Hon"
    dataToFix$Blood.Collection.Experimenter[dataToFix$Blood.Collection.Experimenter == "SC"] <- "Steve Ciciotte"
    dataToFix$Blood.Collection.Experimenter[dataToFix$Blood.Collection.Experimenter == "PF"] <- "Pamelia Fraungruber"
    dataToFix$Blood.Analysis.Experimenter.ID[dataToFix$Blood.Analysis.Experimenter.ID == "SC"] <- "Steve Ciciotte"
  }
  # Urianalysis
  if (dataName == "Urianalysis" | dataName == "EKG" ) {
    dataToFix$Experimenter.ID.[dataToFix$Experimenter.ID. == "OH"] <- "Olivia Hon"
    dataToFix$Experimenter.ID.[dataToFix$Experimenter.ID. == "PF"] <- "Pamelia Fraungruber"
  }

  # ClinicalBloodChemistry

  if (dataName == "ClinicalBloodChemistry") {
    dataToFix$Analyst.ID[dataToFix$Analyst.ID == "SC"] <- "Steve Ciciotte"
    dataToFix$Experimenter.ID[dataToFix$Experimenter.ID == "OH"] <- "Olivia Hon"
  }


  return(dataToFix)

}

#' predict body weight
#' Predict BW
#' updated version v4
#' changed "Date.of.test" to "Date.of.test.New" on line 15
#' Jan 19: Updated the bw to correctly work the the date formats.
#' Jan 19: Also added a progress status, which shows %done
#' Jan 27: Corrected the code to select the closest body weight. Lines 75 - 76 contains the code that fixes the bug
#' Jan 31: Added the folllowing code to remove outliers from predicetedBW
#' predbwMean <- mean(predbw,na.rm = TRUE)
#' predbwSD <- sd(predbw,na.rm = TRUE)
#' uclpredbw <- predbwMean + 2*predbwSD
#' lclpredbw <- predbwMean - 2*predbwSD
#' data$PredictedBW[data$predbw < lclpredbw | data$predbw > uclpredbw] <- NA
#' PredictedBW now contains NA for outliers
#' @param mouse.id - a mouse.id
#' @param predict.age - a predict.age
#' @param avg_aot - a avg_aot
#' @return - predicted body weight
#'
predic.bw <- function(mouse.id,predict.age=16,avg_aot) {
  data=bw.komp
  data.MID <- data[which(data$Mouse.ID == mouse.id),]
  poly.fit <- lm(BW~age+I(age^2)+I(age^3),data=data.MID)
  pred <- predict(poly.fit,newdata=data.frame(age=predict.age))
  names(pred) <- predict.age
  return(pred)
}

#' add body weight
#' Added a new parameter to the fucntion sd_thresh. This is for determining the numbers of sd to be considered as an outlier
#' avg_aot is the average age at which the testing was done so for OFA at 9 weeks its 63days.
#' @param data - data
#' @param bw.data - bw.data
#' @param prefix - prefix
#' @param sd_thresh - a sd_thresh
#' @param avg_aot - avg_aot
#' @return - dataframe with body weight added
#'
addBW <- function(data,bw.data,prefix,sd_thresh=4,avg_aot) { ### SD Threshold set to 4, avg_aot units are days
  mid.col <- which(colnames(data) == "Mouse.ID")
  dob.col <- which(colnames(data) == "Date.of.Birth")
  dot.col <- which(colnames(data) == "Date.of.test.New")
  zeros <- rep(0,dim(data)[1])
  closestbw <- NULL
  data$closestbw <- zeros
  predbw <- NULL
  data$predbw <- zeros
  n <- dim(data)[1]
  cat("Starting\n")
  #####Begin pred BW Section
  for(i in 1:dim(data)[1]) {

    if (i %% 1000 == 0) {
      pct_done = (100 * i / n)
      pct_done = floor(pct_done)
      pct_done = paste(pct_done,"%",sep="")
      cat(pct_done," done\n")
    }

    mid <- as.character(data[i,mid.col])
    dob <- data[i,dob.col]
    dob <- as.Date(strsplit(as.character(dob),' ')[[1]][1])

    dot <- data[i,dot.col]

    if(is.na(dot) || is.na(dob))
    {
      aot <- avg_aot
    }else{
      dot <- as.Date(strsplit(as.character(dot),' ')[[1]][1])
      aot <- as.numeric(dot - dob)
    }
    data.MID <- bw.data[which(bw.data$Mouse.ID == mid),]
    if(dim(data.MID)[1] < 10) {			#minimum N for a line
      predbw <- append(predbw,"NA")
    }
    else {

      poly.fit <- lm(BW~age+I(age^2)+I(age^3),data=data.MID)
      pred <- predict(poly.fit,newdata=data.frame(age=aot))
      if(pred > 75) {
        print(dim(data.MID))
        sp <- smooth.spline(data.MID$age,data.MID$BW)
        pred <- predict(sp,aot)$y
      }
      predbw <- append(predbw,pred)
    }

    data.MID$cia <- data.MID$age
    data.MID$cia_abs <- data.MID$age
    if(dim(data.MID)[1] == 0) {
      closestbw <- append(closestbw,"NA")
    }
    else {
      data.MID$cia <- (data.MID$age - aot)
      data.MID$cia_abs <- abs(data.MID$age - aot)
      min.cia <- which(data.MID$cia_abs == min(data.MID$cia_abs))
      if(length(min.cia) > 1) {
        closest <- mean(data.MID$BW[min.cia]) ### If there are two mice with same value for age, then average them.
      }else {
        closest <- data.MID$BW[min.cia]
      }
      closestbw <- append(closestbw,closest)
    }
    #cat("SrNo =",i,"Done Mouse ID =",mid,"\n")
  }
  data$predbw <- predbw
  data$predbw[which(data$predbw < 0)] <- "NA"
  data$closestbw <- closestbw
  cat("100% Done")
  #colNamepred <- paste(prefix,"predbw",sep="_") # this was originally designed my VP to say "OFA_predbw", my code adds OFA prefix, I'm removing this here
  #colNameclosest <- paste(prefix,"closestbw",sep="_")  # this was originally designed my VP to say "OFA_predbw", my code adds OFA prefix, I'm removing this here
  colNamepred <- ("PredictedBW")
  colNameclosest <- ("ClosestBW")

  colnames(data)[grep(colnames(data),pattern = "predbw")] <- colNamepred
  colnames(data)[grep(colnames(data),pattern = "closestbw")] <- colNameclosest

  data$PredictedBW <- as.numeric(data$PredictedBW)
  predbwMean <- mean(data$PredictedBW,na.rm = TRUE)
  predbwSD <- sd(data$PredictedBW,na.rm = TRUE)
  uclpredbw <- predbwMean + sd_thresh*predbwSD
  lclpredbw <- predbwMean - sd_thresh*predbwSD

  data$PredictedBW_outs <- as.numeric(data$PredictedBW) ## PredictedBW_outs is being used for modeling

  ids_lcl <- which(data$PredictedBW < lclpredbw)
  ids_ucl <- which(data$PredictedBW > uclpredbw)

  if (length(ids_lcl) >0  && length(ids_ucl) >0) {
    ids2rm <- rbind(ids_lcl,ids_ucl)
    data$PredictedBW_outs[ids2rm] <- NA
  }
  if (length(ids_lcl) >0  && length(ids_ucl) == 0) {
    ids2rm <- ids_lcl
    data$PredictedBW_outs[ids2rm] <- NA
  }
  if (length(ids_lcl) == 0  && length(ids_ucl) > 0) {
    ids2rm <- ids_ucl
    data$PredictedBW_outs[ids2rm] <- NA
  }
  data$PredictedBW_outs[is.na(data$PredictedBW_outs)] <- data$ClosestBW[is.na(data$PredictedBW_outs)]
  return(data)
}

#' plotPhenoFunction
#'
#' @param dirName - dirName
#'
plotPhenoFunction <- function(dirName) {
  directory = paste(workDir, dirName, "RCodePlotsDiagnostics.r", sep = "/")
  source(directory, chdir = TRUE)
  source(plotAllFile, chdir = TRUE)
}


#' plotPhenoFunction
#' for phenotypes with no time of test like ECT
#' @param dirName - dirName
#'
plotPhenoFunction2 <- function(dirName) {
  directory = paste(workDir, dirName, "RCodePlotsDiagnostics.r", sep = "/")
  source(directory, chdir = TRUE)
  source(plotAllFile2, chdir = TRUE)
}
