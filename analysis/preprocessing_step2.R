library(ggplot2)

batch = "2022-05-13" #to change for different version of data
projectDir = "/projects/compsci/vmp/USERS/heh/KOMP2022/" #project dir

path.data <- paste0(projectDir, "output/", batch, "_out/")
path.util <- paste0(projectDir, "code/")

source(paste0(path.util,"util_10082019.R"))

reformat.meta <- FALSE
save.meta <- FALSE
save.figure <- FALSE

setwd(path.data)
## it is already formated to KOMP_metadata_03132017.csv
## just read the file KOMP_metadata_03132017.csv
meta.data.file <- paste0(projectDir, "data/KOMP_metadata_03132017.csv")
meta.data <- read.table(meta.data.file, sep=",", header=TRUE)
#####################################################################################
## Reformat all XXOut.csv files (e.g. OFAOut.csv)
## 1. make a new data file containing only variables used in the downstream analyses
## 2. make variable names consistent across different domains
#####################################################################################

domain <- meta.data$pheno.group[!(meta.data$pheno.group %in% c("Urianalysis", "EKG", "ERG"))] #17_domains: The data for Urianalysis is not ready yet.

for(i in domain){
  ##i <- domain[-c(2,8)][4]

  ## Use the code below for domain[c(2,8)]
  ## data <- read.table(paste0(path.data, i ,"out.csv"), sep=",", header=TRUE)

  ## Use the code below for other domains except domain[c(2,8)]
  load(paste0(path.data, i, "out.RData"))
  data <- get(paste0(i))

  col.names <- gsub("OFA_|HB_|LD_|TST_|SLEEP_|RR_|ECT_|PPI_|ABR_|GRIP_|Hematology_|GTT_|HeartWeight_|Insulin_|Urianalysis_|EKG_|ClinicalBloodChemistry_|BodyComp_","", colnames(data))
  col.names <- gsub("Experimenter.ID|Experimeter", "Experimenter", col.names)
  col.names <- gsub("Experimenter", "Experimenter.ID", col.names)
  col.names <- gsub("Arena.ID", "Arena", col.names)
  #    col.names <- gsub("Arena", "Arena.ID", col.names)
  col.names <- gsub("timeCat", "Time.categ", col.names)
  col.names <- gsub("GTT.Data.Fields.Body.weight..g.|Weights.Body.Weight", "BW", col.names)
  col.names <- gsub("PredictedBW_outs", "BW", col.names)
  colnames(data) <- col.names

  pheno.list <- as.character(meta.data[meta.data$pheno.group==i,]$pheno.list)
  pheno.list <- unlist(strsplit(pheno.list, " "))
  pheno.list <- pheno.list[pheno.list%in%colnames(data)] ## remove phenotypes not in pheno.data

  predictor.list <- as.character(meta.data[meta.data$pheno.group==i,]$predictor.list)
  predictor.list <- unlist(strsplit(predictor.list, " "))
  predictor.list <- predictor.list[predictor.list%in%colnames(data)] ## remove predictors not in pheno.data
  predictor.list

  colnames(data)

  ##predictor.list%in%colnames(data)
  ##pheno.list%in%colnames(data)

  new.col <- c(pheno.list, "Mouse.ID","Strain.Name", "GenotypeSymbol", predictor.list)
  print(sum(!new.col%in%colnames(data)))

  data <- data[,new.col]
  data$Date.of.test.New <- factor(MakeDateFormatSame(data$Date.of.test.New))
  levels(data$GenotypeSymbol)[levels(data$GenotypeSymbol)==""] <- "+/+"

  if("Arena"%in%predictor.list){
    data$Arena <- factor(data$Arena)
  }
  levels(data$Strain.Name)[levels(data$Strain.Name)=="C57BL/6NJ "] <- "C57BL/6NJ"

  print(i)
  #print(summary(data))
  #print(table(data$Date.of.test.New))

  write.table(summary(data), file=paste0(path.data, i, "_summary_log.csv"), quote=FALSE, sep=",", row.names = FALSE, col.names=TRUE)
  #write.table(xtabs(~ data$Strain.Name + data$Sex), file=paste0(path.data, i, "_count_per_mutantline.csv"), quote=FALSE, sep=",", row.names = TRUE, col.names=TRUE)
  write.table(data, file=paste0(path.data, i, "Out3.csv"), quote=FALSE, sep=",", row.names = FALSE, col.names=TRUE)

  ##################################
  ## Generate exploratory figures
  ##################################
  if(save.figure){
    pdf(file=paste0(path.data.figures,i,"_date_vs_genotype.pdf"), height=25)
    print(ggplot(data, aes(Date.of.test.New)) + geom_bar(aes(fill=GenotypeSymbol)) + coord_flip())
    dev.off()

    pdf(file=paste0(path.data.figures,i,"_mutantline_vs_sex.pdf"), height=100)
    print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Sex)) + coord_flip())
    dev.off()

    pdf(file=paste0(path.data.figures,i,"_mutantline_vs_genotype.pdf"), height=100)
    print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=GenotypeSymbol)) + coord_flip())
    dev.off()

    pdf(file=paste0(path.data.figures,i,"_mutantline_vs_date.pdf"), height=100, width=30)
    print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Date.of.test.New), colour="black") + coord_flip())
    dev.off()

    if("Room.origin"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_room.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Room.origin)) + coord_flip())
      dev.off()
    }

    if("Experimenter.ID"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_experimenter.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Experimenter.ID)) + coord_flip())
      dev.off()
    }

    if("Time.categ"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_timecateg.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Time.categ)) + coord_flip())
      dev.off()
    }

    if("DayOfWeek"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_dayofweek.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=DayOfWeek)) + coord_flip())
      dev.off()
    }

    if("Month"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_month.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Month)) + coord_flip())
      dev.off()
    }

    if("Arena"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_arena.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(Strain.Name)) + geom_bar(aes(fill=Arena)) + coord_flip())
      dev.off()
    }

    if("BW"%in%predictor.list){
      pdf(file=paste0(path.data.figures,i,"_mutantline_vs_bw_by_sex.pdf"), height=100)
      print(ggplot(subset(data, Strain.Name!="C57BL/6NJ"), aes(x=Strain.Name, y=BW)) + geom_boxplot(aes(fill=Sex)) + coord_flip())
      dev.off()
    }

    pdf(file=paste0(path.data.figures,i,"_date_scatter_by_genotype.pdf"), width=30)
    for(j in pheno.list){
      print(ggplot(data, aes_string(x="Date.of.test.New", y=j)) + geom_point(aes(colour=factor(GenotypeSymbol)))+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    }
    dev.off()

    pdf(file=paste0(path.data.figures,i,"_date_scatter_by_sex.pdf"), width=30)
    for(j in pheno.list){
      print(ggplot(data, aes_string(x="Date.of.test.New", y=j)) + geom_point(aes(colour=Sex))+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    }
    dev.off()

    sink(paste0(path.data.figures,"data_summary_log.txt"), split=TRUE, append=TRUE)
    cat("Phenotype Domain:", i, "\n")
    cat("# of control animals:", nrow(subset(data, Strain.Name=="C57BL/6NJ")), "\n")
    cat("# of mutant animals:", nrow(subset(data, Strain.Name!="C57BL/6NJ")), "\n")
    cat("# of mutant lines:", length(levels(data$Strain.Name))-1, "\n")
    cat("avg # of mutant animals per line:", mean(table(subset(data, Strain.Name!="C57BL/6NJ")$Strain.Name)), "\n")
    cat("\n\n\n")
    sink()
  }
}

