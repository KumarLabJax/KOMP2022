rm(list=ls())

options(warn=1)
args <- commandArgs()
args

str.look <- "-phenogroup="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
pg <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
pg

str.look <- "-modeltype="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
model.type <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
model.type

str.look <- "-filename="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
filename.str <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
filename.str

str.look <- "-wingsize="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
wing.size <- as.integer(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
wing.size

str.look <- "-verbose="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
verbose <- as.integer(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
verbose

if(FALSE){
    rm(list=ls())
    pg <- "OFA"
    ##model.type <- "IMPC"
    ##model.type <- "REGRES"
    model.type <- "LASSO"
    var.equal <- FALSE
    filename.str <- "test"
    wing.size <- 1
    verbose <- FALSE
}

library(lubridate)
library(glmnet) #ridge/lasso/elastic net regression
library(nlme)
##library(lme4)

batch = "2022-05-13" #to change for different version of data

path.wd <- "/projects/compsci/vmp/USERS/heh/KOMP2022/"#project dir
path.meta <- paste0(path.wd,"data/")
path.data <- paste0(path.wd, "output/", batch, "_out/")
path.out <- paste0(path.wd, "output/", batch, "_AT/")
if (!dir.exists(path.out))
{ dir.create(path.out) } #create folder

source(paste0(path.wd, "code/Utils_07182018.R"))
source(paste0(path.wd, "code/AT_Utils_07182018.R"))

filename.str <- paste0(model.type, "_", filename.str)
filename.str

min.num.mutant <- 7  # minimum number of mutant animals per test
min.num.control <- 7 # minimum number of control animals per test

#################
## Read Meta Data
#################
meta.data.file <- paste0(path.meta, "KOMP2018_Meta_Table_v1.2.csv")
meta.data <- read.csv(meta.data.file, header=TRUE)
meta.data <- meta.data[meta.data$Domain != "EKG",]

##Extract phenotype group list
pheno.group.list <- unique(meta.data$Domain)

cat(paste(pg," Started\n"))
cat("\n")

##Load phenotype data
pheno.data <- read.table(paste0(path.data, pg, "Out3.csv"), sep=",", header=TRUE, quote=NULL, comment='',fill=TRUE)

if(verbose){
    cat("\n### Summary of pheno.data ###\n")
    print(summary(pheno.data))
}

##Extract phenotype list
pheno.list <- as.character(subset(meta.data, Domain==pg)$Phenotype)
pheno.list <- pheno.list[pheno.list%in%colnames(pheno.data)] ## remove phenotypes not in pheno.data

##Extract predictor list
pred.list <- as.character(subset(meta.data, Domain==pg)$Predictor_List[1])
pred.list <- unlist(strsplit(pred.list, " "))
pred.list <- replace(pred.list, pred.list=="StrainGeno", "GenotypeSymbol") ## replace "StrainGeno" with "GenotypeSymbol"
pred.list <- pred.list[pred.list%in%colnames(pheno.data)] ## remove predictors not in pheno.data

if(verbose){
    cat("\n### pred.list ###\n")
    print(pred.list)
    #lapply(pheno.data[,pheno.list],class)
}

##-------------------------------------
## Make pheno.data and pred.list ready
##-------------------------------------
for(p in pheno.list){
    pheno.data[,p] <- as.numeric(pheno.data[,p])
    ##pheno.data[,p] <- scale(as.numeric(pheno.data[,p]))
}
if(pg=="RR"){
    pheno.data$Chamber.ID <- as.factor(pheno.data$Chamber.ID)
    pred.list <- pred.list[!pred.list%in%"Arena"] ## Arena has too many NAs ~60%
}
if(pg=="EKG"){
    pheno.list <- pheno.list[!pheno.list%in%"EKG.pNN50...6ms."] ## too low variance
}
if(pg=="LD"){
    pheno.list <- pheno.list[!pheno.list%in%"Collected.Values.Fecal.Boli"] ## too low variance
}
if("Time.categ"%in%pred.list){
    pheno.data$Time.categ <- factor(pheno.data$Time.categ, levels=c("EarlyMorning","Morning","Afternoon","LateAfternoon"))
    levels(pheno.data$Time.categ) <- c("Morning","Morning","Afternoon","Afternoon")
}
if("Month"%in%pred.list){
    pheno.data$Month <- factor(pheno.data$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul",
                                                          "Aug","Sep","Oct","Nov","Dec"))
    levels(pheno.data$Month) <- c("Winter","Winter","Spring","Spring","Spring","Summer","Summer","Summer","Fall","Fall","Fall", "Winter")
}
if("Arena"%in%pred.list){
    pheno.data$Arena <- factor(pheno.data$Arena)
}


##---------------------
## Extract Control data
##---------------------
all.control <- pheno.data[pheno.data$Strain.Name == "C57BL/6NJ",]
all.control$Date.of.test.New <- as.character(all.control$Date.of.test.New)

##----------------------
##  Generate Mutant List
##----------------------
mutant.list <- unique(pheno.data[pheno.data$Strain.Name !="C57BL/6NJ",]$Strain.Name)

## For output
all.out.df <- NULL
all.pred.ind.df <- NULL
for(ml in mutant.list){
    ## ml <- mutant.list[1]
    cat("\n\n\n####################################### ")
    cat(paste("[",pg,"]",ml," Started\n"))
    cat(paste0(ml," (",match(ml,mutant.list),"/",length(mutant.list),")","\n"))
    cat("\n")
    cat("\n")

    mutant <- subset(pheno.data, Strain.Name == ml)
    mutant$Date.of.test.New <- factor(mutant$Date.of.test.New)
    mutant.date.list <- levels(mutant$Date.of.test.New)

    ## Determine control animals based on new control selection rule
    if(wing.size!=0){
        control.date.list <- NULL
        for(ii in mutant.date.list){
            ##ii <- mutant.date.list[1]
            start.d <- as.Date(ii) %m-% months(1)
            end.d <- as.Date(ii) %m+% months(1)
            control.dates <- unique(subset(all.control, Date.of.test.New >= start.d &
                                                    Date.of.test.New <= end.d)$Date.of.test.New)
            control.date.list <- c(control.date.list, control.dates)
        }
        control.date.list <- unique(control.date.list)
        control <- subset(all.control, Date.of.test.New %in% control.date.list)
    } else { ## Use only control animals phenotyped on mutant animals' test dates
        control <- subset(all.control, Date.of.test.New %in% mutant.date.list)
    }
    control$Date.of.test.New <- factor(control$Date.of.test.New)

    ## Convert Genotype Symbols (+/+, -/+, -/-) to Integer (0, 1, 2)
    genosym.case <- factor(mutant$GenotypeSymbol)
    levels(genosym.case) <- list("0"="+/+","1"="-/+","2"="-/-")
    mutant$GenotypeSymbol <- as.integer(levels(genosym.case))[genosym.case]
    genosym.control <- factor(control$GenotypeSymbol)
    levels(genosym.control) <- list("0"="+/+","1"="-/+","2"="-/-")
    control$GenotypeSymbol <- as.integer(levels(genosym.control))[genosym.control]


    num.mutant <- nrow(mutant)
    num.control <- nrow(control)


    if(!(num.mutant >= min.num.mutant&
         num.control >= min.num.control)){
        cat("\n##################\n")
        cat(paste0(ml," skipped \n"))
        cat("##################\n")
        next
    }

    for(pl in pheno.list){
        ##pl <- pheno.list[1]
        if(verbose){
            cat("\n\n\n##################### ")
            cat(paste0(pl," Started!!\n"))
            cat(paste0(pl," (",match(pl, pheno.list),"/",length(pheno.list),")","\n"))
            cat("\n")
            cat(paste("\nCASE: num of NAs in", pl, ":", sum(is.na(mutant[,pl])),"/", nrow(mutant),"\n"))
            cat(paste("CONTROL: num of NAs in", pl, ":", sum(is.na(control[,pl])),"/", nrow(control),"\n"))
            cat(paste("ALL: num of NAs in", pl, ":", sum(is.na(mutant[,pl]))+sum(is.na(control[,pl]))
                     ,"/", nrow(mutant)+nrow(control),"\n"))

            cat("\n### Print pred.list ###\n")
            cat(pred.list,"\n")
        }

        ##Generate test.data, remove samples with NA and drop unused levels.
        test.data <- rbind(mutant, control)
        #pred.list <- pred.list[pred.list != names(which((sapply(test.data, function(x)all(is.na(x))))))]
        #test.data <- test.data[, colnames(test.data) != names(which((sapply(test.data, function(x)all(is.na(x))))))]
        test.data <- MakeTestData(test.data, pl, pred.list, verbose=verbose)

        num.case <- nrow(subset(test.data, GenotypeSymbol!=0))
        num.control <- nrow(subset(test.data, GenotypeSymbol==0))

        num.male.case <- nrow(subset(test.data, GenotypeSymbol!=0&Sex=="Male"))
        num.female.case <- nrow(subset(test.data, GenotypeSymbol!=0&Sex=="Female"))

        num.male.control <- nrow(subset(test.data, GenotypeSymbol==0&Sex=="Male"))
        num.female.control <- nrow(subset(test.data, GenotypeSymbol==0&Sex=="Female"))

        if(verbose){
            cat("\nNum of cases: ",num.case,"\n")
            cat("Num of controls: ",num.control,"\n")
        }

        if(num.case >= min.num.mutant && num.control >= min.num.control){
            if(model.type == "LASSO"){
                out <- RunLASSO.AT(test.data, pl, pred.list, verbose=verbose)
                stat.vec <- out[[1]]
                var.equal <- out[[2]]
                error.code <- out[[3]]
                pred.ind <- out[[4]]
                surv.env.pred.list <- out[[5]]
                out.df <- data.frame(domain=pg, pheno=pl, mutant=ml,
                                     num.case=num.case, num.control=num.control,
                                     num.male.case=num.male.case, num.female.case=num.female.case,
                                     num.male.control=num.male.control, num.female.control=num.female.control,
                                     geno.std.beta=stat.vec[1],
                                     geno.beta=stat.vec[2], geno.stde=stat.vec[3],
                                     geno.tstat=stat.vec[4], geno.pval=stat.vec[5],
                                     var.equal=var.equal, error.code=error.code,
                                     surv.env.pred.list)
                pred.ind.df <- cbind(data.frame(domain=pg, pheno=pl, mutant=ml), data.frame(t(pred.ind)))
                all.pred.ind.df <- rbind(all.pred.ind.df, pred.ind.df)

            }#if model selection
            all.out.df <- rbind(all.out.df, out.df)
        }#if
    }## END for (pl in pheno.list)
}## END for(ml in mutant.list)

rownames(all.out.df) <- NULL
write.table(all.out.df, paste0(path.out, pg, "_stats_", filename.str,".csv")
          , row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",", append=FALSE)

if(model.type=="LASSO"){
    rownames(all.pred.ind.df)<-NULL
    write.table(all.pred.ind.df, paste0(path.out, pg, "_lasso_pred_ind_", filename.str,".csv")
              , row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",", append=FALSE)
}




