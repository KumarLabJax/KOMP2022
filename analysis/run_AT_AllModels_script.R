rm(list=ls())

batch = "2022-05-13" #to change for different version of data

path.wd <- "/projects/compsci/vmp/USERS/heh/KOMP2022/"#project dir
path.meta <- paste0(path.wd, "data/")
path.data <- paste0(path.wd, "output/", batch, "_out/")
path.out <- paste0(path.wd, "output/", batch, "_AT/")
path.ro <- paste0(path.wd, "output/", batch, "_AT/Rout/")

if (!dir.exists(path.data))
{ dir.create(path.data) } #create folder

if (!dir.exists(path.out))
{ dir.create(path.out) } #create folder

if (!dir.exists(path.ro))
{ dir.create(path.ro) } #create folder

R.file <- paste0(path.wd, "analysis/AT_AllModels_script.R")

###############
## R  Params ##
###############

##model.type <- "IMPC"
##model.type <- "REGRES"
model.type <- "LASSO"

date.str <- batch ## add more info here if needed.
wing.size <- 1   #months
verbose <- 0     #1:TRUE 0:FALSE

filename.str <- paste0(date.str,"_wing",wing.size,"months")
filename.str

##Read meta data
meta.data.file <- paste0(path.meta,"KOMP2018_Meta_Table_v1.2.csv")
meta.data <- read.csv(meta.data.file, header=TRUE)
meta.data <- meta.data[meta.data$Domain != "EKG",]

##Extract phenotype group list
pheno.group.list <- unique(meta.data$Domain)

job <- paste0(path.wd, "analysis/run_AT_AllModels_script.sh")
write("sleep 3",job, append = F)
for(pg in pheno.group.list){
        prefix <- paste0("KOMP_",pg,"_",model.type,"_",filename.str)
        Rout.file <- paste0(path.ro, prefix,".Rout")
        write(paste0("R CMD BATCH --no-restore --no-save --no-readline",
                                      " -phenogroup=", pg,
                                      " -modeltype=", model.type,
                                      " -filename=", filename.str,
                                      " -wingsize=", wing.size,
                                      " -verbose=", verbose,
                                      " ", R.file ,           # R file
                                      " ", Rout.file,
                                      " &" ),      # Rout file
              job, append = T)
}

