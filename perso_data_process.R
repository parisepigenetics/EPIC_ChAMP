
# Load libraries 
library(ChAMP)
library(dplyr)
library(yaml)


# Test data_850K : load_pD pheno_data
pD_EPIC <- read.csv("config/pD_MD.txt", header=T,sep=",", stringsAsFactors=F)
head(pD_EPIC)

#View(pD_EPIC)
pD_EPIC$Sample_Name
colnames(pD_EPIC)

#Format pD_files EPIC for ChAMP (change col names)
Sentrix_ID=unlist(lapply(basename(as.character(pD_EPIC$Basename)), function(x) strsplit(gsub("_R",":R",x),":")[[1]][[1]]))
Sentrix_Position=unlist(lapply(basename(as.character(pD_EPIC$Basename)), function(x) strsplit(gsub("_R",":R",x),":")[[1]][[2]]))
head(pD_EPIC)

ChAMP_csv <- data.frame(
  Sample_Name=pD_EPIC$Sample_Name,
  Sample_Plate=pD_EPIC$Chips,
  Sample_Group=pD_EPIC$mutation,
  Sample_Group_2=pD_EPIC$Sex,
  Sample_GSE_ID="",
  Sample_Status=pD_EPIC$Sample_status,
  Mutation_Status=pD_EPIC$Mutation_status,
  Pool_ID="",
  Project="",
  Sample_Well="",
  Sentrix_ID=Sentrix_ID,
  Sentrix_Position=Sentrix_Position,
  Basename=pD_EPIC$Basename,
  Tissue=pD_EPIC$tissue,
  stringsAsFactors=F)

head(ChAMP_csv)

#Filters pD
# keeping all Sample_Status == CTL / Discovery / Validation & Mutation_Status == LOF / abs
# keeping samples with tissue == Blood DNA or Cell type
factor(ChAMP_csv$Sample_Status)
ChAMP_f_csv <- subset(ChAMP_csv, Tissue=="LCL" & Sample_Status %in% c("CTL","Sotos","ATRX","ORC") & Mutation_Status %in% c("LOF","abs"))
head(ChAMP_f_csv)
dim(ChAMP_f_csv)

#CHAMP PROCEDURES_EPIC

# Load data, setwd() in data dir with all idat to create csv file of pD define in Samples Selections
write.table(ChAMP_f_csv,"pD_ChAMP.csv", row.names=F, quote=F, sep=",")






