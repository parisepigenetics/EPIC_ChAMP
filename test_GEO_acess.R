
library(GEOquery)
library(Biobase)
library(ChAMP)

getGEOSuppFiles('GSE191276')
GSE <- getGEO(GEO = 'GSE191276')

pd <- pData(phenoData(GSE[[1]]))

untar(tarfile = 'GSE191276/GSE191276_RAW.tar')
idat_files <- list.files(pattern = 'idat.gz')

for(i in 1:length(idat_files)){
  gunzip(filename = idat_files[i], destname = gsub("[.]gz$", "", idat_files[i]))
}

list_files <- data.frame(list.files(pattern = "idat"))
colnames(list_files) <- "Basename"
list_files_p <- NULL
list_files_p <- unlist(lapply(basename(as.character(list_files$Basename)), function(x) (strsplit(x, "_Gr|_Re")[[1]][1])))
list_files_p <- data.frame(list_files_p)
colnames(list_files_p) <- "Basename"

ChAMP_txt <- data.frame(
  Sample_Name=pd$title, 
  Chips=pd$channel_count, 
  mutation=pd$organism_ch1,
  Sex=pd$'Sex:ch1',
  Sample_status=pd$'case control status:ch1',
  Mutation_status=pd$organism_ch1,
  Basename=list_files_p,
  tissue=pd$source_name_ch1,
  stringsAsFactors=F)

write.table(ChAMP_txt,"pD_ChAMP_test.txt", row.names=F, quote=F, sep=",")

###############################################################################################################################
# Test data_850K : load_pD pheno_data
pD_EPIC <- read.csv("pD_ChAMP_test.txt", header=T,sep=",", stringsAsFactors=F)
head(pD_EPIC)

#View(pD_EPIC)
pD_EPIC$Sample_Name
colnames(pD_EPIC)

#Format pD_files EPIC for ChAMP (change col names)
Sentrix_ID=unlist(lapply(basename(as.character(pD_EPIC$Basename)), function(x) strsplit(gsub("_R",":R",x),":")[[1]][[1]]))
Sentrix_Position=unlist(lapply(basename(as.character(list_files$Basename)), function(x) (strsplit(x, "_")[[1]][3])))
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

# a revoir 
#ChAMP_csv <- na.omit(ChAMP_csv)

write.table(ChAMP_csv,"ChAMP.csv", row.names=F, quote=F, sep=",")




