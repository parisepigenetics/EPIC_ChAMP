# Load libraries 
library(ChAMP)
library(dplyr)
library(yaml)

# Parse arguments :  
args <- commandArgs(TRUE)
result_path <- args[1]

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

idat_path <- paste0( getwd(), "/my_bank")

myLoad <- champ.load(directory = idat_path ,
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="EPIC")

head(myLoad)

## QC procedure
champ.QC(beta = myLoad$beta,
         pheno=myLoad$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=T,
         PDFplot=T,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir=paste0(result_path, "/CHAMP_QC/"))
#dev.off()

## Normalization, by default, method = BMIQ ("PBC","BMIQ","SWAN")
myNorm <- champ.norm(beta=myLoad$beta, arraytype="EPIC", method="BMIQ", resultsDir=paste0(result_path, "/CHAMP_NORM/"))
head(myNorm)
dim(myNorm)
#getwd()
write.table(myNorm, paste0(result_path,"/beta_EPIC_MD"), row.names=F, quote=F, sep="\t")

## Combat correction for bacth effect (450K vs EPIC Samples, samples prep, date etc...)
#myNorm <- champ.runCombat(variablename = "Sample_Status", batchname=c("Sample_Plate"))

## SVD
champ.SVD(beta=myNorm %>% as.data.frame(),pd=myLoad$pd, resultsDir=paste0(result_path, "/CHAMP_SVD/"))

##Cell Heterogeneity (if tissue=Blood)
myRefBase <- champ.refbase(beta=myNorm %>% as.data.frame() ,arraytype="EPIC")
myRefBase$CellFraction
head(myRefBase$CorrectedBeta)

### MDS plot on probes for cell composition
beta_cell <- myNorm[rownames(myNorm) %in% rownames(CellTypeMeans450K),]
#pdf(paste0(QCDir,"/MDS_plot_cpg_cell_compo.pdf"),onefile=T, width=25,height=20)
champ.QC(beta=beta_cell,PDFplot=F, dendrogram=F, densityPlot=T, resultsDir= paste0(result_path,  "/CHAMP_QC_CELL/"))
SVD <- svd(beta_cell)
#dev.off()

# MDS plot with Normalization
champ.QC(beta = myNorm,
         pheno=myLoad$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=T,
         PDFplot=T,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir=paste0(result_path, "/CHAMP_QC_NORM/"))
#dev.off()

## Calling DMPs
myDMP<- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Status, compare.group=c("ORC","CTL"),arraytype="EPIC")

#myDMP <- champ.DMP(beta = myRefBase$CorrectedBeta,pheno=myLoad$pd$Sample_Group, arraytype="EPIC")
#myDMP <- champ.DMP(beta = myNorm,arraytype="EPIC")


## Calling DMR
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter")
head(myDMR)

## GSEA
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450K")
myebayGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
head(myGSEA$DMP)
head(myGSEA$DMR)

# dendogram_beta_850K
head(myNorm)
dim(myNorm)

#dendogram
df<-t(as.matrix(myNorm))
dim(df)
row.names(df)

# Compute distances and hierarchical clustering
dd <- dist(df, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

# Default plot
plot(hc)



# save the R Session in RData folder
save.image(file =paste0(result_path, "/ChAMP.RData")) 




