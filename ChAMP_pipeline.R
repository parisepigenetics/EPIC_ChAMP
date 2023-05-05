# Load libraries 
library(ChAMP)
library(dplyr)
library(yaml)

# Parse arguments :  
args <- commandArgs(TRUE)
result_path <- args[1]

# load the config file
yaml.file <- yaml.load_file('config.yml')

# extract the information from the yaml file
RAW_DATA_PATH <- yaml.file$RAW_DATA_PATH
GSE_NUM <- yaml.file$GSE_NUM
GEO <- yaml.file$GEO
FULL <- yam.file$FULL

if(GEO == TRUE){
  idat_path <- paste0(getwd(), "/", GSE_NUM)
}else{
  idat_path <- paste0(getwd(), "/", RAW_DATA_PATH)
}

if(FULL == TRUE){
  
  champ.process(directory = idat_path)
  
}else{
  

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

} 


# save the R Session in RData folder
save.image(file =paste0(result_path, "/ChAMP.RData")) 




