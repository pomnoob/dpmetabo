
#Normalization using QC-RLSC method
library(tidyverse)
library(NormalizeMets)

#Import data for sample list
sampdata=read_csv("normpeno_pos.csv") %>% 
  remove_rownames() %>% column_to_rownames(var="sample") %>% 
  as.data.frame()

#Import feature data
fdata <- read_csv("normpro_pos.csv") %>% 
  remove_rownames() %>% column_to_rownames(var = "X1") %>%
  as.data.frame()

#Normalization
rlsct <- NormQcsamples(fdata,sampdata,method = c("rlsc"),span=0, lg=FALSE,saveoutput = TRUE,outputname = "QCRSLC_norm")

#Try out statTarget
library(statTarget)

shiftCor("statTargetlist_pos.csv","statTargetfeature_pos.csv", Frule = 0.9, MLmethod = "QCRLSC", QCspan = 0,imputeM = "KNN")

#read normalizaed data
qcrlsc_pos <- read.csv("statTarget/shiftCor/After_shiftCor/shift_sample_cor.csv",header=FALSE,check.names = FALSE, stringsAsFactors = FALSE)

colnames(qcrlsc_pos) <- qcrlsc_pos[1,]

qcrlsc_pos <- as.data.frame(qcrlsc_pos[-1,])

#MSTUS to eliminate possbile osmolity effects
qcrlsc_pos$rowsum <- rowSums(qcrlsc_pos[3:1853],na.rm = FALSE)
for (i in 3:1853){qcrlsc_pos[,i] <- qcrlsc_pos[,i]/qcrlsc_pos[,"rowsum"]}

qc_pos_ac <- qcrlsc_pos[qcrlsc_pos$class==1|qcrlsc_pos$class==3,] %>% dplyr::arrange(sample) %>% dplyr::slice(-27,-28) %>% dplyr::select(-rowsum)

#MetaboAnalystR
qc_pos_ac$class <- c(1:16,-1:-16)
write.csv(qc_pos_ac,"qc_ac.csv",row.names = FALSE)

#MUVR
qc_pos_acm <- qc_pos_ac
for(i in 1:16){qc_pos_acm[i+32,3:1853] <- qc_pos_ac[i,3:1853]-qc_pos_ac[i+16,3:1853]}

muvr_ac <- qc_pos_acm %>% dplyr::slice(33:48) %>% dplyr::select(-class)

muvr_ac$sample <- c(1:16)

# Set method parameters
nCore=detectCores()-1   # Number of processor threads to use
nRep=nCore              # Number of MUVR repetitions
nOuter=8                # Number of outer cross-validation segments
varRatio=0.8            # Proportion of variables kept per iteration 
method='RF'             # Selected core modelling algorithm
# Set up parallel processing
cl=makeCluster(nCore)   
registerDoParallel(cl)

# Perform modelling
MLModel = MUVR(X=muvr_ac, ML=TRUE, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
# 1.0 mins using 7 threads on a Mac Powerbook Pro mid-2015 with 2,8 GHz Intel Core i7.
# Stop parallel processing
stopCluster(cl)
# Examine model performance and output
MLModel$nVar                   # Number of variables for min, mid and max models
MLModel$miss 

plotVAL(MLModel)

plotMV(MLModel, model='min')

getVIP(MLModel, model='min')




