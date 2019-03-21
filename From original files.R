library(tidyverse)

# Load injection order file
inj_order <- read.csv("data/inject_order.csv", stringsAsFactors = F)

# Load positive mode data file
pos_ori <- read.csv("data/pos-ori.csv", stringsAsFactors = F)

pos_ms2 <- dplyr::select(pos_ori, id, MS2.name)  #select metabolites id and MS2 names
pos_ms2$MS2.name <- as.character(pos_ms2$MS2.name)
pos_ms2$MS2.name[is.na(pos_ms2$MS2.name)] <- "unknown"

pos_pro <- dplyr::select(pos_ori, -2, -3, -4, -5, -6, -7, -8)  #select feature profiles for each sample

sampid <- names(pos_pro) %>% str_replace("[0-9][0-9]_pos\\.[ABCD]|[0-9][0-9]_pos\\.QC", "")
names(pos_pro) <- sampid
pos_pro <- dplyr::rename(pos_pro, QC10 = QC0010)

# Get files ready for statTarget sample list file
samp_list <- inj_order
samp_list$batch <- rep(c(1), time(78))
samp_list$class[samp_list$class == "q"] <- "NA"
samp_list <- dplyr::arrange(samp_list, order)
write.csv(samp_list, "data/samp_list.csv", row.names = F)
# Feature file
samp_fea <- dplyr::rename(pos_pro, name = id)
write.csv(samp_fea, "data/samp_fea.csv", row.names = F)

# QC-RLSC normalization
library(statTarget)
shiftCor("data/samp_list.csv", "data/samp_fea.csv", Frule = 0.9, MLmethod = "QCRLSC", QCspan = 0, imputeM = "KNN")
qcnorm_pos <- read.csv("statTarget/shiftCor/After_shiftCor/shift_all_cor.csv", header = F, stringsAsFactors = F)
colnames(qcnorm_pos) <- qcnorm_pos[1, ]
qcnorm_pos <- as.data.frame(qcnorm_pos[-1, ])

# MSTUS to eliminate osmolity effect
qcnorm_pos$rowsum <- rowSums(qcnorm_pos[3:1853], na.rm = FALSE)
for (i in 3:1853) {
    qcnorm_pos[, i] <- qcnorm_pos[, i]/qcnorm_pos[, "rowsum"]
}
qcnorm_pos <- as.data.frame(qcnorm_pos[, -1854])

#data for A and D, data for A,B and C

norm_pos <- qcnorm_pos[!is.na(qcnorm_pos$class),]#remove qc samples

normAD_pos <-norm_pos[norm_pos$class=="a"|norm_pos$class=="d",]%>% dplyr::arrange(sample) %>% dplyr::slice(-27,-28)

normABC_pos<-norm_pos[norm_pos$class=="a"|norm_pos$class=="b"|norm_pos$class=="c",]%>% dplyr::arrange(sample)%>% dplyr::slice(-43,-44)

#log transformation and Pareto scaling
library(statTarget)
write.csv(normAD_pos,"data/normAD_pos.csv",row.names = F)
statAnalysis("data/normAD_pos.csv",Frule = 0.8, normM = "NONE", imputeM = "KNN", glog = TRUE,scaling = "Pareto")
file.copy("statTarget/statAnalysis/scaleData_Pareto/ProcessedTable.csv","data/ts_AD.csv") #copy transformed and scaled data 
tsAD_pos <- read.csv("data/ts_AD.csv",stringsAsFactors = F,header = F)
colnames(tsAD_pos) <- tsAD_pos[1,]
tsAD_pos <- as.data.frame(tsAD_pos[-1,])
tsAD_col <- names(tsAD_pos) %>% str_replace("X","")
names(tsAD_pos) <- tsAD_col

tsAD_pos <- dplyr::rename(tsAD_pos,sample="" )


