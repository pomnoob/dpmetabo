library(tidyverse)

# Load injection order file
inj_order <- read.csv("data/inject_order.csv", stringsAsFactors = F)

# Load negitive mode data file
neg_ori <- read.csv("data/neg-ori.csv", stringsAsFactors = F)

neg_ms2 <- dplyr::select(neg_ori, id, MS2.name)  #select metabolites id and MS2 names
neg_ms2$MS2.name <- as.character(neg_ms2$MS2.name)


neg_pro <- dplyr::select(neg_ori, -2, -3, -4, -5, -6, -7, -8)  #select feature profiles for each sample

sampid <- names(neg_pro) %>% str_replace("[0-9][0-9]_neg\\.[ABCD]|[0-9][0-9]_neg\\.QC", "")
names(neg_pro) <- sampid
neg_pro <- dplyr::rename(neg_pro, QC10 = QC0010)

# Get files ready for statTarget sample list file
samp_list <- inj_order
samp_list$batch <- rep(c(1), time(78))
samp_list$class[samp_list$class == "q"] <- "NA"
samp_list <- dplyr::arrange(samp_list, order)
write.csv(samp_list, "data/samp_list.csv", row.names = F)
# Feature file
samp_fea_neg <- dplyr::rename(neg_pro, name = id)
write.csv(samp_fea_neg, "data/samp_fea_neg.csv", row.names = F)

# QC-RLSC normalization
library(statTarget)
shiftCor("data/samp_list.csv", "data/samp_fea_neg.csv", Frule = 0.9, MLmethod = "QCRLSC", QCspan = 0, imputeM = "KNN")

file.copy("statTarget/shiftCor/After_shiftCor/shift_all_cor.csv","data/qcrlsc-all-neg.csv")
qcnorm_neg <- read.csv("data/qcrlsc-all-neg.csv", header = F, stringsAsFactors = F)
colnames(qcnorm_neg) <- qcnorm_neg[1, ]
qcnorm_neg <- as.data.frame(qcnorm_neg[-1, ])

# MSTUS to eliminate osmolity effect
qcnorm_neg$rowsum <- rowSums(qcnorm_neg[3:2222], na.rm = FALSE)
for (i in 3:2222) {
  qcnorm_neg[, i] <- qcnorm_neg[, i]/qcnorm_neg[, "rowsum"]
}
qcnorm_neg <- as.data.frame(qcnorm_neg[, -2223])

#data for A and D, data for A,B and C
library(tidyverse)

norm_neg <- qcnorm_neg[!is.na(qcnorm_neg$class),]#remove qc samples

normAD_neg <-norm_neg[norm_neg$class=="a"|norm_neg$class=="d",]%>% dplyr::arrange(sample) %>% dplyr::slice(-27,-28)

normABC_neg<-norm_neg[norm_neg$class=="a"|norm_neg$class=="b"|norm_neg$class=="c",]%>% dplyr::arrange(sample)%>% dplyr::slice(-43,-44)

library(statTarget)
write.csv(norm_neg,"data/mstus_neg.csv",row.names = F)
statAnalysis("data/mstus_neg.csv",Frule = 0.8, normM = "NONE", imputeM = "KNN", glog = TRUE,scaling = "Pareto")
file.copy("statTarget/statAnalysis/scaleData_Pareto/ProcessedTable.csv","data/ts_all_neg.csv") #copy transformed and scaled data 


ts_neg <- read.csv("data/ts_all_neg.csv",stringsAsFactors = F,header = F)
colnames(ts_neg) <- ts_neg[1,]
ts_neg <- as.data.frame(ts_neg[-1,])
ts_col <- names(ts_neg) %>% str_replace("X","")
names(ts_neg) <- ts_col
colnames(ts_neg) <- neg_ms2$MS2.name[match(names(ts_neg),neg_ms2$id)]

col12 <- names(ts_neg)
col12[is.na(col12)] <- "sample"
names(ts_neg) <- col12
ts_neg <- ts_neg %>% dplyr::arrange(sample) %>% dplyr::slice(-43,-44,-61,-62)
ts_neg$class <- rep(c("a","b","c","d"),times=c(16,16,16,16))
ts_neg <- dplyr::select(ts_neg,sample,class,everything())



ts_negMS2 <- ts_neg %>% dplyr::select(-contains("unknown"))

#Write data file to csv: After QCRLSC, MSTUS, log transformation and Pareto scaling

write.csv(ts_neg,"data/ts all with MS2 names_neg.csv",row.names = F)

write.csv(ts_negMS2,"data/ts all with known MS2_neg.csv",row.names = F)


ab_negMS2 <- dplyr::slice(ts_negMS2,1:32)

ab_neg <- dplyr::slice(ts_neg,1:32)

ac_negMS2 <- dplyr::slice(ts_negMS2,1:16,33:48)


ac_neg <- dplyr::slice(ts_neg,1:16,33:48)

ad_negMS2 <- dplyr::slice(ts_negMS2,1:16,49:64)

ad_neg<- dplyr::slice(ts_neg,1:16,49:64)

bc_negMS2 <- dplyr::slice(ts_negMS2,17:48)

bc_neg <- dplyr::slice(ts_neg,17:48)

bd_negMS2 <- dplyr::slice(ts_negMS2,17:32,49:64)

bd_neg <- dplyr::slice(ts_neg,17:32,49:64)

cd_negMS2 <- dplyr::slice(ts_negMS2,33:64)

cd_neg <- dplyr::slice(ts_neg,33:64)

abc_negMS2 <- dplyr::slice(ts_negMS2,1:48)

abc_neg <- dplyr::slice(ts_neg,1:48)

abc_negMS2$time <- rep(c(1:3),times=c(16,16,16))
abc_negMS2$subject <- rep(c(1:16),times=3)
abc_negMS2 <- abc_negMS2 %>% dplyr::select(sample,time,subject,everything()) %>% dplyr::arrange(subject,time) %>% dplyr::select(-class)

write.csv(abc_negMS2,"data/time series ABC_neg.csv",row.names = F)


library(MetaboAnalystR)


mSet<-InitDataObjects("conc", "ts", FALSE)
mSet<-SetDesignType(mSet, "time0")
mSet<-Read.TextData(mSet, "data/time series ABC_neg.csv", "rowts", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
timese<-performMB(mSet, 10)
timese_result <- as.data.frame(timese$analSet)
write.csv(timese_result,"data/time series-ABC-results_neg.csv")

library(MetaboAnalystR)

ac_negMS2$Label=c(1:16,-1:-16)
ac_negMS2 <- ac_negMS2 %>% dplyr::select(sample,Label,everything()) %>% dplyr::select(-class)

write.csv(ac_negMS2,"data/ac_negMS2.csv",row.names = F)

ac_pls<-InitDataObjects("conc", "stat", TRUE)
ac_pls<-Read.TextData(ac_pls, "data/ac_negMS2.csv", "rowp", "disc");
ac_pls<-SanityCheckData(ac_pls)
ac_pls<-ReplaceMin(ac_pls);
ac_pls<-PreparePrenormData(ac_pls)
ac_pls<-Normalization(ac_pls, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

ac_pls<-PLSR.Anal(ac_pls, reg=TRUE)

ac_pls<-PlotPLSPairSummary(ac_pls, "figure/MetaboAnalystR_PLSDA/ac pls summary_neg", "png", 72, width=NA, 5)

ac_pls<-PlotPLS2DScore(ac_pls, "figure/MetaboAnalystR_PLSDA/ac pls 2D score plot_neg", "png", 72, width=NA, 1,2,0.95,0,0)

ac_pls<-PlotPLSLoading(ac_pls, "figure/MetaboAnalystR_PLSDA/ac pls loading_neg", "png", 72, width=NA, 1, 2);

ac_pls<-PLSDA.CV(ac_pls, "T",5, "Q2")
ac_pls<-PlotPLS.Classification(ac_pls, "figure/MetaboAnalystR_PLSDA/ac pls cross validation_neg", "png", 72, width=NA)

ac_pls_vip<-PlotPLS.Imp(ac_pls, "figure/MetaboAnalystR_PLSDA/ac pls_vip_neg", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)

ac_vip <- ac_pls_vip$analSet$plsda$vip.mat

ac_vip_d <- as.data.frame(ac_vip)
ac_vip_d <- dplyr::select(ac_vip_d,ac_vip_c1='Comp. 1')

write.csv(ac_vip_d,"data/ac plsda vip_neg.csv")


library(MetaboAnalystR)

ab_negMS2$Label=c(1:16,-1:-16)
ab_negMS2 <- ab_negMS2 %>% dplyr::select(sample,Label,everything()) %>% dplyr::select(-class)

write.csv(ab_negMS2,"data/ab_negMS2.csv",row.names = F)

ab_pls<-InitDataObjects("conc", "stat", TRUE)
ab_pls<-Read.TextData(ab_pls, "data/ab_negMS2.csv", "rowp", "disc");
ab_pls<-SanityCheckData(ab_pls)
ab_pls<-ReplaceMin(ab_pls);
ab_pls<-PreparePrenormData(ab_pls)
ab_pls<-Normalization(ab_pls, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

ab_pls<-PLSR.Anal(ab_pls, reg=TRUE)

ab_pls<-PlotPLSPairSummary(ab_pls, "figure/MetaboAnalystR_PLSDA/ab pls summary_neg", "png", 72, width=NA, 5)

ab_pls<-PlotPLS2DScore(ab_pls, "figure/MetaboAnalystR_PLSDA/ab pls 2D score plot_neg", "png", 72, width=NA, 1,2,0.95,0,0)

ab_pls<-PlotPLSLoading(ab_pls, "figure/MetaboAnalystR_PLSDA/ab pls loading_neg", "png", 72, width=NA, 1, 2);

ab_pls<-PLSDA.CV(ab_pls, "T",5, "Q2")
ab_pls<-PlotPLS.Classification(ab_pls, "figure/MetaboAnalystR_PLSDA/ab pls cross validation_neg", "png", 72, width=NA)

ab_pls_vip<-PlotPLS.Imp(ab_pls, "figure/MetaboAnalystR_PLSDA/ab pls_vip_neg", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)

ab_vip <- ab_pls_vip$analSet$plsda$vip.mat

ab_vip_d <- as.data.frame(ab_vip)
ab_vip_d <- dplyr::select(ab_vip_d,ab_vip_c1='Comp. 1')

write.csv(ab_vip_d,"data/ab plsda vip_neg.csv")




library(MetaboAnalystR)

ad_negMS2$Label=c(1:16,-1:-16)
ad_negMS2 <- ad_negMS2 %>% dplyr::select(sample,Label,everything()) %>% dplyr::select(-class)

write.csv(ad_negMS2,"data/ad_negMS2.csv",row.names = F)

ad_pls<-InitDataObjects("conc", "stat", TRUE)
ad_pls<-Read.TextData(ad_pls, "data/ad_negMS2.csv", "rowp", "disc");
ad_pls<-SanityCheckData(ad_pls)
ad_pls<-ReplaceMin(ad_pls);
ad_pls<-PreparePrenormData(ad_pls)
ad_pls<-Normalization(ad_pls, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

ad_pls<-PLSR.Anal(ad_pls, reg=TRUE)

ad_pls<-PlotPLSPairSummary(ad_pls, "figure/MetaboAnalystR_PLSDA/ad pls summary_neg", "png", 72, width=NA, 5)

ad_pls<-PlotPLS2DScore(ad_pls, "figure/MetaboAnalystR_PLSDA/ad pls 2D score plot_neg", "png", 72, width=NA, 1,2,0.95,0,0)

ad_pls<-PlotPLSLoading(ad_pls, "figure/MetaboAnalystR_PLSDA/ad pls loading_neg", "png", 72, width=NA, 1, 2);

ad_pls<-PLSDA.CV(ad_pls, "T",5, "Q2")
ad_pls<-PlotPLS.Classification(ad_pls, "figure/MetaboAnalystR_PLSDA/ad pls cross validation_neg", "png", 72, width=NA)

ad_pls_vip<-PlotPLS.Imp(ad_pls, "figure/MetaboAnalystR_PLSDA/ad pls_vip_neg", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)

ad_vip <- ad_pls_vip$analSet$plsda$vip.mat

ad_vip_d <- as.data.frame(ad_vip)
ad_vip_d <- dplyr::select(ad_vip_d,ad_vip_c1='Comp. 1')

write.csv(ad_vip_d,"data/ad plsda vip_neg.csv")
```

```{r PLS analysis of BC using MetbcoAnalystR,include=FALSE}
library(MetaboAnalystR)

bc_negMS2$Label=c(1:16,-1:-16)
bc_negMS2 <- bc_negMS2 %>% dplyr::select(sample,Label,everything()) %>% dplyr::select(-class)

write.csv(bc_negMS2,"data/bc_negMS2.csv",row.names = F)

bc_pls<-InitDataObjects("conc", "stat", TRUE)
bc_pls<-Read.TextData(bc_pls, "data/bc_negMS2.csv", "rowp", "disc");
bc_pls<-SanityCheckData(bc_pls)
bc_pls<-ReplaceMin(bc_pls);
bc_pls<-PreparePrenormData(bc_pls)
bc_pls<-Normalization(bc_pls, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

bc_pls<-PLSR.Anal(bc_pls, reg=TRUE)

bc_pls<-PlotPLSPairSummary(bc_pls, "figure/MetaboAnalystR_PLSDA/bc pls summary_neg", "png", 72, width=NA, 5)

bc_pls<-PlotPLS2DScore(bc_pls, "figure/MetaboAnalystR_PLSDA/bc pls 2D score plot_neg", "png", 72, width=NA, 1,2,0.95,0,0)

bc_pls<-PlotPLSLoading(bc_pls, "figure/MetaboAnalystR_PLSDA/bc pls loading_neg", "png", 72, width=NA, 1, 2);

bc_pls<-PLSDA.CV(bc_pls, "T",5, "Q2")
bc_pls<-PlotPLS.Classification(bc_pls, "figure/MetaboAnalystR_PLSDA/bc pls cross validation_neg", "png", 72, width=NA)

bc_pls_vip<-PlotPLS.Imp(bc_pls, "figure/MetaboAnalystR_PLSDA/bc pls_vip_neg", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)

bc_vip <- bc_pls_vip$analSet$plsda$vip.mat

bc_vip_d <- as.data.frame(bc_vip)
bc_vip_d <- dplyr::select(bc_vip_d,bc_vip_c1='Comp. 1')

write.csv(bc_vip_d,"data/bc plsda vip_neg.csv")
```



```{r PLSDA and time series based variable selection, include=F}
library(tidyverse)
timese_result <- timese_result %>% dplyr::select(-type) %>% add_rownames("mid")
ac_vip_d <-rownames_to_column(ac_vip_d,"mid")
ab_vip_d <- rownames_to_column(ab_vip_d,"mid")
bc_vip_d <- rownames_to_column(bc_vip_d,"mid")
pls_times_result_neg <- left_join(timese_result,ac_vip_d,by="mid")
pls_times_result_neg <- left_join(pls_times_result_neg,ab_vip_d,by="mid")
pls_times_result_neg <- left_join(pls_times_result_neg,bc_vip_d,by="mid")
write.csv(pls_times_result_neg,"data/PLSDA and time series results_neg.csv")
pls.time.vip <- dplyr::filter(pls_times_result_neg,ac_vip_c1>2|ab_vip_c1>2|bc_vip_c1>2)

pls.time.vip2 <- dplyr::filter(pls.time.vip,Hotelling.T2>99)
write.csv(pls.time.vip,"data/metabolites VIP over 2.csv",row.names = F)
```



```{r data files for plotting, include=F}

mstus_neg <- read.csv("data/mstus_neg.csv",stringsAsFactors = F,header = F)
colnames(mstus_neg) <- mstus_neg[1,]
mstus_neg <- as.data.frame(mstus_neg[-1,])
mstus_col <- names(mstus_neg) %>% str_replace("X","")
names(mstus_neg) <- mstus_col
mstus_neg <- dplyr::select(mstus_neg,-class)

colnames(mstus_neg) <- neg_ms2$MS2.name[match(names(mstus_neg),neg_ms2$id)]


colmstus <- names(mstus_neg)
colmstus[is.na(colmstus)] <- "sample"
names(mstus_neg) <- colmstus
mstus_neg <- mstus_neg %>% dplyr::arrange(sample)%>% dplyr::slice(-43,-44,-61,-62)
mstus_neg$class <- rep(c("A diet","B diet","C diet","D diet"),times=c(16,16,16,16))
mstus_neg <- dplyr::select(mstus_neg,sample,class,everything())

mstusABC_neg <- dplyr::slice(mstus_neg,1:48)

mid <- pls.time.vip$mid
sampid <- mstusABC_neg$sample
classid <- mstusABC_neg$class

vip.file.neg <- mstusABC_neg[,colnames(mstusABC_neg)%in% mid]
vip.file.neg$sample <- sampid
vip.file.neg$class <- classid



vip.p.neg <- gather(vip.file.neg,key="mid",value = "Abundance",-sample,-class)

```

library(tidyverse)
p1 <- ggplot(data=vip.p.neg,mapping = aes(x=class,y=Abundance))+geom_boxplot()+facet_wrap(~mid,nrow = 4)
p1+ylim(0,0.002)+xlab("Dietary pattern")+ggtitle("尿液代谢产物在不同膳食模式下的峰面积变化图1")+theme(plot.title=element_text(size=rel(1.5), lineheight=.9, hjust=0.5,family="Times"))
