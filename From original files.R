library(tidyverse)

# Load injection order file
inj_order <- read.csv("data/inject_order.csv", stringsAsFactors = F)

# Load positive mode data file
pos_ori <- read.csv("data/pos-ori.csv", stringsAsFactors = F)

pos_ms2 <- dplyr::select(pos_ori, id, MS2.name)  #select metabolites id and MS2 names
pos_ms2$MS2.name <- as.character(pos_ms2$MS2.name)


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

#mixOMics
library(mixOmics)
data(vac18)
gene <- vac18$genes
sample <- vac18$sample
design <- data.frame(sample)
stimulation <- vac18$stimulation

ml <- data.frame(rep(c(1:16),times=2))
ml <- dplyr::rename(ml,sample=1)
group <- rep(c("a","d"),times=c(16,16))

tsAD_pos_m <- tsAD_pos[,-1]
rownames(tsAD_pos_m) <- tsAD_pos[,1]

## pls-da AD
plsda.ad<- plsda(tsAD_pos_m, Y = group, multilevel = ml, ncomp = 10)

plotIndiv(plsda.ad, ind.names = ml,
          group = group, 
          legend = TRUE,
          title = 'AD multilevel PLS-DA, comp 1 - 2')

background = background.predict(plsda.ad, comp.predicted=2, dist = "max.dist") 

plotIndiv(plsda.ad, comp = 1:2,
          group = group, ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.ad <- perf(plsda.ad, validation = "Mfold", folds = 5, 
                         progressBar = FALSE, auc = TRUE, nrepeat = 10) 

plot(perf.plsda.ad, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

auc.plsda.ad = auroc(plsda.ad, roc.comp = 2)

##Tuning spls-da AD
list.keepX <- c(1:10,  seq(20, 300, 10))

tune.splsda.ad <- tune.splsda(tsAD_pos_m, Y = group, multilevel = ml, ncomp = 2, validation = 'Mfold', folds = 5, 
                                 progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                 test.keepX = list.keepX, nrepeat = 10, cpus = 2)

tune.splsda.ad7 <- tune.splsda(tsAD_pos_m, Y = group, multilevel = ml, ncomp = 7, validation = 'Mfold', folds = 5, 
                              progressBar = TRUE, dist = 'max.dist', measure = "BER",
                              test.keepX = list.keepX, nrepeat = 10, cpus = 2)

error <- tune.splsda.ad7$error.rate 

ncomp <- tune.splsda.ad7$choice.ncomp$ncomp


select.keepX <- tune.splsda.ad7$choice.keepX[1:2]  # optimal number of variables to select

splsda.ad <- splsda(tsAD_pos_m, Y = group, multilevel = ml,ncomp = 2, keepX = select.keepX) 

plotIndiv(splsda.ad, 
          group = group, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')

auc.splsda.ad = auroc(splsda.ad, roc.comp = 1)
auc.splsda.ad2 = auroc(splsda.ad, roc.comp = 2)


## pls-da
plsda.abc<- plsda(tsABC_pos_m, Y = group3, multilevel = ml3, ncomp = 10)

plotIndiv(plsda.abc, ind.names = ml3,
          group = group3, 
          legend = TRUE,
          title = 'ABC multilevel PLS-DA, comp 1 - 2')

background = background.predict(plsda.abc, comp.predicted=2, dist = "max.dist") 

plotIndiv(plsda.abc, comp = 1:2,
          group = group3, ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.abc <- perf(plsda.abc, validation = "Mfold", folds = 5, 
                      progressBar = FALSE, auc = TRUE, nrepeat = 10) 

plot(perf.plsda.abc, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

perf.plsda.abc$choice.ncomp

auc.plsda.ad = auroc(plsda.abc, roc.comp = 3)


##Tuning spls-da ABC
list.keepX <- c(1:10,  seq(20, 300, 10))

tune.splsda.abc <- tune.splsda(tsABC_pos_m, Y = group3, multilevel = ml3, ncomp = 3, validation = 'Mfold', folds = 5, 
                              progressBar = TRUE, dist = 'max.dist', measure = "BER",
                              test.keepX = list.keepX, nrepeat = 10, cpus = 2)


error_abc <- tune.splsda.abc$error.rate 

ncomp_abc <- tune.splsda.abc$choice.ncomp$ncomp


select.keepX_abc <- tune.splsda.abc$choice.keepX[1:ncomp_abc]  # optimal number of variables to select

splsda.abc <- splsda(tsABC_pos_m, Y = group3, multilevel = ml3,ncomp = ncomp_abc, keepX = select.keepX_abc) 

plotIndiv(splsda.abc, 
          group = group3, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on ABC, comp 1 & 2')

auc.splsda.abc = auroc(splsda.abc, roc.comp = ncomp_abc)


set.seed(40) # for reproducibility, only when the `cpus' argument is not used
# takes about 1 min to run
perf.abc <- perf(splsda.abc, validation = "Mfold", folds = 5,
                   dist = 'max.dist', nrepeat = 10,
                   progressBar = FALSE) 

# perf.srbct  # lists the different outputs
perf.abc$error.rate

plot(perf.abc, col = color.mixo(5))

par(mfrow=c(1,2))
plot(perf.abc$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.abc$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)

par(mfrow=c(1,1))

# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda.abc, comp = 1)$name, 
                  names(perf.abc$features$stable[[1]]))

#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.abc$features$stable[[1]][ind.match])

data.frame(selectVar(splsda.abc, comp = 1)$value, Freq)
