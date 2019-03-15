
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