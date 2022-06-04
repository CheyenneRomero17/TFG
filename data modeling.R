### R code for TFG
### Cheyenne Romero Freijo

setwd("~/Uni/MatCAD/4t/TFG/data")

############################# UPLOADING THE DATA #############################
data <- read.table('data.txt', header = TRUE)
r1 <- read.table('Risk_A1.txt', header = TRUE)
r2 <- read.table('Risk_A2GO.txt', header = TRUE)
r3 <- read.table('Risk_A3.txt', header = TRUE)
p1 <- read.table('Proxy_A1.txt', header = TRUE)
p2 <- read.table('Proxy_A2GO.txt', header = TRUE)
p3 <- read.table('Proxy_A3.txt', header = TRUE)
volumes <- read.table('added_volumes.txt', header = TRUE)

############################# ANALYSING THE DATA #############################

## Outliers with outlier library
library(outliers)
# Intracranvial volum outliers
idx <- which(outlier(data$IntraCranialVol, logical = TRUE) == TRUE)
patient <- unique(data$RID[idx]) # patient with RID 954, with an IntraCranialVol of 971180.

data2 <- data[-idx,]

# Hippocampal subregions volumes outliers
idx <- NULL
for (i in colnames(volumes)[2:ncol(volumes)]){
  idx <- c(idx, which(outlier(volumes[,i], logical = TRUE) == TRUE))
}

patient <- unique(volumes$Subject[idx]) # 419 1293 558 2394 5296

############################## GATHERING THE DATA ##############################
risk <- merge(r1, r2, all = TRUE)
risk <- merge(risk, r3, all = TRUE)

proxy <- merge(p1, p2, all = TRUE)
proxy <- merge(proxy, p3, all = TRUE)

# To free memory
rm(data); rm(r1); rm(r2); rm(r3); rm(p1); rm(p2); rm(p3); rm(idx); rm(patient)

################################## FUNCTIONS ##################################
data_gather <- function(data, snp_data, volumes, snp, sr){
  idx <-unique(data$RID[data$RID %in% snp_data$IID[!is.na(snp_data[,snps[i]])]])
  data <- data[data$RID %in% idx,]
  data$SNP <- NA
  data$h_subregion <- NA
  
  for (id in idx){
    data$SNP[data$RID==id] <- snp_data[,snp][snp_data$IID==id]
    data$h_subregion[data$RID==id] <- volumes[,sr][volumes$Subject==id]
  }
  
  return(data)
}
