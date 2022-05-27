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
idx <- which(outlier(data$IntraCranialVol, logical = TRUE) == TRUE)
patient <- unique(data$RID[idx]) # patient with RID 954, with an IntraCranialVol of 971180.

data2 <- data[-(which(data$RID == patient)),]
volumes <- volumes[-(which(volumes$Subject == patient)),]

############################## GATHERING THE DATA ##############################
risk <- merge(r1, r2, all = TRUE)
risk <- merge(risk, r3, all = TRUE)

proxy <- merge(p1, p2, all = TRUE)
proxy <- merge(proxy, p3, all = TRUE)

risk <- risk[-(which(risk$IID == patient)),]
proxy <- proxy[-(which(proxy$IID == patient)),]

# To free memory
rm(data); rm(r1); rm(r2); rm(r3); rm(p1); rm(p2); rm(p3)

################################## FUNCTIONS ##################################
# install.packages('nlme')
library(nlme)

data_gather <- function(data, snp_data, volumes, snp, sr){
  data <- data[data$RID %in% snp_data$IID,c(1,2,7,8,9,10,11,12)]
  data$SNP <- NA
  data$h_subregion <- NA
  
  data <- data[data$RID %in% which(!is.na(snp_data[,snp])),]
  
  for (id in data$RID){
    data$SNP[data$RID==id] <- snp_data[,snp][snp_data$IID==id]
    data$h_subregion[data$RID==id] <- volumes[,sr][volumes$Subject==id]
  }
  
  return(data)
}
