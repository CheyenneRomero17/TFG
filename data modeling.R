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

############################## GATHERING THE DATA ##############################
risk <- merge(r1, r2, all = TRUE)
risk <- merge(risk, r3, all = TRUE)

proxy <- merge(p1, p2, all = TRUE)
proxy <- merge(proxy, p3, all = TRUE)

risk <- risk[-(which(risk$IID == patient)),]
proxy <- proxy[-(which(proxy$IID == patient)),]

#################################### MODELS ####################################
# install.packages('lme4')
# library(lme4)
# install.packages('nlme')
library(nlme)

tmp <- data2[data2$RID %in% risk$IID, c(1,2,7,8,9,10,11,12)]
tmp$SNP <- NA
tmp$h_subregion <- NA

##### Simple Model #####
simple_model_fitting <- function(volumes, sr, snp_data, tmp, snp){
  print(snp)
  for (id in snp_data$IID){
    tmp$SNP[tmp$RID==id] <- risk[,snp][snp_data$IID==id]
    tmp$h_subregion[tmp$RID==id] <- volumes[,sr][volumes$Subject==id]
  }
  
  tmp <- tmp[tmp$RID %in% which(!is.na(snp_data[,snp])),]
  fit <- lme(h_subregion~SNP+AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol+VISCODE2, random = ~1|RID, data = tmp, control = lmeControl(opt = "optim"))

  return(list(summary(fit)))
}

fit <- vector(mode = "list", length = 0)
for (snp in (colnames(risk)[7:ncol(risk)])){
  model <- simple_model_fitting(volumes, 'hippocampal_tail', risk, tmp, snp)
  fit <- c(fit, model)
}

##### Stratified Model #####
tmp$DX.bl <- NA
stratified_model_fitting <- function(volumes, h_subregion, snp_data, dem_data, group, tmp, snp){
  print(snp)
  for (id in snp_data$IID){
    tmp$SNP[tmp$RID==id] <- risk[,snp][snp_data$IID==id]
    tmp$h_subregion[tmp$RID==id] <- volumes[,sr][volumes$Subject==id]
    tmp$DX.bl[tmp$RID == id] <- unique(dem_data$DX.bl[dem_data$RID == id])
  }
  
  
  fit <- lme(h_subregion~SNP, random = ~1|ID, data = tmp_data[tmp_data$DX.bl==group,])
  return(list(summary(fit)))
}
fit <- stratified_model_fitting(volumes, 'hippocampal_tail', risk, data2, 'AD')
