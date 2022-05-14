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

#################################### MODELS ####################################
# install.packages('lme4')
# library(lme4)
# install.packages('nlme')
library(nlme)

tmp_data <- data.frame(volumes$Subject[volumes$Subject %in% risk$IID] ,volumes$hippocampal_tail[volumes$Subject %in% risk$IID])
colnames(tmp_data) <- c('ID', 'hippocampal_tail')
tmp_data$SNP <- NA

for (i in risk$IID){
  tmp_data$SNP[tmp_data$ID == i] <- risk$rs4575098_A[risk$IID == i]
}

fit <- lme(hippocampal_tail~SNP, random = ~(1|ID), data = tmp_data)
