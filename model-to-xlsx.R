library(nlme)
library(writexl)

# Selection of data

volume = names(volumes[,2:ncol(volumes)])


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP risk with APOE4
###
#################################

snps = names(risk[,7:ncol(risk)])
data <- data2[data2$RID %in% risk$IID,c(1,2,7,8,9,10,11,12)]

coefs=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs)=c("DateTime","Outcome","Determinant","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP+AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol+APOE4"
    
    # Gathering the needed data
    data <- data_gather(data2, risk, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]=ff
    coefs[count,5]=3230
    coefs[count,6:8]=res[2,c(1,2,4)]
    coefs[count,9]=as.numeric(format(res[2,5], format="e", digits=2))
    count=count+1
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,10,7,8,9,6)]

write_xlsx(coefs, path="risk_apoe.xlsx")


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP risk without APOE4
###
#################################

snps = names(risk[,7:ncol(risk)])
data <- data2[data2$RID %in% risk$IID,c(1,2,7,8,9,10,11,12)]

coefs=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs)=c("DateTime","Outcome","Determinant","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP+AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol"
    
    # Gathering the needed data
    data <- data_gather(data2, risk, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]=ff
    coefs[count,5]=3230
    coefs[count,6:8]=res[2,c(1,2,4)]
    coefs[count,9]=as.numeric(format(res[2,5], format="e", digits=2))
    count=count+1
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,10,7,8,9,6)]

write_xlsx(coefs, path="risk.xlsx")


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP proxy with APOE4
###
#################################

snps = names(proxy[,7:ncol(proxy)])
data <- data2[data2$RID %in% proxy$IID,c(1,2,7,8,9,10,11,12)]

coefs=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs)=c("DateTime","Outcome","Determinant","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP+AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol+APOE4"
    
    # Gathering the needed data
    data <- data_gather(data2, proxy, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]=ff
    coefs[count,5]=3230
    coefs[count,6:8]=res[2,c(1,2,4)]
    coefs[count,9]=as.numeric(format(res[2,5], format="e", digits=2))
    count=count+1
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,10,7,8,9,6)]

write_xlsx(coefs, path="proxy_apoe.xlsx")

##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP proxy without APOE4
###
#################################

snps = names(proxy[,7:ncol(proxy)])
data <- data2[data2$RID %in% proxy$IID,c(1,2,7,8,9,10,11,12)]

coefs=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs)=c("DateTime","Outcome","Determinant","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP+AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol"
    
    # Gathering the needed data
    data <- data_gather(data2, proxy, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]=ff
    coefs[count,5]=3230
    coefs[count,6:8]=res[2,c(1,2,4)]
    coefs[count,9]=as.numeric(format(res[2,5], format="e", digits=2))
    count=count+1
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,10,7,8,9,6)]

write_xlsx(coefs, path="proxy.xlsx")


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP risk stratified by diagnosis
###
#################################

snps = names(risk[,7:ncol(risk)])
data <- data2[data2$RID %in% risk$IID,c(1,2,7,8,9,10,11,12)]
data$DX.bl[data$DX.bl == 'CN'] <- 'ACN'

coefs=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs)=c("DateTime","Outcome","Determinant","Type","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP*DX.bl+AGE+PTGENDER+PTEDUCAT+IntraCranialVol"
    
    # Gathering the needed data
    data <- data_gather(data2, risk, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    # AD
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]='AD'
    coefs[count,5]=ff
    coefs[count,6]=sum(data$DX.bl=='AD')
    coefs[count,7:9]=res[10,c(1,2,4)]
    coefs[count,10]=as.numeric(format(res[10,5], format="e", digits=2))
    # MCI
    coefs[count+1,1]=gsub(" ","_",Sys.time())
    coefs[count+1,2]=volume[j]
    coefs[count+1,3]=snps[i]
    coefs[count+1,4]='MCI'
    coefs[count+1,5]=ff
    coefs[count+1,6]=sum(data$DX.bl=='MCI')
    coefs[count+1,7:9]=res[11,c(1,2,4)]
    coefs[count+1,10]=as.numeric(format(res[11,5], format="e", digits=2))
    # SMC
    coefs[count+2,1]=gsub(" ","_",Sys.time())
    coefs[count+2,2]=volume[j]
    coefs[count+2,3]=snps[i]
    coefs[count+2,4]='SMC'
    coefs[count+2,5]=ff
    coefs[count+2,6]=sum(data$DX.bl=='SMC')
    coefs[count+2,7:9]=res[12,c(1,2,4)]
    coefs[count+2,10]=as.numeric(format(res[12,5], format="e", digits=2))
    
    count=count+3
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,6,11,8,9,10,7)]

write_xlsx(coefs, path="risk_strata.xlsx")


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP proxy stratified by diagnosis
###
#################################

snps = names(proxy[,7:ncol(proxy)])
data <- data2[data2$RID %in% proxy$IID,c(1,2,7,8,9,10,11,12)]
data$DX.bl[data$DX.bl == 'CN'] <- 'ACN'

coefs=as.data.frame(matrix(NA,nrow=1,ncol=10))
colnames(coefs)=c("DateTime","Outcome","Determinant","Type","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP*DX.bl+AGE+PTGENDER+PTEDUCAT+IntraCranialVol"
    
    # Gathering the needed data
    data <- data_gather(data2, proxy, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    # AD
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]='AD'
    coefs[count,5]=ff
    coefs[count,6]=sum(data$DX.bl=='AD')
    coefs[count,7:9]=res[10,c(1,2,4)]
    coefs[count,10]=as.numeric(format(res[10,5], format="e", digits=2))
    # MCI
    coefs[count+1,1]=gsub(" ","_",Sys.time())
    coefs[count+1,2]=volume[j]
    coefs[count+1,3]=snps[i]
    coefs[count+1,4]='MCI'
    coefs[count+1,5]=ff
    coefs[count+1,6]=sum(data$DX.bl=='MCI')
    coefs[count+1,7:9]=res[11,c(1,2,4)]
    coefs[count+1,10]=as.numeric(format(res[11,5], format="e", digits=2))
    # SMC
    coefs[count+2,1]=gsub(" ","_",Sys.time())
    coefs[count+2,2]=volume[j]
    coefs[count+2,3]=snps[i]
    coefs[count+2,4]='SMC'
    coefs[count+2,5]=ff
    coefs[count+2,6]=sum(data$DX.bl=='SMC')
    coefs[count+2,7:9]=res[12,c(1,2,4)]
    coefs[count+2,10]=as.numeric(format(res[12,5], format="e", digits=2))
    
    count=count+3
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,6,11,8,9,10,7)]

write_xlsx(coefs, path="proxy_strata.xlsx")


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP risk interaction with time
###
#################################

snps = names(risk[,7:ncol(risk)])
data <- data2[data2$RID %in% risk$IID,c(1,2,7,8,9,10,11,12)]

coefs=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs)=c("DateTime","Outcome","Determinant","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP*AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol"
    
    # Gathering the needed data
    data <- data_gather(data2, risk, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]=ff
    coefs[count,5]=3230
    coefs[count,6:8]=res[2,c(1,2,4)]
    coefs[count,9]=as.numeric(format(res[2,5], format="e", digits=2))
    count=count+1
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,10,7,8,9,6)]

write_xlsx(coefs, path="risk_interact.xlsx")


##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP proxy interaction with time
###
#################################

snps = names(proxy[,7:ncol(proxy)])
data <- data2[data2$RID %in% proxy$IID,c(1,2,7,8,9,10,11,12)]

coefs=as.data.frame(matrix(NA,nrow=1,ncol=9))
colnames(coefs)=c("DateTime","Outcome","Determinant","Model","N","BetaU","SE","T","P")
count=1
for(j in 1:length(volume)) {
  print(volume[j])
  results=NULL
  for(i in 1:length(snps)) {
    print(snps[i])
    # Generic formula
    ff = "h_subregion~SNP*AGE+PTGENDER+DX.bl+PTEDUCAT+IntraCranialVol"
    
    # Gathering the needed data
    data <- data_gather(data2, proxy, volumes, snps[i], volume[j])
    
    # Fitting Linear Mixed model
    fit <- lme(fixed=as.formula(ff), random=list(RID=pdDiag(form=~AGE)),data=data, method="ML", control = lmeControl(opt="optim"))
    
    #Save analysis details and coefficients
    s = summary(fit)
    res = s$tTable
    
    #Save analysis details and coefficients
    coefs[count,1]=gsub(" ","_",Sys.time())
    coefs[count,2]=volume[j]
    coefs[count,3]=snps[i]
    coefs[count,4]=ff
    coefs[count,5]=3230
    coefs[count,6:8]=res[2,c(1,2,4)]
    coefs[count,9]=as.numeric(format(res[2,5], format="e", digits=2))
    count=count+1
  }
}

coefs$Beta = coefs$BetaU/coefs$SE
coefs = coefs[,c(1,2,3,4,5,10,7,8,9,6)]

write_xlsx(coefs, path="proxy_interact.xlsx")
