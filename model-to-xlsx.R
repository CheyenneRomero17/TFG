
# Selection of data

volume = names(volumes[,2:ncol(volumes)])

##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP risk with APOE4
###
#################################

snps = names(risk[,7:ncol(risk)])

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

library(writexl)
write_xlsx(coefs,file="risk_apoe.xlsx")

##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP risk without APOE4
###
#################################

snps = names(risk[,7:ncol(risk)])

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

library(writexl)
write_xlsx(coefs,file="risk.xlsx")

##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP proxy with APOE4
###
#################################

snps = names(proxy[,7:ncol(proxy)])

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

library(writexl)
write_xlsx(coefs,file="proxy_apoe.xlsx")

##################################
###
### Longitudinal Mixed Models and Diagnostic plots
### Main longitudinal effect SNP proxy without APOE4
###
#################################

snps = names(proxy[,7:ncol(proxy)])

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

library(writexl)
write_xlsx(coefs,file="proxy.xlsx")