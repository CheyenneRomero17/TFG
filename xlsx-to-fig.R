library(corrplot)
library(readxl)

corr <- function(tis, path, strata = FALSE, group = NULL){
  #Create matrices
  tisouts=unique(tis$Outcome)
  tisdets=unique(tis$Determinant)
  
  mat=as.data.frame(matrix(NA,nrow=length(tisdets),ncol=length(tisouts)))
  rownames(mat)=tisdets
  colnames(mat)=tisouts
  
  matb=mat
  matt=mat
  matp=mat
  for(o in tisouts){for(d in tisdets){
  matb[d,o]=tis[tis$Determinant==d & tis$Outcome==o,"Beta"]
  matt[d,o]=tis[tis$Determinant==d & tis$Outcome==o,"T"]
  matp[d,o]=tis[tis$Determinant==d & tis$Outcome==o,"P"]
  }}
  
  
  #Get most significant variants
  selO=grep("t |_",tisouts[rev(order(tisouts))],value=T,invert=T)
  
  best=as.data.frame(matrix(NA,nrow=1,ncol=2))
  c=1
  for(d in grep("^rs",tisdets,value=T)){
  best[c,1]=d
  best[c,2]=min(tis[tis$Determinant==d & tis$Outcome %in% selO,"P"])
  c=c+1
  }
  
  selD=best[order(best$V2),1]
  
  
  col=colorRampPalette(c("red", "white", "blue"))(20)
  
  matsel=as.matrix(matt)
  maxsel=max(abs(range(matsel)))
  psel=as.matrix(matp)
  
  if(strata){
    png(height=8, width=17, file=paste("figs/",path,"/SNP_strata_",group,".png",sep=''),units="in",res=500)
  }else{
    png(height=8, width=17, file=paste("figs/",path,"/SNP_",i,".png",sep=''),units="in",res=500)
  }
  
  #allow plotting (text) outside window!!
  par(xpd=NA)
  
  corrplot(matsel,
  mar=c(2,2,1,2)+.1,
  is.corr=FALSE,
  tl.col="black", tl.srt=60,
  cl.lim=c(-6,6),
  cl.align.text = 'l',
  method="color", tl.cex=1
  )
  
  
  p_format <- function(x)
  {ifelse(x<0.0001,"***",ifelse(x<0.003,"**", ifelse(x<0.01,"*", "")))
  }
  
  pos=expand.grid(nrow(psel):1,1:ncol(psel))
  pos=cbind(pos[,2],pos[,1])
  text(pos, p_format(psel), col ='white')
  
  dev.off()
}

## Make a pretty QQ plot of p-values
qq = function(pvector, ci.alpha, ...) {
	if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
	pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
	o = -log10(sort(pvector,decreasing=F))
	o[which(o > 12)]=12
	e = -log10( ppoints(length(pvector) ))

	lambda<-signif(qchisq(median(pvector),1,lower.tail=F)/qchisq(.5,1),digits=4);
	plot(-1,-1, xlab=paste("expected -logP","  (Lambda = ", lambda, ")", sep=''), ylab=paste("observed -logP",sep=''), xlim=c(0,max(o,e)), ylim=c(0,max(o,e)),...)
	
	cil=-log10(qbeta(1-(ci.alpha/2), 1:length(pvector), length(pvector)-1:length(pvector)))
	ciu=-log10(qbeta(ci.alpha/2, 1:length(pvector), length(pvector)-1:length(pvector)))
	polygon(c(e,rev(e)),c(cil,rev(ciu)),col=adjustcolor("red",alpha.f=0.5),border=F)

	abline(0,1,col="red")

	points(e,o,pch=19,cex=1)

}

##########
## Main ##
##########

files <- c('risk','proxy','risk_apoe','proxy_apoe','risk_interact',
           'proxy_interact','risk_strata','proxy_strata')

### Correlation Heatmap ###

for (i in files){
  if (i %in% c('risk_strata','proxy_strata')){
    coefs <- read_xlsx(paste("~/Uni/MatCAD/4t/TFG/data/",i,".xlsx", sep = ''))
    if (i=='proxy_strata'){
      coefs <- coefs[-(which(coefs$Determinant == 'rs114360492_T')),]
      coefs <- coefs[-(which(coefs$Determinant == 'rs184384746_T')),]
    }
    tis=coefs
    for (j in c('AD', 'MCI', 'SMC')){
      tis=coefs[coefs$Type==j,]
      corr(tis, paste(i,'/',j,sep=''), TRUE, j)
    }
  }else{
    coefs <- read_xlsx(paste("~/Uni/MatCAD/4t/TFG/data/",i,".xlsx", sep = ''))
    tis=coefs
    corr(tis, i)
  }
}

### QQ plots ###

for (i in files){
  if (i %in% c('risk_strata','proxy_strata')){
    coefs <- read_xlsx(paste("~/Uni/MatCAD/4t/TFG/data/",i,".xlsx", sep = ''))
    if (i=='proxy_strata'){
      coefs <- coefs[-(which(coefs$Determinant == 'rs114360492_T')),]
      coefs <- coefs[-(which(coefs$Determinant == 'rs184384746_T')),]
    }
    tis=coefs
    tisouts=unique(tis$Outcome)
    tisdets=unique(tis$Determinant)
    for (j in c('AD', 'MCI', 'SMC')){
      tis=coefs[coefs$Type==j,]

      singlevar=grep("^rs",tisdets,value=T)
      for(o in tisouts){
        p=tis[tis$Outcome==o & tis$Determinant %in% singlevar,"P"]
        
        #Make QQ Plot
        png(paste("figs/",i,"/",j,"/MSxTissue_",o,"_qq.png", sep=''),w=2,h=2,units="in",res=1500)
        par(oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.5,0.5,0),cex=0.6)
        qq(p$P,0.05,main=o)
        dev.off()
      }
    }
  }else{
    coefs <- read_xlsx(paste("~/Uni/MatCAD/4t/TFG/data/",i,".xlsx", sep = ''))
    tis=coefs
    tisouts=unique(tis$Outcome)
    tisdets=unique(tis$Determinant)
    
    singlevar=grep("^rs",tisdets,value=T)
    for(o in tisouts){
      p=tis[tis$Outcome==o & tis$Determinant %in% singlevar,"P"]
      
      #Make QQ Plot
      png(paste("figs/",i,"/MSxTissue_",o,"_qq.png", sep=''),w=2,h=2,units="in",res=1500)
      par(oma=c(0,0,0,0),mar=c(2.5,2.5,2.5,1),mgp=c(1.5,0.5,0),cex=0.6)
      qq(p$P,0.05,main=o)
      dev.off()
    }
  }
}


