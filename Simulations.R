#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%                 Haplotype Association Analysis            %%%%%%%%%%%
#%%%%%%%%                    Brunilda-Sophia Balliu                 %%%%%%%%%%%
#%%%%%%%%                     Leiden January 2013                   %%%%%%%%%%%
#%%%%%%%%          Functions to perform simulation study            %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




args<-commandArgs(TRUE)
Ncases = as.numeric(args[1])
k<- as.numeric(args[2])
i<-as.numeric(args[3])
LD= as.numeric(args[4])
PHASED=as.numeric(args[5])
STAND=as.numeric(args[6])
Ncontrols=3000
################################################################################
########             Packages
################################################################################

library(haplo.stats)


################################################################################
########             Functions needed
################################################################################
setwd('~/Dropbox/03_Haplotype_I/01_Rscripts/')
source("HierarchHaploLikeli.R")
source('RgenericAll.R');
source('Rgenetics.R');
source('RcomputeResources.R') ; 

#mod = activateModule('mod');
################################################################################
########             Simulating the data 
################################################################################

# True parameter values  
theta3.cases=structure(c(seq(0.05,0.5,length.out=k) ,  rep(.01, ncol(combn(k,2))), rep(0,ncol(combn(k,3))+ncol(combn(k,4)))),.Names = deltaNames(k))
#theta3.cases=structure(c(rep(0.5,times=k),rep(0,2^k-k-1)),.Names = deltaNames(k))
theta1.cases=theta321(theta3.cases,k)$theta1 

theta3.controls=structure(c(seq(0.05,0.5,length.out=k) ,  rep(.01, ncol(combn(k,2))), rep(0,ncol(combn(k,3))+ncol(combn(k,4)))),.Names = deltaNames(k))
#theta3.controls=structure(c(rep(0.5,times=k),rep(0,2^k-k-1)),.Names = deltaNames(k))
theta1.controls=theta321(theta3.controls,k)$theta1 

# Use true parameters to simulate data form a multinomial distribution
nr11=rmultinom(Ncases,size=1,prob=theta1.cases)==1
D=matrix(rep(rownames(nr11),Ncases),byrow=F,nrow=Ncases)[nr11]
first_allele_cases=t(sapply(1:length(D),function(i)   as.numeric(unlist(strsplit(D[i],"")))))

nr12=rmultinom(Ncases,size=1,prob=theta1.cases)==1
D=matrix(rep(rownames(nr12),Ncases),byrow=F,nrow=Ncases)[nr12]
second_allele_cases=t(sapply(1:length(D),function(i)   as.numeric(unlist(strsplit(D[i],"")))))

nr21=rmultinom(Ncontrols,size=1,prob=theta1.controls)==1
D=matrix(rep(rownames(nr21),Ncontrols),byrow=F,nrow=Ncontrols)[nr21]
first_allele_controls=t(sapply(1:length(D),function(i)   as.numeric(unlist(strsplit(D[i],"")))))

nr22=rmultinom(Ncontrols,size=1,prob=theta1.controls)==1
D=matrix(rep(rownames(nr22),Ncontrols),byrow=F,nrow=Ncontrols)[nr22]
second_allele_controls=t(sapply(1:length(D),function(i)   as.numeric(unlist(strsplit(D[i],"")))))

# Check: Are the frequencies in the data close to the real values? 
paternal=apply(first_allele_cases[ , 1:k] , 1 , paste , collapse = "" )
maternal=apply(second_allele_cases[ , 1:k] , 1 , paste , collapse = "")
observed=(table(c(paternal,maternal))/(2*Ncases))[eta1Names(k)]
expected=theta1.cases
p1=sum(round((observed-expected),1))

# Save the data
counts=cbind(apply(nr11,1,sum),apply(nr12,1,sum),apply(nr21,1,sum),apply(nr22,1,sum)) ; colnames(counts)=c("nr11","nr12","nr21","nr22")
#for(j in 1:dim(counts)[2]) write.table(t(counts[,j]),file=paste(colnames(counts)[j],"txt",sep="."),quote = FALSE, sep = ";", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE)

# Construct genotype data from haplotype data
gts.cases=first_allele_cases+second_allele_cases
gts.controls=first_allele_controls+second_allele_controls


# Estimate haplotype frequencies from haplotype data
if(k==3){
  geno.cases=1+cbind(first_allele_cases[,1], second_allele_cases[,1],first_allele_cases[,2], second_allele_cases[,2],first_allele_cases[,3], second_allele_cases[,3])
  geno.controls=1+cbind(first_allele_controls[,1], second_allele_controls[,1],first_allele_controls[,2], second_allele_controls[,2],first_allele_controls[,3], second_allele_controls[,3])
}

if(k==4){
  geno.cases=1+cbind(first_allele_cases[,1], second_allele_cases[,1],first_allele_cases[,2], second_allele_cases[,2],first_allele_cases[,3], second_allele_cases[,3],first_allele_cases[,4], second_allele_cases[,4])
  geno.controls=1+cbind(first_allele_controls[,1], second_allele_controls[,1],first_allele_controls[,2], second_allele_controls[,2],first_allele_controls[,3], second_allele_controls[,3],first_allele_controls[,4], second_allele_controls[,4])
}

geno=rbind(geno.cases,geno.controls)

estimated=haplo.em(geno)$hap.prob
names(estimated)=names(theta1.cases)
p2=sum(round(expected-estimated,1))



if(k==3){
  ###############################################################################
  ##########               Simultaneous Optimization              ###############
  ###############################################################################
  if(FALSE){
    ptm <- proc.time()[3]
    # Null Hypothesis
    o=k
    op.NULL=OptimizerSimultaneousRaw(ind.pooled=1:o,InitValues=theta123(estimated,k)$theta3, gts.cases=gts.cases,gts.controls=gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.NULL)[2]) write.table(t(op.NULL[,j]),file=paste(k,o,"a",colnames(op.NULL)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    ptm <- proc.time()[3]-ptm
    
    # Full Alternative Hypothesis
    o=0
    op.ALTFULL=OptimizerSimultaneousRaw(ind.pooled=0,InitValues=c(rbind(theta3.cases,theta3.controls)),gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTFULL)[2]) write.table(t(op.ALTFULL[,j]),file=paste(k,o,"a",colnames(op.ALTFULL)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTFULL=2^k-1
    chisqALTFULL=-2*(op.ALTFULL$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTFULL=pchisq(chisqALTFULL,df=dfALTFULL,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 1 : Testing only allele frequencies
    ind.pooled=2:k
    InitValues=AdjustInitialValuesForPooledOptim(k=k,ind.pooled=ind.pooled, par.cases=theta3.cases,par.controls=theta3.controls)
    op.ALTRES1=OptimizerSimultaneousRaw(ind.pooled=ind.pooled,InitValues=InitValues,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTRES1)[2]) write.table(t(op.ALTRES1[,j]),file=paste(k,o,"a",colnames(op.ALTRES1)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTRES1=4
    chisqALTRES1=-2*(op.ALTRES1$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTRES1=pchisq(chisqALTRES1,df=dfALTRES1,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 2 : Testing allele frequencies and 2nd order LD parameters
    ind.pooled=3:k
    InitValues=AdjustInitialValuesForPooledOptim(k=k,ind.pooled=ind.pooled, par.cases=theta3.cases,par.controls=theta3.controls)
    op.ALTRES2=OptimizerSimultaneousRaw(ind.pooled=ind.pooled,InitValues=InitValues,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTRES2)[2]) write.table(t(op.ALTRES2[,j]),file=paste(k,o,"a",colnames(op.ALTRES2)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTRES2=10
    chisqALTRES2=-2*(op.ALTRES2$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTRES2=pchisq(chisqALTRES2,df=dfALTRES2,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 3 : Testing allele frequencies and 2nd and 3rd order LD parameters
    ind.pooled=k
    InitValues=AdjustInitialValuesForPooledOptim(k=k,ind.pooled=ind.pooled, par.cases=theta3.cases,par.controls=theta3.controls)
    op.ALTRES3=OptimizerSimultaneousRaw(ind.pooled=ind.pooled,InitValues=InitValues,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTRES3)[2]) write.table(t(op.ALTRES3[,j]),file=paste(k,o,"b",colnames(op.ALTRES3)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTRES3=14
    chisqALTRES3=-2*(op.ALTRES3$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTRES3=pchisq(chisqALTRES3,df=dfALTRES3,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 4 : Testing only for 4th order LD parameters
    dfALTRES4=1
    chisqALTRES4=-2*(op.ALTFULL$MinusLogLikelihood-op.ALTRES3$MinusLogLikelihood)
    pvalALTRES4=pchisq(chisqALTRES4,df=dfALTRES3,lower.tail = F)
    
    # Restricted Alternative Hypothesis 5 : Testing for 3rd and 4th order LD parameters
    dfALTRES5=5
    chisqALTRES5=-2*(op.ALTFULL$MinusLogLikelihood-op.ALTRES2$MinusLogLikelihood)
    pvalALTRES5=pchisq(chisqALTRES5,df=dfALTRES5,lower.tail = F)
    
    # Restricted Alternative Hypothesis 6 : Testing for 2nd, 3rd and 4th order LD parameters
    dfALTRES6=11
    chisqALTRES6=-2*(op.ALTFULL$MinusLogLikelihood-op.ALTRES1$MinusLogLikelihood)
    pvalALTRES6=pchisq(chisqALTRES6,df=dfALTRES6,lower.tail = F)
    
    
    write.table(t(c(i,chisqALTFULL, chisqALTRES1 , chisqALTRES2, chisqALTRES3, chisqALTRES4,chisqSING)),file=paste('chisq',i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,pvalALTFULL,pvalALTRES1,pvalALTRES2,pvalALTRES3,pvalALTRES4,pvalSING)),file=paste('pvals',i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,(proc.time()-ptm)[3])),file="time.txt",quote = FALSE, sep = ";", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE)
    
    ptm <- proc.time()[3]-ptm
  }
  
  ###############################################################################
  ##########               Hierarchical Optimization              ###############
  ###############################################################################
  if(TRUE){
    ptm <- proc.time()[3]
    # Null Hypothesis
    ind.pooled=1:k
    op.NULL=OptimizerHierarchical.theta2(ind.pooled,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
  
  
    # Full Alternative Hypothesis
    ind.pooled=0
    op.ALTFULL=OptimizerHierarchical.theta2(ind.pooled,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    dfALTFULL=2^k-1
    chisqALTFULL=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.NULL[2^k-1,6])
    pvalALTFULL=pchisq(chisqALTFULL,df=dfALTFULL,lower.tail = F)
     
    # Single SNP Hypothesis
    dfSING=1
    chisqSING=-2*(op.ALTFULL[1:k,5]+op.ALTFULL[1:k,6] - op.NULL[1:k,6])
    pvalSING=pchisq(chisqSING,df=rep(dfSING,k),lower.tail = F)
  

    # Restricted Alternative Hypothesis 1 : Testing allele frequencies
    ind.pooled=2:k
    op.ALTRES1=OptimizerHierarchical.theta2(ind.pooled,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    dfALTRES1=k
    chisqALTRES1=-2*(sum(op.ALTRES1[2^k-1,6])-op.NULL[2^k-1,6])
    pvalALTRES1=pchisq(chisqALTRES1,df=dfALTRES1,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 2 : Testing allele frequencies and 2nd order LD parameters
    ind.pooled=3:k
    op.ALTRES2=OptimizerHierarchical.theta2.RestHypothesis(ind.pooled=ind.pooled,ParAlt=op.ALTFULL,gts.cases=gts.cases,gts.controls=gts.controls,PHASED=PHASED,STAND=STAND)    
    dfALTRES2=k+dim(combn(k,2))[2]
    chisqALTRES2=-2*(sum(op.ALTRES2[2^k-1,6])-op.NULL[2^k-1,6])
    pvalALTRES2=pchisq(chisqALTRES2,df=dfALTRES2,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 3 : Testing 3rd order LD parameters
    dfALTRES5=1
    chisqALTRES5=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.ALTRES2[2^k-1,6])
    pvalALTRES5=pchisq(chisqALTRES5,df=dfALTRES5,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 4 : Testing 2nd and 3rd order LD parameters
    dfALTRES6=1+dim(combn(k,2))[2]
    chisqALTRES6=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.ALTRES1[2^k-1,6])
    pvalALTRES6=pchisq(chisqALTRES6,df=dfALTRES6,lower.tail = F)
  
  
  Time <- proc.time()[3]-ptm
  
  theta2.cases=theta123(theta1.cases,k)$theta2
  theta2.controls=theta123(theta1.controls,k)$theta2
  
  write.table(t(c(i,op.NULL[,1]-theta2.cases)),file=paste(i,'BiasNull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,c(op.ALTFULL[,1:2]-cbind(theta2.cases,theta2.controls)))),file=paste(i,'BiasAltFull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,c((op.ALTRES1[,1:2]-cbind(theta2.cases,theta2.controls))[1:3,],op.ALTRES1[-c(1:3),1]-theta2.cases[-c(1:3)]))),file=paste(i,'BiasAltRes1.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,c((op.ALTRES2[,1:2]-cbind(theta2.cases,theta2.controls))[1:6,],op.ALTRES2[-c(1:6),1]-theta2.cases[-c(1:6)]))),file=paste(i,'BiasAltRes2.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  
  write.table(t(c(i,op.NULL[,3])),file=paste(i,'SENull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,c(op.ALTFULL[,3:4]))),file=paste(i,'SEAltFull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,c((op.ALTRES1[1:3,3:4])),op.ALTRES1[-c(1:3),3])),file=paste(i,'SEAltRes1.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,c(op.ALTRES2[1:6,3:4]),op.ALTRES2[-c(1:6),3])),file=paste(i,'SEAltRes2.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  
  write.table(t(c(i,chisqALTFULL, chisqALTRES1 , chisqALTRES2, chisqALTRES5, chisqALTRES6,chisqSING)),file=paste(i,'ChiSqStat.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(t(c(i,pvalALTFULL,pvalALTRES1,pvalALTRES2, pvalALTRES5,pvalALTRES6,pvalSING)),file=paste(i,'Pvals.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(Time,file="Time.txt",quote = FALSE, sep = ";", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE)
}

}

if(k==4){
  ###############################################################################
  ##########               Simultaneous Optimization              ###############
  ###############################################################################
  if(FALSE){
    ptm <- proc.time()[3]
    # Null Hypothesis
    o=k
    op.NULL=OptimizerSimultaneousRaw(ind.pooled=1:o,InitValues=theta123(estimated,k)$theta3, gts.cases=gts.cases,gts.controls=gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.NULL)[2]) write.table(t(op.NULL[,j]),file=paste(k,o,"a",colnames(op.NULL)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    ptm <- proc.time()[3]-ptm
    
    # Full Alternative Hypothesis
    o=0
    op.ALTFULL=OptimizerSimultaneousRaw(ind.pooled=0,InitValues=c(rbind(theta3.cases,theta3.controls)),gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTFULL)[2]) write.table(t(op.ALTFULL[,j]),file=paste(k,o,"a",colnames(op.ALTFULL)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTFULL=2^k-1
    chisqALTFULL=-2*(op.ALTFULL$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTFULL=pchisq(chisqALTFULL,df=dfALTFULL,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 1 : Testing only allele frequencies
    ind.pooled=2:k
    InitValues=AdjustInitialValuesForPooledOptim(k=k,ind.pooled=ind.pooled, par.cases=theta3.cases,par.controls=theta3.controls)
    op.ALTRES1=OptimizerSimultaneousRaw(ind.pooled=ind.pooled,InitValues=InitValues,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTRES1)[2]) write.table(t(op.ALTRES1[,j]),file=paste(k,o,"a",colnames(op.ALTRES1)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTRES1=4
    chisqALTRES1=-2*(op.ALTRES1$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTRES1=pchisq(chisqALTRES1,df=dfALTRES1,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 2 : Testing allele frequencies and 2nd order LD parameters
    ind.pooled=3:k
    InitValues=AdjustInitialValuesForPooledOptim(k=k,ind.pooled=ind.pooled, par.cases=theta3.cases,par.controls=theta3.controls)
    op.ALTRES2=OptimizerSimultaneousRaw(ind.pooled=ind.pooled,InitValues=InitValues,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTRES2)[2]) write.table(t(op.ALTRES2[,j]),file=paste(k,o,"a",colnames(op.ALTRES2)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTRES2=10
    chisqALTRES2=-2*(op.ALTRES2$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTRES2=pchisq(chisqALTRES2,df=dfALTRES2,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 3 : Testing allele frequencies and 2nd and 3rd order LD parameters
    ind.pooled=k
    InitValues=AdjustInitialValuesForPooledOptim(k=k,ind.pooled=ind.pooled, par.cases=theta3.cases,par.controls=theta3.controls)
    op.ALTRES3=OptimizerSimultaneousRaw(ind.pooled=ind.pooled,InitValues=InitValues,gts.cases,gts.controls,PHASED=PHASED,STAND=STAND)    
    #for(j in 1:dim(op.ALTRES3)[2]) write.table(t(op.ALTRES3[,j]),file=paste(k,o,"b",colnames(op.ALTRES3)[j],i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    dfALTRES3=14
    chisqALTRES3=-2*(op.ALTRES3$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
    pvalALTRES3=pchisq(chisqALTRES3,df=dfALTRES3,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 4 : Testing only for 4th order LD parameters
    dfALTRES4=1
    chisqALTRES4=-2*(op.ALTFULL$MinusLogLikelihood-op.ALTRES3$MinusLogLikelihood)
    pvalALTRES4=pchisq(chisqALTRES4,df=dfALTRES3,lower.tail = F)
    
    # Restricted Alternative Hypothesis 5 : Testing for 3rd and 4th order LD parameters
    dfALTRES5=5
    chisqALTRES5=-2*(op.ALTFULL$MinusLogLikelihood-op.ALTRES2$MinusLogLikelihood)
    pvalALTRES5=pchisq(chisqALTRES5,df=dfALTRES5,lower.tail = F)
    
    # Restricted Alternative Hypothesis 6 : Testing for 2nd, 3rd and 4th order LD parameters
    dfALTRES6=11
    chisqALTRES6=-2*(op.ALTFULL$MinusLogLikelihood-op.ALTRES1$MinusLogLikelihood)
    pvalALTRES6=pchisq(chisqALTRES6,df=dfALTRES6,lower.tail = F)
    
    
    write.table(t(c(i,chisqALTFULL, chisqALTRES1 , chisqALTRES2, chisqALTRES3, chisqALTRES4,chisqSING)),file=paste('chisq',i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,pvalALTFULL,pvalALTRES1,pvalALTRES2,pvalALTRES3,pvalALTRES4,pvalSING)),file=paste('pvals',i,"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,(proc.time()-ptm)[3])),file="time.txt",quote = FALSE, sep = ";", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE)
    
    ptm <- proc.time()[3]-ptm
  }
  
  ###############################################################################
  ##########               Hierarchical Optimization              ###############
  ###############################################################################
  if(TRUE){
    ptm <- proc.time()[3]
    # Null Hypothesis
    ind.pooled=1:k
    op.NULL=OptimizerHierarchical.theta2(ind.pooled,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED=PHASED,STAND=STAND)    
    
    # Full Alternative Hypothesis
    ind.pooled=0
    op.ALTFULL=OptimizerHierarchical.theta2(ind.pooled,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED=PHASED,STAND=STAND)    
    dfALTFULL=2^k-1
    chisqALTFULL=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.NULL[2^k-1,6])
    pvalALTFULL=pchisq(chisqALTFULL,df=dfALTFULL,lower.tail = F)
    
    # Single SNP Hypothesis
    dfSING=1
    chisqSING=-2*(op.ALTFULL[1:k,5]+op.ALTFULL[1:k,6] - op.NULL[1:k,6])
    pvalSING=pchisq(chisqSING,df=rep(dfSING,k),lower.tail = F)
    
    # Restricted Alternative Hypothesis 1 : Testing allele frequencies
    # Ho: allele frequencies pooled, all else free
    # H1: allele frequencies free, all else free
    ind.pooled=1
    op.ALTRES1=OptimizerHierarchical.theta2.RestHypothesis(ind.pooled=ind.pooled,ParNull=op.NULL,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED=PHASED,STAND=STAND)    
    dfALTRES1=k
    chisqALTRES1=-2*(sum(op.ALTFULL[2^k-1,5:6])-sum(op.ALTRES1[2^k-1,5:6]))
    pvalALTRES1=pchisq(chisqALTRES1,df=dfALTRES1,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 2 : Testing allele frequencies and 2nd order LD parameters
    # Ho: allele frequencies and 2nd order LD pooled, all else free
    # H1: allele frequencies and 2nd order LD free, all else free
    ind.pooled=1:2
    op.ALTRES2=OptimizerHierarchical.theta2.RestHypothesis(ind.pooled=ind.pooled,ParNull=op.NULL,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED=PHASED,STAND=STAND)    
    dfALTRES2=k+dim(combn(k,2))[2]
    chisqALTRES2=-2*(sum(op.ALTFULL[2^k-1,5:6])-sum(op.ALTRES2[2^k-1,5:6]))
    pvalALTRES2=pchisq(chisqALTRES2,df=dfALTRES2,lower.tail = F)
    
    
    # Restricted Alternative Hypothesis 3 : Testing allele frequencies, 2nd and 3rd order LD parameters
    # Ho: allele frequencies, 2nd and 3rd order LD pooled, all else free
    # H1: allele frequencies, 2nd and 3nd order LD free, all else free
    ind.pooled=1:3
    op.ALTRES3=OptimizerHierarchical.theta2.RestHypothesis(ind.pooled=ind.pooled,ParNull=op.NULL,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED=PHASED,STAND=STAND)    
    dfALTRES3=k+dim(combn(k,2))[2]+dim(combn(k,3))[2]
    chisqALTRES3=-2*(sum(op.ALTFULL[2^k-1,5:6])-sum(op.ALTRES3[2^k-1,5:6]))
    pvalALTRES3=pchisq(chisqALTRES3,df=dfALTRES3,lower.tail = F)
    
   if(FALSE) {
     # Restricted Alternative Hypothesis 4 : Testing 4th order LD parameters
     dfALTRES4=1
     chisqALTRES4=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.ALTRES3[2^k-1,6])
     pvalALTRES4=pchisq(chisqALTRES4,df=dfALTRES4,lower.tail = F)
     
     # Restricted Alternative Hypothesis 5 : Testing 3rd and 4th order LD parameters
     dfALTRES5=1+dim(combn(k,3))[2]
     chisqALTRES5=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.ALTRES2[2^k-1,6])
     pvalALTRES5=pchisq(chisqALTRES5,df=dfALTRES5,lower.tail = F)
     
     
     # Restricted Alternative Hypothesis 6 : Testing 2nd, 3rd and 4th order LD parameters
     dfALTRES6=1+dim(combn(k,2))[2]+dim(combn(k,3))[2]
     chisqALTRES6=-2*(sum(op.ALTFULL[2^k-1,5:6])-op.ALTRES1[2^k-1,6])
     pvalALTRES6=pchisq(chisqALTRES6,df=dfALTRES6,lower.tail = F)
   }    
    Time <- proc.time()[3]-ptm
    
    
    theta2.cases=theta123(theta1.cases,k)$theta2
    theta2.controls=theta123(theta1.controls,k)$theta2
    
    write.table(t(c(i,op.NULL[,1]-theta2.cases)),file=paste(i,'BiasNull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c(op.ALTFULL[,1:2]-cbind(theta2.cases,theta2.controls)))),file=paste(i,'BiasAltFull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c((op.ALTRES1[,1:2]-cbind(theta2.cases,theta2.controls))[1:4,],op.ALTRES1[-c(1:4),1]-theta2.cases[-c(1:4)]))),file=paste(i,'BiasAltRes1.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c((op.ALTRES2[,1:2]-cbind(theta2.cases,theta2.controls))[1:10,],op.ALTRES2[-c(1:10),1]-theta2.cases[-c(1:10)]))),file=paste(i,'BiasAltRes2.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c((op.ALTRES3[,1:2]-cbind(theta2.cases,theta2.controls))[1:14,],op.ALTRES3[-c(1:14),1]-theta2.cases[-c(1:14)]))),file=paste(i,'BiasAltRes3.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    
    write.table(t(c(i,op.NULL[,3])),file=paste(i,'SENull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c(op.ALTFULL[,3:4]))),file=paste(i,'SEAltFull.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c((op.ALTRES1[1:4,3:4])),op.ALTRES1[-c(1:4),3])),file=paste(i,'SEAltRes1.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c(op.ALTRES2[1:10,3:4]),op.ALTRES2[-c(1:10),3])),file=paste(i,'SEAltRes2.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,c(op.ALTRES3[1:14,3:4]),op.ALTRES3[-c(1:14),3])),file=paste(i,'SEAltRes3.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    
    write.table(t(c(i,chisqALTFULL, chisqALTRES1 , chisqALTRES2, chisqALTRES3,chisqSING)),file=paste(i,'ChiSqStat.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(t(c(i,pvalALTFULL,pvalALTRES1,pvalALTRES2,pvalALTRES3,pvalSING)),file=paste(i,'Pvals.txt',sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
    write.table(Time,file="Time.txt",quote = FALSE, sep = ";", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE)
  }
  
}
