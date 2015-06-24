#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%                 Haplotype Association Analysis            %%%%%%%%%%%
#%%%%%%%%                    Brunilda-Sophia Balliu                 %%%%%%%%%%%
#%%%%%%%%                     Leiden January 2013                   %%%%%%%%%%%
#%%%%%%%%  Reparametrization and Likelihood Optimization Functions  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(sets) 
library(partitions)
library(numDeriv)   
library(snpStats)
#source('~/Dropbox/03_Haplotype_I/01_Rscripts/RgenericAll.R');
#source('~/Dropbox/03_Haplotype_I/01_Rscripts/Rgenetics.R');
#source('~/Dropbox/03_Haplotype_I/01_Rscripts/RcomputeResources.R') ; 
#parallelize_setEnable(F);

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Function to compute the possible haplotypes for k SNPs in binary notation
eta1Names=function(k) {  
  sapply(0:(2^k - 1), function(j) paste(rev(ord2bin(j, k)), collapse = '')) 
}

# Function to compute all possible SNP sets from k SNPs 
deltaNames=function(k){
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  return(sapply(1:length(Pk),function(i) paste(unlist(Pk[i]),collapse="")))
}

Hap=function(k) {
  to<-vector("list",k)
  ti<-lapply(to,function(x) c(0,1))
  h<-expand.grid(ti)
  return(h)
}

# Function to compute coefficients needed for the higher order cumulants.
Coef=function(v)(-1)^(v-1)*factorial(v-1)      

split.f=function(a) return(split(seq_along(a),a))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#Theta 1 --> Theta 2 --> Theta 3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Function to go from haplotype frequencies to LD parameters
theta123=function(theta1,k) {
  theta2=theta3=rep(NA,(2^k)-1)
  names(theta3)=names(theta2)=deltaNames(k)
  pos.hap=Hap(k) ; pos.hap=pos.hap[order(apply(pos.hap,1,sum),decreasing=T),]
  #theta1-->theta2	 
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  for(j in 1:(2^k-1)) {
    p=unlist(Pk[j])
    c2=sapply(1:2^k,function(i) all(pos.hap[i,p]==1)) 
    bp=pos.hap[c2,]
    bp=sapply(1:dim(bp)[1],function(i) paste(bp[i,],collapse=""))                 
    bp=apply(sapply(1:length(bp), function(i) bp[i]==names(theta1)),1,any)
    theta2[names(theta2)==paste(p,collapse="")]=sum(theta1[bp])
  }
  
  #theta2-->theta3	 
  theta3[1:k]=theta2[1:k] 
  for(j in (k+1):(2^k-1)) {
    ind=as.numeric(unlist(strsplit(names(theta2[j]),split="")))
    partition=t(setparts(length(ind)))
    theta3[names(theta3)==paste(ind,collapse="")]= sum(sapply (1:dim(partition)[1],function(i){
      a=partition[i,]; s=split(seq_along(a),a); l.s=length(s); coeff=Coef(l.s)
      eta=sapply(1:l.s ,function(m) {
        index=ind[unlist(s[m])]
        theta2[names(theta2)==paste(index,collapse="")]
      })
      return(coeff*prod(eta)) })) 
  }
  return(list(theta2=theta2,theta3=theta3))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#Theta 3 --> Theta 2 --> Theta 1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Function to go from LD parameters to haplotype frequencies
theta321=function(theta3,k) {
  # Defining initial values for parametes 
  theta2=rep(NA,(2^k)-1);theta2[1:k]=theta3[1:k] 
  theta1=rep(NA,2^k)
  
  #Names of parameters in theta3 
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  names(theta3)=names(theta2)=sapply(1:length(Pk),function(i) paste(unlist(Pk[i]),collapse=""))
  pos.hap=Hap(k) ; pos.hap=pos.hap[order(apply(pos.hap,1,sum),decreasing=T),]
  names(theta1)=sapply(1:(2^k),function(i) paste(pos.hap[i,],collapse=""))
  
  #theta3-->theta2	
  for(j in (k+1):(2^k-1)) {
    ind=as.numeric(unlist(strsplit(names(theta2[j]),split="")))
    partition=t(setparts(length(ind)))
    d=sum(sapply (2:dim(partition)[1],function(i){
      a=partition[i,]; s=split(seq_along(a),a); l.s=length(s); coeff=Coef(l.s)
      eta=sapply(1:l.s ,function(m) {
        index=ind[unlist(s[m])]
                    		return(theta2[names(theta3)==paste(index,collapse="")])
                    		})
         		     return(coeff*prod(eta)) 
                  })) 
          	a=partition[1,]; s=split(seq_along(a),a)
          	theta2[names(theta3)==paste(ind[unlist(s)],collapse="")]=theta3[names(theta3)==paste(ind[unlist(s)],collapse="")]-d
	      }  
    

   #theta2-->theta1	 
   theta1[1]=theta2[2^k-1]  
   for(j in 1:(2^k-2)) {
      	p=unlist(Pk[j])
      	t=apply(pos.hap,1,sum) ; 
      	c1=(t==k-length(p)) ; c2=sapply(1:2^k,function(i) all(pos.hap[i,p]==0)) 
      	pu=pos.hap[c1&c2,]  
      	pl=(1:k)[-p]	       
      	c3=sapply(1:2^k,function(i) all(pos.hap[i,pl]==1)) ; c4=(t>k-length(p))
      	bp=pos.hap[c3&c4,]
      	bp=sapply(1:dim(bp)[1],function(i) paste(bp[i,],collapse=""))
      	bp=apply(sapply(1:length(bp), function(i) bp[i]==names(theta1)),1,any)
      	theta1[names(theta1)==paste(pu,collapse="")]=theta2[names(theta2)==paste(pl,collapse="")] - sum(theta1[bp])
           }
  theta1[2^k]=1-sum(theta1[1:(2^k-1)])
  return(list(theta1=theta1[eta1Names(k)],theta2=theta2))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#Theta 2 --> Theta 1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#Function to go from "marginal" frequencies of '1'-allele haplotypes to haplotype frequencies
theta21=function(theta2,k) {
  #Defining initial values for parametes 
  theta1=rep(NA,2^k) 
  #Names of parameters in theta3 
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  names(theta2)=sapply(1:length(Pk),function(i) paste(unlist(Pk[i]),collapse=""))
  pos.hap=Hap(k) ; pos.hap=pos.hap[order(apply(pos.hap,1,sum),decreasing=T),]
  names(theta1)=sapply(1:(2^k),function(i) paste(pos.hap[i,],collapse=""))  
  #theta2-->theta1	 
  theta1[1]=theta2[2^k-1]  
  for(j in 1:(2^k-2)) {
    p=unlist(Pk[j])
    t=apply(pos.hap,1,sum) ; 
    c1=(t==k-length(p)) ; c2=sapply(1:2^k,function(i) all(pos.hap[i,p]==0)) 
    pu=pos.hap[c1&c2,]  
    pl=(1:k)[-p]	       
    c3=sapply(1:2^k,function(i) all(pos.hap[i,pl]==1)) ; c4=(t>k-length(p))
    bp=pos.hap[c3&c4,]
    bp=sapply(1:dim(bp)[1],function(i) paste(bp[i,],collapse=""))
    bp=apply(sapply(1:length(bp), function(i) bp[i]==names(theta1)),1,any)
    theta1[names(theta1)==paste(pu,collapse="")]=theta2[names(theta2)==paste(pl,collapse="")] - sum(theta1[bp])
  }
  theta1[2^k]=1-sum(theta1[1:(2^k-1)])
  return(theta1[eta1Names(k)])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#D' = D / Dmax   Standardized LD parameters      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


######################
#For pairwise LD
######################
D.prime=function(final.estimates,k) {
  n=names(final.estimates)[-(1:k)]
  final.d.prime=rep(NA,length(n)) ; names(final.d.prime)=n
  comb.o=combn(k,2)
  to<-vector("list",2)
  ti<-lapply(to,function(x) c(1,2))
  gts<-expand.grid(ti)
  d.prime=sapply (1:dim(comb.o)[2], function(i){
    d=final.estimates[names(final.estimates)==paste(comb.o[,i],collapse="")]
    af=final.estimates[comb.o[,i]]
    to<-vector("list",2)
    ti<-lapply(to,function(x) c(1,2))
    gts<-expand.grid(ti)
    gts=gts[-c(1,dim(gts)[1]),]
    if(d<=0) {dmax=max(-prod(af),-prod(1-af))} else {dmax=min(apply(af*(gts==1)+(1-af)*(gts==2),2,prod) )}
    d.prime=sign(d)*(d/dmax)
    names(d.prime)=paste(comb.o[,i],collapse="")
    return(d.prime)
  })
  bp=sapply(1:length(final.d.prime), function(i) any(names(d.prime)==names(final.d.prime)[i]))
  final.d.prime[bp]=d.prime
  return(final.d.prime)
}



######################
#For all orders of LD
######################
theta32.prime=function(theta3,k) {
  #Defining initial values for parametes 
  theta2=dprime=rep(NA,(2^k)-1);theta2[1:k]=dprime[1:k]=theta3[1:k] 
  
  #Names of parameters in theta3 
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  names(theta3)=names(theta2)=names(dprime)=sapply(1:length(Pk),function(i) paste(unlist(Pk[i]),collapse=""))
  pos.hap=Hap(k) ; pos.hap=pos.hap[order(apply(pos.hap,1,sum),decreasing=T),]
  
  #theta3-->theta2  
  for(j in (k+1):(k+dim(combn(k,2))[2])) {
    ind=as.numeric(unlist(strsplit(names(theta2[j]),split="")))
    partition=t(setparts(length(ind)))
    d=sum(sapply (2:dim(partition)[1],function(i){
      a=partition[i,]; s=split(seq_along(a),a); l.s=length(s); coeff=Coef(l.s)
      eta=sapply(1:l.s ,function(m) {
        index=ind[unlist(s[m])]
        return(theta2[names(theta3)==paste(index,collapse="")])
      })
      return(coeff*prod(eta)) 
    })) 
    a=partition[1,]; s=split(seq_along(a),a)
    theta2[names(theta3)==paste(ind[unlist(s)],collapse="")]=theta3[names(theta3)==paste(ind[unlist(s)],collapse="")]-d
  }  
  
  dprime[(k+1):(k+dim(combn(k,2))[2])]=D.prime(theta3[1:(k+dim(combn(k,2))[2])],k)
  
  for(j in (k+dim(combn(k,2))[2]+1):(2^k-1)) {
    ind=as.numeric(unlist(strsplit(names(theta2[j]),split="")))
    partition=t(setparts(length(ind)))
    d=sapply(2:dim(partition)[1],function(i){
      a=partition[i,]; s=split(seq_along(a),a); l.s=length(s); coeff=Coef(l.s)
      eta=sapply(1:l.s ,function(m) {
        index=ind[unlist(s[m])]
        return(theta2[names(theta3)==paste(index,collapse="")])
      })
      return(c(coeff*prod(eta),min(eta)))}) 
    a=partition[1,]; s=split(seq_along(a),a)
    theta2[names(theta3)==paste(ind[unlist(s)],collapse="")]=theta3[names(theta3)==paste(ind[unlist(s)],collapse="")]-sum(d[1,])
    
    part2=t(combn(ind,length(ind)-1))
    part1=sapply(1:length(ind), function(i) setdiff(ind,part2[i,]))
    bpp=matrix(theta2[sapply(1:dim(part2)[1],function(i) c(which(names(theta2)==part1[i]),which(names(theta2)==apply(part2,1,paste,collapse="")[i])))],ncol=2,byrow=T)
    max.eta=max(0, bpp[,1]-(1-bpp[,2]))
    dmax=(sum(min(d[2,]),d[1,]))*(theta3[j]>0)+(sum(max.eta,d[1,]))*(theta3[j]<0)
    dprime[j]=ifelse((theta3[j]!=0),sign(theta3[j])*(theta3[j]/dmax),0)
  }
  return(list(theta2=theta2,dprime=dprime))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
######%%%%%%%%%%%%%%%%%%%%%%%%%              D = D' * Dmax                    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Function to go from standardized LD parameters to unstandardized ones.

######################
#For pairwise LD
######################
dprime.theta3.p=function(dprime,k) {
  n=names(dprime)[-(1:k)]
  final.d=rep(NA,length(n)) ; names(final.d)=n
  comb.o=combn(k,2)
  to<-vector("list",2)
  ti<-lapply(to,function(x) c(1,2))
  gts<-expand.grid(ti)
  final.d=sapply (1:dim(comb.o)[2], function(i){
    d=dprime[names(dprime)==paste(comb.o[,i],collapse="")]
    af=dprime[comb.o[,i]]
    to<-vector("list",2)
    ti<-lapply(to,function(x) c(1,2))
    gts<-expand.grid(ti)
    gts=gts[-c(1,dim(gts)[1]),]
    if(d<=0) {dmax=max(-prod(af),-prod(1-af))} else {dmax=min(apply(af*(gts==1)+(1-af)*(gts==2),2,prod) )}
    d.final=sign(d)*d*dmax
    names(d.final)=paste(comb.o[,i],collapse="")
    return(d.final)
  })
  return(final.d)
}



######################
# For all orders of LD
######################
dprime.theta3=function(dprime,k) {
  theta2=theta3=rep(NA,(2^k)-1);theta2[1:k]=theta3[1:k]=dprime[1:k] 
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  names(theta3)=names(theta2)=names(dprime)=sapply(1:length(Pk),function(i) paste(unlist(Pk[i]),collapse=""))
  pos.hap=Hap(k) ; pos.hap=pos.hap[order(apply(pos.hap,1,sum),decreasing=T),]
  
  theta3[(k+1):(k+dim(combn(k,2))[2])]=dprime.theta3.p(dprime[1:(k+dim(combn(k,2))[2])],k)
  for(j in (k+1):(k+dim(combn(k,2))[2])) {
    ind=as.numeric(unlist(strsplit(names(theta2[j]),split="")))
    partition=t(setparts(length(ind)))
    d=sum(sapply (2:dim(partition)[1],function(i){
      a=partition[i,]; s=split(seq_along(a),a); l.s=length(s); coeff=Coef(l.s)
      eta=sapply(1:l.s ,function(m) {
        index=ind[unlist(s[m])]
        return(theta2[names(theta3)==paste(index,collapse="")])
      })
      return(coeff*prod(eta)) 
    })) 
    a=partition[1,]; s=split(seq_along(a),a)
    theta2[names(theta3)==paste(ind[unlist(s)],collapse="")]=theta3[names(theta3)==paste(ind[unlist(s)],collapse="")]-d
  }  
  if(length(dprime)>3) {
    for(j in (k+dim(combn(k,2))[2]+1):(2^k-1)) {
      ind=as.numeric(unlist(strsplit(names(theta2[j]),split="")))
      partition=t(setparts(length(ind)))
      d=sapply(2:dim(partition)[1],function(i){
        a=partition[i,]; s=split(seq_along(a),a); l.s=length(s); coeff=Coef(l.s)
        eta=sapply(1:l.s ,function(m) {
          index=ind[unlist(s[m])]
          return(theta2[names(theta3)==paste(index,collapse="")])
        })
        return(c(coeff*prod(eta),min(eta)))}) 
      a=partition[1,]; s=split(seq_along(a),a) 
      part2=t(combn(ind,length(ind)-1))
      part1=sapply(1:length(ind), function(i) setdiff(ind,part2[i,]))
      bpp=matrix(theta2[sapply(1:dim(part2)[1],function(i) c(which(names(theta2)==part1[i]),which(names(theta2)==apply(part2,1,paste,collapse="")[i])))],ncol=2,byrow=T)
      max.eta=max(0, bpp[,1]-(1-bpp[,2]))
      dmax=(sum(min(d[2,]),d[1,]))*(dprime[j]>0)+(sum(max.eta,d[1,]))*(dprime[j]<0)
      theta3[j]=ifelse((dprime[j]!=0),sign(dprime[j])*dprime[j]*dmax,0)
      theta2[names(theta3)==paste(ind[unlist(s)],collapse="")]=theta3[names(theta3)==paste(ind[unlist(s)],collapse="")]-sum(d[1,])
    } }
  return(list(theta2=theta2,theta3=theta3))
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Function to compute the allele frequencies for each marker using observed counts.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
AlleleFreq=function(gts) { 
  sapply(1:ncol(gts),function(i) {
    tab=table(factor(gts[,i],levels=0:2))
    return((2*tab['2']+tab['1']) / (2*nrow(gts)))
  })
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#-logLik(theta1; SNPs) - Used Hierarchical Optimization     
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
MinusLogLikSNPs=function(theta3,D) {
  # theta3 = frequecy of '1' allele at a SINGLE locus\
  # D = genotypes at a SINGLE locus for ALL N individuals, \
  # D should be a vector of length N and for each i D[i] \in {0,1,2}
  # This function will return -log L(D;theta1)
  theta1=c(theta3,1-theta3)
  names(theta1)=c("1","0")
  D=cvtGts122(t(D))$genotypes
  return(-sum(apply(t(sapply(1:dim(D)[1],function(i) {
    c(log(theta1[names(theta1)==D[i,1]]),log(theta1[names(theta1)==D[i,2]]))
  }
    )),1,sum)))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# logLik(theta1; DPL) For Both Hierarchical and Simultaneous Optimization 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Loglikelihood.theta1=function(theta1,DPL,PHASED){
  #theta1 is the vector with the haplotype frequencies
  #For UNPHASED data DPL is a TRUE/FALSE matrix containing the possible diplotypes for ALL individuals. 
  #The dimension of DPL is N rows x 2*(2^k). 
  #Row corresponds to an individual. Columns correspond to 2*(2^k) possible diplotypes 
  if(PHASED==T) { 
    return(sum(apply(t(sapply(1:dim(DPL)[1],function(i) c(log(theta1[names(theta1)==DPL[i,1]]),log(theta1[names(theta1)==DPL[i,2]])))),1,sum)))
  } else {
    etaH1ANDetaH2=as.matrix(expand.grid(theta1,theta1))
    Names=as.matrix(expand.grid(names(theta1),names(theta1)))
    rownames(etaH1ANDetaH2)=paste(Names[,1], Names[,2], sep='/')
    etaH1xetaH2=apply(etaH1ANDetaH2,1,prod)   
    test=all(colnames(DPL)==rownames(etaH1ANDetaH2))
    EtaH1xEtaH2=matrix(rep(etaH1xetaH2,dim(DPL)[1]),byrow=T,nrow=dim(DPL)[1])
    ProbDPL=DPL*EtaH1xEtaH2
    Loglikelihood=sum(log(apply(ProbDPL,1,sum)))
    return(Loglikelihood)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# logLik(theta2; DPL)  For Hierarchical Optimization     
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Loglikelihood.Theta2.Unpooled=function(theta,prev.par,k,D,STAND,PHASED) {
  if(STAND==TRUE){
    dprime=c(prev.par,theta)
    theta3=dprime.theta3(dprime =dprime ,k =k )$theta3
    theta1=theta321(theta3 = theta3 ,k = k )$theta1
  } else{
    theta2=c(prev.par,theta)
    theta1=theta21(theta2,k)
  }
  if((any(theta1<0)|any(theta1>1)|((round(sum(theta1),5)-1)>0.00001))) {l=-log(1e-2400) } else {l=-Loglikelihood.theta1(theta1,D,PHASED)} 
  if(((theta<(-1))|(theta > 1))&(STAND==T)) print(sprintf('Warning: Tested value : %.4f outside (-1,1) interval', theta))
  if(((theta<(0))|(theta > 1))&(STAND==F)) print(sprintf('Warning: Tested value : %.4f outside (0,1) interval', theta))
  #print(sprintf('theta3/lik : %s', paste(theta,round(l,3),sep='/')))
  return(l)
}  

Loglikelihood.Theta2.Pooled=function(theta,prev.par.cases,prev.par.controls,k,D_controls,D_cases,STAND,PHASED) {
  l=Loglikelihood.Theta2.Unpooled(theta = theta,prev.par = prev.par.cases,k = k,D = D_cases,STAND = STAND,PHASED = PHASED) + 
    Loglikelihood.Theta2.Unpooled(theta = theta,prev.par = prev.par.controls,k = k,D = D_controls,STAND = STAND,PHASED = PHASED)
  print(l)
  return(l)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# logLik(theta3; DPL)  For Hierarchical Optimization  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Loglikelihood.theta3.Unpooled=function(theta,prev.par,k,D,interval.LD,STAND,PHASED) {
  if(STAND==T){
    dprime=c(prev.par,theta)   
    theta3=dprime.theta3(dprime,k)$theta3
  } else {theta3=c(prev.par,theta)}
  theta1=theta321(theta3,k)$theta1
  if((any(theta1<0)|any(theta1>1)|((round(sum(theta1),5)-1)>0.00001))) {l=-log(1e-2400) } else {l=-Loglikelihood.theta1(theta1,D,PHASED)} 
  if((theta<interval.LD[1])|(theta>interval.LD[2])) print(sprintf('Warning: Tested value : %.4f outside interval', theta))
  #print(sprintf('theta3/lik : %s', paste(theta,round(l,3),sep='/')))
  write.table(t(c(theta,l)),file='Values.txt',sep=';',append=T,col.names=F,row.names=F,quote=F)
  return(l)
}  

Loglikelihood.theta3.Pooled=function(theta,prev.par.cases,prev.par.controls,k,D_controls,D_cases,interval.LD,STAND,PHASED) {
  l=Loglikelihood.theta3.Unpooled(theta,prev.par.cases,k,D_cases,interval.LD,STAND,PHASED) + Loglikelihood.theta3.Unpooled(theta,prev.par.controls,k,D_controls,interval.LD,STAND,PHASED)
  return(l)}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# logLik(theta3; DPL)  For Simultaneous Optimization   
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For unpooled estimates     
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Loglikelihood.theta3.Unpooled.SimultaneousOptim=function(theta3,k,D,STAND,PHASED) {
    if(STAND==T) theta=dprime.theta3(theta3,k)$theta3
    theta1=theta321(theta3,k)$theta1
    #theta1=structure(par$cumu2multinomial(c(0,theta3)[c(1,4,3,7,2,6,5,8)]),.Names=eta1Names(k))
    if(any(theta1<0)&(min(theta1-0)<1e-8)){ theta1=(theta1+1e-8)/sum(theta1+1e-8)}
    if((any(theta1<0)|any(theta1>1)|((round(sum(theta1),5)-1)>0.00001))) {l=-log(1e-2400) } else {l=-Loglikelihood.theta1(theta1,D,PHASED)} 
    #print(l)
    return(l)
    }  
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For pooled estimates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AdjustInitialValuesForPooledOptim=function(k,ind.pooled,par.cases,par.controls) {
  bin.pooled=unlist(sapply (1:k, function(i) if(any(ind.pooled==i)) rep(TRUE,dim(combn(1:k,i))[2]) else rep(FALSE,dim(combn(1:k,i))[2])))
  NAMES=NULL
  Names=names(par.cases)
  a=sapply(1:length(bin.pooled), function(i) if(bin.pooled[i]) return(F) else return(T))
  InitValues=NULL
  for (i in 1:length(bin.pooled)) {
    if(a[i]) { 
      InitValues=c(InitValues,par.cases[i],par.controls[i]) 
      NAMES=c(NAMES,paste('ca',Names[i],sep='_'),paste('co',Names[i],sep='_'))
    } else { 
      InitValues=c(InitValues,par.cases[i])
      NAMES=c(NAMES,paste('p',Names[i],sep='_'))
    }
  }
  names(InitValues)=NAMES
  return(InitValues)
}

ProfileLoglikelihood.theta3.Unpooled.SimultaneousOptim=function(theta3, theta3.nuisance, k , ind.pooled, D_cases,D_controls,STAND,PHASED) {
  Pk=union(set_power(1:k),set_power(1:k))[-1]
  theta3.cases=theta3.controls=structure(rep(NA,2^k-1), .Names=deltaNames(k))
  ind.not.pooled=setdiff(1:k,ind.pooled)
  Boolean.Index.Not.Pooled=sapply(1:length(theta3.cases),function(j) length(unlist(strsplit(names(theta3.cases)[j],split='')))==ind.not.pooled)
  theta3.cases[Boolean.Index.Not.Pooled]=theta3[1:(length(theta3)/2)]
  theta3.controls[Boolean.Index.Not.Pooled]=theta3[((length(theta3)/2)+1):length(theta3)]
  theta3.cases[!Boolean.Index.Not.Pooled]=theta3.nuisance[1:(length(theta3.nuisance)/2)]
  
  ##Reconstruct parameter vector for cases and controls
  for(i in 1:length(theta)) {
    NAMES=unlist(strsplit(names(theta)[i],split='_'))
    if(NAMES[1]=='p') {theta.cases[names(theta.cases)==NAMES[2]]=theta.controls[names(theta.controls)==NAMES[2]]=theta[i]} else {
      if(NAMES[1]=='ca') {theta.cases[names(theta.cases)==NAMES[2]]=theta[i]} else { theta.controls[names(theta.controls)==NAMES[2]]=theta[i]} } }
  l=Loglikelihood.theta3.Unpooled.SimultaneousOptim(theta=theta.cases,k=k,D=D_cases,STAND=STAND,PHASED=PHASED) + Loglikelihood.theta3.Unpooled.SimultaneousOptim(theta=theta.controls,k,D=D_controls,STAND=STAND,PHASED=PHASED)
  return(l)
}

Loglikelihood.theta3.Pooled.SimultaneousOptim=function(theta,k,D_cases,D_controls,STAND,PHASED) {
  theta.cases=theta.controls=structure(rep(NA,2^k-1), .Names=deltaNames(k))
  ##Reconstruct parameter vector for cases and controls
  for(i in 1:length(theta)) {
    NAMES=unlist(strsplit(names(theta)[i],split='_'))
    if(NAMES[1]=='p') {theta.cases[names(theta.cases)==NAMES[2]]=theta.controls[names(theta.controls)==NAMES[2]]=theta[i]} else {
      if(NAMES[1]=='ca') {theta.cases[names(theta.cases)==NAMES[2]]=theta[i]} else { theta.controls[names(theta.controls)==NAMES[2]]=theta[i]} } }
  l=Loglikelihood.theta3.Unpooled.SimultaneousOptim(theta=theta.cases,k=k,D=D_cases,STAND=STAND,PHASED=PHASED) + Loglikelihood.theta3.Unpooled.SimultaneousOptim(theta=theta.controls,k,D=D_controls,STAND=STAND,PHASED=PHASED)
  return(l)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions Needed to Speed Optimization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This function takes as input a matrix of genotypes and returns a boolean matrix with all possible haplotypes. 
DPL_TrueFalse=function(gts){ 
  k=ncol(gts)
  pos.dipl=as.matrix(expand.grid(eta1Names(k),eta1Names(k)))
  pos.dipl=paste(pos.dipl[,1],pos.dipl[,2],sep='/')
  DPL=matrix(NA,ncol=length(pos.dipl),nrow=nrow(gts)); colnames(DPL)=pos.dipl  
  GTS=matrix(as.vector(t(cvtGts122(t(gts))$genotypes)),ncol = ncol(gts) * 2, byrow = T)
  DPL_ord=gts2dts(GTS);
  for(j in 1:length(DPL_ord)){
    DPL_bin=lapply(DPL_ord[j], function(hts) sapply(hts, function(h) paste(ord2bin(h, digits = ncol(gts)), collapse = ''))); 
    LengthDPL_bin=length(unlist(DPL_bin))
    if(LengthDPL_bin==2) { 
      DPL[j,which(pos.dipl==paste(unlist(DPL_bin),collapse='/'))]=TRUE
      DPL[j,-(which(pos.dipl==paste(unlist(DPL_bin),collapse='/')))]=FALSE
    } else{ 
      INDEX=sapply(seq(1,LengthDPL_bin,by=2), function(jk) which(pos.dipl==paste(unlist(DPL_bin)[jk:(jk+1)],collapse='/')))
      DPL[j,INDEX]=TRUE
      DPL[j,-INDEX]=FALSE
    }
  }
  return(DPL)
}

# Functions to compute intervals for LD parameters

correctFREQS=function(LD,k,prev.par,names.comb){
  theta3=structure(c(prev.par,LD), .Names = c(names(prev.par),names.comb))
  theta1=theta321(theta3,k)$theta1
  return(all(theta1>0)&all(theta1<1))
}

CI.optim=function(k,N=500,prev.par,names.comb, MaxLoops = 5){
  cond=TRUE
  for (iteration in 1:MaxLoops) {
    #print(paste('Attmpt nr. ', iteration, ' to determine the interval',sep=""))
    if(k>2){
      INDEX=strsplit(names(prev.par),split=NULL)
      LENGTH=sapply(1:length(INDEX), function(i) length(unlist(INDEX[i])))
      Limits=max(abs(prev.par[LENGTH==(k-1)]))
      LD=seq(-Limits-Limits*.1,Limits+Limits*.1,length.out=N)
    } else {LD=seq(-1,1,length.out=N)}
    
    CORRECTfreqs=unlist(lapply(LD,correctFREQS,k,prev.par,names.comb))
    N=N+500
    if (!((length(LD[CORRECTfreqs])==1)|(length(LD[CORRECTfreqs])==0)) ) break;
  }
  if (((length(LD[CORRECTfreqs])==1)|(length(LD[CORRECTfreqs])==0)) ) stop('no interval found');
  return(c(min(LD[CORRECTfreqs]),max(LD[CORRECTfreqs])))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
# Hierarchical Optimization Based on Theta 2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
OptimizerHierarchicalRaw.theta2=function(ind.pooled,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED,STAND,dprime=NULL,dprime.cases=NULL, dprime.controls=NULL) {
  gts.cases=first_allele_cases+second_allele_cases
  gts.controls=first_allele_controls+second_allele_controls
  k=dim(gts.cases)[2] ; 
  Ncases=dim(gts.cases)[1]; 
  Ncontrols=dim(gts.controls)[1] ;
  final.estimates.cases=final.variances.cases=final.likelihood.cases=final.estimates.controls=final.variances.controls=final.likelihood.controls=rep(NA,2^k-1)  
  names(final.estimates.cases)=names(final.estimates.controls)=names(final.variances.cases)=names(final.variances.controls)=names(final.likelihood.cases)=names(final.likelihood.controls)=deltaNames(k)
  
  # For allele frequencies 
  if(any(ind.pooled==1)) {
    final.estimates.cases[1:k]=final.estimates.controls[1:k]= AlleleFreq(rbind(gts.cases,gts.controls))  
    final.variances.cases[1:k]=final.variances.controls[1:k]=sqrt((final.estimates.cases[1:k]*(1-final.estimates.cases[1:k]))/(Ncases-1))
    final.likelihood.cases[1:k]=final.likelihood.controls[1:k]=sapply(1:k, function(i) MinusLogLikSNPs(theta3=final.estimates.cases[i],D=c(gts.cases[,i],gts.controls[,i]))) ;
  } else {
    final.estimates.cases[1:k]=AlleleFreq(rbind(gts.cases))  ; final.estimates.controls[1:k]=AlleleFreq(rbind(gts.controls))  
    final.variances.cases[1:k]=sqrt((final.estimates.cases[1:k]*(1-final.estimates.cases[1:k]))/(Ncases-1)) ; final.variances.controls[1:k]=sqrt((final.estimates.controls[1:k]*(1-final.estimates.controls[1:k]))/(Ncases-1))
    final.likelihood.cases[1:k]=sapply(1:k, function(i) MinusLogLikSNPs(theta3=final.estimates.cases[i],D=gts.cases[,i])) ; final.likelihood.controls[1:k]=sapply(1:k, function(i) MinusLogLikSNPs(theta3=final.estimates.controls[i],D=gts.controls[,i])) ;
  }
  
  # For LD parameters 
  for (o in 2:k){
    comb.o=combn(k,o)
    par.est=sapply (1:dim(comb.o)[2], function(i){ 
      print(c(o,i))
      if(PHASED==T) {   
        D_cases=cbind(sapply(1:dim(first_allele_cases)[1],function(j) paste(first_allele_cases[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_cases)[1],function(j) paste(second_allele_cases[j,comb.o[,i]],collapse="")) )
        D_controls=cbind(sapply(1:dim(first_allele_controls)[1],function(j) paste(first_allele_controls[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_controls)[1],function(j) paste(second_allele_controls[j,comb.o[,i]],collapse="")) )
        D=rbind(D_controls,D_cases)
      } else {
        D_cases=DPL_TrueFalse(gts.cases[,comb.o[,i]])
        D_controls=DPL_TrueFalse(gts.controls[,comb.o[,i]])
      }
      names.comb=paste(comb.o[,i],collapse="")
      Pk=union(set_power(comb.o[,i]),set_power(comb.o[,i]))[-1]
      prev.par.names=sapply(1:(length(Pk)-1),function(i) paste(unlist(Pk[i]),collapse=""))
      bp=sapply(1:length(final.estimates.cases), function(i) any(prev.par.names==names(final.estimates.cases)[i]))
      prev.par.cases= final.estimates.cases[bp]
      prev.par.controls= final.estimates.controls[bp]	    
      if(any(o==ind.pooled)) {   
        if(STAND){ #dprime[names.comb]-dprime[names.comb]*.001,dprime[names.comb]+dprime[names.comb]*.001
          op=optimize(f=Loglikelihood.Theta2.Pooled,interval=c(-1,1),prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,STAND=STAND,PHASED=PHASED)
        } else {
          op=optimize(f=Loglikelihood.Theta2.Pooled,interval=c(0,min(prev.par.cases,prev.par.controls)),prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,STAND=STAND,PHASED=PHASED)
        }
        var=hessian(Loglikelihood.Theta2.Pooled,op$minimum,prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,STAND=STAND,PHASED=PHASED)   
        return(c(op$minimum,var,op$objective,names.comb)) 
      } else {
        if(STAND){#dprime.controls[names.comb]-dprime.controls[names.comb]*.001,dprime.controls[names.comb]+dprime.controls[names.comb]*.001
          op.controls=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(-1,1),prev.par=prev.par.controls,k=o,D=D_controls,STAND=STAND,PHASED=PHASED)
          op.cases=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(-1,1),prev.par=prev.par.cases,k=o,D=D_cases,STAND=STAND,PHASED=PHASED,tol=.Machine$double.eps^2)         
        } else {
          op.controls=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(0,min(prev.par.controls)),prev.par=prev.par.controls,k=o,D=D_controls,STAND=STAND,PHASED=PHASED)
          op.cases=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(0,min(prev.par.cases)),prev.par=prev.par.cases,k=o,D=D_cases,STAND=STAND,PHASED=PHASED)
        }       
        var.con=hessian(Loglikelihood.Theta2.Unpooled,op.controls$minimum,prev.par=prev.par.controls,k=o,D=D_controls,STAND=STAND,PHASED=PHASED)
        var.cas=hessian(Loglikelihood.Theta2.Unpooled,op.cases$minimum,prev.par=prev.par.cases,k=o,D=D_cases,STAND=STAND,PHASED=PHASED)
        return(cbind(c(op.cases$minimum,var.cas,op.cases$objective,names.comb), c(op.controls$minimum,var.con,op.controls$objective,names.comb))) 
      }
    })
    par.names=par.est[4,]
    bp=sapply(1:length(final.estimates.cases), function(i) any(par.names==names(final.estimates.cases)[i]))
    if(any(o==ind.pooled)){
      par.est=matrix(as.numeric(par.est[-c(4),]),nrow=3)
      final.estimates.cases[bp]=final.estimates.controls[bp]=par.est[1,] 
      final.variances.cases[bp]=final.variances.controls[bp]=sqrt(1/par.est[2,])
      final.likelihood.cases[bp]=final.likelihood.controls[bp]=par.est[3,]
    } else {
      par.est=matrix(as.numeric(par.est[-c(4,8),]),nrow=6)
      final.estimates.cases[bp]=par.est[1,]     ; final.variances.cases[bp]=sqrt(1/par.est[2,])    ; final.likelihood.cases[bp]=par.est[3,]
      final.estimates.controls[bp]=par.est[4,]  ; final.variances.controls[bp]=sqrt(1/par.est[5,]) ; final.likelihood.controls[bp]=par.est[6,]
    }
  }
  return(cbind(final.estimates.cases,final.estimates.controls,final.variances.cases,final.variances.controls,final.likelihood.cases,final.likelihood.controls))
}

OptimizerHierarchical.theta2=function(ind.pooled,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED,STAND,dprime=NULL,dprime.cases=NULL, dprime.controls=NULL) {
  r = try(OptimizerHierarchicalRaw.theta2(ind.pooled,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED,STAND,dprime=dprime,dprime.cases=dprime.cases, dprime.controls=dprime.controls));
}

OptimizerHierarchicalRaw.theta2.RestHypothesis=function(ind.pooled,ParNull,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED,STAND) {
  gts.cases=first_allele_cases+second_allele_cases
  gts.controls=first_allele_controls+second_allele_controls
  k=dim(gts.cases)[2] ; 
  NCases=dim(gts.cases)[1]; 
  NControls=dim(gts.controls)[1] ;
  final.estimates.cases = ParNull[,'final.estimates.cases']
  final.variances.cases = ParNull[,'final.variances.cases'] 
  final.likelihood.cases = ParNull[,'final.likelihood.cases'] 
  final.estimates.controls = ParNull[,'final.estimates.controls'] 
  final.variances.controls = ParNull[,'final.variances.controls'] 
  final.likelihood.controls = ParNull[,'final.likelihood.controls'] 
  # For LD parameters 
  for (o in (max(ind.pooled)+1):k){
    comb.o=combn(k,o)
    par.est=sapply (1:dim(comb.o)[2], function(i){ 
      print(c(o,i))
      if(PHASED==T) {
        D_cases=cbind(sapply(1:dim(first_allele_cases)[1],function(j) paste(first_allele_cases[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_cases)[1],function(j) paste(second_allele_cases[j,comb.o[,i]],collapse="")) )
        D_controls=cbind(sapply(1:dim(first_allele_controls)[1],function(j) paste(first_allele_controls[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_controls)[1],function(j) paste(second_allele_controls[j,comb.o[,i]],collapse="")) )
        D=rbind(D_controls,D_cases)
      } else {
        D_cases=DPL_TrueFalse(gts.cases[,comb.o[,i]])
        D_controls=DPL_TrueFalse(gts.controls[,comb.o[,i]])
      }
      names.comb=paste(comb.o[,i],collapse="")
      Pk=union(set_power(comb.o[,i]),set_power(comb.o[,i]))[-1]
      prev.par.names=sapply(1:(length(Pk)-1),function(i) paste(unlist(Pk[i]),collapse=""))
      bp=sapply(1:length(final.estimates.cases), function(i) any(prev.par.names==names(final.estimates.cases)[i]))
      prev.par.cases= final.estimates.cases[bp]
      prev.par.controls= final.estimates.controls[bp]      
      
      if(any(o==ind.pooled)) {   
        if(STAND){ #dprime[names.comb]-dprime[names.comb]*.001,dprime[names.comb]+dprime[names.comb]*.001
          op=optimize(f=Loglikelihood.Theta2.Pooled,interval=c(-1,1),prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,STAND=STAND,PHASED=PHASED)
        } else {
          op=optimize(f=Loglikelihood.Theta2.Pooled,interval=c(0,min(prev.par.cases,prev.par.controls)),prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,STAND=STAND,PHASED=PHASED)
        }
        var=hessian(Loglikelihood.Theta2.Pooled,op$minimum,prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,STAND=STAND,PHASED=PHASED)   
        return(c(op$minimum,var,op$objective,names.comb)) 
      } else {
        if(STAND){#dprime.controls[names.comb]-dprime.controls[names.comb]*.001,dprime.controls[names.comb]+dprime.controls[names.comb]*.001
          op.controls=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(-1,1),prev.par=prev.par.controls,k=o,D=D_controls,STAND=STAND,PHASED=PHASED)
          op.cases=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(-1,1),prev.par=prev.par.cases,k=o,D=D_cases,STAND=STAND,PHASED=PHASED,tol=.Machine$double.eps^2)         
        } else {
          op.controls=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(0,min(prev.par.controls)),prev.par=prev.par.controls,k=o,D=D_controls,STAND=STAND,PHASED=PHASED)
          op.cases=optimize(f=Loglikelihood.Theta2.Unpooled,interval=c(0,min(prev.par.cases)),prev.par=prev.par.cases,k=o,D=D_cases,STAND=STAND,PHASED=PHASED)
        }       
        var.con=hessian(Loglikelihood.Theta2.Unpooled,op.controls$minimum,prev.par=prev.par.controls,k=o,D=D_controls,STAND=STAND,PHASED=PHASED)
        var.cas=hessian(Loglikelihood.Theta2.Unpooled,op.cases$minimum,prev.par=prev.par.cases,k=o,D=D_cases,STAND=STAND,PHASED=PHASED)
        return(cbind(c(op.cases$minimum,var.cas,op.cases$objective,names.comb), c(op.controls$minimum,var.con,op.controls$objective,names.comb))) 
      }
    })   
    par.names=par.est[4,]
    bp=sapply(1:length(final.estimates.cases), function(i) any(par.names==names(final.estimates.cases)[i]))
    
    if(any(o==ind.pooled)){
      par.est=matrix(as.numeric(par.est[-c(4),]),nrow=3)
      final.estimates.cases[bp]=final.estimates.controls[bp]=par.est[1,] 
      final.variances.cases[bp]=final.variances.controls[bp]=sqrt(1/par.est[2,])
      final.likelihood.cases[bp]=final.likelihood.controls[bp]=par.est[3,]
    } else {
      par.est=matrix(as.numeric(par.est[-c(4,8),]),nrow=6)
      final.estimates.cases[bp]=par.est[1,]     ; final.variances.cases[bp]=sqrt(1/par.est[2,])    ; final.likelihood.cases[bp]=par.est[3,]
      final.estimates.controls[bp]=par.est[4,]  ; final.variances.controls[bp]=sqrt(1/par.est[5,]) ; final.likelihood.controls[bp]=par.est[6,]
    }
  }
  return(cbind(final.estimates.cases,final.estimates.controls,final.variances.cases,final.variances.controls,final.likelihood.cases,final.likelihood.controls))
}

OptimizerHierarchical.theta2.RestHypothesis=function(ind.pooled,ParNull,first_allele_cases,second_allele_cases,first_allele_controls,second_allele_controls,PHASED,STAND) {
  r = try(OptimizerHierarchicalRaw.theta2.RestHypothesis(ind.pooled=ind.pooled,ParNull=ParNull,first_allele_cases=first_allele_cases,second_allele_cases=second_allele_cases,first_allele_controls=first_allele_controls,second_allele_controls=second_allele_controls,PHASED=PHASED,STAND=STAND));
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
# Hierarchical Optimization Based Theta 3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####

OptimizerHierarchicalRaw.theta3=function(ind.pooled,gts.cases,gts.controls,PHASED,STAND) {
  k=dim(gts.cases)[2] ; 
  NCases=dim(gts.cases)[1]; 
  NControls=dim(gts.controls)[1] ;
  final.estimates.cases=final.variances.cases=final.likelihood.cases=final.estimates.controls=final.variances.controls=final.likelihood.controls=rep(NA,2^k-1)  
  names(final.estimates.cases)=names(final.estimates.controls)=names(final.variances.cases)=names(final.variances.controls)=names(final.likelihood.cases)=names(final.likelihood.controls)=deltaNames(k)
  
  # For allele frequencies  
  if(any(ind.pooled==1)) {
    final.estimates.cases[1:k]=final.estimates.controls[1:k]= AlleleFreq(rbind(gts.cases,gts.controls))  
    final.variances.cases[1:k]=final.variances.controls[1:k]=sqrt(final.estimates.cases[1:k]*(1-final.estimates.cases[1:k]) ) 
    final.likelihood.cases[1:k]=final.likelihood.controls[1:k]=sapply(1:k, function(i) MinusLogLikSNPs(final.estimates.cases[i],D=c(gts.cases[,i],gts.controls[,i]))) ;
  } else {
    final.estimates.cases[1:k]=AlleleFreq(rbind(gts.cases))  ; final.estimates.controls[1:k]=AlleleFreq(rbind(gts.controls))  
    final.variances.cases[1:k]=sqrt(final.estimates.cases[1:k]*(1-final.estimates.cases[1:k])) ; final.variances.controls[1:k]=sqrt(final.estimates.controls[1:k]*(1-final.estimates.controls[1:k]))
    final.likelihood.cases[1:k]=sapply(1:k, function(i) MinusLogLikSNPs(final.estimates.cases[i],D=gts.cases[,i])) ; final.likelihood.controls[1:k]=sapply(1:k, function(i) MinusLogLikSNPs(final.estimates.controls[i],D=gts.controls[,i])) ;
  }
  
  # For LD parameters 
  for (o in 2:k){
    comb.o=combn(k,o)
    par.est=sapply (1:dim(comb.o)[2], function(i){ 
      print(c(o,i))
      if(PHASED==T) {
        D_cases=cbind(sapply(1:dim(first_allele_cases)[1],function(j) paste(first_allele_cases[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_cases)[1],function(j) paste(second_allele_cases[j,comb.o[,i]],collapse="")) )
        D_controls=cbind(sapply(1:dim(first_allele_controls)[1],function(j) paste(first_allele_controls[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_controls)[1],function(j) paste(second_allele_controls[j,comb.o[,i]],collapse="")) )
        D=rbind(D_controls,D_cases)
      } else {
        D_cases=DPL_TrueFalse(gts.cases[,comb.o[,i]])
        D_controls=DPL_TrueFalse(gts.controls[,comb.o[,i]])
      }
      names.comb=paste(comb.o[,i],collapse="")
      Pk=union(set_power(comb.o[,i]),set_power(comb.o[,i]))[-1]
      prev.par.names=sapply(1:(length(Pk)-1),function(i) paste(unlist(Pk[i]),collapse=""))
      bp=sapply(1:length(final.estimates.cases), function(i) any(prev.par.names==names(final.estimates.cases)[i]))
      prev.par.cases= final.estimates.cases[bp]
      prev.par.controls= final.estimates.controls[bp]	    
      
      if(any(o==ind.pooled)) {   
        if(STAND==T) { interval.LD = c(-1,1) } else {
          interval.LD.cases=CI.optim(o,N=500,prev.par.cases,names.comb) ; 
          interval.LD.controls=CI.optim(o,N=500,prev.par.controls,names.comb)       
          interval.LD=c(min(interval.LD.cases,interval.LD.controls),max(interval.LD.cases,interval.LD.controls))
        }
        op=optimize(f=Loglikelihood.theta3.Pooled,interval=interval.LD,prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,interval.LD=interval.LD,STAND=STAND,PHASED=PHASED)
        var=hessian(Loglikelihood.theta3.Pooled,op$minimum,prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,interval.LD=interval.LD,STAND=STAND,PHASED=PHASED)   
        if((op$minimum==interval.LD[1])|(op$minimum==interval.LD[2])) print('Warning! Optimal value is at the border of optimization interval!')
        return(c(op$minimum,var,op$objective,names.comb)) 
      } else {
        if(STAND==T) {interval.LD.cases=interval.LD.controls= c(-1,1)} else {
          interval.LD.cases=CI.optim(o,N=500,prev.par.cases,names.comb) ; interval.LD.controls=CI.optim(o,N=500,prev.par.controls,names.comb)}
        op.controls=optimize(f=Loglikelihood.theta3.Unpooled,interval=interval.LD.controls,prev.par=prev.par.controls,k=o,D=D_controls,interval.LD=interval.LD.controls,STAND=STAND,PHASED=PHASED)
        var.con=hessian(Loglikelihood.theta3.Unpooled,op.controls$minimum,prev.par=prev.par.controls,k=o,D=D_controls,interval.LD=interval.LD.controls,STAND=STAND,PHASED=PHASED)
        op.cases=optimize(f=Loglikelihood.theta3.Unpooled,interval=interval.LD.cases,prev.par=prev.par.cases,k=o,D=D_cases,interval.LD=interval.LD.cases,STAND=STAND,PHASED=PHASED)
        var.cas=hessian(Loglikelihood.theta3.Unpooled,op.cases$minimum,prev.par=prev.par.cases,k=o,D=D_cases,interval.LD=interval.LD.cases,STAND=STAND,PHASED=PHASED)
        return(cbind(c(op.cases$minimum,var.cas,op.cases$objective,names.comb), c(op.controls$minimum,var.con,op.controls$objective,names.comb))) 
      }
    })
    
    par.names=par.est[4,]
    bp=sapply(1:length(final.estimates.cases), function(i) any(par.names==names(final.estimates.cases)[i]))
    
    if(any(o==ind.pooled)){
      par.est=matrix(as.numeric(par.est[-c(4),]),nrow=3)
      final.estimates.cases[bp]=final.estimates.controls[bp]=par.est[1,] 
      final.variances.cases[bp]=final.variances.controls[bp]=sqrt(1/par.est[2,])
      final.likelihood.cases[bp]=final.likelihood.controls[bp]=par.est[3,]
    } else {
      par.est=matrix(as.numeric(par.est[-c(4,8),]),nrow=6)
      final.estimates.cases[bp]=par.est[1,]     ; final.variances.cases[bp]=sqrt(1/par.est[2,])    ; final.likelihood.cases[bp]=par.est[3,]
      final.estimates.controls[bp]=par.est[4,]  ; final.variances.controls[bp]=sqrt(1/par.est[5,]) ; final.likelihood.controls[bp]=par.est[6,]
    }
  }
  return(cbind(final.estimates.cases,final.estimates.controls,final.variances.cases,final.variances.controls,final.likelihood.cases,final.likelihood.controls))
}

OptimizerHierarchical.theta3=function(ind.pooled,gts.cases,gts.controls,PHASED,STAND) {
  r = try(OptimizerHierarchicalRaw.theta3(ind.pooled,gts.cases,gts.controls,PHASED,STAND));
}

# Optimizing Likelihood for Restricted Alternatives for pooled orders>2
OptimizerHierarchicalRaw.theta3.RestHypothesis=function(ind.pooled,ParAlt,gts.cases,gts.controls,PHASED,STAND) {
  k=dim(gts.cases)[2] ; 
  NCases=dim(gts.cases)[1]; 
  NControls=dim(gts.controls)[1] ;
  final.estimates.cases = ParAlt[,'final.estimates.cases']
  final.variances.cases = ParAlt[,'final.variances.cases'] 
  final.likelihood.cases = ParAlt[,'final.likelihood.cases'] 
  final.estimates.controls = ParAlt[,'final.estimates.controls'] 
  final.variances.controls = ParAlt[,'final.variances.controls'] 
  final.likelihood.controls = ParAlt[,'final.likelihood.controls'] 
  
  
  # For LD parameters  
  for (o in ind.pooled){
    comb.o=combn(k,o)
    par.est=sapply (1:dim(comb.o)[2], function(i){ 
      print(c(o,i))
      if(PHASED==T) {
        D_cases=cbind(sapply(1:dim(first_allele_cases)[1],function(j) paste(first_allele_cases[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_cases)[1],function(j) paste(second_allele_cases[j,comb.o[,i]],collapse="")) )
        D_controls=cbind(sapply(1:dim(first_allele_controls)[1],function(j) paste(first_allele_controls[j,comb.o[,i]],collapse="")), sapply(1:dim(second_allele_controls)[1],function(j) paste(second_allele_controls[j,comb.o[,i]],collapse="")) )
        D=rbind(D_controls,D_cases)
      } else {
        D_cases=DPL_TrueFalse(gts.cases[,comb.o[,i]])
        D_controls=DPL_TrueFalse(gts.controls[,comb.o[,i]])
      }
      names.comb=paste(comb.o[,i],collapse="")
      Pk=union(set_power(comb.o[,i]),set_power(comb.o[,i]))[-1]
      prev.par.names=sapply(1:(length(Pk)-1),function(i) paste(unlist(Pk[i]),collapse=""))
      bp=sapply(1:length(final.estimates.cases), function(i) any(prev.par.names==names(final.estimates.cases)[i]))
      prev.par.cases= final.estimates.cases[bp]
      prev.par.controls= final.estimates.controls[bp]	    
      
      if(any(o==ind.pooled)) {   
        if(STAND==T) { interval.LD = c(-1,1) } else {
          interval.LD.cases=CI.optim(o,N=500,prev.par.cases,names.comb) ; 
          interval.LD.controls=CI.optim(o,N=500,prev.par.controls,names.comb)       
          interval.LD=c(min(interval.LD.cases,interval.LD.controls),max(interval.LD.cases,interval.LD.controls))
        }
        op=optimize(f=Loglikelihood.theta3.Pooled,interval=interval.LD,prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,interval.LD=interval.LD,STAND=STAND,PHASED=PHASED)
        var=hessian(Loglikelihood.theta3.Pooled,op$minimum,prev.par.cases=prev.par.cases,prev.par.controls=prev.par.controls, k=o,D_controls=D_controls,D_cases=D_cases,interval.LD=interval.LD,STAND=STAND,PHASED=PHASED)   
        if((op$minimum==interval.LD[1])|(op$minimum==interval.LD[2])) print('Warning! Optimal value is at the border of optimization interval!')
        return(c(op$minimum,var,op$objective,names.comb)) 
      } else {
        if(STAND==T) {interval.LD.cases=interval.LD.controls= c(-1,1)} else {
          interval.LD.cases=CI.optim(o,N=500,prev.par.cases,names.comb) ; interval.LD.controls=CI.optim(o,N=500,prev.par.controls,names.comb)}
        op.controls=optimize(f=Loglikelihood.theta3.Unpooled,interval=interval.LD.controls,prev.par=prev.par.controls,k=o,D=D_controls,interval.LD=interval.LD.controls,STAND=STAND,PHASED=PHASED)
        var.con=hessian(Loglikelihood.theta3.Unpooled,op.controls$minimum,prev.par=prev.par.controls,k=o,D=D_controls,interval.LD=interval.LD.controls,STAND=STAND,PHASED=PHASED)
        op.cases=optimize(f=Loglikelihood.theta3.Unpooled,interval=interval.LD.cases,prev.par=prev.par.cases,k=o,D=D_cases,interval.LD=interval.LD.cases,STAND=STAND,PHASED=PHASED)
        var.cas=hessian(Loglikelihood.theta3.Unpooled,op.cases$minimum,prev.par=prev.par.cases,k=o,D=D_cases,interval.LD=interval.LD.cases,STAND=STAND,PHASED=PHASED)
        return(cbind(c(op.cases$minimum,var.cas,op.cases$objective,names.comb), c(op.controls$minimum,var.con,op.controls$objective,names.comb))) 
      }
    })
    par.names=par.est[4,]
    bp=sapply(1:length(final.estimates.cases), function(i) any(par.names==names(final.estimates.cases)[i]))
    if(any(o==ind.pooled)){
      par.est=matrix(as.numeric(par.est[-c(4),]),nrow=3)
      final.estimates.cases[bp]=final.estimates.controls[bp]=par.est[1,] 
      final.variances.cases[bp]=final.variances.controls[bp]=sqrt(1/par.est[2,])
      final.likelihood.cases[bp]=final.likelihood.controls[bp]=par.est[3,]
    } else {
      par.est=matrix(as.numeric(par.est[-c(4,8),]),nrow=6)
      final.estimates.cases[bp]=par.est[1,]     ; final.variances.cases[bp]=sqrt(1/par.est[2,])    ; final.likelihood.cases[bp]=par.est[3,]
      final.estimates.controls[bp]=par.est[4,]  ; final.variances.controls[bp]=sqrt(1/par.est[5,]) ; final.likelihood.controls[bp]=par.est[6,]
    }
  }
  return(cbind(final.estimates.cases,final.estimates.controls,final.variances.cases,final.variances.controls,final.likelihood.cases,final.likelihood.controls))
}

OptimizerHierarchical.theta3.RestHypothesis=function(ind.pooled,ParAlt,gts.cases,gts.controls,PHASED,STAND) {
  r = try(OptimizerHierarchicalRaw.theta3.RestHypothesis(ind.pooled=ind.pooled,ParAlt=ParAlt,gts.cases=gts.cases,gts.controls=gts.controls,PHASED=PHASED,STAND=STAND));
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Simultaneous Optimization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OptimizerSimultaneousRaw=function(ind.pooled, InitValues, gts.cases,gts.controls,PHASED=FALSE,STAND=FALSE, METHOD="Nelder-Mead",HESSIAN=FALSE) {
  k=dim(gts.cases)[2] ; 
  NCases=dim(gts.cases)[1] ; NControls=dim(gts.controls)[1]
  # Construct Diplotype Matrix Phased or Unphased Data
  if(PHASED==T) { ## Need to change this for phased data !!
    D_cases= NA
    D_controls= NA
  } else {
    D_cases=DPL_TrueFalse(gts.cases)
    D_controls=DPL_TrueFalse(gts.controls)
  }
  ### Optimizing under Null Hypothesis 
  if (length(ind.pooled)==k) { 
    print('Optimizing under the Null Hypothesis')
    op=optim(par=InitValues ,fn=Loglikelihood.theta3.Unpooled.SimultaneousOptim,gr = NULL,k=k,D=rbind(D_cases,D_controls),STAND=STAND,PHASED=PHASED, method = METHOD, hessian = HESSIAN,control = list(maxit = 30000))
    RETURN=list(Estimates=op$par, Covariance=op$hessian, Convergence=op$convergence, MinusLogLikelihood=op$value)
  } else {
    ### Optimizing under Un-restricted Alternative Hypothesis
    if(all(ind.pooled==0)) { 
      print('Optimizing under the Unrestricted Alternative Hypothesis')
      op.cases=optim(par=InitValues[seq(1,2*(2^k-1),by=2)] ,fn=Loglikelihood.theta3.Unpooled.SimultaneousOptim,gr = NULL,k=k,D=D_cases,STAND=STAND,PHASED=PHASED,method = METHOD, hessian = HESSIAN,control = list(maxit = 30000))
      op.controls=optim(par=InitValues[seq(2,2*(2^k-1),by=2)] ,fn=Loglikelihood.theta3.Unpooled.SimultaneousOptim,gr = NULL,k=k,D=D_controls,STAND=STAND,PHASED=PHASED,method = METHOD, hessian = HESSIAN,control = list(maxit = 30000))
      Estimates=cbind(op.cases$par,op.controls$par)
      Covariance=cbind(op.cases$hessian,op.controls$hessian)
      RETURN=list(Estimates=Estimates,Covariance=Covariance, Convergence=c(op.cases$convergence,op.controls$convergence),MinusLogLikelihood=op.cases$value+op.controls$value)
    } else {
      print('Optimizing under the Restricted Alternative Hypothesis')
      ### Optimizing under Restricted Alternatives Hypotheses
      op=optim(par=InitValues,fn=Loglikelihood.theta3.Pooled.SimultaneousOptim,gr = NULL,k=k,D_cases=D_cases,D_controls=D_controls,STAND=STAND,PHASED=PHASED, method = METHOD, hessian = HESSIAN,control = list(maxit = 30000))
      Estimates=matrix(rep(NA,2*2^k-2),ncol=2); rownames(Estimates)=deltaNames(k)
      for(i in 1:length(op$par)) {
        NAMES=unlist(strsplit(names(op$par)[i],split='_'))
        if(NAMES[1]=='p') {Estimates[rownames(Estimates)==NAMES[2],1]=Estimates[rownames(Estimates)==NAMES[2],2]=op$par[i]} else {
          if(NAMES[1]=='ca') {Estimates[rownames(Estimates)==NAMES[2],1]=op$par[i]} else {Estimates[rownames(Estimates)==NAMES[2],2]=op$par[i]} } }
      RETURN=list(Estimates=Estimates, Covariance=op$hessian, Convergence=op$convergence,MinusLogLikelihood=op$value)
      
    } 
  }
  return(RETURN)
}


OptimizerSimultaneous=function(ind.pooled, InitValues, gts.cases,gts.controls,PHASED,STAND) {
  r = try(OptimizerSimultaneousRaw(ind.pooled, InitValues, gts.cases,gts.controls,PHASED,STAND));
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# ADITIONAL FUNCTIONS - Not needed for Hierarchical LD, only for real data 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# From genotype counts to pairs of alleles 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

genexpand<-function(snpcounts, coding=NULL){
  p <- ncol(snpcounts)
  if (is.null(coding)) coding<-cbind(rep(1,p),rep(2,p))
  m<-matrix(ncol=2*p, nrow=nrow(snpcounts))
  for(i in 1:p){
    m[,2*i-1] <- coding[i,1+(snpcounts[,i]>0)]
    m[,2*i] <- coding[i,1+(snpcounts[,i]>1)]
  }
  nms <- colnames(snpcounts)
  rownames(m)<-rownames(snpcounts)
  if (!is.null(nms))
    colnames(m)<-as.vector(t(outer(nms,c(1,2),paste,sep="_")))
  m
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#       Functions to Read the Data from Binary Plink Format    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
readPlinkBinary = function(prefix) {
	      p = read.plink(sprintf("%s.bed", prefix), sprintf("%s.bim",prefix), sprintf("%s.fam", prefix));
	      p@.Data }

ExtractSNPS=function(data,t){
        index=unlist(sapply(1:length(t), function(i) which(colnames(data)==t[i])))
	gts=sapply(index,function(i) ifelse(as.integer(data[,i])== 0, NA, as.integer(data[,i])- 1));
	colnames(gts)=colnames(data)[index] ; rownames(gts)=rownames(data)
	return(gts)
	}
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Functions to Create and Optimize the Real Data      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
createDataWTCCC=function(NAMES_SNPs=NULL,NEIGHBOURS=1,DATA_CASES=NULL,DATA_CONTROLS=NULL,markEXCL=NULL,Fac=20) {
  ## NAMES_SNPs is the SNP we are interested in
  ## NEIGHBOURS is the number of neighbouring SNPs we want to include

  Fac=20
  
  for (NAMES_SNPs in NAMES_MODERATE){ 
    t=which(colnames(DATA_CASES)==NAMES_SNPs)
    t=colnames(DATA_CASES)[(t-(Fac*NEIGHBOURS)):(t+(Fac*NEIGHBOURS))]
    write.table(t, paste(NAMES_SNPs,'txt',sep='.'), quote=F, row.names=F, col.names=F)
    print(NAMES_SNPs)
  }
  
  
  t=which(colnames(DATA_CASES)==NAMES_SNPs)
  if(is.null(markEXCL)){
    t=colnames(DATA_CASES)[(t-NEIGHBOURS):(t+NEIGHBOURS)]
    gts.controls=ExtractSNPS(DATA_CONTROLS,t) ; 
    gts.cases=ExtractSNPS(DATA_CASES,t)
  } else { 
    t=colnames(DATA_CASES)[(t-(Fac*NEIGHBOURS)):(t+(Fac*NEIGHBOURS))]
    gts.controls=ExtractSNPS(DATA_CONTROLS,t) ; gts.cases=ExtractSNPS(DATA_CASES,t)
    
  DIMS=!sapply(1:dim(gts.cases)[2],function(i) (any(colnames(gts.cases)[i]==markEXCL)) |(dim(table(gts.cases[,i]))<3)|(dim(table(gts.controls[,i]))<3))
	LEFT=(1:(Fac))[DIMS[1:(Fac)]] ; LEFT=LEFT[(length(LEFT)-NEIGHBOURS+1):length(LEFT)]
	RIGHT=((Fac+2):(2*Fac+1))[DIMS[(Fac+2):(2*Fac+1)]] ; RIGHT=RIGHT[1:NEIGHBOURS]
	gts.cases=gts.cases[,c(LEFT,Fac+1,RIGHT)]
	gts.controls=gts.controls[,c(LEFT,Fac+1,RIGHT)]
	}
	k=dim(gts.cases)[2]
	print(paste('SNP ',NAMES_SNPs,' WITH ', NEIGHBOURS, ' NEIGHBOURS DONE'))
	r = list(gts.cases = gts.cases, gts.controls = gts.controls, k = k, NAMES_SNPs=NAMES_SNPs, NEIGHBOURS=NEIGHBOURS );
 	save(gts.cases, gts.controls , k , NAMES_SNPs, NEIGHBOURS,file=paste(NAMES_SNPs,'RData',sep='.'));
	#return(r)
	}

	

analayzeWTCCCregion = function(region) with(region, {
		    ptm <- proc.time()
		    ## Null Hypothesis   
		    ind.pooled=1:k
		    # Initial values for cases and controls
		    par.cases=c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    par.controls= c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    # Adjust Initial Values 
		    InitValues= AdjustInitialValuesForPooledOptim(ind.pooled,par.cases,par.controls)
		    op.NULL=OptimizerSimultaneous(ind.pooled, InitValues, gts.cases,gts.controls,PHASED=FALSE,STAND=TRUE) 
   
		    ## Full Alternative Hypothesis
		    ind.pooled=0
		    # Initial values for cases and controls
		    par.cases=c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    par.controls= c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    # Adjust Initial Values 
		    InitValues= AdjustInitialValuesForPooledOptim(ind.pooled,par.cases,par.controls)
		    op.ALTFULL=OptimizerSimultaneous(ind.pooled, InitValues,gts.cases,gts.controls,PHASED=FALSE,STAND=TRUE)    
		    dfALTFULL=2^k-1
		    chisqALTFULL=-2*(op.ALTFULL$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
		    pvalALTFULL=pchisq(chisqALTFULL,df=dfALTFULL,lower.tail = F)*(sum(op.NULL$Convergence,op.ALTFULL$Convergence)==0)+ 1*(sum(op.NULL$Convergence,op.ALTFULL$Convergence)!=0)
		    
		    ## Single SNP with Bonferroni Correction
		    ## STILL NEED TO CORRECT
		    dfSING=k
		    chisqSING=-2*(op.ALTFULL[1:k,5]+op.ALTFULL[1:k,6] -op.NULL[1:k,6])
		    pvalSING=pchisq(chisqSING,df=rep(dfSING,k),lower.tail = F)


		    ## Restricted Alternative Hypothesis 1 : Pooling allele frequencies
		    ind.pooled=1
		    # Initial values for cases and controls
		    par.cases=c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    par.controls= c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    # Adjust Initial Values 
		    InitValues= AdjustInitialValuesForPooledOptim(ind.pooled,par.cases,par.controls)
		    op.ALTRES1=OptimizerSimultaneous(ind.pooled, InitValues,gts.cases,gts.controls,PHASED=FALSE,STAND=TRUE)   
		    dfALTRES1=4
		    chisqALTRES1=-2*(op.ALTRES1$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
		    pvalALTRES1=pchisq(chisqALTRES1,df=dfALTRES1,lower.tail = F)*(sum(op.NULL$Convergence,op.ALTRES1$Convergence)==0)+ 1*(sum(op.NULL$Convergence,op.ALTRES1$Convergence)!=0)

		    ## Restricted Alternative Hypothesis 2 : Pooling allele frequencies and 2nd order LD parameters
		    ind.pooled=1:2
		    # Initial values for cases and controls
		    par.cases=c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    par.controls= c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    # Adjust Initial Values 
		    InitValues= AdjustInitialValuesForPooledOptim(ind.pooled,par.cases,par.controls)
		    op.ALTRES2=OptimizerSimultaneous(ind.pooled, InitValues,gts.cases,gts.controls,PHASED=FALSE,STAND=TRUE)   
		    dfALTRES2=1
		    chisqALTRES2=-2*(op.ALTRES2$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
		    pvalALTRES2=pchisq(chisqALTRES2,df=dfALTRES2,lower.tail = F)*(sum(op.NULL$Convergence,op.ALTRES2$Convergence)==0)+ 1*(sum(op.NULL$Convergence,op.ALTRES2$Convergence)!=0)


		    ## Restricted Alternative Hypothesis 3 : Pooling 2nd and 3rd order LD parameters
		    ind.pooled=2:3
		    # Initial values for cases and controls
		    par.cases=c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    par.controls= c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    # Adjust Initial Values 
		    InitValues= AdjustInitialValuesForPooledOptim(ind.pooled,par.cases,par.controls)
		    op.ALTRES3=OptimizerSimultaneous(ind.pooled, InitValues,gts.cases,gts.controls,PHASED=FALSE,STAND=TRUE)   
		    dfALTRES3=k
		    chisqALTRES3=-2*(op.ALTRES3$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
		    pvalALTRES3=pchisq(chisqALTRES3,df=dfALTRES3,lower.tail = F)*(sum(op.NULL$Convergence,op.ALTRES3$Convergence)==0)+ 1*(sum(op.NULL$Convergence,op.ALTRES3$Convergence)!=0)


		    ## Restricted Alternative Hypothesis 4 : Pooling 3rd order LD parameters
		    ind.pooled=3
		    # Initial values for cases and controls
		    par.cases=c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    par.controls= c(rep(0.5,k),unlist(sapply(2:k,function(i) rep(0,dim(combn(1:k,i))[2])))) 
		    # Adjust Initial Values 
		    InitValues= AdjustInitialValuesForPooledOptim(ind.pooled,par.cases,par.controls)
		    op.ALTRES4=OptimizerSimultaneous(ind.pooled, InitValues,gts.cases,gts.controls,PHASED=FALSE,STAND=TRUE)   
		    dfALTRES4=2*k
		    chisqALTRES4=-2*(op.ALTRES4$MinusLogLikelihood-op.NULL$MinusLogLikelihood)
		    pvalALTRES4=pchisq(chisqALTRES4,df=dfALTRES4,lower.tail = F)*(sum(op.NULL$Convergence,op.ALTRES4$Convergence)==0)+ 1*(sum(op.NULL$Convergence,op.ALTRES4$Convergence)!=0)


		    if((class(op.NULL) != 'try-error')&(class(op.ALTFULL) != 'try-error')&(class(op.ALTRES1) != 'try-error')&(class(op.ALTRES2) != 'try-error')&(class(op.ALTRES3) != 'try-error')&(class(op.ALTRES4) != 'try-error')){
		    write.table(rbind(op.NULL,op.ALTFULL,op.ALTRES1,op.ALTRES2,op.ALTRES3,op.ALTRES4),file=paste(NAMES_SNPs,'ResultsOptim',"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = TRUE)
		    write.table(t(c(NAMES_SNPs,chisqALTFULL, chisqALTRES1 , chisqALTRES2, chisqALTRES3, chisqALTRES4,chisqSING)),file=paste(NAMES_SNPs,'chisq',"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
		    write.table(t(c(NAMES_SNPs,pvalALTFULL,pvalALTRES1,pvalALTRES2,pvalALTRES3,pvalALTRES4,pvalSING)),file=paste(NAMES_SNPs,'pvals',"txt",sep="."),quote = FALSE, sep = ";", append = FALSE, dec = ".", row.names = FALSE, col.names = FALSE)
		    }
		    write.table(t(c(NAMES_SNPs,(proc.time()-ptm)[3])),file="time.txt",quote = FALSE, sep = ";", append = TRUE, dec = ".", row.names = FALSE, col.names = FALSE)
})


