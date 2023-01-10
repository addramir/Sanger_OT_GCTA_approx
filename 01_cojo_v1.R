
#SNP				freq				b		se		p				n				freq_geno	bJ					bJ_se		pJ				LD_r
#1:109262261:C:G  0.0833874		-0.0887	0.0095	9.92501e-21		70467.1			0.09235		0.0551508			0.0103368   9.53374e-08	 -0.562439
#1:109275908:C:T. 0.778479			0.1602	0.0044	3.05679e-290		144450			0.7839		0.170274				0.00480653	6.89357e-275		0

setwd("~/Projects/GCTA/01_approx/")
source("01_utils.R")
library(data.table)

ref=fread("text_2snps.raw",data.table=F)

ref=ref[,7:8]
########
beta=c(-0.0887,-0.1602)
S=c(0.0095,0.0044)
n=c(94595,94595)
f=c(0.0833874,1-0.778479)
f_ref=apply(ref,MAR=2,FUN=function(x){mean(x)/2})
m_ref=10000
LD=cor(ref)

joint_version_0(beta=beta,S=S,n=n,f=f,ref=ref)
joint_version_1(beta=beta,S=S,n=n,f=f,ref=ref)
rm("ref")
joint_version_2(beta=beta,S=S,n=n,f=f,f_ref=f_ref,m_ref=10000,LD=LD)
joint_version_3(beta=beta,S=S,n=n,f=f,LD=LD)

###############

beta1 = -0.0887
beta2 = -0.1602
beta=c(beta1,beta2)
S=c(0.0095,0.0044)
N=length(beta)
m=nrow(ref)
n=c(94595,94595)

#f=apply(ref,2,FUN=function(x){mean(x)/2})
f=c(0.0833874,1-0.778479)

W=apply(ref,2,FUN=function(x){x-mean(x,na.rm=T)})

DW=diag(N)
diag(DW)=diag(t(W)%*%W)

D0=diag(N)
diag(D0)=2*f*(1-f)*n

#ys_y=mean(diag(D0)*S^2*(n-1)+diag(D0)*beta^2)
ys_y=median(2*f*(1-f)*n*S^2+2*f*(1-f)*beta^2)
Neff=ys_y/(S^2*(2*f*(1-f)))-(beta^2/S^2)+1

#Neff=c(70467.1,144450)




#Neff=(median(diag(D)*S^2*(n-1)+diag(D)*beta^2))/(diag(D)*S^2)-(beta^2/S^2)+1


B0=sqrt(D0)%*%solve(sqrt(DW))%*%t(W)%*%W%*%solve(sqrt(DW))%*%sqrt(D0)
beta_j0=solve(B0)%*%D0%*%beta

j=1
k=1
B=array(NA,c(N,N))
for (j in 1:N){
  k=1
  for (k in 1:N){
    B[j,k]=2*min(Neff[j],Neff[k])*(sqrt(f[k]*(1-f[k])*f[j]*(1-f[j]))/sqrt(sum(W[,j]^2)*sum(W[,k]^2)))*(t(W)%*%W)[j,k]
    #B[j,k]=min(Neff[j],Neff[k])*(sqrt(var(ref[,j])*var(ref[,k]))/sqrt(sum(W[,j]^2)*sum(W[,k]^2)))*(t(W)%*%W)[j,k]
  }
}


D=diag(N)
diag(D)=2*f*(1-f)*Neff
#diag(D)=diag(var(ref))*Neff

beta_j=solve(B)%*%D%*%beta

sigma_j=as.vector((ys_y-t(beta_j)%*%D%*%beta)/(min(Neff)-N))
varb=diag(sigma_j*solve(B))
se_b=sqrt(varb)

###############
.b=function(response,pred,S){
  out <- solve(S[pred,pred])%*%S[response,pred]
  return(out)
}

beta1 = -0.0887
beta2 = -0.1602

f1=mean(ref[,1])/2
f2=mean(ref[,2])/2
varg1= 2*f1*(1-f1)
varg2=2*f2*(1-f2)
r=cor(ref[,1],ref[,2])
covg=sqrt(varg1*varg2)*r
covy1=beta1*varg1
covy2=beta2*varg2

S=array(NA,c(3,3))
S[1,1]=1
S[2,2]=varg1
S[3,3]=varg2
S[1,2]=covy1
S[1,3]=covy2
S[2,3]=covg
S[lower.tri(S)]=S[upper.tri(S)]
.b(response=1,pred = 2:3,S=S)



