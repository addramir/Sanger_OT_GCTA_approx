#Chr SNP bp  refA    freq    b   se  p   n   freq_geno   bJ  bJ_se   pJ  LD_r
#1   1:109113473:A:G 109113473   G   0.000672269 0.0017  0.0093  0.854958    2.80301e+06 0.0833  0.00396038  0.00930021  0.670226    -0.0520098
#1   1:109275908:C:T 109275908   T   0.778479    0.1602  0.0044  3.05679e-290    47458.9 0.7839  0.160213    0.00446108  1.8847e-282 0


setwd("~/Projects/GCTA/01_approx/")
source("01_utils.R")
library(data.table)

ref=fread("text_2snps_v2.raw",data.table=F)
ref=ref[,7:8]

beta=c(0.0017,-0.1602)
S=c(0.0093,0.0044)
n=c(94595,94595)
f=c(0.000672269,1-0.778479)
f_ref=apply(ref,MAR=2,FUN=function(x){mean(x)/2})
m_ref=10000
LD=cor(ref)


joint_version_0(beta=beta,S=S,n=n,f=f,ref=ref)
joint_version_1(beta=beta,S=S,n=n,f=f,ref=ref)
rm("ref")
joint_version_2(beta=beta,S=S,n=n,f=f,f_ref=f_ref,m_ref=10000,LD=LD)
joint_version_3(beta=beta,S=S,n=n,f=f,LD=LD)
