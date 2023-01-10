######

joint_version_0=function(beta,S,n,f,ref){

  N=length(beta)
  m=nrow(ref)
  
  W=apply(ref,2,FUN=function(x){x-mean(x,na.rm=T)})
  
  DW=diag(N)
  diag(DW)=diag(t(W)%*%W)
  
  D0=diag(N)
  diag(D0)=2*f*(1-f)*n
  
  
  
  B0=sqrt(D0)%*%solve(sqrt(DW))%*%t(W)%*%W%*%solve(sqrt(DW))%*%sqrt(D0)
  beta_j0=solve(B0)%*%D0%*%beta
  
  ys_y=median(2*f*(1-f)*n*S^2+2*f*(1-f)*beta^2)
  Neff=ys_y/(S^2*(2*f*(1-f)))-(beta^2/S^2)+1
  
  sigma_j=ys_y
  varb=diag(sigma_j*solve(B0))
  se_b=sqrt(varb)
  
  out=list(b=as.vector(beta_j0),se=as.vector(se_b),Neff=Neff,ys_y=ys_y)
  
  out
}


######

joint_version_1=function(beta,S,n,f,ref){
  
  N=length(beta)
  m=nrow(ref)
  
  ys_y=median(2*f*(1-f)*n*S^2+2*f*(1-f)*beta^2)
  Neff=ys_y/(S^2*(2*f*(1-f)))-(beta^2/S^2)+1
  
  W=apply(ref,2,FUN=function(x){x-mean(x,na.rm=T)})
  
  DW=diag(N)
  diag(DW)=diag(t(W)%*%W)
  
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
  
  D_var=diag(N)
  diag(D_var)=2*f*(1-f)
  
  #sigma_j=as.vector((ys_y-t(beta_j)%*%D_var%*%beta)) #/(min(Neff)-N)
  sigma_j=ys_y
  varb=diag(sigma_j*solve(B))
  se_b=sqrt(varb)
  
  out=list(b=as.vector(beta_j),se=as.vector(se_b),Neff=Neff,ys_y=ys_y)
  
  out
}

######

joint_version_2=function(beta,S,n,f,f_ref,m_ref,LD){
  
  N=length(beta)
  
  ys_y=median(2*f*(1-f)*n*S^2+2*f*(1-f)*beta^2)
  Neff=ys_y/(S^2*(2*f*(1-f)))-(beta^2/S^2)+1
  
  varg_ref=2*f_ref*(1-f_ref)
  W_m=LD*(sqrt(varg_ref%o%varg_ref))*(m_ref-1)
  j=1
  k=1
  B=array(NA,c(N,N))
  for (j in 1:N){
    k=1
    for (k in 1:N){
      B[j,k]=2*min(Neff[j],Neff[k])*(sqrt(f[k]*(1-f[k])*f[j]*(1-f[j]))/sqrt(W_m[j,j]*W_m[k,k]))*W_m[j,k]
    }
  }
  
  
  D=diag(N)
  diag(D)=2*f*(1-f)*Neff
  #diag(D)=diag(var(ref))*Neff
  
  beta_j=solve(B)%*%D%*%beta
  
  D_var=diag(N)
  diag(D_var)=2*f*(1-f)
  
  #sigma_j=as.vector((ys_y-t(beta_j)%*%D_var%*%beta)) #/(min(Neff)-N)
  sigma_j=ys_y
  varb=diag(sigma_j*solve(B))
  se_b=sqrt(varb)
  
  out=list(b=as.vector(beta_j),se=as.vector(se_b),Neff=Neff,ys_y=ys_y)
  
  out
}

######
joint_version_3=function(beta,S,n,f,LD){
  
  N=length(beta)
  
  ys_y=median(2*f*(1-f)*n*S^2+2*f*(1-f)*beta^2)
  Neff=ys_y/(S^2*(2*f*(1-f)))-(beta^2/S^2)+1
  
  minNeff=array(NA,c(N,N))
  for (j in 1:N){
    k=1
    for (k in 1:N){
      minNeff[j,k]=min(Neff[j],Neff[k])
    }
  }
  
  varg=2*f*(1-f)
  VG=sqrt(varg%o%varg)
  
  B=minNeff*LD*VG
  
  D=diag(N)
  diag(D)=2*f*(1-f)*Neff
  #diag(D)=diag(var(ref))*Neff
  
  beta_j=solve(B)%*%D%*%beta
  
  D_var=diag(N)
  diag(D_var)=2*f*(1-f)
  
  #sigma_j=as.vector((ys_y-t(beta_j)%*%D_var%*%beta)) #/(min(Neff)-N)
  sigma_j=ys_y
  varb=diag(sigma_j*solve(B))
  se_b=sqrt(varb)
  
  out=list(b=as.vector(beta_j),se=as.vector(se_b),Neff=Neff,ys_y=ys_y)
  
  out
}
