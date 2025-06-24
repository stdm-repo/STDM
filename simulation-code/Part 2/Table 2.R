library(splines)
library(glmnet)
library(KRLS)
library(MASS)

################# Case 1 a #######################

################
## parameters ##
################
N=24 # number of locations
Time=220 # time

p=2
q=1
d=1
C=1
r=max(p,q) # we will discard first r time points
pqd=p+q+d

tseq=c(3,89,155,221) # true time change points


###############
## functions ##
###############
fjlt=function(j,l,t,Time){
  # this function is for fjl(t)
  if(l==1){
    return(sin(t/Time*pi)*2-2)
  }
  else if(l==2){
    return(sin(t/Time*pi)*2+2)
  }
  else{
    return(sin(t/Time*pi)*2-2)
  }
}

cjms=function(j,m,u,v){
  # this function is for cjm(s)
  return(0.05*u*v+sign(u*v)*0.1+0.02*j/4*(-1)^j)
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  error=0
  a_cand=sort(c(change_point,cand_item))
  LL=length(a_cand)
  N=length(XX)
  for(posN in 1:N){
    ZXi=XX[[posN]][,1:((1+pqd)*K)]
    for (i in 1:(LL+1)) {
      if(i==1){
        theta=ginv(t(ZXi[1:(a_cand[1]-r),])%*%ZXi[1:(a_cand[1]-r),])%*%
          t(ZXi[1:(a_cand[1]-r),])%*%Yt[posN,(r+1):a_cand[1]]
        error=error+sum((Yt[posN,(r+1):a_cand[1]]-t(theta)%*%t(ZXi[1:(a_cand[1]-r),]))^2)
      } else if(i==(LL+1)){
        theta=ginv(t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%
          t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%Yt[posN,(a_cand[LL]+1):ncol(Yt)]
        error=error+sum((Yt[posN,(a_cand[LL]+1):ncol(Yt)]-t(theta)%*%t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),]))^2)
      } else {
        theta=ginv(t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
                     ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
          t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%Yt[posN,(a_cand[i-1]+1):a_cand[i]]
        error=error+sum((Yt[posN,(a_cand[i-1]+1):a_cand[i]]-
                           t(theta)%*%t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),]))^2)
      }
      
    }
  }
  return(error)
}

#################
## simulations ##
#################
set.seed(12345)
change_list=list()
error_seq1_list=list()
error_seq2_list=list()
error_seq3_list=list()
for (rep in 1:200) {
  ### locations
  # each quadrant has 6 points
  Sj=matrix(runif(N*2,0,1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1]=-Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W=gausskernel(X = Sj,sigma=1)
  W=W-diag(N)
  W=t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X=matrix(rnorm(Time*N,0,1),nrow = N,ncol = Time)
  
  ### Initiate Y and Ysl
  Yt=matrix(0,nrow=N,ncol=2)
  Ytsl=W%*%Yt
  truebeta=list() # record the true beta for each time
  
  for(t in 3:Time){
    l=lfind(t,tseq)
    #print(l)
    newYt=matrix(0,nrow=N)
    betat=matrix(0,N,pqd+1)
    for(i in 1:N){
      betat[i,1]=0.1+sign(Sj[i,1]*Sj[i,2])*0.0075
      betat[i,2]=fjlt(1,l,t,Time)*cjms(1,1,Sj[i,1],Sj[i,2])
      betat[i,3]=fjlt(2,l,t,Time)*cjms(2,1,Sj[i,1],Sj[i,2])
      betat[i,4]=fjlt(3,l,t,Time)*cjms(3,1,Sj[i,1],Sj[i,2])
      betat[i,5]=fjlt(4,l,t,Time)*cjms(4,1,Sj[i,1],Sj[i,2])
      
      newYt[i]=sum(c(1,Ytsl[i,t-1],Ytsl[i,t-2],Yt[i,t-1],X[i,t])*betat[i,])+rnorm(1,0,1)*0.05
    }
    
    truebeta[[t]]=betat
    Yt=cbind(Yt,newYt)
    Ytsl=cbind(Ytsl,W%*%newYt)
  }
  
  ################
  ## Estimation ##
  ################
  ## cubic B-spline
  knots=c(74,148)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,221/Time))
  K=ncol(BS)
  
  result=list() # record the change points for each location
  XX=list() # record the X(si) for each location
  error_seq1=0
  
  for (posN in 1:N) {
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K*(Time-r)) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      for (j in 1:(t-r)) {
        Xsi[t-r,((j-1)*(1+pqd)*K+1):(j*(1+pqd)*K)]=Zjsi
      }
    }
    
    #best_lambda=cv.glmnet(Xsi/sqrt(Time-r), Yt[posN,(r+1):Time]/sqrt(Time-r), alpha = 1
    #,intercept = FALSE)$lambda.min
    best_lambda=cv.glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1
                          ,intercept = FALSE)$lambda.1se
    #best_model=glmnet(Xsi/sqrt(Time-r), Yt[posN,(r+1):Time]/sqrt(Time-r), alpha = 1,
    #                  lambda = 0.0042,intercept = FALSE)
    
    
    best_model=glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1,
                      lambda =  best_lambda, intercept = FALSE)
    hatgamma=coef(best_model)[-1]
    
    # reshape estimated gamma into matrix, each row represents a time point
    hatgamma2=matrix(hatgamma,nrow=Time-r,byrow = TRUE)
    aa=rowSums(hatgamma2!=0) # count time with non-zero parameters
    result[[posN]]=which(aa!=0)+r
    XX[[posN]]=Xsi
    
    for (t in (r+1):Time) {
      tmp <- hatgamma2[t-r, ]
      BS_tmp <- BS[t-r, ]
      tr <- truebeta[[t]][posN, ]
      
      # Reshape tmp and multiply
      m <- matrix(tmp, ncol = length(BS_tmp), byrow = TRUE)
      est_par <- m %*% BS_tmp
      
      # Update the error
      error_seq1 <- error_seq1 + sum((est_par - tr)^2)
    }
  }
  error_seq1_list[[rep]]=error_seq1/5/N/(Time-2)
  cand=sort(unique(unlist(result)))
  cand=cand[cand>(r+C*(1+pqd)*K)]
  cand=cand[cand<(Time-C*(1+pqd)*K)]
  
  ############
  ## Step 2 ##
  ############
  min_error=0
  change_point=c()
  
  
  error_seq=c()
  for (cand_item in cand) {
    error_seq=c(error_seq,detect(cand_item,c(),XX,Yt))
  }
  error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*(1+pqd)*K*N
  a1=cand[which.min(error_seq)]
  min_error=min(error_seq)
  change_point=c(change_point,a1)
  
  old_min_error=min_error
  LLL=2
  
  while(LLL){
    error_seq=c()
    cand=setdiff(cand,(a1-C*(1+pqd)*K):(a1+C*(1+pqd)*K))
    for (cand_item in cand) {
      error_seq=c(error_seq,detect(cand_item,change_point,XX,Yt))
    }
    error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*LLL*(1+pqd)*K*N
    if(min(error_seq)<old_min_error){
      old_min_error=min(error_seq)
      LLL=LLL+1
      a1=cand[which.min(error_seq)]
      change_point=sort(c(change_point,a1))
    } else{
      break
    }
  }
  
  change_list[[rep]]=change_point
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(74,148)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,221/Time))
  K=ncol(BS)
  
  
  change_point_new=c(r,change_point,Time)
  theta_re_est=list()
  error_seq2=0
  error_seq3=0
  
  for (posN in 1:N) {
    theta_re_est1=matrix(0,nrow=length(change_point_new)-1,ncol=(1+pqd)*K)
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      Xsi[t-r,]=Zjsi
    }
    
    for (t in 1:(length(change_point_new)-1)) {
      Xnew=Xsi[(change_point_new[t]+1-r):(change_point_new[t+1]-r),]
      Ytnew=Yt[posN,(change_point_new[t]+1):change_point_new[t+1]]
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 1, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 1,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
    
    # our method
    for (t in (r+1):Time) {
      l=lfind(t,change_point_new+1)
      tmp=theta_re_est1[l,]
      BS_tmp=BS[t-r,]
      tr=truebeta[[t]][posN,]
      est_par <- sapply(
        seq(1, length(tmp), by = length(BS_tmp)),
        function(i) sum(BS_tmp * tmp[i:(i + 2)])
      )
      error_seq2=error_seq2+sum((est_par-tr)^2)
    }
    
    # direct regression
    Xsi=matrix(0,nrow=Time-r,ncol=pqd) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=pqd)
      Zjsi=c(Ytsl[posN,t-1],Ytsl[posN,t-2],Yt[posN,t-1],X[posN,t])
      
      Xsi[t-r,]=Zjsi
    }
    est_par2=lm(matrix(t(Yt[posN,(r+1):Time]),nrow=Time-r,ncol=1)~Xsi)
    for (t in (r+1):Time) {
      error_seq3=error_seq3+sum((coef(est_par2)-truebeta[[t]][posN,])^2)
    }
  }
  
  error_seq2_list[[rep]]=error_seq2/5/N/(Time-2)
  error_seq3_list[[rep]]=error_seq3/5/N/(Time-2)
  
  print(rep)
}

saveRDS(error_seq1_list, file = "24_220_No_Structure.rds")
saveRDS(error_seq2_list, file = "24_220_OUR.rds")
saveRDS(error_seq3_list, file = "24_220_No_Dynamic.rds")

################# Case 1 b #######################

################
## parameters ##
################
N=48 # number of locations
Time=220 # time

p=2
q=1
d=1
C=1
r=max(p,q) # we will discard first r time points
pqd=p+q+d

tseq=c(3,89,155,221) # true time change points


###############
## functions ##
###############
fjlt=function(j,l,t,Time){
  # this function is for fjl(t)
  if(l==1){
    return(sin(t/Time*pi)*2-2)
  }
  else if(l==2){
    return(sin(t/Time*pi)*2+2)
  }
  else{
    return(sin(t/Time*pi)*2-2)
  }
}

cjms=function(j,m,u,v){
  # this function is for cjm(s)
  return(0.05*u*v+sign(u*v)*0.1+0.02*j/4*(-1)^j)
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  error=0
  a_cand=sort(c(change_point,cand_item))
  LL=length(a_cand)
  N=length(XX)
  for(posN in 1:N){
    ZXi=XX[[posN]][,1:((1+pqd)*K)]
    for (i in 1:(LL+1)) {
      if(i==1){
        theta=ginv(t(ZXi[1:(a_cand[1]-r),])%*%ZXi[1:(a_cand[1]-r),])%*%
          t(ZXi[1:(a_cand[1]-r),])%*%Yt[posN,(r+1):a_cand[1]]
        error=error+sum((Yt[posN,(r+1):a_cand[1]]-t(theta)%*%t(ZXi[1:(a_cand[1]-r),]))^2)
      } else if(i==(LL+1)){
        theta=ginv(t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%
          t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%Yt[posN,(a_cand[LL]+1):ncol(Yt)]
        error=error+sum((Yt[posN,(a_cand[LL]+1):ncol(Yt)]-t(theta)%*%t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),]))^2)
      } else {
        theta=ginv(t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
                     ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
          t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%Yt[posN,(a_cand[i-1]+1):a_cand[i]]
        error=error+sum((Yt[posN,(a_cand[i-1]+1):a_cand[i]]-
                           t(theta)%*%t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),]))^2)
      }
      
    }
  }
  return(error)
}


#################
## simulations ##
#################
set.seed(12345)
change_list=list()
error_seq1_list=list()
error_seq2_list=list()
error_seq3_list=list()
for (rep in 1:200) {
  ### locations
  # each quadrant has 6 points
  Sj=matrix(runif(N*2,0,1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1]=-Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W=gausskernel(X = Sj,sigma=1)
  W=W-diag(N)
  W=t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X=matrix(rnorm(Time*N,0,1),nrow = N,ncol = Time)
  
  ### Initiate Y and Ysl
  Yt=matrix(0,nrow=N,ncol=2)
  Ytsl=W%*%Yt
  truebeta=list() # record the true beta for each time
  
  for(t in 3:Time){
    l=lfind(t,tseq)
    #print(l)
    newYt=matrix(0,nrow=N)
    betat=matrix(0,N,pqd+1)
    for(i in 1:N){
      betat[i,1]=0.1+sign(Sj[i,1]*Sj[i,2])*0.0075
      betat[i,2]=fjlt(1,l,t,Time)*cjms(1,1,Sj[i,1],Sj[i,2])
      betat[i,3]=fjlt(2,l,t,Time)*cjms(2,1,Sj[i,1],Sj[i,2])
      betat[i,4]=fjlt(3,l,t,Time)*cjms(3,1,Sj[i,1],Sj[i,2])
      betat[i,5]=fjlt(4,l,t,Time)*cjms(4,1,Sj[i,1],Sj[i,2])
      
      newYt[i]=sum(c(1,Ytsl[i,t-1],Ytsl[i,t-2],Yt[i,t-1],X[i,t])*betat[i,])+rnorm(1,0,1)*0.05
    }
    
    truebeta[[t]]=betat
    Yt=cbind(Yt,newYt)
    Ytsl=cbind(Ytsl,W%*%newYt)
  }
  
  ################
  ## Estimation ##
  ################
  ## cubic B-spline
  knots=c(74,148)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,221/Time))
  K=ncol(BS)
  
  result=list() # record the change points for each location
  XX=list() # record the X(si) for each location
  error_seq1=0
  
  for (posN in 1:N) {
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K*(Time-r)) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      for (j in 1:(t-r)) {
        Xsi[t-r,((j-1)*(1+pqd)*K+1):(j*(1+pqd)*K)]=Zjsi
      }
    }
    
    #best_lambda=cv.glmnet(Xsi/sqrt(Time-r), Yt[posN,(r+1):Time]/sqrt(Time-r), alpha = 1
    #,intercept = FALSE)$lambda.min
    best_lambda=cv.glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1
                          ,intercept = FALSE)$lambda.1se
    #best_model=glmnet(Xsi/sqrt(Time-r), Yt[posN,(r+1):Time]/sqrt(Time-r), alpha = 1,
    #                  lambda = 0.0042,intercept = FALSE)
    
    
    best_model=glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1,
                      lambda =  best_lambda, intercept = FALSE)
    hatgamma=coef(best_model)[-1]
    
    # reshape estimated gamma into matrix, each row represents a time point
    hatgamma2=matrix(hatgamma,nrow=Time-r,byrow = TRUE)
    aa=rowSums(hatgamma2!=0) # count time with non-zero parameters
    result[[posN]]=which(aa!=0)+r
    XX[[posN]]=Xsi
    
    for (t in (r+1):Time) {
      tmp <- hatgamma2[t-r, ]
      BS_tmp <- BS[t-r, ]
      tr <- truebeta[[t]][posN, ]
      
      # Reshape tmp and multiply
      m <- matrix(tmp, ncol = length(BS_tmp), byrow = TRUE)
      est_par <- m %*% BS_tmp
      
      # Update the error
      error_seq1 <- error_seq1 + sum((est_par - tr)^2)
    }
  }
  error_seq1_list[[rep]]=error_seq1/5/N/(Time-2)
  cand=sort(unique(unlist(result)))
  cand=cand[cand>(r+C*(1+pqd)*K)]
  cand=cand[cand<(Time-C*(1+pqd)*K)]
  
  ############
  ## Step 2 ##
  ############
  min_error=0
  change_point=c()
  
  
  error_seq=c()
  for (cand_item in cand) {
    error_seq=c(error_seq,detect(cand_item,c(),XX,Yt))
  }
  error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*(1+pqd)*K*N
  a1=cand[which.min(error_seq)]
  min_error=min(error_seq)
  change_point=c(change_point,a1)
  
  old_min_error=min_error
  LLL=2
  
  while(LLL){
    error_seq=c()
    cand=setdiff(cand,(a1-C*(1+pqd)*K):(a1+C*(1+pqd)*K))
    for (cand_item in cand) {
      error_seq=c(error_seq,detect(cand_item,change_point,XX,Yt))
    }
    error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*LLL*(1+pqd)*K*N
    if(min(error_seq)<old_min_error){
      old_min_error=min(error_seq)
      LLL=LLL+1
      a1=cand[which.min(error_seq)]
      change_point=sort(c(change_point,a1))
    } else{
      break
    }
  }
  
  change_list[[rep]]=change_point
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(74,148)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,221/Time))
  K=ncol(BS)
  
  
  change_point_new=c(r,change_point,Time)
  theta_re_est=list()
  error_seq2=0
  error_seq3=0
  
  for (posN in 1:N) {
    theta_re_est1=matrix(0,nrow=length(change_point_new)-1,ncol=(1+pqd)*K)
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      Xsi[t-r,]=Zjsi
    }
    
    for (t in 1:(length(change_point_new)-1)) {
      Xnew=Xsi[(change_point_new[t]+1-r):(change_point_new[t+1]-r),]
      Ytnew=Yt[posN,(change_point_new[t]+1):change_point_new[t+1]]
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 1, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 1,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
    
    # our method
    for (t in (r+1):Time) {
      l=lfind(t,change_point_new+1)
      tmp=theta_re_est1[l,]
      BS_tmp=BS[t-r,]
      tr=truebeta[[t]][posN,]
      est_par <- sapply(
        seq(1, length(tmp), by = length(BS_tmp)),
        function(i) sum(BS_tmp * tmp[i:(i + 2)])
      )
      error_seq2=error_seq2+sum((est_par-tr)^2)
    }
    
    # direct regression
    Xsi=matrix(0,nrow=Time-r,ncol=pqd) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=pqd)
      Zjsi=c(Ytsl[posN,t-1],Ytsl[posN,t-2],Yt[posN,t-1],X[posN,t])
      
      Xsi[t-r,]=Zjsi
    }
    est_par2=lm(matrix(t(Yt[posN,(r+1):Time]),nrow=Time-r,ncol=1)~Xsi)
    for (t in (r+1):Time) {
      error_seq3=error_seq3+sum((coef(est_par2)-truebeta[[t]][posN,])^2)
    }
  }
  
  error_seq2_list[[rep]]=error_seq2/5/N/(Time-2)
  error_seq3_list[[rep]]=error_seq3/5/N/(Time-2)
  
  print(rep)
}

saveRDS(error_seq1_list, file = "48_220_No_Structure.rds")
saveRDS(error_seq2_list, file = "48_220_OUR.rds")
saveRDS(error_seq3_list, file = "48_220_No_Dynamic.rds")


################# Case 1 c #######################

################
## parameters ##
################
N=24 # number of locations
Time=440 # time

p=2
q=1
d=1
C=1
r=max(p,q) # we will discard first r time points
pqd=p+q+d

tseq=c(3,149,295,441) # true time change points


###############
## functions ##
###############
fjlt=function(j,l,t,Time){
  # this function is for fjl(t)
  if(l==1){
    return(sin(t/Time*pi)*2-2)
  }
  else if(l==2){
    return(sin(t/Time*pi)*2+2)
  }
  else{
    return(sin(t/Time*pi)*2-2)
  }
}

cjms=function(j,m,u,v){
  # this function is for cjm(s)
  return(0.05*u*v+sign(u*v)*0.1+0.02*j/4*(-1)^j)
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  error=0
  a_cand=sort(c(change_point,cand_item))
  LL=length(a_cand)
  N=length(XX)
  for(posN in 1:N){
    ZXi=XX[[posN]][,1:((1+pqd)*K)]
    for (i in 1:(LL+1)) {
      if(i==1){
        theta=ginv(t(ZXi[1:(a_cand[1]-r),])%*%ZXi[1:(a_cand[1]-r),])%*%
          t(ZXi[1:(a_cand[1]-r),])%*%Yt[posN,(r+1):a_cand[1]]
        error=error+sum((Yt[posN,(r+1):a_cand[1]]-t(theta)%*%t(ZXi[1:(a_cand[1]-r),]))^2)
      } else if(i==(LL+1)){
        theta=ginv(t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%
          t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%Yt[posN,(a_cand[LL]+1):ncol(Yt)]
        error=error+sum((Yt[posN,(a_cand[LL]+1):ncol(Yt)]-t(theta)%*%t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),]))^2)
      } else {
        theta=ginv(t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
                     ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
          t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%Yt[posN,(a_cand[i-1]+1):a_cand[i]]
        error=error+sum((Yt[posN,(a_cand[i-1]+1):a_cand[i]]-
                           t(theta)%*%t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),]))^2)
      }
      
    }
  }
  return(error)
}


#################
## simulations ##
#################
set.seed(12345)
change_list=list()
error_seq1_list=list()
error_seq2_list=list()
error_seq3_list=list()
for (rep in 1:200) {
  ### locations
  # each quadrant has 6 points
  Sj=matrix(runif(N*2,0,1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1]=-Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W=gausskernel(X = Sj,sigma=1)
  W=W-diag(N)
  W=t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X=matrix(rnorm(Time*N,0,1),nrow = N,ncol = Time)
  
  ### Initiate Y and Ysl
  Yt=matrix(0,nrow=N,ncol=2)
  Ytsl=W%*%Yt
  truebeta=list() # record the true beta for each time
  
  for(t in 3:Time){
    l=lfind(t,tseq)
    newYt=matrix(0,nrow=N)
    betat=matrix(0,N,pqd+1)
    for(i in 1:N){
      betat[i,1]=0.1+sign(Sj[i,1]*Sj[i,2])*0.0075
      betat[i,2]=fjlt(1,l,t,Time)*cjms(1,1,Sj[i,1],Sj[i,2])
      betat[i,3]=fjlt(2,l,t,Time)*cjms(2,1,Sj[i,1],Sj[i,2])
      betat[i,4]=fjlt(3,l,t,Time)*cjms(3,1,Sj[i,1],Sj[i,2])
      betat[i,5]=fjlt(4,l,t,Time)*cjms(4,1,Sj[i,1],Sj[i,2])
      
      newYt[i]=sum(c(1,Ytsl[i,t-1],Ytsl[i,t-2],Yt[i,t-1],X[i,t])*betat[i,])+rnorm(1,0,1)*0.05
    }
    
    truebeta[[t]]=betat
    Yt=cbind(Yt,newYt)
    Ytsl=cbind(Ytsl,W%*%newYt)
  }
  
  ################
  ## Estimation ##
  ################
  ## cubic B-spline
  knots=c(146,292)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,441/Time))
  K=ncol(BS)
  
  result=list() # record the change points for each location
  XX=list() # record the X(si) for each location
  error_seq1=0
  
  for (posN in 1:N) {
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K*(Time-r)) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      for (j in 1:(t-r)) {
        Xsi[t-r,((j-1)*(1+pqd)*K+1):(j*(1+pqd)*K)]=Zjsi
      }
    }
    
    best_lambda=cv.glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1
                          ,intercept = FALSE)$lambda.1se
    best_model=glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1,
                      lambda =  best_lambda, intercept = FALSE)
    hatgamma=coef(best_model)[-1]
    
    # reshape estimated gamma into matrix, each row represents a time point
    hatgamma2=matrix(hatgamma,nrow=Time-r,byrow = TRUE)
    aa=rowSums(hatgamma2!=0) # count time with non-zero parameters
    result[[posN]]=which(aa!=0)+r
    XX[[posN]]=Xsi
    
    for (t in (r+1):Time) {
      tmp <- hatgamma2[t-r, ]
      BS_tmp <- BS[t-r, ]
      tr <- truebeta[[t]][posN, ]
      
      # Reshape tmp and multiply
      m <- matrix(tmp, ncol = length(BS_tmp), byrow = TRUE)
      est_par <- m %*% BS_tmp
      
      # Update the error
      error_seq1 <- error_seq1 + sum((est_par - tr)^2)
    }
  }
  error_seq1_list[[rep]]=error_seq1/5/N/(Time-2)
  cand=sort(unique(unlist(result)))
  cand=cand[cand>(r+C*(1+pqd)*K)]
  cand=cand[cand<(Time-C*(1+pqd)*K)]
  
  ############
  ## Step 2 ##
  ############
  min_error=0
  change_point=c()
  
  
  error_seq=c()
  for (cand_item in cand) {
    error_seq=c(error_seq,detect(cand_item,c(),XX,Yt))
  }
  error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*(1+pqd)*K*N
  a1=cand[which.min(error_seq)]
  min_error=min(error_seq)
  change_point=c(change_point,a1)
  
  old_min_error=min_error
  LLL=2
  
  while(LLL){
    error_seq=c()
    cand=setdiff(cand,(a1-C*(1+pqd)*K):(a1+C*(1+pqd)*K))
    for (cand_item in cand) {
      error_seq=c(error_seq,detect(cand_item,change_point,XX,Yt))
    }
    error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*LLL*(1+pqd)*K*N
    if(min(error_seq)<old_min_error){
      old_min_error=min(error_seq)
      LLL=LLL+1
      a1=cand[which.min(error_seq)]
      change_point=sort(c(change_point,a1))
    } else{
      break
    }
  }
  
  change_list[[rep]]=change_point
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(146,292)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,441/Time))
  K=ncol(BS)
  
  
  change_point_new=c(r,change_point,Time)
  theta_re_est=list()
  error_seq2=0
  error_seq3=0
  
  for (posN in 1:N) {
    theta_re_est1=matrix(0,nrow=length(change_point_new)-1,ncol=(1+pqd)*K)
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      Xsi[t-r,]=Zjsi
    }
    
    for (t in 1:(length(change_point_new)-1)) {
      Xnew=Xsi[(change_point_new[t]+1-r):(change_point_new[t+1]-r),]
      Ytnew=Yt[posN,(change_point_new[t]+1):change_point_new[t+1]]
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 1, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 1,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
    
    # our method
    for (t in (r+1):Time) {
      l=lfind(t,change_point_new+1)
      tmp=theta_re_est1[l,]
      BS_tmp=BS[t-r,]
      tr=truebeta[[t]][posN,]
      est_par <- sapply(
        seq(1, length(tmp), by = length(BS_tmp)),
        function(i) sum(BS_tmp * tmp[i:(i + 2)])
      )
      error_seq2=error_seq2+sum((est_par-tr)^2)
    }
    
    # direct regression
    Xsi=matrix(0,nrow=Time-r,ncol=pqd) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=pqd)
      Zjsi=c(Ytsl[posN,t-1],Ytsl[posN,t-2],Yt[posN,t-1],X[posN,t])
      
      Xsi[t-r,]=Zjsi
    }
    est_par2=lm(matrix(t(Yt[posN,(r+1):Time]),nrow=Time-r,ncol=1)~Xsi)
    for (t in (r+1):Time) {
      error_seq3=error_seq3+sum((coef(est_par2)-truebeta[[t]][posN,])^2)
    }
  }
  
  error_seq2_list[[rep]]=error_seq2/5/N/(Time-2)
  error_seq3_list[[rep]]=error_seq3/5/N/(Time-2)
  print(rep)
}

saveRDS(error_seq1_list, file = "24_440_No_Structure.rds")
saveRDS(error_seq2_list, file = "24_440_OUR.rds")
saveRDS(error_seq3_list, file = "24_440_No_Dynamic.rds")



################# Case 1 d #######################

################
## parameters ##
################
N=48 # number of locations
Time=440 # time

p=2
q=1
d=1
C=1
r=max(p,q) # we will discard first r time points
pqd=p+q+d

tseq=c(3,149,295,441) # true time change points


###############
## functions ##
###############
fjlt=function(j,l,t,Time){
  # this function is for fjl(t)
  if(l==1){
    return(sin(t/Time*pi)*2-2)
  }
  else if(l==2){
    return(sin(t/Time*pi)*2+2)
  }
  else{
    return(sin(t/Time*pi)*2-2)
  }
}

cjms=function(j,m,u,v){
  # this function is for cjm(s)
  return(0.05*u*v+sign(u*v)*0.1+0.02*j/4*(-1)^j)
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  error=0
  a_cand=sort(c(change_point,cand_item))
  LL=length(a_cand)
  N=length(XX)
  for(posN in 1:N){
    ZXi=XX[[posN]][,1:((1+pqd)*K)]
    for (i in 1:(LL+1)) {
      if(i==1){
        theta=ginv(t(ZXi[1:(a_cand[1]-r),])%*%ZXi[1:(a_cand[1]-r),])%*%
          t(ZXi[1:(a_cand[1]-r),])%*%Yt[posN,(r+1):a_cand[1]]
        error=error+sum((Yt[posN,(r+1):a_cand[1]]-t(theta)%*%t(ZXi[1:(a_cand[1]-r),]))^2)
      } else if(i==(LL+1)){
        theta=ginv(t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%
          t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),])%*%Yt[posN,(a_cand[LL]+1):ncol(Yt)]
        error=error+sum((Yt[posN,(a_cand[LL]+1):ncol(Yt)]-t(theta)%*%t(ZXi[(a_cand[LL]-r+1):nrow(ZXi),]))^2)
      } else {
        theta=ginv(t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
                     ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%
          t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),])%*%Yt[posN,(a_cand[i-1]+1):a_cand[i]]
        error=error+sum((Yt[posN,(a_cand[i-1]+1):a_cand[i]]-
                           t(theta)%*%t(ZXi[(a_cand[i-1]-r+1):(a_cand[i]-r),]))^2)
      }
      
    }
  }
  return(error)
}


#################
## simulations ##
#################
set.seed(12345)
change_list=list()
error_seq1_list=list()
error_seq2_list=list()
error_seq3_list=list()
for (rep in 1:200) {
  ### locations
  # each quadrant has 6 points
  Sj=matrix(runif(N*2,0,1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1]=-Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W=gausskernel(X = Sj,sigma=1)
  W=W-diag(N)
  W=t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X=matrix(rnorm(Time*N,0,1),nrow = N,ncol = Time)
  
  ### Initiate Y and Ysl
  Yt=matrix(0,nrow=N,ncol=2)
  Ytsl=W%*%Yt
  truebeta=list() # record the true beta for each time
  
  for(t in 3:Time){
    l=lfind(t,tseq)
    newYt=matrix(0,nrow=N)
    betat=matrix(0,N,pqd+1)
    for(i in 1:N){
      betat[i,1]=0.1+sign(Sj[i,1]*Sj[i,2])*0.0075
      betat[i,2]=fjlt(1,l,t,Time)*cjms(1,1,Sj[i,1],Sj[i,2])
      betat[i,3]=fjlt(2,l,t,Time)*cjms(2,1,Sj[i,1],Sj[i,2])
      betat[i,4]=fjlt(3,l,t,Time)*cjms(3,1,Sj[i,1],Sj[i,2])
      betat[i,5]=fjlt(4,l,t,Time)*cjms(4,1,Sj[i,1],Sj[i,2])
      
      newYt[i]=sum(c(1,Ytsl[i,t-1],Ytsl[i,t-2],Yt[i,t-1],X[i,t])*betat[i,])+rnorm(1,0,1)*0.05
    }
    
    truebeta[[t]]=betat
    Yt=cbind(Yt,newYt)
    Ytsl=cbind(Ytsl,W%*%newYt)
  }
  
  ################
  ## Estimation ##
  ################
  ## cubic B-spline
  knots=c(146,292)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,441/Time))
  K=ncol(BS)
  
  result=list() # record the change points for each location
  XX=list() # record the X(si) for each location
  error_seq1=0
  
  for (posN in 1:N) {
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K*(Time-r)) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      for (j in 1:(t-r)) {
        Xsi[t-r,((j-1)*(1+pqd)*K+1):(j*(1+pqd)*K)]=Zjsi
      }
    }
    
    best_lambda=cv.glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1
                          ,intercept = FALSE)$lambda.1se
    best_model=glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1,
                      lambda =  best_lambda, intercept = FALSE)
    hatgamma=coef(best_model)[-1]
    
    # reshape estimated gamma into matrix, each row represents a time point
    hatgamma2=matrix(hatgamma,nrow=Time-r,byrow = TRUE)
    aa=rowSums(hatgamma2!=0) # count time with non-zero parameters
    result[[posN]]=which(aa!=0)+r
    XX[[posN]]=Xsi
    
    for (t in (r+1):Time) {
      tmp <- hatgamma2[t-r, ]
      BS_tmp <- BS[t-r, ]
      tr <- truebeta[[t]][posN, ]
      
      # Reshape tmp and multiply
      m <- matrix(tmp, ncol = length(BS_tmp), byrow = TRUE)
      est_par <- m %*% BS_tmp
      
      # Update the error
      error_seq1 <- error_seq1 + sum((est_par - tr)^2)
    }
  }
  error_seq1_list[[rep]]=error_seq1/5/N/(Time-2)
  cand=sort(unique(unlist(result)))
  cand=cand[cand>(r+C*(1+pqd)*K)]
  cand=cand[cand<(Time-C*(1+pqd)*K)]
  
  ############
  ## Step 2 ##
  ############
  min_error=0
  change_point=c()
  
  
  error_seq=c()
  for (cand_item in cand) {
    error_seq=c(error_seq,detect(cand_item,c(),XX,Yt))
  }
  error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*(1+pqd)*K*N
  a1=cand[which.min(error_seq)]
  min_error=min(error_seq)
  change_point=c(change_point,a1)
  
  old_min_error=min_error
  LLL=2
  
  while(LLL){
    error_seq=c()
    cand=setdiff(cand,(a1-C*(1+pqd)*K):(a1+C*(1+pqd)*K))
    for (cand_item in cand) {
      error_seq=c(error_seq,detect(cand_item,change_point,XX,Yt))
    }
    error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*LLL*(1+pqd)*K*N
    if(min(error_seq)<old_min_error){
      old_min_error=min(error_seq)
      LLL=LLL+1
      a1=cand[which.min(error_seq)]
      change_point=sort(c(change_point,a1))
    } else{
      break
    }
  }
  
  change_list[[rep]]=change_point
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(146,292)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,441/Time))
  K=ncol(BS)
  
  
  change_point_new=c(r,change_point,Time)
  theta_re_est=list()
  error_seq2=0
  error_seq3=0
  
  for (posN in 1:N) {
    theta_re_est1=matrix(0,nrow=length(change_point_new)-1,ncol=(1+pqd)*K)
    Xsi=matrix(0,nrow=Time-r,ncol=(1+pqd)*K) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=(1+pqd)*K)
      Zjsi[1:K]=BS[t-r,]
      Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
      Zjsi[(2*K+1):(3*K)]=Ytsl[posN,t-2]*BS[t-r,]
      Zjsi[(3*K+1):(4*K)]=Yt[posN,t-1]*BS[t-r,]
      Zjsi[(4*K+1):(5*K)]=X[posN,t]*BS[t-r,]
      
      Xsi[t-r,]=Zjsi
    }
    
    for (t in 1:(length(change_point_new)-1)) {
      Xnew=Xsi[(change_point_new[t]+1-r):(change_point_new[t+1]-r),]
      Ytnew=Yt[posN,(change_point_new[t]+1):change_point_new[t+1]]
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 1, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 1,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
    
    # our method
    for (t in (r+1):Time) {
      l=lfind(t,change_point_new+1)
      tmp=theta_re_est1[l,]
      BS_tmp=BS[t-r,]
      tr=truebeta[[t]][posN,]
      est_par <- sapply(
        seq(1, length(tmp), by = length(BS_tmp)),
        function(i) sum(BS_tmp * tmp[i:(i + 2)])
      )
      error_seq2=error_seq2+sum((est_par-tr)^2)
    }
    
    # direct regression
    Xsi=matrix(0,nrow=Time-r,ncol=pqd) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=pqd)
      Zjsi=c(Ytsl[posN,t-1],Ytsl[posN,t-2],Yt[posN,t-1],X[posN,t])
      
      Xsi[t-r,]=Zjsi
    }
    est_par2=lm(matrix(t(Yt[posN,(r+1):Time]),nrow=Time-r,ncol=1)~Xsi)
    for (t in (r+1):Time) {
      error_seq3=error_seq3+sum((coef(est_par2)-truebeta[[t]][posN,])^2)
    }
  }
  
  error_seq2_list[[rep]]=error_seq2/5/N/(Time-2)
  error_seq3_list[[rep]]=error_seq3/5/N/(Time-2)
  print(rep)
}

saveRDS(error_seq1_list, file = "48_440_No_Structure.rds")
saveRDS(error_seq2_list, file = "48_440_OUR.rds")
saveRDS(error_seq3_list, file = "48_440_No_Dynamic.rds")



################# Table 2 #######################
error_seq1_list1=readRDS("24_220_OUR.rds")
error_seq1_list2=readRDS("24_220_No_Dynamic.rds")
error_seq1_list3=readRDS("24_220_No_Structure.rds")
error_seq2_list1=readRDS("48_220_OUR.rds")
error_seq2_list2=readRDS("48_220_No_Dynamic.rds")
error_seq2_list3=readRDS("48_220_No_Structure.rds")
error_seq3_list1=readRDS("24_440_OUR.rds")
error_seq3_list2=readRDS("24_440_No_Dynamic.rds")
error_seq3_list3=readRDS("24_440_No_Structure.rds")
error_seq4_list1=readRDS("48_440_OUR.rds")
error_seq4_list2=readRDS("48_440_No_Dynamic.rds")
error_seq4_list3=readRDS("48_440_No_Structure.rds")

table_res=matrix(0,nrow=4,ncol=6)
colnames(table_res)=c("OUR","OUR_SD",
                      "No-Dynamic","No-Dynamic_SD",
                      "No-Structure","No-Structure_SD")

## Case 1(a)
table_res[1,c(1,2)]=c(mean(unlist(error_seq1_list1)),
                      sd(unlist(error_seq1_list1)))
table_res[1,c(3,4)]=c(mean(unlist(error_seq1_list2)),
                      sd(unlist(error_seq1_list2)))
table_res[1,c(5,6)]=c(mean(unlist(error_seq1_list3)),
                      sd(unlist(error_seq1_list3)))


## Case 1(b)
table_res[2,c(1,2)]=c(mean(unlist(error_seq2_list1)),
                      sd(unlist(error_seq2_list1)))
table_res[2,c(3,4)]=c(mean(unlist(error_seq2_list2)),
                      sd(unlist(error_seq2_list2)))
table_res[2,c(5,6)]=c(mean(unlist(error_seq2_list3)),
                      sd(unlist(error_seq2_list3)))


## Case 1(c)
table_res[3,c(1,2)]=c(mean(unlist(error_seq3_list1)),
                      sd(unlist(error_seq3_list1)))
table_res[3,c(3,4)]=c(mean(unlist(error_seq3_list2)),
                      sd(unlist(error_seq3_list2)))
table_res[3,c(5,6)]=c(mean(unlist(error_seq3_list3)),
                      sd(unlist(error_seq3_list3)))


## Case 1(d)
table_res[4,c(1,2)]=c(mean(unlist(error_seq4_list1)),
                      sd(unlist(error_seq4_list1)))
table_res[4,c(3,4)]=c(mean(unlist(error_seq4_list2)),
                      sd(unlist(error_seq4_list2)))
table_res[4,c(5,6)]=c(mean(unlist(error_seq4_list3)),
                      sd(unlist(error_seq4_list3)))

write.csv(table_res,file="table2.csv",quote=F,row.names = F)



