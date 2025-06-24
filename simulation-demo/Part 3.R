library(splines)
library(glmnet)
library(KRLS)
library(MASS)
library(mosum)
library(fossil)

################# Case 2 a #######################

################
## parameters ##
################
N=40 # number of locations
Time=220 # time
replication_time=1

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
  if(u>0 & v>0){
    return(0.04+0.02*j/4*(-1)^j)
  }
  if(u<0 & v>0){
    return(0.15+0.02*j/4*(-1)^j)
  }
  if(u<0 & v<0){
    return(-0.15+0.02*j/4*(-1)^j)
  }
  if(u>0 & v<0){
    return(-0.04+0.02*j/4*(-1)^j)
  }
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

#################
## simulations ##
#################
set.seed(12345)
location_list=list()
Dis_list=list()
Sj_list=list()
rand_index=c()
run_time=c()
for (rep in 1:replication_time) {
  ## Record start time
  start_time <- Sys.time()
  
  ### locations
  # each quadrant has 6 points
  Sj=matrix(rnorm(N*2,0.4,0.1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(3*N/M+1):N,1]=-Sj[(3*N/M+1):N,1]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj_list[[rep]]=Sj
  
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
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(89,155)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,221/Time))
  K=ncol(BS)
  
  change_point_new=c(r,89,155,Time)
  theta_re_est=list()
  
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
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 0, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 0,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
  }
  
  ## discrepancy measure
  Dis=matrix(0,nrow = N)
  thred=2/sqrt(N)
  for (posN in 1:N) {
    ds=Sj[posN,]
    cand1=c()
    cand2=c()
    cand3=c()
    cand4=c()
    A1=0
    A2=0
    A3=0
    A4=0
    for (posN1 in setdiff(1:N,posN)) {
      if(Sj[posN1,1]-ds[1]<=thred & Sj[posN1,1]-ds[1]>=0){
        if(Sj[posN1,2]-ds[2]<=thred & Sj[posN1,2]-ds[2]>=0){
          cand1=c(cand1,posN1)
          A1=A1+theta_re_est[[posN1]]
        }
        if(ds[2]-Sj[posN1,2]<=thred & ds[2]-Sj[posN1,2]>=0){
          cand3=c(cand3,posN1)
          A3=A3+theta_re_est[[posN1]]
        }
      }
      
      if(ds[1]-Sj[posN1,1]<=thred & ds[1]-Sj[posN1,1]>=0){
        if(Sj[posN1,2]-ds[2]<=thred & Sj[posN1,2]-ds[2]>=0){
          cand2=c(cand2,posN1)
          A2=A2+theta_re_est[[posN1]]
        }
        if(ds[2]-Sj[posN1,2]<=thred & ds[2]-Sj[posN1,2]>=0){
          cand4=c(cand4,posN1)
          A4=A4+theta_re_est[[posN1]]
        }
      }
    }
    
    A1=A1/ifelse(length(cand1)==0, 1, length(cand1))
    A2=A2/ifelse(length(cand2)==0, 1, length(cand2))
    A3=A3/ifelse(length(cand3)==0, 1, length(cand3))
    A4=A4/ifelse(length(cand4)==0, 1, length(cand4))
    
    Dis[posN]=sum((A1-A2)^2)+sum((A2-A3)^2)+sum((A3-A4)^2)+sum((A4-A1)^2)
  }
  
  location_list[[rep]]=mosum(as.vector(Dis),10)$cpts
  cpt.sort=location_list[[rep]]
  Dis_list[[rep]]=Dis
  
  true_label1=c(rep(1,N/M),rep(2,N/M),rep(3,N/M),rep(4,N/M))
  est_label=rep(1,N)
  if(length(cpt.sort)==0){
    rand_index=c(rand_index,rand.index(true_label1,est_label))
  }
  else{
    for (i in 1:(length(cpt.sort)+1)) {
      if(i==1){
        est_label[1:cpt.sort[1]]=i
      }
      else if(i==length(cpt.sort)+1){
        est_label[(cpt.sort[length(cpt.sort)]+1):N]=i
      }
      else{
        est_label[(cpt.sort[i-1]+1):cpt.sort[i]]=i
      }
    }
    rand_index=c(rand_index,rand.index(true_label1,est_label))
  }
  
  ## Record end time and compute duration
  end_time <- Sys.time()
  run_time <- c(run_time, as.numeric(difftime(end_time, start_time, units = "secs")))
  
  cat(sprintf("Dis Values: %.2f\n", Dis_list[[rep]]))
  cat(sprintf("Running time: %.2f seconds\n", run_time))
}


################# Case 2 b #######################

################
## parameters ##
################
N=80 # number of locations
Time=220 # time
replication_time=1

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
  if(u>0 & v>0){
    return(0.04+0.02*j/4*(-1)^j)
  }
  if(u<0 & v>0){
    return(0.15+0.02*j/4*(-1)^j)
  }
  if(u<0 & v<0){
    return(-0.15+0.02*j/4*(-1)^j)
  }
  if(u>0 & v<0){
    return(-0.04+0.02*j/4*(-1)^j)
  }
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

#################
## simulations ##
#################
set.seed(12345)
location_list=list()
Dis_list=list()
Sj_list=list()
rand_index=c()
run_time=c()
for (rep in 1:replication_time) {
  ## Record start time
  start_time <- Sys.time()
  
  ### locations
  # each quadrant has 6 points
  Sj=matrix(rnorm(N*2,0.4,0.1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(3*N/M+1):N,1]=-Sj[(3*N/M+1):N,1]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj_list[[rep]]=Sj
  
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
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(89,155)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,221/Time))
  K=ncol(BS)
  
  change_point_new=c(r,89,155,Time)
  theta_re_est=list()
  
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
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 0, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 0,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
  }
  
  ## discrepancy measure
  Dis=matrix(0,nrow = N)
  thred=2/sqrt(N)
  for (posN in 1:N) {
    ds=Sj[posN,]
    cand1=c()
    cand2=c()
    cand3=c()
    cand4=c()
    A1=0
    A2=0
    A3=0
    A4=0
    for (posN1 in setdiff(1:N,posN)) {
      #for (posN1 in 1:N) {
      if(Sj[posN1,1]-ds[1]<=thred & Sj[posN1,1]-ds[1]>=0){
        if(Sj[posN1,2]-ds[2]<=thred & Sj[posN1,2]-ds[2]>=0){
          cand1=c(cand1,posN1)
          A1=A1+theta_re_est[[posN1]]
        }
        if(ds[2]-Sj[posN1,2]<=thred & ds[2]-Sj[posN1,2]>=0){
          cand3=c(cand3,posN1)
          A3=A3+theta_re_est[[posN1]]
        }
      }
      
      if(ds[1]-Sj[posN1,1]<=thred & ds[1]-Sj[posN1,1]>=0){
        if(Sj[posN1,2]-ds[2]<=thred & Sj[posN1,2]-ds[2]>=0){
          cand2=c(cand2,posN1)
          A2=A2+theta_re_est[[posN1]]
        }
        if(ds[2]-Sj[posN1,2]<=thred & ds[2]-Sj[posN1,2]>=0){
          cand4=c(cand4,posN1)
          A4=A4+theta_re_est[[posN1]]
        }
      }
    }
    
    A1=A1/ifelse(length(cand1)==0, 1, length(cand1))
    A2=A2/ifelse(length(cand2)==0, 1, length(cand2))
    A3=A3/ifelse(length(cand3)==0, 1, length(cand3))
    A4=A4/ifelse(length(cand4)==0, 1, length(cand4))
    
    Dis[posN]=sum((A1-A2)^2)+sum((A2-A3)^2)+sum((A3-A4)^2)+sum((A4-A1)^2)
  }
  
  location_list[[rep]]=mosum(as.vector(Dis),14)$cpts
  cpt.sort=location_list[[rep]]
  Dis_list[[rep]]=Dis
  
  true_label1=c(rep(1,N/M),rep(2,N/M),rep(3,N/M),rep(4,N/M))
  est_label=rep(1,N)
  if(length(cpt.sort)==0){
    rand_index=c(rand_index,rand.index(true_label1,est_label))
  }
  else{
    for (i in 1:(length(cpt.sort)+1)) {
      if(i==1){
        est_label[1:cpt.sort[1]]=i
      }
      else if(i==length(cpt.sort)+1){
        est_label[(cpt.sort[length(cpt.sort)]+1):N]=i
      }
      else{
        est_label[(cpt.sort[i-1]+1):cpt.sort[i]]=i
      }
    }
    rand_index=c(rand_index,rand.index(true_label1,est_label))
  }
  
  ## Record end time and compute duration
  end_time <- Sys.time()
  run_time <- c(run_time, as.numeric(difftime(end_time, start_time, units = "secs")))
  
  cat(sprintf("Dis Values: %.2f\n", Dis_list[[rep]]))
  cat(sprintf("Running time: %.2f seconds\n", run_time))
}


################# Case 2 c #######################

################
## parameters ##
################
N=40 # number of locations
Time=440 # time
replication_time=1

p=2
q=1
d=1
C=1
r=max(p,q)+0 # we will discard first r time points
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
  if(u>0 & v>0){
    return(0.04+0.02*j/4*(-1)^j)
  }
  if(u<0 & v>0){
    return(0.15+0.02*j/4*(-1)^j)
  }
  if(u<0 & v<0){
    return(-0.15+0.02*j/4*(-1)^j)
  }
  if(u>0 & v<0){
    return(-0.04+0.02*j/4*(-1)^j)
  }
}

lfind=function(t,tseq){
  # this function is to find which breaks the time t belongs to
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

#################
## simulations ##
#################
set.seed(12345)
location_list=list()
Dis_list=list()
Sj_list=list()
rand_index=c()
run_time=c()
for (rep in 1:replication_time) {
  ## Record start time
  start_time <- Sys.time()
  
  ### locations
  # each quadrant has 6 points
  Sj=matrix(rnorm(N*2,0.4,0.1),nrow=N,ncol=2)
  M=4
  Sj[(N/M+1):(2*N/M),1]=-Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2]=-Sj[(3*N/M+1):N,2]
  Sj[(3*N/M+1):N,1]=-Sj[(3*N/M+1):N,1]
  Sj[(2*N/M+1):(3*N/M),2]=-Sj[(2*N/M+1):(3*N/M),2]
  Sj_list[[rep]]=Sj
  
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
  
  ############
  ## Step 3 ##
  ############
  ## re-estimate
  knots=c(146,292)/Time # knots for cubic spline
  BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,441/Time))
  K=ncol(BS)
  
  change_point_new=c(r,149,295,Time)
  theta_re_est=list()
  
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
      best_lambda <- cv.glmnet(Xnew, Ytnew, alpha = 0, 
                               intercept = FALSE)$lambda.1se
      best_model=glmnet(Xnew, Ytnew, alpha = 0,
                        lambda =  best_lambda, intercept = FALSE)
      theta_re_est1[t,]=coef(best_model)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)]=0
    theta_re_est[[posN]]=theta_re_est1
  }
  
  ## discrepancy measure
  Dis=matrix(0,nrow = N)
  thred=2/sqrt(N)
  for (posN in 1:N) {
    ds=Sj[posN,]
    cand1=c()
    cand2=c()
    cand3=c()
    cand4=c()
    A1=0
    A2=0
    A3=0
    A4=0
    for (posN1 in setdiff(1:N,posN)) {
      if(Sj[posN1,1]-ds[1]<=thred & Sj[posN1,1]-ds[1]>=0){
        if(Sj[posN1,2]-ds[2]<=thred & Sj[posN1,2]-ds[2]>=0){
          cand1=c(cand1,posN1)
          A1=A1+theta_re_est[[posN1]]
        }
        if(ds[2]-Sj[posN1,2]<=thred & ds[2]-Sj[posN1,2]>=0){
          cand3=c(cand3,posN1)
          A3=A3+theta_re_est[[posN1]]
        }
      }
      
      if(ds[1]-Sj[posN1,1]<=thred & ds[1]-Sj[posN1,1]>=0){
        if(Sj[posN1,2]-ds[2]<=thred & Sj[posN1,2]-ds[2]>=0){
          cand2=c(cand2,posN1)
          A2=A2+theta_re_est[[posN1]]
        }
        if(ds[2]-Sj[posN1,2]<=thred & ds[2]-Sj[posN1,2]>=0){
          cand4=c(cand4,posN1)
          A4=A4+theta_re_est[[posN1]]
        }
      }
    }
    
    A1=A1/ifelse(length(cand1)==0, 1, length(cand1))
    A2=A2/ifelse(length(cand2)==0, 1, length(cand2))
    A3=A3/ifelse(length(cand3)==0, 1, length(cand3))
    A4=A4/ifelse(length(cand4)==0, 1, length(cand4))
    
    Dis[posN]=sum((A1-A2)^2)+sum((A2-A3)^2)+sum((A3-A4)^2)+sum((A4-A1)^2)
  }
  
  location_list[[rep]]=mosum(as.vector(Dis),10)$cpts
  cpt.sort=location_list[[rep]]
  Dis_list[[rep]]=Dis
  
  true_label1=c(rep(1,N/M),rep(2,N/M),rep(3,N/M),rep(4,N/M))
  est_label=rep(1,N)
  if(length(cpt.sort)==0){
    rand_index=c(rand_index,rand.index(true_label1,est_label))
  }
  else{
    for (i in 1:(length(cpt.sort)+1)) {
      if(i==1){
        est_label[1:cpt.sort[1]]=i
      }
      else if(i==length(cpt.sort)+1){
        est_label[(cpt.sort[length(cpt.sort)]+1):N]=i
      }
      else{
        est_label[(cpt.sort[i-1]+1):cpt.sort[i]]=i
      }
    }
    rand_index=c(rand_index,rand.index(true_label1,est_label))
  }
  
  ## Record end time and compute duration
  end_time <- Sys.time()
  run_time <- c(run_time, as.numeric(difftime(end_time, start_time, units = "secs")))
  
  cat(sprintf("Dis Values: %.2f\n", Dis_list[[rep]]))
  cat(sprintf("Running time: %.2f seconds\n", run_time))
}

