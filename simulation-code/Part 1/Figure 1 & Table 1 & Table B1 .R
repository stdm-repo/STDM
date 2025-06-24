library(splines)
library(glmnet)
library(KRLS)
library(MASS)
library(cluster)
library(fossil)
library(doParallel)
library(foreach)

num_cores <- parallel::detectCores() - 10  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

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
r=max(p,q) 
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
  # find which time segment t belongs to in tseq
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  # evaluate candidate change points by computing model error
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
change_list <- foreach(rep = 1:200, .packages = c("splines","glmnet","KRLS","MASS","fossil")) %dopar% {
  ## Set a reproducible seed per iteration
  set.seed(1234 + rep)
  
  ## Record start time
  start_time <- Sys.time()
  
  ### locations
  # each quadrant has 6 points
  Sj = matrix(runif(N*2,0,1), nrow = N, ncol = 2)
  M = 4
  Sj[(N/M+1):(2*N/M),1] = -Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2] = -Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2] = -Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1] = -Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W = gausskernel(X = Sj, sigma = 1)
  W = W - diag(N)
  W = t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X = matrix(rnorm(Time*N, 0, 1), nrow = N, ncol = Time)
  
  ### Initiate Y and Ysl
  Yt = matrix(0, nrow = N, ncol = 2)
  Ytsl = W %*% Yt
  truebeta = list() # record the true beta for each time
  
  for (t in 3:Time) {
    l = lfind(t, tseq)
    newYt = matrix(0, nrow = N)
    betat = matrix(0, N, pqd + 1)
    for (i in 1:N) {
      betat[i,1] = 0.1 + sign(Sj[i,1] * Sj[i,2]) * 0.0075
      betat[i,2] = fjlt(1, l, t, Time) * cjms(1, 1, Sj[i,1], Sj[i,2])
      betat[i,3] = fjlt(2, l, t, Time) * cjms(2, 1, Sj[i,1], Sj[i,2])
      betat[i,4] = fjlt(3, l, t, Time) * cjms(3, 1, Sj[i,1], Sj[i,2])
      betat[i,5] = fjlt(4, l, t, Time) * cjms(4, 1, Sj[i,1], Sj[i,2])
      
      newYt[i] = sum(c(1, Ytsl[i, t-1], Ytsl[i, t-2], Yt[i, t-1], X[i, t]) * betat[i,]) +
        rnorm(1, 0, 1) * 0.05
    }
    
    truebeta[[t]] = betat
    Yt = cbind(Yt, newYt)
    Ytsl = cbind(Ytsl, W %*% newYt)
  }
  
  ## Estimation using cubic B-splines
  knots = c(74, 148) / Time
  BS = ns(((r+1):Time)/Time, knots = knots, Boundary.knots = c(r/Time, 221/Time))
  K = ncol(BS)
  
  result = list()
  XX = list()
  
  for (posN in 1:N) {
    Xsi = matrix(0, nrow = Time-r, ncol = (1+pqd)*K*(Time-r))
    for (t in (r+1):Time) {
      Zjsi = matrix(0, nrow = 1, ncol = (1+pqd)*K)
      Zjsi[1:K] = BS[t-r,]
      Zjsi[(K+1):(2*K)] = Ytsl[posN,t-1] * BS[t-r,]
      Zjsi[(2*K+1):(3*K)] = Ytsl[posN,t-2] * BS[t-r,]
      Zjsi[(3*K+1):(4*K)] = Yt[posN,t-1] * BS[t-r,]
      Zjsi[(4*K+1):(5*K)] = X[posN,t] * BS[t-r,]
      
      for (j in 1:(t-r)) {
        Xsi[t-r, ((j-1)*(1+pqd)*K + 1):(j*(1+pqd)*K)] = Zjsi
      }
    }
    
    best_lambda = cv.glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1, intercept = FALSE)$lambda.1se
    best_model = glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1, lambda = best_lambda, intercept = FALSE)
    hatgamma = coef(best_model)[-1]
    hatgamma2 = matrix(hatgamma, nrow = Time-r, byrow = TRUE)
    aa = rowSums(hatgamma2 != 0)
    result[[posN]] = which(aa != 0) + r
    XX[[posN]] = Xsi
  }
  
  cand = sort(unique(unlist(result)))
  cand = cand[cand > (r + C*(1+pqd)*K)]
  cand = cand[cand < (Time - C*(1+pqd)*K)]
  
  ## Step 2: Select change points
  min_error = 0
  change_point = c()
  error_seq = c()
  for (cand_item in cand) {
    error_seq = c(error_seq, detect(cand_item, c(), XX, Yt))
  }
  error_seq = N*(Time-r)*log(error_seq/N/(Time-r)) + 2*(1+pqd)*K*N
  a1 = cand[which.min(error_seq)]
  min_error = min(error_seq)
  change_point = c(change_point, a1)
  old_min_error = min_error
  LLL = 2
  
  while (LLL) {
    error_seq = c()
    cand = setdiff(cand, (a1-C*(1+pqd)*K):(a1+C*(1+pqd)*K))
    for (cand_item in cand) {
      error_seq = c(error_seq, detect(cand_item, change_point, XX, Yt))
    }
    error_seq = N*(Time-r)*log(error_seq/N/(Time-r)) + 2*LLL*(1+pqd)*K*N
    if (min(error_seq) < old_min_error) {
      old_min_error = min(error_seq)
      LLL = LLL + 1
      a1 = cand[which.min(error_seq)]
      change_point = sort(c(change_point, a1))
    } else {
      break
    }
  }
  
  ## Compute Rand Index
  true_label1 = c(rep(1,89), rep(2,66), rep(3,65))
  est_label = rep(1,220)
  for (i in 1:(length(change_point)+1)) {
    if (i == 1) {
      est_label[1:change_point[1]] = i
    } else if (i == length(change_point)+1) {
      est_label[(change_point[length(change_point)]+1):220] = i
    } else {
      est_label[(change_point[i-1]+1):change_point[i]] = i
    }
  }
  rand_index = rand.index(true_label1, est_label)
  
  ## Record end time and compute duration
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ## Return structured result
  list(
    change_point = change_point,
    rand_index = rand_index,
    run_time_seconds = run_time
  )
}                                                    

saveRDS(change_list, file = "24_220.rds")


################# Case 1 c #######################

################
## parameters ##
################
N=24 # number of locations
Time=440 # time

p=2
q=1
d=1
C=2
r=max(p,q) 
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
  # find which time segment t belongs to in tseq
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  # evaluate candidate change points by computing model error
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
change_list <- foreach(rep = 1:200, .packages = c("splines", "glmnet", "KRLS", "MASS", "fossil")) %dopar% {
  ## Set a reproducible seed
  set.seed(1234 + rep)
  
  ## Record start time
  start_time <- Sys.time()
  
  ### locations
  Sj = matrix(runif(N*2, 0, 1), nrow = N, ncol = 2)
  M = 4
  Sj[(N/M+1):(2*N/M),1] = -Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2] = -Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2] = -Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1] = -Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W = gausskernel(X = Sj, sigma = 1)
  W = W - diag(N)
  W = t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X = matrix(rnorm(Time * N, 0, 1), nrow = N, ncol = Time)
  
  ### Initialize Y and Ysl
  Yt = matrix(0, nrow = N, ncol = 2)
  Ytsl = W %*% Yt
  truebeta = list()
  
  for (t in 3:Time) {
    l = lfind(t, tseq)
    newYt = matrix(0, nrow = N)
    betat = matrix(0, N, pqd + 1)
    for (i in 1:N) {
      betat[i,1] = 0.1 + sign(Sj[i,1] * Sj[i,2]) * 0.0075
      betat[i,2] = fjlt(1, l, t, Time) * cjms(1, 1, Sj[i,1], Sj[i,2])
      betat[i,3] = fjlt(2, l, t, Time) * cjms(2, 1, Sj[i,1], Sj[i,2])
      betat[i,4] = fjlt(3, l, t, Time) * cjms(3, 1, Sj[i,1], Sj[i,2])
      betat[i,5] = fjlt(4, l, t, Time) * cjms(4, 1, Sj[i,1], Sj[i,2])
      
      newYt[i] = sum(c(1, Ytsl[i,t-1], Ytsl[i,t-2], Yt[i,t-1], X[i,t]) * betat[i,]) +
        rnorm(1, 0, 1) * 0.05
    }
    truebeta[[t]] = betat
    Yt = cbind(Yt, newYt)
    Ytsl = cbind(Ytsl, W %*% newYt)
  }
  
  ## Estimation
  knots = c(146, 292) / Time
  BS = ns(((r+1):Time)/Time, knots = knots, Boundary.knots = c(r/Time, 441/Time))
  K = ncol(BS)
  
  result = list()
  XX = list()
  
  for (posN in 1:N) {
    Xsi = matrix(0, nrow = Time - r, ncol = (1 + pqd) * K * (Time - r))
    for (t in (r+1):Time) {
      Zjsi = matrix(0, nrow = 1, ncol = (1 + pqd) * K)
      Zjsi[1:K] = BS[t-r,]
      Zjsi[(K+1):(2*K)] = Ytsl[posN, t-1] * BS[t-r,]
      Zjsi[(2*K+1):(3*K)] = Ytsl[posN, t-2] * BS[t-r,]
      Zjsi[(3*K+1):(4*K)] = Yt[posN, t-1] * BS[t-r,]
      Zjsi[(4*K+1):(5*K)] = X[posN, t] * BS[t-r,]
      
      for (j in 1:(t - r)) {
        Xsi[t - r, ((j - 1) * (1 + pqd) * K + 1):(j * (1 + pqd) * K)] = Zjsi
      }
    }
    
    best_lambda = cv.glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1, intercept = FALSE)$lambda.1se
    best_model = glmnet(Xsi, Yt[posN,(r+1):Time], alpha = 1, lambda = best_lambda, intercept = FALSE)
    hatgamma = coef(best_model)[-1]
    hatgamma2 = matrix(hatgamma, nrow = Time - r, byrow = TRUE)
    aa = rowSums(hatgamma2 != 0)
    result[[posN]] = which(aa != 0) + r
    XX[[posN]] = Xsi
  }
  
  cand = sort(unique(unlist(result)))
  cand = cand[cand > (r + C * (1 + pqd) * K)]
  cand = cand[cand < (Time - C * (1 + pqd) * K)]
  
  ## Step 2
  min_error = 0
  change_point = c()
  error_seq = sapply(cand, function(cand_item) detect(cand_item, c(), XX, Yt))
  error_seq = N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * (1 + pqd) * K * N
  a1 = cand[which.min(error_seq)]
  min_error = min(error_seq)
  change_point = c(change_point, a1)
  old_min_error = min_error
  LLL = 2
  
  while (LLL) {
    cand = setdiff(cand, (a1 - C * (1 + pqd) * K):(a1 + C * (1 + pqd) * K))
    error_seq = sapply(cand, function(cand_item) detect(cand_item, change_point, XX, Yt))
    error_seq = N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * LLL * (1 + pqd) * K * N
    if (min(error_seq) < old_min_error) {
      old_min_error = min(error_seq)
      a1 = cand[which.min(error_seq)]
      change_point = sort(c(change_point, a1))
      LLL = LLL + 1
    } else {
      break
    }
  }
  
  ## Rand Index
  true_label1 = c(rep(1, 149), rep(2, 146), rep(3, 145))
  est_label = rep(1, 440)
  for (i in 1:(length(change_point) + 1)) {
    if (i == 1) {
      est_label[1:change_point[1]] = i
    } else if (i == length(change_point) + 1) {
      est_label[(change_point[length(change_point)] + 1):440] = i
    } else {
      est_label[(change_point[i - 1] + 1):change_point[i]] = i
    }
  }
  rand_index = rand.index(true_label1, est_label)
  
  ## Record end time and compute runtime
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ## Return result
  list(
    change_point = change_point,
    rand_index = rand_index,
    run_time_seconds = run_time
  )
}


saveRDS(change_list, file = "24_440.rds")


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
  # find which time segment t belongs to in tseq
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  # evaluate candidate change points by computing model error
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
change_list <- foreach(rep = 1:200, .packages = c("splines", "glmnet", "KRLS", "MASS", "fossil")) %dopar% {
  set.seed(1234 + rep)  # unique seed for each repetition
  start_time <- Sys.time()
  
  ### locations
  Sj = matrix(runif(N * 2, 0, 1), nrow = N, ncol = 2)
  M = 4
  Sj[(N/M+1):(2*N/M),1] = -Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2] = -Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2] = -Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1] = -Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W = gausskernel(X = Sj, sigma = 1)
  W = W - diag(N)
  W = t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X = matrix(rnorm(Time * N, 0, 1), nrow = N, ncol = Time)
  
  ### Initialize Y and Ysl
  Yt = matrix(0, nrow = N, ncol = 2)
  Ytsl = W %*% Yt
  truebeta = list()
  
  for (t in 3:Time) {
    l = lfind(t, tseq)
    newYt = matrix(0, nrow = N)
    betat = matrix(0, N, pqd + 1)
    for (i in 1:N) {
      betat[i,1] = 0.1 + sign(Sj[i,1] * Sj[i,2]) * 0.0075
      betat[i,2] = fjlt(1, l, t, Time) * cjms(1, 1, Sj[i,1], Sj[i,2])
      betat[i,3] = fjlt(2, l, t, Time) * cjms(2, 1, Sj[i,1], Sj[i,2])
      betat[i,4] = fjlt(3, l, t, Time) * cjms(3, 1, Sj[i,1], Sj[i,2])
      betat[i,5] = fjlt(4, l, t, Time) * cjms(4, 1, Sj[i,1], Sj[i,2])
      
      newYt[i] = sum(c(1, Ytsl[i,t-1], Ytsl[i,t-2], Yt[i,t-1], X[i,t]) * betat[i,]) +
        rnorm(1, 0, 1) * 0.05
    }
    truebeta[[t]] = betat
    Yt = cbind(Yt, newYt)
    Ytsl = cbind(Ytsl, W %*% newYt)
  }
  
  ### Estimation
  knots = c(74, 148) / Time
  BS = ns(((r+1):Time) / Time, knots = knots, Boundary.knots = c(r / Time, 221 / Time))
  K = ncol(BS)
  
  result = list()
  XX = list()
  
  for (posN in 1:N) {
    Xsi = matrix(0, nrow = Time - r, ncol = (1 + pqd) * K * (Time - r))
    for (t in (r+1):Time) {
      Zjsi = matrix(0, nrow = 1, ncol = (1 + pqd) * K)
      Zjsi[1:K] = BS[t - r,]
      Zjsi[(K+1):(2*K)] = Ytsl[posN, t-1] * BS[t - r,]
      Zjsi[(2*K+1):(3*K)] = Ytsl[posN, t-2] * BS[t - r,]
      Zjsi[(3*K+1):(4*K)] = Yt[posN, t-1] * BS[t - r,]
      Zjsi[(4*K+1):(5*K)] = X[posN, t] * BS[t - r,]
      
      for (j in 1:(t - r)) {
        Xsi[t - r, ((j - 1) * (1 + pqd) * K + 1):(j * (1 + pqd) * K)] = Zjsi
      }
    }
    
    best_lambda = cv.glmnet(Xsi, Yt[posN, (r+1):Time], alpha = 1, intercept = FALSE)$lambda.1se
    best_model = glmnet(Xsi, Yt[posN, (r+1):Time], alpha = 1, lambda = best_lambda, intercept = FALSE)
    hatgamma = coef(best_model)[-1]
    
    hatgamma2 = matrix(hatgamma, nrow = Time - r, byrow = TRUE)
    aa = rowSums(hatgamma2 != 0)
    result[[posN]] = which(aa != 0) + r
    XX[[posN]] = Xsi
  }
  
  cand = sort(unique(unlist(result)))
  cand = cand[cand > (r + C * (1 + pqd) * K)]
  cand = cand[cand < (Time - C * (1 + pqd) * K)]
  
  ## Step 2
  min_error = 0
  change_point = c()
  
  error_seq = sapply(cand, function(cand_item) detect(cand_item, c(), XX, Yt))
  error_seq = N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * (1 + pqd) * K * N
  a1 = cand[which.min(error_seq)]
  min_error = min(error_seq)
  change_point = c(change_point, a1)
  
  old_min_error = min_error
  LLL = 2
  while (LLL) {
    cand = setdiff(cand, (a1 - C * (1 + pqd) * K):(a1 + C * (1 + pqd) * K))
    error_seq = sapply(cand, function(cand_item) detect(cand_item, change_point, XX, Yt))
    error_seq = N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * LLL * (1 + pqd) * K * N
    if (min(error_seq) < old_min_error) {
      old_min_error = min(error_seq)
      a1 = cand[which.min(error_seq)]
      change_point = sort(c(change_point, a1))
      LLL = LLL + 1
    } else {
      break
    }
  }
  
  true_label1 = c(rep(1, 89), rep(2, 66), rep(3, 65))
  est_label = rep(1, 220)
  for (i in 1:(length(change_point) + 1)) {
    if (i == 1) {
      est_label[1:change_point[1]] = i
    } else if (i == length(change_point) + 1) {
      est_label[(change_point[length(change_point)] + 1):220] = i
    } else {
      est_label[(change_point[i - 1] + 1):change_point[i]] = i
    }
  }
  rand_index = rand.index(true_label1, est_label)
  
  ## End time and result
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(
    change_point = change_point,
    rand_index = rand_index,
    run_time_seconds = run_time
  )
}


saveRDS(change_list, file = "48_220.rds")


################# Case 1 d #######################

################
## parameters ##
################
N=48 # number of locations
Time=440 # time

p=2
q=1
d=1
C=2
r=max(p,q)
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
  # find which time segment t belongs to in tseq
  i=1
  while(t>=tseq[i]){
    i=i+1
  }
  return(i-1)
}

detect=function(cand_item,change_point,XX,Yt){
  # evaluate candidate change points by computing model error
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
change_list <- foreach(rep = 1:200, .packages = c("splines", "glmnet", "KRLS", "MASS", "fossil")) %dopar% {
  
  # Set unique seed for reproducibility
  set.seed(1234 + rep)
  start_time <- Sys.time()
  
  ### locations
  Sj = matrix(runif(N * 2, 0, 1), nrow = N, ncol = 2)
  M = 4
  Sj[(N/M+1):(2*N/M),1] = -Sj[(N/M+1):(2*N/M),1]
  Sj[(3*N/M+1):N,2] = -Sj[(3*N/M+1):N,2]
  Sj[(2*N/M+1):(3*N/M),2] = -Sj[(2*N/M+1):(3*N/M),2]
  Sj[(2*N/M+1):(3*N/M),1] = -Sj[(2*N/M+1):(3*N/M),1]
  
  ### pre-defined W
  W = gausskernel(X = Sj, sigma = 1)
  W = W - diag(N)
  W = t(apply(W, 1, function(X) X / sum(X)))
  
  ### X
  X = matrix(rnorm(Time * N, 0, 1), nrow = N, ncol = Time)
  
  ### Initialize Y and Ysl
  Yt = matrix(0, nrow = N, ncol = 2)
  Ytsl = W %*% Yt
  truebeta = list()
  
  for (t in 3:Time) {
    l = lfind(t, tseq)
    newYt = matrix(0, nrow = N)
    betat = matrix(0, N, pqd + 1)
    for (i in 1:N) {
      betat[i,1] = 0.1 + sign(Sj[i,1] * Sj[i,2]) * 0.0075
      betat[i,2] = fjlt(1, l, t, Time) * cjms(1, 1, Sj[i,1], Sj[i,2])
      betat[i,3] = fjlt(2, l, t, Time) * cjms(2, 1, Sj[i,1], Sj[i,2])
      betat[i,4] = fjlt(3, l, t, Time) * cjms(3, 1, Sj[i,1], Sj[i,2])
      betat[i,5] = fjlt(4, l, t, Time) * cjms(4, 1, Sj[i,1], Sj[i,2])
      
      newYt[i] = sum(c(1, Ytsl[i,t-1], Ytsl[i,t-2], Yt[i,t-1], X[i,t]) * betat[i,]) +
        rnorm(1, 0, 1) * 0.05
    }
    truebeta[[t]] = betat
    Yt = cbind(Yt, newYt)
    Ytsl = cbind(Ytsl, W %*% newYt)
  }
  
  ### Estimation
  knots = c(146, 292) / Time
  BS = ns(((r + 1):Time) / Time, knots = knots, Boundary.knots = c(r / Time, 441 / Time))
  K = ncol(BS)
  
  result = list()
  XX = list()
  
  for (posN in 1:N) {
    Xsi = matrix(0, nrow = Time - r, ncol = (1 + pqd) * K * (Time - r))
    for (t in (r+1):Time) {
      Zjsi = matrix(0, nrow = 1, ncol = (1 + pqd) * K)
      Zjsi[1:K] = BS[t - r,]
      Zjsi[(K+1):(2*K)] = Ytsl[posN, t-1] * BS[t - r,]
      Zjsi[(2*K+1):(3*K)] = Ytsl[posN, t-2] * BS[t - r,]
      Zjsi[(3*K+1):(4*K)] = Yt[posN, t-1] * BS[t - r,]
      Zjsi[(4*K+1):(5*K)] = X[posN, t] * BS[t - r,]
      
      for (j in 1:(t - r)) {
        Xsi[t - r, ((j - 1) * (1 + pqd) * K + 1):(j * (1 + pqd) * K)] = Zjsi
      }
    }
    
    best_lambda = cv.glmnet(Xsi, Yt[posN, (r+1):Time], alpha = 1, intercept = FALSE)$lambda.1se
    best_model = glmnet(Xsi, Yt[posN, (r+1):Time], alpha = 1, lambda = best_lambda, intercept = FALSE)
    hatgamma = coef(best_model)[-1]
    
    hatgamma2 = matrix(hatgamma, nrow = Time - r, byrow = TRUE)
    aa = rowSums(hatgamma2 != 0)
    result[[posN]] = which(aa != 0) + r
    XX[[posN]] = Xsi
  }
  
  cand = sort(unique(unlist(result)))
  cand = cand[cand > (r + C * (1 + pqd) * K)]
  cand = cand[cand < (Time - C * (1 + pqd) * K)]
  
  ### Step 2
  min_error = 0
  change_point = c()
  
  error_seq = sapply(cand, function(cand_item) detect(cand_item, c(), XX, Yt))
  error_seq = N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * (1 + pqd) * K * N
  a1 = cand[which.min(error_seq)]
  min_error = min(error_seq)
  change_point = c(change_point, a1)
  old_min_error = min_error
  LLL = 2
  
  while (LLL) {
    cand = setdiff(cand, (a1 - C * (1 + pqd) * K):(a1 + C * (1 + pqd) * K))
    error_seq = sapply(cand, function(cand_item) detect(cand_item, change_point, XX, Yt))
    error_seq = N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * LLL * (1 + pqd) * K * N
    if (min(error_seq) < old_min_error) {
      old_min_error = min(error_seq)
      a1 = cand[which.min(error_seq)]
      change_point = sort(c(change_point, a1))
      LLL = LLL + 1
    } else {
      break
    }
  }
  
  true_label1 = c(rep(1, 149), rep(2, 146), rep(3, 145))
  est_label = rep(1, 440)
  for (i in 1:(length(change_point) + 1)) {
    if (i == 1) {
      est_label[1:change_point[1]] = i
    } else if (i == length(change_point) + 1) {
      est_label[(change_point[length(change_point)] + 1):440] = i
    } else {
      est_label[(change_point[i - 1] + 1):change_point[i]] = i
    }
  }
  rand_index = rand.index(true_label1, est_label)
  
  ## Record run time
  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ## Return structured result
  list(
    change_point = change_point,
    rand_index = rand_index,
    run_time_seconds = run_time
  )
}


saveRDS(change_list, file = "48_440.rds")

stopCluster(cl)


################# Plot of Figure 1 #######################
change_list1=readRDS("24_220.rds")
change_list2=readRDS("24_440.rds")
change_list3=readRDS("48_220.rds")
change_list4=readRDS("48_440.rds")

## Case 1(a)
change_points <- lapply(change_list1, function(x) x$change_point)
table1=table(unlist(change_points))
data <- data.frame(table1)
colnames(data)=c("Number","Count")
complete_range <- data.frame(Number = 18:220)
complete_data <- merge(complete_range, data, by = "Number", all.x = TRUE)
complete_data$Count[is.na(complete_data$Count)] <- 0
table_format <- t(complete_data$Count)
colnames(table_format) <- complete_data$Number
table_matrix <- as.matrix(table_format)

png("Figure 1 a.png", width = 800, height = 600, res = 120)
barplot(table_matrix,main = "Case 1(a)",cex.main = 2,cex.axis = 1.2)
abline(v = 85,  col = "red", lwd = 2)
abline(v = 163, col = "red", lwd = 2)
dev.off()

## Case 1(c)
change_points <- lapply(change_list2, function(x) x$change_point)
table1=table(unlist(change_points))
data <- data.frame(table1)
colnames(data)=c("Number","Count")

complete_range <- data.frame(Number = 33:440)
complete_data <- merge(complete_range, data, by = "Number", all.x = TRUE)
complete_data$Count[is.na(complete_data$Count)] <- 0
table_format <- t(complete_data$Count)
colnames(table_format) <- complete_data$Number
table_matrix <- as.matrix(table_format)

png("Figure 1 c.png", width = 800, height = 600, res = 120)
barplot(table_matrix,main = "Case 1(c)",cex.main = 2,cex.axis = 1.2)
abline(v = 140,  col = "red", lwd = 2)
abline(v = 313, col = "red", lwd = 2)
dev.off()

## Case 1(b)
change_points <- lapply(change_list3, function(x) x$change_point)
table1=table(unlist(change_points))
data <- data.frame(table1)
colnames(data)=c("Number","Count")

complete_range <- data.frame(Number = 18:220)
complete_data <- merge(complete_range, data, by = "Number", all.x = TRUE)
complete_data$Count[is.na(complete_data$Count)] <- 0
table_format <- t(complete_data$Count)
colnames(table_format) <- complete_data$Number
table_matrix <- as.matrix(table_format)

png("Figure 1 b.png", width = 800, height = 600, res = 120)
barplot(table_matrix,main = "Case 1(b)",cex.main = 2,cex.axis = 1.2)
abline(v = 85,  col = "red", lwd = 2)
abline(v = 163, col = "red", lwd = 2)
dev.off()

## Case 1(d)
change_points <- lapply(change_list4, function(x) x$change_point)
table1=table(unlist(change_points))
data <- data.frame(table1)
colnames(data)=c("Number","Count")

complete_range <- data.frame(Number = 33:440)
complete_data <- merge(complete_range, data, by = "Number", all.x = TRUE)
complete_data$Count[is.na(complete_data$Count)] <- 0
table_format <- t(complete_data$Count)
colnames(table_format) <- complete_data$Number
table_matrix <- as.matrix(table_format)

png("Figure 1 d.png", width = 800, height = 600, res = 120)
barplot(table_matrix,main = "Case 1(d)",cex.main = 2,cex.axis = 1.2)
abline(v = 140,  col = "red", lwd = 2)
abline(v = 313, col = "red", lwd = 2)
dev.off()

################# Table 1 #######################
change_list1=readRDS("24_220.rds")
change_list2=readRDS("24_440.rds")
change_list3=readRDS("48_220.rds")
change_list4=readRDS("48_440.rds")

table_res=matrix(0,nrow=4,ncol=6)
colnames(table_res)=c("t1","t1+-3","t2","t2+-3","ARI","ARI_SD")

## Case 1(a)
check_l_r=function(x,bound){
  a1=88-bound
  b1=88+bound
  a2=154-bound
  b2=154+bound
  
  y1 <- any(x >= a1 & x <= b1)
  y2 <- any(x >= a2 & x <= b2)
  
  return(c(y1,y2))
}

change_points <- lapply(change_list1, function(x) x$change_point)
rd <- lapply(change_list1, function(x) x$rand_index)
res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,0))
}
table_res[1,c(1,3)]=colSums(res)/100

res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,3))
}
table_res[1,c(2,4)]=colSums(res)/100

table_res[1,c(5,6)]=c(mean(unlist(rd)),sd(unlist(rd)))

## Case 1(c)
check_l_r=function(x,bound){
  a1=148-bound
  b1=148+bound
  a2=294-bound
  b2=294+bound
  
  y1 <- any(x >= a1 & x <= b1)
  y2 <- any(x >= a2 & x <= b2)
  
  return(c(y1,y2))
}

change_points <- lapply(change_list2, function(x) x$change_point)
rd <- lapply(change_list2, function(x) x$rand_index)
res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,0))
}
table_res[3,c(1,3)]=colSums(res)/100

res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,3))
}
table_res[3,c(2,4)]=colSums(res)/100

table_res[3,c(5,6)]=c(mean(unlist(rd)),sd(unlist(rd)))

## Case 1(b)
check_l_r=function(x,bound){
  a1=88-bound
  b1=88+bound
  a2=154-bound
  b2=154+bound
  
  y1 <- any(x >= a1 & x <= b1)
  y2 <- any(x >= a2 & x <= b2)
  
  return(c(y1,y2))
}

change_points <- lapply(change_list3, function(x) x$change_point)
rd <- lapply(change_list3, function(x) x$rand_index)
res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,0))
}
table_res[2,c(1,3)]=colSums(res)/100

res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,3))
}
table_res[2,c(2,4)]=colSums(res)/100

table_res[2,c(5,6)]=c(mean(unlist(rd)),sd(unlist(rd)))

## Case 1(d)
check_l_r=function(x,bound){
  a1=148-bound
  b1=148+bound
  a2=294-bound
  b2=294+bound
  
  y1 <- any(x >= a1 & x <= b1)
  y2 <- any(x >= a2 & x <= b2)
  
  return(c(y1,y2))
}

change_points <- lapply(change_list4, function(x) x$change_point)
rd <- lapply(change_list4, function(x) x$rand_index)
res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,0))
}
table_res[4,c(1,3)]=colSums(res)/100

res=c()
for(i in 1:100){
  seq=change_points[[i]]
  res=rbind(res,check_l_r(seq,3))
}
table_res[4,c(2,4)]=colSums(res)/100

table_res[4,c(5,6)]=c(mean(unlist(rd)),sd(unlist(rd)))

write.csv(table_res,file="table1.csv",quote=F,row.names = F)

################# Appendix B table 1 #######################
change_list1=readRDS("24_220.rds")
change_list2=readRDS("24_440.rds")
change_list3=readRDS("48_220.rds")
change_list4=readRDS("48_440.rds")

table_res=matrix(0,nrow=1,ncol=8)
colnames(table_res)=c("Case 1a","Case 1a SE",
                      "Case 1b","Case 1b SE",
                      "Case 1c","Case 1c SE",
                      "Case 1d","Case 1d SE")

time1 <- lapply(change_list1, function(x) x$run_time_seconds)
table_res[1,1]=mean(unlist(time1))
table_res[1,2]=sd(unlist(time1))
time2 <- lapply(change_list3, function(x) x$run_time_seconds)
table_res[1,3]=mean(unlist(time2))
table_res[1,4]=sd(unlist(time2))
time3 <- lapply(change_list2, function(x) x$run_time_seconds)
table_res[1,5]=mean(unlist(time3))
table_res[1,6]=sd(unlist(time3))
time4 <- lapply(change_list4, function(x) x$run_time_seconds)
table_res[1,7]=mean(unlist(time4))
table_res[1,8]=sd(unlist(time4))

write.csv(table_res,file="tableB1.csv",quote=F,row.names = F)
