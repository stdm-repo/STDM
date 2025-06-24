library(dplyr)
library(splines)
library(tidyr)
library(glmnet)
library(MASS)
library(KRLS)

##########################
detect=function(cand_item,change_point,K,p,q,r,XX,Yt){
  error=0
  a_cand=sort(c(change_point,cand_item))
  LL=length(a_cand)
  N=length(XX)
  
  estimate_theta <- function(X_seg, Y_seg) {
    XtX_inv <- ginv(t(X_seg) %*% X_seg)
    theta <- XtX_inv %*% t(X_seg) %*% Y_seg
    residuals <- Y_seg - t(theta) %*% t(X_seg)
    return(sum(residuals^2))
  }
  
  
  
  for (posN in 1:N) {
    ZXi <- XX[[posN]][, 1:((p + q + 3) * K)]
    for (i in 1:(LL + 1)) {
      
      if (i == 1) {
        idx_X <- 1:(a_cand[1] - r)
        idx_Y <- (r + 1):a_cand[1]
        
      } else if (i == (LL + 1)) {
        idx_X <- (a_cand[LL] - r + 1):nrow(ZXi)
        idx_Y <- (a_cand[LL] + 1):ncol(Yt)
        
      } else {
        idx_X <- (a_cand[i - 1] - r + 1):(a_cand[i] - r)
        idx_Y <- (a_cand[i - 1] + 1):a_cand[i]
      }
      
      error <- error + estimate_theta(ZXi[idx_X, ], Yt[posN, idx_Y])
    }
  }
  
  return(error)
}


compute_Dis <- function(Sj, theta_re_est, N) {
  Dis <- numeric(N)
  thred <- 2 / sqrt(N)
  
  for (posN in 1:N) {
    ds <- Sj[posN, ]
    A1 <- A2 <- A3 <- A4 <- 0
    cand1 <- cand2 <- cand3 <- cand4 <- c()
    
    for (posN1 in setdiff(1:N, posN)) {
      if (Sj[posN1, 1] - ds[1] <= thred && Sj[posN1, 1] - ds[1] >= 0) {
        if (Sj[posN1, 2] - ds[2] <= thred && Sj[posN1, 2] - ds[2] >= 0) {
          cand1 <- c(cand1, posN1)
          A1 <- A1 + theta_re_est[[posN1]]
        }
        if (ds[2] - Sj[posN1, 2] <= thred && ds[2] - Sj[posN1, 2] >= 0) {
          cand3 <- c(cand3, posN1)
          A3 <- A3 + theta_re_est[[posN1]]
        }
      }
      if (ds[1] - Sj[posN1, 1] <= thred && ds[1] - Sj[posN1, 1] >= 0) {
        if (Sj[posN1, 2] - ds[2] <= thred && Sj[posN1, 2] - ds[2] >= 0) {
          cand2 <- c(cand2, posN1)
          A2 <- A2 + theta_re_est[[posN1]]
        }
        if (ds[2] - Sj[posN1, 2] <= thred && ds[2] - Sj[posN1, 2] >= 0) {
          cand4 <- c(cand4, posN1)
          A4 <- A4 + theta_re_est[[posN1]]
        }
      }
    }
    
    A1 <- A1 / ifelse(length(cand1) == 0, 1, length(cand1))
    A2 <- A2 / ifelse(length(cand2) == 0, 1, length(cand2))
    A3 <- A3 / ifelse(length(cand3) == 0, 1, length(cand3))
    A4 <- A4 / ifelse(length(cand4) == 0, 1, length(cand4))
    
    Dis[posN] <- sum((A1 - A2)^2) + sum((A2 - A3)^2) + sum((A3 - A4)^2) + sum((A4 - A1)^2)
  }
  
  return(Dis)
}


est_dis <- function(p, q, r, Y, Sj, X1, X2) {
  # Step 1: Construct kernel weights and transform Y
  W <- gausskernel(X = Sj, sigma = 1)
  W <- W - diag(N)
  W <- t(apply(W, 1, function(X) X / sum(X)))  # Row-normalize the weight matrix
  Ytsl <- W %*% Y  # Transformed response using spatial smoothing
  
  # Step 2: Construct B-spline basis for time
  knots <- c(36, 72) / Time
  BS <- ns(((r + 1):Time) / Time, knots = knots, Boundary.knots = c(r / Time, (Time + 1) / Time))
  K <- ncol(BS)
  
  # Step 3: Fit model for each location
  result <- list()
  XX <- list()
  est_par_list <- list()
  est_par_list_partial <- list()
  for (posN in 1:N) {
    Xsi <- matrix(0, nrow = Time - r, ncol = (p + q + 3) * K * (Time - r))
    
    for (t in (r + 1):Time) {
      Zjsi <- numeric((p + q + 3) * K)
      Zjsi[1:K] <- BS[t - r, ]
      
      for (pt in 1:p) {
        Zjsi[(pt * K + 1):((pt + 1) * K)] <- Ytsl[posN, t - pt] * BS[t - r, ]
      }
      for (qt in 1:q) {
        Zjsi[((p + qt) * K + 1):((p + qt + 1) * K)] <- Y[posN, t - qt] * BS[t - r, ]
      }
      
      Zjsi[((p + q + 1) * K + 1):((p + q + 2) * K)] <- X1[posN, t - 1] * BS[t - r, ]
      Zjsi[((p + q + 2) * K + 1):((p + q + 3) * K)] <- X2[posN, t - 1] * BS[t - r, ]
      
      for (j in 1:(t - r)) {
        Xsi[t - r, ((j - 1) * (p + q + 3) * K + 1):(j * (p + q + 3) * K)] <- Zjsi
      }
    }
    
    best_lambda <- cv.glmnet(Xsi, Y[posN, (r + 1):Time], alpha = 1, nfolds = 5, intercept = FALSE)$lambda.1se
    best_model <- glmnet(Xsi, Y[posN, (r + 1):Time], alpha = 1, lambda = best_lambda, intercept = FALSE)
    hatgamma <- coef(best_model)[-1]
    hatgamma2 <- matrix(hatgamma, nrow = Time - r, byrow = TRUE)
    
    aa <- rowSums(hatgamma2 != 0)
    result[[posN]] <- which(aa != 0) + r
    XX[[posN]] <- Xsi
    
    est_par1 <- NULL
    for (t in (r + 1):Time) {
      tmp <- hatgamma2[t - r, ]
      m <- matrix(tmp, ncol = length(BS[t - r, ]), byrow = TRUE)
      est_par <- m %*% BS[t - r, ]
      est_par1 <- cbind(est_par1, est_par)
    }
    est_par_list[[posN]] <- est_par1
    est_par_list_partial[[posN]] <- est_par1[,(15-r):(32-r)]
  }
  
  ini_Dis <- compute_Dis(Sj, est_par_list, N)
  ini_Dis_partial <- compute_Dis(Sj, est_par_list_partial, N)
  
  # Step 4: Change point detection using AIC
  cand <- sort(unique(unlist(result)))
  cand <- cand[cand > (r + 3 * K) & cand < (Time - 3 * K)]
  error_seq <- sapply(cand, function(cand_item) detect(cand_item, c(), K, p, q, r, XX, Y))
  error_seq <- N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * (p + q + 3) * K * N
  a1 <- cand[which.min(error_seq)]
  min_error <- min(error_seq)
  change_point <- a1
  
  old_min_error <- min_error
  LLL <- 1
  
  while (LLL) {
    cand <- setdiff(cand, (a1 - (p + q + 3) * K):(a1 + (p + q + 3) * K))
    error_seq <- sapply(cand, function(cand_item) detect(cand_item, change_point, K, p, q, r, XX, Y))
    error_seq <- N * (Time - r) * log(error_seq / N / (Time - r)) + 2 * (LLL + 1) * (p + q + 3) * K * N
    
    if (min(error_seq) < old_min_error) {
      old_min_error <- min(error_seq)
      LLL <- LLL + 1
      a1 <- cand[which.min(error_seq)]
      change_point <- sort(c(change_point, a1))
    } else {
      break
    }
  }
  
  # Step 5: Re-estimate parameters segment-wise
  change_point_new <- c(r, change_point, Time)
  theta_re_est <- list()
  theta_direct <- list()
  theta_re_est_partial <- list()
  theta_direct_partial <- list()
  
  for (posN in 1:N) {
    theta_re_est1 <- matrix(0, nrow = length(change_point_new) - 1, ncol = (p + q + 3) * K)
    Xsi <- matrix(0, nrow = Time - r, ncol = (p + q + 3) * K)
    
    for (t in (r + 1):Time) {
      Zjsi <- numeric((p + q + 3) * K)
      Zjsi[1:K] <- BS[t - r, ]
      
      for (pt in 1:p) {
        Zjsi[(pt * K + 1):((pt + 1) * K)] <- Ytsl[posN, t - pt] * BS[t - r, ]
      }
      for (qt in 1:q) {
        Zjsi[((p + qt) * K + 1):((p + qt + 1) * K)] <- Y[posN, t - qt] * BS[t - r, ]
      }
      
      Zjsi[((p + q + 1) * K + 1):((p + q + 2) * K)] <- X1[posN, t - 1] * BS[t - r, ]
      Zjsi[((p + q + 2) * K + 1):((p + q + 3) * K)] <- X2[posN, t - 1] * BS[t - r, ]
      
      Xsi[t - r, ] <- Zjsi
    }
    
    for (t in 1:(length(change_point_new) - 1)) {
      Xnew <- Xsi[(change_point_new[t] + 1 - r):(change_point_new[t + 1] - r), ]
      Ytnew <- as.numeric(Y[posN, (change_point_new[t] + 1):change_point_new[t + 1]])
      best_lambda2 <- cv.glmnet(Xnew, Ytnew, alpha = 0, intercept = FALSE)$lambda.1se
      best_model2 <- glmnet(Xnew, t(Ytnew), alpha = 0, lambda = best_lambda2, intercept = FALSE)
      theta_re_est1[t, ] <- coef(best_model2)[-1]
    }
    theta_re_est1[is.na(theta_re_est1)] <- 0
    theta_re_est[[posN]] <- theta_re_est1
    theta_re_est_partial[[posN]] <- theta_re_est1[2,]
    
    # direct regression
    Xsi=matrix(0,nrow=Time-r,ncol=p+q+2) # initial X(si)
    for (t in (r+1):Time) {
      Zjsi=matrix(0,nrow=1,ncol=p+q+2)
      Zjsi=c(Ytsl[posN,(t-1):(t-p)],Y[posN,(t-1):(t-q)],X1[posN,t],X2[posN,t])
      
      Xsi[t-r,]=Zjsi
    }
    est_par2=lm(matrix(t(Y[posN,(r+1):Time]),nrow=Time-r,ncol=1)~Xsi)
    theta_direct[[posN]] <- coef(est_par2)[-1]
    
  }
  
  Dis_direct <- compute_Dis(Sj, theta_direct, N)
  Dis <- compute_Dis(Sj, theta_re_est, N)
  Dis_partial <- compute_Dis(Sj, theta_re_est_partial, N)
  
  
  return(list(ini_Dis=ini_Dis, ini_Dis_partial=ini_Dis_partial,
              Dis_direct=Dis_direct,Dis_partial=Dis_partial,
              time_change_point=change_point,Dis=Dis))
}


##########################
data=read.csv("data/combined_data.csv")
id_list=read.csv("data/location.csv")
set.seed(123)
N=nrow(id_list) # number of locations
Time=12*9 # time

Sj=id_list[,2:3]
Sj=scale(Sj)
X1=data[,c(1,5,7)]
X2=data[,c(1,6,7)]
Y1=data[,c(1,4,7)]
Y=Y1 %>%
  pivot_wider( names_from = time, 
               values_from = count,
               values_fill = list(count = 0))
X1=X1 %>%
  pivot_wider( names_from = time, 
               values_from = avg_customer,
               values_fill = list(avg_customer = 0))
X2=X2 %>%
  pivot_wider( names_from = time, 
               values_from = avg_mile,
               values_fill = list(avg_mile = 0))
Y=as.matrix(Y[,2:ncol(Y)])
X1=as.matrix(X1)
X2=as.matrix(X2)

tmp=est_dis(1,1,1,Y,Sj,X1,X2)

write.csv(tmp$Dis,"Dis_our.csv", row.names = FALSE)
write.csv(tmp$Dis_direct,"Dis_structure.csv", row.names = FALSE)
write.csv(tmp$ini_Dis,"Dis_dynamic.csv", row.names = FALSE)
