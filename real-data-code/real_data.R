library(arrow)
library(dplyr)
library(KRLS)
library(splines)
library(tidyr)
library(glmnet)
library(MASS)

##########################
detect=function(cand_item,change_point,XX,Yt){
  error=0
  a_cand=sort(c(change_point,cand_item))
  LL=length(a_cand)
  N=length(XX)
  for(posN in 1:N){
    ZXi=XX[[posN]][,1:(5*K)]
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

data=read.csv("data/combined_data.csv")
id_list=read.csv("data/location.csv")

start_time=Sys.time()
N=nrow(id_list) # number of locations
Time=12*9 # time
r=1
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

### pre-defined W
W=gausskernel(X = Sj,sigma=1)
W=W-diag(N)
W=t(apply(W, 1, function(X) X / sum(X)))
Ytsl=W%*%Y

## cubic B-spline
knots=c(36,72)/Time # knots for cubic spline
BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,(Time+1)/Time))
K=ncol(BS)
result=list()
XX=list()
est_par_list=list()
for (posN in 1:N) {
  Xsi=matrix(0,nrow=Time-r,ncol=5*K*(Time-r)) # initial X(si)
  est_par1=c()
  for (t in (r+1):Time) {
    Zjsi=matrix(0,nrow=1,ncol=5*K)
    Zjsi[1:K]=BS[t-r,]
    Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
    Zjsi[(2*K+1):(3*K)]=Y[posN,t-1]*BS[t-r,]
    Zjsi[(3*K+1):(4*K)]=X1[posN,t-1]*BS[t-r,]
    Zjsi[(4*K+1):(5*K)]=X2[posN,t-1]*BS[t-r,]
    
    for (j in 1:(t-r)) {
      Xsi[t-r,((j-1)*5*K+1):(j*5*K)]=Zjsi
    }
  }
  
  best_lambda=cv.glmnet(Xsi, Y[posN,(r+1):Time], alpha = 1
                        ,intercept = FALSE)$lambda.1se
  
  
  best_model=glmnet(Xsi, Y[posN,(r+1):Time], alpha = 1,
                    lambda =  best_lambda, intercept = FALSE)
  hatgamma=coef(best_model)[-1]
  hatgamma2=matrix(hatgamma,nrow=Time-r,byrow = TRUE)
  aa=rowSums(hatgamma2!=0) # count time with non-zero parameters
  
  result[[posN]]=which(aa!=0)+r
  XX[[posN]]=Xsi
  
  for (t in (r+1):Time) {
    tmp <- hatgamma2[t-r, ]
    BS_tmp <- BS[t-r, ]
    
    # Reshape tmp and multiply
    m <- matrix(tmp, ncol = length(BS_tmp), byrow = TRUE)
    est_par <- m %*% BS_tmp
    est_par1=cbind(est_par1,est_par)
  }
  
  est_par_list[[posN]]=est_par1
}

C=1
cand=sort(unique(unlist(result)))
cand=cand[cand>(r+C*3*K)]
cand=cand[cand<(Time-C*3*K)]

min_error=0
change_point=c()


error_seq=c()
for (cand_item in cand) {
  error_seq=c(error_seq,detect(cand_item,c(),XX,Y))
}
error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*5*K*N
a1=cand[which.min(error_seq)]
min_error=min(error_seq)
change_point=c(change_point,a1)

old_min_error=min_error
LLL=2

while(LLL){
  error_seq=c()
  cand=setdiff(cand,(a1-C*5*K):(a1+C*5*K))
  for (cand_item in cand) {
    error_seq=c(error_seq,detect(cand_item,change_point,XX,Y))
  }
  error_seq=N*(Time-r)*log(error_seq/N/(Time-r))+2*LLL*5*K*N
  if(min(error_seq)<old_min_error){
    old_min_error=min(error_seq)
    LLL=LLL+1
    a1=cand[which.min(error_seq)]
    change_point=sort(c(change_point,a1))
  } else{
    break
  }
}

## re-estimate
knots=c(36,72)/Time # knots for cubic spline
BS=ns(((r+1):Time)/Time,knots=knots, Boundary.knots = c(r/Time,(Time+1)/Time))
K=ncol(BS)

change_point_new=c(r,change_point,Time)
theta_re_est=list()

for (posN in 1:N) {
  theta_re_est1=matrix(0,nrow=length(change_point_new)-1,ncol=5*K)
  Xsi=matrix(0,nrow=Time-r,ncol=5*K) # initial X(si)
  for (t in (r+1):Time) {
    Zjsi=matrix(0,nrow=1,ncol=5*K)
    Zjsi[1:K]=BS[t-r,]
    Zjsi[(K+1):(2*K)]=Ytsl[posN,t-1]*BS[t-r,]
    Zjsi[(2*K+1):(3*K)]=Y[posN,t-1]*BS[t-r,]
    Zjsi[(3*K+1):(4*K)]=X1[posN,t-1]*BS[t-r,]
    Zjsi[(4*K+1):(5*K)]=X2[posN,t-1]*BS[t-r,]
    
    Xsi[t-r,]=Zjsi
  }
  
  for (t in 1:(length(change_point_new)-1)) {
    Xnew=Xsi[(change_point_new[t]+1-r):(change_point_new[t+1]-r),]
    Ytnew=Y[posN,(change_point_new[t]+1):change_point_new[t+1]]
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

end_time=Sys.time()
run_time=end_time-start_time
print(run_time)

################ Save data for Figure 5 and 6 ###################
write.csv(Dis,"Dis.csv", row.names = FALSE)

Y11=apply(Y[,(r+1):32],1,mean)
Y13=apply(Y[,33:108],1,mean)
Y14=cbind(Y11,Y13)
write.csv(Y14,"Dis2.csv", row.names = FALSE)

################# Plot for Figure 4 ########################
result1=matrix(0,nrow = N, ncol = 5)
result2=matrix(0,nrow = N, ncol = 5)
result3=matrix(0,nrow = N, ncol = 5)
for(posN in 1:N){
  tmp=est_par_list[[posN]]
  tmp1=tmp[,1:(15-r)]
  result1[posN,]=apply(tmp1, 1, mean)
  tmp2=tmp[,(16-r):(32-r)]
  result2[posN,]=apply(tmp2, 1, mean)
  tmp3=tmp[,(33-r):(Time-r)]
  result3[posN,]=apply(tmp3, 1, mean)
}

normalize_matrix <- function(mat) {
  (mat - mean(mat)) / sd(mat)
}

result1 <- normalize_matrix(result1)
result2 <- normalize_matrix(result2)
result3 <- normalize_matrix(result3)

library(ggplot2)
library(reshape2)
library(patchwork)

melt_data1 <- melt(result1)
melt_data2 <- melt(result2)
melt_data3 <- melt(result3)

melt_data1$label1="Jan 2011-March 2012"
melt_data2$label1="April 2012-August 2013"
melt_data3$label1="September 2013-December 2019"
melt_data=rbind(melt_data1,melt_data2,melt_data3)
melt_data=as.data.frame(melt_data)

melt_data[melt_data$value<0,]$value=-1
melt_data[melt_data$value>5,]$value=6
melt_data$label1 <- factor(melt_data$label1, levels = c("Jan 2011-March 2012", "April 2012-August 2013", "September 2013-December 2019"))

ggplot(melt_data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  facet_wrap(~ label1, scales = "free_y") +
  scale_fill_gradient2(low="#edf8b1",mid="#7fcdbb", high="#2c7fb8", midpoint=2.5,breaks=c(0,2.5,5),
                       labels=c("0","2.5",
                                "5")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "5 Factors", y = "Locations")

ggsave("Figure 4.png", width = 8, height = 6, dpi = 300)
