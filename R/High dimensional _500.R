#Loading packages
library(quantreg)
library(gam)
library(mvtnorm)
require(Matrix) 
require(metaSEM)
require(matrixcalc)
require(expm)
require(evmix)
require(lme4)
require(REEMtree) #this package needs a session restart to work
library(caTools)
require(MixRF)
library(LongituRF)
library(lqmm)
library(quantregForest)


pc<- proc.time()


# ------------------------------------------------------------------
# Scenario specification
# ------------------------------------------------------------------
# 1 = Normal intercept / Normal error
# 2 = Normal intercept / ALD error
# 3 = Normal intercept / t error
# 4 = ALD intercept / Normal error
# 5 = ALD intercept / ALD error
# 6 = ALD intercept / t error

scenario_id <- 1   # change to 1,2,3,4,5, or 6

scenario_list <- list(
  list(intercept = "normal", error = "normal"),
  list(intercept = "normal", error = "ald"),
  list(intercept = "normal", error = "t"),
  list(intercept = "ald",    error = "normal"),
  list(intercept = "ald",    error = "ald"),
  list(intercept = "ald",    error = "t")
)

current_scenario <- scenario_list[[scenario_id]]

draw_random_intercept <- function(n, dist) {
  if (dist == "normal") {
    return(rnorm(n, mean = 0, sd = 1))
  }
  if (dist == "ald") {
    return(ald::rALD(n, mu = 0, sigma = 1, p = 0.5))
  }
  stop("Unsupported intercept distribution: ", dist)
}

draw_error <- function(n, dist) {
  if (dist == "normal") {
    return(rnorm(n, mean = 0, sd = 1))
  }
  if (dist == "ald") {
    return(ald::rALD(n, mu = 0, sigma = 1, p = 0.5))
  }
  if (dist == "t") {
    return(rt(n, df = 4))
  }
  stop("Unsupported error distribution: ", dist)
}

r = 1000 #no.of iterations
n = 500 #no.of subjects
m = 6 #no.of timepoints
p=0.8 #proportion of train subset
tao = 0.5 #quantile
t=5 #no.of variables to start with the trees
m1=list()
fx_eff_lmm=matrix(nrow = r,ncol=25)
ranef_lmm=list()
MSE1=matrix(nrow=r,ncol = 3)
colnames(MSE1)=c("Train","Test","Test2")
predict_LMM=matrix(nrow=r,ncol = n*round((p*m)))
predict2_LMM=matrix(nrow=r,ncol = n*round(p*m))

R1=matrix(nrow=r,ncol=2)
REEMresult=list()
ranef_REEM=list()
predict_REEM=matrix(nrow=r,ncol=n*round(p*m))
predict2_REEM=matrix(nrow=r,ncol=n*round(p*m))

MSE2=matrix(nrow=r,ncol = 3)
colnames(MSE2)=c("Train","Test","Test2")
prediction_error1=matrix(nrow=r,ncol = n*round(p*m))
ME=list()
ranef_ME=list()
MSE4=matrix(nrow=r,ncol = 3)
colnames(MSE4)=c("Train","Test","Test2")

predict_ME=matrix(nrow=r,ncol =n*round(p*m))
predict2_ME=matrix(nrow=r,ncol =n* round(p*m))

Reemforest=list()
ranef_Reemforest=list()
predict_Reemforest=matrix(nrow=r,ncol = n*round(p*m))
predict2_Reemforest=matrix(nrow=r,ncol = n*round(p*m))

MSE5=matrix(nrow=r,ncol = 3)
colnames(MSE5)=c("Train","Test","Test2")

qlmm=list()
ranef_qlmm=list()
predict_lqmm=matrix(nrow=r,ncol = n*round(p*m))
predict2_lqmm=matrix(nrow=r,ncol = n*round(p*m))

MSE7=matrix(nrow=r,ncol = 3)
colnames(MSE7)=c("Train","Test","Test2")

qf=list()
predict_qf=matrix(nrow=r,ncol = n*round(p*m))
predict2_qf=matrix(nrow=r,ncol = n*round(p*m))

MSE9=matrix(nrow=r,ncol = 3)
colnames(MSE9)=c("Train","Test","Test2")

ID = rep(1:n, each = m)
X1 = rep(c(0,1), each = m*n/2, length.out = m*n) ##treatment group
## X2 = rep(rgamma(100, 30, 1), length.out = 10) ##BMI averaged
#setseed to fix the vectors but each time we should run the set.seed line together with the vector
X2 = rgamma(m*n, 30, 1) ##BMI averaged yet varied across the same subject
X3 = rep(1:m, length.out = m*n) ##timepoints
X4 = 0.25*X2+rnorm(m*n, mean = 0, sd = 1) ##correlated with X2
X5 = rt(m*n,10,0)
X6 = rnorm(m*n,10,5)
X7 = runif(m*n,5,10)
X8 = rt(m*n,25,5)
X9 = rnorm(n*m,100,20)
X10 = rchisq(n*m,15,5) 
X11 = cos(X7)+rnorm(m*n,mean=0,sd=1)
X12 = rlnorm(m*n,mean=0,sd=1)
Z_1 = runif(m*n, 0, 1) 
Z_2 = rnorm(m*n, mean = 30, sd = 4.2)
Z_3 = 2*Z_1 + rnorm(m*n, mean = 0, sd = 2) ##correlated with Z_1
Z_4 = log(Z_2)+rnorm(m*n, mean = 0, sd = 1) ##correlated with Z_2
Z_5 = rchisq(m*n,10,0)
Z_6 = rnorm(m*n, mean = 0, sd = 1)
Z_7 = rt(m*n,20,0)
Z_8 = rexp(m*n,12)  
Z_9 = Z_5+0.25*Z_6+rnorm(m*n, mean = 0, sd = 1)
Z_10 = 2*Z_1 +rnorm(m*n, mean = 0, sd = 2)
Z_11 = X5+2*X6+rnorm(m*n, mean = 0, sd = 1)
Z_12 = log(Z_9)+runif(m*n,0,1)
d.f = data.frame(ID,1,X1,X2,X3,X4,Z_1,Z_2,Z_3,Z_4,
                 X5,X6,X7,X8,X9,X10,X11,X12,Z_5,Z_6,Z_7,Z_8,Z_9,Z_10,Z_11,Z_12)
b = as.vector(c(10,5,3.5,1.25,5,7.5,2,-2.5,-0.25,
                1,4,8,6,4.2,-3,-0.25,2.5,1.75,4.25,-4,-1,1.5,-0.45,10,1.75))

v = 2.94 #v: variance of errors (homoscedastic)
# Defining the covariance matrix for subject i
sigma = v

nu = 4
# mu=20
# mu0=15

mat=Matrix::Diagonal(m)

sdvec=10
rho=0.1
#Defining the covariance matrix for all subjects
sigma_whole = metaSEM::bdiagRep(v*mat,n)
# Generate random effects first

u_i<- draw_random_intercept(n, current_scenario$intercept)
u=rep(u_i,each=m)
X = d.f[,-1]
data=as.data.frame(cbind(X,ID))

results=list()
test_Y=list()
test2_Y=list()
test=list()
train=list()
train_Y=list()
datalist_train=list()
datalist_test=list()
predict_LMM1=matrix(nrow=r,ncol = n*round((p*m)))
predict_REEM1=matrix(nrow=r,ncol=n*round(p*m))
predict_ME1=matrix(nrow=r,ncol =n*round(p*m))
predict_lqmm1=matrix(nrow=r,ncol = n*round(p*m))
predict_qf1=matrix(nrow=r,ncol = n*round(p*m))
NN1=matrix(ncol = 5, nrow = r)

ranef_REEM1=list()
ranef_ME1=list()
ranef_qlmm1=list()
ranef_lmm1=list()

for (i in 1:r) {
  drop(data)
  data=as.data.frame(cbind(X,ID))
  
  print(sprintf("Iteration: %d", i)) ######### display number of iteration
  flush.console()
  Sys.sleep(1)
  
  #Data 
  options(scipen = 999)
  
  e <- draw_error(n * m, current_scenario$error)
  e1 = matrix(t(e), ncol = 1, byrow = TRUE)##t() to get the transpose of the matrix
  
  
  # Generating y|r conditional distribution
  Y = as.matrix(X) %*% b + u + e1
  data_=cbind(data,u,e)
  e11=subset(data_$e,X3<=round(p*m))
  u11=subset(data_$u,X3<=round(p*m))
  #is.positive.definite(sigma_whole)
  data=as.data.frame(cbind(data,Y))
  names(data)[names(data) == "X1.1"] <- "Intercept"
  data=data[-1]
  
  train_ID  <- subset(data$ID, X3<=round(p*m))
  test_ID   <- subset(data$ID, X3>round(p*m))
  
  train_Y[[i]]<-subset(data$Y, X3<=round(p*m))
  test_Y[[i]]<-subset(data$Y, X3>round(p*m))
  
  
  test[[i]]=subset(data,X3>round(p*m))
  train[[i]]=subset(data,X3<=round(p*m))
  #test2=as.data.frame(cbind(test2_X1,test2_X2,test2_X3,test2_Z1,test2_Z2,test2_ID,test2_Y[[i]]))
  #train2=as.data.frame(cbind(train2_X1,train2_X2,train2_X3,train2_Z1,train2_Z2,train2_ID,train2_Y))
  
  #LMM
  m1[[i]]=lmer(Y~1+X1+X2+X3+X4+Z_1+Z_2+Z_3+Z_4+
                X5+X6+X7+X8+X9+X10+X11+X12+Z_5+Z_6+Z_7+Z_8+Z_9+Z_10+Z_11+Z_12 
               +(1|ID), data = train[[i]], 
               REML = TRUE)
  #m2[[i]]=lmer(Y~1+X1+X2+X3+Z_1+(1|ID), data = train2, REML = TRUE)
  # summary(m1) #package doesn't compute p-values
  
  fx_eff_lmm[i,]=fixef(m1[[i]])
  
  ranef_lmm[[i]]=as.matrix(unlist(ranef(m1[[i]])$ID))
  predict_LMM[i,]=predict(m1[[i]],newdata=test[[i]])
  #predict2_LMM[i,]=predict(m2[[i]],newdata=test2)
  
  MSE1[i,1]=mean(residuals(m1[[i]])^2)
  MSE1[i,2]=mean((predict_LMM[i,]-test_Y[[i]])^2)
  #MSE1[i,3]=mean((predict2_LMM[i,]-test2_Y[[i]])^2)
  
  R1[i,]=MuMIn::r.squaredGLMM(m1[[i]]) #extract rsquared m:marginal and c:conditional
  
  
  #REEM tree
  REEMresult[[i]]<-REEMtree::REEMtree(Y~1+X1+X2+X3+X4+Z_1+Z_2+Z_3+Z_4+
                                        X5+X6+X7+X8+X9+X10+X11+X12+
                                        Z_5+Z_6+Z_7+Z_8+Z_9+Z_10+Z_11+Z_12, 
                                      data=train[[i]], 
                                      random= ~1|ID)
  # print(REEMresult)
  # plot(REEMresult)
  ranef_REEM[[i]]=ranef.REEMtree(REEMresult[[i]])
  test_n = purrr::map_dfr(seq_len(nrow(train[[i]])/nrow(test[[i]])),
                          ~test[[i]]) 
  predict_REEM[i,]=predict.REEMtree(REEMresult[[i]],
                                    newdata = test_n) #the number of observations should be
  #predict2_REEM[i,]=predict.REEMtree(REEMresult[[i]],newdata = test2) #the number of observations should be
  
  #the same in the training and testing datasets
  MSE2[i,1] = mean((REEMresult[[1]][["residuals"]])^2)
  #residuals_test = as.matrix(predict_REEM[i,]-test_Y[[i]])
  MSE2[i,2] =mean((predict_REEM[i,]-test_Y[[i]])^2)
  
  #MSE2[i,3] =mean((predict2_REEM[i,]-test2_Y[[i]])^2)
  
  prediction_error1[i,]=
    (1/((m/2)*n))*(test_Y[[i]]-predict_REEM[i,])^2
  
  data=as.data.frame(cbind(X,ID,Y))
  names(data)[names(data) == "X1.1"] <- "Intercept"
  test_X=subset(data[1:25],X3>round(p*m))
  train_X=subset(data[1:25],X3<=round(p*m))
  
  train_Y[[i]]<-subset(data$Y, X3<=round(p*m))
  test_Y[[i]]<-subset(data$Y, X3>round(p*m))
  
  sample = sample.split(data, SplitRatio = p)
  train2_X = subset(data[1:25], sample == TRUE)
  test2_X  = subset(data[1:25], sample == FALSE)
  
  train2_ID  <- subset(data$ID, sample == TRUE)
  test2_ID   <- subset(data$ID, sample == FALSE)
  
  train2_Y<-subset(data$Y, sample == TRUE)
  test2_Y[[i]]<-subset(data$Y, sample == FALSE)
  #library(LongituRF)
  
  data_train=as.data.frame(cbind(train_X,train_ID,
                                 train_Y[[i]]))
  names(data_train)[names(data_train) == "X1.1"] <- "Intercept"
  datalist_train[[i]]=list(as.matrix(train_X[2:25]),
                           train_ID,
                           as.matrix(train_X[1]),
                           train_Y[[i]])
  
  data_test=as.data.frame(cbind(test_X,test_ID,test_Y[[i]]))
  names(data_test)[names(data_test) == "X1.1"] <- "Intercept"
  datalist_test[[i]]=list(as.matrix(test_X[2:25]),
                          test_ID,as.matrix(test_X[1]),test_Y[[i]])
  
  data_test2=as.data.frame(cbind(test2_X,test2_ID,test2_Y[[i]]))
  names(data_test2)[names(data_test2) == "X1.1"] <- "Intercept"
  datalist_test2=list(as.matrix(test2_X[2:25]),
                      test2_ID,as.matrix(test2_X[1]),test2_Y[[i]])
  
  
  #ME random forest
  ME[[i]]=MERF(X=datalist_train[[i]][[1]],Y=datalist_train[[i]][[4]],
               id=datalist_train[[i]][[2]],Z=datalist_train[[i]][[3]],
               time=datalist_train[[i]][[2]], sto="none",mtry=t)
  
  
  ranef_ME[[i]]=ME[[i]]$random_effects
  
  predict_ME[i,]=predict(ME[[i]],X=datalist_test[[i]][[1]],Y=datalist_test[[i]][[4]],
                         id=datalist_test[[i]][[2]],Z=datalist_test[[i]][[3]],
                         time=datalist_test[[i]][[2]])
  # predict2_ME[i,]=predict(ME[[i]],X=datalist_test2[[1]],Y=datalist_test2[[4]],
  #                         id=datalist_test2[[2]],Z=datalist_test2[[3]],
  #                         time=datalist_test2[[2]])
  MSE4[i,1]=ME[[i]]$forest$mse[ME[[i]]$forest$ntree]
  MSE4[i,2]=mean((predict_ME[i,]-test_Y[[i]])^2)
  #MSE4[i,3]=mean((predict2_ME[i,]-test2_Y[[i]])^2)
  
  #Quantile regression LMM 
  qlmm[[i]]=lqmm(Y~1+X1+X2+X3+X4+Z_1+Z_2+Z_3+Z_4+
                   X5+X6+X7+X8+X9+X10+X11+X12+Z_5+Z_6+Z_7+Z_8+Z_9+Z_10+Z_11+Z_12,~1,group=ID,
                 data=train[[i]],tau=0.5)
  #summary(qlmm)
  #qlmm$theta_x
  predict_lqmm[i,]=predict(qlmm[[i]],newdata = test[[i]])
  #predict2_lqmm[i,]=predict(qlmm[[i]],newdata = test2)
  
  ranef_qlmm[[i]]=ranef.lqmm(qlmm[[i]]) #random effects
  # MSE6=mean(resid(qlmm)[,1]^2) #0.1 quantile
  # MSE6
  MSE7[i,1]=mean(resid(qlmm[[i]])^2) #0.5 quantile
  MSE7[i,2]=mean((predict_lqmm[i,]-test_Y[[i]])^2)
  #MSE7[i,3]=mean((predict2_lqmm[i,]-test2_Y[[i]])^2)
  
  # MSE7
  # MSE8=mean(resid(qlmm)[,3]^2) #0.9 quantile
  # MSE8
  
  #Quantile Regression Forests
  
  qf[[i]]=quantregForest(x=train[[i]][1:24],y=train[[i]]$Y ,
                         data=train[[i]],
                         keep.inbag = T,mtry=5)
  # qf$rsq
  predict_qf[i,]=predict(qf[[i]],newdata = test[[i]],what=0.5)
  #predict2_qf[i,]=predict(qf[[i]],newdata = test2,what=0.5)
  
  MSE9[i,1] =qf[[i]]$mse[qf[[i]]$forest$ntree]
  MSE9[i,2]=mean((predict_qf[i,]-test_Y[[i]])^2)
  #MSE9[i,3]=mean((predict2_qf[i,]-test2_Y[[i]])^2)
  
  # MSE9
  
  #bias_qf=(test_Y-predict_qf)^2
  
  results[[i]]=matrix(nrow=1,ncol=7)
  colnames(results[[i]])=c("LMM","REEMtree",
                           "Mixed effect random forest",
                           "REEM forest","QLMM",
                           "Quantile Random Forest",
                           "Quantile Gradient boost")
  
  rownames(results[[i]])=c("MSE")
  results[[i]][1,]=c(MSE1[i,2],MSE2[i,2],
                     MSE4[i,2],"",
                     MSE7[i,2],MSE9[i,2],"")
  
  predict_LMM1[i,]=predict(m1[[i]],newdata=train[[i]])
  predict_REEM1[i,]=predict.REEMtree(REEMresult[[i]],
                                     newdata = train[[i]])
  predict_ME1[i,]=predict(ME[[i]],X=datalist_train[[i]][[1]],Y=datalist_train[[i]][[4]],
                          id=datalist_train[[i]][[2]],Z=datalist_train[[i]][[3]],
                          time=datalist_train[[i]][[2]])
  predict_lqmm1[i,]=predict(qlmm[[i]],newdata = train[[i]])
  predict_qf1[i,]=predict(qf[[i]],newdata = train[[i]],what=0.5)
  
  ranef_lmm1[[i]]= matrix(rep(as.numeric(t(ranef_lmm[[i]])), each = round(p*m)), 
                          ncol = ncol(ranef_lmm[[i]]), byrow = TRUE)
  ranef_REEM1[[i]]=matrix(rep(as.numeric(t(ranef_REEM[[i]])), 
                              each = round(p*m)), 
                          ncol = ncol(ranef_REEM[[i]]), byrow = TRUE)
  ranef_ME1[[i]]=matrix(rep(as.numeric(t(ranef_ME[[i]])), 
                            each = round(p*m)), ncol = ncol(ranef_ME[[i]]), byrow = TRUE)
  ranef_qlmm1[[i]]=matrix(rep(as.numeric(t(ranef_qlmm[[i]])), 
                              each = round(p*m)), ncol = ncol(ranef_qlmm[[i]]), byrow = TRUE)
  
}

proc.time() - pc
NN=matrix(ncol = 5, nrow = r)
for (i in 1:r){
  for (j in 1:(round((1-p)*m)*n)){
    NN[i,1] =sum((predict_LMM[i,j]-test_Y[[i]][j])^2)
    NN[i,2] =sum((predict_REEM[i,j]-test_Y[[i]][j])^2)
    NN[i,3] =sum((predict_ME[i,j]-test_Y[[i]][j])^2)
    NN[i,4] =sum((predict_lqmm[i,j]-test_Y[[i]][j])^2)
    NN[i,5] =sum((predict_qf[i,j]-test_Y[[i]][j])^2)
    
  }
}
NN=as.data.frame(NN)
colnames(NN)=c("LMM","REEM tree","MERF","QLMM","QRF")
MSE_results=cbind(mean(NN[,1]),mean(NN[,2]),mean(NN[,3]),mean(NN[,4]),mean(NN[,5]))
for (i in 1:r){
  for (j in 1:(round(p*m)*n)){
    fixed[i,1] =sum(((predict_LMM[i,j]-ranef_lmm1[[i]][j])-
                     (train_Y[[i]][j]-e11[j]-u11[j]))^2)
    fixed[i,2] =sum(((predict_REEM[i,j]-ranef_REEM1[[i]][j,])-
                     (train_Y[[i]][j]-e11[j]-u11[j]))^2)
    fixed[i,3] =sum(((predict_ME[i,j]-ranef_ME1[[i]][j])-
                     (train_Y[[i]][j]-e11[j]-u11[j]))^2)
    fixed[i,4] =sum(((predict_lqmm[i,j]-ranef_qlmm1[[i]][j,])-
                     (train_Y[[i]][j]-e11[j]-u11[j]))^2)
    fixed[i,5] =sum(((predict_qf[i,j])-
                     (train_Y[[i]][j]-e11[j]-u11[j]))^2)
    
  }}

fixed=as.data.frame(fixed)
colnames(fixed)=c("LMM","REEM tree","MERF","QLMM","QRF")
fixedMSE_results=cbind(mean(fixed[,1]),mean(fixed[,2]),mean(fixed[,3]),mean(fixed[,4]),mean(fixed[,5]))

abss=matrix(ncol = 5, nrow = r)
for (i in 1:r){
  for (j in 1:(round((1-p)*m)*n)){
    abss[i,1] =sum(abs(predict_LMM[i,j]-test_Y[[i]][j]))
    abss[i,2] =sum(abs(predict_REEM[i,j]-test_Y[[i]][j]))
    abss[i,3] =sum(abs(predict_ME[i,j]-test_Y[[i]][j]))
    abss[i,4] =sum(abs(predict_lqmm[i,j]-test_Y[[i]][j]))
    abss[i,5] =sum(abs(predict_qf[i,j]-test_Y[[i]][j]))
    
  }
}
abss=as.data.frame(abss)
colnames(abss)=c("LMM","REEM tree","MERF","QLMM","QRF")
MAE_results=cbind(mean(abss[,1]),mean(abss[,2]),mean(abss[,3]),mean(abss[,4]),mean(abss[,5]))

write.csv(MSE_results,"scenario1_h500_6_MSE_results.csv")
write.csv(fixedMSE_results,"scenario1_h500_6_fixedMSE_results.csv")
write.csv(MAE_results,"scenario1_h500_6_MAE_results.csv")
colnames(MAE_results)=c("LMM","REEM tree","MERF","QLMM","QRF")
