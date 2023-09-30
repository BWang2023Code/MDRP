#Code for project "Predicting absolute risk for a person with missing risk factors"
#load population dataset, age-specific competing mortality rates and log of hazard ratio
load("~MDRP dataload.RData")
#load functions
source("~MDRP functions Github.r")
seed=1
set.seed(seed)
#dimension
p=6
n=500
#number of imputation
n.imp=70
#population size
n_popu=dim(ref_X_population)[1]
#smaple a complete random sample from populaiton
ID.sample=sample(1:n_popu,n,replace = FALSE)
X_sample=ref_X_population[ID.sample,2:(1+p)]
######################################################################
#The  breast cancer hazard used in the model 
lambda=0.45*10^(-4)
gamma=2.15
time=mort_inc[,1]
lambda2_time=mort_inc[,2]
#calculate the true relative risk and the absolute risk of women with complete infromation
RR_complete=sapply(1:n,relativerisk,data_X=X_sample)
AR_sample.complete=AR_rr(RR_complete,mort_inc)
#specify missing mechanism and generate missing values under different assumptions: MAR and NMAR 1 to 3

#MAR setting
perc.missing = 0.8
idx.incomplete.MAR=c(1,1,0,0,0,0)
idx.covariates.MAR=c(0,0,1,1,1,1)
idx.weights.MAR=c(0,0,1,1,1,1)
mechanism="MAR"
X_miss_MAR=MARmissing(X_sample,idx.incomplete=idx.incomplete.MAR,idx.covariates=idx.covariates.MAR,idx.weights=idx.weights.MAR,perc.missing,seed=seed)
X_sample.miss=X_miss_MAR
##################################################
#NMAR_id=1 to 3 refer to setting NMAR 1 to 3
library(SimCorMultRes)
NMAR_id=1
perc.missing = 0.8
idx.incomplete.NMAR=c(1,1,0,0,0,0)
idx.covariates.NMAR=c(1,1,1,0,0,0)
idx.weights.NMAR=ifelse(NMAR_id==1,c(1,1,1,0,0,0),ifelse(NMAR_id==2,c(1,1,-1,0,0,0),c(2,-1,-1,0,0,0)))
X_miss_NMAR=NMARmissing(X_sample,idx.incomplete=idx.incomplete.NMAR,idx.covariates=idx.covariates.NMAR,idx.weights=idx.weights.NMAR,perc.missing,seed=seed)
X_sample.miss=X_miss_NMAR[["X_sample"]]

#sample of reference sample with complete data from population
n_ref=5000
X_ref_data=ref_X_population[sample(1:n_popu,n_ref,replace = FALSE),]
rownames(X_ref_data)=1:n_ref

####################################################################################################  
####################################################################################################  
####################################################################################################  
####################################################################################################  
library(dplyr)
#If under MAR assumption
ID.miss=which(rowSums(is.na(X_sample.miss))>0)
CAR_sample.miss_MAR_popu=sapply(1:n,CondAR_MAR,X_miss=X_sample.miss,X_popu=ref_X_population[,2:(1+p)],mort_inc=mort_inc)
##################################################  
CARhat_sample.miss_MAR_m=sapply(1:n,CondARhat_MAR,X_miss=X_sample.miss,X_ref=X_ref_data[,2:(1+p)],mort_inc=mort_inc,n.imp=n.imp)
##################################################
#Gail's imputation_lowest risk group
X.imputed.Gail=t(sapply(1:n,LowestLevelImputation, X_miss=X_sample.miss))
RR_Gail=sapply(1:n,relativerisk,data_X=X.imputed.Gail)
AR_hat_Gail=AR_rr(RR_Gail,mort_inc)
##################################################
#Risk score (iCARE and others, assuming MAR)
ConAR_hat_RS=Reduce(rbind,lapply(1:n,CondAR_RiskScore,X_miss=X_sample.miss,X_ref=X_ref_data[,2:(1+p)],mort_inc=mort_inc,n.imp=n.imp))
ConAR_hat_CERS=ConAR_hat_RS[,1]
ConAR_hat_HDRS_m=ConAR_hat_RS[,2]
ConAR_hat_KNNRS_m=ConAR_hat_RS[,3]
ConAR_hat_CDRS_m=ConAR_hat_RS[,4]
################################################
#multiple imputation based on reference dataset with all new patients added
library(mice)
X.miss_mice=rbind(X_sample.miss,X_ref_data[,2:(1+p)])
init = mice(X.miss_mice, maxit=0) 
meth = init[["method"]]
predM = init[["predictorMatrix"]]
meth[c("famhist","parity")]="logreg"
imputed_Data_MICE <- mice(X.miss_mice, m=n.imp, maxit = 20, method=meth, predictorMatrix=predM, seed = seed)
AR_hat_mice_m=NULL
mice.imputed.data=list()
for (m in 1:n.imp) {
  X.imputed.m=complete(imputed_Data_MICE,m)[1:n,]
  mice.imputed.data[[m]]=X.imputed.m
  RR_hat_m=sapply(1:n, relativerisk, data=X.imputed.m)
  AR_hat_mice_m=cbind(AR_hat_mice_m,AR_rr(RR=RR_hat_m,mort_inc=mort_inc))
  
}
AR_miss_hat_MICE=rowMeans(AR_hat_mice_m)
####################################################################################################
#multiple imputation method using MICE based on reference dataset with only new patient added
MICEbyi_result=sapply(1:n, MICEbyi,X_miss=X_sample.miss,X_ref=X_ref_data[,2:(1+p)],mort_inc=mort_inc,n.imp=n.imp)
####################################################################################################  
#results summary
true_pred=cbind(AR_sample.complete,CAR_sample.miss_MAR_popu,AR_hat_Gail,ConAR_hat_CERS,AR_miss_hat_MICE,MICEbyi_result,ConAR_hat_HDRS_m,ConAR_hat_KNNRS_m,ConAR_hat_CDRS_m,CARhat_sample.miss_MAR_m)
colnames(true_pred)=c("AR","CAR","LRL","CERS","MICE","MICEi","HDRS","KNNRS","CDRS","CARhat")
result=round(t(apply(cbind(AR_hat_Gail,ConAR_hat_CERS,AR_miss_hat_MICE,MICEbyi_result,ConAR_hat_HDRS_m,ConAR_hat_KNNRS_m,ConAR_hat_CDRS_m,CARhat_sample.miss_MAR_m),2,function(x){return(c(mean((x-CAR_sample.miss_MAR_popu)),mean(abs(x-CAR_sample.miss_MAR_popu)),mean((x-CAR_sample.miss_MAR_popu)^2)))})),2)
rownames(result)=c("LRL","CERS","MICE","MICEi","HDRS","KNNRS","CDRS","CARhat") 
colnames(result)=c("Bias","aBias","MSE") 

####################################################################################################  
####################################################################################################  
####################################################################################################


####################################################################################################  
####################################################################################################  
#If under any of the NMAR assumption
ID.miss=which(rowSums(is.na(X_sample.miss))>0)
#joint distribution is used to calculate the true absolute risk under NMAR assumption
p_pred_NMAR=sum(idx.covariates.NMAR==1)
name.pred_NMAR=colnames(X_sample[idx.covariates.NMAR==1])
weight_NMAR=X_miss_NMAR[["weight"]]
names(weight_NMAR)=name.pred_NMAR
BETA_NMAR=X_miss_NMAR[["beta_NMAR"]]
X.com_NMAR=ref_X_population[,-1][,idx.covariates.NMAR==1]
X.com_NMAR=apply(X.com_NMAR, 2, function(x){as.numeric(x)})
colnames(X.com_NMAR)=name.pred_NMAR
latent_correlation_matrix <- toeplitz(c(1, 0.2))
x_NMAR <-rep(X.com_NMAR%*%as.matrix(weight_NMAR),each=sum(idx.incomplete.NMAR==1))
simulated_binary_dataset_NMAR <- rbin(clsize = sum(idx.incomplete.NMAR==1), intercepts = 0,
                                      betas = BETA_NMAR, xformula = ~x_NMAR, 
                                      cor.matrix = latent_correlation_matrix,link = "probit")
res_NMAR=simulated_binary_dataset_NMAR[["Ysim"]]
colnames(res_NMAR)=c("V1","V2")
X.com_NMAR=as.data.frame(cbind(X.com_NMAR,res_NMAR))
X.com_NMAR[,"M"]=ifelse(X.com_NMAR[,"V1"]==1 & X.com_NMAR[,"V2"]==1,1,ifelse(X.com_NMAR[,"V1"]==1 & X.com_NMAR[,"V2"]==0,2,ifelse(X.com_NMAR[,"V1"]==0 & X.com_NMAR[,"V2"]==1,3,4)))
X.com_NMAR=X.com_NMAR[,-c(4:5)]
joint_miss_NMAR=aggregate(M~famhist+parity+age_1stb, data=X.com_NMAR, FUN=function(x)  c(table(x), count=length(x)))
joint_miss_NMAR=cbind(joint_miss_NMAR[,1:sum(idx.covariates.NMAR)],joint_miss_NMAR[,"M"][,-5]/joint_miss_NMAR[,"M"][,"count"])
#################################################
ID.miss=which(rowSums(is.na(X_sample.miss))>0)
CAR_sample.miss_NMAR_popu=sapply(1:n,CondAR_NMAR,X_miss=X_sample.miss,X_popu=ref_X_population[,2:(1+p)],idx.covariates=idx.covariates.NMAR,idx.incomplete=idx.incomplete.NMAR,weight=weight_NMAR,mort_inc=mort_inc,joint_miss_NMAR=joint_miss_NMAR)
#################################################
CARhat_sample.miss_MAR=sapply(1:n,CondARhat_MAR,X_miss=X_sample.miss,X_ref=X_ref_data[,2:(1+p)],mort_inc=mort_inc,n.imp=n.imp)
##################################################
#lowest risk group method
X.imputed.Gail=t(sapply(1:n,LowestLevelImputation, X_miss=X_sample.miss))
RR_Gail=sapply(1:n,relativerisk,data_X=X.imputed.Gail)
AR_hat_Gail=AR_rr(RR_Gail,mort_inc)
##################################################
#Risk score (iCARE and others, assuming MAR)
ConAR_hat_RS=Reduce(rbind,lapply(1:n,CondAR_RiskScore,X_miss=X_sample.miss,X_ref=X_ref_data[,2:(1+p)],mort_inc=mort_inc,n.imp=n.imp))
ConAR_hat_CERS=ConAR_hat_RS[,1]
ConAR_hat_HDRS_m=ConAR_hat_RS[,2]
ConAR_hat_KNNRS_m=ConAR_hat_RS[,3]
ConAR_hat_CDRS_m=ConAR_hat_RS[,4]
##################################################
#multiple imputation based on reference dataset with all new patients added
library(mice)
X.miss_mice=rbind(X_sample.miss,X_ref_data[,2:(1+p)])
init = mice(X.miss_mice, maxit=0) 
meth = init[["method"]]
predM = init[["predictorMatrix"]]
meth[c("famhist","parity")]="logreg"
imputed_Data_MICE <- mice(X.miss_mice, m=n.imp, maxit = 20, method=meth, predictorMatrix=predM, seed = seed)
AR_hat_mice_m=NULL
mice.imputed.data=list()
for (m in 1:n.imp) {
  X.imputed.m=complete(imputed_Data_MICE,m)[1:n,]
  mice.imputed.data[[m]]=X.imputed.m
  RR_hat_m=sapply(1:n, relativerisk, data=X.imputed.m)
  AR_hat_mice_m=cbind(AR_hat_mice_m,AR_rr(RR=RR_hat_m,mort_inc=mort_inc))
  AR_miss_hat_MICE=rowMeans(AR_hat_mice_m)
  
}
####################################################################################################
#multiple imputation method using MICE based on reference dataset with only new patient added
MICEbyi_result=sapply(1:n, MICEbyi,X_miss=X_sample.miss,X_ref=X_ref_data[,2:(1+p)],mort_inc=mort_inc,n.imp=n.imp)
####################################################################################################  
#results summary
true_pred=cbind(AR_sample.complete,CAR_sample.miss_NMAR_popu,AR_hat_Gail,ConAR_hat_CERS,AR_miss_hat_MICE,MICEbyi_result,ConAR_hat_HDRS_m,ConAR_hat_KNNRS_m,ConAR_hat_CDRS_m,CARhat_sample.miss_MAR_m)
colnames(true_pred)=colnames(true_pred)=c("AR","CAR","LRL","CERS","MICE","MICEi","HDRS","KNNRS","CDRS","CARhat")
result=round(t(apply(cbind(AR_hat_Gail,ConAR_hat_CERS,AR_miss_hat_MICE,MICEbyi_result,ConAR_hat_HDRS_m,ConAR_hat_KNNRS_m,ConAR_hat_CDRS_m,CARhat_sample.miss_MAR_m),2,function(x){return(c(mean((x-CAR_sample.miss_NMAR_popu)),mean(abs(x-CAR_sample.miss_NMAR_popu)),mean((x-CAR_sample.miss_NMAR_popu)^2)))})),2)
rownames(result)=c("LRL","CERS","MICE","MICEi","HDRS","KNNRS","CDRS","CARhat") 
colnames(result)=c("Bias","aBias","MSE") 

