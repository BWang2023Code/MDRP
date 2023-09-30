#simulate missing values under MAR assumption
MARmissing=function(X_sample,idx.incomplete,idx.covariates,idx.weights,perc.missing,seed){
  n=dim(X_sample)[1]
  p_pred=sum(idx.covariates==1)
  findBeta_MAR=function(beta_MAR,X_sample,p_pred,idx.weights,idx.covariates){
    weight=as.matrix(rep(1/p_pred,p_pred)*idx.weights[idx.weights!=0],p_pred,1)*beta_MAR
    X.com=X_sample[,idx.covariates==1]
    X.com=apply(as.data.frame(X.com),2, function(x){as.numeric(x)})
    eta=exp(X.com%*%weight)
    p.miss=eta/(1+eta)
    indm_MAR_j=NULL
    for (j in 1:sum(idx.incomplete)) {
      indm_MAR_j=cbind(indm_MAR_j,rbinom(n, 1, p.miss))
    }
    return(mean(rowSums(indm_MAR_j)>0))
  }
  B_MAR=seq(-10,10,0.01)
  beta_MAR=mean(B_MAR[which(abs(sapply(B_MAR, findBeta_MAR,X_sample=X_sample,p_pred=p_pred,idx.weights=idx.weights,idx.covariates=idx.covariates)-perc.missing)<0.025)])
  weight=as.matrix(rep(1/p_pred,p_pred)*idx.weights[idx.weights!=0],p_pred,1)*beta_MAR
  X.com=X_sample[,idx.covariates==1]
  X.com=apply(as.data.frame(X.com),2, function(x){as.numeric(x)})
  eta=exp(X.com%*%weight)
  p.miss=eta/(1+eta)
  if(sum(idx.incomplete)==1){
    X_sample[which(rbinom(n, 1, p.miss)==1),idx.incomplete==1]=NA
  }else{
    for (j in 1:sum(idx.incomplete)) {
      indm_MAR_j=which(rbinom(n, 1, p.miss)==1)
      X_sample[indm_MAR_j,idx.incomplete==1][,j]=NA
    }
  }
  return(X_sample)
}
#simulate missing values under NMAR assumptions
NMARmissing=function(X_sample,idx.incomplete,idx.covariates,idx.weights,perc.missing,seed){
  set.seed(seed)
  n=dim(X_sample)[1]
  p_pred=sum(idx.covariates==1)
  findBeta_NMAR=function(beta_NMAR,X_sample,p_pred,idx.weights,idx.covariates){
    weight=rep(1/p_pred,p_pred)*idx.weights[idx.weights!=0]
    x.simu_mnar=apply(X_sample[,idx.covariates==1], 2, function(x){as.numeric(x)})
    colnames(x.simu_mnar)=colnames(X_sample[,idx.covariates==1])
    latent_correlation_matrix <- toeplitz(c(1, 0.2))
    x <-rep(x.simu_mnar%*%as.matrix(weight),each=sum(idx.incomplete==1))
    simulated_binary_dataset <- rbin(clsize = sum(idx.incomplete==1), intercepts = 0,
                                     betas = beta_NMAR, xformula = ~x, 
                                     cor.matrix = latent_correlation_matrix,link = "probit")
    res=simulated_binary_dataset[["Ysim"]]
    return(mean(rowSums(res)>0))
  }
  B_NMAR=seq(-3,3,0.01)
  beta_NMAR=mean(B_NMAR[which(abs(sapply(B_NMAR, findBeta_NMAR,X_sample=X_sample,p_pred=p_pred,idx.weights=idx.weights,idx.covariates=idx.covariates)-perc.missing)<0.05)])
  weight=rep(1/p_pred,p_pred)*idx.weights[idx.weights!=0]
  x.simu_mnar=apply(X_sample[,idx.covariates==1], 2, function(x){as.numeric(x)})
  colnames(x.simu_mnar)=colnames(X_sample[,idx.covariates==1])
  latent_correlation_matrix <- toeplitz(c(1, 0.2))
  x <-rep(x.simu_mnar%*%as.matrix(weight),each=sum(idx.incomplete==1))
  simulated_binary_dataset <- rbin(clsize = sum(idx.incomplete==1), intercepts = 0,
                                   betas = beta_NMAR, xformula = ~x, 
                                   cor.matrix = latent_correlation_matrix,link = "probit")
  res=simulated_binary_dataset[["Ysim"]]
  if(sum(idx.incomplete)==1){
    X_sample[which(res==1),idx.incomplete==1]=NA
  }else{
    for (j in 1:sum(idx.incomplete)) {
      X_sample[which(res[,j]==1),idx.incomplete==1][,j]=NA
    }
  }
  l=list()
  l[[1]]=X_sample
  l[[2]]=weight
  l[[3]]=beta_NMAR
  names(l)=c("X_sample","weight","beta_NMAR")
  return(l)
}
#Calculate relative risk
relativerisk=function(i,data_X){
  #LogHR = c(1.4, 1.2, 0.18, 0.36, 0.1, 0.2, -0.05, -0.075, -0.15, -0.55, -0.45, 0.2, 0.325, 0.3, 0.15, 0.275, 0.22)
  #names(LogHR)=c("famhist2", "parity2", "age_1stb2", "age_1stb3", "height2", "height3", "bmi2", "bmi3", "bmi4", "HRT2", "HRT3", "bmi2HRT2", "bmi3HRT2", "bmi4HRT2", "bmi2HRT3", "bmi3HRT3", "bmi4HRT3" )
  woman_i=data_X[i,]
  if(dim(data_X)[2]==1){names(woman_i)=names(data_X)}
  para_i=paste0(names(woman_i),woman_i)
  para_i=c(para_i,paste0(para_i[which(substr(names(woman_i),1,3)=="bmi")],para_i[which(substr(names(woman_i),1,3)=="HRT")]))
  rr_i=exp(sum(LogHR[names(LogHR)%in%para_i]))
  return(rr_i)
}
#Calculate absolute risk, projection interval of length 10
AR_rr=function(RR,mort_inc){
  a=10
  time=mort_inc[,1]
  lambda2time=mort_inc[,2]
  lambda1time=t(outer((lambda*gamma*time^(gamma-1)),RR))
  IntgrlLngth=1
  AbsRisk <- rep(NA,length(RR))
  for (i in 1:length(RR)){
    RskWrk <- 0
    Cumlambda <- 0
    StrtIntvl=which(time>50)[1]-1
    lambda1itime <-lambda1time[i,]
    NumbrIntvl=which(time>50+a)[1]-which(time>50)[1]
    rr_i=RR[i]
    for (j in 1:NumbrIntvl){
      jintvl <- StrtIntvl+j-1
      lambdaij <- lambda1itime[jintvl]+lambda2time[jintvl]
      Cumlambda <- Cumlambda+lambdaij*IntgrlLngth
      PIj <- ((lambda1itime[jintvl]/lambdaij)*exp(-Cumlambda))*(1-exp(-lambdaij*IntgrlLngth))
      RskWrk <- RskWrk+PIj 
    } 
    AbsRisk[i] <- 100*RskWrk 
  }
  return(AbsRisk)
}
#Calculate conditional absolute risk under MAR setting
CondAR_MAR=function(i,X_miss,X_popu,mort_inc){
  p=6
  if(sum(is.na(X_miss[i,]))==0){
    AR_i=AR_rr(relativerisk(1,X_miss[i,]),mort_inc)
  }else{
    X_i=X_miss[i,]
    miss_pattern_i=which(is.na(X_i))
    obs_pattern_i=which(!is.na(X_i))
    Xobs_i=X_i[obs_pattern_i]
    Xobs_ref_i=X_popu[,obs_pattern_i]
    Xmis_ref_i=X_popu[,miss_pattern_i]
    id_cond_i=which(apply(Xobs_ref_i, 1, function(x){identical(as.numeric(x),as.numeric(Xobs_i))}))
    ncond_i=length(id_cond_i)
    if(length(miss_pattern_i)==1){
      tt=table(Xmis_ref_i[id_cond_i])
      Cond_i=cbind(as.numeric(names(tt)),as.numeric(tt),tt/sum(tt))
      colnames(Cond_i)=c(colnames(X_popu)[miss_pattern_i],"n","prob")
    }else{
      Cond_i=Xmis_ref_i[id_cond_i,] %>% group_by_all %>% count
      Cond_i[,"prob"]=Cond_i[,"n"]/ncond_i
    }
    Cond_i=as.data.frame(cbind(Xobs_i %>% slice(rep(1:n(), each = dim(Cond_i)[1])),data.frame(Cond_i)))
    Cond_i[,"RR"]=sapply(1:dim(Cond_i)[1],relativerisk,data=Cond_i)
    Cond_i[,"AR"]=AR_rr(Cond_i[,"RR"],mort_inc)
    AR_i=sum(Cond_i[,"prob"]*Cond_i[,"AR"])
  }
  return(AR_i)
}

#estimate conditional absolute risk under MAR setting based on reference sample
CondARhat_MAR=function(i,X_miss,X_ref,mort_inc,n.imp){
  p=6
  if(sum(is.na(X_miss[i,]))==0){
    AR_i=AR_rr(relativerisk(1,X_miss[i,]),mort_inc)
  }else{
    X_i=X_miss[i,]
    n_ref=dim(X_ref)[1]
    miss_pattern_i=which(is.na(X_i))
    obs_pattern_i=which(!is.na(X_i))
    Xobs_i=X_i[obs_pattern_i]
    Xobs_ref_i=X_ref[,obs_pattern_i]
    Xmis_ref_i=X_ref[,miss_pattern_i]
    id_cond_i=which(apply(Xobs_ref_i, 1, function(x){identical(as.numeric(x),as.numeric(Xobs_i))}))
    if(length(id_cond_i)==0){
      R_obsi=log(relativerisk(1,Xobs_i))
      R_obsi_ref=log(sapply(1:n_ref, relativerisk, data=Xobs_ref_i))
      if(length(miss_pattern_i)==1){
        Xmis_ref_i=as.data.frame(Xmis_ref_i)
        colnames(Xmis_ref_i)=names(X_ref)[miss_pattern_i]
      }
      R_mis_ref=log(sapply(1:n_ref, relativerisk, data=Xmis_ref_i))
      R_i_ref=cbind(R_obsi_ref,R_mis_ref)
      R_i_ref_knn=cbind(R_obsi_ref,R_mis_ref,dist=abs(R_obsi_ref-R_obsi))
      R_i_ref_knn=R_i_ref_knn[order(R_i_ref_knn[,"dist"],decreasing=FALSE),]
      R_i_ref_knn=R_i_ref_knn[1:n.imp,]
      rr_hat_knn=exp(R_i_ref_knn[,"R_mis_ref"]+R_obsi)
      AR_i=mean(sapply(rr_hat_knn,AR_rr,mort_inc=mort_inc))
    }else{
      ncond_i=length(id_cond_i)
      if(length(miss_pattern_i)==1){
        tt=table(Xmis_ref_i[id_cond_i])
        Cond_i=cbind(as.numeric(names(tt)),as.numeric(tt),tt/sum(tt))
        colnames(Cond_i)=c(colnames(X_ref)[miss_pattern_i],"n","prob")
      }else{
        Cond_i=Xmis_ref_i[id_cond_i,] %>% group_by_all %>% count
        Cond_i[,"prob"]=Cond_i[,"n"]/ncond_i
      }
      Cond_i=as.data.frame(cbind(Xobs_i %>% slice(rep(1:n(), each = dim(Cond_i)[1])),data.frame(Cond_i)))
      Cond_i[,"RR"]=sapply(1:dim(Cond_i)[1],relativerisk,data=Cond_i)
      Cond_i[,"AR"]=AR_rr(Cond_i[,"RR"],mort_inc)
      AR_i=sum(Cond_i[,"prob"]*Cond_i[,"AR"])
    }
  }
  return(AR_i)
}


#Calculate conditional absolute risk under NMAR setting
CondAR_NMAR=function(i,X_miss,X_popu,idx.covariates,idx.incomplete,weight,mort_inc,joint_miss_NMAR){
  p=6
  n_popu=dim(X_popu)[1]
  if(sum(is.na(X_miss[i,]))==0){
    AR_i=AR_rr(relativerisk(1,X_miss[i,]),mort_inc)
  }else{
    X_i=X_miss[i,]
    name_miss=colnames(X_i[idx.incomplete==1])
    name_pred=colnames(X_i[idx.covariates==1])
    miss_pattern_i=which(is.na(X_i))
    obs_pattern_i=which(!is.na(X_i))
    Xobs_i=X_i[obs_pattern_i]
    Xobs_ref_i=X_popu[,obs_pattern_i]
    Xmis_ref_i=X_popu[,miss_pattern_i]
    id_cond_i=which(apply(as.data.frame(Xobs_ref_i), 1, function(x){identical(as.numeric(x),as.numeric(Xobs_i))}))
    
    if(length(miss_pattern_i)==1){
      tt=table(Xmis_ref_i[id_cond_i])
      Cond_i=cbind(as.numeric(names(tt)),as.numeric(tt),tt/n_popu)
      colnames(Cond_i)=c(colnames(X_popu)[miss_pattern_i],"n","prob")
    }else{
      Cond_i=Xmis_ref_i[id_cond_i,] %>% group_by_all %>% count
      Cond_i[,"prob"]=Cond_i[,"n"]/n_popu
    }
    Cond_i=as.data.frame(cbind(Xobs_i %>% slice(rep(1:n(), each = dim(Cond_i)[1])),data.frame(Cond_i)))
    Cond_i[,"RR"]=sapply(1:dim(Cond_i)[1],relativerisk,data=Cond_i)
    Cond_i[,"AR"]=AR_rr(Cond_i[,"RR"],mort_inc)
    joint_m_i=ifelse(identical(as.numeric(miss_pattern_i),c(1,2)),1,ifelse(as.numeric(miss_pattern_i)==1,2,ifelse(as.numeric(miss_pattern_i)==2,3,4)))
    Cond_i[,"P.miss"]=merge(Cond_i[,name_pred],joint_miss_NMAR,by=name_pred)[,as.character(joint_m_i)]
    Cond_i[,"CondP"]=Cond_i[,"prob"]*Cond_i[,"P.miss"]/sum(Cond_i[,"prob"]*Cond_i[,"P.miss"])
    AR_i=sum(Cond_i[,"CondP"]*Cond_i[,"AR"])
  }
  return(AR_i)
}


#Using risk score method to estimate conditional absolute risk under MAR setting based on reference sample
CondAR_RiskScore=function(i,X_miss,X_ref,mort_inc,n.imp){
  p=6
  X_i=X_miss[i,]
  n_ref=dim(X_ref)[1]
  miss_pattern_i=which(is.na(X_i))
  obs_pattern_i=which(!is.na(X_i))
  Xobs_i=X_i[obs_pattern_i]
  Xobs_ref_i=X_ref[,obs_pattern_i]
  Xmis_ref_i=X_ref[,miss_pattern_i]
  
  R_obsi=log(relativerisk(1,Xobs_i))
  R_obsi_ref=log(sapply(1:n_ref, relativerisk, data=Xobs_ref_i))
  
  
  if(length(miss_pattern_i)==1){
    Xmis_ref_i=as.data.frame(Xmis_ref_i)
    colnames(Xmis_ref_i)=names(X_ref)[miss_pattern_i]
  }
  R_mis_ref=log(sapply(1:n_ref, relativerisk, data=Xmis_ref_i))
  R_i_ref=cbind(R_obsi_ref,R_mis_ref)
  
  R_i_ref=R_i_ref[order(R_i_ref[,"R_obsi_ref"],decreasing=FALSE),]
  Qstrata=quantile(as.numeric(R_i_ref[,"R_obsi_ref"]),probs = seq(0,1,0.1))
  Strata_ind_i=ifelse(R_obsi>=max(Qstrata),length(Qstrata)-1,ifelse(R_obsi<=min(Qstrata),1,which(R_obsi<=Qstrata)[1]-1))
  Rmis_srata_i=as.numeric(R_i_ref[which(R_i_ref[,"R_obsi_ref"]>=Qstrata[Strata_ind_i] & R_i_ref[,"R_obsi_ref"]<=Qstrata[Strata_ind_i+1]),"R_mis_ref"])
  rr_hat_HD=exp(Rmis_srata_i[sample(1:length(Rmis_srata_i),n.imp,replace = T)]+R_obsi)
  #########################################################################
  R_i_ref_knn=cbind(R_obsi_ref,R_mis_ref,dist=abs(R_obsi_ref-R_obsi))
  R_i_ref_knn=R_i_ref_knn[order(R_i_ref_knn[,"dist"],decreasing=FALSE),]
  R_i_ref_knn=R_i_ref_knn[1:n.imp,]
  rr_hat_knn=exp(R_i_ref_knn[,"R_mis_ref"]+R_obsi)
  #########################################################################
  lm.Ri <- lm(R_i_ref[,"R_mis_ref"] ~ R_i_ref[,"R_obsi_ref"])
  rr_hat_CERS=exp(as.numeric(unlist(lm.Ri["coefficients"])%*%c(1,R_obsi)+R_obsi))
  rr_hat_CDRS=exp(as.numeric(rnorm(n.imp,unlist(lm.Ri["coefficients"])%*%c(1,R_obsi),sigma(lm.Ri))+R_obsi))
  #########################################################################
  CAR_hat_HotD_Rs_i=mean(sapply(rr_hat_HD,AR_rr,mort_inc=mort_inc))
  CAR_hat_KNN_Rs_i=mean(sapply(rr_hat_knn,AR_rr,mort_inc=mort_inc))
  CAR_hat_condexp_CERS_i=AR_rr(rr_hat_CERS,mort_inc)
  CAR_hat_condexp_CDRS_i=mean(sapply(rr_hat_CDRS,AR_rr,mort_inc=mort_inc))
  l=c(CAR_hat_condexp_CERS_i,CAR_hat_HotD_Rs_i,CAR_hat_KNN_Rs_i,CAR_hat_condexp_CDRS_i)
  return(l)
  
}
#Using Lowest Level method to estimate conditional absolute risk
LowestLevelImputation=function(i,X_miss){
  p=6
  Risk_factor_names=c("famhist","parity","age_1stb","height","bmi","HRT")
  Risk_factor_names_1=Risk_factor_names[1:4]
  Lowestrisklevel_1=c(1,1,1,1)
  Risk_factor_names_2=Risk_factor_names[5:6]
  X_i=X_miss[i,]
  miss_pattern_i=which(is.na(X_i))
  if(length(miss_pattern_i)==1){
    if(sum(Risk_factor_names[miss_pattern_i]%in%Risk_factor_names_2)>0){
      X_i[miss_pattern_i]=ifelse(X_i[,'HRT']==1,4,ifelse(X_i[,'HRT']==2,1,1))
    }else{
      X_i[miss_pattern_i]=Lowestrisklevel_1[miss_pattern_i]
    }
  }else{
    X_i[which(1:4%in%miss_pattern_i)]=Lowestrisklevel_1[which(1:4%in%miss_pattern_i)]
    if(sum(is.na(X_i))>0){
      X_i[is.na(X_i)]=ifelse(X_i[,'HRT']==1,4,ifelse(X_i[,'HRT']==2,1,1))
    }
  }
  return(as.factor(X_i))
}
#multiple imputation method using MICE based on reference dataset with only new patient added
MICEbyi=function(i,X_miss,X_ref,mort_inc,n.imp){
  if(sum(is.na(X_miss[i,]))==0){
    AR_miss_hat_MICE_i20=AR_miss_hat_MICE_i50=AR_miss_hat_MICE_i=AR_rr(RR=relativerisk(1,X_miss[i,]),mort_inc=mort_inc)
  }else{
    M=n.imp
    X.miss_mice_i=rbind(X_miss[i,],X_ref)
    init_i = mice(X.miss_mice_i, maxit=0) 
    meth_i = init_i[["method"]]
    predM_i = init_i[["predictorMatrix"]]
    meth_i[c("famhist","parity")]="logreg"
    imputed_Data_MICE_i <- mice(X.miss_mice_i, m=M, maxit = 20, method=meth_i, predictorMatrix=predM_i, seed = sample(1:1000,1))
    AR_hat_mice_m_i=NULL
    mice.imputed.data_i=list()
    for (m in 1:M) {
      X.imputed.m_i=complete(imputed_Data_MICE_i,m)[1,]
      mice.imputed.data_i[[m]]=X.imputed.m_i
      RR_hat_m_i=sapply(1, relativerisk, data=X.imputed.m_i)
      AR_hat_mice_m_i=cbind(AR_hat_mice_m_i,AR_rr(RR=RR_hat_m_i,mort_inc=mort_inc))
      AR_miss_hat_MICE_i=mean(AR_hat_mice_m_i)
      
    }

  }
  return(AR_miss_hat_MICE_i)
}
