DML_HDMA=function(x,m,y,z,K){
  
  n=nrow(m)
  p=ncol(m)
  q=ncol(z)
  fold=split(1:n,cut(1:n,K))
  ###############################################################################################################
  f1=c()
  f2=c()
  f3=matrix(NA,n,p)
  for (k in 1:K) {
    
    va=fold[[k]]
    tr=c(1:n)[-va]
    
    f1[va]=ML_SS(z[tr,],y[tr],z[va,],tree=500)
    f2[va]=ML_SS(z[tr,],x[tr],z[va,],tree=500)
    for (j in 1:p) { f3[va,j]=ML_SS(z[tr,],m[tr,j],z[va,],tree=500) }
  }
  ###############################################################################
  Y=y-f1
  X=x-f2
  M=m-f3
  
  mar_y=c()
  for (j in 1:p) {
    mar_y[j]=lm(Y~X+M[,j]-1)$coefficients[2]
  }
  mar_m=lm(M~X-1)$coefficients
  
  sis=floor(n/log(n))
  p=p*(sis>=p)+sis*(sis<p)
  index1=order(abs(mar_y*mar_m),decreasing=T)[1:p]
  M=M[,index1]
  ######################################################################################
  bb=lm(Y~X+M-1)$coefficients[-1]
  gg=lm(M~X-1)$coefficients
  weight=1/(abs(bb*gg))
  
  cv=99999
  for (tun in seq(0.01,2,by=0.01)) {
    fit=cv.glmnet(x=cbind(X,M),y=Y,lambda=exp(seq(-10,2,length.out=100)),penalty.factor=c(0,weight^tun),intercept=F)
    if(min(fit$cvm)<cv){ beta=coef(fit,s=fit$lambda.min);cv=min(fit$cvm) }
  }
  index2=which(round(beta[-c(1,2)],3)!=0)
  M=M[,index2]
  ##############################################################################################
  be=lm(Y~X+M-1)
  be=summary(be)$coefficients
  ###############################################################################################
  alpha_est=be[1,1]
  alpha_sd=be[1,2]
  alpha_l=alpha_est-qnorm(0.975)*alpha_sd
  alpha_r=alpha_est+qnorm(0.975)*alpha_sd
  alpha_p=(1-pnorm(abs(alpha_est)/alpha_sd))*2
  ################################################################################################
  beta_est=be[-1,1]
  beta_sd=be[-1,2]
  beta_l=beta_est-qnorm(0.975)*beta_sd
  beta_r=beta_est+qnorm(0.975)*beta_sd
  Tb=abs(beta_est)/beta_sd
  beta_p=(1-pnorm(Tb))*2
  ##################################################################################################
  ga=na.omit(matrix(NA,1,4))
  GA=lm(M~X-1)
  for (j in 1:length(index2)) { ga=rbind(ga,summary(GA)[[j]]$coefficients) }
  
  gamma_est=ga[,1]
  gamma_sd=ga[,2]
  gamma_l=gamma_est-qnorm(0.975)*gamma_sd
  gamma_r=gamma_est+qnorm(0.975)*gamma_sd
  Tg=abs(gamma_est)/gamma_sd
  gamma_p=(1-pnorm(Tg))*2
  ###################################################################################################
  med_est=beta_est*gamma_est
  med_sd=sqrt(beta_est^2*gamma_sd^2+gamma_est^2*beta_sd^2)
  T_JS=pmax(Tb,Tg)
  T_med=abs(med_est)/med_sd
  
  med_Sobel=(T_JS>(sqrt(n)/log(n)))*(1-pnorm(T_med))*2+(T_JS<(sqrt(n)/log(n)))*(1-pnorm(T_med,sd=0.5))*2
  med_p=p.adjust(med_Sobel,method="bonferroni")
  med_l=(T_JS>(sqrt(n)/log(n)))*(med_est-qnorm(0.975)*med_sd)+(T_JS<(sqrt(n)/log(n)))*(med_est-qnorm(0.975,sd=0.5)*med_sd)
  med_r=(T_JS>(sqrt(n)/log(n)))*(med_est+qnorm(0.975)*med_sd)+(T_JS<(sqrt(n)/log(n)))*(med_est+qnorm(0.975,sd=0.5)*med_sd)
  ########################################################################################################
  result=data.frame(alpha_est,alpha_sd,alpha_l,alpha_r,alpha_p,
                    beta_est,beta_sd,beta_l,beta_r,beta_p,
                    gamma_est,gamma_sd,gamma_l,gamma_r,gamma_p,med_est,med_sd,med_l,med_r,med_Sobel,med_p)
  rownames(result)=index1[index2]
  return(result)
}