#--------------------------------------------------
#----------------- Simulation Study
#--------------------- exsiting menthods including: LASSO SCAD MCP
#--------------------- written in R 3.6.1
#--------------------------------------------------

library(snow)  # for parallelization 
library(snowfall) # for parallelization 
set.seed(0713)

sfInit(parallel=TRUE, cpus=10, type="SOCK")
sfLibrary(snow)
sfLibrary(snowfall)

Xtype = 1 # type of X {1,2}
corr = 0.75 # correlation {0.25, 0.75}
sp = 1 # value of s for sparsity {0, 1}
code = 'jia_example'
Iteration = 50 # number of replication 

everyloop<- function(iter){
  
  Xtype = 1
  corr = 0.75
  sp = 1
  code = 'jia_example'
  K = 10
  Sys.sleep(1)
  library(ncvreg) # for LASSO
  library(glmnet) # for SCAD and MCP
  p = 1000
  n = 50
  option=c(1,2,5,8) # true active set
  signal = c(2.1, -2.4, 2.7, -3.0) # true signal
  beta_true=rep(0, p)
  beta_true[option]=signal
  beta_true=as.matrix(beta_true) # true coefficiants
  # initialize vectors to save the results
  count_scad=rep(0,7)
  count_lasso=rep(0,7)
  count_mcp=rep(0,7)
  #--------------------------------------------------
  #------------------ read X and Y generated in net_sp_worker.jl
  #--------------------------------------------------
  file_path = paste('/simulation_data/',code,'_X_K='
                    ,K, '_n=',n,'_d=',p,'_Xtype=',Xtype,'_corr=',corr,'_sp=',sp,'_loop=',iter,'.csv',sep='')
  X_feature = data.matrix(read.table(file_path, quote="\"", comment.char=""))
  file_path = paste('/simulation_data/',code,'_Y_K='
                    ,K, '_n=',n,'_d=',p,'_Xtype=',Xtype,'_corr=',corr,'_sp=',sp,'_loop=',iter,'.csv',sep='')
  Y_net = data.matrix(read.table(file_path, quote="\"", comment.char=""))
  
  X = matrix(0, nrow = n*(n-1)/2, ncol = p+n)
  Y = rep(0, n*(n-1)/2)
  tmp = 1;
  for (i in 2:n){
    for (j in 1:(i-1)){
      vec = rep(0, n)
      vec[i] = 1
      vec[j] = 1
      X[tmp,1:p] = X_feature[i,] * X_feature[j,]
      X[tmp,(p+1):(p+n)] = vec
      Y[tmp] = Y_net[i, j]
      tmp = tmp + 1
    }
  }
  # standardize X
  X_std  =  apply(X, MARGIN=2, FUN= function(x) { (x - mean(x))/sd(x) })
  #--------------------------------------------------
  #------------------ get values for criteria:
  #------------------  
  #------------------   1: whether the true models is selected
  #------------------   2: whether the true models is contained in the selected set
  #------------------   3: true positive rate (TPR)
  #------------------   4: true positive rate (TPR) given model size <=4
  #------------------   5: false discovery rate (FDR)
  #------------------   6: model error
  #------------------   7: selected model size (LASSO only)
  #--------------------------------------------------
  
  #--------------------------------------------------
  #------------------  lasso bic
  #--------------------------------------------------
  fit_lasso = ncvreg(X_std, Y, family = "binomial", penalty="lasso", penalty.factor=c(rep(1, p), rep(0, n)), max.iter = 50000)
  # model selected by BIC
  lam = fit_lasso$lambda[which.min(BIC(fit_lasso))]
  beta_hat = coef(fit_lasso, lambda=lam)[2:(p+1)] 
  beta_hat_lasso = coef(fit_lasso, lambda=lam)[2:(p+n+1)] # save the LASSO estimator
  rank= which(beta_hat!=0)
  #------------------ 
  if ((sum(option %in% rank[1:length(option)])==length(option)) && (length(rank)==length(option))) {
    count_lasso[1] = 1
    } # whether the true models is selected
  if (sum(option %in% rank)==length(option)) {
    count_lasso[2] = 1
    } # whether the true models is contained in the selected set
  count_lasso[3] = length(intersect(rank, option)) / length(option) #TPR
  count_lasso[5] = length(setdiff(rank, option)) /length(rank) # FDR
  count_lasso[6] = t(beta_true-beta_hat)%*%(beta_true-beta_hat) # model error
  count_lasso[7] = length(rank) # save model size

  # TPR given model size <= 4
  rank = vector()
  i = 0
  while (length(rank) < 4){
   i = i + 1
   lam = fit_lasso$lambda[i]
   beta_hat = coef(fit_lasso, lambda=lam)[2:(p+1)]
   rank= which(beta_hat!=0)
  }
  count_lasso[4] = length(intersect(rank, option)) / length(option)
  #--------------------------------------------------
  #------------------  utility functions
  #--------------------------------------------------
  dePenalty <- function(z,lamb,penalty='SCAD',a_scad=3.7,a_mcp=2){
    
    if (penalty =='SCAD'){
      return(1*(z<=lamb)+pmax((a_scad*lamb-z),0)/((a_scad-1)*lamb)*(lamb<z))
    }
    if (penalty=='MCP'){
      return(pmax((1-z/(a_mcp*lamb)),0))
    }
  }

  sigmoid <- function(z){1/(1+exp(-z))}

  # loglikelihood for logistic reg
  logLikelihood_LG <- function(theta, X, y){
    h <- sigmoid(X %*% theta)
    J <- (t(y)%*%log(h)+t(1-y)%*%log(1-h))
    return(J)
  }

  LLA = function(X, y, lambda, penalty){
    
    # Step 1 using Lasso
    w = c(rep(1, p), rep(0, n))
    beta_int = coef(glmnet(X , y, family = 'binomial', alpha=1,
                           lambda = lambda, penalty.factor = w, intercept = FALSE))[-1]
    # Step 2 using local linear approximation of SCAD
    if (penalty!='LASSO'){
      w = c(rep(1, p), rep(0, n))
      for(j in 1:p){
        w[j] = dePenalty(abs(beta_int[j]), lambda , penalty= penalty, a_scad=3.7, a_mcp=2)
      }
      beta = coef(glmnet(X, y, family = 'binomial', alpha=1,
                         lambda = lambda, penalty.factor = w, intercept = FALSE))[-1]
    }else{
      beta = beta_int
    }
    return(beta)
  }

  BIC_calc = function(X, y, beta){
    n = dim(X)[1]
    num_active_feature = length(which(beta[1:p]!=0))
    k = num_active_feature + n
    BIC = k * log(n) - 2 * logLikelihood_LG(beta, X, y)
    # BIC = k*ln(n) - 2ln(L) 
    return(list(BIC = BIC, beta = beta, rank = num_active_feature))
  }

  glm_BIC <- function(X, y, lambda_list, penalty){
    
    BIC_list = c()
    rank_list = c()
    for( i in 1: length(lambda_list) ){
      lambda = lambda_list[i]
      beta = LLA(X, y, lambda, penalty)
      BIC_list[i]= BIC_calc(X, y, beta)$BIC
      rank_list[i] = BIC_calc(X, y, beta)$rank
    }
    return(list(BIC = BIC_list, rank = rank_list))
  }

  #----------------------------------------------------------------------------------
  #------------------  SCAD BIC one-step estimator by local linear approximation algorithm
  #----------------------------------------------------------------------------------
  penalty = 'SCAD'
  lambda_list = round(exp(seq(log(0.01), log(0.2), length.out = 100)), digits = 10)
  BIC_result_list = glm_BIC(X_std, Y, lambda_list, penalty = penalty)
  # Find minimum BIC score's corresponding lambda
  id = which(BIC_result_list$BIC==min(BIC_result_list$BIC, na.rm = TRUE))
  id = tail(id,1)
  lambda_BIC = lambda_list[id]
  beta_BIC = LLA(X_std, Y, lambda_BIC, penalty = penalty)
  rank = which(beta_BIC[1:p]!=0)
  beta_hat = beta_BIC[1:p]
  #------------------ 
  if ((sum(option %in% rank[1:length(option)])==length(option)) && (length(rank)==length(option))) {
    count_scad[1]=1
    } # whether the true models is selected
  if (sum(option %in% rank)==length(option)) {
    count_scad[2]=1
    } # whether the true models is contained in the selected set
  count_scad[3]=length(intersect(rank, option))  / length(option) # TPR
  count_scad[5]=length(setdiff(rank, option))/ length(rank) # FDR
  count_scad[6]=t(beta_true-beta_hat)%*%(beta_true-beta_hat) # model error
  
  # TPR given model size <= 4
  id = min(which(BIC_result_list$rank<=4), na.rm = TRUE)
  lambda_given_size = lambda_list[id]
  beta_given_size = LLA(X_std, Y, lambda_given_size, penalty = penalty)
  rank = which(beta_given_size[1:p]!=0)
  count_scad[4]=length(intersect(rank, option)) / length(option)
  #----------------------------------------------------------------------------------
  #------------------  MCP BIC one-step estimator by local linear approximation algorithm
  #----------------------------------------------------------------------------------
  penalty = 'MCP'
  lambda_list = round(exp(seq(log(0.01), log(0.2), length.out = 100)), digits = 10)
  BIC_result_list = glm_BIC(X_std, Y, lambda_list, penalty = penalty)
  # Find minimum BIC score's corresponding lambda
  id = which(BIC_result_list$BIC==min(BIC_result_list$BIC, na.rm = TRUE))
  id = tail(id,1)
  lambda_BIC = lambda_list[id]
  beta_BIC = LLA(X_std, Y, lambda_BIC, penalty = penalty)
  rank = which(beta_BIC[1:p]!=0)
  beta_hat = beta_BIC[1:p]
  #------------------ 
  if ((sum(option %in% rank[1:length(option)])==length(option)) && (length(rank)==length(option))) {
    count_mcp[1]=1
    }  # whether the true models is selected
  if (sum(option %in% rank)==length(option)) {
    count_mcp[2]=1
    } # whether the true models is contained in the selected set
  count_mcp[3]=length(intersect(rank, option)) / length(option) # TPR
  count_mcp[5]=length(setdiff(rank, option)) /length(rank) # FDR
  count_mcp[6]=t(beta_true-beta_hat)%*%(beta_true-beta_hat) # model error
  # TPR given model size <= 4
  id = min(which(BIC_result_list$rank<=4), na.rm = TRUE)
  lambda_given_size = lambda_list[id]
  beta_given_size = LLA(X_std, Y, lambda_given_size, penalty = penalty)
  rank = which(beta_given_size[1:p]!=0)
  count_mcp[4] = length(intersect(rank, option))/ length(option)
  #--------------------------------------------------
  #------------------  combine results
  #--------------------------------------------------
  result=rbind(count_lasso, count_scad, count_mcp)
  return(result)
}
#--------------------------------------------------
#------------------ write results 
#--------------------------------------------------
findata<-sfClusterApplyLB(1:Iteration, everyloop)
write.csv(findata,file = paste('/simulation_result_raw/', code, '_LLA_om_raw_Xtype='
                               ,Xtype,'_corr=',corr,'_sp=',sp,'.csv',sep=''))
data=read.csv(paste('/simulation_result_raw/', code, '_LLA_om_raw_Xtype='
                    ,Xtype,'_corr=',corr,'_sp=',sp,'.csv',sep=''))
data=as.matrix(data)
data=data[,2:(Iteration*7+1)]
data=as.numeric(data)
data=matrix(data, ncol=Iteration*7, nrow=3)
K_list = data[1,7]
count = 1
for (i in 1:(Iteration-1)){
  if (sum(is.na(data[,(i*7+1):(i*7+7)])) == 0){
    count = count + 1
    data[,1:7]=data[,1:7]+data[,(i*7+1):(i*7+7)]
}
  K_list = c(K_list, data[1,(i*7+7)])
}
result=data[,1:7]/count
rownames(result) <-c("LASSO", "SCAD", "MCP")
# write.csv(K_list,file = paste('/simulation_result/',code,'_lasso_K_Xtype=',Xtype,'_corr=',corr,'_sp=',sp,'.csv',sep=''))
write.csv(result,file = paste('/simulation_result/',code,'_LLA_om_Xtype=',Xtype,'_corr=',corr,'_sp=',sp,'.csv',sep=''))
