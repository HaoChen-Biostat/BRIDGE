rm(list=ls())


library(sparcl)
library(mvtnorm)
library(MASS)
library(huge)
library(Matrix, warn.conflicts = FALSE)
library(expm, warn.conflicts = FALSE)
library(mnormt)
library(clime)
library(glmnet)
library(foreach)
library(doParallel)
library(e1071)
library(iterators)
library(parallel)
library(DensParcorr)
library(sparsegl)
library(lpSolve)
library(gplots, warn.conflicts = FALSE)
library(igraph, warn.conflicts = FALSE)

#set.seed()


#Load the functions required for multivariate sparse regression (theta estimation) and module identification: functions.R
# and small-world network generation:SW function.R

source("functions.R")
source("SW function.R")




# -------------------------------------------------------------
# Simulation parameter settings
# -------------------------------------------------------------

l<-500
p<-100
q=50
N<-150
times<-10
boot_strap <- 100
N_tr<-0.7*N
N_ts<-0.3*N


# Initialize containers to store simulation results

var.iden <- list()
error.rate <- array()
Y_pre.vote<-array()
Y_pre <- matrix(nrow=boot_strap,ncol=2*N_ts)
result_list<-list()

B_image_list<-list()
B_image_scale_list<-list()

# Generate temporal covariance matrices


sigmaT = "band" #autreg

if(sigmaT=="autreg"){
  sigmaT1 <- 0.4^abs(outer(1:q,1:q,"-"))
  sigmaT2 <- 0.5^abs(outer(1:q,1:q,"-"))
} else if(sigmaT=="band"){
  sigmaT1 <- sigmaT2 <- 1/(abs(outer(1:q,1:q,"-"))+1)
  sigmaT1[abs(row(sigmaT1)-col(sigmaT1))>4] <- 0
  sigmaT2[abs(row(sigmaT2)-col(sigmaT2))>6] <- 0
} 


# Generate spatial covariance matrices (small-world network or hub network)



G1 <- createS(n=10,p=p,topology="small-world",m=10,banded.n=6,precision=T)
Theta1 <- (G1!=0)*1


# Precision matrix 

omega1.total = Theta1 * sample(c(-1,1), p * p, replace = TRUE) * runif(p * p, 0.3,0.5)
omega1.total[lower.tri(omega1.total, diag = FALSE)] <- 0
omega1.total <- omega1.total + t(omega1.total)
diag(omega1.total) = abs(min(eigen(omega1.total)$values)) + 0.5
sigma1.total = solve(omega1.total)

omega1.total = solve(sigma1.total)
omega1.total[abs(omega1.total) < 10 ^ -4] <- 0

omega2.total <- omega1.total
omega2.total[1:(p/10*1),1:(p/10*1)] <- -1*omega2.total[1:(p/10*1),1:(p/10*1)]

diag(omega2.total) <- diag(omega1.total)
sigma2.total <- solve(omega2.total)

# True difference network

delta <- omega1.total-omega2.total
delta <- as.matrix(delta)
delta_vec <- as.vector(delta[upper.tri(delta, diag = FALSE)])

omega1_N <- list()
omega2_N <- list()

# Add noise to generate individual-level precision matrices
for(i in 1:N){ 
  omega1_N[[i]] <- omega1.total + Theta1 * matrix(rnorm(p*p, mean = 0, sd = 0.2),nrow=p,ncol=p)
  omega2_N[[i]] <- omega2.total + Theta1 * matrix(rnorm(p*p, mean = 0, sd = 0.2),nrow=p,ncol=p)
}  


# Parallel computing setup

cores<-10
cl <- makeCluster(cores)
registerDoParallel(cl)


#  Main loop: times repetitions

for (r in 1:times) {
  

  # Generate brain imaging data

  X1_w <- list()
  X2_w <- list()
  X1_omega <- list()
  X2_omega <- list()
  W1 <- list()
  W2 <- list()
  
  X1_w <- foreach(i = 1:N, .errorhandling = 'pass', .packages = c("Matrix", "mnormt","expm")) %dopar% {
    #control group
    omega1 <- as.matrix(omega1_N[[i]])
    omega1[lower.tri(omega1, diag = FALSE)] <- 0
    omega1 <- omega1 + t(omega1)
    diag(omega1) = abs(min(eigen(omega1)$values)) + 1
    sigma1 = solve(omega1)
    # 
    omega1 = solve(sigma1)
    omega1[abs(omega1) < 10 ^ -4] <- 0
    # 
    Z = rmnorm(n=p, mean=0, sigmaT1)
    SQS = sqrtm(sigma1)  
    X1_w = SQS%*%Z 
    rm(Z)
    gc()
    rm(SQS)
    gc()
    return(X1_w)
  }
  X1_omega <- foreach(i = 1:N, .errorhandling = 'pass', .packages = "DensParcorr") %dopar% {
    X1.clime <- DensParcorr(t(X1_w[[i]]),dens.level =.5,select=TRUE)
    X1_omega <- X1.clime$selected.precision
    rm(X1.clime)
    gc()
    return(X1_omega)
  }
  print(r)
  
  X2_w <- foreach(i = 1:N, .errorhandling = 'pass', .packages = c("Matrix", "mnormt","expm")) %dopar% {
    #case group
    omega2 <- omega2_N[[i]]
    omega2[lower.tri(omega2, diag = FALSE)] <- 0
    omega2 <- omega2 + t(omega2)
    diag(omega2) = abs(min(eigen(omega2)$values)) + 1
    sigma2 = solve(omega2)
    
    omega2 = solve(sigma2)
    omega2[abs(omega2) < 10 ^ -4] <- 0
    
    Z = rmnorm(n=p, mean=0, sigmaT2)
    SQS = sqrtm(sigma2)
    X2_w = SQS%*%Z
    rm(Z)
    gc()
    rm(SQS)
    gc()
    return(X2_w)
  }
  X2_omega <- foreach(i = 1:N, .errorhandling = 'pass', .packages = "DensParcorr") %dopar% {
    X2.clime <- DensParcorr(t(X2_w[[i]]),dens.level =.5,select=TRUE)
    X2_omega <- X2.clime$selected.precision
    rm(X2.clime)
    gc()
    return(X2_omega)
  }
  print(r)
  
  
  # Fisher transformation and vectorization
  X1_vec<-matrix(nrow=N,ncol=p*(p-1)/2)
  W1_vec <- array()
  X2_vec<-matrix(nrow=N,ncol=p*(p-1)/2)
  W2_vec <- array()
  

  mut_in_f<-function(x){
    1/2*log((1+x)/(1-x))
  }
  
  W1<-lapply(lapply(X1_omega,cov2cor),mut_in_f)
  W2<-lapply(lapply(X2_omega,cov2cor),mut_in_f)
  
  for (i in 1:N) {
    W1_vec <- as.vector(W1[[i]][upper.tri(W1[[i]], diag = FALSE)])
    X1_vec[i,] <- c(W1_vec)
    
    W2_vec <- as.vector(W2[[i]][upper.tri(W2[[i]], diag = FALSE)])
    X2_vec[i,] <- c(W2_vec)
  }
  
  ## Combine and standardize data
  B_image<-rbind(X1_vec,X2_vec)
  B_image_scale_pri<-scale(B_image,center = TRUE,scale = TRUE)
  B_image_list[[r]]<-B_image
  B_image_scale_list[[r]]<-B_image_scale_pri
  
  
  # Feature selection: Top 500 by Pearson correlation
  pBB<-p*(p-1)/2
  data_y<-array(c(rep(1,N),rep(0,N)),dim = c(2*N,1))
  pearson_list<-vector()
  for (i in 1:pBB) {
    pearson_list[i]<-cor (B_image_scale_pri[,i], data_y, method="pearson")
  }
  pearson_list<-abs(pearson_list)
  order_image<-order(pearson_list,decreasing = TRUE)[1:500]
  order_image_number<-order_image
  delta_image<-which(delta_vec!=0)
  B_image_scale<-B_image_scale_pri[,order_image_number]

  
  #Construct theta matrix 
  pB<-dim(B_image_scale)[2]
  pG<-l
  n<-2*N

  ## Build 10 regulatory modules (simulating true gene-imaging associations)
  #block1
  block1<-cbind(matrix(rep(c(rnorm(7, -0.6,0.1), rep(0, pB-7)),10 ), ncol=10, nrow=pB), matrix(0, ncol=pG-10, nrow=pB))
  
 
  #block2
  block2<-cbind(matrix(0, ncol=80, nrow=pB),
                matrix(rep( c(rep(0,80)  ,rnorm(15, 1.1,0.1), rep(0, pB-95)),10 ),  ncol=10, nrow=pB),
                matrix(0, ncol=pG-90, nrow=pB))
  theta2<-block1+block
  #block3
  block3<-cbind(matrix(0, ncol=85, nrow=pB),
                matrix(rep( c(rep(0,90)  ,rnorm(10, -0.8,0.1), rep(0, pB-100)),15 ),  ncol=15, nrow=pB),
                matrix(0, ncol=pG-100, nrow=pB))
  theta3<-block3+theta2

  #block4
  block4<-cbind(matrix(0, ncol=180, nrow=pB),
                matrix(rep( c(rep(0,180)  ,rnorm(18,0.8,0.1), rep(0, pB-198)),10 ),  ncol=10, nrow=pB),
                matrix(0, ncol=pG-190, nrow=pB))
  theta4<-theta3+block4

  #block5
  block5<-cbind(matrix(0, ncol=230, nrow=pB),
                matrix(rep( c(rep(0,230)  ,rnorm(8,-1.0,0.1), rep(0, pB-238)),14 ),  ncol=14, nrow=pB),
                matrix(0, ncol=pG-244, nrow=pB))
  theta5<-theta4+block5

  #block6
  block6<-cbind(matrix(0, ncol=300, nrow=pB),
                matrix(rep( c(rep(0,300)  ,rnorm(12,0.9,0.1), rep(0, pB-312)),16 ),  ncol=16, nrow=pB),
                matrix(0, ncol=pG-316, nrow=pB))
  theta6<-theta5+block6

  #block7
  block7<-cbind(matrix(0, ncol=360, nrow=pB),
                matrix(rep( c(rep(0,360)  ,rnorm(20,0.7,0.1), rep(0, pB-380)),10 ),  ncol=10, nrow=pB),
                matrix(0, ncol=pG-370, nrow=pB))
  theta7<-theta6+block7

  #block8
  block8<-cbind(matrix(0, ncol=400, nrow=pB),
                matrix(rep( c(rep(0,400)  ,rnorm(10,-1.2,0.1), rep(0, pB-410)),10 ),  ncol=10, nrow=pB),
                matrix(0, ncol=pG-410, nrow=pB))
  theta8<-theta7+block8

  #block9
  block9<-cbind(matrix(0, ncol=440, nrow=pB),
                matrix(rep( c(rep(0,440)  ,rnorm(15,1.3,0.1), rep(0, pB-455)),15 ),  ncol=15, nrow=pB),
                matrix(0, ncol=pG-455, nrow=pB))
  theta9<-theta8+block9

  #block10
  block10<-cbind(matrix(0, ncol=480, nrow=pB),
                 matrix(rep( c(rep(0,475)  ,rnorm(15,-1.4,0.1), rep(0, pB-490)),10 ),  ncol=10, nrow=pB),
                 matrix(0, ncol=pG-490, nrow=pB))
  theta<-theta9+block10
  image(1:nrow(theta), 1:ncol(theta),theta)

  # generate gene expression data 
  epsilon<-matrix(rnorm(pG*n), ncol=pG)
  G_scale<-matrix(0,n,l)
  G_scale<-B_image_scale%*%theta + epsilon
  
  
  

  
  

  G<-G_scale
  B<-B_image_scale
  n<-dim(G)[[1]]
  pG<-dim(G)[[2]]
  pB<-dim(B)[[2]]
  a <- matrix(nrow=boot_strap,ncol=(pB+pG))
  
  #multivariate sparse regression to estimate theta
  alpha = 1
  nfolds= 5 
  theta_hat0<-cv.multi_elnet(G,B,alpha,nfolds)
  theta_update<-theta_hat0$theta_esti
  image(1:nrow(theta_update), 1:ncol(theta_update),theta_update)
  
  #Iterative biclustering to identify modules
  theta_list<-lay_list<-list()
  S=10
  for(s in 1:S){ 
    theta_s<-theta_update
    lay_list[[s]]<-testbic(theta_update) ## sparse clustering 
    theta_update<-theta_list[[s]]<-update_bic(theta_s,lay_list[[s]]) ## subtract identified module
  }
  
  bicluster_module<-matrix(0,pB,pG)
  for (i in 1:length(lay_list)) {
    bicluster_module<-bicluster_module+lay_list[[i]]
  }
  image(1:nrow(bicluster_module), 1:ncol(bicluster_module),bicluster_module)
  
  ###PCA dimension reduction within module
  COV_list<-s_d_list<-s_w_list<-list()
  idx_list<-idz_list<-list()
  pc_list<-list()
  s_cutoff<-c()
  c=0.9
  for(s in 1:S){
    slay<-lay_list[[s]]
    idG<-which(apply(slay,2,sum)!=0)
    idB<-which(apply(slay,1,sum)!=0)
    COV<-cov(cbind(G[,idG],B[,idB]))
    temp<-svd(COV)
    w<-temp$u
    d<-temp$d
    ind_d<-order(d,decreasing = T)
    s_cutoff[s]<-which.min(cumsum(sort(d/sum(d),decreasing = T))<c)
    index_pc<-ind_d[1:s_cutoff[s]] 
    pc_list[[s]]<-cbind(cbind(G[,idG],B[,idB]))%*%w[,index_pc] 
  }
  
  # Combine all module PCs
  X<- Reduce(cbind,pc_list)
  # Independent variables (not in any module)
  lay_reduce<-Reduce("+",lay_list)
  Z<-cbind(B[,which(apply(lay_reduce,1,sum)==0)],G[,which(apply(lay_reduce,2,sum)==0)])
  
  
  
  #------------------------- Module information---------------------------------------------
  main<-cbind(X,Z)
  group_ind<-cbind(s_cutoff[which(s_cutoff>=1)],1:length(s_cutoff[which(s_cutoff>=1)]))
  group_g<-unlist(apply(group_ind,1,function(x){rep(x[2],x[1])}))
  group_s<-seq(from=max(group_ind[,2])+1,to=max(group_ind[,2])+dim(Z)[2])
  group<-c(group_g,group_s)

  image_number<-which(apply(lay_reduce,1,sum)==0)
  gene_number<-which(apply(lay_reduce,2,sum)==0)
  delta_vec_number<-match(delta_image,order_image_number)
  delta_gene<-unique(which(theta[delta_vec_number,]!=0,arr.ind = TRUE)[,2])

  module_number<-matrix(0,40,S)
  for (i in 1:S) {
    module_number[1:length(which(apply(lay_list[[i]],1,sum)!=0)),i]<-which(apply(lay_list[[i]],1,sum)!=0)
  }
  module_number_gene<-matrix(0,40,S)
  for (i in 1:S) {
    module_number_gene[1:length(which(apply(lay_list[[i]],2,sum)!=0)),i]<-which(apply(lay_list[[i]],2,sum)!=0)
  }

  delta_module_intersect<-list()
  delta_module_intersect_loacation<-array()
  for (i in 1:S) {
    delta_module_intersect[i]<-length(intersect(delta_vec_number,module_number[,i]))
    
  }
  delta_module_intersect_gene<-list()
  delta_module_intersect_loacation_gene<-array()
  for (i in 1:S) {
    delta_module_intersect_gene[i]<-length(intersect(delta_gene,module_number_gene[,i]))
    
  }

  delta_module_intersect_loacation<-which(delta_module_intersect!=0)
  delta_module_intersect_loacation_gene<-which(delta_module_intersect_gene!=0)

  single_delta_intersect<-intersect(delta_vec_number,image_number)
  single_delta_intersect_gene<-intersect(delta_gene,gene_number)

  X_number<-group_ind[,1]

  PCA_number_location<-matrix(0,max(X_number),S)
  PCA_number_location[1:X_number[1],1]<-1:sum(X_number[1:1])
  for (i in 2:S) {
    PCA_number_location[1:X_number[i],i]<-(sum(X_number[1:i-1])+1):(sum(X_number[1:i-1])+X_number[i])
  }

  single_coef<-matrix(0,length(image_number),2)
  single_coef[,1]<-image_number
  single_coef[,2]<-(dim(X)[2]+1):(dim(X)[2]+dim(B[,image_number])[2])
  
  single_coef_gene<-matrix(0,length(gene_number),2)
  single_coef_gene[,1]<-gene_number
  single_coef_gene[,2]<-((dim(X)[2]+dim(B[,image_number])[2])+1):((dim(X)[2]+dim(B[,image_number])[2])+dim(G[,gene_number])[2])#后

  delta_gene<-unique(which(theta[delta_vec_number,]!=0,arr.ind = TRUE)[,2])+pB
  delta_gene_image<-c(delta_vec_number,delta_gene)

  delta_gene_image_0<-setdiff(c(1:(pB+pG)),delta_gene_image)

  #Ensemble Penalized Logistic Regression
  data_x<-main
  groups<-group
  data_y<-array(c(rep(1,N),rep(0,N)),dim = c(2*N,1))
  data_x_AD<-data_x[1:N,]
  data_x_CN<-data_x[(N+1):(2*N),]
  data_y_AD<-data_y[c(1:N)]
  data_y_CN<-data_y[c((N+1):(2*N))]

  tr <- sample(1:N,N_tr,replace = F)
  #test data
  data_x_test <- rbind(data_x_AD[-tr,],data_x_CN[-tr,])
  data_y_test <- array(c(rep(1,nrow(data_x_AD[-tr,])),rep(0,nrow(data_x_CN[-tr,]))),dim = c(nrow(data_x_test),1))

  alpha_lambda_cvm<-list()
  cvm<-array()
  gl_alpha<-array()
  gl_lambda<-array()
  grouplasso<-list()
  grouplasso.cvfit_boot_strap<-foreach(j = 1:boot_strap, .errorhandling = 'pass', .packages = "sparsegl") %dopar% {
    ids <- sample(tr, N_tr, replace=T)
    #training data
    data_x_train<- rbind(data_x_AD[ids,],data_x_CN[ids,])
    data_y_train <- array(c(rep(1,N_tr),rep(0,N_tr)),dim = c(N_tr*2,1))
    
    for (ii in 1:19) {
      gl_alpha<-ii*0.05
      grouplasso.cvfit<-cv.sparsegl(data_x_train, data_y_train, group = groups,family = 'binomial',pred.loss = "misclass",nfolds = 10,asparse =gl_alpha,pf_sparse=c(rep(0,length(group_g)),rep(1,length(group_s))),pf_group=c(sqrt(s_cutoff),rep(0,length(group_s))))
      lambda.min_bx<-grouplasso.cvfit$lambda.min
      cvm_bx<-grouplasso.cvfit[["cvm"]][which(grouplasso.cvfit$lambda==grouplasso.cvfit$lambda.min)]
      alpha_lambda_cvm[[ii]]<-list(gl_alpha,lambda.min_bx,cvm_bx,grouplasso.cvfit)
    }
    
    for (i in 1:19) {
      cvm[i]<-alpha_lambda_cvm[[i]][[3]]
      gl_alpha[i]<-alpha_lambda_cvm[[i]][[1]]
      gl_lambda[i]<-alpha_lambda_cvm[[i]][[2]]
      grouplasso[[i]]<-alpha_lambda_cvm[[i]][[4]]
    }
    min<-which.min(cvm)
    cvm_min<-cvm[min]
    gl_alpha_min<-gl_alpha[min]
    gl_lambda_min<-gl_lambda[min]
    grouplasso.cvfit_boot_strap<-grouplasso[[min]]
    
    
    return(grouplasso.cvfit_boot_strap)
  }
  print(r)
  
  
  # # Collect results from each bootstrap iteration
  for (j in 1:boot_strap) {
    grouplasso.cvfit<-grouplasso.cvfit_boot_strap[[j]]
    Y_pre[j,] <- as.vector(as.numeric(predict(grouplasso.cvfit,  newx = data_x_test, s = grouplasso.cvfit$lambda.min,type='class')))
    coef <- coef(grouplasso.cvfit, s =grouplasso.cvfit$lambda.min )
    coef_number<-(which(coef!=0)-1)[-1]
    PCA_coef_intersect<-list()
    PCA_coef_intersect_loacation<-array()
    for (i in 1:S) {
      PCA_coef_intersect[i]<-length(intersect(coef_number,PCA_number_location[,i]))
    }
    PCA_coef_intersect_loacation<-which(PCA_coef_intersect!=0)
    module_image_location<-intersect(c(1:pB),module_number[,PCA_coef_intersect_loacation])
    module_gene_location<-intersect(c(1:pG),module_number_gene[,PCA_coef_intersect_loacation])
    
    single_coef_location<-intersect(single_coef[,2],coef_number)-dim(X)[2]
    single_coef_location_gene<-intersect(single_coef_gene[,2],coef_number)-(dim(X)[2]+dim(B[,image_number])[2])
    single_image_location<-single_coef[c(single_coef_location),1]
    single_gene_location<-single_coef_gene[c(single_coef_location_gene),1]
    
    all_image_location<-c(module_image_location,single_image_location)
    all_gene_location<-c(module_gene_location,single_gene_location)+pB
    all_image_gene_location<-c(all_image_location,all_gene_location)
    
    all_image_gene_location_0<-setdiff(c(1:(pB+pG)),all_image_gene_location)
    #------------------------------------------------------------------------------------------
    for (a_t in 1:(pB+pG)) {
      if (a_t %in% all_image_gene_location ){
        a[j,a_t]<-1                        
      }else{
        a[j,a_t]<-0
      }
    }
    
    
  }
  
  
  
  var.iden[[r]] <- a
  #classification error vote
  for (i in 1:(2*N_ts)) {
    if(sum(Y_pre[,i])>(boot_strap/2)){
      Y_pre.vote[i]<-1
    } else if(sum(Y_pre[,i])<(boot_strap/2)){
      Y_pre.vote[i]<-0
    } else if(sum(Y_pre[,i])==(boot_strap/2)){
      Y_pre.vote[i]<-sample(c(0,1),1)
    }
  }
  
  error.rate[r] <- sum((data_y_test - Y_pre.vote)^2) / length(data_y_test)
  
  #save
  result_list[[r]]<-list(B_image_scale,G_scale,theta,theta_hat0$theta_esti,bicluster_module,lay_list,group_g,group_s,s_cutoff,main,group,tr,delta_image,delta_vec,delta_vec_number,delta_gene)
  names(result_list[[r]])<-c("B_image_scale","G_scale","theta","theta_hat0$theta_esti","bicluster_module","lay_list","group_g","group_s","s_cutoff","main","group","tr","delta_image","delta_vec","delta_vec_number","delta_gene")


  
}




parallel::stopCluster(cl)
end_time <- Sys.time()
cat("计算时间：", end_time - start_time, "\n")



# Classification error

error.rate_final <- mean(error.rate)
sd_final <- sd(error.rate)  



#---------------------tpr,tdr,roc----------------------------
a_sum1<-array(dim=(pB+pG))
a_1<-array(dim=(pB+pG))
precision<-matrix(nrow=times,ncol=boot_strap)
tpr<-matrix(nrow=times,ncol=boot_strap)
tdr<-matrix(nrow=times,ncol=boot_strap)
tnr<-matrix(nrow=times,ncol=boot_strap)
spec<-matrix(nrow=times,ncol=boot_strap)


for (r in 1:times) {
  a_sum <- t(as.data.frame(colSums(var.iden[[r]])))
  for (l in 1:(boot_strap)) {
    for (i in 1:((pB+pG))) {
      if (a_sum[i]>=l) {
        a_sum1[i]=1
      } else {
        a_sum1[i]=0
      }
    }
    for (b in 1:((pB+pG))) {
      if (a_sum1[b]==0 & (b %in% delta_gene_image_0))
      {
        a_1[b]<-1
      } 
      else if(a_sum1[b]!=0 & (b %in% delta_gene_image)){
        a_1[b]<-2
      }
      else if(a_sum1[b]==0 & (b %in% delta_gene_image)){
        a_1[b]<-3
      }
      else if(a_sum1[b]!=0 & (b %in% delta_gene_image_0)){
        a_1[b]<-4
      }
      else{
        a_1[b]<-0
      }
    }
    precision[r,l] <- sum(a_1==2) / sum(a_sum1!=0)
    tpr[r,l] <- sum(a_1==2) / length(delta_gene_image)
    
    tnr[r,l] <- sum(a_1==1) / length(delta_gene_image_0)
    spec[r,l] <- sum(a_1==1)/ length(delta_gene_image_0)
    tdr[r,l] <- sum(a_1==2)/(sum(a_1==4)+sum(a_1==2)) 
  }
}
P<-colMeans(precision)
R <- colMeans(tpr)
F<-1-colMeans(spec)

specificity<-colMeans(spec)



plot(F, R, col=3,lwd=2, type="l",xlab="1-Spec",ylab="Recall", main="ROC Curve")




#------------------save_result----------------------------------
save.image("output_20.RData")



