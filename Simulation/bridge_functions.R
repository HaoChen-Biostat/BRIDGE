############################################################
# bridge_functions.R:Functions required to run BRIDGE
############################################################

############################################################
# testbic():apply sparse biclustering to the matrix and identify one candidate module.
# Inputs: Coefficient matrix X to be clustered.
# Output: A binary bicluster indicator matrix with the same dimensions as X.
############################################################
testbic<-function(X){ ## X is scaled for the later process
  
  XX=X
  
  n=dim(X)[1]
  p=dim(X)[2]
  
  final_Cs=rep(0,n)
  final_Cv=rep(0,p)
  final_w=matrix(0,n,p)
  final_bicluster=matrix(0,n,p)
  
  result<-KMeansSparseCluster(X,2,sqrt(p))
  
  K1=which(result[[1]]$Cs==1) ## cluster 1
  K2=which(result[[1]]$Cs==2) ## cluster 2
  n1=length(K1)
  n2=length(K2)
  
  ### btw cluster sum of squares for p features
  a=matrix(0,1,p)
  for (i in 1:(n-1)) {
    for (ii in (i+1):n){
      a=a+(X[i,]-X[ii,])^2
    }
  }
  a=a*2/n
  
  w_estimation=result[[1]]$ws
  
  w_expect=0
  
  B=100
  
  for (i in 1:B){
    X_per=X[sample(n),]
    b_per=matrix(0,2,p)
    
    for (i in 1:(n1-1)){
      for (ii in (i+1):n1){
        b_per[1,]=b_per[1,]+(X_per[K1[i],]-X_per[K1[ii],])^2
      }
    }
    b_per[1,]=2*b_per[1,]/n1
    
    for (i in 1:(n2-1)){
      for (ii in (i+1):n2){
        b_per[2,]=b_per[2,]+(X_per[K2[i],]-X_per[K2[ii],])^2
      }
    }
    b_per[2,]=2*b_per[2,]/n2
    
    
    w_per=a-colSums(b_per)
    
    w_per=w_per/base::norm(w_per,'2')
    
    w_per_order=sort(w_per)
    
    w_expect=w_expect+w_per_order
  }
  
  w_expect=w_expect/B
  
  
  ks_test=ks.test(w_estimation,w_expect)
  
  
  
  if (ks_test$p.value<0.05){
    
    w_true=w_expect
    w_order=sort(w_estimation)
    w_t_order=sort(w_true)
    diff_v=w_order-w_t_order
    id=p-which.max(diff(w_order-w_t_order))
    select_v=sort(w_estimation,decreasing = T,index.return=T)$ix[1:id]
    
    
    
    sK=which.min(c(length(K1),length(K2)))
    
    id_S=which(result[[1]]$Cs==sK)
    #id_S2=setdiff(1:n,id_S)
    
    id_V=select_v#which(select_v[which.max(objective),]==1)
    print(length(id_V)) 
    
    if ((length(id_S)>1) && (length(id_V))>1){
      
      final_Cs[id_S]=1
      final_Cv[id_V]=1
      
      final_bicluster[id_S,id_V]=1
      
    } else {
      print("no cluster identified 1")
    }
  } else {
    print("no cluster identified 2")
  }
  return(final_bicluster)
  
}

############################################################
# update_bic(): remove the mean signal of an identified bicluster before searching for the next module.
############################################################
update_bic<-function(theta_prev,bcluster){
  theta_update<- theta_prev
  for(j in 1:pG){
    ind<-which(bcluster[,j]==1)
    mean_bic<-mean(theta_prev[ind,j])
    mean_nbic<-mean(theta_prev[-ind,j])
    for(i in ind){
      if(bcluster[i,j]==1){
        theta_update[i,j] = theta_prev[i,j]-mean_bic+mean_nbic
      }
    }
  }
  return(theta_update)
}

############################################################
# lambda.max_j()
# Purpose: Compute the maximum elastic-net lambda candidate for one response.
############################################################
lambda.max_j<-function(z,xj,alpha){
  n<-length(xj)
  # note that z is the matrix as the n by pz design matrix
  # and xj is a n by 1 vector as the response
  return(max(t(z)%*%xj)/(n*alpha))
}

############################################################
#cv.multi_elnet()
# Purpose: estimate the coefficient matrix.
############################################################
cv.multi_elnet<-function(x,z,alpha,nfolds){
  # type of measure is mean squared error
  # nfolds = 10 
  n<-dim(x)[[1]]
  px<-dim(x)[[2]]
  pz<-dim(z)[[2]]
  # find lambda max for all j = 1, ..., px
  lambda_max<-double(px)
  for(j in 1:px){
    lambda_max[j]<-lambda.max_j(z,x[,j],alpha)
  }
  lambda.max<-max(lambda_max)
  
  # given epsilon=0.001 and the length of sequence 100
  # construct the lambda sequence for cross-validation
  lambda.seq = exp(seq(log(0.0001*lambda.max),log(lambda.max),length=100))
  
  nf = nfolds
  foldid = sample(rep(seq(nfolds), length = n))
  lambda.seq=sort(lambda.seq,decreasing = T)
  
  l_mse<-matrix(1e+20,px,length(lambda.seq))
  
  fit_result=list()
  
  for(i in 1:px){
    xi<-x[,i]
    fit<-cv.glmnet(z,xi,foldid=foldid,alpha=alpha,lambda=lambda.seq)
    fit_result[[i]]=fit$glmnet.fit
    l_mse[i,1:length(fit$cvm)]=fit$cvm
  }
  
  lambda_idx=which.min(colSums(l_mse))
  
  lchoose=lambda.seq[lambda_idx]
  
  
  theta_esti<-matrix(0,nrow=pz,ncol=px)
  for(i in 1:px){
    xi<-x[,i]
    theta_esti[,i]<-as.numeric(fit_result[[i]]$beta[,lambda_idx])
    indx<-which(theta_esti[,i]!=0)
  }
  
  return(list(lambda.seq=lambda.seq,
              lambda.minmse=lchoose,
              theta_esti=theta_esti))
}

############################################################
# build_temporal_cov()
# Purpose: Construct group-specific temporal covariance matrices under the AR or banded setting.
############################################################
build_temporal_cov <- function(q, sigmaT = c("autreg", "band")){
  sigmaT <- match.arg(sigmaT)
  
  if(sigmaT == "autreg"){
    sigmaT1 <- 0.4^abs(outer(1:q, 1:q, "-"))
    sigmaT2 <- 0.5^abs(outer(1:q, 1:q, "-"))
  } else {
    sigmaT1 <- sigmaT2 <- 1 / (abs(outer(1:q, 1:q, "-")) + 1)
    sigmaT1[abs(row(sigmaT1) - col(sigmaT1)) > 4] <- 0
    sigmaT2[abs(row(sigmaT2) - col(sigmaT2)) > 6] <- 0
  }
  
  list(sigmaT1 = sigmaT1, sigmaT2 = sigmaT2)
}

############################################################
# build_hub_spatial_base()
# Purpose: Generate the hub spatial network, group population matrices, and differential-network.
############################################################
build_hub_spatial_base <- function(p, hub_g = 5, flip_block_ratio = 1/8){
  G1 <- huge.generator(
    n = 10,
    d = p,
    graph = "hub",
    g = hub_g,
    prob = NULL
  )
  
  Theta1 <- G1$theta
  
  omega1.total <- Theta1 * sample(c(-1, 1), p * p, replace = TRUE) * runif(p * p, 0.3, 0.5)
  omega1.total[lower.tri(omega1.total, diag = FALSE)] <- 0
  omega1.total <- omega1.total + t(omega1.total)
  diag(omega1.total) <- abs(min(eigen(omega1.total)$values)) + 0.5
  
  sigma1.total <- solve(omega1.total)
  
  omega1.total <- solve(sigma1.total)
  omega1.total[abs(omega1.total) < 10^-4] <- 0
  
  omega2.total <- omega1.total
  block_end <- floor(p * flip_block_ratio)
  omega2.total[1:block_end, 1:block_end] <- -1 * omega2.total[1:block_end, 1:block_end]
  diag(omega2.total) <- diag(omega1.total)
  sigma2.total <- solve(omega2.total)
  
  delta <- omega1.total - omega2.total
  delta <- as.matrix(delta)
  delta_vec <- as.vector(delta[upper.tri(delta, diag = FALSE)])
  
  list(
    G1 = G1,
    Theta1 = Theta1,
    omega1.total = omega1.total,
    omega2.total = omega2.total,
    sigma1.total = sigma1.total,
    sigma2.total = sigma2.total,
    delta = delta,
    delta_vec = delta_vec
  )
}

############################################################
# ：build_smallworld_spatial_base()
# Purpose: Generate the small-world spatial network, group population matrices, and differential-network.
############################################################
build_smallworld_spatial_base <- function(p,
                                          flip_block_ratio = 1/10,
                                          sw_m = 10,
                                          sw_banded_n = 6,
                                          sw_source = "SW function.R"){
  source(sw_source)
  
  G1 <- createS(
    n = 10,
    p = p,
    topology = "small-world",
    m = sw_m,
    banded.n = sw_banded_n,
    precision = TRUE
  )
  
  Theta1 <- (G1 != 0) * 1
  
  omega1.total <- Theta1 * sample(c(-1, 1), p * p, replace = TRUE) * runif(p * p, 0.3, 0.5)
  omega1.total[lower.tri(omega1.total, diag = FALSE)] <- 0
  omega1.total <- omega1.total + t(omega1.total)
  diag(omega1.total) <- abs(min(eigen(omega1.total)$values)) + 0.5
  
  sigma1.total <- solve(omega1.total)
  
  omega1.total <- solve(sigma1.total)
  omega1.total[abs(omega1.total) < 10^-4] <- 0
  
  omega2.total <- omega1.total
  block_end <- floor(p * flip_block_ratio)
  omega2.total[1:block_end, 1:block_end] <- -1 * omega2.total[1:block_end, 1:block_end]
  diag(omega2.total) <- diag(omega1.total)
  sigma2.total <- solve(omega2.total)
  
  delta <- omega1.total - omega2.total
  delta <- as.matrix(delta)
  delta_vec <- as.vector(delta[upper.tri(delta, diag = FALSE)])
  
  list(
    G1 = G1,
    Theta1 = Theta1,
    omega1.total = omega1.total,
    omega2.total = omega2.total,
    sigma1.total = sigma1.total,
    sigma2.total = sigma2.total,
    delta = delta,
    delta_vec = delta_vec
  )
}

############################################################
# build_spatial_base()
# Purpose: Dispatch to the hub or small-world spatial generator according to spatial_type.
############################################################
build_spatial_base <- function(p,
                               spatial_type = c("hub", "smallworld"),
                               hub_g = 5,
                               flip_block_ratio = NULL,
                               sw_m = 10,
                               sw_banded_n = 6,
                               sw_source = "SW function.R"){
  spatial_type <- match.arg(spatial_type)
  
  if(spatial_type == "hub"){
    if(is.null(flip_block_ratio)) flip_block_ratio <- 1/8
    out <- build_hub_spatial_base(
      p = p,
      hub_g = hub_g,
      flip_block_ratio = flip_block_ratio
    )
  } else {
    if(is.null(flip_block_ratio)) flip_block_ratio <- 1/10
    out <- build_smallworld_spatial_base(
      p = p,
      flip_block_ratio = flip_block_ratio,
      sw_m = sw_m,
      sw_banded_n = sw_banded_n,
      sw_source = sw_source
    )
  }
  
  out
}

############################################################
# generate_subject_precision()
# Purpose: Add subject-level noise to the population precision matrices.
############################################################
generate_subject_precision <- function(N, p, Theta1, omega1.total, omega2.total, noise_sd = 0.2){
  omega1_N <- vector("list", N)
  omega2_N <- vector("list", N)
  
  for(i in 1:N){
    omega1_N[[i]] <- omega1.total + Theta1 * matrix(rnorm(p * p, mean = 0, sd = noise_sd), nrow = p, ncol = p)
    omega2_N[[i]] <- omega2.total + Theta1 * matrix(rnorm(p * p, mean = 0, sd = noise_sd), nrow = p, ncol = p)
  }
  
  list(omega1_N = omega1_N, omega2_N = omega2_N)
}

############################################################
# simulate_group_X()
# Purpose: Generate subject-level fMRI matrices from spatial and temporal covariance structures.
############################################################
simulate_group_X <- function(omega_list, sigmaT, N, p){
  foreach(i = 1:N, .errorhandling = "pass", .packages = c("Matrix", "mnormt", "expm")) %dopar% {
    omega <- as.matrix(omega_list[[i]])
    omega[lower.tri(omega, diag = FALSE)] <- 0
    omega <- omega + t(omega)
    diag(omega) <- abs(min(eigen(omega)$values)) + 1
    sigma <- solve(omega)
    
    omega <- solve(sigma)
    omega[abs(omega) < 10^-4] <- 0
    
    Z <- rmnorm(n = p, mean = 0, varcov = sigmaT)
    SQS <- sqrtm(sigma)
    X_w <- SQS %*% Z
    
    rm(Z, SQS)
    gc()
    return(X_w)
  }
}

############################################################
# estimate_group_precision()
# Purpose: Estimate a sparse precision matrix from each subject-level fMRI matrix.
############################################################
estimate_group_precision <- function(X_list, N){
  foreach(i = 1:N, .errorhandling = "pass", .packages = "DensParcorr") %dopar% {
    fit <- DensParcorr(t(X_list[[i]]), dens.level = .5, select = TRUE)
    out <- fit$selected.precision
    rm(fit)
    gc()
    return(out)
  }
}

############################################################
# vectorize_fisher_features()
# Purpose: Convert precision matrices to correlations, apply the Fisher transform, and vectorize upper-triangular edges.
############################################################
vectorize_fisher_features <- function(omega_list, N, p){
  mut_in_f <- function(x){
    1/2 * log((1 + x) / (1 - x))
  }
  
  W <- lapply(lapply(omega_list, cov2cor), mut_in_f)
  X_vec <- matrix(nrow = N, ncol = p * (p - 1) / 2)
  
  for (i in 1:N){
    W_vec <- as.vector(W[[i]][upper.tri(W[[i]], diag = FALSE)])
    X_vec[i, ] <- c(W_vec)
  }
  
  list(W = W, X_vec = X_vec)
}

############################################################
# generate_theta_blocks_10()
# Purpose: Generate the coefficient matrix under the 10-module setting.
############################################################
generate_theta_blocks_10 <- function(pB, pG){
  block1 <- cbind(
    matrix(rep(c(rnorm(7, -0.6, 0.1), rep(0, pB - 7)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 10, nrow = pB)
  )
  
  block2 <- cbind(
    matrix(0, ncol = 80, nrow = pB),
    matrix(rep(c(rep(0, 80), rnorm(15, 1.1, 0.1), rep(0, pB - 95)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 90, nrow = pB)
  )
  theta2 <- block1 + block2
  
  block3 <- cbind(
    matrix(0, ncol = 85, nrow = pB),
    matrix(rep(c(rep(0, 90), rnorm(10, -0.8, 0.1), rep(0, pB - 100)), 15), ncol = 15, nrow = pB),
    matrix(0, ncol = pG - 100, nrow = pB)
  )
  theta3 <- block3 + theta2
  
  block4 <- cbind(
    matrix(0, ncol = 180, nrow = pB),
    matrix(rep(c(rep(0, 180), rnorm(18, 0.8, 0.1), rep(0, pB - 198)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 190, nrow = pB)
  )
  theta4 <- theta3 + block4
  
  block5 <- cbind(
    matrix(0, ncol = 230, nrow = pB),
    matrix(rep(c(rep(0, 230), rnorm(8, -1.0, 0.1), rep(0, pB - 238)), 14), ncol = 14, nrow = pB),
    matrix(0, ncol = pG - 244, nrow = pB)
  )
  theta5 <- theta4 + block5
  
  block6 <- cbind(
    matrix(0, ncol = 300, nrow = pB),
    matrix(rep(c(rep(0, 300), rnorm(12, 0.9, 0.1), rep(0, pB - 312)), 16), ncol = 16, nrow = pB),
    matrix(0, ncol = pG - 316, nrow = pB)
  )
  theta6 <- theta5 + block6
  
  block7 <- cbind(
    matrix(0, ncol = 360, nrow = pB),
    matrix(rep(c(rep(0, 360), rnorm(20, 0.7, 0.1), rep(0, pB - 380)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 370, nrow = pB)
  )
  theta7 <- theta6 + block7
  
  block8 <- cbind(
    matrix(0, ncol = 400, nrow = pB),
    matrix(rep(c(rep(0, 400), rnorm(10, -1.2, 0.1), rep(0, pB - 410)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 410, nrow = pB)
  )
  theta8 <- theta7 + block8
  
  block9 <- cbind(
    matrix(0, ncol = 440, nrow = pB),
    matrix(rep(c(rep(0, 440), rnorm(15, 1.3, 0.1), rep(0, pB - 455)), 15), ncol = 15, nrow = pB),
    matrix(0, ncol = pG - 455, nrow = pB)
  )
  theta9 <- theta8 + block9
  
  block10 <- cbind(
    matrix(0, ncol = 480, nrow = pB),
    matrix(rep(c(rep(0, 475), rnorm(15, -1.4, 0.1), rep(0, pB - 490)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 490, nrow = pB)
  )
  
  theta <- theta9 + block10
  return(theta)
}

############################################################
# generate_theta_blocks_15()
# Purpose: Generate the coefficient matrix under the 15-module setting.
############################################################
generate_theta_blocks_15 <- function(pB, pG){
  block1 <- cbind(
    matrix(rep(c(rnorm(7, -0.6, 0.1), rep(0, pB - 7)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 10, nrow = pB)
  )
  
  block2 <- cbind(
    matrix(0, ncol = 80, nrow = pB),
    matrix(rep(c(rep(0, 80), rnorm(12, 1.1, 0.1), rep(0, pB - 92)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 90, nrow = pB)
  )
  theta2 <- block1 + block2
  
  block3 <- cbind(
    matrix(0, ncol = 86, nrow = pB),
    matrix(rep(c(rep(0, 88), rnorm(10, -0.8, 0.1), rep(0, pB - 98)), 12), ncol = 12, nrow = pB),
    matrix(0, ncol = pG - 98, nrow = pB)
  )
  theta3 <- theta2 + block3
  
  block4 <- cbind(
    matrix(0, ncol = 110, nrow = pB),
    matrix(rep(c(rep(0, 110), rnorm(8, 0.8, 0.1), rep(0, pB - 118)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 120, nrow = pB)
  )
  theta4 <- theta3 + block4
  
  block5 <- cbind(
    matrix(0, ncol = 130, nrow = pB),
    matrix(rep(c(rep(0, 130), rnorm(12, -1.0, 0.1), rep(0, pB - 142)), 8), ncol = 8, nrow = pB),
    matrix(0, ncol = pG - 138, nrow = pB)
  )
  theta5 <- theta4 + block5
  
  block6 <- cbind(
    matrix(0, ncol = 160, nrow = pB),
    matrix(rep(c(rep(0, 160), rnorm(12, 0.9, 0.1), rep(0, pB - 172)), 6), ncol = 6, nrow = pB),
    matrix(0, ncol = pG - 166, nrow = pB)
  )
  theta6 <- theta5 + block6
  
  block7 <- cbind(
    matrix(0, ncol = 200, nrow = pB),
    matrix(rep(c(rep(0, 200), rnorm(10, 0.7, 0.1), rep(0, pB - 210)), 6), ncol = 6, nrow = pB),
    matrix(0, ncol = pG - 206, nrow = pB)
  )
  theta7 <- theta6 + block7
  
  block8 <- cbind(
    matrix(0, ncol = 230, nrow = pB),
    matrix(rep(c(rep(0, 230), rnorm(8, -1.2, 0.1), rep(0, pB - 238)), 8), ncol = 8, nrow = pB),
    matrix(0, ncol = pG - 238, nrow = pB)
  )
  theta8 <- theta7 + block8
  
  block9 <- cbind(
    matrix(0, ncol = 270, nrow = pB),
    matrix(rep(c(rep(0, 270), rnorm(12, 1.3, 0.1), rep(0, pB - 282)), 12), ncol = 12, nrow = pB),
    matrix(0, ncol = pG - 282, nrow = pB)
  )
  theta9 <- theta8 + block9
  
  block10 <- cbind(
    matrix(0, ncol = 300, nrow = pB),
    matrix(rep(c(rep(0, 300), rnorm(9, -1.4, 0.1), rep(0, pB - 309)), 9), ncol = 9, nrow = pB),
    matrix(0, ncol = pG - 309, nrow = pB)
  )
  theta10 <- theta9 + block10
  
  block11 <- cbind(
    matrix(0, ncol = 350, nrow = pB),
    matrix(rep(c(rep(0, 350), rnorm(10, 1.5, 0.1), rep(0, pB - 360)), 7), ncol = 7, nrow = pB),
    matrix(0, ncol = pG - 357, nrow = pB)
  )
  theta11 <- theta10 + block11
  
  block12 <- cbind(
    matrix(0, ncol = 380, nrow = pB),
    matrix(rep(c(rep(0, 380), rnorm(12, -1.6, 0.1), rep(0, pB - 392)), 6), ncol = 6, nrow = pB),
    matrix(0, ncol = pG - 386, nrow = pB)
  )
  theta12 <- theta11 + block12
  
  block13 <- cbind(
    matrix(0, ncol = 410, nrow = pB),
    matrix(rep(c(rep(0, 410), rnorm(10, 0.6, 0.1), rep(0, pB - 420)), 10), ncol = 10, nrow = pB),
    matrix(0, ncol = pG - 420, nrow = pB)
  )
  theta13 <- theta12 + block13
  
  block14 <- cbind(
    matrix(0, ncol = 440, nrow = pB),
    matrix(rep(c(rep(0, 440), rnorm(8, -0.7, 0.1), rep(0, pB - 448)), 8), ncol = 8, nrow = pB),
    matrix(0, ncol = pG - 448, nrow = pB)
  )
  theta14 <- theta13 + block14
  
  block15 <- cbind(
    matrix(0, ncol = 470, nrow = pB),
    matrix(rep(c(rep(0, 470), rnorm(6, 1.4, 0.1), rep(0, pB - 476)), 6), ncol = 6, nrow = pB),
    matrix(0, ncol = pG - 476, nrow = pB)
  )
  
  theta <- theta14 + block15
  return(theta)
}

############################################################
# generate_theta_blocks()
# Purpose: Dispatch to the 10-module or 15-module theta generator.
############################################################
generate_theta_blocks <- function(pB, pG, block_setting = c("block10", "block15")){
  block_setting <- match.arg(block_setting)
  
  if(block_setting == "block10"){
    theta <- generate_theta_blocks_10(pB = pB, pG = pG)
  } else {
    theta <- generate_theta_blocks_15(pB = pB, pG = pG)
  }
  
  theta
}

############################################################
# get_pca_cutoff()
# Purpose: Determine the number of module PCs needed to reach a cumulative explained-variance threshold.
############################################################
get_pca_cutoff <- function(d, threshold = 0.9){
  prop <- cumsum(sort(d / sum(d), decreasing = TRUE))
  which(prop >= threshold)[1]
}

###########################Functions for variable extraction and indexing#################################

get_delta_vec_number <- function(delta_image, order_image_number){
  idx <- match(delta_image, order_image_number)
  as.integer(na.omit(idx))
}


get_module_indices <- function(lay_list, S){
  module_number <- vector("list", S)
  module_number_gene <- vector("list", S)
  
  for (i in 1:S){
    module_number[[i]] <- which(apply(lay_list[[i]], 1, sum) != 0)
    module_number_gene[[i]] <- which(apply(lay_list[[i]], 2, sum) != 0)
  }
  
  list(module_number = module_number,
       module_number_gene = module_number_gene)
}


make_group_info <- function(s_cutoff, Z){
  valid_idx <- which(s_cutoff >= 1)
  
  if(length(valid_idx) == 0){
    group_ind <- matrix(numeric(0), nrow = 0, ncol = 2)
    group_g <- numeric(0)
    group_s <- seq_len(ncol(Z))
    group <- group_s
  } else {
    group_ind <- cbind(s_cutoff[valid_idx], 1:length(valid_idx))
    group_g <- unlist(apply(group_ind, 1, function(x){rep(x[2], x[1])}))
    group_s <- seq(from = max(group_ind[, 2]) + 1,
                   to = max(group_ind[, 2]) + ncol(Z))
    group <- c(group_g, group_s)
  }
  
  list(group_ind = group_ind,
       group_g = group_g,
       group_s = group_s,
       group = group)
}


make_component_location <- function(group_sizes){
  if(length(group_sizes) == 0) return(list())
  
  out <- vector("list", length(group_sizes))
  start <- 1
  
  for(i in 1:length(group_sizes)){
    if(group_sizes[i] <= 0){
      out[[i]] <- integer(0)
    } else {
      out[[i]] <- start:(start + group_sizes[i] - 1)
      start <- start + group_sizes[i]
    }
  }
  out
}


build_single_coef_maps <- function(X_ncol, image_number, gene_number, B_sub, G_sub){
  if(length(image_number) > 0){
    single_coef <- matrix(0, length(image_number), 2)
    single_coef[, 1] <- image_number
    single_coef[, 2] <- (X_ncol + 1):(X_ncol + ncol(B_sub))
  } else {
    single_coef <- matrix(0, 0, 2)
  }
  
  if(length(gene_number) > 0){
    single_coef_gene <- matrix(0, length(gene_number), 2)
    single_coef_gene[, 1] <- gene_number
    single_coef_gene[, 2] <- ((X_ncol + ncol(B_sub)) + 1):((X_ncol + ncol(B_sub)) + ncol(G_sub))
  } else {
    single_coef_gene <- matrix(0, 0, 2)
  }
  
  list(single_coef = single_coef,
       single_coef_gene = single_coef_gene)
}


extract_selected_locations <- function(coef_number,
                                       component_location_list,
                                       module_number,
                                       module_number_gene,
                                       single_coef,
                                       single_coef_gene,
                                       X_ncol,
                                       B_single_ncol,
                                       pB){
  selected_component <- which(sapply(component_location_list, function(idx){
    length(intersect(coef_number, idx)) > 0
  }))
  
  if(length(selected_component) > 0){
    module_image_location <- intersect(c(1:pB), unique(unlist(module_number[selected_component])))
    module_gene_location_raw <- unique(unlist(module_number_gene[selected_component]))
  } else {
    module_image_location <- integer(0)
    module_gene_location_raw <- integer(0)
  }
  
  if(nrow(single_coef) > 0){
    single_coef_location <- intersect(single_coef[, 2], coef_number) - X_ncol
    if(length(single_coef_location) > 0){
      single_image_location <- single_coef[single_coef_location, 1]
    } else {
      single_image_location <- integer(0)
    }
  } else {
    single_image_location <- integer(0)
  }
  
  if(nrow(single_coef_gene) > 0){
    single_coef_location_gene <- intersect(single_coef_gene[, 2], coef_number) - (X_ncol + B_single_ncol)
    if(length(single_coef_location_gene) > 0){
      single_gene_location_raw <- single_coef_gene[single_coef_location_gene, 1]
    } else {
      single_gene_location_raw <- integer(0)
    }
  } else {
    single_gene_location_raw <- integer(0)
  }
  
  all_image_location <- unique(c(module_image_location, single_image_location))
  all_gene_location <- unique(c(module_gene_location_raw, single_gene_location_raw)) + pB
  all_image_gene_location <- unique(c(all_image_location, all_gene_location))
  
  list(
    module_image_location = module_image_location,
    module_gene_location_raw = module_gene_location_raw,
    single_image_location = single_image_location,
    single_gene_location_raw = single_gene_location_raw,
    all_image_gene_location = all_image_gene_location
  )
}

############################################################
# tune_sparsegl_grouped()
# Purpose: Tune sparse-group logistic regression over the alpha grid using cross-validation.
############################################################
tune_sparsegl_grouped <- function(data_x_train, data_y_train,
                                  group, group_g, group_s, s_cutoff,
                                  alpha_grid = seq(0.05, 0.95, by = 0.05),
                                  nfolds = 10){
  alpha_lambda_cvm <- vector("list", length(alpha_grid))
  
  for (ii in seq_along(alpha_grid)) {
    gl_alpha <- alpha_grid[ii]
    
    fit <- tryCatch(
      cv.sparsegl(
        data_x_train,
        as.vector(data_y_train),
        group = group,
        family = "binomial",
        pred.loss = "misclass",
        nfolds = nfolds,
        asparse = gl_alpha,
        pf_sparse = c(rep(0, length(group_g)), rep(1, length(group_s))),
        pf_group = c(sqrt(s_cutoff), rep(0, length(group_s)))
      ),
      error = function(e) e
    )
    
    if(inherits(fit, "error")){
      alpha_lambda_cvm[[ii]] <- list(gl_alpha, NA, Inf, fit)
    } else {
      lambda.min_bx <- fit$lambda.min
      cvm_bx <- fit[["cvm"]][which(fit$lambda == fit$lambda.min)]
      alpha_lambda_cvm[[ii]] <- list(gl_alpha, lambda.min_bx, cvm_bx, fit)
    }
  }
  
  cvm <- sapply(alpha_lambda_cvm, function(x) x[[3]])
  gl_alpha <- sapply(alpha_lambda_cvm, function(x) x[[1]])
  gl_lambda <- sapply(alpha_lambda_cvm, function(x) x[[2]])
  fit_list <- lapply(alpha_lambda_cvm, function(x) x[[4]])
  
  min_id <- which.min(cvm)
  best_fit <- fit_list[[min_id]]
  
  list(
    alpha_lambda_cvm = alpha_lambda_cvm,
    cvm = cvm,
    gl_alpha = gl_alpha,
    gl_lambda = gl_lambda,
    min_id = min_id,
    best_fit = best_fit
  )
}

############################################################
# fit_train_standardizer()
# Purpose: Estimate column means and standard deviations from the outer training set only.
############################################################
fit_train_standardizer <- function(X_train){
  X_train <- as.matrix(X_train)
  center <- colMeans(X_train, na.rm = TRUE)
  center[!is.finite(center)] <- 0

  scale_value <- apply(X_train, 2, stats::sd, na.rm = TRUE)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1

  list(center = center, scale = scale_value)
}

############################################################
# apply_train_standardizer()
# Purpose: Apply training-derived centering and scaling to a matrix with matching columns.
############################################################
apply_train_standardizer <- function(X, prep){
  X <- as.matrix(X)
  out <- sweep(X, 2, prep$center, FUN = "-")
  out <- sweep(out, 2, prep$scale, FUN = "/")
  out[!is.finite(out)] <- 0
  as.matrix(out)
}

############################################################
# screen_image_train_test()
# Purpose: Standardize imaging features and screen them by label correlation using only outer-training information.
############################################################
screen_image_train_test <- function(B_train_raw, B_test_raw,
                                    y_train, top_k = 500){
  B_train_raw <- as.matrix(B_train_raw)
  B_test_raw <- as.matrix(B_test_raw)
  y_train <- as.numeric(y_train)

  if(nrow(B_train_raw) != length(y_train)){
    stop("length(y_train) must equal nrow(B_train_raw).")
  }
  if(ncol(B_train_raw) != ncol(B_test_raw)){
    stop("Training and test imaging matrices must have identical columns.")
  }

  prep <- fit_train_standardizer(B_train_raw)
  B_train_std_all <- apply_train_standardizer(B_train_raw, prep)
  B_test_std_all <- apply_train_standardizer(B_test_raw, prep)

  marginal_score <- vapply(seq_len(ncol(B_train_std_all)), function(j){
    z <- suppressWarnings(stats::cor(B_train_std_all[, j], y_train,
                                     method = "pearson",
                                     use = "pairwise.complete.obs"))
    if(is.finite(z)) abs(z) else 0
  }, numeric(1))

  top_k <- min(as.integer(top_k), ncol(B_train_std_all))
  selected_index <- order(marginal_score, decreasing = TRUE)[seq_len(top_k)]

  B_train <- B_train_std_all[, selected_index, drop = FALSE]
  B_test <- B_test_std_all[, selected_index, drop = FALSE]
  colnames(B_train) <- paste0("Img", seq_len(ncol(B_train)))
  colnames(B_test) <- colnames(B_train)

  list(
    B_train = B_train,
    B_test = B_test,
    selected_index = selected_index,
    marginal_score = marginal_score,
    standardizer = prep
  )
}

############################################################
# standardize_gene_train_test()
# Purpose: Estimate gene standardization on outer training data and transform both splits.
############################################################
standardize_gene_train_test <- function(G_train_raw, G_test_raw){
  prep <- fit_train_standardizer(G_train_raw)
  G_train <- apply_train_standardizer(G_train_raw, prep)
  G_test <- apply_train_standardizer(G_test_raw, prep)
  colnames(G_train) <- paste0("Gene", seq_len(ncol(G_train)))
  colnames(G_test) <- colnames(G_train)

  list(G_train = G_train, G_test = G_test, standardizer = prep)
}

############################################################
# majority_vote_with_ties()
# Purpose: Aggregate bootstrap class predictions and resolve ties using the mean probability.
############################################################
majority_vote_with_ties <- function(pred_mat, mean_prob = NULL) {
  pred_mat <- as.matrix(pred_mat)
  out <- rep(NA_real_, ncol(pred_mat))
  
  for (j in seq_len(ncol(pred_mat))) {
    cur <- pred_mat[, j]
    cur <- cur[is.finite(cur)]
    
    if (length(cur) == 0) {
      out[j] <- NA_real_
    } else if (sum(cur == 1) > sum(cur == 0)) {
      out[j] <- 1
    } else if (sum(cur == 1) < sum(cur == 0)) {
      out[j] <- 0
    } else {
      if (!is.null(mean_prob) && length(mean_prob) >= j && is.finite(mean_prob[j])) {
        out[j] <- as.numeric(mean_prob[j] >= 0.5)
      } else {
        out[j] <- NA_real_
      }
    }
  }
  out
}

############################################################
# run_module_discovery_once()
# Purpose: Estimate theta on outer training data and perform S sequential biclustering steps.
############################################################
run_module_discovery_once <- function(B_train, G_train,
                                      S = 10, alpha = 1, nfolds = 5) {
  # IMPORTANT: B_train and G_train must contain outer-training subjects only.
  G <- as.matrix(G_train)
  B <- as.matrix(B_train)
  
  pG_local <- ncol(G)
  old_pG_exists <- exists("pG", envir = .GlobalEnv)
  if (old_pG_exists) old_pG <- get("pG", envir = .GlobalEnv)
  assign("pG", pG_local, envir = .GlobalEnv)
  
  on.exit({
    if (old_pG_exists) {
      assign("pG", old_pG, envir = .GlobalEnv)
    } else if (exists("pG", envir = .GlobalEnv)) {
      rm("pG", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  theta_hat0 <- cv.multi_elnet(G, B, alpha, nfolds)
  theta_update <- theta_hat0$theta_esti
  
  theta_list <- vector("list", S)
  lay_list <- vector("list", S)
  
  for (s in seq_len(S)) {
    theta_s <- theta_update
    tmp <- suppressWarnings(suppressMessages(testbic(theta_update)))
    lay_list[[s]] <- tmp
    theta_update <- theta_list[[s]] <- update_bic(theta_s, lay_list[[s]])
  }
  
  bicluster_module <- matrix(0, nrow = ncol(B), ncol = ncol(G))
  for (i in seq_along(lay_list)) {
    bicluster_module <- bicluster_module + lay_list[[i]]
  }
  
  modules <- get_module_indices(lay_list, S)
  module_number <- modules$module_number
  module_number_gene <- modules$module_number_gene
  lay_reduce <- Reduce("+", lay_list)
  
  list(
    theta_hat0 = theta_hat0,
    theta_hat = theta_hat0$theta_esti,
    theta_update = theta_update,
    theta_list = theta_list,
    lay_list = lay_list,
    bicluster_module = bicluster_module,
    lay_reduce = lay_reduce,
    module_number = module_number,
    module_number_gene = module_number_gene
  )
}

############################################################
# build_bridge_method_data()
# Purpose: Train module-wise PCA and construct the grouped BRIDGE classification design.
############################################################
build_bridge_method_data <- function(B_train, G_train, B_test, G_test,
                                     discovery_obj, var_explained = 0.9) {
  B_train <- as.matrix(B_train)
  G_train <- as.matrix(G_train)
  B_test <- as.matrix(B_test)
  G_test <- as.matrix(G_test)

  lay_reduce <- discovery_obj$lay_reduce
  module_number <- discovery_obj$module_number
  module_number_gene <- discovery_obj$module_number_gene
  S <- length(module_number)
  
  pc_train_list <- vector("list", S)
  pc_test_list <- vector("list", S)
  pca_models <- vector("list", S)
  s_cutoff <- rep(0, S)
  
  for (s in seq_len(S)) {
    idG <- module_number_gene[[s]]
    idB <- module_number[[s]]
    
    if (length(idG) + length(idB) <= 1) {
      s_cutoff[s] <- 0
      pc_train_list[[s]] <- NULL
      pc_test_list[[s]] <- NULL
      pca_models[[s]] <- NULL
      next
    }
    
    H_train <- cbind(G_train[, idG, drop = FALSE],
                     B_train[, idB, drop = FALSE])
    H_test <- cbind(G_test[, idG, drop = FALSE],
                    B_test[, idB, drop = FALSE])

    # PCA center, covariance, loadings and retained dimension are estimated
    # exclusively from the outer training subjects.
    pca_center <- colMeans(H_train)
    H_train_centered <- sweep(H_train, 2, pca_center, FUN = "-")
    H_test_centered <- sweep(H_test, 2, pca_center, FUN = "-")

    COV <- stats::cov(H_train_centered)
    temp <- svd(COV)
    w <- temp$u
    d <- temp$d
    
    s_cutoff[s] <- get_pca_cutoff(d, threshold = var_explained)
    index_pc <- seq_len(s_cutoff[s])
    rotation <- w[, index_pc, drop = FALSE]

    pc_train_list[[s]] <- H_train_centered %*% rotation
    pc_test_list[[s]] <- H_test_centered %*% rotation
    pca_models[[s]] <- list(
      gene_index = idG,
      image_index = idB,
      center = pca_center,
      rotation = rotation,
      retained_components = s_cutoff[s],
      singular_values = d,
      variance_threshold = var_explained
    )
  }
  
  valid_pc <- !sapply(pc_train_list, is.null)
  if (!any(valid_pc)) {
    X_train <- matrix(0, nrow = nrow(B_train), ncol = 0)
    X_test <- matrix(0, nrow = nrow(B_test), ncol = 0)
  } else {
    X_train <- Reduce(cbind, pc_train_list[valid_pc])
    X_test <- Reduce(cbind, pc_test_list[valid_pc])
  }
  
  image_number <- which(apply(lay_reduce, 1, sum) == 0)
  gene_number <- which(apply(lay_reduce, 2, sum) == 0)
  
  Z_train <- cbind(
    B_train[, image_number, drop = FALSE],
    G_train[, gene_number, drop = FALSE]
  )
  Z_test <- cbind(
    B_test[, image_number, drop = FALSE],
    G_test[, gene_number, drop = FALSE]
  )
  
  group_info <- make_group_info(s_cutoff = s_cutoff, Z = Z_train)
  group_ind <- group_info$group_ind
  group_g <- group_info$group_g
  group_s <- group_info$group_s
  group <- group_info$group
  
  X_number <- if (nrow(group_ind) > 0) group_ind[, 1] else numeric(0)
  component_location_list <- make_component_location(X_number)
  
  maps <- build_single_coef_maps(
    X_ncol = ncol(X_train),
    image_number = image_number,
    gene_number = gene_number,
    B_sub = B_train[, image_number, drop = FALSE],
    G_sub = G_train[, gene_number, drop = FALSE]
  )
  
  list(
    type = "bridge",
    pB = ncol(B_train),
    pG = ncol(G_train),
    p_total = ncol(B_train) + ncol(G_train),
    main_train = cbind(X_train, Z_train),
    main_test = cbind(X_test, Z_test),
    X = X_train,
    X_train = X_train,
    X_test = X_test,
    Z_train = Z_train,
    Z_test = Z_test,
    s_cutoff = s_cutoff,
    group_ind = group_ind,
    group_g = group_g,
    group_s = group_s,
    group = group,
    image_number = image_number,
    gene_number = gene_number,
    component_location_list = component_location_list,
    single_coef = maps$single_coef,
    single_coef_gene = maps$single_coef_gene,
    module_number = module_number,
    module_number_gene = module_number_gene,
    pca_models = pca_models
  )
}

############################################################
# extract_selected_bridge_like()。
# Purpose: Map nonzero BRIDGE regression coefficients back to screened imaging and gene variables.
############################################################
extract_selected_bridge_like <- function(coef_number, method_data) {
  locs <- extract_selected_locations(
    coef_number = coef_number,
    component_location_list = method_data$component_location_list,
    module_number = method_data$module_number,
    module_number_gene = method_data$module_number_gene,
    single_coef = method_data$single_coef,
    single_coef_gene = method_data$single_coef_gene,
    X_ncol = ncol(method_data$X),
    B_single_ncol = length(method_data$image_number),
    pB = method_data$pB
  )
  sort(unique(locs$all_image_gene_location))
}

############################################################
# run_bridge_bootstrap_method()
# Purpose: Run bootstrap fitting, sparse-group tuning, test prediction, ensembling, and feature selection.
############################################################
run_bridge_bootstrap_method <- function(method_data, y_train, y_test,
                                        boot_strap = 10,
                                        true_active,
                                        p_total,
                                        tune_nfolds = 10,
                                        seed = 123) {
  set.seed(seed)
  
  # PCA/modules in method_data were learned from the outer training set.
  # main_test was produced only by applying those stored transformations.
  data_x_train_outer <- as.matrix(method_data$main_train)
  data_x_test <- as.matrix(method_data$main_test)
  groups <- method_data$group
  group_g <- method_data$group_g
  group_s <- method_data$group_s
  s_cutoff <- method_data$s_cutoff[which(method_data$s_cutoff >= 1)]

  y_train_num <- if (is.factor(y_train) || is.character(y_train)) {
    as.numeric(as.character(y_train) == "AD")
  } else as.numeric(y_train)
  y_test_num <- if (is.factor(y_test) || is.character(y_test)) {
    as.numeric(as.character(y_test) == "AD")
  } else as.numeric(y_test)
  y_test_fac <- factor(ifelse(y_test_num == 1, "AD", "HC"),
                       levels = c("HC", "AD"))

  idx_ad <- which(y_train_num == 1)
  idx_hc <- which(y_train_num == 0)
  if (length(idx_ad) == 0 || length(idx_hc) == 0) {
    stop("Both classes must occur in the outer training set.")
  }
  
  alpha_grid_fixed <- seq(0.05, 0.95, 0.05)
  
  boot_res <- foreach::foreach(
    j = 1:boot_strap,
    .errorhandling = "pass",
    .packages = c("sparsegl")
  ) %dopar% {
    set.seed(seed + j)
    
    boot_ad <- sample(idx_ad, length(idx_ad), replace = TRUE)
    boot_hc <- sample(idx_hc, length(idx_hc), replace = TRUE)
    boot_ids <- c(boot_ad, boot_hc)
    data_x_train <- data_x_train_outer[boot_ids, , drop = FALSE]
    data_y_train_num <- y_train_num[boot_ids]
    
    tmp_fit <- tryCatch(
      tune_sparsegl_grouped(
        data_x_train = data_x_train,
        data_y_train = data_y_train_num,
        group = groups,
        group_g = group_g,
        group_s = group_s,
        s_cutoff = s_cutoff,
        alpha_grid = alpha_grid_fixed,
        nfolds = tune_nfolds
      ),
      error = function(e) NULL
    )
    
    if (is.null(tmp_fit) || is.null(tmp_fit$best_fit)) {
      return(list(
        pred_class = rep(NA_real_, nrow(data_x_test)),
        pred_prob = rep(NA_real_, nrow(data_x_test)),
        selected_all = integer(0),
        fit_summary = NULL
      ))
    }
    
    fit_obj <- tmp_fit$best_fit
    
    fit_summary <- list(
      alpha_grid_result = tmp_fit$alpha_lambda_cvm,
      best_index = tmp_fit$min_id,
      alpha = tmp_fit$gl_alpha[tmp_fit$min_id],
      lambda = tmp_fit$gl_lambda[tmp_fit$min_id]
    )
    
    pred_class <- tryCatch({
      as.vector(as.numeric(predict(fit_obj, newx = data_x_test,
                                   s = fit_obj$lambda.min, type = "class")))
    }, error = function(e) rep(NA_real_, nrow(data_x_test)))
    
    pred_prob <- tryCatch({
      as.vector(as.numeric(predict(fit_obj, newx = data_x_test,
                                   s = fit_obj$lambda.min, type = "response")))
    }, error = function(e) rep(NA_real_, nrow(data_x_test)))
    
    coef_fit <- tryCatch(coef(fit_obj, s = fit_obj$lambda.min), error = function(e) NULL)
    
    selected_all <- integer(0)
    if (!is.null(coef_fit)) {
      coef_number <- (which(coef_fit != 0) - 1)[-1]
      selected_all <- extract_selected_bridge_like(
        coef_number = coef_number,
        method_data = method_data
      )
      selected_all <- selected_all[selected_all >= 1 & selected_all <= p_total]
    }
    
    list(
      pred_class = pred_class,
      pred_prob = pred_prob,
      selected_all = selected_all,
      fit_summary = fit_summary
    )
  }
  
  boot_res <- lapply(boot_res, function(x) {
    if (inherits(x, "error") || is.null(x)) {
      return(list(
        pred_class = rep(NA_real_, length(y_test_num)),
        pred_prob = rep(NA_real_, length(y_test_num)),
        selected_all = integer(0),
        fit_summary = NULL
      ))
    }
    x
  })
  
  pred_class_mat <- do.call(rbind, lapply(boot_res, function(x) x$pred_class))
  pred_prob_mat  <- do.call(rbind, lapply(boot_res, function(x) x$pred_prob))
  selected_sets <- lapply(boot_res, function(x) x$selected_all)
  fit_summaries <- lapply(boot_res, function(x) x$fit_summary)
  
  sel_mat_full <- matrix(0, nrow = boot_strap, ncol = p_total)
  for (j in seq_len(boot_strap)) {
    idx <- selected_sets[[j]]
    if (length(idx) > 0) {
      sel_mat_full[j, idx] <- 1
    }
  }
  
  mean_prob <- colMeans(pred_prob_mat, na.rm = TRUE)
  vote_class <- majority_vote_with_ties(pred_class_mat, mean_prob = mean_prob)
  
  acc <- mean(vote_class == y_test_num, na.rm = TRUE)
  auc <- safe_auc(y_test_fac, mean_prob, positive = "AD")
  
  sel_eval <- bridge_selection_from_bootstrap(
    sel_mat = sel_mat_full,
    true_active = true_active,
    p_total = p_total
  )
  
  list(
    method = "BRIDGE",
    train_acc = NA_real_,
    test_acc = acc,
    train_auc = NA_real_,
    test_auc = auc,
    pred_class_mat = pred_class_mat,
    pred_prob_mat = pred_prob_mat,
    vote_class = vote_class,
    mean_prob = mean_prob,
    acc_rule = "majority_vote_tie_by_mean_prob",
    auc_rule = "mean_probability",
    selection_eval = sel_eval,
    selected_all = sel_eval$selected_all,
    tpr = sel_eval$tpr,
    fpr = sel_eval$fpr,
    selected_sets = selected_sets,
    fit_summaries = fit_summaries
  )
}
