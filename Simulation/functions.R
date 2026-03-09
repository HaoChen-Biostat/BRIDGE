library(sparcl)
library(MASS)
library(glmnet)

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

lambda.max_j<-function(z,xj,alpha){
  n<-length(xj)
  # note that z is the matrix as the n by pz design matrix
  # and xj is a n by 1 vector as the response
  return(max(t(z)%*%xj)/(n*alpha))
}


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
    # print(i)
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


#### function for updating and biclustering 
update_xtilda<-function(s,theta_hat,xtilda_old,z){
  pz<-dim(z)[2]
  px<-dim(theta_hat)[2]
  shat<-matrix(0,nrow=pz,ncol=px)
  for(j in 1:px){
    ind<-which(s[,j]==1)
    mean1<-mean(theta_hat[ind,j])
    mean2<-mean(theta_hat[-ind,j])
    shat[which(s[,j]!=0),j]<- mean1-mean2
  }
  xtilda <- xtilda_old - z%*%shat
  return(xtilda)
}

### soft_thresholding
sthres<-function(a,b){
  if(abs(a)<=b){
    return(0) 
  }else{
    return(sign(a)*(abs(a)-b))
  }
}
### group lasso update for each group 
glasso_update<-function(b,x_tt,y_tt,lambda){
  b0<-b
  d=1
  while(d>10^-4){
    jj<-dim(x_tt)[2]
    for(j in 1:jj){
      theta_rest<-b[-j]
      if(length(theta_rest)>1){
        yc <- y_tt - x_tt[,-j]%*%theta_rest
      }else{
        yc <- y_tt - x_tt[,-j]*theta_rest
      }
      fr<-function(theta_j){1/2*sum((yc-x_tt[,j]*theta_j)^2)+lambda*sqrt(jj)*sqrt(theta_j^2)}
      fit<-optimize(fr,interval = c(-2,2),tol = 0.0001,maximum = F)
      b[j]<-fit$minimum
    }
    d<-sqrt(sum((b0-b)^2))
    # print(d)
    b0<-b
    # no need to update yc or y_tt since we update b itself
  }
  return(b) ## beta_m_update[[i]]
}


### estimated y hat
gf<-function(cbeta_env,cbeta_main,cbeta_inter,
             cenviro,cmain,cinter){
  main_ef<-cenviro%*%cbeta_env+ cmain%*%cbeta_main
  L = length(cinter)
  inter_eff<-lapply(1:L,function(i,a,b){a[[i]]%*%(b[i,]*cbeta_main)},a=cinter,b=cbeta_inter)
  return(main_ef+Reduce("+",inter_eff))
}

#### objective function 
qf<-function(cbeta_main,cbeta_inter, ## vector
             gf_value,y,
             lam1,lam2){
  mse_sq<-sum((y-gf_value)^2)/n
  
  pen<-double(length(count_group))
  for(i in 1:length(count_group)){
    ind_g<-which(group==i)
    if(length(ind_g)>1){
      pen[i]<-lam1*sqrt(length(ind_g))*sqrt(sum(cbeta_main[ind_g]^2))+lam1*sqrt(length(ind_g))*sum(apply(cbeta_inter[,ind_g],1,function(x){sqrt(sum(x^2))}))
    }else{
      pen[i]<-lam2*abs(cbeta_main[ind_g])+lam2*sum(abs(cbeta_inter[,ind_g]))
    }
  }
  return(sum(pen)+mse_sq)
}

# joint interaction model
joint_model<-function(env,cmain,cinter,y,lam1,lam2){
  #### starting values 
  beta_env_old<-rep(0,L)
  beta_main_old<-rep(0,p)
  beta_inter_old<-matrix(0,nrow=L,ncol=p)
  yfit<-rep(0,n)
  q<-c()
  loop=0
  d=1
  while(d>10^-4 && loop<201){
    loop=loop+1
    q_old<-qf(beta_main_old,beta_inter_old,yfit,y,lam1,lam2)
    ## update alpha
    ytilda<-y-yfit+env%*%beta_env_old
    beta_env_old<-solve(t(env)%*%env)%*%t(env)%*%ytilda
    yfit<-y-ytilda+env%*%beta_env_old
    ## update beta_m i.e. coefficients for group main effects 
    ## for ind_groups
    xtilda<-cmain[,ind_groups]+Reduce("+",lapply(1:L,function(i,a,b){a[[i]]%*%diag(b[i,])},a=cinter,b=beta_inter_old))[,ind_groups]
    beta_g<-beta_main_old[ind_groups]
    for(i in groups_id){
      ytilda<-y-yfit+xtilda%*%beta_g
      id_g<-which(group==i)
      noi<-ytilda-xtilda[,-id_g]%*%beta_g[-id_g]
      temp<-t(xtilda[,id_g])%*%noi
      if(sqrt(sum(temp^2))<lam1*sqrt(length(id_g))){
        # print(id_g)
        beta_g[id_g]<-rep(0,length(id_g))
      }else{
        beta_g[id_g]<-glasso_update(beta_g[id_g],xtilda[,id_g],noi,lam1)
      }
      yfit<-y-ytilda+xtilda%*%beta_g
    }
    beta_main_old[ind_groups]<-beta_g
    ## update gamma_k
    xtilda<-cmain[,ind_single]+Reduce("+",lapply(1:L,function(i,a,b){a[[i]]%*%diag(b[i,])},a=cinter,b=beta_inter_old))[,ind_single]
    beta_s<-beta_main_old[ind_single]
    for(i in 1:length(ind_single)){
      ytilda<-y-yfit+xtilda%*%beta_s
      noi<-ytilda-xtilda[,-i]%*%beta_s[-i]
      beta_s[i]<-sthres((t(xtilda[,i])%*%xtilda[,i])^(-1)*xtilda[,i]%*%noi,(t(xtilda[,i])%*%xtilda[,i])^(-1)*lam2)
      yfit<-y-ytilda+xtilda%*%beta_s
    }
    
    beta_main_old[ind_single]<-beta_s
    ## update eta_ml
    xtilda<-lapply(cinter,function(x){x[,ind_groups]%*%diag(beta_g)})
    eta<-beta_inter_old[,ind_groups]
    groups_id_nonz<-unique(group[which(beta_g!=0)])
    for(l in 1:L){
      eta_nol<-eta[-l,]
      xtilda_nol<-xtilda[-l]
      for(j in groups_id_nonz){
        ytilda<-y-yfit+Reduce("+",lapply(1:L,function(i,a,b){a[[i]]%*%b[i,]},a=xtilda,b=eta))
        id_g<-which(group==j)
        noi<-ytilda-Reduce("+",lapply(1:(L-1),function(i,a,b){a[[i]][,-id_g]%*%b[i,-id_g]},a=xtilda_nol,b=eta_nol))
        temp<-t(xtilda[[l]][,id_g])%*%noi
        if(sqrt(sum(temp^2))<lam1*sqrt(length(id_g))){
          eta[l,id_g]<-rep(0,length(id_g))
        }else{
          eta[l,id_g]<-glasso_update(eta[l,id_g],xtilda[[l]][,id_g],noi,lam1)
        }
        yfit<-y-ytilda+Reduce("+",lapply(1:L,function(i,a,b){a[[i]]%*%b[i,]},a=xtilda,b=eta))
      }
    }
    beta_inter_old[,ind_groups]<-eta
    
    ## update tau_kl
    xtilda<-lapply(cinter,function(x){x[,ind_single]%*%diag(beta_s)})
    tau<-beta_inter_old[,ind_single]
    single_id_nonz<-which(beta_s!=0)
    for(l in 1:L){
      tau_nol<-tau[-l,]
      xtilda_nol<-xtilda[-l]
      for(j in single_id_nonz){
        ytilda<-y-yfit+Reduce("+",lapply(1:L,function(i,a,b){a[[i]]%*%b[i,]},a=xtilda,b=tau))
        
        noi<-ytilda-Reduce("+",lapply(1:(L-1),function(i,a,b){a[[i]]%*%b[i,]},a=xtilda_nol,b=tau_nol))
        tau[l,j]<-sthres((t(xtilda[[l]][,j])%*%xtilda[[l]][,j])^(-1)*xtilda[[l]][,j]%*%noi,(t(xtilda[[l]][,j])%*%xtilda[[l]][,j])^(-1)*lam2)
        yfit<-y-ytilda+Reduce("+",lapply(1:L,function(i,a,b){a[[i]]%*%b[i,]},a=xtilda,b=tau))
        
      }
    }
    beta_inter_old[,ind_single]<-tau
    q_new<-qf(beta_main_old,beta_inter_old,yfit,y,lam1,lam2)
    d<-abs(q_new-q_old)/abs(q_old)
    q[loop]<-q_new
    print(d)
  }
  
  ## cut to zero if the estimation is smaller than 10^-4
  beta_main_old[which(abs(beta_main_old)<10^-4)]=0
  beta_inter_old<-t(apply(beta_inter_old,1,function(x)
  {s = which(abs(x)<10^-4)
  x[s]=rep(0,length(s))
  return(x)}))
  ## return the estimated para and seq of objective 
  return(list(beta_main=beta_main_old,beta_inter=beta_inter_old,
              beta_env=beta_env_old,yhat=yfit,obj=q))
  
}


