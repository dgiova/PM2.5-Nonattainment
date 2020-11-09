
## code for makephi function, from HEIfunctions package
makephi <- function(coords, scale= 10){
  -log(0.05)/(max(dist(coords))/scale)
}



mvpsmod <- function(formula, data, trt, coords, nsamp, thin,
                    tuning = list(A = 0.1, psi = 0.2, theta = 1), 
                    prior = list(KIG = rep(.01,2), psi = rep(.01,2), theta1 = rep(.6,2), theta2 = rep(10,2)),
                    starting = list(B=NULL, A = NULL, psi = NULL, theta = NULL), outputbinsize = NULL, outputfilename = NULL){
  
  require(spBayes)
  require(corpcor)
  require(MASS)
  
  ncov <- dim(model.matrix(formula, data))[2]
  
  print("Note: function does not accept missing data in covariates.")
  logit	= function(theta, a, b){
    return(log((theta-a)/(b-theta)))
  }
  
  
  logitInv	= function(z, a, b){
    return ( b-(b-a)/(1+exp(z)) );
  }
  
  ###########################
  ##### Format the Data #####
  ###########################
  
  n 	<- dim(data)[1]
  q	<- 2
  
  names  <- colnames(model.matrix(formula, data))
  data$a <- data[,trt]
  data$y <- data[,all.vars(formula)[1]]
  data$ytemp 	<- 1
  nterm 	<- length(labels(terms(formula)))
  terms 	<- NULL
  for(ii in 1:nterm){
    if(ii == 1) terms <- labels(terms(formula))[ii]
    if(ii > 1) terms 	<- paste(terms, "+", labels(terms(formula))[ii], sep=" ")
  }
  formulatemp <- as.formula(paste("ytemp ~ ", terms))
  covs 	 	<- model.matrix(formulatemp, data)
  
  makeXmat	= function(Xvals){
    #Xvals	= cbind(rep(1,dim(Xvals)[[1]]), Xvals)
    p	= dim(Xvals)[[2]]
    Xout	= matrix(0, n*q, p*q)
    Xout[seq(1,n*q,q),1:p]		= Xvals
    Xout[seq(2,n*q,q),(p+1):(2*p)]= Xvals
    return(Xout)
  }
  
  X	= makeXmat(covs)
  p	= ncol(X)
  
  y0			= rep(NA,n)
  y1			= rep(NA,n)
  y0[data$a==0]	= data$y[data$a==0]
  y1[data$a==1]	= data$y[data$a==1]
  
  Y			= rep(-999,n*q)
  Y[seq(1,n*q,q)]	= y0
  Y[seq(2,n*q,q)]	= y1
  ismissy		= (is.na(Y))
  nwithmissing	= sum(ismissy)
  
  a			= rep(NA, n*q)
  a[seq(1,n*q,q)]	= data$a
  a[seq(2,n*q,q)]	= data$a
  
  
  
  ###########################
  ##### Select Priors #######
  ###########################
  
  ####Priors for AA' (Inverse Wishart)
  #####With q=2, this is just two inverse gamma priors
  KIG_a	= prior$KIG
  KIG_b	= prior$KIG
  
  ####Priors for Psi (Inverse gamma)
  psiig_a	= prior$psi
  psiig_b	= prior$psi
  
  ####Priors for theta (uniform)
  thetaunif_a	= prior$theta1 #miles
  thetaunif_b	= prior$theta2 #miles
  
  ##################################
  ##### Select Starting Values #####
  ##################################
  
  my0	= summary(lm(y0~data$ps))
  my1	= summary(lm(y1~data$ps))
  
  missy1	= ismissy[seq(1,n*q,q)]
  missy2	= ismissy[seq(2,n*q,q)]
  missingindy	= list(missy1, missy2)
  
  Y[seq(1,n*q,q)][missy1==1]	= rnorm(sum(missy1==1), mean(Y[seq(1,n*q,q)], na.rm=T), 
                                     sd(Y[seq(1,n*q,q)], na.rm=T))
  Y[seq(2,n*q,q)][missy2==1]	= rnorm(sum(missy2==1), mean(Y[seq(2,n*q,q)], na.rm=T), 
                                     sd(Y[seq(2,n*1,q)], na.rm=T))	
  
  #Regression Coefficients
  if(is.null(starting$B) == F){
    B	= starting$B
  }
  if(is.null(starting$B) == T){
    B 	= rep(0, ncov*2)
  }
  
  #Spatial variance parameters, log transform for proposal
  if(is.null(starting$A) == T){
    initA1	= log(my0$sigma^2)
    initA2	= log(my1$sigma^2)
  }else{
    initA1	= starting$A
    initA2 	= starting$A
  }
  
  #Vector of log of residual variances (Psi)
  if(is.null(starting$psi) == T){
    initpsis	= log(c(my0$sigma^2, my1$sigma^2)*.1)
  }else{
    initpsis 	= log(rep(starting$psi, 2))
  }
  
  #Spatial Range parameters
  if(is.null(starting$theta) == T){
    inittheta	= c(1,1) # original scale
    inittheta	= logit(inittheta,thetaunif_a, thetaunif_b) #logit scale for proposals
  }else{
    inittheta 	= logit(starting$theta,thetaunif_a, thetaunif_b)
  }
  
  ##################################
  ##### Select Tuning Parameters ###
  ##################################
  
  ####Proposal SDS for A1 and A2
  A1propsds	= tuning$A
  A2propsds	= tuning$A
  
  ####Proposal SDS for log(diag(Psi))
  psipropsds	= rep(tuning$psi, 2)
  
  ###Proposal SDS for logit(theta) (aka, phi)
  thetapropsds	= rep(tuning$theta, 2)
  
  
  #################################
  ##### Define Some Functions #####
  #################################
  
  
  loglike_sp	= function(paramvals,Yvals){
    llike	= 0
    K1	= exp(paramvals[Aindex[1]])
    K2	= exp(paramvals[Aindex[2]])
    
    K	= createspatialsig(K1,K2,rho)
    
    if (!is.na(K)[[1]]){
      ### Untransform theta to the original scale
      thetatemp	= logitInv(paramvals[thetaindex], thetaunif_a, thetaunif_b)
      
      det	= mvCovInvLogDet(coords=coords, cov.model="exponential",
                           V=K, Psi=diag(exp(paramvals[psiindex]), nrow=q), 
                           theta=thetatemp, modified.pp=FALSE, SWM=TRUE)
      
      ## q inverse gamma priors for spatial variance
      ###Jacobian for log transformation: \sum log(sigmasq)
      llike	= llike+ sum(-(KIG_a+1)*paramvals[Aindex] -KIG_b/exp(paramvals[Aindex]) + paramvals[Aindex])
      
      ## q inverse gamma priors
      ###Jacobian for log transformation: \sum log(sigmasq)
      llike	= llike+ sum(-(psiig_a+1)*paramvals[psiindex] -psiig_b/exp(paramvals[psiindex]) + paramvals[psiindex])
      
      ## Uniform part for theta- this is the jacobian for the LogitInv transformation = (theta-a)(b-theta) / (b-a) 
      llike	= llike+ sum(log(thetatemp - thetaunif_a) + log(thetaunif_b - thetatemp))
      
      outlist		= list(llike, det$C, det$C.inv, det$log.det)
      names(outlist)	= c("ll","C","Cinv", "logdet")
    }# !isna(K)
    
    if (is.na(K)[[1]]){
      outlist		= list(-Inf, 1)
      names(outlist)	= c("ll", "notpd")
    }
    
    return(outlist)
  }	
  
  loglike_norm	= function(Yvals,Bvals,Cinvval,logdetval){
    ##Normal Part
    YXB	= (Yvals-X%*%Bvals)
    llike	= -.5*logdetval -.5*t(YXB)%*%Cinvval%*%YXB
    return(llike)
  }
  
  
  MHstep_sp	= function(index, paramvals, Bvals, Yvals, Cval, Cinvval, logdetval,currentllsp, currentllnorm){
    accept	= 0
    notpd		= 1
    llretsp	= currentllsp
    llretnorm	= currentllnorm
    currentll	= currentllsp+currentllnorm
    Cret		= Cval
    Cinvret	= Cinvval
    logdetret	= logdetval
    currentparams	= paramvals
    props			= currentparams
    props[index]	= rnorm(length(index), paramvals[index], propsds[index])
    llanddet		= loglike_sp(props,Yvals)
    Cprop		= llanddet$C
    Cinvprop	= llanddet$Cinv
    logdetprop	= llanddet$logdet
    llpropsp	= llanddet$ll
    
    if (llanddet$ll != -Inf){
      notpd	= 0
      llpropnorm	= loglike_norm(Yvals,Bvals,Cinvprop,logdetprop)
      llprop	= llpropsp+llpropnorm
      
      ratio		= exp(llprop-currentll)
      ratio[ratio>1] = 1
      
      if (runif(1)<=ratio){
        currentparams[index]	= props[index]
        llretsp	= llpropsp
        llretnorm	= llpropnorm
        Cret		= Cprop
        Cinvret	= Cinvprop
        logdetret	= logdetprop
        accept	= 1
      }
    }## llanddet$ll !=-Inf
    
    outlist		= list(currentparams, llretsp,llretnorm, accept, Cret, Cinvret, logdetret, notpd)
    names(outlist) 	= c("params", "llsp","llnorm", "accepted", "C", "Cinv", "logdet", "notpd")
    return(outlist)
  }
  
  updateBeta	= function(Cinvval, Yvals){
    S_beta	= solve(t(X)%*%Cinvval%*%X)
    Mu_beta  	= S_beta%*%t(X)%*%Cinvval%*%Yvals
    return(mvrnorm(1,Mu_beta,S_beta))
  }
  
  createspatialsig	= function(K1mat, K2mat, rho){
    sig		= matrix(0,q,q)
    sig[1,1]	= K1mat
    sig[2,2]	= K2mat
    sig[1,2]	= rho*sqrt(K1mat*K2mat)
    sig[2,1]	= sig[1,2]
    
    if(!is.positive.definite(sig)){
      sig	= NA
    }
    return(sig)
  }
  
  
  ########################################
  ##### Initialize Sampling Matrices #####
  ########################################
  
  params	= c(initA1,initA2, initpsis, inittheta)  #Everything here is on the transformed scale for proposals
  nparams	= length(params)
  propsds	= c(A1propsds, A2propsds, psipropsds, thetapropsds)
  Aindex	= 1:2
  psiindex	= 3:4
  thetaindex	= 5:6
  accepted 	= rep(0, length(propsds))
  
  binsize 	= 1000
  rho		= 0
  notpd		= 0
  
  saveindex = seq(1, nsamp, by = thin)
  samples	= matrix(NA, nrow=length(saveindex), ncol=nparams+1+length(B))
  dimnames(samples)[[2]]	= c("K[1,1]", "K[2,1]", "K[2,2]", paste("Psi", 1:q, sep=""), paste("Theta", 1:2, sep=""), 
                             paste("B0", 0:((p/2)-1), sep=""), paste("B1", 0:((p/2)-1), sep=""))
  ysims		= matrix(NA, nrow=length(saveindex), ncol=length(Y))
  
  ##Calculate initial log likelihood	
  initsp	= loglike_sp(params,Y)
  Cinv		= initsp$Cinv
  C		= initsp$C
  logdet 	= initsp$logdet
  llsp		= initsp$ll
  llnorm	= loglike_norm(Y,B,Cinv,logdet) 
  ll		= llsp+llnorm
  
  iterno	= 1
  donesampling= FALSE
  kk 		= 1
  
  while (donesampling==FALSE){
    B 	= updateBeta(Cinv,Y)
    llnorm= loglike_norm(Y,B,Cinv,logdet)
    ll	= llsp+llnorm
    
    ######## PROPOSE #########
    ###### Metropolis steps for spatial parameters #########
    #a0 parameters
    paramindex	= c(Aindex[1], psiindex[1], thetaindex[1])
    mhstep	= MHstep_sp(paramindex,params,B,Y, C, Cinv, logdet, llsp, llnorm)
    params	= mhstep$params
    llsp		= mhstep$llsp
    llnorm	= mhstep$llnorm
    Cinv		= mhstep$Cinv
    C		= mhstep$C
    logdet	= mhstep$logdet
    accepted[paramindex] = accepted[paramindex]+mhstep$accepted
    notpd		= notpd+mhstep$notpd
    
    #a1parameters
    paramindex	= c(Aindex[2], psiindex[2], thetaindex[2])
    mhstep	= MHstep_sp(paramindex,params,B,Y, C, Cinv, logdet, llsp, llnorm)
    params	= mhstep$params
    llsp		= mhstep$llsp
    llnorm	= mhstep$llnorm
    Cinv		= mhstep$Cinv
    C		= mhstep$C
    logdet	= mhstep$logdet
    accepted[paramindex]	= accepted[paramindex]+mhstep$accepted
    notpd		= notpd+mhstep$notpd
    
    ll 		= llsp+llnorm	
    
    #Simulate Missing Y
    Sy_ystar	= C[!ismissy, ismissy]
    Sy_y		= C[!ismissy,!ismissy]
    Sy_y_inv	= solve(Sy_y)
    Systar_ystar= C[ismissy,ismissy]
    
    m_pred	= X[ismissy,]%*%B + t(Sy_ystar)%*%Sy_y_inv%*%(Y[!ismissy] - X[!ismissy,]%*%B)
    S_pred	= Systar_ystar - t(Sy_ystar)%*%Sy_y_inv%*%Sy_ystar
    
    Y[ismissy]	= mvrnorm(1, m_pred,S_pred)
    
    llnorm= loglike_norm(Y,B,Cinv,logdet)  	 	
    ### Create output variables on the original scale
    
    K1out	= exp(params[Aindex[1]])
    K2out	= exp(params[Aindex[2]])
    Kout	= createspatialsig(K1out, K2out,rho)
    
    thetaout		= logitInv(params[thetaindex], thetaunif_a, thetaunif_b)
    
    if(iterno%%thin==0){
      samples[kk,]	= c(Kout[lower.tri(diag(1,q), TRUE)], exp(params[psiindex]), thetaout, B)
      ysims[kk,]		= Y
      kk = kk+1
    }
    
    if (!is.null(outputbinsize) & iterno%%outputbinsize==0){
      out  	<- list()
      out$samples <- samples[1:(kk-1), ]
      out$y0	<- ysims[1:(kk-1),seq(1,n*q,q)]
      out$y1 	<- ysims[1:(kk-1),seq(2,n*q,q)]
      out$coords 	<- coords
      out$trt	<- data$a
      out$formula <- formula
      save(out, file=outputfilename)
      rm(out)
    }
    
    if (iterno%%binsize==0){
      print(paste("Iteration:", iterno))
      print(c("Accepted:", round(accepted/iterno,3)))
      print(c("Number of not PD:", notpd))
    }
    
    if (iterno>=nsamp){donesampling=TRUE}	
    
    iterno = iterno+1
    
  }#while donesampling==FALSE
  
  out		<- list()
  out$samples <- samples
  out$y0	<- ysims[,seq(1,n*q,q)]
  out$y1 	<- ysims[,seq(2,n*q,q)]
  out$coords 	<- coords
  out$trt	<- data$a
  return(out)
}
