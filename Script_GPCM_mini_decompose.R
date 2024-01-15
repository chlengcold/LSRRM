library(psych)
library(pscl)
library(dplyr)
library(mvtnorm)
library(foreach)
library(doParallel)
library(ranger)
library(tidyverse)

seqlast = function(from, to, by){
  vec = do.call(seq, list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

SNIRT_D_GPCM = function(y, J, nchain=2, epoch=300, nbatch=min(c(nrow(y),30)), iters=(epoch*N/nbatch)+1, burnin=round((iters-1)/3), thin=round((iters-burnin-1)/1000)){
  
  # MCMC
  sb = function(x,a,b){
    new_b = rgamma(1,2.5+a,x*b/(x+b))
    return(new_b)
  }
  
  ss3 = function(x,a,b){
    I = length(x)
    s = pscl::rigamma(1,a+I/2,b+sum(x^2)/2)
    return(s)
  }
  
  ss2 = function(x,a,b){
    I = length(x)
    s = pscl::rigamma(1,a+I,b+sum(x))
    return(s)
  }
  
  ptau = function(y,a,tp,tf,tau,St){
    
    N = length(tp)
    J = length(tau) +1
    
    tpf = (a*tp + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf) - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob * m_y,1,sum)
    
    (-sum(tau^2)/(2*St))+sum(log(sum_Prob))
  }
  
  pa = function(y,a,tp,tf,tau,Sa){
    
    N = length(tp)
    J = length(tau)+1
    
    tpf = (a*tp + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf) - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob = apply(Prob * m_y,1,sum)
    
    -a/Sa+sum(log(sum_Prob))
  }
  
  pthetaP = function(y,a,tp,tf,tau,tfn,STT){
    
    N = length(tf) 
    J = ncol(tau) +1
    
    ## Tp
    tpf = (a*tp + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    mc_tau = t(apply(tau,1, cumsum))
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf) - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = NULL
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob_p = apply(Prob * m_y,1,sum)
    
    tmp = (tp^2*STT[2,2] - 2 * tp * tfn * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    (-0.5*(tmp)+sum(log(sum_Prob_p)))
  }
  
  pthetaF = function(y,a,tp,tf,tau,tpn,STT){
    
    N = length(tp) 
    J = length(tau) +1
    
    ## Tf
    tpf = (a*tp + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf) - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = NULL
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob_f = apply(Prob * m_y,1,sum)
    
    tmp = ( tf^2*STT[1,1] - 2 * tpn * tf * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    (-0.5*(tmp)+sum(log(sum_Prob_f)))
  }
  
  pla = function(la,tp,tf,v){
    STT = matrix(c(1, rep(-la,2), 1),2,2)
    tmp = (tf^2*STT[1,1] + tp^2*STT[2,2] - 2 * tp * tf * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    
    Sl = v*(1-la^2)
    
    (-0.5*sum(la^2/Sl)
      -(N/2)*log(det(STT))-0.5*sum(tmp))
  }  
  
  # Estimate
  eap = function(x,node){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,node,c]/nchain 
    }
    out = apply(tmp, 1, mean)
    return(out)
  } 
  
  eap1d = function(x,node){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,node,c]/nchain 
    }
    out = mean(tmp)
    return(out)
  } 
  
  MCMC = function(x){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,,c]/nchain 
    }
    return(tmp)
  } 
  
  RR_s = function(x){
    len = dim(x)[1]
    nchain = dim(x)[2]
    B = NULL
    W = NULL
    B_tmp = NULL
    W_tmp = NULL
    for(c in 1:nchain){
      B_tmp = c(B_tmp, mean(x[,c]))
      W_tmp = c(W_tmp, var(x[,c]))
    }
    B = c(B, len*var(B_tmp))
    W = c(W, mean(W_tmp))
    
    Var = (len-1)*W/len + B/len
    sqrt(Var/W)
  }
  
  RR = function(x){
    items = dim(x)[1]
    len = dim(x)[2]
    nchain = dim(x)[3]
    B = NULL
    W = NULL
    for(i in 1:items){
      B_tmp = NULL
      W_tmp = NULL
      for(c in 1:nchain){
        B_tmp = c(B_tmp, mean(x[i,,c]))
        W_tmp = c(W_tmp, var(x[i,,c]))
      }
      B = c(B, len*var(B_tmp))
      W = c(W, mean(W_tmp))
    }
    Var = (len-1)*W/len + B/len
    sqrt(Var/W)
  }
  
  RR_m = function(x){
    items = dim(x)[1]
    Dim = dim(x)[2]
    len = dim(x)[3]
    nchain = dim(x)[4]
    B = NULL
    W = NULL
    for(d in 1:Dim){
      B_tmp2 = NULL
      W_tmp2 = NULL
      for(i in 1:items){
        B_tmp = NULL
        W_tmp = NULL
        for(c in 1:nchain){
          B_tmp = c(B_tmp, mean(x[i,d,,c]))
          W_tmp = c(W_tmp, var(x[i,d,,c]))
        }
        B_tmp2 = c(B_tmp2, len*var(B_tmp))
        W_tmp2 = c(W_tmp2, mean(W_tmp))
      }
      B = cbind(B,B_tmp2)
      W = cbind(W,W_tmp2)
    }
    Var = (len-1)*W/len + B/len
    Out = sqrt(Var/W)
    colnames(Out) = paste('Dim',1:Dim)
    
    return(Out)
  }
  
  r = 0
  N = nrow(y)
  K = J
  Batch = seqlast(1,N+1,by=nbatch)
  BatchLen = length(Batch)-1
  
  ## Initial value
  {
    ### Chain
    a = array(1, dim = c(N, iters, nchain))
    Sa = 0.5
    tp = array(scale(apply(y,1,sum)),dim = c(N, iters, nchain))
    tf = array(scale(apply(y,2,sum)),dim = c(N, iters, nchain))
    tau = array(matrix(seq(-2,2,length.out=J-1),N,J-1,byrow=T),dim = c(N, J-1, iters, nchain))
    St = rep(0.5, N)
    A = array(rep(0,1),dim=c(1,iters,nchain))
  }
  
  {
    Ja = rep(0.5,N)
    Jtau = matrix(rep(0.1,J-2),nrow=N,ncol=J-2,byrow=T)
    JA = 0.2
    Jtf = rep(0.5,N)
    Jtp = rep(0.5,N)
  }

  # Parallel backend
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(min(nchain,n.cores), type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  ## Estimation
  tmp = foreach(c = 1:nchain, .combine = list, .multicombine=TRUE) %dopar%{
    iter = 1
    for(e in 1:epoch){
      for(b in 1:BatchLen){
        posi = sample(1:N, replace = F)
        batch_seq = Batch[b]:(Batch[b+1]-1) 
        batch = posi[batch_seq]
      
        Sr = ss3(A[,iter,c],3,10)
        
        ### r
        ### A
        proposal = rnorm(1,A[,iter,c],JA)
        R = exp(pla(proposal,tp[,iter,c],tf[,iter,c],5)-
              pla(A[,iter,c],tp[,iter,c],tf[,iter,c],5))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        A[,iter+1,c] = ifelse(accept, proposal, A[,iter,c])

        JA = sd(A[,1:(iter+1),c])+0.00001
        
        STT = matrix(c(1, rep(-A[,iter+1,c],2), 1),2,2)
        
        for(n in 1:N){
          n_theta = unique(c(n, batch))[-1]
          {
            ### a
            proposal = rlnorm(1,log(a[n,iter,c]),Ja[n])
            R = exp(log(proposal)+pa(y[n_theta,n],proposal,tp[n_theta,iter,c],tf[n,iter,c],tau[n,,iter,c],Sa)
                    -log(a[n,iter,c])-pa(y[n_theta,n],a[n,iter,c],tp[n_theta,iter,c],tf[n,iter,c],tau[n,,iter,c],Sa) )
            R[is.na(R)] = 0
            accept = rbinom(1,1,min(1,R))
            a[n,iter+1,c] = ifelse(accept,proposal,a[n,iter,c])
            
            Ja[n] = sd(log(a[n,1:(iter+1),c])) + 0.0001
            
            ### tau
            proposal = mvtnorm::rmvnorm(1,tau[n,1:(J-2),iter,c],diag(Jtau[n,],J-2))
            proposal = sort(c(proposal,-sum(proposal)))
            R = exp(ptau(y[n_theta,n],a[n,iter+1,c], tp[n_theta,iter,c],tf[n,iter,c],proposal,St[n])
                   -ptau(y[n_theta,n],a[n,iter+1,c], tp[n_theta,iter,c],tf[n,iter,c],tau[n,,iter,c],St[n]) )
            if(is.na(R)==T){R=0}
            accept = rbinom(1,1,min(1,R))
            tau[n,,iter+1,c] = ifelse(rep(accept,K-1),proposal,tau[n,,iter,c])
            
            Jtau[n,] = apply(tau[n,,1:(iter+1),c],1,var)[1:(J-2)] + rep(0.0001,J-2)
            St[n] = ss3(tau[n,,iter+1,c],3,10)
            
            ### thetaF
            proposal = rnorm(1,tf[n,iter,c],Jtf[n])
            R = exp(pthetaF(y[n_theta,n],a[n,iter+1,c], tp[n_theta,iter,c],proposal,tau[n,,iter+1,c],tp[n,iter,c],STT)-
                    pthetaF(y[n_theta,n],a[n,iter+1,c], tp[n_theta,iter,c],tf[n,iter,c],tau[n,,iter+1,c],tp[n,iter,c],STT))
            if(is.na(R)==T){R=0}
            accept = rbinom(1,1,min(1,R))
            tf[n,iter+1,c] = ifelse(accept,proposal,tf[n,iter,c])
            
            Jtf[n] = sd(tf[n,1:(iter+1),c]) + 0.0001
          }
        }
        Sa = ss2((a[,iter+1,c]),3,8)
        
        for(n in 1:N){
          n_theta = unique(c(n, batch))[-1]
          {
            ### thetaP
            proposal = rnorm(1,tp[n,iter,c],Jtp[n])
            R = exp(pthetaP(y[n,n_theta],a[n_theta,iter+1,c],proposal,tf[n_theta,iter+1,c],tau[n_theta,,iter+1,c],tf[n,iter+1,c],STT)-
                pthetaP(y[n,n_theta],a[n_theta,iter+1,c],tp[n,iter,c],tf[n_theta,iter+1,c],tau[n_theta,,iter+1,c],tf[n,iter+1,c],STT))
            if(is.na(R)==T){R=0}
            accept = rbinom(1,1,min(1,R))
            tp[n,iter+1,c] = ifelse(accept,proposal,tp[n,iter,c])
            
            Jtp[n] = sd(tp[n,1:(iter+1),c]) + 0.0001
          }
        }
      
        iter = iter + 1
      }
    }
    out = list(tau = tau[,,,c],
               A = A[,,c],
               a = a[,,c],
               tp = tp[,,c],
               tf = tf[,,c])
    return(out)
    
  }
  
  stopCluster(my.cluster)
  
  for(i in 1:nchain){
    tp[,,i] = tmp[[i]]$tp
    tf[,,i] = tmp[[i]]$tf
    a[,,i] = tmp[[i]]$a
    A[,,i] = tmp[[i]]$A
    tau[,,,i] = tmp[[i]]$tau
  }
  
  # EAP
  node = seq(burnin+thin,iters,by=thin)
  
  ID = row.names(y)
  a_out = eap(a, node); names(a_out) = ID 
  tp_out = eap(tp, node); names(tp_out) = ID
  tf_out = eap(tf, node); names(tf_out) = ID 
  tau_out = matrix(0,nrow=N,ncol=(J-1))
  for(j in 1:(J-1)){
    tau_out[,j] = eap(tau[,j,,], node) 
  }
  A_out = eap1d(A, node)
  Rho_out = -A_out
  
  # tp_MCMC = MCMC(tp)
  # tf_MCMC = MCMC(tf)
  # tau_MCMC = MCMC(tau)
  
  tp_RR = RR(tp[,node,])
  tf_RR = RR(tf[,node,])
  tau_RR = RR_m(tau[,,node,])
  a_RR = RR(a[,node,])
  A_RR = RR_s(A[,node,])
  
  out = list(
    Estimate = list(
      a = a_out,
      ThetaP = tp_out,
      ThetaF = tf_out,
      Tau = tau_out,
      Rho = Rho_out
    ),
    MCMC = list(
      a = a,
      ThetaP = tp,
      ThetaF = tf,
      Tau = tau,
      L = A
    ),
    RSquare = list(
      a = a_RR,
      ThetaP = tp_RR,
      ThetaF = tf_RR,
      Tau = tau_RR,
      L = A_RR
    ),
    Node = node
  )
  return(out)
}