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

SNIRT_D_RSM = function(y, J, nchain=2, Update=F, mod=NULL, epoch=300, nbatch=min(c(nrow(y),30)), iters=(epoch*N/nbatch)+1, burnin=round((iters-1)/3), thin=max(1,round((iters-burnin-1)/1000))){
  
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
  
  ptau = function(y,tp,tf,tau,St){
    
    N = length(tp)
    J = length(tau) +1
    
    prob = array(1, dim=c(N,N,J))
    m_tp = matrix(tp,N,N, byrow=F) 
    m_tf = matrix(tf,N,N, byrow=T)
    m_t = array(m_tp + m_tf, dim=c(N,N,J-1))
    
    m_I = array(rep(1:(J-1),each=N^2), dim=c(N,N,J-1))
    mc_tau = array(rep(cumsum(tau),each = N^2),dim=c(N,N,J-1) )
    prob[,,2:J] = exp((m_I *m_t - mc_tau) ) 
    c_Prob = array(rowSums(prob, dims = 2), dim=c(N,N,J))
    
    m_I = array(rep(1:J,each=N^2), dim=c(N,N,J))  
    m_y = array(y, dim=c(N,N,J))
    Prob = (m_y==m_I)*prob/c_Prob
    sum_Prob = rowSums(Prob, dims=2)
    diag(sum_Prob) = 1
    
    (-sum(tau^2)/(2*St))+sum(log(sum_Prob))
  }
  
  
  pthetaP = function(y,tp,tf,tau,tfn,STT){
    
    N = length(tf) 
    J = length(tau) +1
    
    ## Tp
    tpf = (tp + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf - mc_tau))
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == (m_I+1)
    sum_Prob_p = apply(Prob * m_y,1,sum)
    
    tmp = (tp^2*STT[2,2] - 2 * tp * tfn * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    (-0.5*(tmp)+sum(log(sum_Prob_p)))
  }
  
  pthetaF = function(y,tp,tf,tau,tpn,STT){
    
    N = length(tp) 
    J = length(tau) +1
    
    ## Tf
    tpf = (tp + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf - mc_tau))
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == (m_I+1)
    sum_Prob_f = apply(Prob * m_y,1,sum)
    
    tmp = (tf^2*STT[1,1] - 2 * tpn * tf * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    (-0.5*(tmp)+sum(log(sum_Prob_f)))
  }
  
  plambda = function(l,tp,tf,la,v0,d0){
    
    N = length(tf)
    
    STT = matrix(c(l[1], -la*l[1], -la*l[1], l[1]*la^2 + l[2]),2,2)
    tmp = (tf^2*STT[1,1] + tp^2*STT[2,2] - 2 * tp * tf * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    
    (sum((-v0-1)*log(l)-d0/l)+
        -(N/2)*log(det(STT))-0.5*sum(tmp))
  }
  
  pla = function(la,tp,tf,l,v){
    
    Sl = v*rep(l[2:length(l)],1:(length(l)-1))
    
    STT = matrix(c(l[1], -la*l[1], -la*l[1], l[1]*la^2 + l[2]),2,2)
    tmp = (tf^2*STT[1,1] + tp^2*STT[2,2] - 2 * tp * tf * STT[1,2])/(STT[1,1]*STT[2,2]-STT[1,2]^2)
    
    (-0.5*sum(la^2/Sl)
      -(N/2)*log(det(STT))-0.5*sum(tmp))
  }  
  
  Sigmat = function(l, la){
    STT = matrix(c(l[1], -la*l[1], -la*l[1], l[1]*la^2 + l[2]),2,2)
    return(STT)
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
  
  r = 0
  N = nrow(y)
  K = J
  Batch = seqlast(1,N+1,by=nbatch)
  BatchLen = length(Batch)-1
  
  ## Initial value
  {
    ### Chain
    tp = array(scale(apply(y,1,sum)),dim = c(N, iters, nchain))
    tf = array(scale(apply(y,2,sum)),dim = c(N, iters, nchain))
    tau = array(seq(-2,2,length.out=J-1),dim = c(J-1, iters, nchain))
    Stp = array(1,dim = c(1, iters, nchain))
    Stf = array(1,dim = c(1, iters, nchain))
    St = 0.5
    A = array(rep(0,1),dim=c(1,iters,nchain))
    L = array(rep(1,2),dim=c(2,iters,nchain))
  }
  
  {
    Jtau = rep(0.1,J-2)
    JL = rep(0.2,2)
    JA = 0.2
    Jtf = rep(0.5,N)
    Jtp = rep(0.5,N)
  }
  
  time1 = Sys.time()
  
  # Parallel backend
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(min(n.cores,nchain), type = "PSOCK")
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
        
        ### tau
        {
          proposal = mvtnorm::rmvnorm(1,tau[1:(J-2),iter,c],diag(Jtau,J-2))
          proposal = sort(c(proposal,-sum(proposal)))
          R = exp(ptau(y[batch,batch],tp[batch,iter,c],tf[batch,iter,c],proposal,St)
                 -ptau(y[batch,batch],tp[batch,iter,c],tf[batch,iter,c],tau[,iter,c],St) )
          if(is.na(R)==T){R=0}
          accept = rbinom(1,1,min(1,R))
          tau[,iter+1,c] = ifelse(rep(accept,K-1),proposal,tau[,iter,c])
        }
        Jtau = apply(tau[,1:(iter+1),c],1,var)[1:(J-2)] + rep(0.0001,J-2)
        
        ### var
        {
          St = ss3(tau[,iter+1,c],3,10)
        }
        
        ### L
        proposal2 = L[,iter,c]
        for(l in 1:2){
          ### Lambda1, 3
          b3 = rgamma(1, shape=5.5, rate = 2*L[l,iter,c]/(2+L[l,iter,c]))
          temp = proposal2
          proposal = rlnorm(1,log(L[l,iter,c]),JL[l])
          proposal2[l] = proposal
          R = exp(plambda(proposal2,tp[,iter,c],tf[,iter,c],A[,iter,c],3,b3)+log(proposal2[l])-
                    (plambda(temp,tp[,iter,c],tf[,iter,c],A[,iter,c],3,b3)+log(temp[l])))
          R[is.na(R)] = 0
          accept = rbinom(1,1,min(1,R))
          proposal2[l] = ifelse(accept,proposal,L[l,iter,c]) 
        }
        L[,iter+1,c] = proposal2
        JL = apply(log(L[,1:(iter+1),c]),1,sd)+0.00001
        
        ### A
        proposal = rnorm(1,A[,iter,c],JA)
        R = exp(pla(proposal,tp[,iter,c],tf[,iter,c],L[,iter+1,c],5)-
                  pla(A[,iter,c],tp[,iter,c],tf[,iter,c],L[,iter+1,c],5))
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        A[,iter+1,c] = ifelse(accept, proposal, A[,iter,c])
        
        JA = sd(A[,1:(iter+1),c])+0.00001
        
        STT = Sigmat(L[,iter+1,c],A[,iter+1,c])
        
        for(n in 1:N){
          n_theta = unique(c(n, batch))[-1]
          {
            ### thetaF
            proposal = rnorm(1,tf[n,iter,c],Jtf[n])
            R = exp(pthetaF(y[n_theta,n],tp[n_theta,iter,c],proposal,tau[,iter+1,c],tp[n,iter,c],STT)-
                    pthetaF(y[n_theta,n],tp[n_theta,iter,c],tf[n,iter,c],tau[,iter+1,c],tp[n,iter,c],STT))
            if(is.na(R)==T){R=0}
            accept = rbinom(1,1,min(1,R))
            tf[n,iter+1,c] = ifelse(accept,proposal,tf[n,iter,c])
            
            Jtf[n] = sd(tf[n,1:(iter+1),c]) + 0.0001
          }
        }
        
        for(n in 1:N){
          n_theta = unique(c(n, batch))[-1]
          {
            ### thetaP
            proposal = rnorm(1,tp[n,iter,c],Jtp[n])
            R = exp(pthetaP(y[n,n_theta],proposal,tf[n_theta,iter+1,c],tau[,iter+1,c],tf[n,iter+1,c],STT)-
                    pthetaP(y[n,n_theta],tp[n,iter,c],tf[n_theta,iter+1,c],tau[,iter+1,c],tf[n,iter+1,c],STT))
            if(is.na(R)==T){R=0}
            accept = rbinom(1,1,min(1,R))
            tp[n,iter+1,c] = ifelse(accept,proposal,tp[n,iter,c])
            
            Jtp[n] = sd(tp[n,1:(iter+1),c]) + 0.0001
          }
        }
        iter = iter + 1

      }
    }
    out = list(tau = tau[,,c],
               tp = tp[,,c],
               tf = tf[,,c], 
               A = A[,,c],
               L = L[,,c])
    return(out)
    
  }
  
  stopCluster(my.cluster)
  
  for(i in 1:nchain){
    tp[,,i] = tmp[[i]]$tp
    tf[,,i] = tmp[[i]]$tf
    tau[,,i] = tmp[[i]]$tau
    A[,,i] = tmp[[i]]$A
    L[,,i] = tmp[[i]]$L
  }
  
  # EAP
  node = seq(burnin+thin,iters,by=thin)
  
  ID = row.names(y)
  
  tp_out = eap(tp, node); names(tp_out) = ID
  tf_out = eap(tf, node); names(tf_out) = ID 
  tau_out = eap(tau, node) 
  L_out = eap(L, node) 
  A_out = eap1d(A, node)
  STT = Sigmat(L_out,A_out)
  Stp_out = STT[1,1]
  Stf_out = STT[2,2]
  Rho_out = STT[1,2] /sqrt(Stp_out * Stf_out)
  
  tp_RR = RR(tp[,node,])
  tf_RR = RR(tf[,node,])
  tau_RR = RR(tau[,node,])
  L_RR = RR(L[,node,])
  A_RR = RR_s(A[,node,])
  
  out = list(
    Estimate = list(
      ThetaP = tp_out,
      ThetaF = tf_out,
      Tau = tau_out,
      Rho = Rho_out,
      Stp = Stp_out,
      Stf = Stf_out
    ),
    MCMC = list(
      ThetaP = tp,
      ThetaF = tf,
      Tau = tau,
      D = L,
      L = A
    ),
    RSquare = list(
      ThetaP = tp_RR,
      ThetaF = tf_RR,
      Tau = tau_RR,
      D = L_RR, 
      L = A_RR
    ),
    Node = node
  )
  return(out)
}

