##### Estimation #####
ProcrustesMatching = function(chain_output, Args){
  
  N=Args[['N']]
  V=Args[['V']]
  nchain=Args[['nchain']]
  K=Args[['K']]
  
  post_samples = NULL
  for(c in 1:nchain){
    post_samples = rbind(post_samples, as.matrix(chain_output[[c]][[1]]$samples))
  }
  n_var = colnames(post_samples)
  sel_Xi = grep('xi', n_var)
  
  len = dim(post_samples)[1]
  
  est_logL = NULL
  for(c in 1:nchain){
    tmp_logL = chain_output[[c]][[1]]$samples[,grep('logL_raw',n_var)] %>% array(dim=c(len/nchain, N,N,K))
    est_logL = c(est_logL, rowSums(tmp_logL,dims=1))
  }
  
  ref_point = which.max(est_logL) 
  post_Xi = post_samples[,sel_Xi] %>% array(.,dim=c(len,N,V))
  
  post_Xi_ref = post_Xi[ref_point,,]
  
  post_Xi_arr = array(0, dim=c(N, V, len))
  for(i in 1:len){
    Xi_sample = post_Xi[i,,] # c(len, N, V)
    proc = MCMCpack::procrustes(Xi_sample, post_Xi_ref, translation = T, dilation = T)
    post_Xi_arr[,,i] = proc$X.new # c(N, V, len)
  }
  
  est_Xi = rowMeans(post_Xi_arr, dim=2)
  
  return(est_Xi)
}

##### Model Fit #####
logL = function(y,thetaS, thetaR, tauj, lambda, xi){
  logL_outpput = 0
  cats = length(tauj)
  N = length(thetaS)
  r = matrix(0,N,N)
  for(i in 1: N){
    for(j in (1:N)[-i]){
      # Retrivel from knowledge
      ## Probability of y
      logit = thetaS[i] + thetaR[j] - lambda * sqrt(sum((xi[i,]- xi[j,])^2))- tauj
      py = exp(cumsum(c(0, logit)))
      py = py / sum(py)
      
      # Reponse of person n
      for(k in 1:(cats+1)){
        logL_outpput = logL_outpput + (y[i,j]==k)*log(py[k])
      }
    }
  }
  
  return(logL_outpput)
}

logL_InnerProd = function(y,thetaS, thetaR, tauj, lambda, xi){
  logL_outpput = 0
  cats = length(tauj)
  N = length(thetaS)
  r = matrix(0,N,N)
  for(i in 1: N){
    for(j in (1:N)[-i]){
      # Retrivel from knowledge
      ## Probability of y
      Dist = sum(xi[i,] * xi[j,])
      logit = thetaS[i] + thetaR[j] + lambda * Dist - tauj
      py = exp(cumsum(c(0, logit)))
      py = py / sum(py)
      
      # Reponse of person n
      for(k in 1:(cats+1)){
        logL_outpput = logL_outpput + (y[i,j]==k)*log(py[k])
      }
    }
  }
  
  return(logL_outpput)
}

logL_Proj = function(y,thetaS, thetaR, tauj, lambda, xi){
  logL_outpput = 0
  cats = length(tauj)
  N = length(thetaS)
  r = matrix(0,N,N)
  for(i in 1: N){
    for(j in (1:N)[-i]){
      # Retrivel from knowledge
      ## Probability of y
      Dist = sum(xi[i,] * xi[j,]) / sqrt(sum(xi[j,]^2))
      logit = thetaS[i] + thetaR[j] + lambda * Dist - tauj
      py = exp(cumsum(c(0, logit)))
      py = py / sum(py)
      
      # Reponse of person n
      for(k in 1:(cats+1)){
        logL_outpput = logL_outpput + (y[i,j]==k)*log(py[k])
      }
    }
  }
  
  return(logL_outpput)
}

# DIC
DIC_LSRRM = function(y, nchain=2, K=5, v=2, Dist='Dist'){
  # ppDIC
  DIC = LogL = NULL
  {
    N = nrow(y)
  }
  
  len = dim(post_samples)[1]
  ppdic = post_samples[,grep('logL', n_var_est)] %>% as.matrix %>% apply(.,1,sum) %>% sum()
  ppdic = ppdic/len
  
  if(Dist=='Dist'){
    logl = logL(y, Est_Theta[,1], Est_Theta[,2], Est_Tau, Est_Lambda, Est_Xi)
  }else if(Dist=='InnerProd'){
    logl = logL_InnerProd(y, Est_Theta[,1], Est_Theta[,2], Est_Tau, Est_Lambda, Est_Xi)
  }else if(Dist=='Proj'){
    logl = logL_Proj(y, Est_Theta[,1], Est_Theta[,2], Est_Tau, Est_Lambda, Est_Xi)
  }else{cat('Model is not specified correctly.\n'); break}
  
  pdic = 2*(logl-ppdic)
  DIC = c(DIC, -2*logl+2*pdic)
  LogL = c(LogL, -2*logl)
  
  out = list(LogL = LogL, DIC=DIC)
  return(out)
}

##### Sender Fit & Receiver Fit #####
# Sender Fit
Fit_s = function(y, thetaS, thetaR, tauj, lambda, xiS, xiR, Mod='Dist'){
  
  ProbS = function(thetaS, thetaR, tauj, lambda, xiS, xiR, Mod='Dist'){
    cats = length(tauj)
    N = length(thetaR)
    Prob = matrix(0,N,cats+1)
    for(i in 1: N){
      if(Mod=='Dist'){
        Dist = -sqrt(sum((xiS - xiR[i,])^2))
      }else if(Mod=='InnerProd'){
        Dist = sum(xiS * xiR[i,])
      }else if(Mod=='Proj'){
        Dist = sum(xiS * xiR[i,]) / sqrt(sum(xiR[i,]^2))
      }else{
        cat('Model is not specified correctly.\n')
        break
      }
      logit = thetaS + thetaR[i] + lambda * Dist - tauj
      py = exp(cumsum(c(0, logit)))
      py = py / sum(py)
      
      Prob[i,] = py
    }
    return(Prob)
  }
  
  K = length(tauj)+1
  N = length(thetaR) 
  
  Prob_S = ProbS(thetaS, thetaR, tauj, lambda, xiS, xiR, Mod='Dist') 
  m_y = matrix(0:(K-1), nrow=N, ncol=K, byrow=T)
  E_S = apply(Prob_S * m_y, 1, sum)
  m_E_S = matrix(E_S, nrow=N, ncol=K, byrow=F)
  V_S = apply(Prob_S * (m_y - m_E_S)^2,1,sum)
  y_S = y-1
  
  S_Fit = (sum((y_S - E_S)^2)/sum(V_S))
  
  return(S_Fit)
}

# Receiver fit
Fit_r = function(y, thetaS, thetaR, tauj, lambda, xiS, xiR, Mod='Dist'){
  
  ProbR = function(thetaS, thetaR, tauj, lambda, xiS, xiR, Mod='Dist'){
    cats = length(tauj)
    N = length(thetaS)
    Prob = matrix(0,N,cats+1)
    for(i in 1: N){
      if(Mod=='Dist'){
        Dist = -sqrt(sum((xiS[i,] - xiR)^2))
      }else if(Mod=='InnerProd'){
        Dist = sum(xiS[i,] * xiR)
      }else if(Mod=='Proj'){
        Dist = sum(xiS[i,] * xiR) / sqrt(sum(xiR^2))
      }else{
        cat('Model is not specified correctly.\n')
        break
      }
      logit = thetaS[i] + thetaR + lambda * Dist - tauj
      py = exp(cumsum(c(0, logit)))
      py = py / sum(py)
      
      Prob[i,] = py
    }
    return(Prob)
  }
  
  K = length(tauj)+1
  N = length(thetaS) 
  
  Prob_R = ProbR(thetaS, thetaR, tauj, lambda, xiS, xiR, Mod='Dist') 
  m_y = matrix(0:(K-1), nrow=N, ncol=K, byrow=T)
  E_R = apply(Prob_R * m_y, 1, sum)
  m_E_R = matrix(E_R, nrow=N, ncol=K, byrow=F)
  V_R = apply(Prob_R * (m_y - m_E_R)^2,1,sum)
  y_R = y-1
  
  R_Fit = (sum((y_R - E_R)^2)/sum(V_R))
  
  return(R_Fit)
}


data_gen_fit = function(thetaP, thetaF, tauj, lambda, xiS, xiR, SR='s', Mod='Dist'){
  if(SR=='s'){
    cats = length(tauj)
    N = length(thetaF)
    r = matrix(0,1,N)
    for(i in 1:1){
      for(j in 1:N){
        
        if(Mod=='Dist'){
          Dist = -sqrt(sum((xiS - xiR[j,])^2))
        }else if(Mod=='InnerProd'){
          Dist = sum(xiS * xiR[j,])
        }else if(Mod=='Proj'){
          Dist = sum(xiS * xiR[j,]) / sqrt(sum(xiR[j,]^2))
        }else{
          cat('Model is not specified correctly.\n')
          break
        }
        
        logit = thetaP + thetaF[j] + lambda*Dist - tauj
        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        x = rcat(1, py)
        
        r[i,j] = x
      }
    }
  }else if(SR=='r'){
    cats = length(tauj)
    N = length(thetaP)
    r = matrix(0,N,1)
    for(i in 1:N){
      for(j in 1:1){
        
        if(Mod=='Dist'){
          Dist = -sqrt(sum((xiS[i,] - xiR)^2))
        }else if(Mod=='InnerProd'){
          Dist = sum(xiS[i,] * xiR)
        }else if(Mod=='Proj'){
          Dist = sum(xiS[i,] * xiR) / sqrt(sum(xiR^2))
        }else{
          cat('Model is not specified correctly.\n')
          break
        }
        
        logit = thetaP[i] + thetaF + lambda*Dist - tauj
        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        x = rcat(1, py)
        
        r[i,j] = x
      }
    }
  }
  
  return(r)
}
