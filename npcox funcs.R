
beta1 = function(t){t^2}
beta2 = function(t){1-t}

npsimu = function(n, cenpara = NULL){
  
  cva  = data.frame(Z1 = runif(n,0.01,1), Z2 = runif(n,0.01,2))
  Z    = as.matrix(cva)
  if(is.null(cenpara)){ cenpara = 4 }
  
  cen   = runif(n,0,cenpara)
  Ti    = array();X = array()
  for(i in 1:n){
    tem   = log(runif(1,0,1))
    lamt  = function(t){  exp(beta1(t)*Z[i,1] + Z[i,2]*beta2(t) - 1)  }
    Ft    = function(t){  integrate(lamt,0,t)$value + tem  }
    Ti[i] = ifelse(Ft(0)*Ft(10)<0,unlist(uniroot(Ft,c(0,10),tol = 1e-4)$root),10)
    X[i]  = min(Ti[i],cen[i])
  }
  delta   = as.numeric(Ti < cen)
  cenrate = 1 - sum(delta)/n
  
  data = list(cva = cva, delta = delta, obstime = X)
  return(data)
}


spsimu = function(n, cenpara = NULL){
  
  cva  = data.frame(Z1 = runif(n,0.01,3), Z2 = runif(n,0.01,2))
  Z    = as.matrix(cva)
  if(is.null(cenpara)){ cenpara = 6 }
  
  cen   = c()
  for(i in 1:n){cen[i] = min(runif(1,0,cenpara), 2)}
  Ti    = array();X = array()
  for(i in 1:n){
    tem   = log(runif(1,0,1))
    lamt  = function(t){  exp(1*Z[i,1] + Z[i,2]*beta2(t) - 1)  }
    Ft    = function(t){  integrate(lamt,0,t)$value + tem  }
    Ti[i] = ifelse(Ft(0)*Ft(10)<0,unlist(uniroot(Ft,c(0,10),tol = 1e-4)$root),10)
    X[i]  = min(Ti[i],cen[i])
  }
  delta   = as.numeric(Ti < cen)
  cenrate = 1 - sum(delta)/n
  
  data = list(cva = cva, delta = delta, obstime = X)
  return(data)
}


# temp = res; xrange = c(650,3200); CIlevel = 0.95
npplot = function(temp, xrange = NULL, CIlevel = 0.95){
  # print("Please decide the partition of a figure with par(...)")
  temp1 = temp
  ind = sort(unique(c(which(is.na(temp$temporal_coef[,1])),which(is.na(temp$temporal_coef[,1])))))
  if(length(ind)>0){
    temp1$temporal_coef = temp$temporal_coef[-ind,]
    temp1$temporal_coef_SEE = temp$temporal_coef_SEE[-ind,]
    temp1$data = temp$data[-ind,]
    temp  = temp1
  }
  
  coef = as.matrix(temp$temporal_coef)
  see  = as.matrix(temp$temporal_coef_SEE)
  data = as.data.frame(temp$data)
  obs  = data$obstime
  h    = temp$bandwidth
  hmin = max(which(obs<min(obs)+h))
  hmax = min(which(obs>max(obs)-h)) - 1
  
  CIquan = qnorm(1-(1-CIlevel)/2)
  
  r    = ncol(coef)
  n    = nrow(coef)
  if(!is.null(xrange)){
    hmin = max(which(obs<xrange[1]))
    hmax = min(which(obs>xrange[2])) - 1
  }
  cc1  = c(hmin:hmax)
  
  for(j in 1:r){
    uppci = coef[cc1,j]+CIquan*see[cc1,j]
    lowci = coef[cc1,j]-CIquan*see[cc1,j]
    if(sum(is.na(see)) == n*r){
      upp   = max(coef[cc1,j]) + 0.2*abs(max(coef[cc1,j]))
      low   = min(coef[cc1,j]) - 0.2*abs(min(coef[cc1,j]))
    }else{
      upp   = max(uppci) + 0.2*max(uppci)
      low   = min(lowci) - 0.2*min(lowci)
    }
    upp = max(upp,  abs(upp)*0.1)
    low = min(low, -abs(low)*0.1)
    plot(obs[cc1], coef[cc1,j],type = "l",main = colnames(coef)[j], xlab = "Time",
         ylab = "Time-varying coefficients", ylim = c(low, upp))
    lines(obs[cc1], uppci, col = 'blue')
    lines(obs[cc1], lowci, col = 'blue')
    abline(h = 0, col = "red")
  }
}



spcox = function(cva_cons, cva_time, delta, obstime, SE = FALSE, bandwidth = FALSE, resamp = 100){
  
  # PH part: β(t) estimation
  # Bandwidth assignment or selection function
  Kh  = function(h,x){ifelse(abs(x/h)<1, 0.75*(1-(x/h)^2),0)}
  fm2 = function(x){x = array(x,dim = c(length(x),1)); return(x%*%t(x))}
  hmx = function(){ min(which(Xs>max(obstime)-h)) }# floor(.9*n)min(which(Xs>max(Xs)-h)) - 1
  
  PE   = function(bhat,Zs,ds,Xs){
    # prediction error calculation, refer to Lutian(2005)
    bz  = array(0, dim = c(n,n))
    for(i in 1:n){
      for(j in 1:n){
        bz[i,j] = sum(bhat[i,]*Zs[j,])
      }
    }
    pie = array(0, n)
    for(l in 1:n){
      pie[l] = -( bz[l,l] - log( sum(exp(bz[l,l:n])) ) )
    }
    return(sum(pie))
  }
  
  band = function(){
    
    if(bandwidth){
      # print("--------------------------------------------------------------------")
      # print(paste("Note: The bandwidth is predesigned as ", bandwidth, sep = ""))
      h = bandwidth
    }else{
      # Prediction error for bandwidth selection
      # print("--------------------------------------------------------------------")
      # print("Note: No predesigned bandwidth, commencing bandwidth selection procedure")
      
      diff  = Xs
      for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1] }
      diff  = sort(diff,decreasing = T)
      hmi   = diff[5]
      sep   = (quantile(Xs,0.85)/3 - Xs[1])/30
      ha    = seq(max(hmi,Xs[1]), quantile(Xs,0.85), by = sep)
      PErec = array(0, length(ha))
      for(l in 1:length(ha)){
        h    = ha[l]
        #print(paste("Calculating PE for h = ",h, sep = ""))
        bhat = esti(h,Zs,ds,Xs)
        PErec[l] = PE(bhat,Zs,ds,Xs)
      }
      # Bandwidth selection
      h = ha[which(PErec == min(PErec))]
      print(paste("Note: The selected bandwidth is ", h, sep = ""))
    }
    return(h)
  }
  
  cons = function(bhat,h,Zs,ds,Xs){
    s = r
    
    hmin = max(which(Xs<min(Xs)+h))
    hmax = hmx()-2
    hmax1 = floor(hmax*0.95)
    
    bhatcheck = sum(is.na(bhat))
    if(bhatcheck > 0){
      ind1 = which(is.na(bhat[,1]))
      for(l in ind1){
        bhat[l,] = bhat[l-1,]
      }
    }
    
    # esti of the partly constant effect
    S0 = array(0, n)
    S1 = array(0, dim = c(n,r))
    S2 = array(0, dim = c(n,r,r))
    bz = array(0, dim = c(n,n))
    for(j in 1:n){
      for(i in 1:n){
        bz[j,i] = sum(bhat[j,]*Zs[i,])
      }
    }
    for(j in 1:n){
      S0[j] = sum(exp(bz[j,j:n]))
      if(j == n){
        S1[j,]  = exp(bz[j,j:n])*Zs[n,]
        S2[j,,] = exp(bz[j,j:n])*fm2(Zs[n,])
      }else{
        for(i in j:n){
          S1[j,]  = S1[j,]  + exp(bz[j,i])*Zs[i,]
          S2[j,,] = S2[j,,] + exp(bz[j,i])*fm2(Zs[i,])
        }
      }
    }
    Vbt = array(0, dim = c(n,r,r))
    for(j in 1:n){ Vbt[j,,] = S2[j,,]/S0[j] - fm2(S1[j,]/S0[j] ) }
    
    Ibt  = array(0, dim = c(n,r,r))
    Iinv = array(0, dim = c(n,r,r))
    for(j in hmin:hmax1){
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[j]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      for(i in min_i:max_i){
        Ibt[j,,] = Ibt[j,,] + Vbt[i,,]*Ker[i]*ds[i]/n
      }
      Iinv[j,,] = solve(Ibt[j,,])
    }
    Jt    = array(0, dim = c(n,r1,r1))
    for(i in hmin:hmax1){ Jt[i,,] = solve(Iinv[i,1:r1,1:r1]) }# v(c(Jt))
    diff  = Xs
    for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1]}
    Jtint = array(0,dim = c(r1,r1))
    for(i in hmin:hmax1){ Jtint = Jtint + Jt[i,,]*diff[i] }
    wop   = array(0, dim = c(n,r1,r1))
    Jtint_inv = solve(Jtint)
    for(i in hmin:hmax1){ wop[i,,] = Jtint_inv%*%Jt[i,,] }
    beta1 = rep(0,r1)
    for(i in hmin:hmax1){ beta1 = beta1 + wop[i,,]%*%bhat[i,1:r1]*diff[i]}
    beta1 = c(beta1)
    beta1_sd = sqrt(diag(Jtint_inv))
    
    # Ub = function(beta){
    #   temp = array(NA,dim = c(n,s))
    #   S0 = array(n); S0pie = array(n)
    #   for(i in n:1){
    #     S0pie[i] = exp(c(beta%*%Zs[i,]))
    #     if(i == n){ S0[i] = S0pie[i] }else{ S0[i] = S0pie[i] + S0[i+1] } 
    #   }
    #   S1 = array(0,dim = c(n,s)); S1pie = array(0,dim = c(n,s))
    #   for(i in n:1){
    #     S1pie[i,] = S0pie[i]*Zs[i,]
    #     if(i == n){ S1[i,] = S1pie[i,] }else{ S1[i,] = S1pie[i,] + S1[i+1,] }
    #   }
    #   for(i in 1:n){
    #     temp[i,] = ds[i]*(Zs[i,] - S1[i,]/S0[i])
    #   }
    #   return(colSums(temp))
    # }
    # bhat = nleqslv(rep(1,s), Ub)$x
    # 
    # Ubb = function(beta){
    #   temp = array(NA,dim = c(n,s))
    #   S0 = array(n); S0pie = array(n)
    #   for(i in n:1){
    #     S0pie[i] = exp(c(beta%*%Zs[i,]))
    #     if(i == n){ S0[i] = S0pie[i] }else{ S0[i] = S0pie[i] + S0[i+1] } 
    #   }
    #   S1 = array(0,dim = c(n,s)); S1pie = array(0,dim = c(n,s))
    #   for(i in n:1){
    #     S1pie[i,] = S0pie[i]*Zs[i,]
    #     if(i == n){ S1[i,] = S1pie[i,] }else{ S1[i,] = S1pie[i,] + S1[i+1,] }
    #   }
    #   
    #   S2 = array(0,dim = c(n,s,s)); S2pie = array(0,dim = c(n,s,s))
    #   for(j in n:1){
    #     S2pie[j,,] = S0pie[j]*Zs[j,]%*%t(Zs[j,])
    #     if(j == n){S2[j,,] = S2pie[j,,]}else{S2[j,,] = S2pie[j,,] + S2[j+1,,]}
    #   }
    #   ldd = array(0,dim = c(n,s,s)); lddpie = array(0,dim = c(n,s,s))
    #   for(i in 1:n){
    #     lddpie[i,,] = -ds[i]/(S0[i]^2)*(S0[i]*S2[i,,] - S1[i,]%*%t(S1[i,]))
    #     if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
    #   }
    #   
    #   return(ldd[n,,])
    # }
    # SEE = sqrt(diag(solve(-Ubb(bhat))))
    # ,bhat[1:r1]， SEE[1:r1]
    
    return(data.frame(constant_coef = c(beta1), constant_coef_SEE = beta1_sd))
  }
  
  esti = function(h,Zs,ds,Xs){
    
    cri  = 0.001
    nit  = 100 # number of iteration
    conv = array(0,n)
    
    bhat = array(0,dim = c(n,r))
    hmin = max(which(Xs<min(Xs)+h))
    hmax = min(which(Xs>max(Xs)-h)) - 1
    for(t in hmin:hmax){
      brec     = array(0,dim = c(nit,r))
      brec[1,] = rep(.1,r)
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[t]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      
      for(ite in 2:nit){
        
        beta = brec[ite-1,]
        
        G = array(n); Gpie = array(n)
        for(i in n:min_i){
          Gpie[i] = exp(c(beta%*%Zs[i,]))
          if(i == n){ G[i] = Gpie[i] }else{ G[i] = Gpie[i] + G[i+1] }
        }
        Gd = array(0,dim = c(n,r)); Gdpie = array(0,dim = c(n,r))
        for(i in n:min_i){
          Gdpie[i,] = Gpie[i]*Zs[i,]
          if(i == n){ Gd[i,] = Gdpie[i,] }else{ Gd[i,] = Gdpie[i,] + Gd[i+1,] }
        }
        ld = array(0,dim = c(n,r)); ldpie = array(0,dim = c(n,r))
        for(i in non0_i){
          ldpie[i,] = ds[i]*Ker[i]*(Zs[i,] - Gd[i,]/G[i]);
          if(i == 1){ld[i,] = ldpie[i,]}else{ ld[i,] = ld[i-1,] + ldpie[i,] };
        }
        Gmat = array(0,dim = c(n,r,r)); Gmatpie = array(0,dim = c(n,r,r))
        for(j in n:min_i){
          Gmatpie[j,,] = Gpie[j]*Zs[j,]%*%t(Zs[j,])
          if(j == n){Gmat[j,,] = Gmatpie[j,,]}else{Gmat[j,,] = Gmatpie[j,,] + Gmat[j+1,,]}
        }
        ldd = array(0,dim = c(n,r,r)); lddpie = array(0,dim = c(n,r,r))
        for(i in non0_i){
          lddpie[i,,] = -ds[i]*Ker[i]/(G[i]^2)*(G[i]*Gmat[i,,] - Gd[i,]%*%t(Gd[i,]))
          if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
        }
        
        if(abs(det(ldd[max_i,,]))>1e-8){
          temp = c(solve(ldd[max_i,,])%*%ld[max_i,])
          brec[ite,] = beta - temp
        }else{
          brec[ite,] = beta + 0.01
        }
        
        if(sum(abs(brec[ite,]-brec[ite-1,])) < cri) break
        if(max(abs(brec[ite,])) > 10){ conv[t] = 1; brec[ite,] = NA; break }
      }
      bhat[t,] = brec[ite,]
    }
    for(i in 1:(hmin-1)){ bhat[i,] = bhat[hmin,] }
    for(i in (hmax+1):n){ bhat[i,] = bhat[hmax,] }
    return(list(bhat = bhat, conv = conv))
    
  }
  
  esti_sd = function(h,Zs,ds,Xs){
    
    cri   = 0.001
    nit   = 100
    
    hmin = max(which(Xs<min(Xs)+h))
    hmax = min(which(Xs>max(Xs)-h)) - 1
    bsd  = array(0,dim = c(n,r))
    for(t in hmin:hmax){
      brec     = array(0,dim = c(nit,r))
      brec[1,] = rep(.1,r)
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[t]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      
      betaM = array(0,dim = c(M,r))
      conv2 = array(0,M)
      for(k in 1:M){
        Gi = rexp(n)
        for(ite in 2:nit){
          
          beta = brec[ite-1,]
          
          G = array(n); Gpie = array(n)
          for(i in n:min_i){
            Gpie[i] = exp(c(beta%*%Zs[i,]))
            if(i == n){ G[i] = Gpie[i] }else{ G[i] = Gpie[i] + G[i+1] } # View(cbind(G,Gpie))
          }
          Gd = array(0,dim = c(n,r)); Gdpie = array(0,dim = c(n,r))
          for(i in n:min_i){
            Gdpie[i,] = Gpie[i]*Zs[i,]
            if(i == n){ Gd[i,] = Gdpie[i,] }else{ Gd[i,] = Gdpie[i,] + Gd[i+1,] }; # View(cbind(Gd,Gdpie))
          }
          ld = array(0,dim = c(n,r)); ldpie = array(0,dim = c(n,r))
          for(i in non0_i){
            ldpie[i,] =  Gi[i]*ds[i]*Ker[i]*(Zs[i,] - Gd[i,]/G[i]); # View(cbind(ldpie,ld))
            if(i == 1){ld[i,] = ldpie[i,]}else{ ld[i,] = ld[i-1,] + ldpie[i,] }; # View(cbind(ldpie,ds,Ker))
          }
          Gmat = array(0,dim = c(n,r,r)); Gmatpie = array(0,dim = c(n,r,r))
          for(j in n:min_i){
            Gmatpie[j,,] = Gpie[j]*Zs[j,]%*%t(Zs[j,])
            if(j == n){Gmat[j,,] = Gmatpie[j,,]}else{Gmat[j,,] = Gmatpie[j,,] + Gmat[j+1,,]}
          }
          ldd = array(0,dim = c(n,r,r)); lddpie = array(0,dim = c(n,r,r))
          for(i in non0_i){
            lddpie[i,,] = -Gi[i]*ds[i]*Ker[i]/(G[i]^2)*(G[i]*Gmat[i,,] - Gd[i,]%*%t(Gd[i,]))
            if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
          }
          
          if(abs(det(ldd[max_i,,]))>1e-8){
            temp = c(solve(ldd[max_i,,])%*%ld[max_i,])
            brec[ite,] = beta - temp
          }else{
            brec[ite,] = beta + 0.01
          }
          
          if(sum(abs(brec[ite,]-brec[ite-1,])) < cri) break
          if(max(abs(brec[ite,]))>10){ conv2[k] = 1; break }
        }
        betaM[k,] = beta
      }
      if(sum(conv2) > 0 & sum(conv2) < (M-1) ){
        bsd[t,] = apply(betaM[-which(conv2==1),], 2, sd)
      }else if(sum(conv2) >= (M-1)){
        bsd[t,] = NA
      }else{
        bsd[t,] = apply(betaM, 2, sd)
      }
    }
    
    for(i in 1:(hmin-1)){ bsd[i,] = bsd[hmin,] }
    for(i in (hmax+1):n){ bsd[i,] = bsd[hmax,] }
    return(bsd)
  }
  
  # Rename of covar and obstime
  covname = c(colnames(cva_cons), colnames(cva_time))
  Z1  = cva_cons
  Z2  = cva_time
  if(sum(class(Z1) != c('matrix')) == 0) stop("Please transform cva_cons into matrix with specified variable name", call. = FALSE)
  if(sum(class(Z2) != c('matrix')) == 0) stop("Please transform cva_time into matrix with specified variable name", call. = FALSE)
  if(sum( nchar(colnames(Z1)) == 0 )){   stop("Please transform cva_cons into matrix with specified variable name", call. = FALSE)  }
  if(sum( nchar(colnames(Z2)) == 0 )){   stop("Please transform cva_time into matrix with specified variable name", call. = FALSE)  }
  X   = obstime
  M   = resamp
  
  # Some constants
  Z   = cbind(Z1,Z2)
  covname = colnames(Z)
  
  n   = nrow(Z)
  r1  = ncol(Z1)
  r2  = length(ncol(Z2))
  r   = r1 + r2
  n1  = length(delta)
  n2  = length(X)
  ind = sum(is.na(Z)) + sum(is.na(delta)) + sum(is.na(X))
  
  # sort order for X
  XZ  = as.matrix(cbind(X,Z,delta))
  sor = XZ[order(XZ[,1]),]
  Xs  = sor[,1]
  Zs  = sor[,(1:r)+1]
  ds  = sor[,r+2]
  
  if(n != n1 || n != n2 || n1 != n2){
    return("Incorrect length of covariates, censoring indicator or observed time")
  }else if(ind > 0){
    return("data contains NA's")
  }else{
    h = band()
    #print("Commencing estimation of temporal coefficient")
    tem1 = esti(h,Zs,ds,Xs)
    bhat = tem1$bhat
    convind = tem1$conv
    
    if(SE){
      #print("Commencing estimation of standard error")
      bsd = esti_sd(h,Zs,ds,Xs)
    }else{
      bsd = array(NA, dim = c(nrow(bhat), ncol(bhat)) )
    }
    
    #print("Commencing estimation of constant effect")
    conres  = cons(bhat,h,Zs,ds,Xs)
    if(length(convind) > 0){
      bhat[convind,] = NA
      bsd[convind,]  = NA
    }
    data = cbind(Zs, Xs, ds)
    colnames(data) = c(covname,'obstime','delta')
    colnames(bhat) = covname
    colnames(bsd)  = paste(covname,"_SEE", sep = "")
    
    if(r2 > 0){
      ress = list(temporal_coef = bhat[,(1:r2)+r1], temporal_coef_SEE = bsd[,(1:r2)+r1],
                  constant_coef = conres$constant_coef, constant_coef_SEE = conres$constant_coef_SEE,
                  bandwidth = h, data = data,
                  "Time points (not converged)" = paste('obstime = ', Xs[convind], sep = ""))
    }else{
      ress = list(temporal_coef = NULL, temporal_coef_SEE = NULL,
                  constant_coef = conres$constant_coef, constant_coef_SEE = conres$constant_coef_SEE,
                  bandwidth = h, data = data,
                  "Time points (not converged)" = paste('obstime = ', Xs[convind], sep = ""))
    }
    return(ress)
  }
}



npcox = function(cva, delta, obstime, SE = FALSE, bandwidth = FALSE, resamp = 100){
  
  Z   = cva
  if(sum(class(Z) != c('matrix')) == 0) stop("Please transform covariate into matrix with specified variable name", call. = FALSE)
  if(sum( nchar(colnames(Z)) == 0 )){   stop("Please transform covariate into matrix with specified variable name", call. = FALSE) }
  covname = colnames(cva)
  
  # PH part: β(t) estimation
  Kh   = function(h,x){ifelse(abs(x/h)<1, 0.75*(1-(x/h)^2),0)}
  
  PE   = function(bhat,Zs,ds,Xs){
    # prediction error calculation, refer to Lutian(2005)
    bz  = array(0, dim = c(n,n))
    for(i in 1:n){
      for(j in 1:n){
        bz[i,j] = sum(bhat[i,]*Zs[j,])
      }
    }
    pie = array(0, n)
    for(l in 1:n){
      pie[l] = -( bz[l,l] - log( sum(exp(bz[l,l:n])) ) )
    }
    return(sum(pie))
  }
  
  band = function(){
    
    if(bandwidth){
      h = bandwidth
    }else{
      diff  = Xs
      for(i in 2:n){ diff[i] = Xs[i] - Xs[i-1] }
      diff  = sort(diff,decreasing = T)
      hmi   = as.numeric(quantile(Xs,0.05))
      sep   = (as.numeric(quantile(Xs,0.2)) - hmi)/50
      ha    = seq(hmi, as.numeric(quantile(Xs,0.2)), by = sep)
      PErec = array(0, length(ha))
      for(l in 1:length(ha)){
        h    = ha[l]
        bhat = esti(h,Zs,ds,Xs)
        PErec[l] = PE(bhat$bhat,Zs,ds,Xs)
        # print(paste(Sys.time(), '  ', l ))
      }
      # Bandwidth selection
      h = ha[which(PErec == min(PErec[which(!is.na(PErec))]))]
      # cbind(ha, PErec)
      print(paste("Note: The selected bandwidth is ", h, sep = ""))
    }
    return(h)
  }
  
  esti = function(h,Zs,ds,Xs){
    
    cri  = 0.001
    nit  = 100 # number of iteration
    conv = array(0,n)
    
    bhat = array(0,dim = c(n,r))
    hmin = max(which(Xs<min(Xs)+h))
    hmax = min(which(Xs>max(Xs)-h)) - 1
    for(t in hmin:hmax){
      brec     = array(0,dim = c(nit,r))
      brec[1,] = rep(.1,r)
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[t]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      
      for(ite in 2:nit){
        
        beta = brec[ite-1,]
        
        G = array(n); Gpie = array(n)
        for(i in n:min_i){
          Gpie[i] = exp(c(beta%*%Zs[i,]))
          if(i == n){ G[i] = Gpie[i] }else{ G[i] = Gpie[i] + G[i+1] }
        }
        Gd = array(0,dim = c(n,r)); Gdpie = array(0,dim = c(n,r))
        for(i in n:min_i){
          Gdpie[i,] = Gpie[i]*Zs[i,]
          if(i == n){ Gd[i,] = Gdpie[i,] }else{ Gd[i,] = Gdpie[i,] + Gd[i+1,] }
        }
        ld = array(0,dim = c(n,r)); ldpie = array(0,dim = c(n,r))
        for(i in non0_i){
          ldpie[i,] = ds[i]*Ker[i]*(Zs[i,] - Gd[i,]/G[i]);
          if(i == 1){ld[i,] = ldpie[i,]}else{ ld[i,] = ld[i-1,] + ldpie[i,] };
        }
        Gmat = array(0,dim = c(n,r,r)); Gmatpie = array(0,dim = c(n,r,r))
        for(j in n:min_i){
          Gmatpie[j,,] = Gpie[j]*Zs[j,]%*%t(Zs[j,])
          if(j == n){Gmat[j,,] = Gmatpie[j,,]}else{Gmat[j,,] = Gmatpie[j,,] + Gmat[j+1,,]}
        }
        ldd = array(0,dim = c(n,r,r)); lddpie = array(0,dim = c(n,r,r))
        for(i in non0_i){
          lddpie[i,,] = -ds[i]*Ker[i]/(G[i]^2)*(G[i]*Gmat[i,,] - Gd[i,]%*%t(Gd[i,]))
          if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
        }
        
        if(abs(det(ldd[max_i,,]))>1e-8){
          temp = c(solve(ldd[max_i,,])%*%ld[max_i,])
          brec[ite,] = beta - temp
        }else{
          brec[ite,] = beta + 0.01
        }
        
        if(sum(abs(brec[ite,]-brec[ite-1,])) < cri) break
        if(max(abs(brec[ite,])) > 10){ conv[t] = 1; brec[ite,] = NA; break }
      }
      bhat[t,] = brec[ite,]
    }
    for(i in 1:(hmin-1)){ bhat[i,] = bhat[hmin,] }
    for(i in (hmax+1):n){ bhat[i,] = bhat[hmax,] }
    return(list(bhat = bhat, conv = conv))
    
  }
  
  esti_sd = function(h,Zs,ds,Xs){
    
    cri   = 0.001
    nit   = 100
    
    hmin = max(which(Xs<min(Xs)+h))
    hmax = min(which(Xs>max(Xs)-h)) - 1
    bsd  = array(0,dim = c(n,r))
    for(t in hmin:hmax){
      brec     = array(0,dim = c(nit,r))
      brec[1,] = rep(.1,r)
      Ker   = array(n); for(i in 1:n){ Ker[i] = Kh(h,Xs[i]-Xs[t]) }
      non0_i= which(Ker>0)
      min_i = min(non0_i)
      max_i = max(non0_i)
      
      betaM = array(0,dim = c(M,r))
      conv2 = array(0,M)
      for(k in 1:M){
        Gi = rexp(n)
        for(ite in 2:nit){
          
          beta = brec[ite-1,]
          
          G = array(n); Gpie = array(n)
          for(i in n:min_i){
            Gpie[i] = exp(c(beta%*%Zs[i,]))
            if(i == n){ G[i] = Gpie[i] }else{ G[i] = Gpie[i] + G[i+1] } # View(cbind(G,Gpie))
          }
          Gd = array(0,dim = c(n,r)); Gdpie = array(0,dim = c(n,r))
          for(i in n:min_i){
            Gdpie[i,] = Gpie[i]*Zs[i,]
            if(i == n){ Gd[i,] = Gdpie[i,] }else{ Gd[i,] = Gdpie[i,] + Gd[i+1,] }; # View(cbind(Gd,Gdpie))
          }
          ld = array(0,dim = c(n,r)); ldpie = array(0,dim = c(n,r))
          for(i in non0_i){
            ldpie[i,] =  Gi[i]*ds[i]*Ker[i]*(Zs[i,] - Gd[i,]/G[i]); # View(cbind(ldpie,ld))
            if(i == 1){ld[i,] = ldpie[i,]}else{ ld[i,] = ld[i-1,] + ldpie[i,] }; # View(cbind(ldpie,ds,Ker))
          }
          Gmat = array(0,dim = c(n,r,r)); Gmatpie = array(0,dim = c(n,r,r))
          for(j in n:min_i){
            Gmatpie[j,,] = Gpie[j]*Zs[j,]%*%t(Zs[j,])
            if(j == n){Gmat[j,,] = Gmatpie[j,,]}else{Gmat[j,,] = Gmatpie[j,,] + Gmat[j+1,,]}
          }
          ldd = array(0,dim = c(n,r,r)); lddpie = array(0,dim = c(n,r,r))
          for(i in non0_i){
            lddpie[i,,] = -Gi[i]*ds[i]*Ker[i]/(G[i]^2)*(G[i]*Gmat[i,,] - Gd[i,]%*%t(Gd[i,]))
            if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
          }
          
          if(abs(det(ldd[max_i,,]))>1e-8){
            temp = c(solve(ldd[max_i,,])%*%ld[max_i,])
            brec[ite,] = beta - temp
          }else{
            brec[ite,] = beta + 0.01
          }
          
          if(sum(abs(brec[ite,]-brec[ite-1,])) < cri) break
          if(max(abs(brec[ite,]))>10){ conv2[k] = 1; break }
        }
        betaM[k,] = beta
      }
      if(sum(conv2) > 0 & sum(conv2) < (M-1) ){
        bsd[t,] = apply(betaM[-which(conv2==1),], 2, sd)
      }else if(sum(conv2) >= (M-1)){
        bsd[t,] = NA
      }else{
        bsd[t,] = apply(betaM, 2, sd)
      }
    }
    
    for(i in 1:(hmin-1)){ bsd[i,] = bsd[hmin,] }
    for(i in (hmax+1):n){ bsd[i,] = bsd[hmax,] }
    return(bsd)
  }
  
  # Rename of covar and obstime
  X   = obstime
  M   = resamp # if(!resamp){M = 100}else{}
  
  # Some constants
  n   = nrow(Z)
  r   = ncol(Z)
  n1  = length(delta)
  n2  = length(X)
  ind = sum(is.na(Z)) + sum(is.na(delta)) + sum(is.na(X))
  
  # sort order for X
  XZ  = as.matrix(cbind(X,Z,delta))
  sor = XZ[order(XZ[,1]),]
  Xs  = as.numeric(sor[,1])
  Zs  = as.matrix(sor[,(1:r)+1])
  ds  = as.numeric(sor[,r+2])
  
  if(n != n1 || n != n2 || n1 != n2){
    return("Incorrect length of covariates, censoring indicator or observed time")
  }else if(ind > 0){
    return("data contains NA's")
  }else{
    h    = band()
    if(length(h) == 0){
      stop('No available bandwidth can be provided, please specify a bandwidth')
    }else{
      tem1 = esti(h,Zs,ds,Xs)
      bhat = tem1$bhat
      conv = tem1$conv
      
      if(SE){
        #print("Commencing estimation of standard error")
        bsd = esti_sd(h,Zs,ds,Xs)
      }else{
        bsd  = array(NA, dim = c(nrow(bhat), ncol(bhat)) )
      }
      
      ind     = conv
      convind = which(ind == 1)
      if(length(convind) > 0){
        bhat[convind,] = NA
        bsd[convind,]  = NA
      }
      data = cbind(Zs, Xs, ds)
      colnames(data) = c(covname,'obstime','delta')
      colnames(bhat) = covname
      colnames(bsd)  = paste(covname,"_SEE", sep = "")
      
      return(list(temporal_coef = bhat, temporal_coef_SEE = bsd, bandwidth = h, data = data,
                  inconverged = convind,
                  "Time points (not converged)" = paste('obstime = ', Xs[convind], sep = "")))
    }
  }
}



library(nleqslv)
fcox = function(Z,status,time){
  X     = time
  delta = status
  s     = ncol(Z)
  n     = nrow(Z)
  
  # Reorder of all the info
  XZ  = cbind(X,Z,delta)
  sor = XZ[order(XZ[,1]),] 
  Xs  = sor[,1]
  Zs  = sor[,(1:s)+1]
  ds  = sor[,s+2]
  
  Ub = function(beta){
    temp = array(NA,dim = c(n,s))
    S0 = array(n); S0pie = array(n)
    for(i in n:1){
      S0pie[i] = exp(c(beta%*%Zs[i,]))
      if(i == n){ S0[i] = S0pie[i] }else{ S0[i] = S0pie[i] + S0[i+1] } 
    }
    S1 = array(0,dim = c(n,s)); S1pie = array(0,dim = c(n,s))
    for(i in n:1){
      S1pie[i,] = S0pie[i]*Zs[i,]
      if(i == n){ S1[i,] = S1pie[i,] }else{ S1[i,] = S1pie[i,] + S1[i+1,] }
    }
    for(i in 1:n){
      temp[i,] = ds[i]*(Zs[i,] - S1[i,]/S0[i])
    }
    return(colSums(temp))
  }
  bhat = nleqslv(rep(1,s), Ub)$x
  
  Ubb = function(beta){
    temp = array(NA,dim = c(n,s))
    S0 = array(n); S0pie = array(n)
    for(i in n:1){
      S0pie[i] = exp(c(beta%*%Zs[i,]))
      if(i == n){ S0[i] = S0pie[i] }else{ S0[i] = S0pie[i] + S0[i+1] } 
    }
    S1 = array(0,dim = c(n,s)); S1pie = array(0,dim = c(n,s))
    for(i in n:1){
      S1pie[i,] = S0pie[i]*Zs[i,]
      if(i == n){ S1[i,] = S1pie[i,] }else{ S1[i,] = S1pie[i,] + S1[i+1,] }
    }
    
    S2 = array(0,dim = c(n,s,s)); S2pie = array(0,dim = c(n,s,s))
    for(j in n:1){
      S2pie[j,,] = S0pie[j]*Zs[j,]%*%t(Zs[j,])
      if(j == n){S2[j,,] = S2pie[j,,]}else{S2[j,,] = S2pie[j,,] + S2[j+1,,]}
    }
    ldd = array(0,dim = c(n,s,s)); lddpie = array(0,dim = c(n,s,s))
    for(i in 1:n){
      lddpie[i,,] = -ds[i]/(S0[i]^2)*(S0[i]*S2[i,,] - S1[i,]%*%t(S1[i,]))
      if(i == 1){ldd[i,,] = lddpie[i,,]}else{ldd[i,,] = ldd[i-1,,] + lddpie[i,,]}
    }
    
    return(ldd[n,,])
  }
  SEE = sqrt(diag(solve(-Ubb(bhat))))
  
  res = data.frame(beta = bhat, SE = SEE, TS = bhat/SEE)
  rownames(res) = colnames(Z)
  return(res)
}
