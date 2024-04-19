# Huber psi function #
# - The function evaluates the first derivative 
#of Huber's loss function - #

Huber.psi = function(X, C )
{
  psi = ifelse(abs(X) <= C, X, C*sign(X))
  return(psi)
}

#Huber rho function
#Huber's loss function before taking 
# the derivative #
Huber.rho = function(X, C)
{
  rho = ifelse(abs(X) <= C, 0.5*(X^2),
               C*(abs(X) - (C*0.5)))
  return(rho)
}
##############################################
# -- Defining the robust spatial function -- #
##############################################

# In the following function #  
# - site is the location matrix - #
# - y is the vector of observations possibly with outliers - #
# - Residuals truncated outside of Cval*res.Scale to  +/- this limit. - #
# - res.Scale is the MAD adjusted to be unbiased for the standard deviation - #
# - for a normal sample. - #
# - maxIteration is the maximum number of iterations. - #
# - tol is a convergence criterion. - #


spatial.robust.new = function(site, datVec, Cval,
                              fhat.old,
                              range.par, lambda.par, 
                              nu.par, MLE = TRUE,
                              maxIteration = 200, 
                              tol = 1e-7, 
                              verbose = FALSE)
{
  
  # computing the residuals
  res.vec.old = (datVec-fhat.old)

  #Initial Estimates of the covariance parameters
  aRangeOld = range.par 
  lambdaOld = lambda.par
  
  kk = 0 # iteration index
  conv.flag = FALSE
  
  while(kk <= maxIteration) 
  {
    kk = kk + 1
    
    if(verbose){cat(kk, lambdaOld, aRangeOld, fill = TRUE)}
    
    # Computing the Pseudo-data and fitting a spatial model 
    
    PSy = fhat.old +  Huber.psi(res.vec.old, Cval)
    
    obj = spatialProcess(site, PSy,
                           lambda = lambdaOld, 
                           aRange = aRangeOld,
                           cov.args = list(Covariance = "Matern", 
                                           smoothness = nu.par)) 
    
    #Updating the estimates
    fhat.new = predict(obj)
      
    lambdanew = obj$summary["lambda"]
    aRangenew = obj$summary["aRange"]
    
    # Setting up the convergence criterion
    test = mean(abs(fhat.new - fhat.old)) / mean(abs(fhat.old))
    if (test <= tol) 
    {
      conv.flag <- TRUE
      break
    } 
    fhat.old = fhat.new
    res.vec.old = datVec-fhat.old
  }
  
  # Final fit with pseudo data after convergence
  
  aRangeMLE = aRangenew
  lambdaMLE = lambdanew
  
  obj.final = spatialProcess(site, PSy, aRange = aRangeMLE, 
                             lambda = lambdaMLE,
                             cov.args = list(Covariance = "Matern", 
                              smoothness = nu.par))
  
  class(obj.final) = c(class(obj), "robustGP")
  
  # Add some extra components from robust fitting 
  obj.final$yRaw = datVec
  obj.final$Cval = Cval
  obj.final$convergenceInfo<- data.frame(iter= kk, test=test, conv.flag=conv.flag)
  obj.final$call<- match.call()
  return(obj.final)
}

spatial.robust.IRWLS = function(site, datVec, Cval,
                                fhat.old,
                                range.par, lambda.par,
                                nu.par, MLE = FALSE,
                                maxIteration = 200,
                                tol = 1e-7,
                                verbose = FALSE)
{
  
  # computing the residuals
  res.vec.old = (datVec-fhat.old)
  
  #Initial Estimates of the covariance parameters
  aRangeOld = range.par
  lambdaOld = lambda.par
  
  kk = 0 # iteration index
  conv.flag = FALSE
  
  while(kk <= maxIteration)
  {
    kk = kk + 1
    
    if(verbose){cat(kk, lambdaOld, aRangeOld, fill = TRUE)}
    
    # Computing the Pseudo-data and fitting a spatial model
    
    W<- Huber.rho(res.vec.old, Cval)/(res.vec.old^2)
    
    if(!MLE){
      obj = spatialProcess(site, datVec,
                           lambda = lambdaOld,
                           aRange = aRangeOld,
                           weights= W,
                           cov.args = list(Covariance = "Matern",
                                           smoothness = nu.par))
      lambdaNew = obj$summary["lambda"]
      aRangeNew = obj$summary["aRange"]
    }
    else{
      obj = spatialProcess(site, datVec,
                           cov.params.start=
                             list( lambda = lambdaOld,
                                   aRange = aRangeOld),
                           weights= W,
                           cov.args = list(Covariance = "Matern",
                                           smoothness = nu.par))
      lambdaNew = obj$summary["lambda"]
      aRangeNew = obj$summary["aRange"]
      cat( lambdaNew,aRangeNew, fill=TRUE )
      lambdaOld<- lambdaNew
      aRangeOld<- aRangeNew
    }
    
    #Updating the estimates
    fhat.new = predict(obj)
    
    
    # Setting up the convergence criterion
    test = mean(abs(fhat.new - fhat.old)) / mean(abs(fhat.old))
    if (test <= tol)
    {
      conv.flag <- TRUE
      break
    }
    fhat.old = fhat.new
    res.vec.old = datVec-fhat.old
  }
  
  # Final fit with pseudo data after convergence
  
  aRangeMLE = aRangeNew
  lambdaMLE = lambdaNew
  
  obj.final = spatialProcess(site, datVec, aRange = aRangeMLE,
                             lambda = lambdaMLE,
                             weights= W,
                             cov.args = list(Covariance = "Matern",
                                             smoothness = nu.par))
  
  class(obj.final) = c(class(obj), "robustGP")
  
  # Add some extra components from robust fitting
  obj.final$yRaw = datVec
  obj.final$Cval = Cval
  obj.final$convergenceInfo<- c(iter= kk, test=test, conv.flag=conv.flag)
  obj.final$call<- match.call()
  return(obj.final)
}
