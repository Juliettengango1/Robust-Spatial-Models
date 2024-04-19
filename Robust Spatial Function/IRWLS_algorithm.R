RobSpatialProcess = function(){
  # Robust spatial process using M-estimation
  
  # Created on: 10.10.2023
  # Modified on: 02.07.2024
  # Abhijit Mandal, University of Texas at El Paso, USA
  # Joint work with Doug Nychka, Soutir Bandyopadhyay and Juliette Mukangango
  # ____________________________________________________________________________
  library(fields)
  
  # true covariance parameters
  nu = 3/2         #smoothness param.
  sigma2 = 1       #process variance
  tau2 = 0.02      #nugget effect
  theta = 0.4      #aRange param.
  
  # min and max for param = c(theta, sigma2), to be used in optim
  param_min = c(0.01, 0.01); param_max = c(2, 5) 
  
  c = 1.28 # param for Huber function
  
  # site
  X_range = seq(1, 10, length.out = 10)
  X = as.matrix(expand.grid(X_range, X_range)) 
  colnames(X) = c("long", "lat")
  N = nrow(X)                #sample size
  
  
  # covariance matrix
  Dist_site = rdist(X, X)
  covTrue = Matern(Dist_site, range = theta, nu = nu)
  covT = chol(covTrue)
  
  beta = c(1.5, -2)
  # Generating pure data
  z = t(covT) %*% rnorm(N)
  y = X %*% beta + z + sqrt(tau2) * rnorm(N)
  
  # Clasical spatial model
  spmodel  = spatialProcess(x=X, y=y, reltol = 1e-6, #gridN = 10,
                            cov.args = list(Covariance = "Matern",
                                            smoothness = nu)) #,
  # smoothness = nu))
  
  # Outputs from Classical Spatial Process
  cat("\nClassical Spatial Process:\n")
  print(spmodel$summary)        
  
  cat("\nRobust Spatial Process:\n")
  # the last outputs may be used for an initial value, just checking if the true values works
  Rob_spatialProcess(X=X, y=y, z=z, Dist_site=Dist_site, beta=beta, theta=theta,
                     tau2=tau2, sigma2=sigma2, 
                     c=1.28, nu = 3/2, tol=1e-4, maxItr=100,
                     param_min = param_min, param_max = param_max)
  
  
  
}
# ____________________________________________________________________________
#                       Matern Correlation Matrix
# ____________________________________________________________________________
Matern_corr = function(R, theta, nu = 3/2){
  # The special cases of Matern correlation matrix 
  if (nu == 0.5) return(exp(-R/theta))
  if (nu == 3/2) return((1 + sqrt(3)*R/theta) * exp(-sqrt(3)*R/theta))
  if (nu == Inf) return(exp(-0.5 * (R/theta)^2) )
  stop("error in nu")
} 

# ____________________________________________________________________________
#                       Huber function and derivatives
# ____________________________________________________________________________
Huber_rho_with_diff = function(x, c=1.28){
  # Huber rho function with the first two derivatives
  
  is.x.less.c = (abs(x) <= c)
  rho_val = is.x.less.c * (0.5 * x^2) + 
    (1 - is.x.less.c) * (c * abs(x) - 0.5 * c^2)
  
  rho_1st_diff = is.x.less.c * x + 
    (1 - is.x.less.c) * (c * sign(x))
  
  rho_2nd_diff = is.x.less.c + 0
  
  return(list(rho_val=rho_val, rho_1st_diff=rho_1st_diff, rho_2nd_diff=rho_2nd_diff))
}

# ____________________________________________________________________________
#            Pseudo Likelihood for estimating theta and sigma2
# ____________________________________________________________________________
Pseudo_lik = function(X, y, z, Dist_site, beta, theta, tau2, sigma2, 
                      c=1.28, nu = 3/2){
  # Pseudo-likelihood function for estimating theta and tau^2
  
  n = length(y)
  
  epsilon = y - X %*% beta - z
  Huber_outputs = Huber_rho_with_diff(x = epsilon/sqrt(tau2), c=c)
  
  trace_D = sum(Huber_outputs$rho_val)
  D2 = diag(as.vector(Huber_outputs$rho_2nd_diff))
  
  Sigma = sigma2 * Matern_corr(R = Dist_site, theta = theta, nu=nu)
  
  lik_val = (-trace_D - 0.5* sum(z * solve(Sigma, z)) ) - 
    log( sqrt(det(Sigma %*% D2 + tau2 * diag(n)) ) )
}

# ____________________________________________________________________________
#                       Iterative algo for robust estimates
# ____________________________________________________________________________
Rob_spatialProcess = function(X, y, z, Dist_site, beta, theta, tau2, sigma2, 
                              c=1.28, nu = 3/2, tol=1e-4, maxItr=100,
                              param_min = c(0.01, 0.01), param_max = c(2, 5)){
  
  # min and max for param = c(theta, sigma2), to be used in optim
  
  for (r in 1:maxItr){
    
    z_old = z
    param_old = c(theta, sigma2) 
    
    Pseudo_lik_theta = function(x) Pseudo_lik(X=X, y=y, z=z, Dist_site=Dist_site, 
                                              beta=beta, theta=x[1], tau2=tau2, 
                                              sigma2=x[2], c=c, nu=nu)
    
    # Update theta and sigma2
    param = optim(param_old, Pseudo_lik_theta, lower = param_min,
                  upper = param_max, method = "L-BFGS-B")$par
    theta = param[1]; sigma2 = param[2]
    
    # Update tau^2 and beta
    for (i in 1:maxItr){ 
      beta_old = beta
      
      epsilon = y - X %*% beta - z
      tau2 = (mad(epsilon))^2
      eps_by_tau = epsilon/sqrt(tau2)
      Huber_outputs = Huber_rho_with_diff(x = eps_by_tau, c=c)
      W = diag(as.vector(Huber_outputs$rho_1st_diff/eps_by_tau))
      beta = solve(t(X) %*% W %*% X, t(X) %*% W %*% (y - z))
      
      if (sqrt(sum((beta_old - beta)^2) / max(1e-20, sum(beta_old^2))) < tol) break 
      #tau2 may be added
    }
    
    # Update z
    spatial_obj = spatialProcess(x=X, y=y, 
                                 cov.args=list(aRange=theta), lambda= tau2/sigma2,
                                 mKrig.args= list(sigma2=sigma2) ,
                                 smoothness=nu)
    
    z = spatial_obj$fitted.values - predict(spatial_obj, just.fixed=TRUE)
    
    if (sqrt(sum((z_old - z)^2) / max(1e-20, sum(z_old^2))) < tol) break  
    #other parameters may be added
    
    if (r == maxItr) message("Robust estimate did not converge.")
  }
  
  return(list(beta=beta, theta=theta, tau2=tau2, sigma2=sigma2))
}
