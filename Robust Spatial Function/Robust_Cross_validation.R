#function that find the minima and their locations 
#from the fine grid

which.min.image<-function(obj){
  ind.z <- which.max.matrix(-1*obj$z)
  return(list(x = obj$x[ind.z[, 1]], y = obj$y[ind.z[, 2]], 
              z = obj$z[ind.z], ind = ind.z))
}

#The function to compute cross validation

crossV.fn= function(n, k.fold.val, site,
                    Zvec, Cval, grd.mat,
                    nu.par, maxIteration = 200,
                    tol = 1e-4)
{
  
  #Creating K equally sized folds
  set.seed(222)
  indShuffle<- sample(1:n, n, replace = FALSE)
  ind.mat = matrix(indShuffle, nrow = k.fold.val)
  
  #length of the 2-D grid
  ngrd = nrow(grd.mat)
  
  #Initializing the loss matrix
  loss.mat = matrix(NA,ncol = k.fold.val, nrow = ngrd)
  
  #Initializing the convergence matrix
  convInfo<-matrix(NA, ncol = k.fold.val, nrow = ngrd)
  
  for (jj in 1:k.fold.val) {
    testIndex = ind.mat[jj,]
    
    #validation set
    s.test = site[testIndex,]
    Z.test = Zvec[testIndex]
    
    #Training set
    s.train = site[-testIndex,]
    Z.train = Zvec[-testIndex]
    
    fhat.train = rep(median(Z.train), length(Z.train))
    
    for (ll in 1:ngrd) {
      
      #fixed values for the covariance Parameters
      lambdaVal = grd.mat[ll,1]
      rangeVal = grd.mat[ll,2]
      
      #Fitting a robust model
      fit.train = spatial.robust.IRWLS(s.train, Z.train,
                                       Cval, fhat.train,
                                       range.par = rangeVal,
                                       lambda.par = lambdaVal,
                                       nu.par, tol = tol)
      
      #Predicting the robust model on the test sample
      yhat = predict(fit.train, s.test)
      
      #Computing the loss error
      loss.mat[ll,jj] = mean(Huber.rho((Z.test-yhat), Cval))
      
      #Looking the model's convergence
      convInfo[ll,jj]<-fit.train$convergenceInfo[1]

    }
  }
  
  #Fitting a curve on the grid
  look<-Tps(log10(grd.mat), rowMeans(loss.mat), lambda = 0, m = 3)
  
  #Evaluating the curve on a fine set of points
  obj.pred<-predictSurface(look, nx = 80, ny = 80)
  
  minPars<-which.min.image(obj.pred)
 
  OptParam = list(minPars = minPars,
                  loss.mat = loss.mat,
                  parGrid = grd.mat)
  
  return(OptParam)
}
