######################################################
# -- Pure Data Generation from Matern correlation -- #
######################################################

# - range.par is the range parameter 
# - nu.par is the smoothness, which is set at 1 - #

# Location generation #
# - The following function will create a grid of locations - #
# - It will return an n x 2 matrix - #

site.fn = function(domain.size, len)
{
  x.grid = seq(0, domain.size, length.out = len)
  y.grid = seq(0, domain.size, length.out = len)
  site.mat = make.surface.grid(list(x.grid,y.grid))
  
  return(site.mat)
}

# Data generation # 
# - The following function will create an M x n matrix where - #
# - every row of this big data matrix corresponds to one - #
# - simulated data set. We have M such data sets generated on the - #
# - same grid created above - #

dat.gen = function(sim, site, range.par, nu.par, sigma2 = 1)
{
  nn = nrow(site)
  dat.mat = matrix(0, nrow = sim, ncol = nn)
  
  for (j in 1:sim)
  {
    dist.mat = rdist(site, site)
    Sigma = sigma2*Matern(dist.mat, range = range.par, nu = nu.par)  
    Sigma.c = chol(Sigma)   # Cholesky factor in R
    
    # -- Simulated data -- #    
    sim.dat = t(Sigma.c) %*% rnorm(nn)    
    dat.mat[j, ] = sim.dat
  }
    return(dat.mat)  
}

# THE FOLLOWING EXAMPLE CAN BE RUN TO UNDERSTAND THE OUTPUTS #
# Example: 
# loc.mat = site.fn(domain.size = 10, len=10)
# dim(loc.mat)
# plot(loc.mat)
# dat.tmp = dat.gen(sim = 2, site = loc.mat, range.par = 1, nu.par = 1, sigma2 = 1)
# dim(dat.tmp)

# set.panel(1,2)
# quilt.plot(loc.mat, dat.tmp[1,])
# quilt.plot(loc.mat, dat.tmp[2,])





