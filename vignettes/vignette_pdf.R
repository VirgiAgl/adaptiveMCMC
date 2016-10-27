## ----initialisations, eval=FALSE-----------------------------------------
#  
#  d = 8                             #  Dimension of the parameter space
#  n_iter = 20000                    #  Total number of iterations
#  x1 = rep(5,d)                     #  Vector of inital values
#  t_adapt = 2000                    #  When to start adapting
#  adapt = "AM"                      #  Algorithm to run.
#  target = pi_norm_corr             #  Target distribution function
#  cov_estimator="Sample covariance" # Coviance matrix estimator

## ----initialisations2, eval=FALSE----------------------------------------
#  
#  X = mcmc(target = target,
#           n_iter = n_iter,
#           x_1 = x1,
#           adapt=adapt,
#           t_adapt = t_adapt
#           )

