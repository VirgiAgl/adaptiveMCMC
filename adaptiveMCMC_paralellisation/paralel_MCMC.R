library(parallel)

start = Sys.time()

d = 8
n_iter = 200     #  Total number of iterations
x_1 = rep(5,d)      #  Vector of inital values
t_adapt = 50      #  When to start adapting


times = 14

mclapply(1:times,function(t){
  if (t==1) {
    adapt = 'None'
    X_MH = mcmc(target = pi_norm_corr,
                n_iter = n_iter,
                x_1 = x_1,
                adapt=adapt)
    save(X_MH,file='R/data/X_MH.rda')
  }

  if (t==2) {
    adapt = 'AM'
    cov_estimator = 'Sample covariance'
    X_AM = mcmc(target = pi_norm_corr,
                n_iter = n_iter,
                x_1 = x_1,
                adapt=adapt,
                t_adapt = t_adapt,
                cov_estimator=cov_estimator
    )
    save(X_AM,file='R/data/X_AM.rda')
  }

  if (t==3) {
    adapt = 'AM2'
    cov_estimator = 'Sample covariance'
    X_AM2 = mcmc(target = pi_norm_corr,
                n_iter = n_iter,
                x_1 = x_1,
                adapt=adapt,
                t_adapt = t_adapt,
                cov_estimator=cov_estimator
    )
    save(X_AM2,file='R/data/X_AM2.rda')
  }

  if (t==4) {
    adapt = 'AM'
    cov_estimator = 'Shrinkage estimator'
    X_AM_sh = mcmc(target = pi_norm_corr,
                n_iter = n_iter,
                x_1 = x_1,
                adapt=adapt,
                t_adapt = t_adapt,
                cov_estimator=cov_estimator
    )
    save(X_AM_sh,file='R/data/X_AM_sh.rda')
  }

  if (t==5) {
    adapt = 'AM'
    cov_estimator = 'Thresholding estimator'
    X_AM_th = mcmc(target = pi_norm_corr,
                   n_iter = n_iter,
                   x_1 = x_1,
                   adapt=adapt,
                   t_adapt = t_adapt,
                   cov_estimator=cov_estimator
    )
    save(X_AM_th,file='R/data/X_AM_th.rda')
  }

  if (t==6) {
    adapt = 'AM2'
    cov_estimator = 'Shrinkage estimator'
    X_AM2_sh = mcmc(target = pi_norm_corr,
                   n_iter = n_iter,
                   x_1 = x_1,
                   adapt=adapt,
                   t_adapt = t_adapt,
                   cov_estimator=cov_estimator
    )
    save(X_AM2_sh,file='R/data/X_AM2_sh.rda')
  }

  if (t==7) {
    adapt = 'AM2'
    cov_estimator = 'Thresholding estimator'
    X_AM2_th = mcmc(target = pi_norm_corr,
                    n_iter = n_iter,
                    x_1 = x_1,
                    adapt=adapt,
                    t_adapt = t_adapt,
                    cov_estimator=cov_estimator
    )
    save(X_AM2_th,file='R/data/X_AM2_th.rda')
  }

  if (t==8) {
    adapt = 'None'
    X_MH_banana = mcmc(target = pi_banana8,
                n_iter = n_iter,
                x_1 = x_1,
                adapt=adapt)
    save(X_MH_banana,file='R/data/X_MH_banana.rda')
  }

  if (t==9) {
    adapt = 'AM'
    cov_estimator = 'Sample covariance'
    X_AM_banana = mcmc(target = pi_banana8,
                n_iter = n_iter,
                x_1 = x_1,
                adapt=adapt,
                t_adapt = t_adapt,
                cov_estimator=cov_estimator
    )
    save(X_AM_banana,file='R/data/X_AM_banana.rda')
  }

  if (t==10) {
    adapt = 'AM2'
    cov_estimator = 'Sample covariance'
    X_AM2_banana = mcmc(target = pi_banana8,
                 n_iter = n_iter,
                 x_1 = x_1,
                 adapt=adapt,
                 t_adapt = t_adapt,
                 cov_estimator=cov_estimator
    )
    save(X_AM2_banana,file='R/data/X_AM2_banana.rda')
  }

  if (t==11) {
    adapt = 'AM'
    cov_estimator = 'Shrinkage estimator'
    X_AM_sh_banana = mcmc(target = pi_banana8,
                   n_iter = n_iter,
                   x_1 = x_1,
                   adapt=adapt,
                   t_adapt = t_adapt,
                   cov_estimator=cov_estimator
    )
    save(X_AM_sh_banana,file='R/data/X_AM_sh_banana.rda')
  }

  if (t==12) {
    adapt = 'AM'
    cov_estimator = 'Thresholding estimator'
    X_AM_th_banana = mcmc(target = pi_banana8,
                   n_iter = n_iter,
                   x_1 = x_1,
                   adapt=adapt,
                   t_adapt = t_adapt,
                   cov_estimator=cov_estimator
    )
    save(X_AM_th_banana,file='R/data/X_AM_th_banana.rda')
  }

  if (t==13) {
    adapt = 'AM2'
    cov_estimator = 'Shrinkage estimator'
    X_AM2_sh_banana = mcmc(target = pi_banana8,
                    n_iter = n_iter,
                    x_1 = x_1,
                    adapt=adapt,
                    t_adapt = t_adapt,
                    cov_estimator=cov_estimator
    )
    save(X_AM2_sh_banana,file='R/data/X_AM2_sh_banana.rda')
  }

  if (t==14) {
    adapt = 'AM2'
    cov_estimator = 'Thresholding estimator'
    X_AM2_th_banana = mcmc(target = pi_banana8,
                    n_iter = n_iter,
                    x_1 = x_1,
                    adapt=adapt,
                    t_adapt = t_adapt,
                    cov_estimator=cov_estimator
    )
    save(X_AM2_th_banana,file='R/data/X_AM2_th_banana.rda')
  }

}
,mc.cores = 14)

now = Sys.time()
print(now - start)

