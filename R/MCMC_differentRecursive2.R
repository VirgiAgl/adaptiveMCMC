library(MASS)
library(Matrix)
library(mvtnorm)
library(CVTuningCov)
library(corpcor)
library(PDSCE)
library(ggplot2)

############################################################################
# Still to implement
# ------------------
# Second algorithm for adaptation
# Burn in
# Normal MH acceptance probability
# Do something with the optimal alpha
# Suboptimality factor
# Implement and test higher dimensional problem
# Speed up code
# Check updating acceptance forumla for AM3 and 4
############################################################################

mcmc = function(target, n_iter = 100, x_1, adapt="None", cov_estimator="Sample covariance", t_adapt = Inf){

  start = Sys.time()

  #Print detail of run to screen
  cat("Running MCMC targeting ",
      as.character(substitute(target)),
      "with " ,
      n_iter,
      "iterations.")
  cat("\nAdaptation algorithm is: ",
      adapt,
      "."
  )
  cat("\nCovariance matrix estimator is: ",
      cov_estimator,
      "."
  )

  #####################################################################
  # Define and infer some constants
  #####################################################################
  optimal_alpha = 0.234  #  The optimal acceptance rate
  d = length(x_1)  #  Infer dimension of the problem

  #####################################################################
  # Set up some structures for storing the markov chain and its properties
  #####################################################################
  X = matrix(NA, n_iter+1, d)  #  Initialise matrix to store the markov chain
  X[1,] = x_1  #  Store the first step of the markov chain
  acceptances = c(rep(NA,n_iter+1))  #  Initialise counter for accepted states
  acceptances[1] = 0  #  Initialise counter for accepted states
  expected_values = matrix(NA,n_iter+1,d)  #  Initialise counter for accepted states
  expected_values[1,] = sapply(x_1,mean)  #  Initialise counter for accepted states

  #####################################################################
  # Some initialisation needed for AM adaptation
  #####################################################################
  if (adapt=="AM") {
    s_d = (2.4)^2/d  #  The scaling parameter
    epsilon = 0.0001  #  Constant (small compared to support of target)
    sigma = diag(1,d,d)  #  Initialise the covariance matrix
    #sigma = nearPD(sigma)$mat
    mean_X = matrix(NA, n_iter+1, d) #  Initialise matrix to store means
    mean_X[1,] = x_1 #  First mean is first step of chain
    s_d = (2.4)^2/d  #  The scaling parameter
    beta = 0.05
    }
  if (adapt=="None") {
    sigma = diag(1,d,d) #  No correlation in proposal
  }

  if (adapt == "AM2" ) {
    s_d = (2.4)^2/d  #  The scaling parameter
    beta = 0.05
    mean_X = matrix(NA, n_iter+1, d) #  Initialise matrix to store means
    mean_X[1,] = x_1 #  First mean is first step of chain
    sigma = diag(1,d,d)
  }

  #####################################################################
  # Create the chain
  #####################################################################
  for (i in 2:(n_iter+1)) {

    if(i %% 5000 == 0){
      cat("\n","The current timestep count is ", i)
    }

    if (adapt=="AM") {
      Y = q_norm(mu=X[i-1,], sigma=sigma) #proposal distribution sampled at the current point
      numerator = target(Y)
      denom = target(X[i-1,])
    }

    if (adapt=="AM2"){
      if (i <= 2*d) {
        Y = q_norm(mu=X[i-1,], sigma=sigma) #proposal distribution sampled at the current point
      }

      if (i == 2*d+1) {
        if (cov_estimator=="Sample covariance"){
          X_until_now = X[1:(i-1),] #the rows of X that we have filled so far
          sigma = s_d*cov(X_until_now)
          sigma = nearPD(sigma)$mat
          Y = (beta*q_norm(mu=X[i-1,], sigma=diag(d))) + ((1-beta)*q_norm(mu=X[i-1,],sigma=sigma))
        }
        if (cov_estimator=="Shrinkage estimator"){
          X_until_now = X[1:(i-1),]
          covariance=cor.shrink(X_until_now)
          lambda = attr(covariance,"lambda")
          sigma=s_d*covariance
          sigma=sigma[(1:d),(1:d)]
          sigma = nearPD(sigma)$mat
          Y = (beta*q_norm(mu=X[i-1,], sigma=diag(d))) + ((1-beta)*q_norm(mu=X[i-1,],sigma=sigma))
        }
        if (cov_estimator=="Thresholding estimator"){
          X_until_now = X[1:(i-1),]
          covariance=band.chol(X_until_now,k=4)
          sigma=s_d*covariance
          sigma = nearPD(sigma)$mat
          Y = (beta*q_norm(mu=X[i-1,], sigma=diag(d))) + ((1-beta)*q_norm(mu=X[i-1,],sigma=sigma))
        }
      }

      if (i > 2*d+1) {
        if (cov_estimator=="Sample covariance") {
        sigma = ((i-2)/(i-1))*sigma + (s_d/(i-1))*(((i-1)*mean_X[i-2,]%*%t(mean_X[i-2,]))-(i*mean_X[i-1,]%*%t(mean_X[i-1,]))+mean_X[i-1,]%*%t(mean_X[i-1,]))
        sigma = nearPD(sigma)$mat
        Y = (beta * q_norm(mu=X[i-1,], sigma=diag(1,d,d))) + ((1-beta) * q_norm(mu=X[i-1,], sigma=sigma))
        }
        if (cov_estimator=="Shrinkage estimator") {
          sigma = ((i-2)/(i-1))*sigma + (s_d/(i-1)*(1-lambda))*(((i-1)*mean_X[i-2,]%*%t(mean_X[i-2,]))-(i*mean_X[i-1,]%*%t(mean_X[i-1,]))+mean_X[i-1,]%*%t(mean_X[i-1,]))+((lambda*s_d*diag(d))/(i-1))
          sigma = nearPD(sigma)$mat
          Y = (beta * q_norm(mu=X[i-1,], sigma=diag(1,d,d))) + ((1-beta) * q_norm(mu=X[i-1,], sigma=sigma))
        }
        if (cov_estimator=="Thresholding estimator") {
          X_until_now = X[1:(i-1),]
          covariance=band.chol(X_until_now,k=4)
          sigma=s_d*covariance
          sigma = nearPD(sigma)$mat
          Y = (beta * q_norm(mu=X[i-1,], sigma=diag(1,d,d))) + ((1-beta) * q_norm(mu=X[i-1,], sigma=sigma))
        }

      }
      numerator = target(Y)
      denom = target(X[i-1,])
    }


    if (adapt=="None") {
      Y = q_norm(mu=X[i-1,], sigma=diag(d)) #proposal distribution sampled at the current point
      numerator = target(Y)*q_norm(X[i-1,],sigma=diag(d)) #*************NEED TO CHANGE THIS************
      denom = target(X[i-1,])*q_norm(Y,sigma=diag(d))
    }

    alpha = min(1,numerator/denom)

    if(runif(1) < alpha) {
      X[i,] = Y
      acceptances[i] = acceptances[i-1] + 1
      } else{
      X[i,] = X[i-1,]
      acceptances[i] = acceptances[i-1]
    }

    if (adapt=="AM") {
      mean_X[i,] = ((mean_X[i-1,] * (i-1)) + X[i,])/i

      if (i == t_adapt) { #the first adaptation step
         #the rows of X that we have filled so far
        X_until_now = X[1:i,]
        if (cov_estimator=="Sample covariance") {
          sigma = s_d * ( cov(X_until_now) + epsilon * diag(d))
        }
        if (cov_estimator=="Shrinkage estimator") {
          covariance=cor.shrink(X_until_now)
          lambda = attr(covariance,"lambda")
          sigma=s_d*covariance
          sigma=sigma[(1:d),(1:d)]
          sigma = nearPD(sigma)$mat
          }
        if (cov_estimator=="Thresholding estimator") {
          covariance=band.chol(X_until_now,k=4)
          sigma=s_d*covariance
          sigma = nearPD(sigma)$mat
          }
      }
      if (i > t_adapt) { #the later adaptation step
        X_until_now = X[1:i,] #the rows of X that we have filled so far

        if (cov_estimator=="Sample covariance") {
          #recursive sigma
          sigma = ((i-2)/(i-1))*sigma + (s_d/(i-1))*(((i-1)*mean_X[i-2,]%*%t(mean_X[i-2,]))-(i*mean_X[i-1,]%*%t(mean_X[i-1,]))+mean_X[i-1,]%*%t(mean_X[i-1,])+epsilon*diag(d))
          sigma = nearPD(sigma)$mat
        }
        if (cov_estimator=="Shrinkage estimator") {
          sigma = ((i-2)/(i-1))*sigma + (s_d/(i-1)*(1-lambda))*(((i-1)*mean_X[i-2,]%*%t(mean_X[i-2,]))-(i*mean_X[i-1,]%*%t(mean_X[i-1,]))+mean_X[i-1,]%*%t(mean_X[i-1,]))+((lambda*s_d*diag(d))/(i-1))
          sigma = nearPD(sigma)$mat
        }
        if (cov_estimator=="Thresholding estimator") {
          covariance=band.chol(X_until_now,k=4)
          sigma=s_d*covariance
          sigma = nearPD(sigma)$mat
        }
      }
    }

    if (adapt=="AM2") {
      mean_X[i,] = ((mean_X[i-1,] * (i-1)) + X[i,])/i
    }
    X_until_now = X[1:i,]
    expected_values[i,] = apply(X_until_now,2,mean)
    alpha_values<-acceptances/n_iter

  }  #  End of loop

  cat("\nMCMC finished in", Sys.time() - start, "seconds.")
  cat("\nThe acceptance rate was:", acceptances[n_iter]/n_iter)
  results = list(X=X, acceptance_rate=alpha_values, sample_mean=expected_values)
  return(results)
}

############################################################################
# Some common target functions for testing MCMC
############################################################################
#2d banana
pi_banana = function(x, B=0.03) {
  exp(-x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2)
}

#3d banana
pi_banana3 = function(x, B=0.03) {
  exp(-(x[1]^2)/200 - (1/2)*(x[2]+B*x[1]^2-100*B)^2-(1/2)*(x[3]^2))
}

#d-dimensional banana
pi_banana3 = function(x, B=0.03) {
  exp(-(x[1]^2)/200 - (1/2)*(x[2]+B*x[1]^2-100*B)^2-(1/2)*(sum(x[3:20]^2)))
}

# Uncorrelated Gaussian distribution  (d=2)
pi_norm = function(x) {
  d=2
  dmvnorm(x, mean = rep(0, d), sigma = matrix(c(1,0,0,1),2,2), log = FALSE)
}

# Correlated Gaussian distribution  (d=8)
# Generating a random positive-definite matrix
Posdef <- function (n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}
sigmaPD <- Posdef(n=8)

pi_norm_corr= function(x) {
  d=8
  dmvnorm(x, mean = rep(0, d), sigma = sigmaPD, log = FALSE)
}


#plot(pi_norm_corr(cbind(seq(-5,5,0.1),seq(-5,5,0.1))))

#values<-rmvnorm(n = 500,mean=rep(0,2),sigma=matrix(c(1,0.8,0.8,1),2,2))
#plot(values[,1], values[,2])

############################################################################
# Some common proposal distributions to draw from
############################################################################

#define a multivariate normal proposal function

q_norm = function(mu, sigma) {
  mvrnorm(n = 1, mu = mu, Sigma = sigma)
}

############################################################################
# Testing the MCMC function
############################################################################
n_iter = 100000  #  Total number of iterations
t_adapt = 1000  #  Time step at which adaptation begins
x_1 = rep(5,8)   #  Vector of inital values
adapt = "AM"   #  Choose the type of adaptation. "AM" or "None" currently.
cov_estimator="Shrinkage"

X = mcmc(target = pi_norm_corr, n_iter = n_iter, x_1 = x_1, adapt=adapt, cov_estimator, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM markov chain")
#hist(X[,1], xlab="First coordinate of the AM markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold')

plot(iterations,X$acceptance_rate)
plot(iterations,X$sample_mean[,2])
