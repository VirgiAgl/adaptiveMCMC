library(MASS)
library(Matrix)
library(mvtnorm)
library(CVTuningCov)
library(corpcor)
library(PDSCE)
library(ggplot2)
library(gridExtra)


############################################################################
# Still to implement
# ------------------
# Burn in
# Suboptimality factor
# Implement and test higher dimensional problem
# Speed up code
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
  acceptances = c(rep(NA,n_iter+1))  #  Initialise vector to store the number of accepted states
  acceptances[1] = 0  #  Store the first number of accepted states
  expected_values = matrix(NA,n_iter+1,d)  #  Initialise matrix to store the expected valeus
  expected_values[1,] = sapply(x_1,mean)  #  Store the first step expected value

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
      if (i <= t_adapt) {
        Y = q_norm(mu=X[i-1,], sigma=sigma) #proposal distribution sampled at the current point
      }

      if (i == t_adapt+1) {
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

      if (i > t_adapt+1) {
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
      numerator = target(Y)*q_norm(X[i-1,],sigma=diag(d))
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

      if (i == t_adapt) {  #the first adaptation step
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
      if (i > t_adapt) {       #the later adaptation step
        X_until_now = X[1:i,]  #the rows of X that we have filled so far

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
    }  #  End of loop

  alpha_values<-acceptances/seq(1,length(acceptances),by=1)

  cat("\nMCMC finished in", Sys.time() - start, "seconds.")
  cat("\nThe last acceptance rate was:", acceptances[n_iter]/n_iter)
  results = list(X=X, acceptance_rates=alpha_values, sample_mean=expected_values, n_acceptances=acceptances)
  return(results)
}

############################################################################
# Some common target functions for testing MCMC
############################################################################
#2d banana
pi_banana = function(x, B=0.03) {
  exp(-x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2)
}
pi_banana_matrix = function(x, B=0.03) {
  exp(-x[,1]^2/200 - 1/2*(x[,2]+B*x[,1]^2-100*B)^2)
}

#3d banana
pi_banana3 = function(x, B=0.03) {
  exp(-(x[1]^2)/200 - (1/2)*(x[2]+B*x[1]^2-100*B)^2-(1/2)*(x[3]^2))
}

#d-dimensional banana
pi_banana8 = function(x, B=0.03) {
  exp(-(x[1]^2)/200 - (1/2)*(x[2]+B*x[1]^2-100*B)^2-(1/2)*(sum(x[3:8]^2)))
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



############################################################################
# Some common proposal distributions to draw from
############################################################################

#Multivariate Gaussian proposal distribution

q_norm = function(mu, sigma) {
  mvrnorm(n = 1, mu = mu, Sigma = sigma)
}


############################################################################
# Define a function to compute moving average acceptance rates
############################################################################

#moving_acceptance = function(n_moving, n_iter,X){
#  moving_alpha = rep(NA, n_iter/n_moving)
#  moving_alpha[1]= X$acceptances[n_moving]/n_moving
#  for (i in seq(2*n_moving, n_iter, by=n_moving)) {
#    moving_alpha[i/n_moving]=(X$acceptances[i]-X$acceptances[i-n_moving])/n_moving
#  }
#  return(moving_alpha)
#}



############################################################################
# Define functions to produce plots
############################################################################

plotTarget = function(target, d, num_values){
  Z = matrix(cbind(rep(seq(from=(-10), to=10,length.out = 1000),d)), ncol=d)
  functionValues = target(Z)
  plot(seq(from=(-10), to=10,length.out = 1000),functionValues)
}


#plotTarget(pi_norm_corr, 8, 100)
#plotTarget(pi_banana_matrix, 2, 100)


plotIterations = function(X, n_iter, title){
  plot = qplot(seq(from= 1, to=n_iter+1, by=1), X ,geom="point",
            main="",
            xlab="Iteration index", ylab=title, xlim=c(0,n_iter+1), alpha = I(1/1000), size=I(2.5)) +
            theme(axis.title.x = element_text(size = 12), title = element_text(colour='black'),
            axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
            axis.line = element_line(colour="black", size = 1, linetype = "solid"))+
            geom_point(size=2, shape=19, alpha=0.1)
  return(plot)
}


plotComponents = function(X, Y , Xtitle, Ytitle){
  plot = qplot(X, Y ,geom="point", main="", xlab=Xtitle, ylab=Ytitle, alpha = I(1/1000), size=I(2.5)) +
    theme(axis.title.x = element_text(size = 12), title = element_text(colour='black'),
          axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.line = element_line(colour="black", size = 1, linetype = "solid"))+
          geom_point(size=2, shape=19, alpha=0.1)
  return(plot)
}

plotComponents3D = function(X, Y , Xtitle, Ytitle){
  plot = qplot(X, Y ,geom="point", main="", xlab=Xtitle, ylab=Ytitle, alpha = I(1/1000), size=I(2.5)) +
    theme(axis.title.x = element_text(size = 12), title = element_text(colour='black'),
          axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.line = element_line(colour="black", size = 1, linetype = "solid"))+
    geom_point(size=2, shape=19, alpha=0.1)
  return(plot)
}

plotACF = function(X, title){
  plot = acf(X, main=title)
}


#
#
# plot1=plotIterations(X$X[,1], 1000, title="Values of the first MC component")
# plot2=plotIterations(X$X[,2], 1000, title="Values of the second MC component")
# grid.arrange(plot1, plot2, ncol=2)
#
# plot3=plotComponents(X$X[,1], X$X[,2],  Xtitle = "First component of the MC", Ytitle = "Second component of the MC")
# plot3
#
# plot4=plotIterations(X$acceptance_rates, 5000, "Alpha")
# #plot4
#
# plot5=plotIterations(X$sample_mean[,1], 5000, 'Sample expected values')
# #plot5
#
# plotACF(X$X[,2], 'Autocorrelation function for the first MC component')
#


# ############################################################################
# # Testing the MCMC function
# ############################################################################
# n_iter = 100000     #  Total number of iterations
# t_adapt = 1000    #  Time step at which adaptation begins
# x_1 = rep(5,8)    #  Vector of inital values
# adapt = "AM"      #  Choose the type of adaptation. "AM" or "None" currently.
# cov_estimator="Sample covariance"    #  Choose the type of covariance matrix estimator. "Sample covariance", "Shrinkage estimator" or "Thresholding estimator".
#
# X = mcmc(target = pi_norm_corr, n_iter = n_iter, x_1 = x_1, adapt=adapt, cov_estimator, t_adapt = t_adapt)
#
# plot1=plotIterations(X$X[,1], 1000, title="Values of the first MC component")
#
#
#
# iterations<-c(seq(from= 1, to=n_iter+1, by=1))  #costruct the sequence of iteration steps


#moving_accep=plot(moving_acceptance(n_moving, n_iter, X))
#n_moving=10
#
# moving_alpha = rep(NA, n_iter/n_moving)
# #  moving_alpha[1]= X$acceptances[n_moving]/n_moving
# #  for (i in seq(2*n_moving, n_iter, by=n_moving)) {
# #    moving_alpha[i/n_moving]=(X$acceptances[i]-X$acceptances[i-n_moving])/n_moving
# #  }
# apply(X_chain[1:10],2,mean)
# apply(X_chain[2:11],2,mean)


#Mov_mean<-function(X,window=2000, col_num=1) {
#  data_to_consider=as.matrix(X$X[,col_num])
#  chain_row=nrow(data_to_consider)
#  chain_col=ncol(data_to_consider)
#  fin_matrix=matrix(NA,ncol = chain_col, nrow =chain_row-window+1 )
#  for (i in 1:(chain_row-window+1)){
#    ext=(i+(window-1))
#    data=as.matrix(data_to_consider[i:ext,])
#    fin_matrix[i,]=apply(data,2,mean)
#  }
#  return(fin_matrix)
#}
