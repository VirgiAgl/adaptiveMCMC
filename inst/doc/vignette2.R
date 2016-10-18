## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, out.width = '600px', dpi = 200, out.height = '400px', fig.align='center',fig.width=6, fig.height=6, size='normalsize')

## ------------------------------------------------------------------------

library(adaptiveMH)

n_iter = 10000  #  Total number of iterations
t_adapt = 1000  #  Time step at which adaptation begins
x_1 = rep(-10,2) #  Vector of inital values
adapt = "AM"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_norm, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)
color = c(217)
color_transparant = adjustcolor(color, alpha.f=0.3)
plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain", 
     col=color)
hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")


## ------------------------------------------------------------------------

adapt = "None"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_norm, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")


## ------------------------------------------------------------------------

adapt = "AM"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_banana, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
#hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold', xlab="First coordinate of the AM markov chain",  ylab="Second coordinate of the AM Markov chain")

## ------------------------------------------------------------------------

adapt = "None"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_banana, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
#hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold', xlab="First coordinate of the AM markov chain",  ylab="Second coordinate of the AM markov chain")

## ------------------------------------------------------------------------

n_iter = 10000 #  Total number of iterations
x_1 = rep(1,20) #  Vector of inital values
adapt = "AM2"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_banana3, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
#hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold', xlab="First coordinate of the AM markov chain",  ylab="Second coordinate of the AM Markov chain")
library(scatterplot3d)
scatterplot3d(x=X[,2],y=X[,3],z=X[,1], xlab="Second coordinate of the AM2 Markov chain", ylab="Third coordinate of the AM2 Markov chain", zlab="First coordinate of the AM2 Markov chain", color="gold")


