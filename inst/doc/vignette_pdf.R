## ----chunkname-----------------------------------------------------------
library(adaptiveMH)

n_iter = 100  #  Total number of iterations
t_adapt = 50  #  Time step at which adaptation begins
x_1 = rep(-10,2) #  Vector of inital values
adapt = "AM"  #  Choose the type of adaptation.

X = mcmc(target = pi_norm, 
	 n_iter = n_iter, 
	 x_1 = x_1, 
	 adapt=adapt, 
	 t_adapt = t_adapt
	)
plot(X[,1], 
     xlab="Iteration index", 
     ylab="First coordinate of the AM Markov chain"
    )

hist(X[,1], 
     xlab="First coordinate of the AM Markov chain", 
     ylab="Frequency", 
     breaks=50, 
     main=""
    )

library(ggplot2)

iterations<-c(seq(from= 1, to=n_iter+1, by=1))

plot<-qplot(iterations, 
	   X[,2],
	   geom="point", 
	   main="Adaptive MCMC for multidimensional correlated Gaussian distribution", 
	   xlab="Iteration", 
	   ylab="", 
 	   xlim=c(0,10000), 
	   alpha = I(1/50),
	   size=I(2.5)) +  theme(axis.title.x = element_text(size = 15), 
	   title = element_text(colour='black'), 
	   axis.text.x=element_text(colour="black"),
 	   axis.text.y=element_text(colour="black"), 
	   axis.line = element_line(colour="black", 
	   			    size = 1, 
				    linetype = "solid"
				   )
				)+ 
				theme(axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid"))
plot



## ----chunkname2----------------------------------------------------------
adapt = "None"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_norm, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")


## ----chunk3--------------------------------------------------------------
adapt = "AM"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_banana, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
#hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold', xlab="First coordinate of the AM markov chain",  ylab="Second coordinate of the AM Markov chain")


## ----chunk4--------------------------------------------------------------
adapt = "None"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_banana, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
#hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold', xlab="First coordinate of the AM markov chain",  ylab="Second coordinate of the AM markov chain")


## ------------------------------------------------------------------------
n_iter = 1000 #  Total number of iterations
x_1 = rep(1,20) #  Vector of inital values
adapt = "AM2"  #  Choose the type of adaptation. "AM" or "None" currently.

X = mcmc(target = pi_banana3, n_iter = n_iter, x_1 = x_1, adapt=adapt, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM Markov chain")
#hist(X[,1], xlab="First coordinate of the AM Markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold', xlab="First coordinate of the AM markov chain",  ylab="Second coordinate of the AM Markov chain")
library(scatterplot3d)
scatterplot3d(x=X[,2],y=X[,3],z=X[,1], xlab="Second coordinate of the AM2 Markov chain", ylab="Third coordinate of the AM2 Markov chain", zlab="First coordinate of the AM2 Markov chain", color="gold")

