############################################################################
# Multidimensional correlated Guassian distribution d=8
############################################################################
n_iter = 10000     #  Total number of iterations
t_adapt = 2000     #  Time step at which adaptation begins
x_1 = rep(5,8)     #  Vector of inital values
adapt = "None"     #  Choose the type of adaptation. "None", "AM" or "AM2".
cov_estimator=""   #  Choose the covariance matrix estimators. "Sample covariance", "Shrinkage estimator" or "Thresholding estimator".


X = mcmc(target = pi_norm_corr, n_iter = n_iter, x_1 = x_1, adapt=adapt, cov_estimator, t_adapt = t_adapt)

#plot(X[,1],xlab="Iteration index", ylab="First coordinate of the AM markov chain")
#hist(X[,1], xlab="First coordinate of the AM markov chain", ylab="Frequency", breaks=50, main="")
plot(X[,2],X[,1],col='gold')
plot(X[,1],col='gold')

iterations<-c(seq(from= 1, to=n_iter+1, by=1))

plot1<-qplot(iterations, X[,1],geom="point", main="Adaptive MCMC for multidimensional correlated Gaussian distribution", xlab="Iteration", ylab="", xlim=c(0,100002), alpha = I(1/50), size=I(2.5)) +  theme(axis.title.x = element_text(size = 15), title = element_text(colour='black'), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.line = element_line(colour="black", size = 1, linetype = "solid"))+ theme( axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid")) + geom_point()

plot2<-qplot(iterations, X[,2],geom="point", main="Adaptive MCMC for multidimensional correlated Gaussian distribution", xlab="Iteration", ylab="", xlim=c(0,100002), alpha = I(1/50), size=I(2.5)) +  theme(axis.title.x = element_text(size = 15), title = element_text(colour='black'), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.line = element_line(colour="black", size = 1, linetype = "solid"))+ theme( axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid")) + geom_point()

plot3<-qplot(x=X[,1], main="Adaptive MCMC for multidimensional correlated Gaussian distribution", xlab="values", ylab="Frequency", ylim=c(0,100002), alpha = I(1/50), size=I(2.5)) +  theme(axis.title.x = element_text(size = 15), title = element_text(colour='black'), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.line = element_line(colour="black", size = 1, linetype = "solid")) + geom_histogram(bins = 200)

plot4<-qplot(X[,2], main="Adaptive MCMC for multidimensional correlated Gaussian distribution", xlab="Iteration", ylab="", xlim=c(0,100002), alpha = I(1/50), size=I(2.5), binwidth=30) +  theme(axis.title.x = element_text(size = 15), title = element_text(colour='black'), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"), axis.line = element_line(colour="black", size = 1, linetype = "solid"))+ theme( axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid")) + geom_point()
