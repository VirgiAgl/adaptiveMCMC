## ----initialisations, eval=FALSE-----------------------------------------
#  
#  d = 8                 #  Dimension of the parameter space
#  n_iter = 20000        #  Total number of iterations
#  x1 = rep(5,d)         #  Vector of inital values
#  t_adapt = 2000        #  When to start adapting
#  adapt = "AM"          #  Algorithm to run
#  target = pi_norm_corr #  Target distribution function

## ----initialisations2, eval=FALSE----------------------------------------
#  
#  X = mcmc(target = target,
#           n_iter = n_iter,
#           x_1 = x1,
#           adapt=adapt,
#           t_adapt = t_adapt
#           )

## ----mh_on_correlated_Gaussian_plots,out.width='0.7\\linewidth',out.height='0.8\\linewidth',fig.align="center",fig.cap="The first (top) and second (bottom) components of the Markov chain resulting from Metropolis-Hastings (MH) targetting a correlated 8-dimensional Gaussian distribution.", echo=FALSE----

load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_MH.rda')

plot1=plotIterations(
                    X_MH$X[,1],
                    n_iter,
                    title="1st component"
                    )
plot2=plotIterations(
                    X_MH$X[,2],
                    n_iter,
                    title="2nd component"
                    )
grid.arrange(plot1, plot2, nrow=2)


## ----mh_on_correlated_Gaussian_plots2,out.width='0.7\\linewidth',out.height='0.4\\linewidth',fig.align="center",fig.cap="The relationship between the first two components of the Markov chain generated using Metropolis-Hasting targeting the correlated Gaussian distribution.", echo=FALSE----

plot3=plotComponents(
                    X_MH$X[,1],
                    X_MH$X[,2],
                    Xtitle = "1st component",
                    Ytitle = "2nd component"
                    )
plot3

## ----am_on_correlated_Gaussian,out.width='0.7\\linewidth',out.height='0.8\\linewidth',fig.align="center",fig.cap="The first (top) and second (bottom) components of the Markov chain resulting from the AM algorithm targetting a correlated 8-dimensional Gaussian distribution.", echo=FALSE----

load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM.rda')

plot6=plotIterations(
                    X_AM$X[,1],
                    n_iter,
                    title="1st component"
                    )
plot7=plotIterations(
                    X_AM$X[,2],
                    n_iter,
                    title="2nd component"
                    )
grid.arrange(plot6, plot7, nrow=2)

## ----am_on_correlated_Gaussian2,out.width='0.7\\linewidth',out.height='0.4\\linewidth',fig.align="center",fig.cap="The relationship between the first two components of the Markov chain resulting from the adaptive metropolis algorithm (AM) targeting a correlated 8-dimensional Gaussian distribution.", echo=FALSE----

plot8=plotComponents(
                    X_AM$X[,1],
                    X_AM$X[,2],
                    Xtitle = "1st component",
                    Ytitle = "2nd component"
                    )
plot8

## ----mh_am_banana,out.width='0.7\\linewidth',out.height='0.8\\linewidth',fig.align="center",fig.cap="The first and second components of the samples resulting from MH algorithm (left) and AM algorithm (right) targeting an 8-dimensional banana shaped distribution.", echo=FALSE----

load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_MH_banana.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM_banana.rda')

plot11=plotComponents(
                    X_MH_banana$X[,2],
                    X_MH_banana$X[,1],
                    Xtitle = "1st component",
                    Ytitle = "2nd component"
                    )

plot12=plotComponents(
                    X_AM_banana$X[,2],
                    X_AM_banana$X[,1],
                    Xtitle = "1st component",
                    Ytitle = "2nd component"
                    )
grid.arrange(plot11, plot12, ncol=2)



## ----am1_am2_diffCov,out.width='0.9\\linewidth',out.height='0.9\\linewidth',fig.align="center",fig.cap="The acceptance probability for the implemented algorithms with a 8-dimensional correlated Gaussian target distribution. The optimal acceptance probability is also shown.", echo=FALSE----

load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM2.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM_sh.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM_th.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM2_th.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM2_sh.rda')



iteration=seq(from= 1, to=n_iter+1, by=1)


df1<-data.frame(iterations=iteration,acceptance_probability=X_MH$acceptance_rates)
df2<-data.frame(iterations=iteration,acceptance_probability=X_AM$acceptance_rates)
df3<-data.frame(iterations=iteration,acceptance_probability=X_AM2$acceptance_rates)
df4<-data.frame(iterations=iteration,acceptance_probability=X_AM_sh$acceptance_rates)
df5<-data.frame(iterations=iteration,acceptance_probability=X_AM2_sh$acceptance_rates)
df6<-data.frame(iterations=iteration,acceptance_probability=X_AM_th$acceptance_rates)
df7<-data.frame(iterations=iteration,acceptance_probability=X_AM2_th$acceptance_rates)
df8<-data.frame(iterations=iteration,acceptance_probability=rep(0.234,n_iter+1))

ggplot(df1,aes(iteration,acceptance_probability))+geom_line(aes(color="MH"))+
      geom_line(data=df2,aes(color="AM"))+
      geom_line(data=df3,aes(color="AM2"))+
      geom_line(data=df4,aes(color="AM + shrinkage"))+
      geom_line(data=df5,aes(color="AM2 + shrinkage"))+
      geom_line(data=df6,aes(color="AM + thresholding"))+
      geom_line(data=df7,aes(color="AM2 + thresholding"))+
      geom_line(data=df8,aes(color="Optimal acceptance"))+
      labs(color="Alorithm:")


## ----am1_am2_banana_diffCov,out.width='0.9\\linewidth',out.height='0.9\\linewidth',fig.align="center",fig.cap="The acceptance probability for the implemented algorithms with a multidimensional banana shaped target distribution. ", echo=FALSE----


load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM_sh_banana.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM_th_banana.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM2_th_banana.rda')
load(file='/data/toucanet/aglietti/adaptiveMCMC/R/data/X_AM2_sh_banana.rda')



iteration=seq(from= 1, to=n_iter+1, by=1)

df1<-data.frame(iterations=iteration,acceptance_probability=X_MH_banana$acceptance_rates)
df2<-data.frame(iterations=iteration,acceptance_probability=X_AM_banana$acceptance_rates)
df3<-data.frame(iterations=iteration,acceptance_probability=X_AM2_banana$acceptance_rates)
df4<-data.frame(iterations=iteration,acceptance_probability=X_AM_sh_banana$acceptance_rates)
df5<-data.frame(iterations=iteration,acceptance_probability=X_AM2_sh_banana$acceptance_rates)
df6<-data.frame(iterations=iteration,acceptance_probability=X_AM_th_banana$acceptance_rates)
df7<-data.frame(iterations=iteration,acceptance_probability=X_AM2_th_banana$acceptance_rates)
df8<-data.frame(iterations=iteration,acceptance_probability=rep(0.234,n_iter+1))

ggplot(df1,aes(iteration,acceptance_probability))+geom_line(aes(color="MH"))+
      geom_line(data=df2,aes(color="AM"))+
      geom_line(data=df3,aes(color="AM2"))+
      geom_line(data=df4,aes(color="AM + shrinkage"))+
      geom_line(data=df5,aes(color="AM2 + shrinkage"))+
      geom_line(data=df6,aes(color="AM + thresholding"))+
      geom_line(data=df7,aes(color="AM2 + thresholding"))+
      geom_line(data=df8,aes(color="Optimal acceptance"))+
      labs(color="Alorithm:")


