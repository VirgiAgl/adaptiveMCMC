## ----mh_on_correlated_guassian-------------------------------------------

d = 8
n_iter = 200     #  Total number of iterations
x_1 = rep(5,d)      #  Vector of inital values
t_adapt = 50      #  When to start adapting
adapt = "None"      #  Choose the type of adaptation.

X_MH = mcmc(target = pi_norm_corr,
         n_iter = n_iter,
         x_1 = x_1,
         adapt=adapt)

## ----mh_on_correlated_guassian_plots,out.width='0.7\\textwidth',fig.cap="This is the first graph", echo=FALSE----

plot1=plotIterations(
                    X_MH$X[,1],
                    n_iter,
                    title="Values of the first MC component"
                    )
plot2=plotIterations(
                    X_MH$X[,2],
                    n_iter,
                    title="Values of the second MC component"
                    )
grid.arrange(plot1, plot2, ncol=2)


plot3=plotComponents(
                    X_MH$X[,1],
                    X_MH$X[,2],
                    Xtitle = "First component of the MC",
                    Ytitle = "Second component of the MC"
                    )
plot3

plot4=plotIterations(X_MH$acceptance_rates, n_iter, "Alpha")
plot4

plot5=plotIterations(X_MH$sample_mean[,1], n_iter, "Sample expected values")
plot5



## ----am_on_correlated_guassian-------------------------------------------

adapt = "AM"       #  Choose the type of adaptation.
cov_estimator="Sample covariance"    #  Choose the type of covariance matrix estimator.

X_AM = mcmc(target = pi_norm_corr,
         n_iter = n_iter,
         x_1 = x_1,
         adapt=adapt,
         t_adapt = t_adapt,
         cov_estimator=cov_estimator
         )

plot6=plotIterations(
                    X_AM$X[,1],
                    n_iter,
                    title="Values of the first MC component"
                    )
plot7=plotIterations(
                    X_AM$X[,2],
                    n_iter,
                    title="Values of the second MC component"
                    )
grid.arrange(plot6, plot7, ncol=2)


plot8=plotComponents(
                    X_AM$X[,1],
                    X_AM$X[,2],
                    Xtitle = "First component of the MC",
                    Ytitle = "Second component of the MC"
                    )
plot8

plot9=plotIterations(X_AM$acceptance_rates, n_iter, "Alpha")
plot9

plot10=plotIterations(X_AM$sample_mean[,1], n_iter, "Sample expection for the first component of MC")
plot10



## ----mh_am_banana--------------------------------------------------------

x_1 = rep(5,2)      #  Vector of inital values
target = pi_banana


X_MH_banana = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="None"
         )

plot11=plotComponents(
                    X_MH_banana$X[,2],
                    X_MH_banana$X[,1],
                    Xtitle = "First component of the MC",
                    Ytitle = "Second component of the MC"
                    )

X_AM_banana = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         t_adapt = t_adapt,
         adapt="AM"
         )

plot12=plotComponents(
                    X_AM_banana$X[,2],
                    X_AM_banana$X[,1],
                    Xtitle = "First component of the MC",
                    Ytitle = "Second component of the MC"
                    )
grid.arrange(plot11, plot12, ncol=2)



## ----am1_am2_banana, fig.width=5-----------------------------------------

x_1 = rep(5,8)      #  Vector of inital values
target = pi_banana8


X_AM2_banana = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM2"
         )


plot13=plotIterations(
                    X_AM_banana$X[,1],
                    n_iter,
                    title="Values of the first MC component"
                    )
plot14=plotIterations(
                    X_AM2_banana$X[,1],
                    n_iter,
                    title="Values of the first MC component"
                    )

grid.arrange(plot13, plot14, ncol=2)


## ----am1_am2_banana_diffCov, fig.width=5---------------------------------
x_1 = rep(5,8)      #  Vector of inital values
cov_estimator1="Shrinkage estimator"
cov_estimator2="Thresholding estimator"
target = pi_banana8

X_MH_banana = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="None"
         )


X_AM_banana = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         t_adapt = t_adapt,
         adapt="AM"
         )

X_AM_banana_sh= mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM",
         t_adapt = t_adapt,
         cov_estimator=cov_estimator1
         )

X_AM2_banana_sh = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM2",
         cov_estimator=cov_estimator1
         )

X_AM_banana_th= mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM",
         t_adapt = t_adapt,
         cov_estimator=cov_estimator2
         )

X_AM2_banana_th = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM2",
         cov_estimator=cov_estimator2
         )

iteration=seq(from= 1, to=n_iter+1, by=1)
#data=data.frame(iteration_count=iteration, am_alpha=X_AM_banana_sh$acceptance_rates)


df1<-data.frame(iterations=iteration,acceptance_probability=X_MH_banana$acceptance_rates)
df2<-data.frame(iterations=iteration,acceptance_probability=X_AM_banana$acceptance_rates)
df3<-data.frame(iterations=iteration,acceptance_probability=X_AM2_banana$acceptance_rates)
df4<-data.frame(iterations=iteration,acceptance_probability=X_AM_banana_sh$acceptance_rates)
df5<-data.frame(iterations=iteration,acceptance_probability=X_AM2_banana_sh$acceptance_rates)
df6<-data.frame(iterations=iteration,acceptance_probability=X_AM_banana_th$acceptance_rates)
df7<-data.frame(iterations=iteration,acceptance_probability=X_AM2_banana_th$acceptance_rates)
df8<-data.frame(iterations=iteration,acceptance_probability=rep(0.234,n_iter+1))

ggplot(df1,aes(iterations,acceptance_probability))+geom_line(aes(color="MH"))+
      geom_line(data=df2,aes(color="AM"))+
      geom_line(data=df3,aes(color="AM2"))+
      geom_line(data=df4,aes(color="AM + shrinkage"))+
      geom_line(data=df5,aes(color="AM2 + shrinkage"))+
      geom_line(data=df6,aes(color="AM + thresholding"))+
      geom_line(data=df7,aes(color="AM2 + thresholding"))+
      geom_line(data=df8,aes(color="Optimal acceptance"))+
      labs(color="Alorithm:")

# all lines of the mean in the same plot with a legend


## ----mean_diffCov, fig.width=5-------------------------------------------

target=pi_norm_corr


X_MH = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt=adapt)

X_AM = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt=adapt,
         t_adapt = t_adapt,
         cov_estimator=cov_estimator
         )
X_AM2 = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM2"
         )
X_AM_sh= mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM",
         t_adapt = t_adapt,
         cov_estimator=cov_estimator1
         )

X_AM2_sh = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM2",
         cov_estimator=cov_estimator1
         )

X_AM_th= mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM",
         t_adapt = t_adapt,
         cov_estimator=cov_estimator2
         )

X_AM2_th = mcmc(target = target,
         n_iter = n_iter,
         x_1 = x_1,
         adapt="AM2",
         cov_estimator=cov_estimator2
         )



df1<-data.frame(x=iteration,y=X_MH$sample_mean[,1])
df2<-data.frame(x=iteration,y=X_AM$sample_mean[,1])
df3<-data.frame(x=iteration,y=X_AM2$sample_mean[,1])
df4<-data.frame(x=iteration,y=X_AM_sh$sample_mean[,1])
df5<-data.frame(x=iteration,y=X_AM2_sh$sample_mean[,1])
df6<-data.frame(x=iteration,y=X_AM_th$sample_mean[,1])
df7<-data.frame(x=iteration,y=X_AM2_th$sample_mean[,1])
df8<-data.frame(x=iteration, y=rep(0,n_iter+1))

ggplot(df1,aes(x,y))+geom_line(aes(color="First line"))+
      geom_line(data=df2,aes(color="Second line"))+
      geom_line(data=df3,aes(color="Third line"))+
      geom_line(data=df4,aes(color="Forth line"))+
      geom_line(data=df5,aes(color="Fifth line"))+
      geom_line(data=df6,aes(color="Sixth line"))+
      geom_line(data=df7,aes(color="Seventh line"))+
      geom_line(data=df8,aes(color="Black"))+
      labs(color="Legend text")

# all lines of the mean in the same plot with a legend


