################################################################################
############### Installation instructions for the rethinking package ###########
################################################################################

#Unless you are an experience Rstan user, we recommend following the simplified
# installation instructions, by running the following commands:

#First install some auxiliary packages (UNCOMMENT THIS COMMAND)
#install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
#Then install the 'lite' version of the package (UNCOMMENT THIS COMMAND)
#devtools::install_github("rmcelreath/rethinking@slim")

#Once the packages are installed, you can reinstate the comment symbols above (so that 
# you don't re-install the packages every time you run the script)

################################################################################

library(ggplot2)
library(rethinking)



#update df name
# Let's load some hypothetical field data 
# (this matches the example in the tutorial- see Table 6.1 & Figure 6.1A)
# N: 
dx2 <- data.frame('t' = c(0.5,1,2,3),'N' = c(96,93,86,80), 'surv' = c(86,66,38,22))

# For visualisation purposes, we may wish to add confidence intervals here
# Note: this is not necessary for the calculation of the median survival
dx2$ci1 <- NA
dx2$ci2 <- NA
for(i in 1:dim(dx2)[1]){
  dx2$ci1[i] <- binom.test(c(dx2$surv[i], dx2$N[i] - dx2$surv[i]))$conf.int[1]
  dx2$ci2[i] <- binom.test(c(dx2$surv[i], dx2$N[i] - dx2$surv[i]))$conf.int[2]
}  

#Now, let's fit the S-shaped curve (Equation 6.1 in the tutorial),
# using the Bayesian model defined in Equation 6.2

m2 <- quap(
  alist(
    surv ~ dbinom(N, p),
    p <- exp(20 - 20/(1 - (t/L)*(t/L))), # value of k fixed to 20
    L ~ dunif(5.5,20.7) # prior distribution for L
  ),
  data = dx2#, chains = 3, cores = 1, iter = 3600
)
precis(m2)

#Those who have installed the full rethinking package (see Appendix A of the tutorial)

# m2 <- ulam(
#   alist(
#     surv ~ dbinom(N, p),
#     p <- exp(20 - 20/(1 - (t/L)*(t/L))),
#     L ~ dunif(5.5,20.7)
#   ),
#   data = dx2, chains = 3, cores = 1, iter = 3600
# )
# precis(m2)
# stancode(m2) # this function allows us to view the raw stan code for this ulam() model

#Whichever method you've used, we now need to extract posterior samples from the fitted model
#We calculate the value of t_{50} using Equation 6.3

es2 <- extract.samples(m2)
store <- rep(0,5000) # extract 5000 posterior samples
for(i in 1:5000){
  #calculate Eq. 6.3 using each posterior sample for parameter L
  store[i] <- es2$L[i]*sqrt(log(2)/(log(2) + 20))
}
x1 <- mean(store) # Mean estimate of t_50
x2 <- quantile(store, probs = c(0.025,0.5,0.975)) # 95% credible interval for t_50

# This finishes the calculation of t_50. Below we repeat the steps required to reproduce 
# Figure 6.1B (the visualisation of median survival) for this dataset

# To generate the black curve (and the grey shaed area, showing the uncertainty)
# that represents the fitted S-shaped curve, we need to estimate its value at a range of timepoints

#Use the 5000 samples to estimate the curve at 500 time points
storeM <- matrix(0, nrow = 5000, ncol = 500)
for(j in 1:500){
  for(i in 1:5000){
    years <- (j-1)/100 # a value of time (in years)
    storeM[i,j] <- exp(20 - 20/(1 - (years/es2$L[i])*(years/es2$L[i])))
  }
}
MSF_curve <- apply(storeM, 2, mean)
MSF_CrI <- apply(storeM, 2, HPDI , prob=0.95)
# save results to a data frame
df <- data.frame('t' = seq(0,499)/100, 'mean' = MSF_curve,
                 'ci1' = MSF_CrI[1,], 'ci2' = MSF_CrI[2,])

#Let's build up the plot, starting with the black curve
MSF_plot <- ggplot() + geom_line(data = df, aes(x=t,y=mean)) + 
  theme_classic() +  
  ylab('Functional survival of ITNs (proportion)') + 
  xlab('Time since ITNs distributed (years)')
MSF_plot

#Add uncertainty in the curve
MSF_plot <- MSF_plot + geom_ribbon(data = df, aes(x=t,ymin = ci1, ymax = ci2), alpha = .3)
MSF_plot

#Add the data points
MSF_plot <- MSF_plot + 
geom_point(data = dx2, aes(x = t, y = surv/N), color = 'slateblue') + 
  geom_errorbar(data = dx2, aes(x = t, ymin = ci1, ymax = ci2), 
                color = 'slateblue', width = .1)
MSF_plot 

#Add the estimate of t_50
MSF_plot <- MSF_plot + annotate('point', x = x1, y = 0.5, color = 'green') + 
  annotate('segment', x = x2[1], xend = x2[3], y = 0.5, yend = 0.5, color = 'green')
MSF_plot 

#you may wish to annotate the plot with the esimated t_50 value 
MSF_plot <- MSF_plot +
annotate('text', label = expression(paste('Measurement of ',t[50],' =')), x = 3.5, y = 0.9) + 
  annotate('text', label = paste(  round(x1,3),'[',round(x2[1],3),',',round(x2[3],3),']'), x = 3.5, y = 0.83) 
MSF_plot 


################################################################################
###### Simualation procedure, to estimate precision of the t_50 estimate #######
################################################################################

# For this procedure, we are not starting with a dataset. This means we need to think about
# how nets are lost to follow up over time

## We assume constant loss of nets over time (exponential decay)
## Here, 'pc_loss' is the percentage of nets lost by the end of the follow up period
## 'duration' determines the follow-up period (in years)
## This function returns the the 'loss rate': this enables us to estimate the
## percentage of nets remaining at intermediate time points
calculate_r <- function(pc_loss, duration = 3){
  duration/(log(100) - log(100 - pc_loss))
}

#so, if 20% of ITNs are lost over three years...
rc <- calculate_r(pc_loss = 20, duration = 3)

#... we can estimate net loss as a continuous function of time
ggplot() + geom_function(fun = function(x) exp(-x/rc)) + theme_classic() + 
  xlim(c(0,4)) + geom_vline(xintercept = 3, color = 'orange2') + 
  geom_hline(yintercept = 0.8, color = 'orange2') + 
  geom_vline(xintercept = c(0.5,1,2), alpha = .3, color = 'cyan')

#so how many nets would one expect to be remaining after 6/12/24 months?
#(Here we write this as 0.5/1/2 years)
exp(-0.5/rc)
exp(-1/rc)
exp(-2/rc)

#This function 'loss_sim' will do the following tasks: 
# 1) for a given value of loss rate over 3 years ('loss_rate'), it will estimate
#  the number of nets lost at earlier time points (0.5,1,2 years);
# 2) For a user-chosen choice of the true values of the proportion of nets remaining
#  in use at each point, binomial sampling is used to generate a synthetic dataset
#  from a community study
# 3) The Bayesian model is fitted, to estimate L (and thus t_{50}) 

loss_sim <- function(start_samp = 100, loss_rate = 10, verbose = F,
                     surv_prob = c(0.88,.63,.48,.27), full_model = F){
  rc <- calculate_r(pc_loss = loss_rate, duration = 3)
  
  #so how many nets would one expect to be remaining after 6/12/24 mths?
  bray <- round(start_samp*c(exp(-0.5/rc), exp(-1/rc), exp(-2/rc), exp(-3/rc)))
  dxh <- data.frame('t' = c(0.5,1,2,3), 'N' = bray)
  
  dxh$surv <- rbinom(4,dxh$N, surv_prob)
  
  if(verbose == T){
    print(dxh)
  }
  
  if(full_model == T){
    m2 <- ulam(
      alist(
        surv ~ dbinom(N, p),
        p <- exp(20 - 20/(1 - (t/L)*(t/L))),
        L ~ dunif(5.5,20.7)
      ),
      data = dxh, chains = 3, cores = 3, iter = 4000
    )
  }else{
    m2 <- quap(
      alist(
        surv ~ dbinom(N, p),
        p <- exp(20 - 20/(1 - (t/L)*(t/L))),
        L ~ dunif(5.5,20.7)
      ),
      data = dxh#, chains = 3, cores = 3, iter = 4000
    )
  }
  
  if(verbose == T){
    print(precis(m2))    
  }
  
  es2 <- extract.samples(m2)
  store <- rep(0,5000)
  for(i in 1:5000){
    store[i] <- es2$L[i]*sqrt(log(2)/(log(2) + 20))
  }
  
  return(c(mean(store),
           quantile(store, probs = c(0.025))[[1]],quantile(store, probs = c(0.5))[[1]],
           quantile(store, probs = c(0.975))[[1]],
           quantile(store, probs = c(0.975))[[1]] - quantile(store, probs = c(0.025))[[1]],
           start_samp, loss_rate))
}

loss_sim(full_model = F, verbose = T, start_samp = 100, loss_rate = 20)

loss_vec <- seq(5,30,5)
samp_vec <- seq(30,120,10)
storeM <- matrix(0,nrow = length(samp_vec)*length(loss_vec), ncol = 7)
count <- 1
for(i in 1:length(samp_vec)){
  for(j in 1:length(loss_vec)){
    storeM[count,] <- loss_sim(start_samp = samp_vec[i], loss_rate = loss_vec[j],
                               full_model = F)
    print(paste0('count equals: ',count))
    count <- count + 1
  }
}
storeMDF <- as.data.frame(storeM)
colnames(storeMDF) <- c('mean','cr1','med','cr2','prec', 'sample','loss')

ggplot(storeMDF) + geom_point(aes(x = sample, y = prec, color = factor(loss))) + 
  facet_wrap(~loss) + theme_classic()

#this can be quite noisy, largely due to the binomial sampling.
# If we wish, we could average over this, by repeating the process (say, 50 times),
# for each combination of parameter values

loss_vec <- seq(5,20,5)
samp_vec <- seq(100,300,50)
storeM2 <- matrix(0,nrow = length(samp_vec)*length(loss_vec)*50, ncol = 7)
count <- 1
for(i in 1:length(samp_vec)){
  for(j in 1:length(loss_vec)){
    for(k in 1:50){
      storeM2[count,] <- loss_sim(start_samp = samp_vec[i], 
                                  loss_rate = loss_vec[j], full_model = F)
      print(paste0('count equals: ',count))
      count <- count + 1
    }
  }
}
storeM2DF <- as.data.frame(storeM2)
colnames(storeM2DF) <- c('mean','cr1','med','cr2','prec', 'sample','loss')

#Post-process the data
strM <- matrix(0,ncol = 5, nrow = length(samp_vec)*length(loss_vec))
count <- 1
for(i in 1:length(samp_vec)){
  for(j in 1:length(loss_vec)){
    #Select all outputs for this combination of parameter values
    aux <- storeM2DF[storeM2DF$sample==samp_vec[i] & storeM2DF$loss==loss_vec[j],]
    aux1 <- mean(aux$prec)#mean
    aux2 <- quantile(aux$prec, probs = c(0.025,0.975))
    aux3 <- aux2[2] - aux2[1]
    strM[count,] <- c(aux1, aux2[1], aux2[2], samp_vec[i], loss_vec[j])
    count <- count + 1  
  }
}
strMDF <- as.data.frame(strM)
colnames(strMDF) <- c('mean_prec', 'ci1', 'ci2', 'sample', 'loss')

ggplot(strMDF) + 
  theme_classic() + facet_wrap(~loss) + 
  geom_errorbar(aes(x = sample, ymin = ci1, ymax = ci2), alpha = .8, color = 'grey', width = .1) + 
  geom_point(aes(x = sample, y = mean_prec, color = factor(loss))) + 
  labs(color = 'Loss to\nfollow up (%)') + xlab('Starting sample size (no. of ITNs)') + 
  ylab(expression(paste("Precision of measurement of ", t[50],' (years)',sep=""))) + 
  ggtitle(expression(paste('Precision of ', t[50],' measurement'))) + ylim(c(0,NA)) + xlim(c(85,NA)) + 
  theme(legend.position = 'bottom')







