library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one
tunnel_data <- tunnel_sim(npos = 3, n_nets = 30,
                          sigma_net = 0.7, n_tunnels = 5, verbose = T, n_mosq = 50)
head(tunnel_data)
str(tunnel_data)

table(tunnel_data$arm) # two arms, 'A' and 'B'

#Fit GLM. Include fixed effects for day & tunnel. These variables should be 
# factor variables, if they are not already
fit <- glm(
  cbind(tot_dead, total - tot_dead) ~ arm + day + tunnel, #  + compartment. As arms don't move
  family = binomial, data = tunnel_data
)
summary(fit)

#then do NIM 
OR <- exp(coef(summary(fit))['armB',"Estimate"])
OR_lower <- exp(coef(summary(fit))['armB',"Estimate"] - 
                  1.96*coef(summary(fit))['armB','Std. Error'])
OR_upper <- exp(coef(summary(fit))['armB',"Estimate"] + 
                  1.96*coef(summary(fit))['armB','Std. Error'])

sum(tunnel_data[tunnel_data$arm=='A',]$tot_dead)/sum(tunnel_data[tunnel_data$arm=='A',]$total)
FIC_mortality <- sum(tunnel_data[tunnel_data$arm=='A',]$tot_dead)/sum(tunnel_data[tunnel_data$arm=='A',]$total)
non_inf_margin <- ((FIC_mortality - 0.07) / (1- (FIC_mortality - 0.07))) / (FIC_mortality / (1- FIC_mortality)) 
c(OR_lower, OR, OR_upper)

################################################################################
####################  Simulation-based power calculation   #####################
################################################################################

#Here, we will use a function to simulate the synthetic datasets, and another 
# function to carryout the non-inferiority assessment

nsim <- 1000
store_power <- rep(NA, nsim) # a container for the outcome of each simulated study

for(i in 1:nsim){
  tunnel_data <- tunnel_sim(sigma_net = 0.9, n_nets = 30, #reps = 1,
                      mort = 0.3, n_mosq = 45)
  store_power[i] <- tunnel_NIM(data = tunnel_data, verbose = F)
}

#Power estimate:
mean(store_power)
#95% confidence intervals for power estimate
binom.test(table(factor(store_power,c(1,0))))$conf.int

################################################################################
#Now remake Figure 5.2 from the tutorial

dtX2 <- data.frame() # empty dataframe, in which to store results
num_mosq <- seq(10,60,10)
morts <- seq(0.1,0.5,0.1)

for(l in 1:length(num_mosq)){ # mosquito numbers
  for(k in 1:length(morts)){ # mortality
    
    store_sigma_net <- rep(0,9)
    store_sigma_net_ci1 <- rep(0,9)
    store_sigma_net_ci2 <- rep(0,9)
    for(j in 1:length(store_sigma_net)){
      
      store_power <- rep(NA, nsim)
      for(i in 1:nsim){
        test <- tunnel_sim(sigma_net = 0.1*j, n_nets = 30, 
                            mort = morts[k], n_mosq = num_mosq[l])
        #print(head(test))
        #print(i)
        store_power[i] <- tunnel_NIM(data = test, verbose = F)
      }
      #print(paste0('J',j))
      
      store_sigma_net[j] <- mean(store_power)
      store_sigma_net_ci1[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[1]
      store_sigma_net_ci2[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[2]
    }
    aux <- data.frame('sigma_net' = 0.1*(1:9), 'power' = store_sigma_net,
                      'power_ci1' = store_sigma_net_ci1, 
                      'power_ci2' = store_sigma_net_ci2, 
                      'mort' = rep(morts[k], length(store_sigma_net_ci1)),
                      'num_mosq' = rep(num_mosq[l], length(store_sigma_net_ci1))
    )
    dtX2 <- rbind(dtX2, aux)
    print(paste0('k equal: ',k))
  }
  print(paste0('l equal: ',l))
}
head(dtX2)

#VT facet with new labels
f_names <- c(
  `10` = "10 mosquitoes",
  `20` = "20 mosquitoes",
  `30` = "30 mosquitoes",
  `40` = "40 mosquitoes",
  `50` = "50 mosquitoes",
  `60` = "60 mosquitoes"
)

ggplot(dtX2) + geom_line(aes(x = sigma_net, y = 100*power, color = factor(mort))) + 
  facet_wrap(~num_mosq, labeller = as_labeller(f_names)) + theme_classic() + 
  geom_point(aes(x = sigma_net, y = 100*power, color = factor(mort))) + 
  geom_hline(yintercept = 80, color = 'grey') + ylim(c(0,100)) + 
  xlab('Between-net heterogeneity (s.d.)') + 
  ylab('Statistical Power (%)') + labs(color = 'Average\nMosquito\nMortality')

