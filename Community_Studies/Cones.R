library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one
cone_data <- cone_sim(reps = 4, npos = 3, n_nets = 30, 
                       verbose = T, num_mosq = 5, nday = 30)
head(cone_data)
str(cone_data)

fit <- glm(
  cbind(tot_dead, total - tot_dead) ~ arm + day,
  family = binomial, data = cone_data
)
summary(fit)

#then do NIM 
OR <- exp(coef(summary(fit))['armB',"Estimate"])
OR_lower <- exp(coef(summary(fit))['armB',"Estimate"] - 
                   1.96*coef(summary(fit))['armB','Std. Error'])
OR_upper <- exp(coef(summary(fit))['armB',"Estimate"] + 
                   1.96*coef(summary(fit))['armB','Std. Error'])

sum(cone_data[cone_data$arm=='A',]$tot_dead)/sum(cone_data[cone_data$arm=='A',]$total)
FIC_mortality <- sum(cone_data[cone_data$arm=='A',]$tot_dead)/sum(cone_data[cone_data$arm=='A',]$total)
non_inf_margin <- ((FIC_mortality - 0.07) / (1- (FIC_mortality - 0.07))) / (FIC_mortality / (1- FIC_mortality)) 
c(OR_lower, OR, OR_upper)

OR_lower > non_inf_margin #Non-inferiority assessment

################################################################################
####################  Simulation-based power calculation   #####################
################################################################################

#Here, we will use a function to simulate the synthetic datasets, and another 
# function to carryout the non-inferiority assessment

nsim <- 1000
store_power <- rep(NA, nsim) # a container for the outcome of each simulated study

for(i in 1:nsim){
 cone_data <- cone_sim(cone_mort = 0.4, npos = 3, rep = 4, n_nets = 30,
                       nday = 20, verbose = F)
 store_power[i] <- cone_NIM(dataset = cone_data, verbose = F)
}

#Power estimate:
mean(store_power)
#95% confidence intervals for power estimate
binom.test(table(factor(store_power,c(1,0))))$conf.int

################################################################################
#Now remake Figure 5.1 from the tutorial

#Note: this will be quite slow!

nsim <- 1000
dc <- data.frame() # empty data frame to store results
num_rep <- c(3,4,5) # number of replicates to perform
morts <- seq(0.1,0.4,0.1)

for(l in 1:length(num_rep)){
  for(k in 1:length(morts)){ # mortality
    
    #Tip: you could speed up the simulations, by choosing fewer values for sigma_net. E.g. seq(0.1,0.9,0.2)-- 5 values, or c(0.3,0.6,0.9)-- 3 values
    store_sigma_net <- rep(0,9)
    store_sigma_net_ci1 <- rep(0,9)
    store_sigma_net_ci2 <- rep(0,9)
    for(j in 1:length(store_sigma_net)){
      
      store_power <- rep(NA, nsim)
      for(i in 1:nsim){
        test <- cone_sim(sigma_net = 0.1*j, n_nets = 30, reps = num_rep[l], npos = 3,
                         cone_mort = morts[k], num_mosq = 5, verbose = F, nday = 16)
        store_power[i] <- cone_NIM(dataset = test, verbose = F)
      }
      print(paste0('J',j))
      
      store_sigma_net[j] <- mean(store_power)
      store_sigma_net_ci1[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[1]
      store_sigma_net_ci2[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[2]
    }
    aux <- data.frame('sigma_net' = 0.1*(1:9), 'power' = store_sigma_net,
                      'power_ci1' = store_sigma_net_ci1, 
                      'power_ci2' = store_sigma_net_ci2, 
                      'mort' = rep(morts[k], length(store_sigma_net_ci1)),
                      'reps' = rep(num_rep[l], length(store_sigma_net_ci1))
    )
    dc <- rbind(dc, aux)
    print(paste0('k equal: ',k))
  }
  print(paste0('l equal: ',l))
}

head(dc)

#labels for the facets
m_names <- c(
  `0.1` = "Mosquito mortality = 0.1",
  `0.2` = "Mosquito mortality = 0.2",
  `0.3` = "Mosquito mortality = 0.3",
  `0.4` = "Mosquito mortality = 0.4"
)

ggplot(dc) + geom_line(aes(x = sigma_net, y = 100*power, color = factor(reps))) + 
  geom_hline(yintercept = 80) +
  facet_wrap(~mort, labeller = as_labeller(m_names)) + theme_classic() + 
  labs(color = 'Number of\nReplicates') + 
  xlab('Between-net heterogeneity (s.d.)') + ylab('Statistical power (%)')
