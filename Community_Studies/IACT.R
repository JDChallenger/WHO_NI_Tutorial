library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one

#By default, the study has 9 arms and 18 compartments
iact_data <- IACT_sim(n_day = 28, n_mosq = 15, verbose = T,
                      mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25),
                      sigma_net = 0.9, n_nets = 30)
str(iact_data)
head(iact_data)

#Make 'net' a factor variable. Change the levels of the factor

iact_data$net <- as.factor(iact_data$net)
iact_data$net <- relevel(iact_data$net, 'N3')
levels(iact_data$net)

fit <- glm(
  cbind(tot_dead, total - tot_dead) ~ net + day + sleeper, #  + compartment. As arms don't move
  family = binomial, data = iact_data
)
summary(fit)

OR <- exp(coef(summary(fit))['netN6',"Estimate"])
OR_lower <- exp(coef(summary(fit))['netN6',"Estimate"] - 
                   1.96*coef(summary(fit))['netN6','Std. Error'])
OR_upper <- exp(coef(summary(fit))['netN6',"Estimate"] + 
                   1.96*coef(summary(fit))['netN6','Std. Error'])

sum(iact_data[iact_data$net=='N3',]$tot_dead)/sum(iact_data[iact_data$net=='N3',]$total)
FIC_mortality <- sum(iact_data[iact_data$net=='N3',]$tot_dead)/sum(iact_data[iact_data$net=='N3',]$total)
non_inf_margin <- ((FIC_mortality - 0.07) / (1- (FIC_mortality - 0.07))) / (FIC_mortality / (1- FIC_mortality)) 

#Make non-inferiority assessment
OR_lower > non_inf_margin 

################################################################################
####################  Simulation-based power calculation   #####################
################################################################################

#Here, we will use a function to simulate the synthetic datasets, and another 
# function to carryout the non-inferiority assessment

nsim <- 1000 # the number of simulations to carry out
store_power <- rep(NA, nsim) # a container for the outcome of each simulated study
for(i in 1:nsim){
  iact_data <- IACT_sim(sigma_net = 0.9, n_day = 44, n_mosq = 20, 
                       verbose = F, n_nets = 30)
  store_power[i] <- IACT_NIM(dataset = iact_data, verbose = F)
}

#Power estimate:
mean(store_power)
#95% confidence intervals for power estimate
binom.test(table(factor(store_power,c(1,0))))$conf.int

#Now we will make Figure 4.1 in the tutorial. Note: this figure used a higher number
# of simulations, to further reduce stochastic variation. But note that changing 
# nsim to 5000 will make this very slow to run (by default nsim=1000 here)

nsim <- 1000
di <- data.frame()
dayz <- seq(28,76,16)
num_mosq <- seq(15,25,5)
vblist <- c(T,rep(F,nsim-1)) # list for the verbose option. Hence we will see the details of the study for the first simulation only

for(l in 1:length(num_mosq)){ # no. of mosquitoes
  for(k in 1:length(dayz)){ # days
    
    store_sigma_net <- rep(0,9)
    store_sigma_net_ci1 <- rep(0,9)
    store_sigma_net_ci2 <- rep(0,9)
    for(j in 1:length(store_sigma_net)){
      
      store_power <- rep(NA, nsim)
      for(i in 1:nsim){
        mosdata <- IACT_sim(sigma_net = 0.1*j, n_day = dayz[k], n_mosq = num_mosq[l], 
                             verbose = vblist[i], n_nets = 30)
        store_power[i] <- IACT_NIM(dataset = mosdata, verbose = vblist[i])
      }
      print(paste0('J',j))
      
      store_sigma_net[j] <- mean(store_power)
      store_sigma_net_ci1[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[1]
      store_sigma_net_ci2[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[2]
    }
    aux <- data.frame('sigma_net' = 0.1*(1:9), 'power' = store_sigma_net,
                      'power_ci1' = store_sigma_net_ci1, 
                      'power_ci2' = store_sigma_net_ci2, 
                      'days' = rep(dayz[k], length(store_sigma_net_ci1)),
                      'num_mosq' = rep(num_mosq[l], length(store_sigma_net_ci1))
    )
    di <- rbind(di, aux)
    print(paste0('k equal: ',k))
  }
  print(paste0('l equal: ',l))
}
head(di)
date()  

f_names <- c(
  `15` = "15 mosquitoes",
  `20` = "20 mosquitoes",
  `25` = "25 mosquitoes"
)

ggplot(di) + geom_hline(yintercept = 80, alpha = .25) + 
  geom_line(aes(x = sigma_net, y = 100*power, color = factor(days))) +
  labs(color = 'No. of\nStudy Days', fill = 'No. of\nStudy Days') +
  #scale_x_continuous(breaks = seq(0,0.9,0.3)) +
  ylab('Statistical Power (%)') + xlab('Between ITN heterogeneity (s.d.)') +
  theme_classic() + facet_wrap(~num_mosq, labeller = as_labeller(f_names)) + 
  ylim(c(50,100)) + geom_ribbon(aes(x = sigma_net, 
      ymin = 100*power_ci1, ymax = 100*power_ci2, fill = factor(days)),alpha = .3)
