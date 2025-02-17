library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one

#Update function name
iact_data <- IACT_sim2()
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
  iact_data <- IACT_sim2(sigma_net = 0.9, n_day = 44, n_mosq = 20, 
                       verbose = F, n_nets = 30)
  store_power[i] <- IACT_NIM(dataset = iact_data, verbose = F)
}

#Power estimate:
mean(store_power)
#95% confidence intervals for power estimate
binom.test(table(factor(store_power,c(1,0))))$conf.int
