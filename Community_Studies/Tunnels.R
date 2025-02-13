library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one
tunnel_data <- tunnel_sim2(verbose = T)
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

#Here, we will use a function to simulate the synethetic datasets, and another 
# function to carryout the non-inferiority assessment

nsim <- 1000
store_power <- rep(NA, nsim) # a container for the outcome of each simulated study

for(i in 1:nsim){
  tunnel_data <- tunnel_sim2(sigma_net = 0.9, n_nets = 30, #reps = 1,
                      mort = 0.3, n_mosq = 45)
  store_power[i] <- tunnel_NIM(data = tunnel_data, verbose = F)
}

#Power estimate:
mean(store_power)
#95% confidence intervals for power estimate
binom.test(table(factor(store_power,c(1,0))))$conf.int
