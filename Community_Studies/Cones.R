library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one
cone_data <- cone_sim2()
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
 cone_data <- cone_sim2(cone_mort = 0.4, npos = 4, rep = 4, nday = 20, verbose = F) # check the number of replicates required 
 store_power[i] <- cone_NIM(dataset = cone_data, verbose = F)
}

#Power estimate:
mean(store_power)
#95% confidence intervals for power estimate
binom.test(table(factor(store_power,c(1,0))))$conf.int