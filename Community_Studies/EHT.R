library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one

#Update function name
mosdata <- EHT_sim2(n_arms = 9, meanMos = 10)
str(mosdata)
head(mosdata)

table(mosdata$net)
table(mosdata$day)

levels(mosdata$net)
mosdata$net <- relevel(mosdata$net, 'E2')
levels(mosdata$net) # view levels again

#Fit regression model. day, hut & sleeper included as fixed effects

fit <- glm(
  cbind(tot_dead, n - tot_dead) ~ net + day + sleeper + hut, 
  family = binomial, data = mosdata
)
summary(fit)

OR <- exp(coef(summary(fit))['netE5',"Estimate"])
OR_lower <- exp(coef(summary(fit))['netE5',"Estimate"] - 
                   1.96*coef(summary(fit))['netE5','Std. Error'])
OR_upper <- exp(coef(summary(fit))['netE5',"Estimate"] + 
                   1.96*coef(summary(fit))['netE5','Std. Error'])

FIC_mortality <- sum(mosdata[mosdata$net=='E2',]$tot_dead)/sum(mosdata[mosdata$net=='E2',]$n)
non_inf_margin <- ((FIC_mortality - 0.07) / (1- (FIC_mortality - 0.07))) / (FIC_mortality / (1- FIC_mortality)) 

OR_lower > non_inf_margin # Is the non-inferiority criterion satisfied? This should return 'TRUE', if it is


################################################################################
####################  Simulation-based power calculation   #####################
################################################################################

#Here, we will use a function to simulate the synthetic datasets, and another 
# function to carryout the non-inferiority assessment

#full details of the arguments of each function

nsim <- 1000 # the number of simulations to carry out
store_power <- rep(NA, nsim) # a container for the outcome of each simulated study
for(i in 1:nsim){
  mosdata <- EHT_sim(n_arms = 9, npr = 9, mos_det = 0,
                     meanMos = 10, rotations = 1,
                     mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25),
                     sigma_net = 0.5)
  store_power[i] <- EHT_NIM(dataset = mosdata, verbose = F)
}

mean(store_power)
binom.test(table(factor(store_power,c(1,0))))$conf.int


