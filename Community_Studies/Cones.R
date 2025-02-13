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


################################################################################
####################  Simulation-based power calculation   #####################
################################################################################

#Here, we will use a function to simulate the synethetic datasets, and another 
# function to carryout the non-inferiority assessment