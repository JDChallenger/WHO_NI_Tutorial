library(ggplot2)

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



################################################################################
####################  Simulation-based power calculation   #####################
################################################################################

#Here, we will use a function to simulate the synethetic datasets, and another 
# function to carryout the non-inferiority assessment