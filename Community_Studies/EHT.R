library(ggplot2)
source('Community_Studies/useful_functions_CS.R')

################################################################################
#########################   Analysing a dataset   ##############################
################################################################################

#Here, we will simulate a dataset, and then analyse it.
# The simulated dataset could be replaced with a real one

#Update function name
mosdata <- EHT_sim(n_arms = 9, meanMos = 10)
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

################################################################################
#Now we re-create Figure 3.1 from the tutorial
# This will be very slow to run!
# This figure requires three different study durations: we'll simulate these separately

nsim <- 1000
########## 9-arm trial; ONE ROTATION

dgm <- data.frame()
num_mosq <- seq(5,20,5)
#morts <- seq(0.1,0.5,0.1) vary mortality too??

for(k in 1:length(num_mosq)){ # mosquito numbers
  store_sigma_net <- rep(0,9) # change these names? 
  store_sigma_net_ci1 <- rep(0,9)
  store_sigma_net_ci2 <- rep(0,9)
  for(j in 1:length(store_sigma_net)){
    
    store_power <- rep(NA, nsim)
    for(i in 1:nsim){
      
      mosdata <- EHT_sim(n_arms = 9, npr = 9, mos_det = 0, meanMos = num_mosq[k], rotations = 1, 
                          mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25),
                          sigma_net = 0.1*j)
      
      store_power[i] <- EHT_NIM(dataset = mosdata, verbose = F)
      
    }
    print(paste0('J',j))
    
    store_sigma_net[j] <- mean(store_power)
    store_sigma_net_ci1[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[1]
    store_sigma_net_ci2[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[2]
  }
  aux <- data.frame('sigma_net' = 0.1*(1:9), 'power' = store_sigma_net,
                    'power_ci1' = store_sigma_net_ci1, 
                    'power_ci2' = store_sigma_net_ci2, 
                    # 'mort' = rep(morts[k], length(store_sigma_net_ci1)),
                    'num_mosq' = rep(num_mosq[k], length(store_sigma_net_ci1)) )
  dgm <- rbind(dgm, aux)
  print(paste0('k equal: ', k))
}  
date()


########## 9-arm trial; 1.5 ROTATIONS
dgma <- data.frame()
num_mosq <- seq(5,20,5)
#morts <- seq(0.1,0.5,0.1) vary mortality too??

for(k in 1:length(num_mosq)){ # mosquito numbers
  store_sigma_net <- rep(0,9) # change these names? 
  store_sigma_net_ci1 <- rep(0,9)
  store_sigma_net_ci2 <- rep(0,9)
  for(j in 1:length(store_sigma_net)){
    
    store_power <- rep(NA, nsim)
    for(i in 1:nsim){
      
      mosdata <- EHT_sim(n_arms = 9, npr = 9, mos_det = 0, meanMos = num_mosq[k], rotations = 1.5, 
                          mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25),
                          sigma_net = 0.1*j)
      
      store_power[i] <- EHT_NIM(dataset = mosdata, verbose = F)
      
    }
    #print(paste0('J',j))
    
    store_sigma_net[j] <- mean(store_power)
    store_sigma_net_ci1[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[1]
    store_sigma_net_ci2[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[2]
  }
  aux <- data.frame('sigma_net' = 0.1*(1:9), 'power' = store_sigma_net,
                    'power_ci1' = store_sigma_net_ci1, 
                    'power_ci2' = store_sigma_net_ci2, 
                    # 'mort' = rep(morts[k], length(store_sigma_net_ci1)),
                    'num_mosq' = rep(num_mosq[k], length(store_sigma_net_ci1)) )
  dgma <- rbind(dgma, aux)
  print(paste0('k equal: ', k))
}  
date()


########## 12-arm trial; ONE ROTATION (field-aged nets in two arms) 

dgmb <- data.frame()
num_mosq <- seq(5,20,5)

for(k in 1:length(num_mosq)){ # mosquito numbers
  store_sigma_net <- rep(0,9) # change these names? 
  store_sigma_net_ci1 <- rep(0,9)
  store_sigma_net_ci2 <- rep(0,9)
  for(j in 1:length(store_sigma_net)){
    
    store_power <- rep(NA, nsim)
    for(i in 1:nsim){
      
      mosdata <- EHT_sim(n_arms = 12, npr = 12, mos_det = 0, meanMos = 20, rotations = 1, 
                          mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25,0.5,0.35,0.25),
                          sigma_net = 0.1*j)
      
      #Re label arms 10-12
      mosdata[mosdata$net=='E9',]$net <- 'E3'
      mosdata[mosdata$net=='E10',]$net <- 'E4'
      mosdata[mosdata$net=='E11',]$net <- 'E5'
      
      store_power[i] <- EHT_NIM(dataset = mosdata, verbose = F)
      
    }
    #print(paste0('J',j))
    
    store_sigma_net[j] <- mean(store_power)
    store_sigma_net_ci1[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[1]
    store_sigma_net_ci2[j] <- binom.test(table(factor(store_power,c(1,0))))$conf.int[2]
  }
  aux <- data.frame('sigma_net' = 0.1*(1:9), 'power' = store_sigma_net,
                    'power_ci1' = store_sigma_net_ci1, 
                    'power_ci2' = store_sigma_net_ci2, 
                    # 'mort' = rep(morts[k], length(store_sigma_net_ci1)),
                    'num_mosq' = rep(num_mosq[k], length(store_sigma_net_ci1)) )
  dgmb <- rbind(dgmb, aux)
  print(paste0('k equal: ', k))
}  
date()

ggplot() + geom_line(data = dgm, aes(x = sigma_net, y = 100*power, color = factor(num_mosq),
                                     linetype = 'a')) + 
  geom_line(data = dgma, aes(x = sigma_net, y = 100*power, color = factor(num_mosq),
                             linetype = 'b')) + 
  geom_line(data = dgmb, aes(x = sigma_net, y = 100*power, color = factor(num_mosq),
                             linetype = 'c')) + 
  theme_classic() + ylab('Statistical Power (%)') + facet_wrap(~num_mosq) + 
  geom_hline(yintercept = 80, color = 'grey') + 
  labs(color = 'Mean no.\nmosquitoes\nper night') + 
  scale_linetype_manual(values = c('solid','dashed','dotted'), name = 'Study duration',
                        labels = c('9-arm; npr=9;\n1 rotation','9-arm; npr=9;\n1.5 rotations','12-arm; npr=12;\n1rotation')) +
  xlab('Between-net heterogeneity (s.d.)') + ylim(c(0,100)) +
  theme(legend.key.size = unit(1.75, 'lines')) + theme(legend.spacing.y = unit(1.0, 'cm'))

