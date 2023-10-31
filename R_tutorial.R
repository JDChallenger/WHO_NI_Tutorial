#if you are missing a package, you can install with the following command:
#install.packages("ggplot2")
library(ggplot2)
library(lme4)
library(cowplot)
library(tidyr)
library(stringr)
source('useful_functions_tutorial.R')

df <- read.csv("example_dataset.csv")
str(df)
head(df)

table(df$treatment)
table(df$hut)
table(df$sleeper)
table(df$day)

#This user-defined function will give the (unadjusted) mortalities, or 
#blood feeding proportion in each trial arm
summm(df, vec = df$treatment, td = 'tot_dead', tot = 'total')
summm(df, vec = df$treatment, td = 'tot_bf', tot = 'total')

#Let's save these:
tab_mortality <- summm(df, vec = df$treatment, td = 'tot_dead', tot = 'total', table = 1)
tab_bf <- summm(df, vec = df$treatment, td = 'tot_bf', tot = 'total', table = 1)

#And do the same, this time stratifying by ITN, instead of trial arm
tab_mortality_ITN <- summm(df, vec = df$ITN, td = 'tot_dead', tot = 'total', table = 1)
tab_bf_ITN <- summm(df, vec = df$ITN, td = 'tot_bf', tot = 'total', table = 1)

#These variables should be factor variables in R
df$hut <- as.factor(df$hut)
df$sleeper <- as.factor(df$sleeper)
df$day <- as.factor(df$day)
df$treatment <- as.factor(df$treatment)
df$ITN <- as.factor(df$ITN)

#########################################################################
#####          1. Mosquito mortality (unwashed ITNS)                #####
#########################################################################

#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Active_comparator_unwashed') 
levels(df$treatment)

fit1 <- 
  glm(
    cbind(tot_dead, total - tot_dead) ~
      treatment + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit1)

OR1 <- exp(coef(summary(fit1))['treatmentCandidate_unwashed',"Estimate"])
OR1_lower <- exp(coef(summary(fit1))['treatmentCandidate_unwashed',"Estimate"] - 
             1.96*coef(summary(fit1))['treatmentCandidate_unwashed','Std. Error'])
OR1_upper <- exp(coef(summary(fit1))['treatmentCandidate_unwashed',"Estimate"] + 
             1.96*coef(summary(fit1))['treatmentCandidate_unwashed','Std. Error'])

#### What should the non-inferiority margin be???
#First work out the FIC mortality- we'll use the value 
#taken directly from the data for this

#First convert percentage into a proportion
FIC_mortality1 <- tab_mortality[tab_mortality$Arm=='Active_comparator_unwashed',]$Percentage / 100
non_inf_margin1 <- ((FIC_mortality1 - 0.07) / (1- (FIC_mortality1 - 0.07))) / (FIC_mortality1 / (1- FIC_mortality1)) 

NI_1 <- plot_NI_OR(OR = OR1, ORl = OR1_lower, ORu = OR1_upper, mortality = 1,
           NIM = non_inf_margin1, precision = 3, title = 'Candidate vs. Active Comparator (unwashed)')

#Now prepare a plot of the estimated mortalities (not required for the non-inferiority assessment)
mFE(model = fit1, vec = df$treatment, intercept = 'Active_comparator_unwashed', bfi = 0, name = 'treatment')
ofs1 <- new_median_FE(model = fit1, FE = c('hut','sleeper','day'))
mk1 <- mFE(model = fit1, vec = df$treatment, intercept = 'Active_comparator_unwashed', bfi = 0, 
           name = "treatment", offset = ofs1)
summm(df, vec = df$treatment, td = 'tot_dead', tot = 'total')
mk1a <- mk1[-grep(" washed", mk1$Arm),]
mk1a$ord <- c(1,3,2,4) 

p1 <- ggplot(data = mk1a) + 
  geom_errorbarh(aes(y = ord, xmin = Lower_95pc_CI, xmax = Upper_95pc_CI), height = 0) + 
  geom_point(aes(y=ord, x=Mortality, colour = Arm), size = 3) +
  xlim(c(0,1)) + xlab('Proportion of mosquitoes blood fed') +
  theme_classic() + ylab('') + theme(axis.line.y = element_blank(),
                                     axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_discrete(breaks = c('Candidate unwashed','Active comparator unwashed',
                                  'Standard comparator unwashed','Control')) +
  theme(legend.position = c(0.8,0.3)) + labs(color = '') + # add washed status to labs??
  ggtitle('Mosquito mortality (unwashed ITNs)')
p1
cowplot::plot_grid(p1,NI_1, nrow = 1, rel_widths = c(0.6,0.4), labels = c('A','B'))
ggsave('Assessment1.pdf', height = 5.5, width = 9)

# Now also check that the candidate net is superior to the standard comparator
#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Standard_comparator_unwashed') 
levels(df$treatment)

fit1a <- 
  glm(
    cbind(tot_dead, total - tot_dead) ~
      treatment + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit1a)

coef(summary(fit1a))['treatmentCandidate_unwashed',"Pr(>|z|)"]
if(coef(summary(fit1a))['treatmentCandidate_unwashed',"Pr(>|z|)"] < 0.05 & 
   coef(summary(fit1a))['treatmentCandidate_unwashed',"Estimate"] > 0){
  print('Candidate superior to standard comparator (mosquito mortality, unwashed nets)')
}else{
  print('Candidate NOT superior to standard comparator (mosquito mortality, unwashed nets)')
}

# For the non-inferiority plot, we will now show an alternative way of presenting
# the same information. This uses the function 'variable_NIM', which 
# highlights the fact that the non-inferiority margin is variable 
# (i.e. it depends on the performance of the first-in-class product)
variable_NIM(OR = OR1, ORl = OR1_lower, ORu = OR1_upper, 
             FIC = FIC_mortality1, mortality = 1, ymin = 0.2, ymax = 0.5)
ggsave('variable_NIM1.pdf', width = 6, height = 5)

#########################################################################
######         2. Mosquito mortality (washed ITNS)                 ######
#########################################################################
#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Active_comparator_washed') 
levels(df$treatment)

fit2 <- 
  glm(
    cbind(tot_dead, total - tot_dead) ~
      treatment + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit2)

OR2 <- exp(coef(summary(fit2))['treatmentCandidate_washed',"Estimate"])
OR2_lower <- exp(coef(summary(fit2))['treatmentCandidate_washed',"Estimate"] - 
                   1.96*coef(summary(fit2))['treatmentCandidate_washed','Std. Error'])
OR2_upper <- exp(coef(summary(fit2))['treatmentCandidate_washed',"Estimate"] + 
                   1.96*coef(summary(fit2))['treatmentCandidate_washed','Std. Error'])

FIC_mortality2 <- tab_mortality[tab_mortality$Arm=='Active_comparator_washed',]$Percentage / 100
non_inf_margin2 <- ((FIC_mortality2 - 0.07) / (1- (FIC_mortality2 - 0.07))) / (FIC_mortality2 / (1- FIC_mortality2)) 

plot_NI_OR(OR = OR2, ORl = OR2_lower, ORu = OR2_upper, mortality = 1,
           NIM = non_inf_margin2, precision = 3, title = 'Candidate vs. Active Comparator (washed)')

mFE(model = fit2, vec = df$treatment, intercept = 'Active_comparator_washed', bfi = 0, name = "treatment")
ofs2 <- new_median_FE(model = fit2, FE = c('hut','sleeper','day'))
mk2 <- mFE(model = fit2, vec = df$treatment, intercept = 'Active_comparator_washed', bfi = 0, name = "treatment", offset = ofs2)
summm(df, vec = df$treatment, td = 'tot_dead', tot = 'total')
mk2a <- mk2[-grep("unwashed", mk2$Arm),]
mk2a$ord <- c(1,3,2,4) 

p2 <- ggplot(data = mk2a) + 
  geom_errorbarh(aes(y = ord, xmin = Lower_95pc_CI, xmax = Upper_95pc_CI), height = 0) + 
  geom_point(aes(y=ord, x=Mortality, colour = Arm), size = 3) +
  xlim(c(0,1)) + xlab('Proportion of mosquitoes blood fed') +
  theme_classic() + ylab('') + theme(axis.line.y = element_blank(),
                                     axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_discrete(breaks = c('Candidate washed','Active comparator washed','Standard comparator washed','Control')) +
  theme(legend.position = c(0.8,0.3)) + labs(color = '') + # add washed status to labs??
  ggtitle('Mosquito Mortality (washed ITNs)')
p2

# Now also check that the candidate net is superior to the standard comparator
#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Standard_comparator_washed') 
levels(df$treatment)

fit2a <- 
  glm(
    cbind(tot_dead, total - tot_dead) ~
      treatment + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit2a)

coef(summary(fit2a))['treatmentCandidate_washed',"Pr(>|z|)"]
if(coef(summary(fit2a))['treatmentCandidate_washed',"Pr(>|z|)"] < 0.05 & 
   coef(summary(fit2a))['treatmentCandidate_washed',"Estimate"] > 0){
  print('Candidate superior to standard comparator (mosquito mortality, washed nets)')
}else{
  print('Candidate NOT superior to standard comparator (mosquito mortality, washed nets)')
}

#########################################################################
######      3. Mosquito mortality (unwashed & washed combined)     ######
#########################################################################

df$ITN <- relevel(df$ITN, 'Active_comparator') 
levels(df$ITN)
fit3 <- 
  glm(
    cbind(tot_dead, total - tot_dead) ~
      ITN + hut + sleeper + day + wash,
    family = binomial, data = df)
summary(fit3)

OR3 <- exp(coef(summary(fit3))['ITNCandidate',"Estimate"])
OR3_lower <- exp(coef(summary(fit3))['ITNCandidate',"Estimate"] - 
                   1.96*coef(summary(fit3))['ITNCandidate','Std. Error'])
OR3_upper <- exp(coef(summary(fit3))['ITNCandidate',"Estimate"] + 
                   1.96*coef(summary(fit3))['ITNCandidate','Std. Error'])

FIC_mortality3 <- tab_mortality_ITN[tab_mortality_ITN$Arm=='Active_comparator',]$Percentage / 100
non_inf_margin3 <- ((FIC_mortality3 - 0.07) / (1- (FIC_mortality3 - 0.07))) / (FIC_mortality3 / (1- FIC_mortality3)) 

plot_NI_OR(OR = OR3, ORl = OR3_lower, ORu = OR3_upper, mortality = 1,
           NIM = non_inf_margin3, precision = 3, title = 'Candidate vs. Active Comparator (combined)')

mFE(model = fit3, vec = df$ITN, intercept = 'Active_comparator', bfi = 0, name = "ITN")
ofs3 <- new_median_FE(model = fit3, FE = c('hut','sleeper','day'))
mk3 <- mFE(model = fit3, vec = df$ITN, intercept = 'Active_comparator', bfi = 0, name = "ITN", offset = ofs3)
summm(df, vec = df$ITN, td = 'tot_dead', tot = 'total')
mk3$ord <- c(1,3,2,4) 

p3 <- ggplot(data = mk3) + 
  geom_errorbarh(aes(y = ord, xmin = Lower_95pc_CI, xmax = Upper_95pc_CI), height = 0) +
  geom_point(aes(y = ord, x = Mortality, colour = Arm), size = 3) +
  xlim(c(0,1)) + xlab('Proportion of mosquitoes blood fed') +
  theme_classic() + ylab('') + theme(axis.line.y = element_blank(),
                                     axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_discrete(breaks = c('Candidate','Active comparator','Standard comparator','Control')) +
  theme(legend.position = c(0.8,0.3)) + labs(color = '') + # add washed status to labs??
  ggtitle('Mosquito Mortality (combined analysis)')
p3

# Now also check that the candidate net is superior to the standard comparator
#Change baseline treatment category
df$ITN <- relevel(df$ITN, 'Standard_comparator') 
levels(df$ITN)

fit3a <- 
  glm(
    cbind(tot_dead, total - tot_dead) ~
      ITN + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit3a)

coef(summary(fit3a))['ITNCandidate',"Pr(>|z|)"]
if(coef(summary(fit3a))['ITNCandidate',"Pr(>|z|)"] < 0.05 & 
   coef(summary(fit3a))['ITNCandidate',"Estimate"] > 0){
  print('Candidate superior to standard comparator (mosquito mortality, combined analysis)')
}else{
  print('Candidate NOT superior to standard comparator (mosquito mortality, combined analysis)')
}

#########################################################################
######            4. Blood Feeding (unwashed ITNS)                 ######
#########################################################################
#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Active_comparator_unwashed') 
levels(df$treatment)

fit4 <- 
  glm(
    cbind(tot_bf, total - tot_bf) ~
      treatment + hut + sleeper + day,
    family = binomial, data = df)
summary(fit4)

OR4 <- exp(coef(summary(fit4))['treatmentCandidate_unwashed',"Estimate"])
OR4_lower <- exp(coef(summary(fit4))['treatmentCandidate_unwashed',"Estimate"] - 
                   1.96*coef(summary(fit4))['treatmentCandidate_unwashed','Std. Error'])
OR4_upper <- exp(coef(summary(fit4))['treatmentCandidate_unwashed',"Estimate"] + 
                   1.96*coef(summary(fit4))['treatmentCandidate_unwashed','Std. Error'])

FIC_bf4 <- tab_bf[tab_bf$Arm=='Active_comparator_unwashed',]$Percentage / 100
non_inf_margin4 <- ((FIC_bf4 + 0.07) / (1- (FIC_bf4 + 0.07))) / (FIC_bf4 / (1- FIC_bf4)) 

plot_NI_OR(OR = OR4, ORl = OR4_lower, ORu = OR4_upper, mortality = 0,
           NIM = non_inf_margin4, precision = 3, title = 'Candidate vs. Active Comparator (unwashed)')

mFE(model = fit4, vec = df$treatment, intercept = 'Active_comparator_unwashed', bfi = 1, name = "treatment")
ofs4 <- new_median_FE(model = fit4, FE = c('hut','sleeper','day'))
mk4 <- mFE(model = fit4, vec = df$treatment, intercept = 'Active_comparator_unwashed', bfi = 1, name = "treatment", offset = ofs4)
summm(df, vec = df$treatment, td = 'tot_bf', tot = 'total')
mk4a <- mk4[-grep(" washed", mk4$Arm),]
mk4a$ord <- c(1,3,2,4) 

p4 <- ggplot(data = mk4a) + 
  geom_errorbarh(aes(y = ord, xmin = Lower_95pc_CI, xmax = Upper_95pc_CI), height = 0) + 
  geom_point(aes(y = ord, x = Blood.Feeding, colour = Arm), size = 3) +
  xlim(c(0,1)) + xlab('Proportion of mosquitoes blood fed') +
  theme_classic() + ylab('') + theme(axis.line.y = element_blank(),
                                     axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_discrete(breaks = c('Candidate unwashed','Active comparator unwashed',
                                  'Standard comparator unwashed','Control')) +
  theme(legend.position = c(0.8,0.3)) + labs(color = '') + # add washed status to labs??
  ggtitle('Blood Feeding (unwashed ITNs)')
p4

# With Assessment 1 & assessment4, you can build a useful summary of the whole trial, using this function:
summary_output <- tidy_blf_FE(data = df, model_fit = fit1, name = "treatment", vec = df$treatment,
            model_fit_blf = fit4, offset = c(ofs1,ofs4),
            intercept = 'Active_comparator_unwashed', first_cat = 'Control')
summary_output
#save this? E.g.
write.csv(summary_output, 'trial_summary.csv', row.names = F, col.names = T)

# Now also check that the candidate net is superior to the standard comparator
#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Standard_comparator_unwashed') 
levels(df$treatment)

fit4a <- 
  glm(
    cbind(tot_bf, total - tot_bf) ~
      treatment + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit4a)

coef(summary(fit4a))['treatmentCandidate_unwashed',"Pr(>|z|)"]
if(coef(summary(fit4a))['treatmentCandidate_unwashed',"Pr(>|z|)"] < 0.05 & 
   coef(summary(fit4a))['treatmentCandidate_unwashed',"Estimate"] < 0){
  print('Candidate superior to standard comparator (blood feeding, unwashed nets)')
}else{
  print('Candidate NOT superior to standard comparator (blood feeding, unwashed nets)')
}

# For the non-inferiority plot, we will now show an alternative way of presenting
# the same information. This uses the function 'variable_NIM', which 
# highlights the fact that the non-inferiority margin is variable 
# (i.e. it depends on the performance of the first-in-class product)
variable_NIM(OR = OR4, ORl = OR4_lower, ORu = OR4_upper, 
             FIC = FIC_bf4, mortality = 0, ymin = 0.09, ymax = 0.5, xmax = 2)

#########################################################################
######              5. Blood Feeding (washed ITNS)                 ######
#########################################################################
df$treatment <- relevel(df$treatment, 'Active_comparator_washed') 
levels(df$treatment)

fit5 <- 
  glm(
    cbind(tot_bf, total - tot_bf) ~
      treatment + hut + sleeper + day,
    family = binomial, data = df)
summary(fit5)

OR5 <- exp(coef(summary(fit5))['treatmentCandidate_washed',"Estimate"])
OR5_lower <- exp(coef(summary(fit5))['treatmentCandidate_washed',"Estimate"] - 
                   1.96*coef(summary(fit5))['treatmentCandidate_washed','Std. Error'])
OR5_upper <- exp(coef(summary(fit5))['treatmentCandidate_washed',"Estimate"] + 
                   1.96*coef(summary(fit5))['treatmentCandidate_washed','Std. Error'])

FIC_bf5 <- tab_bf[tab_bf$Arm=='Active_comparator_washed',]$Percentage / 100
non_inf_margin5 <- ((FIC_bf5 + 0.07) / (1- (FIC_bf5 + 0.07))) / (FIC_bf5 / (1- FIC_bf5)) 

plot_NI_OR(OR = OR5, ORl = OR5_lower, ORu = OR5_upper, mortality = 0,
           NIM = non_inf_margin5, precision = 3, title = 'Candidate vs. Active Comparator (washed)')

mFE(model = fit5, vec = df$treatment, intercept = 'Active_comparator_washed', bfi = 1, name = "treatment")
ofs5 <- new_median_FE(model = fit5, FE = c('hut','sleeper','day'))
mk5 <- mFE(model = fit5, vec = df$treatment, intercept = 'Active_comparator_washed', bfi = 1, name = "treatment", offset = ofs5)
summm(df, vec = df$treatment, td = 'tot_bf', tot = 'total')
mk5a <- mk5[-grep("unwashed", mk5$Arm),]
mk5a$ord <- c(1,3,2,4) 

p5 <- ggplot(data = mk5a) + 
  geom_errorbarh(aes(y = ord, xmin = Lower_95pc_CI, xmax = Upper_95pc_CI), height = 0) + 
  geom_point(aes(y = ord, x = Blood.Feeding, colour = Arm), size = 3) +
  xlim(c(0,1)) + xlab('Proportion of mosquitoes blood fed') +
  theme_classic() + ylab('') + theme(axis.line.y = element_blank(),
                                     axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_discrete(breaks = c('Candidate washed','Active comparator washed',
                                  'Standard comparator washed','Control')) +
  theme(legend.position = c(0.8,0.3)) + labs(color = '') + # add washed status to labs??
  ggtitle('Blood Feeding (washed ITNs)')
p5

# Now also check that the candidate net is superior to the standard comparator
#Change baseline treatment category
df$treatment <- relevel(df$treatment, 'Standard_comparator_washed') 
levels(df$treatment)

fit5a <- 
  glm(
    cbind(tot_bf, total - tot_bf) ~
      treatment + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit5a)

coef(summary(fit5a))['treatmentCandidate_washed',"Pr(>|z|)"]
if(coef(summary(fit5a))['treatmentCandidate_washed',"Pr(>|z|)"] < 0.05 & 
   coef(summary(fit5a))['treatmentCandidate_washed',"Estimate"] < 0){
  print('Candidate superior to standard comparator (blood feeding, washed nets)')
}else{
  print('Candidate NOT superior to standard comparator (blood feeding, washed nets)')
}


#########################################################################
######          6. Blood Feeding (unwashed & washed combined)      ######
#########################################################################

levels(df$ITN)
fit6 <- 
  glm(
    cbind(tot_bf, total - tot_bf) ~
      ITN + hut + sleeper + day + wash,
    family = binomial, data = df)
summary(fit6)

OR6 <- exp(coef(summary(fit6))['ITNCandidate',"Estimate"])
OR6_lower <- exp(coef(summary(fit6))['ITNCandidate',"Estimate"] - 
                   1.96*coef(summary(fit6))['ITNCandidate','Std. Error'])
OR6_upper <- exp(coef(summary(fit6))['ITNCandidate',"Estimate"] + 
                   1.96*coef(summary(fit6))['ITNCandidate','Std. Error'])

FIC_bf6 <- tab_bf_ITN[tab_bf_ITN$Arm=='Active_comparator',]$Percentage / 100
non_inf_margin6 <- ((FIC_bf6 + 0.07) / (1- (FIC_bf6 + 0.07))) / (FIC_bf6 / (1- FIC_bf6)) 

plot_NI_OR(OR = OR6, ORl = OR6_lower, ORu = OR6_upper, mortality = 0,
           NIM = non_inf_margin6, precision = 3, title = 'Candidate vs. Active Comparator (combined)')

mFE(model = fit6, vec = df$ITN, intercept = 'Active_comparator', bfi = 1, name = "ITN")
ofs6 <- new_median_FE(model = fit6, FE = c('hut','sleeper','day'))
mk6 <- mFE(model = fit6, vec = df$ITN, intercept = 'Active_comparator', bfi = 1, name = "ITN", offset = ofs6)
summm(df, vec = df$treatment, td = 'tot_bf', tot = 'total')
mk6$ord <- c(1,3,2,4) 

p6 <- ggplot(data = mk6) + 
  geom_errorbarh(aes(y = ord, xmin = Lower_95pc_CI, xmax = Upper_95pc_CI), height = 0) + 
  geom_point(aes(y = ord, x = Blood.Feeding, colour = Arm), size = 3) +
  xlim(c(0,1)) + xlab('Proportion of mosquitoes blood fed') +
  theme_classic() + ylab('') + theme(axis.line.y = element_blank(),
                                     axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_discrete(breaks = c('Candidate','Active comparator','Standard comparator','Control')) +
  theme(legend.position = c(0.8,0.3)) + labs(color = '') + # add washed status to labs??
  ggtitle('Blood Feeding (combined analysis)')
p6

# Now also check that the candidate net is superior to the standard comparator
#Change baseline treatment category
df$ITN <- relevel(df$ITN, 'Standard_comparator') 
levels(df$ITN)

fit6a <- 
  glm(
    cbind(tot_bf, total - tot_bf) ~
      ITN + hut + sleeper + day,# + wash,
    family = binomial, data = df)
summary(fit6a)

coef(summary(fit6a))['ITNCandidate',"Pr(>|z|)"]
if(coef(summary(fit6a))['ITNCandidate',"Pr(>|z|)"] < 0.05 & 
   coef(summary(fit6a))['ITNCandidate',"Estimate"] < 0){
  print('Candidate superior to standard comparator (blood feeding, combined analysis)')
}else{
  print('Candidate NOT superior to standard comparator (blood feeding, combined analysis)')
}

