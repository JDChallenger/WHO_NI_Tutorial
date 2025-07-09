# code adapted from '/...2023/RT_auto_day40_EMM.R'

library(ggplot2)
library(lme4)
library(readxl)
library(emmeans)

#We need a dataset to look at. I have uploaded 12 simulated datasets for users
# to practise with. You can replace these with your own data, but ensure that the 
# variable names are consistent with the simulated datasets.

# In the Excel sheet, the datasets are saved in separate sheets ('sim1', 'sim2', 'sim3', etc.)

dx <-  readxl::read_xlsx('Regeneration_Time/twelve_datasets_bioassay.xlsx', sheet = 'sim2')
head(dx)
#For plotting, it is more useful to have 'day' as a numeric variable
# For the GLM, it needs to be a factor variable. So we make a 2nd version
dx$day <- as.numeric(dx$day) # numeric 'day'
dx$day_f <- as.factor(dx$day) # factor variable for 'day'
table(dx$day)
days <- unique(dx$day)

#This function calculates the raw (unadjusted) mortality estimates from the dataset
average_mortality_ci <- function(data = df, days_v = days){
  av_mort <- rep(0,length(days_v))
  av_mort_c1 <- rep(0,length(days_v))
  av_mort_c2 <- rep(0,length(days_v))
  for(i in 1:length(days_v)){
    denom <- sum(data[data$day==days_v[i],]$total)
    num <- sum(data[data$day==days_v[i],]$tot_dead)
    av_mort[i] <- num/denom 
    xr <- c(rep(1,num),rep(0,denom - num))
    av_mort_c1[i] <- binom.test(table(factor(xr,c(1,0))))$conf.int[1]
    av_mort_c2[i] <- binom.test(table(factor(xr,c(1,0))))$conf.int[2]
  }
  dfm <- data.frame('day' = days_v, 'mort' = av_mort,
                    'mort_c1' = av_mort_c1,'mort_c2' = av_mort_c2)
  return(dfm)
}

#We can then estimate these, and plot:
mm <- average_mortality_ci(data = dx, days_v = days)
ggplot(mm) + 
  geom_point(aes(x=day, y=mort)) + theme_classic() + 
  geom_errorbar(aes(x=day, ymin = mort_c1, ymax = mort_c2)) + 
  xlab('Day') + ylab('Mortality') + ylim(c(0,1))

#However, our main interest is the predictions from the GLM
# remember to use the factor variable for day
fit <- glm(
  cbind(tot_dead, total - tot_dead) ~
    day_f + batch, 
  family = binomial, data = dx)
summary(fit)
#use the emmeans package to generate predicted mortalities, averaging over batch effects
pred <- emmeans(fit, "day_f", type = "response")
#convert to data frame
pred2 <- as.data.frame(pred)
#Add a numeric variable for day to this data frame (careful converting from a factor!)
pred2$day <- as.numeric(as.character(pred2$day_f))

# extract a vector of the model predictions
xv2 <- pred2$prob
ld <- length(xv2) - 1

# Now we loop over the model-estimated mortalities, and select the first one
# that fits the criteria for being the regeneration time (no later mortality estimate
# with a mortality more than 5% higher)

count <- 0
for(j in 1:ld){
  bb <- all(xv2[j] + 0.05 > xv2[(j+1):(ld+1)])
  if(bb==T){
    rt <- days[j]
    print(paste0('Regeneration time is: ',rt))
    count <- 1
    rtG3 <- j
    break # finish the for loop when the RT has been identified
  }
  if(j==ld & count==0){
    print('Routine has identified the regeneration time as being the final measured value. Inspect data & model predictions, to check this makes sense')
    rt <- days[j+1]
    print(paste0('Regeneration time is: ',rt))
    rtG3 <- j + 1 #i.e. the last value
  }
}

# Plot: model-adjusted mortality; unadjusted mortality; and identify the RT via our method
# Open black circles: unadjusted mortality estimates (error bars are the 95% CIs)
# Closed purple circles: model-adjusted mortality estimates
# Solid pink horizontal line: model-adjusted estimate of mortality at the time point we establish as the R.T.
# Dashed pink horizontal line: mortality 5% higher than the established R.T.
# The purple box indicates the estimated R.T. 
ggplot() + geom_hline(yintercept = 
      c(xv2[rtG3],xv2[rtG3] + 0.05), color = 'maroon2', linetype = c('solid','dashed')) + 
  geom_point(data = mm, aes(x=day, y=mort, color = 'a'), shape = 1) + theme_classic() + 
  geom_errorbar(data = mm, aes(x=day, ymin = mort_c1, ymax = mort_c2), width = .5) + 
  geom_point(data = pred2, aes(x = day, y = prob, color = 'b')) + 
  geom_errorbar(data = pred2, aes(x=day, ymin = asymp.LCL, ymax = asymp.UCL), 
                color = 'purple') +
  xlab('Day') + ylab('Mortality') + ylim(c(0,max(1,xv2[rtG3] + 0.05))) + 
  annotate('rect', xmin = days[rtG3]-1.5, xmax = days[rtG3] + 1.5,
           ymin = xv2[rtG3] - 0.05, ymax = xv2[rtG3] + 0.05,
           color = 'purple', fill = 'purple', alpha = .15) + 
  scale_color_manual(values = c('black','purple'), 
                     labels = c('Unadjusted estimate','Adjusted estimate'),
                     name = 'Mosquito mortality') + theme(legend.position = c(0.78,0.2),
                legend.background = element_blank(),
                legend.box.background = element_rect(colour = "black")) + 
  ggtitle('Estimation of regeneration time')




