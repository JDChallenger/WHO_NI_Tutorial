library(ggplot2)
library(lme4)
library(readxl)
library(emmeans)

#We need a dataset to look at. I have uploaded 12 simulated datasets for users
# to practise with. You can replace these with your own data, but ensure that the 
# variable names are consistent with the simulated datasets.

#days <- c(1,2,3,5,7,10,15,20,30,40) #do we need day0? # remove?
dat1 <-  readxl::read_xlsx('Regeneration_Time/twelve_datasets_chem.xlsx', sheet = 'sim1')
str(dat1)

#For plotting, it is more useful to have 'day' as a numeric variable
# For the GLM, it needs to be a factor variable. So we make a 2nd version
dat1$day <- as.numeric(dat1$day) # numeric 'day'
dat1$day_f <- as.factor(dat1$day) # factor variable for 'day'

days <- unique(dat1$day)

#function to obtain (unadjusted) average concentration measured on each day in the dataset
av_chem <- function(data = df, days_v = days){
  av_ch <- rep(0,length(days_v))

  for(i in 1:length(days_v)){
    av_ch[i] <- mean(data[data$day==days_v[i],]$conc)
  }
  dfm <- data.frame('day' = days_v, 'av_ch' = av_ch)
  return(dfm)
}
mc <- av_chem(data = dat1)
#Then can plot
ggplot(mc) + geom_point(aes(x = day, y = av_ch)) + theme_classic() + 
  ylab('Surface concentration [g/kg]') + xlab('Day')

#Now fit linear model to data, with fixed-effects for day (as a factor variable) & batch
fitC <- lm(conc ~ day_f + batch, data = dat1)
summary(fitC)
rc <- emmeans(fitC, "day_f")
pred_c <- as.data.frame(rc)
#Add a numeric variable for day to this data frame (careful converting from a factor!)
pred_c$day <- as.numeric(as.character(pred_c$day_f))

#this variable gives the predicted chemical concentrations
mc2 <- pred_c$emmean

# Now we loop over the model-estimated chemical concentrations, and select the first one
# that fits the criteria for being the regeneration time (no later concentration estimate
# more than 10% higher)
ld <- length(mc2) - 1
count <- 0
rt <- 0
for(j in 1:ld){
  bb <- all(mc2[j]*1.1 > mc2[(j+1):(ld+1)])
  if(bb==T){
    print(paste0('Regeneration time is: ', days[j]))
    rt <- j
    rt_c <- days[j]
    count <- 1
    break # end for loop
  }
  if(j==ld & count==0){
    print('Routine has identified the regeneration time as being the final measured value. Inspect data, to check this makes sense')
    rt <- j
    rt_c <- days[j]
  }
}

# Plot: model-adjusted chemical conc; unadjusted concentrations; and identify the RT via our method
# Open black circles: unadjusted concentration estimates
# Closed orange circles: model-adjusted mortality estimates (with 95% CIs)
# Solid pink horizontal line: model-adjusted estimate of mortality at the time point we establish as the R.T.
# Dashed pink horizontal line: concentration 10% higher than measured at the established R.T.
# The orange box indicates the estimated R.T. 

ggplot() + geom_hline(yintercept = c(mc2[rt],mc2[rt]*1.1), color = 'maroon2',
                        linetype = c('solid','dashed')) + 
  geom_point(data = mc, aes(x = day, y = av_ch, color = 'a'), shape = 1) +
  geom_point(data = pred_c, aes(x = day, y = emmean, color = 'b')) + 
  theme_classic() + ylim(c(0.66*min(mc2),NA)) + 
  geom_errorbar(data = pred_c, aes(x = day, ymin = lower.CL, ymax = upper.CL), color = 'orange') + 
  scale_color_manual(values = c('black','orange'), 
                       labels = c('Unadjusted estimate','Adjusted estimate'),
                       name = 'Mosquito mortality') + 
    theme(legend.position = c(0.78,0.2), legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black")) + 
  annotate('rect', xmin = days[rt] - 0.75, xmax = days[rt] + 0.75,
           ymin = mc2[rt]*0.95, ymax = mc2[rt]*1.05,
           color = 'orange', fill = 'orange', alpha = .15) +
  ylab('Surface concentration [g/kg]') + xlab('Day') + 
  ggtitle('Estimation of regeneration time - Chemical analysis')


