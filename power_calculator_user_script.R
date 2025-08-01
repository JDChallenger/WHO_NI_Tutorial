
library(lme4)
#devtools::install_github("pcdjohnson/GLMMmisc")
library(GLMMmisc) #this package has a non-standard installation. Use the command above to install 
#library(ggplot2)
library(parallel)
source('power_calculator_functions_FE.R')

#####################################################
# 1. Info required

# Does the hut trial involve insecticide-treated nets, or indoor residual spraying?
# The code we've developed is designed to test for either 
# superiority or non-inferiority (of either mortality or blood-feeding inhibition)
# We can make this test between two trial arms. 
# However, in some ITN trials, both washed & unwashed nets of the same type are included.
# Therefore, we allow the option of combining data from washed & unwashed net of the 
# same type before hypothesis testing. The variable 'trial' determines the primary study
# question; it can take different values (described below)

####################### EHTs involving ITNs ####################### 
#### MEASURING MORTALITY

# trial=1. Superiority between two trial arms
# trial=2. Superiority between two ITNs (combining data from washed & unwashed ITNs of two types i.e. involving 4 trial arms)
# trial=3. Non-inferiority between two trial arms
# trial=4. Non-inferiority between two ITNs (combining data from washed & unwashed ITNs of two types i.e. involving 4 trial arms)

#### MEASURING BLOOD-FEEDING INHIBITION

# trial=5. Superiority between two trial arms
# trial=6. Superiority between two ITNs (combining data from washed & unwashed ITNs of two types i.e. involving 4 trial arms)
# trial=7. Non-inferiority between two trial arms
# trial=8. Non-inferiority between two ITNs (combining data from washed & unwashed ITNs of two types i.e. involving 4 trial arms)

#How many trial arms in total? (n_arms)

## 'npr', or 'Nights per round': How many nights should an ITN stay in a hut before the nets are rotated? (This parameter used to be called 'npr', or nights per week)

#Expected behaviour in each arm (either for mosquito mortality, or blood-feeding inhibition)
mortalities <- c(0.05, 0.2, 0.15, 0.25, 0.15, 0.30, 0.2) 
blood_feeding <- c(0.50, 0.30, 0.30, 0.25, 0.30, 0.30, 0.25)
#Note: the length of this list should equal the number of trial arms
#length(mortalities)==n_arms
#length(blood_feeding)==n_arms

#Trial may contain multiple products. We'll need to specify which arms are
# involved in the hypothesis testing (aoi) and what type of trial we're conducting
#trial = 1; aoi <- c(4,6) # We'll test whether the latter arm is superior to the former
#trial = 2; aoi <- c(4,5,6,7)# We'll test whether ITN2 is superior to ITN1. In order, these should be ITN1 (unwashed), ITN1 (washed), ITN2 (unwashed), ITN2 (washed). 
#trial = 3; aoi <- c(4,6) # We'll test whether the latter arm is non-inferior to the former
#trial = 4; aoi <- c(4,5,6,7) # We'll test whether ITN2 is non-inferior to ITN1. # In order, these should be ITN1 (unwashed), ITN1 (washed), ITN2 (unwashed), #ITN2 (washed). 

# Length of trial (rotations). 
# Described in terms of number of complete 'rotations' of the trial. 
# For example, if you have 6 huts and 6 trial arms, one rotation would take 6 weeks
# to complete (each net spends one week in a single hut).
# Needs to be at least 1 at the moment, but you can enter things like '1.25' or '1.5'

####################### End of ITN-specific parameters ##################

#How many mosquitoes per night per hut? First we specify the mean number (meanMos) 
#Should the mosquito counts be a constant value ('deterministic', det=1), or be sampled from a
#negative binomial distn (det=0) with the given mean and
#dispersion parameter (dispMos) ?
#If you are drawing from a negative-binomial distribution, this code can
#help you visualise the distribution's shape:
meanMos <- 10
dispMos <- 1
hist(rnbinom(1000, mu = meanMos, size = dispMos))

#Variability present in the assay..Here this is described by a single 
#observation-level random effect (varO), that represents overdispersion.
#This is for simplicity: when calculating power, the regression model will be 
#adjusted for day/hut/sleeper. From datasets of past EHTs, the value of varO can
#be seen to vary widely (see Supplementary Table 1 of this paper:
# https://doi.org/10.1016/j.crpvbd.2023.100115)
# The default value of 0.9 (this is the variance of the random effect) is 
# a reasonable value to use, if you're unsure. 

# Before calculating power, it is recommended to run a single trial with the 
# relevant parameter values. In this way, you can explore the simulated dataset
# in the same way you would for a real hut trial. This allows you to check that 
# the study design you have specified in the R function matches your desired design.

xc <- simulate_trial_ITN(n_arms = 6, npr = 6, 
                         responses = c(0.1,0.5,0.5,0.425,0.55,0.33),
               varO = 0.9, rotations = 1, mos_det = 0, meanMos = 20, dispMos = 1)
dim(xc)
head(xc)
table(xc$net)
table(xc$hut)
table(xc$sleeper)
hist(xc$n) # The distribution of mosquito counts
#We can calculate the total number of mosquitoes in each arm like this: 
tapply(xc$n, xc$net, mean, na.rm = F)
#... and the total number of dead mosquitoes in each arm like this:
tapply(xc$response, xc$net, mean, na.rm = F)

#The sleepers should rotate round the different huts and arms
#Note: might not be exactly equal for trials with incomplete rotations
table(xc[xc$net=='C',]$sleeper)
table(xc[xc$net=='E2',]$sleeper)
table(xc[xc$hut==2,]$sleeper)

#Check how many days this trial took (note: doesn't include rest days)
table(xc$day) 

#Note: there is also a variable called 'night' in the dataset- this just 
#denotes the day in a given round i.e. it takes a value between 1 and 'npr'


### Another function performs the hypothesis testing. We have to provide the function
# with a dataset, and tell it what to test for.
#Here's an example, using the dataset we've just generated above (xc)

#Is Arm #6 superior to Arm #4, in terms of mosquito mortality?
#This function returns a value of 1 if the null hypothesis is rejected; otherwise,
# it returns zero
hypothesis_test(trial = 1, aoi = c(4,6), dataset = xc)

#Another example, this time for non-inferiority
#Is Arm #6 non-inferior to Arm #4, in terms of mosquito mortality?
#There are two routes to select the non-inferiority margin (NIM):-
#(i) Selected a NIM with a fixed value on the log-odds scale (e.g., for mortality,
#one could select NIM =0.7) ; (ii) Select a variable NIM, based upon the 
#performance of the comparator product. The second option is preferred by WHO,
#and will be used throughout this tutorial. To use option (ii), we set 
# NIMvar=1 to use the variable NIM, and then choose the percentage 
# variation permitted using NIMpc (NIMpc = 7 means that 7% variation is permitted)
hypothesis_test(trial = 3, aoi = c(4,6), NIMvar = 1, NIMpc = 7, dataset = xc)
#This is how you would used the fixed non-inferiority margin:
hypothesis_test(trial = 3, aoi = c(4,6), NIMvar = 0, NIM = 0.7, dataset = xc)

####
# The function below simulates many trials ('nsim' specifies how many).
# For each one, hypothesis testing will be performed (defined by the values chosen for 
# 'aoi' and 'trial' above). Statistical power is given by the percentage of trials 
# for which the null hypothesis is rejected.

# This function can be slow to run, as a large number of trials need to be simulated
# It may be possible to parallelise the code, if you have a computer with multiple cores.
# Running this command (from the 'parallel' package) will tell you how many cores you can use
detectCores()
# To parallelise (i.e. speed up) this process requires different code for Mac & Windows
#computers. Here are the options:

# parallelise = 0; Don't parallise code (use this option if you are in doubt)
# parallelise = 1; Parallise code for Windows
# parallelise = 2; Parallise code for Mac OS

# The function power_calculator_ITN(...) will return (i) A printed statement of the power estimate & 
# associated 95% confidence intervals (ii) a vector of three numbers: the power estimate, the lower 
# confidence interval, and the upper confidence interval. The latter is more useful if you want to calculate
# how power varies with (e.g.) the average number of mosquitoes per hut per night (meanMos). We'll 
# show an example of this in the IRS section below.

#########################################################################################
##############                  Worked Example for ITNs                   ###############
#########################################################################################


# Let's calculate power for an illustrative example. We wish to carry out a non-inferiority trial,
# comparing a washed candidate net to a washed comparator net (this is the comparison required for WHO PQ).
# For both of these, we expect a mosquito mortality of about 20% (in this example).
# We run a 7-arm trial: in the responses vector, we put the washed
# mortality of the comparator net in the third position, and the washed candidate net in the fifth position.
# we use the argument 'aoi' to indicate the arms we wish to use for the non-inferiority assessment, and
# we set 'trial = 3' (see above). We perform one full rotation and run the trial for 6 days a week 
# (npr = 'nights per rotation' = 6). Based on recent observations (at our imagined field site),
# we expect about 10 mosquitoes per hut per night.
# We use the variable NIM outlined above, based on a maximum acceptable 
# difference in mortality of 7% 

t1 <- Sys.time()
power_calculator_ITN(parallelise = 0, trial = 3, npr = 6, rotations = 1, 
     nsim = 500, n_arms = 7, mos_det = 1, meanMos = 10, varO = .9, 
     NIMvar = 1, NIMpc = 7,
     dispMos = 1, aoi = c(3,5), responses = c(0.15,0.5,0.2,0.5,0.2,0.4,0.3))
t2 <- Sys.time()
t2 - t1
#system("say Just finished!") #On Mac, this command alerts you that the function has finished
# I don't think this works on Windows, but you could look at the beepr package and use the function beepr::beep()

# For the default parameters used here, the power is rather low (~65%, although this will vary a little,
# each time you estimate it, due to random variation. Increasing the number of simulations will
# reduce this variation). Two options you could use to increase the power here could be: 
# (i) Choose npr = 7, i.e. allow the sleepers to be rotated around all the huts 
# (which for a 7-hut trial takes 7 days) before rotating the nets; 
# (ii) Run a longer trial, by (e.g.) using 'rotations=1.5', to run an extra half rotation
# In this instance, making both these changes should result in power > 80%. 
# Alternatively, running the study in a location (or a time of year) with higher mosquito densities
# would make the study easier to power.

####################### EHTs involving IRS ############################## 

# n_arms: how many trial arms (including untreated control) 
# rep_IRS: how many huts per IRS product? 
# rep_C: How many huts per untreated control?
# Then the number of huts in the trial will be rep_C + (n_arms - 1 * rep_IRS). This'll be the same as the number of volunteers required

#How many days will the trial last? (trial_days)
#mortalities (or blood-feeding) in each arm
mortalities_IRS <- c(0.10, 0.40, 0.50, 0.47)
#blood_feeding_IRS <- c(0.50, 0.30, 0.30, 0.25)

#### MEASURING MORTALITY

# trial=9. Superiority between two trial arms
# trial=10. Non-inferiority between two trial arms

#### MEASURING BLOOD-FEEDING INHIBITION

# trial=11. Superiority between two trial arms
# trial=12. Non-inferiority between two trial arms

####################### End of IRS-specific parameters ##################

#Before calculating power, let's simulate 1 trial, to check everything looks OK
#As before, we need info on mosquito numbers.

xd <- simulate_trial_IRS(n_arms = 4, rep_IRS = 2, rep_C = 2, 
                         responses = mortalities_IRS,
                  trial_days = 75, varO = 0.8, mos_det = 0, meanMos = 4, 
                  dispMos = 0.5)
dim(xd)
head(xd)
table(xd$net) #Note: for IRS, 'net' really means 'trial arm'. I will try to update these.
table(xd$hut)
table(xd$sleeper)
xd[xd$hut==1,]

#For the hypothesis testing, we can use the same function as for ITNs
#We just need to update the value of 'trial'

#For example: is Arm #4 superior to Arm #2, in terms of mosquito mortality?
#This function returns a value of 1 if the null hypothesis is rejected; otherwise,
# it returns zero
hypothesis_test(trial = 9, aoi = c(2,4), dataset = xd)

########## 
# Now let's simulate a large number of trials, 

# parallelise = 0; Don't parallise code (use this option if you are in doubt)
# parallelise = 1; Parallise code for Windows
# parallelise = 2; Parallise code for Mac OS

# The function power_calculator_IRS(...) will return (i) A printed statement of the power estimate & 
# associated 95% confidence intervals (ii) a vector of three numbers: the power estimate, the lower 
# confidence interval, and the upper confidence interval. The latter is more useful if you want to calculate
# how power varies with (e.g.) the average number of mosquitoes per hut per night (meanMos). We'll 
# show an example of this below.

t1 <- Sys.time()
power_calculator_IRS(parallelise = 0, trial = 9, varO = 0.9, trial_days = 25,
                     rep_C = 2, rep_IRS = 4, nsim = 300, n_arms = 4, mos_det = 1, meanMos = 11, 
                 dispMos = 1.4, aoi = c(2,4), responses = mortalities_IRS)
t2 <- Sys.time()
t2 - t1
#system("say Just finished!") #On Mac, this command alerts you that the function has finished
# I don't think this works on Windows, but you could look at the beepr package and use the function beepr::beep()

#How many simulations is enough? Depends how much precision you need. 

#Finally, let's vary meanMos and see how this influences the power

av_mosquitoes <- seq(3,23,4) # values to use
store_power <- rep(0,length(av_mosquitoes)) # vector to hold power estimates
store_ci1 <- rep(0,length(av_mosquitoes)) # vector to hold lower CI
store_ci2 <- rep(0,length(av_mosquitoes)) # vector to hold upper CI

#Make a loop: each time we vary meanMos
for(j in 1:length(av_mosquitoes)){
  xx <- power_calculator_IRS(parallelise = 0, trial = 9, varO = 0.9, trial_days = 25,
                             rep_C = 2, rep_IRS = 4, nsim = 400, n_arms = 4, mos_det = 0, meanMos = av_mosquitoes[j], 
                             dispMos = 1.4, aoi = c(2,4), responses = mortalities_IRS)
  store_power[j] <- xx[1]
  store_ci1[j] <- xx[2]
  store_ci2[j] <- xx[3]
}

library(ggplot2) # We'll use ggplot2 to look at the output

#Make a data.frame
data_power <- data.frame('mosquitoes' = av_mosquitoes, 'power' = store_power,
                         'ci1' = store_ci1, 'ci2' = store_ci2)
ggplot(data_power) + geom_point(aes(x = mosquitoes, y = power)) + theme_classic() +
  xlab('Mean no. of mosquitoes') + ylab('Study Power (%)') + ylim(c(0,100)) +
  geom_errorbar(aes(x = mosquitoes, ymin = ci1, ymax = ci2)) + 
  geom_hline(yintercept = 80, color = 'orange') # orange line indicates 80% power


