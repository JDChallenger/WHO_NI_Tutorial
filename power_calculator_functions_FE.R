library(lme4)
#devtools::install_github("pcdjohnson/GLMMmisc")
library(GLMMmisc)
#library(ggplot2)
library(parallel)

#user-defined function, to convert from log-odds scale to probability scale
InvLogit <- function(X){
  exp(X)/(1+exp(X))
}
#user-defined function
repeat_fn <- function(sttr, repp, repC){
  bray <- rep(sttr[1],repC)
  l <- length(sttr) - 1
  for(i in 1:l){
    bray <- c(bray,rep(sttr[i+1],repp))
  }
  return(bray)
}

simulate_trial_ITN <- function(n_arms, npw, mos_det = 0, meanMos, dispMos = 1.5,
                           rotations = 1, responses, varO = 0.9){
  #Check length(responses) == n_arms
  if(n_arms != length(responses)){
    print('Operation was not executed. Check the number of arms, and corresponding vector of mosquito mortalities')
    return(-9)
  }
  if(npw > n_arms){
    print('Operation was not executed. Rethink trial design (npw should be less than or equal to n_arms)')
    return(-9)
  }
  if(rotations < 1){ #At the moment, we need at least 1 rotation
    rotations <- 1
    print('Trial simulated for rotations = 1')
  }
  if(varO < 0 ){ #Reject negative variance for the random effect
    print('The variance of a random effect must be positive')
    return(-9)
  }
  if(meanMos < 0 | dispMos < 0 ){ #Mean & dispersion parameters must both be positive
    print('Parameters for mosquito counts should be positive')
    return(-9)
  }
  
  n_volunteer <- n_arms
  n_hut <- n_arms
  
  #Make a vector of trial arms
  aux <- c('C',paste0('E',seq(1,n_arms-1)))
  aux2 <- aux
  for(k in 1:(n_arms-1)){
    aux2 <- rbind(aux2,c(aux[(n_arms-k+1):n_arms],aux[1:(n_arms-k)]))
  }
  
  mosdata <-
    expand.grid(
      hut = factor(1:ncol(aux2)),
      week = 1:nrow(aux2),
      night = 1:npw
    )
  mosdata <- mosdata[order(mosdata$hut, mosdata$week, mosdata$night),]
  
  mosdata$net <- NA
  count <- 1
  for(i in 1:n_arms){
    for(j in 1:n_arms){
      mosdata$net[((npw*(count-1))+1):(npw*count)]  <- rep(aux2[i,j],npw)
      count <- count + 1
    }
  }
  #table(mosdata$net, useNA = 'a')
  mosdata$sleeper <- NA
  aux3 <- sample(1:n_arms) #was '7'- a mistake?
  
  for(i in 1:n_arms){ # number of weeks for one rotation
    aux4 <- c(aux3[c(i:n_arms)],aux3[seq_len(i-1)])
    #print(aux4)
    for(j in 1:npw){ #nights per week
      mosdata[mosdata$week==i & mosdata$night==j,]$sleeper <- 
        c(aux4[c(j:n_arms)],aux4[seq_len(j-1)])
    }
  }
  
  #Do we have rotation > 1? What if rotation > 2??
  if(rotations > 1){
    if(rotations <= 2){
      mosdata$day <- npw*(mosdata$week - 1) + mosdata$night
      extra <- round((rotations - 1)*dim(mosdata)[1] / n_arms)
      mosdataX <- mosdata[mosdata$day <= extra,]
      mosdataX$day <- mosdataX$day + max(mosdata$day)
      mosdataX$week <- mosdataX$week + max(mosdata$week)
      #Now add it on to the existing trial design
      mosdata <- rbind(mosdata,mosdataX)
    }else{
      mosdata$day <- npw*(mosdata$week - 1) + mosdata$night
      intt <- as.integer(rotations - 1)
      #print(paste0('intt: ',intt))
      dec <- rotations - intt - 1
      mosdataY <- mosdata
      mosdataZ <- mosdata
      for(j in 1:intt){
        mosdataY$day <- mosdataZ$day + max(mosdataZ$day)*j
        mosdataY$week <- mosdataZ$week + max(mosdataZ$week)*j
        mosdata <- rbind(mosdata, mosdataY)
      }
      #print(paste0('dim is ',dim(mosdata)[1]))
      extra <- round(dec*dim(mosdataZ)[1] / n_arms)
      #print(paste0('extra: ',extra))
      mosdataX <- mosdataZ[mosdataZ$day <= extra,]
      mosdataX$day <- mosdataX$day + max(mosdata$day)
      mosdataX$week <- mosdataX$week + max(mosdata$week)
      #print(paste0('dimX is ',dim(mosdataX)[1]))
      mosdata <- rbind(mosdata,mosdataX)
    }
  }else{
    mosdata$day <- npw*(mosdata$week - 1) + mosdata$night
  }
  
  #Now make a unique identifier for each data point
  mosdata$observation <- factor(formatC(1:nrow(mosdata), flag="0", width=3))
  mosdata$sleeper <- factor(mosdata$sleeper)
  mosdata$net <- factor(mosdata$net)
  
  mosdata$n <- NA
  if(mos_det==1){
    mosdata$n <- meanMos
  }else{
    l <- dim(mosdata)[1]
    mosdata$n <- rnbinom(l, mu = meanMos, size = dispMos)
  }
  or_vec <- (responses[1:n_arms] / (1-responses[1:n_arms])) / (responses[1] / (1-responses[1]))
  names(or_vec) <- aux
  
  mosdata <- 
    sim.glmm(design.data = mosdata,
             fixed.eff = list(
               intercept = qlogis(responses[1]),
               net = log(or_vec)),
             rand.V = c(observation = varO),
             distribution = "binomial")
  return(mosdata)
}

simulate_trial_IRS <- function(trial_days, n_arms, rep_IRS, rep_C,  mos_det = 0, meanMos,
                               dispMos = 1.5, responses, varO = 0.9){
  if(meanMos < 0 | dispMos < 0 ){ 
    print('Parameters for mosquito counts should be positive')
    return(-9)
  }
  
  nhuts <- (n_arms-1)*rep_IRS + rep_C
  aux <- repeat_fn(sttr = c('C',paste0('E',seq(1,n_arms-1))), repp = rep_IRS, repC = rep_C)
  aux_resp <- repeat_fn(sttr = responses, repp = rep_IRS, repC = rep_C)
  
  mosdata <-
    expand.grid(
      hut = seq(1,nhuts),
      day = 1:trial_days
    )
  mosdata$net <- NA #change from net to something else??
  for(i in 1:trial_days){
    mosdata[mosdata$day == i,]$net = aux
  }
  
  mosdata$sleeper <- NA
  mosdata$sleeper <- ((mosdata$hut + mosdata$day) %% nhuts) + 1
  mosdata$observation <- factor(formatC(1:nrow(mosdata), flag="0", width=3))
  mosdata$sleeper <- factor(mosdata$sleeper)
  mosdata$net <- factor(mosdata$net)
  mosdata$hut <- factor(mosdata$hut)
  
  mosdata$n <- NA
  if(mos_det==1){
    mosdata$n <- meanMos
  }else{
    l <- dim(mosdata)[1]
    mosdata$n <- rnbinom(l, mu = meanMos, size = dispMos)
  }
  or_vec <- (responses[1:n_arms] / (1-responses[1:n_arms])) / (responses[1] / (1-responses[1]))
  names(or_vec) <- c('C',paste0('E',seq(1,n_arms-1)))
  
  mosdata <- 
    sim.glmm(design.data = mosdata,
             fixed.eff = list(
               intercept = qlogis(responses[1]),
               net = log(or_vec)),
             rand.V = c(observation = varO),
             distribution = "binomial")
  return(mosdata)
}  

hypothesis_test <- function(trial,aoi,NIM=0.7, dataset){ 
  count <- 0
  
  l <- length(unique(dataset$net))
  aux <- c('C',paste0('E',seq(1,l-1)))
  
  #MORTALITY TRIALS
  
  if(trial == 1 | trial == 9){
    if(length(aoi) != 2){
      print('Check aoi')
      return(-9)
    }
    #relevel
    dataset$net <- relevel(dataset$net,aux[aoi[1]]) 
    #levels(dataset$net)
    
    fit_model <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper, 
        family = binomial, data = dataset)
    #summary(fit_model)
    
    labl <- paste0('net',aux[aoi[2]])
    
    if(coef(summary(fit_model))[labl, "Pr(>|z|)"] <0.05 & 
       coef(summary(fit_model))[labl,"Estimate"]>0){
      return(1)
    }else{
      return(0)
    }
    
    count <- count + 1
  }
  if(trial == 2){
    if(length(aoi) != 4){
      print('Check aoi')
      return(-8)
    }
    dataset2 <- dataset
    dataset2$net[dataset2$net ==  aux[aoi[2]] ] <- aux[aoi[1]]
    dataset2$net[dataset2$net ==  aux[aoi[4]] ] <-  aux[aoi[3]]
    table(dataset2$net)
    
    dataset2$net <- relevel(dataset2$net,aux[aoi[1]]) 
    levels(dataset2$net)
    
    fit_model <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper,
        family = binomial, data = dataset2)
    summary(fit_model)
    
    labl <- paste0('net',aux[aoi[3]])
    
    if(coef(summary(fit_model))[labl, "Pr(>|z|)"] <0.05 & 
       coef(summary(fit_model))[labl,"Estimate"]>0){
      return(1)
    }else{
      return(0)
    }
    
    count <- count + 1
  }
  if(trial == 3 | trial == 10){
    if(length(aoi) != 2){
      print('Check aoi')
      return(-7)
    }
    #relevel
    dataset$net <- relevel(dataset$net,aux[aoi[1]]) 
    #levels(dataset$net)
    
    fit_n <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper, 
        family = binomial, data = dataset) 
    summary(fit_n)
    
    labl <- paste0('net', aux[aoi[2]])
    
    exp(coef(summary(fit_n))[labl,'Estimate'])
    exp(coef(summary(fit_n))[labl,'Estimate'] - 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    exp(coef(summary(fit_n))[labl,'Estimate'] + 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    
    if(exp(coef(summary(fit_n))[labl,'Estimate'] -
           1.96*coef(summary(fit_n))[labl,'Std. Error']) > NIM){
      return(1)
    }else{
      return(0)
    }
    
    count <- count + 1
  }
  if(trial == 4){
    if(length(aoi) != 4){
      print('Check aoi')
      return(-6)
    }
    
    dataset2 <- dataset
    dataset2$net[dataset2$net ==  aux[aoi[2]] ] <- aux[aoi[1]]
    dataset2$net[dataset2$net ==  aux[aoi[4]] ] <-  aux[aoi[3]]
    table(dataset2$net)
    
    dataset2$net <- relevel(dataset2$net, aux[aoi[1]]) 
    levels(dataset2$net)
    
    fit_n <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper,
        family = binomial, data = dataset2)
    summary(fit_n)
    
    labl <- paste0('net',aux[aoi[3]])
    #exp(coef(summary(fit_n))[labl,'Estimate'])
    #exp(coef(summary(fit_n))[labl,'Estimate'] - 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    #exp(coef(summary(fit_n))[labl,'Estimate'] + 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    
    if(exp(coef(summary(fit_n))[labl,'Estimate'] -
           1.96*coef(summary(fit_n))[labl,'Std. Error']) > NIM){
      return(1)
    }else{
      return(0)
    }
    count <- count + 1
  }
  #BLOOD FEEDING INHIBITION
  if(trial == 5 | trial == 11){
    if(length(aoi) != 2){
      print('Check aoi')
      return(-9)
    }
    #relevel
    dataset$net <- relevel(dataset$net,aux[aoi[1]]) 
    #levels(dataset$net)
    
    fit_model <-
      glm(
        cbind(response, n - response) ~
          net  + day + hut + sleeper,
        family = binomial, data = dataset)
    summary(fit_model)
    
    labl <- paste0('net',aux[aoi[2]])
    
    if(coef(summary(fit_model))[labl, "Pr(>|z|)"] <0.05 & 
       coef(summary(fit_model))[labl,"Estimate"]<0){ #check
      return(1)
    }else{
      return(0)
    }
    count <- count + 1
  }
  if(trial == 6){
    if(length(aoi) != 4){
      print('Check aoi')
      return(-8)
    }
    dataset2 <- dataset
    dataset2$net[dataset2$net ==  aux[aoi[2]] ] <- aux[aoi[1]]
    dataset2$net[dataset2$net ==  aux[aoi[4]] ] <-  aux[aoi[3]]
    
    dataset2$net <- relevel(dataset2$net,aux[aoi[1]]) 
    #levels(dataset2$net)
    
    fit_model <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper,
        family = binomial, data = dataset2)
    #summary(fit_model)
    labl <- paste0('net',aux[aoi[3]])
    
    if(coef(summary(fit_model))[labl, "Pr(>|z|)"] <0.05 & 
       coef(summary(fit_model))[labl,"Estimate"]<0){ #check
      return(1)
    }else{
      return(0)
    }
    count <- count + 1
  }
  if(trial == 7 | trial == 12){
    if(length(aoi) != 2){
      print('Check aoi')
      return(-7)
    }
    if(NIM < 1){
      print('Check NIM: it should be >1 for blood-feeding inhibition')
      return(-8)
    }
    #relevel
    dataset$net <- relevel(dataset$net,aux[aoi[1]]) 
    #levels(dataset$net)
    
    fit_n <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper,
        family = binomial, data = dataset) 
    #summary(fit_n)
    labl <- paste0('net', aux[aoi[2]])
    
    #exp(coef(summary(fit_n))[labl,'Estimate'])
    #exp(coef(summary(fit_n))[labl,'Estimate'] - 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    #exp(coef(summary(fit_n))[labl,'Estimate'] + 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    
    if(exp(coef(summary(fit_n))[labl,'Estimate'] +
           1.96*coef(summary(fit_n))[labl,'Std. Error']) < NIM){
      return(1)
    }else{
      return(0)
    }
    count <- count + 1
  }
  if(trial == 8){
    if(length(aoi) != 4){
      print('Check aoi')
      return(-6)
    }
    if(NIM < 1){
      print('Check NIM: it should be >1 for blood-feeding inhibition')
      return(-8)
    }
    
    dataset2 <- dataset
    dataset2$net[dataset2$net ==  aux[aoi[2]] ] <- aux[aoi[1]]
    dataset2$net[dataset2$net ==  aux[aoi[4]] ] <-  aux[aoi[3]]
    
    dataset2$net <- relevel(dataset2$net, aux[aoi[1]]) 
    #levels(dataset2$net)
    
    fit_n <-
      glm(
        cbind(response, n - response) ~
          net + day + hut + sleeper,
        family = binomial, data = dataset2)
    #summary(fit_n)
    
    labl <- paste0('net',aux[aoi[3]])
    #exp(coef(summary(fit_n))[labl,'Estimate'])
    #exp(coef(summary(fit_n))[labl,'Estimate'] - 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    #exp(coef(summary(fit_n))[labl,'Estimate'] + 1.96*coef(summary(fit_n))[labl,'Std. Error'])
    
    if(exp(coef(summary(fit_n))[labl,'Estimate'] +
           1.96*coef(summary(fit_n))[labl,'Std. Error']) < NIM){
      return(1)
    }else{
      return(0)
    }
    count <- count + 1
  }
  
  if(count<1){
    return('Check value of trial (analyses not carried out)')
  }
}

#Add NIM?
power_calculator_ITN <- function(parallelise = 0, trial, npw, rotations = 1, varO = 0.9, 
                 nsim = 1000, n_arms, mos_det = 0, meanMos, dispMos = 1.5,
                 aoi, responses, NIM = 0.7){
  if(parallelise!=0 & parallelise!=1 & parallelise!=2){
    print('Operation was not executed. Parallelise must take a value of (i) 0 (code not parallelised); (ii) Parallelised for Windows; (iii) Parallilised for Mac. If in doubt, set to zero')
    return(-9)
  }
  if(parallelise==0){
    sz <- lapply(1:nsim, function(...) hypothesis_test(trial = trial, aoi = aoi, 
                                                    NIM = NIM, dataset = simulate_trial_ITN(n_arms = n_arms, npw = npw, 
                                                                                   rotations = rotations, mos_det = mos_det, meanMos = meanMos, 
                                                                                   dispMos = dispMos, responses = responses,
                                                                                   varO = varO)))
  }
  if(parallelise==1){ # Under construction
    ncores <- detectCores() - 1
    
    n_armsY <- n_arms
    trialY <- trial
    aoiY <- aoi
    NIMY <- NIM
    dispMosY <- dispMos
    meanMosY <- meanMos
    varOY <- varO
    responsesY <- responses
    npwY <- npw
    mos_detY <- mos_det
    rotationsY <- rotations
    
    cl <- makeCluster(ncores)
    clusterExport(cl,'npwY', envir = environment())
    clusterExport(cl,'rotationsY', envir = environment())
    clusterExport(cl,'responsesY', envir = environment())
    clusterExport(cl,'varOY', envir = environment())
    clusterExport(cl,'mos_detY', envir = environment())
    clusterExport(cl,'meanMosY', envir = environment())
    clusterExport(cl,'dispMosY', envir = environment())
    clusterExport(cl,'n_armsY', envir = environment())
    clusterExport(cl,'aoiY', envir = environment())
    clusterExport(cl,'NIMY', envir = environment())
    clusterExport(cl,'trialY', envir = environment())
    clusterExport(cl,'hypothesis_test')
    clusterExport(cl,'simulate_trial')
    #clusterExport(cl,'lazy') #replace with name of this function???
    clusterExport(cl,'InvLogit')
    clusterEvalQ(cl, {
      library(lme4)
      library(GLMMmisc)
      #library("optimx")
    })
    sz <- parLapply(cl, 1:nsim, function(...) hypothesis_test(trial = trialY, aoi = aoiY, NIM = NIMY,
                                                              dataset = simulate_trial_ITN(n_arms = n_armsY, npw = npwY,
                                                                                       rotations = rotationsY, mos_det = mos_detY, meanMos = meanMosY,
                                                                                       dispMos = dispMosY, responses = responsesY,
                                                                                       varO = varOY)))
    stopCluster(cl)
  }
  if(parallelise==2){
    ncores <- detectCores() - 1
    sz <- mclapply(1:nsim, function(...) hypothesis_test(trial = trial, aoi = aoi, NIM = NIM,
                                                dataset = simulate_trial_ITN(n_arms = n_arms, npw = npw, 
                                                              rotations = rotations, mos_det = mos_det, meanMos = meanMos, 
                                                                dispMos = dispMos, responses = responses,
                                                                  varO = varO)))
  }
  mm <- unlist(sz)
  print(paste0('Power Estimate: ',100*sum(mm)/length(mm),'%, 95% CI: [',
               100*round(binom.test(table(factor(mm,c(1,0))))$conf.int[1],3),',',
               100*round(binom.test(table(factor(mm,c(1,0))))$conf.int[2],3),']'))
  return(c(100*sum(mm)/length(mm),100*binom.test(table(factor(mm,c(1,0))))$conf.int[1],
           100*binom.test(table(factor(mm,c(1,0))))$conf.int[2]))
}

power_calculator_IRS <- function(parallelise = 0, trial, varO = .9, nsim = 400, 
                             trial_days, n_arms, rep_IRS, rep_C, mos_det = 0, meanMos,
                             dispMos = 1.5, aoi, responses){
  if(parallelise!=0 & parallelise!=1 & parallelise!=2){
    print('Operation was not executed. Parallelise must take a value of (i) 0 (code not parallelised); (ii) Parallelised for Windows; (iii) Parallilised for Mac. If in doubt, set to zero')
    return(-9)
  }
  if(parallelise==0){
    sz <- lapply(1:nsim, function(...) hypothesis_test(trial = trial, aoi = aoi, 
                                          dataset = simulate_trial_IRS(trial_days = trial_days, n_arms = n_arms,
                                             rep_IRS = rep_IRS, rep_C = rep_C, mos_det = mos_det, meanMos = meanMos, 
                                              dispMos = dispMos, responses = responses,varO = varO)))
  }
  if(parallelise==1){ # Under construction
    ncores <- detectCores() - 1
    
    n_armsY <- n_arms
    trial_daysY <- trial_days
    rep_CY <- rep_C
    rep_IRSY <- rep_IRS
    trialY <- trial
    aoiY <- aoi
    dispMosY <- dispMos
    meanMosY <- meanMos
    varOY <- varO
    responsesY <- responses
    mos_detY <- mos_det
    
    cl <- makeCluster(ncores)
    clusterExport(cl,'rep_CY', envir = environment())
    clusterExport(cl,'rep_IRSY', envir = environment())
    clusterExport(cl,'trial_daysY', envir = environment())
    clusterExport(cl,'responsesY', envir = environment())
    clusterExport(cl,'varOY', envir = environment())
    clusterExport(cl,'mos_detY', envir = environment())
    clusterExport(cl,'meanMosY', envir = environment())
    clusterExport(cl,'dispMosY', envir = environment())
    clusterExport(cl,'n_armsY', envir = environment())
    clusterExport(cl,'aoiY', envir = environment())
    clusterExport(cl,'trialY', envir = environment())
    clusterExport(cl,'hypothesis_test')
    clusterExport(cl,'simulate_trial_IRS')
    clusterExport(cl,'repeat_fn')
    clusterExport(cl,'InvLogit')
    clusterEvalQ(cl, {
      library(lme4)
      library(GLMMmisc)
      #library("optimx")
    })
    sz <- parLapply(cl, 1:nsim, function(...) hypothesis_test(trial = trialY, aoi = aoiY,
                                                dataset = simulate_trial_IRS(n_arms = n_armsY, 
                                                       rep_C = rep_CY, rep_IRS = rep_IRSY, trial_days = trial_daysY,
                                                      mos_det = mos_detY, meanMos = meanMosY,
                                                    dispMos = dispMosY, responses = responsesY,
                                                     varO = varOY)))
    stopCluster(cl)
  }
  if(parallelise==2){
    ncores <- detectCores() - 1
    sz <- mclapply(1:nsim, function(...) hypothesis_test(trial = trial, aoi = aoi, 
                                            dataset = simulate_trial_IRS(n_arms = n_arms, 
                                                  rep_C = rep_C, rep_IRS = rep_IRS,
                                              trial_days = trial_days, mos_det = mos_det, meanMos = meanMos, 
                                                dispMos = dispMos, responses = responses,
                                                  varO = varO))) #add mc.cores = ncores?
  }
  mm <- unlist(sz)
  print(paste0('Power Estimate: ',100*round(sum(mm)/length(mm),3),'%, 95% CI: [',
               100*round(binom.test(table(factor(mm,c(1,0))))$conf.int[1],3),',',
               100*round(binom.test(table(factor(mm,c(1,0))))$conf.int[2],3),']'))
  return(c(100*sum(mm)/length(mm),100*binom.test(table(factor(mm,c(1,0))))$conf.int[1],
           100*binom.test(table(factor(mm,c(1,0))))$conf.int[2]))
}
