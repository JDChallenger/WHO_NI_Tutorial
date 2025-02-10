
InvLogit <- function(X){
  exp(X)/(1+exp(X))
}


#TIDY, and change name!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Don't output the LO????
EHT_sim2 <- function(n_arms, npr = 9, mos_det = 0, meanMos, dispMos = 1.5, verbose = F,
                     rotations = 1,sigma_hut = 0.2, sigma_sleep = 0.3, sigma_net = 0.4, sigma_day = 0.1,
                     mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25), ITN_reps = 30){
  #Check length(responses) == n_arms
  if(n_arms != length(mortalities)){
    print('Operation was not executed. Check the number of arms, and corresponding vector of mosquito mortalities')
    return(-9)
  }
  if(npr > n_arms){
    print('Operation was not executed. Rethink trial design (npr should be less than or equal to n_arms)')
    return(-9)
  }
  if(rotations < 1){ #At the moment, we need at least 1 rotation
    rotations <- 1
    print('Trial simulated for rotations = 1')
  }
  #if(varO < 0 ){ #Reject negative variance for the random effect
  # print('The variance of a random effect must be positive')
  #  return(-9)
  #}
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
      night = 1:npr
    )
  mosdata <- mosdata[order(mosdata$hut, mosdata$week, mosdata$night),] # would day be better than night? 
  
  mosdata$net <- NA
  count <- 1
  for(i in 1:n_arms){
    for(j in 1:n_arms){
      mosdata$net[((npr*(count-1))+1):(npr*count)]  <- rep(aux2[i,j],npr)
      count <- count + 1
    }
  }
  #table(mosdata$net, useNA = 'a')
  mosdata$sleeper <- NA
  aux3 <- sample(1:n_arms) #was '7'- a mistake?
  
  for(i in 1:n_arms){ # number of weeks for one rotation
    aux4 <- c(aux3[c(i:n_arms)],aux3[seq_len(i-1)])
    #print(aux4)
    for(j in 1:npr){ #nights per round
      mosdata[mosdata$week==i & mosdata$night==j,]$sleeper <- 
        c(aux4[c(j:n_arms)],aux4[seq_len(j-1)])
    }
  }
  
  #Do we have rotation > 1? What if rotation > 2??
  if(rotations > 1){
    if(rotations <= 2){
      mosdata$day <- npr*(mosdata$week - 1) + mosdata$night
      extra <- round((rotations - 1)*dim(mosdata)[1] / n_arms)
      mosdataX <- mosdata[mosdata$day <= extra,]
      mosdataX$day <- mosdataX$day + max(mosdata$day)
      mosdataX$week <- mosdataX$week + max(mosdata$week)
      #Now add it on to the existing trial design
      mosdata <- rbind(mosdata,mosdataX)
    }else{
      mosdata$day <- npr*(mosdata$week - 1) + mosdata$night
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
    mosdata$day <- npr*(mosdata$week - 1) + mosdata$night
  }
  
  #Now make a unique identifier for each data point
  mosdata$observation <- factor(formatC(1:nrow(mosdata), flag="0", width=3))
  mosdata$sleeper <- factor(mosdata$sleeper)
  mosdata$net <- factor(mosdata$net)
  mosdata$day <- factor(mosdata$day)
  
  mosdata$n <- NA
  if(mos_det==1){
    mosdata$n <- meanMos
  }else{
    l <- dim(mosdata)[1]
    mosdata$n <- rnbinom(l, mu = meanMos, size = dispMos)
  }
  
  mosdata$net_id <- NA # round down, so that the last net has extra data pts
  data_pts_per_ITN <- dim(mosdata[mosdata$net =='C',])[1]/ITN_reps # NOTE: this may now vary by arm, if we're doubling up. Need to do iteratively, for each arm??
  if(verbose==T){
    print(paste0('data pts per replicate: ',data_pts_per_ITN))
  }
  remainder <- dim(mosdata[mosdata$net =='C',])[1] %% ITN_reps # this could be zero
  if(verbose==T){
    print(paste0(remainder, ' replicates get ',ceiling(data_pts_per_ITN),' data points; ',ITN_reps - remainder,' replicates get ',floor(data_pts_per_ITN),' data points'))
  }
  
  #for each arm in turn, generate the entries for net_id. I'm assuming that order doesn't matter here
  lu <- unique(as.character(mosdata$net))
  for(i in 1:length(lu)){
    lizt <- c()
    for(j in 1:remainder){
      for(k in 1:ceiling(data_pts_per_ITN)){
        lizt <- c(lizt, paste(lu[i], j, sep = '_'))
      }
    }
    for(j in 1:(ITN_reps - remainder)){
      for(k in 1:floor(data_pts_per_ITN)){
        lizt <- c(lizt, paste(lu[i], j + remainder, sep = '_'))
      }
    }
    #print(paste0('arm: ',lu[i],', \n',lizt))
    mosdata[mosdata$net==lu[i],]$net_id <- lizt
  }
  
  
  # n_nets <- round(dim(mosdata[mosdata$net =='C',])[1]/4)#ceiling(dim(mosdata[mosdata$net =='C',])[1]/4)#ceiling(table(mosdata$net)[1][[1]]/4)
  # if(verbose == T){
  #   #print(table(mosdata$net))
  #   print('n_nets: ',n_nets)
  # }
  # 
  # 
  # lu <- unique(as.character(mosdata$net))
  # mosdata2 <- data.frame()
  # for(i in 1:length(lu)){ # length(lu)
  #   aux <- mosdata[mosdata$net==lu[i],]
  #   ln <- dim(aux)[1]
  #   #count <- 0
  #   for(j in 1:n_nets){
  #     if(j < n_nets){
  #       aux$net_id[(4*(j-1) + 1) : (4*j)] <- paste(lu[i], j, sep = '_')
  #     }else{
  #       #aux$net_id[(4*(j-1) + 1) : (4*j - 2)] <- paste(lu[i], j, sep = '_')
  #       aux$net_id[(4*(j-1) + 1) : ln] <- paste(lu[i], j, sep = '_')
  #       #print(count)
  #     }
  #     #count <- count + 1
  #   }
  #   mosdata2 <- rbind(mosdata2, aux)
  # }
  # #View(mosdata2)
  # #table(mosdata2$net_id)
  # 
  
  
  #MOrtalities:
  # C: Untreated Control (5% mortality)
  # E1: Unwashed Test ITN (60% mortality)
  # E2: 20x Washed Test ITN (25% mortality)
  # E3: 12mth field-aged ITN (50% mortality)
  # E4: 24mth field-aged ITN (35% mortality)
  # E5: 36mth field-aged ITN (25% mortality)
  # E6: Unwashed positive control (60% mortality)
  # E7: 20x washed positive control (25% mortality)
  # E8: 36mth field-aged positive control (25% mortality)
  
  #Hence:
  vec_arm <- qlogis(mortalities)#qlogis(c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25))
  
  vec_day <- rnorm(length(unique(mosdata$day)),0,sigma_day)#rnorm(n_arms*npr,0,sigma_day)
  vec_hut <- rnorm(n_arms,0,sigma_hut)
  vec_sleep <- rnorm(n_arms,0,sigma_sleep)
  vec_net <- rnorm(ITN_reps*n_arms,0,sigma_net) # n_nets is now a user-supplied argument
  
  ll <- dim(mosdata)[1]
  lu2 <- unique(mosdata$net_id)
  lud <- unique(mosdata$day)
  
  if(verbose==T & length(vec_net)==length(lu2)){
    print('OK')
  }
  
  mosdata$LO <- NA
  for(i in 1:ll){
    s1 <- which((1:n_arms)==mosdata[i,1]) # hut
    s2 <- which(lud==mosdata[i,6]) #day
    s3 <- which((1:n_arms)==mosdata[i,5]) #sleeper
    s4 <- which(lu==mosdata[i,4]) # arm
    s5 <- which(lu2==mosdata[i,9]) # net ID
    
    #print(c(s1,s2,s3,s4,s5))
    
    bray <- vec_hut[s1] + vec_day[s2] + vec_sleep[s3] + vec_arm[s4] + vec_net[s5]
    #print(bray)
    mosdata$LO[i] <- bray
    
  }
  
  #simulate trial
  mosdata$tot_dead <- NA
  for(i in 1:ll){
    aux <- rbinom(1, mosdata$n[i], InvLogit(mosdata$LO[i]))
    mosdata$tot_dead[i] <- aux
  }
  
  #Then analyse data in GLMM
  
  mosdata$sleeper <- as.factor(mosdata$sleeper)
  mosdata$hut <- as.factor(mosdata$hut)
  mosdata$day <- as.factor(mosdata$day)
  
  return(mosdata)
}

#TIDY AND ANNOTATE!!!!!!!!!!! CHANGE NAME????????????
cone_sim2 <- function(cone_mort = 0.5, reps = 4, npos = 4, n_nets = 10, 
                      sigma_net = 0.7, nday = -9, verbose = T, num_mosq = 5){ # add day?
  cone_data <-
    data.frame(expand.grid(
      #period = 3,#c(0,1,2,3),
      net_position = seq(1,npos),
      replicates = seq(1,reps),
      llin_code = sample(c(paste('A',seq(1,n_nets), sep = '_'), paste('B',seq(1,n_nets), sep = '_')))
      #arm = c('A','B')
    ))
  cone_data$arm <- substr(cone_data$llin_code,1,1)
  
  cone_data$total <- num_mosq
  
  cone_data$day <- NA
  if(nday < 1){
    lu <- unique(cone_data$llin_code)
    for(i in 1:(lu/2)){ #number of LLINs
      #cone_data[cone_data$llin_code==paste('A',i, sep = '_'),]$day <- i 
      #cone_data[cone_data$llin_code==paste('B',i, sep = '_'),]$day <- i 
      cone_data[cone_data$llin_code==lu[2*i - 1],]$day <- i 
      cone_data[cone_data$llin_code==lu[2*i],]$day <- i 
      pday <- 2*dim(cone_data[cone_data$llin_code==lu[2*i - 1],]$day)[1]
    }
    if(verbose==T){
      print(paste0('Cone tests per day: ',pday))
    }
    #table(cone_data$day, useNA = 'a')
  }else{
    d1 <- dim(cone_data)[1]
    if(d1%%nday==0){ # is the number of data points an exact multiple of the number of days, or is there a remainder?
      data_points_per_day <- d1/nday
      if(verbose==T){
        print(paste0('Cone tests per day: ',data_points_per_day))
      }
      for(i in 1:nday){
        cone_data[(1 + (i-1)*data_points_per_day):(i*data_points_per_day),]$day <- i
      }
    }else{
      data_points_per_day <- floor(d1/nday)
      if(verbose==T){
        print(paste0('Cone tests per day: ',data_points_per_day))
        print('Note: the final day will have extra cone tests, to finish off the study')
      }
      for(i in 1:(nday)){
        cone_data[(1 + (i-1)*data_points_per_day):(i*data_points_per_day),]$day <- i
      }
      # finishing the study on the final day:
      # You could change the RHS to 'nday + 1', if you prefer to add an extra day to the study to finish off, 
      # instead of increasing the number of cone tests on the final day
      cone_data[(1 + nday*data_points_per_day):d1,]$day <- nday
    }
  }
  cone_data$day <- as.factor(cone_data$day)
  
  lst <- unique(cone_data$llin_code)
  lstA <- unique(cone_data$arm)
  lstD <- unique(cone_data$day)
  
  vec_arm <- qlogis(c(cone_mort,cone_mort))
  vec_llin <- rnorm(length(unique(cone_data$llin_code)),0,sigma_net)
  vec_day <- rnorm(length(unique(cone_data$day)),0,0.2)
  
  #ready to simulate?
  ll <- dim(cone_data)[1]
  #what affects mortality here? Arm & llin_code
  cone_data$LO <- NA
  for(i in 1:ll){
    s1 <- which(lstA==cone_data[i,4]) # arm
    s2 <- which(lst==cone_data[i,3]) #llin code
    s3 <- which(lstD==cone_data[i,6]) #day
    #print(c(s1,s2))
    
    bray <- vec_arm[s1] + vec_llin[s2] + vec_day[s3] #+ vec_arm[s4] + vec_net[s5]
    cone_data$LO[i] <- bray
    
  }
  #simulate trial
  cone_data$tot_dead <- NA
  for(i in 1:ll){
    aux <- rbinom(1, cone_data$total[i], InvLogit(cone_data$LO[i]))
    cone_data$tot_dead[i] <- aux
  }
  
  return(cone_data) 
  
}  