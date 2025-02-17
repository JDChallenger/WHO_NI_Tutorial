# Convert from logit scale to probability (or proportion) scale
InvLogit <- function(X){
  exp(X)/(1+exp(X))
}


#TIDY, and change name!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    for(i in 1:(length(lu)/2)){ #number of LLINs
      #cone_data[cone_data$llin_code==paste('A',i, sep = '_'),]$day <- i 
      #cone_data[cone_data$llin_code==paste('B',i, sep = '_'),]$day <- i 
      cone_data[cone_data$llin_code==lu[2*i - 1],]$day <- i 
      cone_data[cone_data$llin_code==lu[2*i],]$day <- i 
      pday <- 2*dim(cone_data[cone_data$llin_code==lu[2*i - 1],])[1] # '$day taken off'
      #print(dim(cone_data[cone_data$llin_code==lu[2*i - 1],]))
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
  
  return(cone_data[, !(names(cone_data) %in% c('LO'))]) 
  
}

#Function for non-inf assessment for cone assay data
cone_NIM <- function(dataset, NIM_pc = 0.07, int_cat = 'A',  FE_label = 'armB', verbose = T){
  
  FIC_mortality <- sum(dataset[dataset$arm==int_cat,]$tot_dead)/sum(dataset[dataset$arm==int_cat,]$total)
  if(verbose == T){
    print(FIC_mortality)
  }
  non_inf_margin <- ((FIC_mortality - NIM_pc) / (1- (FIC_mortality - NIM_pc))) / (FIC_mortality / (1- FIC_mortality)) 
  if(verbose == T){
    print(non_inf_margin)
  }
  
  fit <- glm(
    cbind(tot_dead, total - tot_dead) ~ arm + day, 
    family = binomial, data = dataset 
  )
  summary(fit)
  
  OR1 <- exp(coef(summary(fit))[FE_label,"Estimate"])
  OR1_lower <- exp(coef(summary(fit))[FE_label,"Estimate"] - 
                     1.96*coef(summary(fit))[FE_label,'Std. Error'])
  OR1_upper <- exp(coef(summary(fit))[FE_label,"Estimate"] + 
                     1.96*coef(summary(fit))[FE_label,'Std. Error'])
  
  aux2 <- -1
  if(OR1_lower > non_inf_margin){
    aux2 <- 1
  }else{
    aux2 <- 0
  }
  #if verbose is T, print OR & 95% CI?
  
  if(aux2 < 0){
    print('Error detected')
  }
  return(aux2)
}

#Clean & annotate!!
tunnel_sim2 <- function(mort = 0.4, #reps = 1, #need reps???
                        npos = 3, n_nets = 30,
                        sigma_net = 0.7, n_tunnels = 5, verbose = F, n_mosq = 50){
  
  #what if this isn't an integer??
  n_day <- ceiling(2*n_nets*npos / n_tunnels) # number of arms * number of nets * number of net positions
  if(verbose==T){
    print(n_day)
  }
  
  tunnel_data <- data.frame(expand.grid(
    day = seq(1,n_day),
    tunnel = seq(1,n_tunnels)#,
  ))
  
  net_list <- c(paste('A',seq(1,n_nets), sep = '_'),
                paste('B',seq(1,n_nets), sep = '_'))
  #shuffle order
  net_list <- sample(net_list)
  tunnel_data$llin_code <- rep(net_list,3)#c(net_list, net_list, net_list) # repeat three times?
  
  pos_list <- c()
  for(i in 1:npos){
    pos_list <- c(pos_list,rep(i,2*n_nets))  # for two arms
  }
  
  tunnel_data$net_position <- pos_list#c(rep(1,2*n_nets), rep(2,2*n_nets), rep(3,2*n_nets))
  
  tunnel_data$arm <- substr(tunnel_data$llin_code,1,1)
  tunnel_data$day <- as.factor(tunnel_data$day)
  tunnel_data$tunnel <- as.factor(tunnel_data$tunnel)
  tunnel_data$total <- n_mosq
  
  lst <- unique(tunnel_data$llin_code)
  lstA <- unique(tunnel_data$arm)
  lstD <- unique(tunnel_data$day)
  lstT <- unique(tunnel_data$tunnel)
  
  vec_arm <- qlogis(c(mort,mort))
  vec_llin <- rnorm(2*n_nets, 0, sigma_net)
  vec_day <- rnorm(length(unique(tunnel_data$day)),0,0.25)
  vec_tunnel <- rnorm(length(unique(tunnel_data$tunnel)),0,0.25)
  
  #print(head(tunnel_data))
  
  ll <- dim(tunnel_data)[1]
  #what affects mortality here? Arm & llin_code
  tunnel_data$LO <- NA
  for(i in 1:ll){
    s1 <- which(lstA==tunnel_data[i,5]) # arm
    s2 <- which(lst==tunnel_data[i,3]) #llin code
    s3 <- which(lstD==tunnel_data[i,1]) #day
    s4 <- which(lstT==tunnel_data[i,2]) #tunnel
    
    #print(c(s1,s2,s3,s4))
    
    bray <- vec_arm[s1] + vec_llin[s2] + vec_day[s3] + vec_tunnel[s4] #+ vec_net[s5]
    tunnel_data$LO[i] <- bray
    
  }
  
  #simulate trial
  tunnel_data$tot_dead <- NA
  for(i in 1:ll){
    aux <- rbinom(1, tunnel_data$total[i], InvLogit(tunnel_data$LO[i]))
    tunnel_data$tot_dead[i] <- aux
  }
  return(tunnel_data[, !(names(tunnel_data) %in% c('LO'))])
}

tunnel_NIM <- function(dataset, NIM_pc = 0.07, int_cat = 'A',  FE_label = 'armB', verbose = T){
  
  FIC_mortality <- sum(dataset[dataset$arm==int_cat,]$tot_dead)/sum(dataset[dataset$arm==int_cat,]$total)
  if(verbose == T){
    print(FIC_mortality)
  }
  non_inf_margin <- ((FIC_mortality - NIM_pc) / (1- (FIC_mortality - NIM_pc))) / (FIC_mortality / (1- FIC_mortality)) 
  if(verbose == T){
    print(non_inf_margin)
  }
  
  fit <- glm(
    cbind(tot_dead, total - tot_dead) ~ arm + day + tunnel, #  + compartment. As arms don't move
    family = binomial, data = dataset # looks weird?
  )
  summary(fit)
  
  OR1 <- exp(coef(summary(fit))[FE_label,"Estimate"])
  OR1_lower <- exp(coef(summary(fit))[FE_label,"Estimate"] - 
                     1.96*coef(summary(fit))[FE_label,'Std. Error'])
  OR1_upper <- exp(coef(summary(fit))[FE_label,"Estimate"] + 
                     1.96*coef(summary(fit))[FE_label,'Std. Error'])
  
  aux2 <- -1
  if(OR1_lower > non_inf_margin){
    aux2 <- 1
  }else{
    aux2 <- 0
  }
  
  if(aux2 < 0){
    print('Error detected')
  }
  return(aux2)
}

#change name, add more annotation. Should the number of days be left empty? For the user to define
IACT_sim2 <- function(n_day = 24, n_comp = 18, n_mosq = 15, n_arm = 9, verbose = F,
                      mortalities = c(0.05,0.6,0.25,0.5,0.35,0.25,0.60,0.25,0.25),
                      sigma_net = 0.4, n_nets = 30){
  if(n_comp%%n_arm !=0){
    print('Results may be unrealiable if the number of compartments is 
         not a multiple of (or equal to) the number of trial arms')
  } 
  if(length(mortalities) != n_arm  ){
    print('The number of arms (n_arm) should match the length of the mortalities vector')
  }
  if(n_nets < 30){
    print('Please select at least 30 nets per arm')
    #break?
  }
  
  mosdata <-
    expand.grid(
      compartment = factor(1:n_comp),
      day = 1:n_day
    )
  
  aux3 <- 1:n_comp #sample(1:18) # number of sleepers must match number of compartments
  
  mosdata$sleeper <- NA
  for(i in 1:n_day){
    bm <- i%%n_comp
    #print(bm)
    if(bm==0){
      #print('bray')
      mosdata[mosdata$day==i,]$sleeper <- c(aux3[c(n_comp:n_comp)],aux3[seq_len(n_comp-1)])
      
    }else{
      #print(c(aux3[c(bm:n_comp)],aux3[seq_len(n_comp-1)]))
      mosdata[mosdata$day==i,]$sleeper <- c(aux3[c(bm:n_comp)],aux3[seq_len(bm-1)])
    }
    #print(i)
  }
  
  aux2 <- paste('N', 1:n_arm, sep = '')
  ratio <- round(n_comp/n_arm)
  aux2b <- rep(aux2,ratio) # This assumes that we have twice as many compartments as arms
  aux2d <- sample(aux2b)#aux2b[order(aux2b)]
  
  mosdata$net <- NA
  for(i in 1:n_day){
    mosdata[mosdata$day==i,]$net <- aux2d
  }
  
  mosdata$netcode <- NA
  bray <- n_day * n_comp / n_arm / n_nets # not necessarily an integer
  #which means
  bm <- (n_day * n_comp / n_arm ) %% n_nets
  lu <- unique(mosdata$net)
  if(bm==0){
    for(i in 1:length(lu)){
      if(verbose==T){
        print(paste0('Each net is used ',bray,' times'))
      }
      
      mosdata[mosdata$net==lu[i],]$netcode <- rep(paste0(lu[i],'_',seq(1,n_nets)),bray)  # if bray an integer!!
    }
  }else{
    flr <- floor(n_day * n_comp / n_arm / n_nets)
    rem <- (n_day * n_comp / n_arm ) %% n_nets
    if(verbose==T){
      print(paste0('In each arm, ',n_nets - rem,' nets are used ',flr,' times; ',rem,' are used ',flr + 1,' times.'))
    }  
    
    lu <- unique(mosdata$net)
    for(i in 1:length(lu)){
      mosdata[mosdata$net==lu[i],]$netcode <- c(rep(paste0(lu[i],'_',seq(1,n_nets)),flr), 
                                                paste0(lu[i],'_',seq(1,rem)))  
    }
  }
  
  #add mosquitoes
  mosdata$total <- n_mosq
  
  vec_arm <- qlogis(mortalities)
  
  #Variability due to... sleepers? Compartments? Night? Individual net??
  #For each data point, work out a bespoke log-odds mortality
  
  vec_day <- rnorm(n_day,0,0.2)
  vec_comp <- rnorm(n_comp,0,0.2)
  vec_sleep <- rnorm(n_comp,0,0.1)
  vec_net <- rnorm(n_nets*n_arm,0,sigma_net) 
  
  ll <- dim(mosdata)[1]
  aux4 <- unique(mosdata$netcode)
  
  mosdata$LO <- NA
  for(i in 1:ll){
    s1 <- which((1:n_comp) == mosdata[i,1]) # comp
    s2 <- which((1:n_day) == mosdata[i,2]) #day
    s3 <- which((1:n_comp) == mosdata[i,3]) #sleeper
    s4 <- which(aux2 == mosdata[i,4])
    s5 <- which(aux4 == mosdata[i,5])
    
    bray <- vec_comp[s1] + vec_day[s2] + vec_sleep[s3] + vec_arm[s4] + vec_net[s5]
    mosdata$LO[i] <- bray
  }
  
  #simulate trial
  mosdata$tot_dead <- NA
  for(i in 1:ll){
    aux <- rbinom(1, mosdata$total[i], InvLogit(mosdata$LO[i]))
    mosdata$tot_dead[i] <- aux
  }
  if(verbose==T){
    print(paste0('Dim. of data: ',dim(mosdata)))
  }
  
  return(mosdata[, !(names(mosdata) %in% c('LO'))]) 
}  

#######################################################################
# function to do non-inf assessment
IACT_NIM <- function(dataset, NIM_pc = 0.07, int_cat = 'N3',  FE_label = 'netN6', verbose = T){
  
  FIC_mortality <- sum(dataset[dataset$net==int_cat,]$tot_dead)/sum(dataset[dataset$net==int_cat,]$total)
  if(verbose == T){
    print(FIC_mortality)
  }
  non_inf_margin <- ((FIC_mortality - NIM_pc) / (1- (FIC_mortality - NIM_pc))) / (FIC_mortality / (1- FIC_mortality)) 
  if(verbose == T){
    print(non_inf_margin)
  }
  
  dataset$sleeper <- as.factor(dataset$sleeper)
  dataset$day <- as.factor(dataset$day)
  dataset$net <- as.factor(dataset$net)
  dataset$net <- relevel(dataset$net, 'N3') 
  levels(dataset$net)
  
  fit <- glm(
    cbind(tot_dead, total - tot_dead) ~ net + day + sleeper, 
    family = binomial, data = dataset 
  )
  summary(fit)
  
  OR1 <- exp(coef(summary(fit))[FE_label,"Estimate"])
  OR1_lower <- exp(coef(summary(fit))[FE_label,"Estimate"] - 
                     1.96*coef(summary(fit))[FE_label,'Std. Error'])
  OR1_upper <- exp(coef(summary(fit))[FE_label,"Estimate"] + 
                     1.96*coef(summary(fit))[FE_label,'Std. Error'])
  
  aux2 <- -1
  if(OR1_lower > non_inf_margin){
    aux2 <- 1
  }else{
    aux2 <- 0
  }
  #if verbose is T, print OR & 95% CI?
  if(verbose == T){
    print(paste0(round(OR1,4),' [',round(OR1_lower,4),', ',round(OR1_upper,4),']'))
  }
  
  if(aux2 < 0){
    print('Error detected')
  }
  return(aux2)
}
