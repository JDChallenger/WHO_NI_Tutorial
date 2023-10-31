InvLogit <- function(X){
  exp(X)/(1+exp(X))
}

NIM_answer <- function(type = 'mort', vecOR = c(1,0.1,1.9), nim = 0.7){
  ans <- 'NOT non-inferior'
  count <- 0
  if(type =='mort'){
    count <- 1
    if(vecOR[2] > nim){
      ans <- 'Non-inferior'
    }
  }
  
  if(type == 'bfi'){
    if(nim < 1){
      'Check non-inferiority margin (should be greater than 1 for blood feeding)'
    }
    count <- 1
    if(vecOR[3] < nim){
      ans <- 'Non-inferior'
    }
  }
  if(count==0){
    'Not carried out'
  }
  return(ans)
}

summm <- function(data, vec, td = 'tot_dead', tot = 'total', table = 0, precision = 3){
  
  ud <- unique(vec)
  lud <- length(ud)
  if(table==0){
    for(i in 1:lud){
      print(paste0(ud[i],': ',round(100*sum(data[vec==ud[i],td])/sum(data[vec==ud[i],tot]),precision),'%'))
    }
  }else{
    dh <- data.frame('Arm' = as.character(), 'Percentage' = as.numeric())
    for(i in 1:lud){
     aux <- data.frame('Arm' = ud[i], 'Percentage' = round(100*sum(data[vec==ud[i],td])/sum(data[vec==ud[i],tot]),precision)) 
    dh <- rbind(dh,aux)
    }
    return(dh)
  }

  
}


##mortality summary for FE-only models
mFE <- function(model, vec, intercept, bfi = 0, name = "treatment", offset = 0, precision = 3){
  
  if(bfi==1){
    dfe <- data.frame('Arm' = as.character(),
                      'Blood Feeding' = as.numeric(),
                      'Lower_95pc_CI' = as.numeric(), 'Upper_95pc_CI' = as.numeric())
  }else{
    dfe <- data.frame('Arm' = as.character(),
                      'Mortality' = as.numeric(),
                      'Lower_95pc_CI' = as.numeric(), 'Upper_95pc_CI' = as.numeric())
  }
  dfe[1,] <- c(intercept,round(InvLogit(coef(summary(model))['(Intercept)','Estimate'] + offset),precision),
               round(InvLogit(coef(summary(model))['(Intercept)','Estimate'] - 
                                1.96*coef(summary(model))['(Intercept)','Std. Error'] + offset),precision),
               round(InvLogit(coef(summary(model))['(Intercept)','Estimate'] + 
                                1.96*coef(summary(model))['(Intercept)','Std. Error'] + offset),precision))
  udd <- unique(as.character(vec)) 
  udd <- udd[udd != intercept]
  for(i in 1:length(udd)){
    j <- i+1
    rho <- vcov(model)[1,j]/(sqrt(vcov(model)[1,1])*sqrt(vcov(model)[j,j]))
    #Standard deviation for the difference in the fixed effects
    sigma <- sqrt(vcov(model)[1,1] + vcov(model)[j,j] + 
                    2 * rho *(sqrt(vcov(model)[1,1]) *(sqrt(vcov(model)[j,j]))))
    central <- coef(summary(model))['(Intercept)','Estimate'] + 
      coef(summary(model))[paste0(name,udd[i]),'Estimate']
    ctl <- round(InvLogit(central + offset),3)
    upp <- round(InvLogit(central + offset + 1.96*sigma),3)
    low <- round(InvLogit(central + offset - 1.96*sigma),3)
    dfe[1+i,] <- c(udd[i],ctl,low,upp)
  }
  if(bfi==1){
    dfe$Blood.Feeding <- as.numeric(dfe$Blood.Feeding)
  }else{
    dfe$Mortality <- as.numeric(dfe$Mortality)
  }
  dfe$Lower_95pc_CI <- as.numeric(dfe$Lower_95pc_CI)
  dfe$Upper_95pc_CI <- as.numeric(dfe$Upper_95pc_CI)
  dfe$Arm <- gsub('_',' ',dfe$Arm)
  if(table(stringr::str_detect(dfe$Arm,'Control'))[2][[1]]==1){
    dfe <- rbind(dfe[dfe$Arm == "Control",], dfe[dfe$Arm != "Control",])
  }
  return(dfe)
}

plot_NI_OR <- function(OR, ORl, ORu, mortality = 1, NIM = 0.7, precision = 2, 
                       title = 'Candidate vs. Active Comparator'){ #washing in title??
  
  if(mortality == 1){
    
    title2 <- paste0(title,'\n - mosquito mortality')
    labl <- as.character(paste0('OR=',
                                format(round(OR,precision),nsmall=precision),' [',
                                format(round(ORl,precision),nsmall=precision),', ',
                                format(round(ORu,precision),nsmall=precision),']'))
    print(labl)
    print(paste0('NIM: ',format(round(NIM,precision),nsmall=precision)))
    print(NIM_answer(type = 'mort', vecOR = c(OR,ORl,ORu), nim = NIM))
    upper <- max(c(2,ORu + 0.35))
    pl <- ggplot() + geom_point(aes(x=1,y=OR)) + theme_classic() + xlim(c(0.94,1.06)) + 
      geom_line(aes(x = c(1,1), y = c(ORl,ORu))) + xlab('') + ylab('Odds Ratio') + 
      geom_line(aes(x=c(0.965,1.06), y = c(NIM,NIM)), color = 'grey') + #linetype = 'dashed'
      geom_hline(yintercept = 1, color = 'grey', linetype = 'dashed') +
      ylim(c(0,NA)) + theme(axis.line.x = element_blank(),
                            axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
      annotate('text', x=0.94, y = NIM, label = 'N.I.M.', color = 'grey40', hjust = 0) + 
      annotate('text', x=1, y = upper - 0.17, label = labl) + 
      geom_rect(aes(xmin=0.94, xmax=1.06, ymin=NIM, ymax=upper), fill="green4", alpha=0.35) + 
      geom_rect(aes(xmin=0.94, xmax=1.06, ymin=0, ymax = NIM), fill="red3", alpha=0.35) +
      annotate('text', x = 1.05, y = 0.99, label = 'Favours candidate >>>',#bquote('Favours candidate'~'\u2192'),
               angle = 90, alpha = 1, vjust = 1) #ylim(0,1.5) +
    pl
  }else{
    
    title2 <- paste0(title,'\n - blood feeding')
    labl <- as.character(paste0('OR=',
                                format(round(OR,precision),nsmall=precision),' [',
                                format(round(ORl,precision),nsmall=precision),', ',
                                format(round(ORu,precision),nsmall=precision),']'))
    print(labl)
    print(paste0('NIM: ',format(round(NIM,precision),nsmall=precision)))
    print(NIM_answer(type = 'bfi', vecOR = c(OR,ORl,ORu), nim = NIM))
    upper <- max(c(2,ORu + 0.35))
    pl <- ggplot() + geom_point(aes(x = 1, y = OR)) + theme_classic() + xlim(c(0.94,1.06)) + 
      geom_line(aes(x = c(1,1), y = c(ORl,ORu))) + xlab('') + ylab('Odds Ratio') + 
      #geom_hline(yintercept = 1/0.7, color = 'grey') + #linetype = 'dashed'
      geom_line(aes(x=c(0.965,1.06), y = c(NIM,NIM)), color = 'grey') + 
      geom_hline(yintercept = 1, color = 'grey', linetype = 'dashed') +
      theme(axis.line.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
      annotate('text', x=0.94, y = NIM, label = 'N.I.M.', color = 'grey40', hjust = 0) + 
      annotate('text', x=1, y = upper - 0.17, label = labl) + 
      
      geom_rect(aes(xmin=0.94, xmax=1.06, ymin=0, ymax=NIM), fill="green4", alpha=0.35) + 
      geom_rect(aes(xmin=0.94, xmax=1.06, ymin=NIM, ymax = upper), fill="red3", alpha=0.35) +
      annotate('text', x = 1.05, y = .9, label = 'Favours candidate >>>',#bquote('Favours candidate'~'\u2192'),
               angle = 270, alpha = 1, vjust = 1, hjust = 0) #+ ylim(c(0,NA))
  }
  
  return(pl + ggtitle(title2))
}

variable_NIM <- function(pc_diff = 0.07, mortality = 1, OR, ORl, ORu,
                         FIC, ymin = 0.3, ymax = 0.7, xmin = 0.0, xmax = 1.5){
  if(mortality==1){
    xx <- seq(pc_diff,0.99,0.01)
    
    or <- rep(0,length(xx))
    for(i in 1:length(xx)){
      or[i] <- ((xx[i]-pc_diff)/(1-xx[i] + pc_diff))/((xx[i])/(1-xx[i]))
    }
    df <- data.frame('mortality' = xx, 'NIM' = or)
    
    ggplot() + geom_path(data = df, aes(y = mortality, x = or)) + theme_classic() + 
      geom_point(aes(y = FIC, x=OR), size = 2.3) + 
      geom_line(aes(y = c(FIC,FIC), x = c(ORl,ORu))) + 
      xlab('Odds Ratio') + ylab('Mortality of the active comparator') + 
      geom_ribbon(data = df, aes(y = mortality, xmin = 0.02, xmax = or), fill="red3", alpha=0.35) + 
      geom_ribbon(data = df, aes(y = mortality, xmin = or, xmax = xmax), fill="green4", alpha=0.35) + 
      geom_vline(xintercept = 1, color = 'grey66', linetype = 'dashed') + 
      annotate('text', x = 1.01, y = ymin + 0.8*(ymax-ymin), label = 'Favours candidate >>>',#bquote('Favours candidate'~'\u2192'),
               alpha = 1, hjust = 0) +
      #annotate('text', x = df$NIM[5] + 0.05, y = df$mortality[5], 
      #         label = 'Non-inferiority margin',angle = 0, hjust = 0) + 
      scale_y_continuous(#breaks = c(0.4,0.5,0.6), 
        limits = c(ymin,ymax)) + 
      scale_x_continuous(breaks = c(0.25,0.5,0.75,1,1.25), limits = c(xmin,xmax)) + 
      annotate('text', x = 0.765, y = ymin + 0.05, label = 'N.I.M.',angle = 0)
  }else{
    xx <- seq(0.05,0.99 - pc_diff - 0.01,0.01)
    
    or <- rep(0,length(xx))
    for(i in 1:length(xx)){
      or[i] <- ((xx[i] + pc_diff)/(1 - xx[i] - pc_diff))/((xx[i])/(1-xx[i]))
    }
    df <- data.frame('bf' = xx, 'NIM' = or)
    ggplot() + geom_path(data = df, aes(y = bf, x = or)) + theme_classic() + 
      geom_point(aes(y = FIC, x=OR), size = 2.3) + 
      geom_line(aes(y = c(FIC,FIC), x = c(ORl,ORu))) + 
      xlab('Odds Ratio') + ylab('Blood-feeding proportion for the active comparator') + 
      geom_ribbon(data = df, aes(y = bf, xmin = 0.02, xmax = or), fill="green4", alpha=0.35) + 
      geom_ribbon(data = df, aes(y = bf, xmin = or, xmax = xmax), fill="red3", alpha=0.35) + 
      geom_vline(xintercept = 1, color = 'grey66', linetype = 'dashed') + 
      annotate('text', x = 1.01, y = ymin + 0.8*(ymax-ymin), label = '<<< Favours candidate',#bquote('Favours candidate'~'\u2192'),
               alpha = 1, hjust = 1) +
      #annotate('text', x = df$NIM[5] + 0.05, y = df$bf[5], 
      #         label = 'Non-inferiority margin',angle = 0, hjust = 0) + 
      scale_y_continuous(#breaks = c(0.4,0.5,0.6), 
        limits = c(ymin,ymax)) + 
      scale_x_continuous(breaks = c(0.25,0.5,0.75,1,1.25), limits = c(xmin,xmax)) + 
      annotate('text', x = 1.45, y = ymin + 0.05, label = 'N.I.M.',angle = 0)
  }
}

new_median_FE <- function(model, FE = c('hut','sleeper','day')){
  l <- length(FE)
  ofs <- 0
  count <- 0
  if(l==3){
    xx <- as.data.frame(tidyr::crossing(names(model$coefficients)[grep(FE[1],names(model$coefficients))],
                                        names(model$coefficients)[grep(FE[2],names(model$coefficients))], 
                                        names(model$coefficients)[grep(FE[3],names(model$coefficients))]))
    colnames(xx) <- FE
    #print(head(xx))
    
    stor <- rep(0,dim(xx)[1])
    for(i in 1:(dim(xx)[1])){
      stor[i] <- coef(summary(model))[xx[i,1],'Estimate'] +
        coef(summary(model))[xx[i,2],'Estimate'] + coef(summary(model))[xx[i,3],'Estimate']
      #print(i)
      
    }
    ofs <- median(stor,na.rm = T)
    count <- count + 1
  }
  
  if(l==2){
    xx <- as.data.frame(tidyr::crossing(names(model$coefficients)[grep(FE[1],names(model$coefficients))],
                                        names(model$coefficients)[grep(FE[2],names(model$coefficients))]))
    colnames(xx) <- FE
    #print(head(xx))
    
    stor <- rep(0,dim(xx)[1])
    for(i in 1:(dim(xx)[1])){
      stor[i] <- coef(summary(model))[xx[i,1],'Estimate'] +
        coef(summary(model))[xx[i,2],'Estimate']# + coef(summary(model))[xx[i,3],'Estimate']
      #print(i)
      
    }
    ofs <- median(stor,na.rm = T)
    count <- count + 1
  }
  
  if(l==1){
    xx <- as.data.frame(names(model$coefficients)[grep(FE[1],names(model$coefficients))])
    colnames(xx) <- FE[1]
    #print(head(xx))
    
    stor <- rep(0,dim(xx)[1])
    for(i in 1:(dim(xx)[1])){
      stor[i] <- coef(summary(model))[xx[i,1],'Estimate'] #+
        #coef(summary(model))[xx[i,2],'Estimate']# + coef(summary(model))[xx[i,3],'Estimate']
      #print(i)
      
    }
    ofs <- median(stor,na.rm = T)
    count <- count + 1
  }
  if(count != 1){
    print('Check number of fixed effects entered (must be 1, 2, or 3)')
  }
  return(ofs)
}

tidy_blf_FE <- function(data, model_fit, vec, intercept, name, first_cat = "Control", model_fit_blf, offset = c(0,0)){
  
  dfx <- mFE(model = model_fit, vec = vec, name = name, intercept = intercept, bfi = 0, offset = offset[1])#mortality_summary(model_fit)
  dfx$summ <- paste0(100*dfx$Mortality,'%, [',100*dfx$Lower_95pc_CI,',',100*dfx$Upper_95pc_CI,']')
  data$Arm <- vec 
  
  dfx$count <- NA
  dfx$av <- NA
  for(i in 1:(dim(dfx)[1])){
    tr <- gsub(' ','_',dfx$Arm[i]) 

    aux <- sum(data[data$Arm==tr,]$total) 
    dfx$count[i] <- aux
    tr2 <- table(data$Arm)[tr][[1]]
    dfx$av[i] <- round(aux/tr2,1)
  }
  dfx
  dfx2 <- rbind(dfx[dfx$Arm == first_cat,],dfx[dfx$Arm!=first_cat,])
  dfx2
  dfx3 <- dplyr::select(dfx2, c('Arm','count','av','summ'))
  colnames(dfx3) <- c('Arm','Total mosquitoes','Mosquitoes per hut per night','Mortality [95% CI]')
  
  dfy <- mFE(model = model_fit_blf, vec = vec, name, intercept = intercept, bfi = 1, offset = offset[2])
  dfy$summ <- paste0(100*dfy$Blood.Feeding,'%, [',100*dfy$Lower_95pc_CI,',',100*dfy$Upper_95pc_CI,']')

  dfy$bfi <- round(100*(1-dfy$Blood.Feeding/dfy[dfy$Arm==first_cat,]$Blood.Feeding),3)

  dfy <- dfy[,c(1,5,6)]
  colnames(dfy) <-c('Arm','Blood Feeding [95% CI]','B.F.I.(%)')
  
  #dfy$Arm[1] <- intercept
  dfx4 <- merge(dfx3,dfy,by='Arm')
  dfx4$`B.F.I.(%)`<-ifelse(dfx4$`B.F.I.(%)`<0,'<0',dfx4$`B.F.I.(%)`)
  dfx5 <- rbind(dfx4[dfx4$Arm == first_cat,],dfx4[dfx4$Arm!=first_cat,])
  dfx5$`B.F.I.(%)`[1] <- '-'
  
  return(dfx5)
}

