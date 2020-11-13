#Input data and WBModel adapted to France

setwd("X:/Luis/Model")
library(dplyr)

rpert <- function(n,a,l,b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
  }
  a+(b-a)*rbeta(n,v,w)
}

WBModel <- function( MaxIterations    = 1,
                     MaxDays          = 365*5,
                     Detailed         = TRUE,
                     InitialOccCells  = 50,
                     MaxMortProbAd    = 1 - ((0.4)^(1/365)),          # Minimum Survival Prob for adults Lange et al. (2015), we recalculate per day 
                     MaxMortProbSAd   = 1 - ((0.4)^(1/365)),          # Minimum Survival Prob for Subadults Lange et al. (2015) we recalculate per day 
                     MaxMortProbPig   = 1 - ((0.1)^(1/365)),          # Minimum Survival Prob for piglets Lange et al. (2015) we recalculate per day
                     SurvivalProbAdF  = 0.9988^365,                   # Yearly Survival Prob for adults Kueling et al. (2013)    
                     SurvivalProbAdM  = 0.9986^365,                   # Yearly Survival Prob for adults Kueling et al. (2013) 
                     SurvivalProbSAdF = 0.9988^365,                   # Yearly Survival Prob for subadults (we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)
                     SurvivalProbSAdM = 0.9986^365,                   # Yearly Survival Prob for subadults(we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)    
                     SurvivalProbPigF = 0.50,                         # Yearly Survival Prob for piglets Lange et al. (2015) 
                     SurvivalProbPigM = 0.45,                         # Yearly Survival Prob for piglets Lange et al. (2015), reduced to consider higher survival of females
                     AgeProbab        = c(0.39,0.24,0.15,0.09,0.06,0.03,0.02,0.01,0.01), # Age probabilities from Lange et al. (2015)
                     ## BreedSizeFactor  =                            # A factor to allow groups to breed as soon as the group size is less than the maximum allowed size.
                     ProbSplitSDS     =  1/7,  #cumsum(rep(1/7,7))    # Probability of short term split if all conditions are satisfyed for a group to split per day during the weekof spliting
                     ProbSplitMSA     = 1/(6*7), #cumsum(rep(1/(6*7),6*7))     # Probability of male spliting per group per day of the week (7days) of the total spliting period (6 weeks)
                     SettlingProb     = c(0.75,1),                    # probability to settle for a male in a habitat cell if it was only good habitate (0.5) or good habitate with females (1.0)
                     FemPropGroup     = mean(c(14.7,18.1,12.8,14.6)), # Proportion of females within a group Merta et al. (2015)  
                     MalPropGroup     = mean(c(9.2,9,5.7,7.9)),       # Proportion of males within a group Merta et al. (2015)
                     SubAPropGroup    = mean(c(33.4,34.3,22.8,24.1)), # Proportion of subadults within a group Merta et al. (2015)
                     PigletPropGroup  = mean(c(42.7,38.6,58.7,53.4)), # Proportion of piglets within a group Merta et al. (2015)
                     
                     ### ASF related parameters
                     WGProbInf        = 1-(1-0.05)^(1/7),             # Direct contact Lange et al. (2015EFSA) and Lange and Thulke (2017)
                     CarcProbInf      = 1-(1-0.2)^(1/7),              # Carcass Lange et al. (2015EFSA) and Lange and Thulke (2017)
                     ProbAccesCarc    = 0.9,                          # Lange and Thulke (2017)
                     ProbDeathInf     = 0.95,                         # Blome (2012)  
                     CarcDaySur       = 28,                           # Lange and Thulke (2017)
                     DieInConEdge     = 0.8,                          # Lange et al. (2015EFSA)
                     TimeSeedInf      = 0,
                     NumbGSeedInf     = 0,
                     MinBCap          = 1,
                     ModeBCap         = 3,
                     MaxBCap          = 6,
                     AREA             = 1,
                     
                     RepProbList      = list(Week=1:26,
                                             Prob=c(0,0,0,0.01, 0.05,0.1,0.2,0.23, 0.28,0.325,0.35,0.325, 0.3,0.28,0.25,0.22, 0.18,0.17,0.16,0.15,
                                                    0.14,0.13,0.12,0.08, 0.04,0.0)),
                     NumOfSprProbList = list(Number=0:8,
                                             Prob=c(0.01,0.07,0.16,0.25,0.25,0.16,0.07,0.02,0.001)),
                     DirectionList    =  list(c(6,7,8), c(8,10,13), c(11,12,13), c(6,9,11)), ## North, ## East, ## South, and ## West
                     
                     HabitatProb      = c(0.75,0.25),         # The probability to select a new pixel while male walking (we model that they prefer to walk in a habitat cell) 
                     FileWildBoarMat = "Inputs/EA_9km/Input_EA9.csv",
                     runID            = "EA9NI"
){
  
  WBMat <- read.table(FileWildBoarMat, sep=",", header = T)
  WBMat <- cbind(WBMat,0)
  WBMat <- WBMat[,-1]
  colnames(WBMat) <- c("ID", "Lon", "Lat", "Forest_cover", "Habitat", "Nb1", "Nb2",
                       "Nb3", "Nb4", "Nb5", "Nb6", "Nb7", "Nb8", "Colors", "Breeding")
  
  
  DayOutInfMat  <- matrix(0,ncol=MaxIterations,nrow=MaxDays)
  DayOutPopMat  <- matrix(0,ncol=MaxIterations,nrow=MaxDays)
  DayOutAniMat  <- matrix(0,ncol=MaxIterations,nrow=MaxDays)
  
  NewInfGroups  <- matrix(numeric(0),ncol=3)
  NewInfAnimals <- matrix(numeric(0),ncol=3)
  
  EpDuration           <- rep(0,MaxIterations)
  FreqRelapse          <- rep(0,MaxIterations)
  cumDeath             <- rep(0,MaxIterations)
  Fpopcount            <- matrix(0, MaxDays, nrow(WBMat))
  Finfcount            <- matrix(0, MaxDays, nrow(WBMat))
  Fcarcount            <- matrix(0, MaxDays, nrow(WBMat))
  Fimmunecount         <- matrix(0, MaxDays, nrow(WBMat))
  
  for(Iter in 1:MaxIterations){
    set.seed(Iter)
    ## Initialize the population matrix
    PopMatWB <- matrix(0,ncol=15,nrow=50000)
    ## Add the habitat capacity
    WBMat[WBMat[,5]==1,15] <- ceiling(rpert(sum(WBMat[,5]==1),MinBCap,ModeBCap,MaxBCap))

    TMPHP <- WBMat[,5]==1
                    
    homePixelsAll   <- sample(WBMat[TMPHP,1],InitialOccCells)
  
    InitSizeMat <- matrix(0,ncol=15)
    
    for(i in 1:InitialOccCells){
      ## The initial distribution of the groups was based on Merta et al. (2015). They were calculated relative to the average
      ## Female percentage. for instance males were average % males divided by average % females. 
      females       <- round(runif(1,3,4))
      males         <- round(females*MalPropGroup/FemPropGroup)
      subAdults     <- round(females*SubAPropGroup/FemPropGroup)
      Piglets       <- round(females*PigletPropGroup/FemPropGroup)
      GroupID       <- rep(i,sum(females,males,subAdults,Piglets))
      Sex           <- c(rep(1,females),rep(0,males),rbinom(subAdults,1,prob=0.5),rbinom(Piglets,1,prob=0.5))
      AgeCat        <- c(rep(3,females),rep(3,males),rep(2,subAdults),rep(1,Piglets))
      Dam           <- c(rep(0,sum(females,males)),sample(1:females,subAdults,rep=T),sample(1:females,Piglets,rep=T))
      # Females deliver between Jan and June, so the age of the piglets will be between 183 and 365 days.
      # The same thing for sub-adults but with 365 days extra.
      tmpAgeSubA    <- sample((183:365)+365,females,rep=T)
      tmpAgePig     <- sample((183:365),females,rep=T)
      Age           <- c(sample(2:10,size=sum(females,males),T,prob=AgeProbab)*365,tmpAgeSubA[Dam[(females+males+1):(females+males+subAdults)]],
                         tmpAgePig[Dam[(females+males+subAdults+1):(females+males+subAdults+Piglets)]])
      # Adjust the Dam number to fit the actual ID of the DAMs
      Dam[(females+males+1):(females+males+subAdults)] <- Dam[(females+males+1):(females+males+subAdults)] + (dim(InitSizeMat)[1]) - 1 # -1 because we have extra row at start
      Dam[(females+males+subAdults+1):(females+males+subAdults+Piglets)] <- Dam[(females+males+subAdults+1):(females+males+subAdults+Piglets)] + (dim(InitSizeMat)[1]) - 1
      Breed         <- rep(0,sum(females,males,subAdults,Piglets))
      HomePixel     <- rep(homePixelsAll[i],sum(females,males,subAdults,Piglets))
      CurrPixel     <- HomePixel
      infectStatus  <- rep(0,sum(females,males,subAdults,Piglets))
      SplitStatus   <- rep(0,sum(females,males,subAdults,Piglets))
      SplitMale     <- rep(0,sum(females,males,subAdults,Piglets))
      TimeDeath     <- rep(0,sum(females,males,subAdults,Piglets))
      TimeToInfect  <- rep(0,sum(females,males,subAdults,Piglets))
      TimeToDeath   <- rep(0,sum(females,males,subAdults,Piglets))
      IDs           <- (max(InitSizeMat[,1])+1):(max(InitSizeMat[,1])+sum(females,males,subAdults,Piglets))
      
      InitMatWBPop  <- cbind(IDs,GroupID,Sex,AgeCat,Age,Breed,HomePixel,CurrPixel,infectStatus,SplitStatus,Dam,SplitMale,TimeDeath,TimeToInfect,TimeToDeath)
      InitSizeMat   <- rbind(InitSizeMat,InitMatWBPop)
      WBMat[unique(HomePixel),15] <- females
    }
    InitSizeMat <- InitSizeMat[-1,] # This is just to remove the first line that includes zeros, necessary for initialization to keep track of IDs.
    PopMatWB[1:dim(InitSizeMat)[1],] <- InitSizeMat
    
    GroupsToSplit        <- matrix(numeric(0),ncol=4)  
    
    colnames(PopMatWB) <- c("IDs", "Group_ID", "Sex", "Age_Cat", "Age_days", "Breed",
                            "Home_pixel", "Current_pixel", "Infect_status", "Split_status", 
                            "Dam", "Split_male", "Time_death","Time_to_infect", "Time_to_death")
    
    ## Store the groups that have split this year
    SplittedGroups  <- numeric(0)
    cumDeathPar     <- 0
    Year            <- 1
    gTime           <- 0
    Criteria        <- TRUE
    TMPOutYesterday <- FALSE
    OnlyOnce        <- TRUE
    
    while(Criteria & (gTime < MaxDays)){
      
     #if(gTime<=TimeSeedInf) Criteria <- TRUE else Criteria <- any(PopMatWB[,9]%in%1:3)
      
      MatSizeExp <- sum(PopMatWB[,1]==0)
      if(MatSizeExp<1000){
        ExtraRows <- matrix(rep(0,15000),ncol=15)
        PopMatWB <- rbind(PopMatWB,ExtraRows)
      } 

      gTime <- gTime + 1
      PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
      
      if(gTime > 365) Year <- ceiling(gTime/365)  ## Count the year 
      PopMatWB[PopMatWB[,5]>0,5] <- PopMatWB[PopMatWB[,5]>0,5] + 1
      PopMatWB[PopMatWB[,5]>=365,4]     <- 2
      PopMatWB[PopMatWB[,5]>=(365*2),4] <- 3
      
      ### initiate variables and matrix that have to be initiated on yearly basis
      if(gTime %in% c(1, (366 * 1:(MaxDays/365)))){
        
        
        ### Mortality probabilities 
        
        
        ProbMortFact <- rnorm(1,1,0.1)
        ProbMortAdF <- 1 - (ProbMortFact*SurvivalProbAdF)^(1/365)
        ProbMortAdF[ProbMortAdF>MaxMortProbAd]   <- MaxMortProbAd
        ProbMortAdM <- 1 - (ProbMortFact*SurvivalProbAdM)^(1/365)
        ProbMortAdM[ProbMortAdM>MaxMortProbAd]   <- MaxMortProbAd
        ProbMortSAdF <- 1 - (ProbMortFact*SurvivalProbSAdF)^(1/365)
        ProbMortSAdF[ProbMortSAdF>MaxMortProbSAd] <- MaxMortProbSAd
        ProbMortSAdM <- 1 - (ProbMortFact*SurvivalProbSAdM)^(1/365)
        ProbMortSAdM[ProbMortSAdM>MaxMortProbSAd] <- MaxMortProbSAd
        ProbMortPigF <- 1 - (ProbMortFact*SurvivalProbPigF)^(1/365)
        ProbMortPigF[ProbMortPigF>MaxMortProbPig] <- MaxMortProbPig
        ProbMortPigM <- 1 - (ProbMortFact*SurvivalProbPigM)^(1/365)
        ProbMortPigM[ProbMortPigM>MaxMortProbPig] <- MaxMortProbPig
        
        
        ### Reproduction
        
        
        PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        AnimalsWInG     <- unlist(sapply(unique(PopMatWB[,2]),function(x) 1:sum(PopMatWB[,2]==x)))
        BreedCapCells   <- c(WBMat[PopMatWB[,7],15],rep(0,sum(PopMatWB[,7]==0)))
        
        #GroupSize       <- cbind(sort(unique(PopMatWB[,2])),tapply(PopMatWB[,1],PopMatWB[,2],length))
        # Animals are allowed to breed when they are females and the group size is less than the maximum allowed size (size that can be maintained by the pixel habitat)
        # and the animals above the breeding capacity based on age.
        AniAllowBreed   <- AnimalsWInG<=BreedCapCells & PopMatWB[,3]==1 &PopMatWB[,9]<2
      
        PopMatWB[AniAllowBreed,6]    <- (sample(RepProbList$Week,sum(AniAllowBreed),rep=T,prob=RepProbList$Prob)*7) + ((Year-1)*365)
       
        
        ### Reset Male Splitting
        
        
        PopMatWB[,12]          <- 0 #!PopMatWB[,2]%in%GroupsToSplit
        IndexRem               <- which(PopMatWB[,2]%in%GroupsToSplit[,1])
        PopMatWB[IndexRem,]    <- 0
       
        PopMatWB               <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        cumDeathPar            <- cumDeathPar + length(IndexRem)
        GroupsToSplit          <- matrix(numeric(0),ncol=4) 
      }
      
      
      ### Model Mortality 
      
      
      ## We exculde infectious and died animals from disease in order to allow the carcass to persist
      AdultsFToDie    <- which(PopMatWB[,4]==3&PopMatWB[,3]==1&!PopMatWB[,9]%in%2:3)
      AdultsFToDie    <- AdultsFToDie[rbinom(length(AdultsFToDie),1,ProbMortAdF)==1]
      AdultsMToDie    <- which(PopMatWB[,4]==3&PopMatWB[,3]==0&!PopMatWB[,9]%in%2:3)
      AdultsMToDie    <- AdultsMToDie[rbinom(length(AdultsMToDie),1,ProbMortAdM)==1]
      
      SubAdultsFToDie <- which(PopMatWB[,4]==2&PopMatWB[,3]==0&!PopMatWB[,9]%in%2:3)
      SubAdultsFToDie <- SubAdultsFToDie[rbinom(length(SubAdultsFToDie),1,ProbMortSAdF)==1]
      SubAdultsMToDie <- which(PopMatWB[,4]==2&PopMatWB[,3]==1&!PopMatWB[,9]%in%2:3)
      SubAdultsMToDie <- SubAdultsMToDie[rbinom(length(SubAdultsMToDie),1,ProbMortSAdM)==1]
      
      PigsFToDie      <- which(PopMatWB[,4]==1&PopMatWB[,3]==0&!PopMatWB[,9]%in%2:3)
      PigsFToDie      <- PigsFToDie[rbinom(length(PigsFToDie),1,ProbMortPigF)==1]
      PigsMToDie      <- which(PopMatWB[,4]==1&PopMatWB[,3]==1&!PopMatWB[,9]%in%2:3)
      PigsMToDie      <- PigsMToDie[rbinom(length(PigsMToDie),1,ProbMortPigM)==1]
      
      TooOldAni       <- which(PopMatWB[,5]==(11*365))
      
      # Piglets that are less than 8 weeks and do not have a mother will die.
      # This is based on Petersen (1999) that showed that the majority of piglets would graze on their own when they reach 8 weeks 
      TooYoungNoMoth  <- which(PopMatWB[,5]< (8*7) & !(PopMatWB[,11]%in%PopMatWB[,1])&!PopMatWB[,9]%in%2:3) 
      
      ToDieNormal <- c(AdultsFToDie,AdultsMToDie,SubAdultsFToDie,SubAdultsMToDie,PigsFToDie,PigsMToDie,TooYoungNoMoth,TooOldAni)
      
      if(length(ToDieNormal)>0){
        PopMatWB[ToDieNormal,] <- 0
        PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        
        cumDeathPar <- cumDeathPar+length(ToDieNormal)
      }
      
      
      ### Model Reproduction
      
      
      ## On daily basis, from Jan to end June, check if there are animals to deliver and make them deliver
        if(gTime%in%PopMatWB[,6]){
        DelIndex        <- which(PopMatWB[,6]%in%gTime&PopMatWB[,9]!=3)
        NumOfSpring     <- sample(NumOfSprProbList$Number,length(DelIndex),rep=T,prob=NumOfSprProbList$Prob)
        DelIndex        <- DelIndex[NumOfSpring>0]
        NumOfSpring     <- NumOfSpring[NumOfSpring>0]
        if(length(NumOfSpring)>0){
          IDNew              <- (max(PopMatWB[,1])+1):(max(PopMatWB[,1])+sum(NumOfSpring))
          IndexLocation      <- (sum(PopMatWB[,1] > 0)+1):(sum(PopMatWB[,1] > 0)+length(IDNew))
          PopMatWB[IndexLocation,1]  <- IDNew
          PopMatWB[IndexLocation,2]  <- rep(PopMatWB[DelIndex,2],NumOfSpring)
          PopMatWB[IndexLocation,3]  <- rbinom(sum(NumOfSpring),1,0.5)
          PopMatWB[IndexLocation,4]  <- 1
          PopMatWB[IndexLocation,5]  <- 1
          PopMatWB[IndexLocation,6]  <- 0
          PopMatWB[IndexLocation,7]  <- rep(PopMatWB[DelIndex,7],NumOfSpring)
          PopMatWB[IndexLocation,8]  <- rep(PopMatWB[DelIndex,7],NumOfSpring)
          # Remember to make infection status dependent on Dam's infection status
          PopMatWB[IndexLocation,9]  <- 0
          PopMatWB[IndexLocation,10] <- 0
          PopMatWB[IndexLocation,11] <- rep(PopMatWB[DelIndex,1],NumOfSpring)
          PopMatWB[IndexLocation,12] <- 0
          PopMatWB[IndexLocation,13] <- 0
          PopMatWB[IndexLocation,14] <- 0
          PopMatWB[IndexLocation,15] <- 0
          
          PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
          
        }
      }
      
      
      ### Model females groups splitting
      
      
      ## This code checks whether splitting of subadult females may happen. 
      ## It occurs only once in week 28 Kramer-Schadt et al. (2009)
      if(gTime%in%((28*7+((Year-1)*365))+0:6)){
        ## Calculate the day in the week
        #DayInWeekSSplit <- (1-(ceiling((gTime-(365*(Year-1)))/7) - ((gTime-(365*(Year-1)))/7)))*7
        PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        AnimalsWInG     <- unlist(sapply(unique(PopMatWB[,2]),function(x) 1:sum(PopMatWB[,2]==x)))
        BreedCapCells   <- c(WBMat[PopMatWB[,7],15],rep(0,sum(PopMatWB[,7]==0)))
        
        ## Identify the animals within the groups that will split
        ## Females subadults > breading capacity, did not split before and are not sick
        GroupSplitNum   <- cbind(sort(unique(PopMatWB[,2])),(tapply((AnimalsWInG>BreedCapCells & PopMatWB[,4]==2 & PopMatWB[,3]==1 & PopMatWB[,10]==0& PopMatWB[,9]<2),  PopMatWB[,2], sum)))
        GroupSplitNum   <- GroupSplitNum[GroupSplitNum[,2]>=2,,drop=FALSE]
        
        ## Idetify groups where there are only males
        ## To allow females to join groups where there are only males.
        PixelMalesOnly  <- cbind(sort(unique(PopMatWB[,2])),(tapply((PopMatWB[,3]==0),PopMatWB[,2], all)))
        PixelMalesOnly  <- PixelMalesOnly[PixelMalesOnly[,2]==1,,drop=FALSE]
        
        GroupsSplitShD  <- sapply(PopMatWB[match(GroupSplitNum[,1],PopMatWB[,2]),7],function(x) {
          Distance <- sqrt((WBMat[x,2]-WBMat[,2])^2 + (WBMat[x,3]-WBMat[,3])^2)/1000
          tmp1     <- which(Distance<=9)
          any(WBMat[tmp1,5]==1 & !(tmp1%in%PopMatWB[!PopMatWB[,2]%in%PixelMalesOnly[,1],7]))})
        GroupSplitNumSD <- cbind(GroupSplitNum,GroupsSplitShD)
        GroupSplitNumSD <- GroupSplitNumSD[GroupSplitNumSD[,3]==1,,drop=FALSE]
        
        ## Groups to split, may split during any day of the week. No need for a probability here the group splits only once a year.
        #if(dim(GroupSplitNumSD)[1]>0) GroupSplitNumSD <- GroupSplitNumSD[rbinom(dim(GroupSplitNumSD)[1],1,prob=ProbSplitSDS)==1,,drop=FALSE]
        PopMatWB[PopMatWB[,2]%in%GroupSplitNumSD[,1],10] <- 1
        
        ## Make short term split to happen
        if(dim(GroupSplitNumSD)[1]>0){
          OriginPixel     <- PopMatWB[match(GroupSplitNumSD[,1],PopMatWB[,2]),7]
          
          TargetPixel <- numeric(0)
          for(xx in OriginPixel){
            Distance <- sqrt((WBMat[xx,2]-WBMat[,2])^2 + (WBMat[xx,3]-WBMat[,3])^2)/1000
            tmp1      <- which(Distance<=9)
            tmp2      <- tmp1[WBMat[tmp1,5]==1 & !(tmp1%in%PopMatWB[,7])]
            if(length(tmp2)>1)  tmp3 <- sample(tmp2,1)
            if(length(tmp2)==1) tmp3 <- tmp2
            if(length(tmp2)==0) tmp3 <- 0
            TargetPixel <- c(TargetPixel,tmp3)
          }
          
          GroupSplitNumSD <- GroupSplitNumSD[TargetPixel>0,,drop=FALSE]
          OriginPixel     <- OriginPixel[TargetPixel>0]
          TargetPixel     <- TargetPixel[TargetPixel>0]
          GroupSplitNumSD <- cbind(GroupSplitNumSD,TargetPixel) 
          
          ## This part of the code to allow movement
          if(dim(GroupSplitNumSD)[1]>0){
            
            ## A list to keep track for the pixels where the pigs have been. 
            ## Make sure that this information is exported on daily basis, because the 
            ## list will be re-initiate every day splitting may happen.
            PixelsMoved <- as.list(matrix(0,ncol=length(TargetPixel)))
            
            
            CurrentPos <- OriginPixel
            for(i in 1:length(TargetPixel)){
              Trail <- 0
              prevEdge <- CurrentPos[i]
              while(CurrentPos[i]!=TargetPixel[i] & Trail < 10){
                Trail      <- Trail + 1
                Edges      <- unlist(WBMat[CurrentPos[i],6:13])
                Edges      <- Edges[Edges>0&Edges!=prevEdge]
                if(length(Edges)>0){
                  DistEdges  <- sqrt((WBMat[TargetPixel[i],2]-WBMat[Edges,2])^2 + (WBMat[TargetPixel[i],3]-WBMat[Edges,3])^2)/1000
                  prevEdge   <- CurrentPos[i]
                  NewPosition<- Edges[DistEdges==min(DistEdges)]
                  if(length(NewPosition)>1) NewPosition<-sample(NewPosition,1)
                  CurrentPos[i] <- NewPosition
                  # Here we keep track of the pixels where the pigs moved to.
                  PixelsMoved[[i]]<- c(PixelsMoved[[i]],NewPosition)
                  #print(c(i,OriginPixel[i],prevEdge,Trail,CurrentPos[i],TargetPixel[i]))
                }
              }
            }
            
            ## Here we assume that groups that did not find the way, did not actually split as the edges was not connecting.
            IndexNotSplit   <- which(CurrentPos!=TargetPixel)
            GroupSplitNumSD <- GroupSplitNumSD[-IndexNotSplit,,drop=FALSE]
            
            ## Make females from different groups moving to a new pixel to form a new group.
            if(dim(GroupSplitNumSD)[1]>0){
              newGroupIDs <- (max(PopMatWB[,2])+1):(max(PopMatWB[,2])+dim(GroupSplitNumSD)[1])
              if(sum(duplicated(GroupSplitNumSD[,4]))>0){
                for(l in unique(GroupSplitNumSD[,4])){
                  TEMP  <- which(GroupSplitNumSD[,4]==l)
                  if(length(TEMP)>1){
                    newGroupIDs[TEMP] <- newGroupIDs[TEMP[1]]
                  }
                }
              }
              
              names(PixelsMoved) <- GroupSplitNumSD[,1]
              for(b in 1:dim(GroupSplitNumSD)[1]){
                
                ## If there is a male group already in the pixel and is not splitting, then make them one group with the new females.
                ## Notice the code above allows them to come into a new pixel only if it is empty or there are only males. 
                ## Dont worry tat IndexMalComb1 is matching with all PopMatWB[,7]
                IndexMalComb1  <- GroupSplitNumSD[b,4]%in%PopMatWB[,7]
                if(IndexMalComb1){
                  IndexMalComb2  <- unique(PopMatWB[PopMatWB[,7]==GroupSplitNumSD[b,4],2])%in%GroupsToSplit[,1]
                  if(!IndexMalComb2){
                    AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & PopMatWB[,3]==1 & PopMatWB[,4]==2 & PopMatWB[,9]<2 & PopMatWB[,2]%in%GroupSplitNumSD[b,1]
                    MalesInPixel     <- PopMatWB[,3]==0 & PopMatWB[,9]!=3 & PopMatWB[,7] == GroupSplitNumSD[b,4] & !PopMatWB[,2]%in% GroupsToSplit[,1]
                    newGroupIDsAn    <- rep(newGroupIDs[b],sum(AnimalsCanSplitF)+sum(MalesInPixel)) 
                    NewHomePixel     <- rep(GroupSplitNumSD[b,4],sum(AnimalsCanSplitF)+sum(MalesInPixel))
                    PopMatWB[AnimalsCanSplitF|MalesInPixel,2]  <- newGroupIDsAn
                    PopMatWB[AnimalsCanSplitF|MalesInPixel,7]  <- NewHomePixel
                    PopMatWB[AnimalsCanSplitF|MalesInPixel,8]  <- NewHomePixel
                  }
                  if(IndexMalComb2){
                    AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & PopMatWB[,3]==1 & PopMatWB[,4]==2 & PopMatWB[,9]<2 & PopMatWB[,2]%in%GroupSplitNumSD[b,1]
                    newGroupIDsAn    <- rep(newGroupIDs[b],sum(AnimalsCanSplitF)) 
                    NewHomePixel     <- rep(GroupSplitNumSD[b,4],sum(AnimalsCanSplitF))
                    PopMatWB[AnimalsCanSplitF,2]  <- newGroupIDsAn
                    PopMatWB[AnimalsCanSplitF,7]  <- NewHomePixel
                    PopMatWB[AnimalsCanSplitF,8]  <- NewHomePixel
                  }   
                }
                if(!IndexMalComb1){
                  AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & PopMatWB[,3]==1 & PopMatWB[,4]==2 & PopMatWB[,9]<2 & PopMatWB[,2]%in%GroupSplitNumSD[b,1]
                  newGroupIDsAn    <- rep(newGroupIDs[b],sum(AnimalsCanSplitF)) 
                  NewHomePixel     <- rep(GroupSplitNumSD[b,4],sum(AnimalsCanSplitF))
                  PopMatWB[AnimalsCanSplitF,2]  <- newGroupIDsAn
                  PopMatWB[AnimalsCanSplitF,7]  <- NewHomePixel
                  PopMatWB[AnimalsCanSplitF,8]  <- NewHomePixel
                }   
              }
            }
          } 
        }
      }#Closes short distance splitting
      
      ## After the end of female splitting, we reset female splitting.
      if(gTime%in%((28*7+((Year-1)*365))+7)){
        PopMatWB[,10]    <- 0
        SplittedGroups   <- numeric(0)
        PixelsMoved      <- as.list(matrix(0,ncol=1))
      }
      
      
      ### Model Male splitting
      
      
      ## First we determine the period of spliting of males, as defined in Lange et al., (2012). 
      ## Males may find a pixel to live in otherwise they will keep Wondering around 
      ## until they either die or find a place to live in
      if(gTime%in%((25*7+((Year-1)*365)):(30*7+((Year-1)*365)))){
        ## Sort the matrix
        PopMatWB    <-  PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        ## Define the males as adults and subadults
        ## Identify the groups where subadult males may split
        ## Male subadults from groups where males have not yet split this year.
        #GroupSplittedThisYear <- cbind(sort(unique(PopMatWB[,2])),tapply(PopMatWB[,12],PopMatWB[,2],sum))
        #GroupSplittedThisYear <- GroupSplittedThisYear[GroupSplittedThisYear[,2]>0,,drop=FALSE]
        
        ## Calculate nuber of adult males
        NumbAdultMInG <-  cbind(sort(unique(PopMatWB[,2])),tapply(PopMatWB[,3]==0&PopMatWB[,4]==3,PopMatWB[,2],sum))
        ## No need for subadult mails to split if there is need for them in the group
        tmpIndexSG    <- NumbAdultMInG[match(PopMatWB[,2],NumbAdultMInG[,1]),2]
        ## Groups that may split must not be a new splilt groups and have more than 2 adult males
        GroupSubAM    <- unique(PopMatWB[!(PopMatWB[,2]%in%SplittedGroups) & PopMatWB[,3]==0 & PopMatWB[,4]==2 & PopMatWB[,12]==0 & tmpIndexSG>2,2])
        if(length(GroupSubAM)>0) 
          GroupSubAM <- GroupSubAM[rbinom(length(GroupSubAM),1,ProbSplitMSA)==1]
        if(length(GroupSubAM)>0){
          
          ## Identify the subadault males that may split.
          ## We must ensure that the groups have adult males before the subadults are split.
          ## No need for the subadults to find new areas if they can nurish and breed in their home
          #TEMPMat <- cbind(sort(unique(PopMatWB[,2])),tapply((PopMatWB[,3]==0 & PopMatWB[,4]==3 & PopMatWB[,9] < 2),  PopMatWB[,2], sum))
          tmpSAM <- which(PopMatWB[,2]%in%GroupSubAM & PopMatWB[,3]==0 & PopMatWB[,4]==2 & PopMatWB[,9] < 2)# & TEMPMat[match(PopMatWB[,2],TEMPMat[,1]),2]>=2)
          ## Check which ones that would split
          if(length(tmpSAM)>0) {
            ## Make subadults to split based on probability from 25 to 75% (Truve et al., 2004).
            tmpSAM           <-tmpSAM[rbinom(length(tmpSAM),1,prob=runif(length(tmpSAM),0.25,0.75))==1]
            AinmFromGToSplit <- cbind(sort(unique(PopMatWB[tmpSAM,2])),tapply(tmpSAM,PopMatWB[tmpSAM,2],length))
            AinmFromGToSplit <- AinmFromGToSplit[AinmFromGToSplit[,2]>1,,drop=FALSE]
            tmpSAM           <-tmpSAM[PopMatWB[tmpSAM,2] %in% AinmFromGToSplit[,1]]
            
            if(length(tmpSAM)>0) {# do not worry, if the dimension is > 0 then there will be at least 2 animals from one group to split ;-) check above.
              PopMatWB[tmpSAM,12] <- 1
              ## Store the group number of the groups that have splitted
              SplittedGroups <- c(SplittedGroups,unique(PopMatWB[tmpSAM,2]))
              ## Number of spliting males per group
              AinmFromGToSplit   <- cbind(AinmFromGToSplit,round(runif(dim(AinmFromGToSplit)[1],1,4)),unique(PopMatWB[tmpSAM,8]))
              ## The third column, informs about the direction
              ## Which should be dipicted from the DirectionList. each group would walk in a specific direction
              ## If they cannot find a connecting edge in that direction, a random cell is then selected.
              ## Make them to prefer to move through habitat cells rather than fields.
              indexsNGroups  <- (max(PopMatWB[,2])+1):(max(PopMatWB[,2])+dim(AinmFromGToSplit)[1])
              indexNGPigs    <- rep(indexsNGroups,AinmFromGToSplit[,2])
              indexNewGNum   <- sort(which(PopMatWB[,2]%in%AinmFromGToSplit[,1] & PopMatWB[,12] == 1),decreasing = T)
              PopMatWB[indexNewGNum,2] <- indexNGPigs
              AinmFromGToSplit[,1] <- indexsNGroups
              GroupsToSplit <- rbind(GroupsToSplit,AinmFromGToSplit)      
            }
          }
        }
      }
      
      if(dim(GroupsToSplit)[1]>0){
        ## First we check that all groups still exist in the population
        Checktmp         <- which(!GroupsToSplit[,1]%in%PopMatWB[,2])
        if(length(Checktmp)>0) GroupsToSplit[-Checktmp,,drop=FALSE]
        ## Determine the number of pixels they may move per day and the direction
        NumbPixMovesTod  <- round(runif(dim(GroupsToSplit)[1],0,1)) ## Lange et al. (2015)
        MovedPixMale <- as.list(matrix(0,ncol=length(NumbPixMovesTod>0)))
        ToRemovejj <- numeric(0)
        if(any(NumbPixMovesTod>0)){
          for(i in 1:max(NumbPixMovesTod)){
            if(any(PopMatWB[,12]==1) & any(NumbPixMovesTod>=i)){
              IndexGroup   <- which(NumbPixMovesTod>=i)
              if(length(IndexGroup)>0){
                for(jj in IndexGroup){
                  #Find the directions for all neighbours
                  NeighboursIndex <- unlist(WBMat[GroupsToSplit[jj,4],DirectionList[[GroupsToSplit[jj,3]]]])
                  #Keep The ones in the selected area :exclude 0's
                  NeighboursIndex <- NeighboursIndex[NeighboursIndex>0]
                  #Check the habitat cat
                  NbHabCat <- WBMat[NeighboursIndex,5]
                  #Keep accessible and suitable ones only
                  SuitableNb <- which(NbHabCat>0 & NbHabCat<3)
                  TSuitableNb <- NeighboursIndex[SuitableNb]
                  
                  # if(length(TSuitableNb)>0){
                  #if only 1 suitable neighbour => go there
                  if(length(TSuitableNb)==1) Destination <- TSuitableNb
                  #if more thqn one select one neighbour randomly
                  if(length(TSuitableNb)> 1) Destination <- sample(TSuitableNb,1,prob=HabitatProb[NbHabCat[SuitableNb]])
                  # }
                  if(length(TSuitableNb)==0){
                    NeighboursIndex <- unlist(WBMat[GroupsToSplit[jj,4],6:13])
                    NeighboursIndex <- NeighboursIndex[NeighboursIndex>0]
                    NbHabCat <- WBMat[NeighboursIndex,5]
                    SuitableNb <- which(NbHabCat>0 & NbHabCat<3)
                    TSuitableNb <- NeighboursIndex[SuitableNb]
                    if(length(TSuitableNb)>1)  Destination <- sample(TSuitableNb,1,prob=HabitatProb[NbHabCat[SuitableNb]])
                    if(length(TSuitableNb)==1) Destination <- TSuitableNb
                  }
                  
                  ## Here we add the list with the moved pixels
                  MovedPixMale[[jj]] <- c(MovedPixMale[[jj]],Destination)
                  ## Is this only a suitable habitat pixel or a suitable habitat and has a female(s) with no or 1 male
                  
                  TMPIndPix <- (WBMat[Destination,5]==1 & sum(PopMatWB[PopMatWB[,7]==Destination&PopMatWB[,9]!=3,3]==1)>0 & sum(PopMatWB[PopMatWB[,7]%in%Destination,9]!=3)>0 &
                                  sum(PopMatWB[PopMatWB[,7]==Destination,3]==0)<2) #& unique(PopMatWB[PopMatWB[,7]%in%tmp5,2])<2
                  # the first part of the c() represent a good habitat cell and not occupied
                  Freepixel      <- which(c(WBMat[Destination,5]==1 & !(Destination%in%PopMatWB[,7]), TMPIndPix))
                  
                  ## We decide whether they will settle in this pixel or not based on a random process
                  ## Notice here that the animals that die during splitting do not affect the splitting. 
                  ## We checked above that all groups in splitting matrix have animals in the Population matrix
                  if(length(Freepixel)>0){
                    Settle <- rbinom(1,1,prob = SettlingProb[max(Freepixel)])==1
                    if(Settle&!TMPIndPix){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      #newGroupIDM             <- rep(max(PopMatWB[,2])+1,length(indexspPigs))
                      #PopMatWB[indexspPigs,2] <- newGroupIDM
                      PopMatWB[indexspPigs,7] <- Destination
                      PopMatWB[indexspPigs,8] <- Destination
                      PopMatWB[indexspPigs,12]<- 0
                      ToRemovejj              <- c(ToRemovejj,jj)
                    }
                    
                    ## If the group will settle with a pre-existing female group, then their group number will be the same as the females
                    if(Settle&TMPIndPix){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      IndexSelPigs            <- which(PopMatWB[,7]%in%Destination & PopMatWB[,9]!=3 & !PopMatWB[,2]%in%GroupsToSplit[jj,1] & !PopMatWB[,2]%in%GroupsToSplit[GroupsToSplit[,4]%in%PopMatWB[,7],1])
                      newGroupIDM             <- rep(unique(PopMatWB[IndexSelPigs,2]),length(indexspPigs)) #this would allow settling in an infected pixel with carcass 
                      PopMatWB[indexspPigs,2] <- newGroupIDM
                      PopMatWB[indexspPigs,7] <- Destination
                      PopMatWB[indexspPigs,8] <- Destination
                      PopMatWB[indexspPigs,12]<- 0
                      ToRemovejj              <- c(ToRemovejj,jj)
                    }
                    if(!Settle){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      PopMatWB[indexspPigs,8] <- Destination
                      GroupsToSplit[jj,4]     <- Destination
                    }
                  }
                  
                  if(length(Freepixel)==0) {
                    indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                    PopMatWB[indexspPigs,7] <- Destination
                    PopMatWB[indexspPigs,8] <- Destination ## we set the home pixel as the current pixel for splitting groups as they do not have home yet
                    GroupsToSplit[jj,4]     <- Destination
                  }
                  
                }
              }
            }
          }
          names(MovedPixMale) <- GroupsToSplit[,1]
          if(length(ToRemovejj)>0) GroupsToSplit  <- GroupsToSplit[-ToRemovejj,,drop=FALSE]
        }
      }
      
      ### Seed infection
      
      
       # if(TimeSeedInf==gTime){
       #   ## Here I seed randomly in 2 groups and I infect 1 animal (make it directly infectious). We make sure that the groups have at least 3 animals.
       #   NumPGroup  <- tapply(PopMatWB[,2],PopMatWB[,2],length)
       #   temp1      <- NumPGroup[NumPGroup>2]
       #   temp2      <- as.numeric(names(temp1))
       #   temp3      <- temp2[temp2>0]
       #   RandSeedG  <- sample(temp3,NumbGSeedInf)
       #   SampSeedAn <- sapply(RandSeedG, function(x) sample(which(PopMatWB[,2]==x),1))
       #   PopMatWB[SampSeedAn,9]  <- 2
       #   PopMatWB[SampSeedAn,14] <- gTime   # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
       #   PopMatWB[SampSeedAn,15] <- PopMatWB[SampSeedAn,14] + round(rpert(length(SampSeedAn),1,5,7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
       # }
       # 
      ## Here we model between group spread
      ## First we model infection between groups within the same area
      ## Here we model within group infection and infection from connecting edges via carcasses
      ## Infectious groups are those with status ==2 and carcases ==3 (0=sus, 1= exp, 2=inf, 3=dead, 4=immune)
      if(sum(PopMatWB[,9])>0){
        infGroups <- unique(PopMatWB[PopMatWB[,9]%in%2:3,2]) 
        for(i in infGroups){
          NumInfPG  <- sum(PopMatWB[PopMatWB[,2]==i ,9]==2)
          TotDead   <- sum(PopMatWB[PopMatWB[,2]==i,9]==3 & (gTime-PopMatWB[PopMatWB[,2]==i,15])<=CarcDaySur)
          ## Prob from infectious individuals
          probWGDCInf <- 1-(1-WGProbInf)^NumInfPG
          ## Prob from infectious carcasses within the group
          ProbWGCInf  <- 1-(1-CarcProbInf)^TotDead
          ## Prob from infectious carcasses in connecting edges (based on current pixel NOT Home pixel as they maybe splitting)
          tmpConEdgs  <- unlist(WBMat[unique(PopMatWB[PopMatWB[,2]==i,8]),6:13])
          tmpConEdgs  <- tmpConEdgs[tmpConEdgs>0]
          tmpDInCE    <- sum(PopMatWB[,8]%in%tmpConEdgs & PopMatWB[,9]==3)
          ProbBGCInf  <- 1-(1-CarcProbInf)^tmpDInCE
          ## Estimate total probability of infection from contact
          ProbWGInf   <- 1-((1-probWGDCInf)*(1-ProbWGCInf)*(1-ProbBGCInf))
          
          ## Suseptible animals
          tmpSus <- which(PopMatWB[,9]==0 & PopMatWB[,2]==i)
          if(length(tmpSus)>0){
            NewInfPG <- tmpSus[rbinom(length(tmpSus),1,prob=ProbWGInf)==1]
            if(length(NewInfPG)>0){
              PopMatWB[NewInfPG,9]  <- 1 #
              PopMatWB[NewInfPG,14] <- gTime + round(rpert(length(NewInfPG),1,5,9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
              PopMatWB[NewInfPG,15] <- PopMatWB[NewInfPG,14] + round(rpert(length(NewInfPG),1,5,7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
            }
          }
      
          ## Infection of neighboring susceptible groups. 
          ## Here we model the susceptible groups that are around infectious herd i, but we also take into
          ## account that the suceptible herd could have risk from other infectious herds 
          ## Or carcasses around it or a carcass in its own pixel
          GroupsAtRisk <- unique(PopMatWB[(PopMatWB[,8]%in%tmpConEdgs & !PopMatWB[,2]%in%infGroups),2])
          if(length(GroupsAtRisk)>0){
            GroupsAtRisk <- GroupsAtRisk[GroupsAtRisk>0]
            PixelsAtRisk <- unique(PopMatWB[PopMatWB[,2]%in%GroupsAtRisk,8])
            NumInfAtRisk <- sapply(PixelsAtRisk,function(x) {tmp1 <- c(unlist(WBMat[x,6:13]),x) # include the possibility of an infected carcass in the pixel from another group
            tmp2 <- PopMatWB[,8]%in%tmp1
            sum(tmp2&PopMatWB[,9]%in%2:3)}) 
            if(length(NumInfAtRisk)>0){
              ProbInfGAR   <- 1-(1-CarcProbInf)^NumInfAtRisk
              StatusGAR    <- GroupsAtRisk[rbinom(length(GroupsAtRisk),1,ProbInfGAR)==1]
              if(length(StatusGAR)>0){
                NumAnimPGAR <- sapply(StatusGAR,function(x)sum(PopMatWB[,2]==x))
                NumInfAnim  <- round(runif(length(StatusGAR),1,3))
                NumInfAnim[NumInfAnim>NumAnimPGAR]  <- 1
                for(ss in 1:length(StatusGAR)){
                  IndexToSelect1 <- which(PopMatWB[,2]==StatusGAR[ss])
                  if(NumAnimPGAR[ss] > NumInfAnim[ss]) IndexToSelect2 <- sample(IndexToSelect1,NumInfAnim[ss]) else IndexToSelect2 <- IndexToSelect1
                  PopMatWB[IndexToSelect2,9]  <- 1
                  PopMatWB[IndexToSelect2,14] <- gTime + round(rpert(length(IndexToSelect2),1,5,9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
                  PopMatWB[IndexToSelect2,15] <- PopMatWB[IndexToSelect2,14] + round(rpert(length(IndexToSelect2),1,5,7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
                  
                }
              }
            }
          }
        }
            
            ### Groups that can be infected via carcass, while moving.
            
            
            ## First male groups can be infected while moving (if susceptibles) or infect others (if infectious)
            GroupsSplit    <- as.numeric(names(MovedPixMale))
            if(length(GroupsSplit)>0){
              MovedGroups    <- unname(unlist(lapply(MovedPixMale,function(x)length(x)>1)))
              MovedGroupsTod <- GroupsSplit[MovedGroups]
              if(length(MovedGroupsTod)>0){
                SusMovMalG     <- sapply(MovedGroupsTod,function(x) all(PopMatWB[,9]==0&PopMatWB[,2]==x & PopMatWB[,12]==1))
                SusMaleGroups  <- MovedGroupsTod[SusMovMalG]
                if(length(SusMaleGroups)>0){
                  for(i in SusMaleGroups){
                    tmp <- unname(unlist(MovedPixMale[names(MovedPixMale)==i]))
                    tmp <- tmp[tmp>0]
                    ## pixels that have carcasses
                    AnyInfPopInPix <- tmp[sapply(tmp,function(x) any(PopMatWB[,8]%in%x & PopMatWB[,9]%in%3))] 
                    if(length(AnyInfPopInPix)>0){
                      ## number of carcasses per pixel
                      NumCarPP     <- sapply(AnyInfPopInPix,function(x) sum(PopMatWB[,8]%in%x & PopMatWB[,9]==3) )
                      probInfSplit <- 1-(1-CarcProbInf)^sum(NumCarPP)
                      NewInfSG     <- rbinom(1,1,probInfSplit)
                      if(NewInfSG==1){
                        SusAnimal   <- which(PopMatWB[,2]==i & PopMatWB[,12]==1)
                        NumInfAnim  <- round(runif(1,1,3))
                        NumInfAnim[NumInfAnim>length(SusAnimal)] <- SusAnimal
                        if(SusAnimal > NumInfAnim) IndexToSelect3 <- sample(SusAnimal,NumInfAnim) else IndexToSelect3 <- SusAnimal
                        PopMatWB[IndexToSelect3,9]  <- 1
                        PopMatWB[IndexToSelect3,14] <- gTime + round(rpert(length(IndexToSelect3),1,5,9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
                        PopMatWB[IndexToSelect3,15] <- PopMatWB[IndexToSelect3,14] + round(rpert(length(IndexToSelect3),1,5,7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
                      }
                    }
                  }
                }
              }
            }
       
            
            ### Groups split females
            
            
            GroupsSplit    <- as.numeric(names(PixelsMoved))
            if(length(GroupsSplit)>0){
              MovedGroups    <- unname(unlist(lapply(PixelsMoved,function(x)length(x)>1)))
              MovedGroupsTod <- GroupsSplit[MovedGroups]
              if(length(MovedGroupsTod)>0){
                SusMovFMalG    <- sapply(MovedGroupsTod,function(x) all(PopMatWB[,9]==0 &PopMatWB[,2]==x & PopMatWB[,10]==1))
                SusFemaleGroups<- MovedGroupsTod[SusMovFMalG]
                if(length(SusFemaleGroups)>0){
                  for(i in SusFemaleGroups){
                    tmp <- unname(unlist(PixelsMoved[names(PixelsMoved)==i]))
                    tmp <- tmp[tmp>0]
                    ## pixels that have carcasses
                    AnyInfPopInPix <- tmp[sapply(tmp,function(x) any(PopMatWB[,8]%in%x & PopMatWB[,9]%in%3))] 
                    if(length(AnyInfPopInPix)>0){
                      ## number of carcasses per pixel
                      NumCarPP     <- sapply(AnyInfPopInPix,function(x) sum(PopMatWB[,8]%in%x & PopMatWB[,9]==3) )
                      probInfSplit <- 1-(1-CarcProbInf)^sum(NumCarPP)
                      NewInfSG     <- rbinom(1,1,probInfSplit)
                      if(NewInfSG==1){
                        SusAnimal   <- which(PopMatWB[,2]==i & PopMatWB[,10]==1)
                        NumInfAnim  <- round(runif(1,1,3))
                        NumInfAnim[NumInfAnim>length(SusAnimal)] <- SusAnimal
                        if(SusAnimal > NumInfAnim) IndexToSelect3 <- sample(SusAnimal,NumInfAnim) else IndexToSelect3 <- SusAnimal
                        PopMatWB[IndexToSelect3,9]  <- 1
                        PopMatWB[IndexToSelect3,14] <- gTime + round(rpert(length(IndexToSelect3),1,5,9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
                        PopMatWB[IndexToSelect3,15] <- PopMatWB[IndexToSelect3,14] + round(rpert(length(IndexToSelect3),1,5,7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
                      }
                    }
                  } 
                }
              }
            }
                  
            ### Update Animal status
            
            
            # Update animals to Infectious status
            IndexToInf <- PopMatWB[,9] == 1 & PopMatWB[,14] == gTime
            if(sum(IndexToInf)>0) PopMatWB[IndexToInf,9] <- 2
            
            # Update animals to Dead/Immune state
            IndexToDOrI <- which(PopMatWB[,9] == 2 & PopMatWB[,15] == gTime)
            if(length(IndexToDOrI)>0){
              tmpToDie <- rbinom(length(IndexToDOrI),1,ProbDeathInf)
              if(any(tmpToDie==1)) {
                PopMatWB[IndexToDOrI[tmpToDie==1],9]  <- 3
                ## dead animals may die within their own home range 20% or any of the connecting edges 80% (splitting male groups are excluded).
                tmpDieOSG   <- IndexToDOrI[tmpToDie==1 & (!PopMatWB[IndexToDOrI,2]%in%GroupsToSplit[,1])]
                DieOutGroup <- tmpDieOSG[rbinom(length(tmpDieOSG),1,DieInConEdge)]
                if(length(DieOutGroup)>0){
                  EdgesDie <- sapply(DieOutGroup,function(x) {
                    tmp1 <- unlist(WBMat[PopMatWB[x,8],6:13])
                    tmp2 <- WBMat[tmp1,5]
                    tmp3 <- tmp1[tmp2%in%1:2]
                    tmp3 <- tmp3[tmp3>0]
                    if(length(tmp3)>1) sample(tmp3,1) else(PopMatWB[x,8])
                  })
                  ### If an animal will die in a connecting edge, its pixel home and current will change but not the group number as the probablity
                  ### of infection from Carcasses within a home-range and the connecting edges is the same according to Lange et al. (2015,EFSA) .
                  PopMatWB[DieOutGroup,7] <- EdgesDie
                  PopMatWB[DieOutGroup,8] <- EdgesDie
                }
              }
              if(any(tmpToDie==0)) PopMatWB[IndexToDOrI[tmpToDie==0],9]  <- 4
            }
            
            ## Remove dead carcasses over 28 days from death from the matrix
            indexDeadOld <- PopMatWB[,9] == 3 & ((gTime-PopMatWB[,15]) > CarcDaySur)
            
            if(sum(indexDeadOld)>0){
              PopMatWB[indexDeadOld,] <- 0
              PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
            } 
            
            ## Dead animals must be removed from the splitting mechanism
            # Identify group number of animals that died today
            IndexDiedTod  <- sort(unique(PopMatWB[PopMatWB[,9]==3 & PopMatWB[,15] == gTime,2]))
            MIndexDiedTod <- IndexDiedTod%in%GroupsToSplit[,1]
            if(sum(MIndexDiedTod)>0){
              IndexToRemSM <- IndexDiedTod[MIndexDiedTod]
              NumRemTodFSM <- tapply((PopMatWB[PopMatWB[,2]%in%IndexToRemSM,9]==3 & PopMatWB[PopMatWB[,2]%in%IndexToRemSM,15] == gTime),PopMatWB[PopMatWB[,2]%in%IndexToRemSM,2],sum)
              GroupsToSplit[match(IndexToRemSM,GroupsToSplit[,1]),2] <- GroupsToSplit[match(IndexToRemSM,GroupsToSplit[,1]),2] - NumRemTodFSM
              if(any(GroupsToSplit[,2]<1)){
                IndexDelFSM   <- which(GroupsToSplit[,2]<1)
                GroupsToSplit <- GroupsToSplit[-IndexDelFSM,,drop=FALSE]
              }
            }
      }
      
      ## Here er re-initiate the daily matrix for males and females movements to insure that there is no bleed from previous days
      MovedPixMale <- as.list(matrix(0,ncol=1))
      PixelsMoved <- as.list(matrix(0,ncol=1))
      
      ## Make the daily summary
      DayOutInfMat[gTime,Iter] <- length(unique(PopMatWB[PopMatWB[,9]%in%2:4,2]))
      DayOutPopMat[gTime,Iter] <- length(unique(PopMatWB[,2]))
      DayOutAniMat[gTime,Iter] <- sum(PopMatWB[,2]>0)
      
      
      InfGroupsD         <- unique(PopMatWB[PopMatWB[,9]%in%2:4,2])
      InfGroupsD         <- InfGroupsD[!InfGroupsD%in%NewInfGroups[NewInfGroups[,1]==Iter,3]]
      InfAnimals         <- which(PopMatWB[,9]==2)
      InfAnimals         <- InfAnimals[!InfAnimals%in%NewInfAnimals[NewInfAnimals[,1]==Iter,3]]
      
      if(length(InfGroupsD)>0){
        NewInfGroups    <- rbind(NewInfGroups,cbind(Iter,gTime,InfGroupsD))
      }
      if(length(InfAnimals)>0){
        NewInfAnimals   <- rbind(NewInfAnimals,cbind(Iter,gTime,InfAnimals))
      }
      ## a variable to count whether the disease faded out this iteration and started up again. this is counted once/iteration despite in while loop   
      TMPOutToday <- TMPOutYesterday & sum(PopMatWB[,9]%in%1:2)>0
      if(TMPOutToday&OnlyOnce){
        FreqRelapse[Iter] <- 1
        OnlyOnce <- FALSE
      }
      TMPOutYesterday <- gTime>TimeSeedInf & sum(PopMatWB[,9]%in%1:2)==0 & sum(PopMatWB[,9]==3)>0
      
      print(c(Iter,gTime))
      
      if(Iter == MaxIterations){
        
        #Daily Population Count
        PopMatWB <- as.data.frame(PopMatWB)
        Popcount <- PopMatWB %>% 
          group_by(Current_pixel) %>% 
          dplyr::summarise(n())
        Popcount <- t(Popcount)
        #Final population count (number of animals) each day (1825 rows) and in each cell (323 columns)
        Fpopcount[gTime, Popcount[1,-1]] = Popcount[2,-1] 
        #-1 because the first column in Popcount are the count of empty cells (0)
        #The observations on 1st row in Popcount that match each column in Fpopcount will be the number on line 2 of popcount, 
        #which is the total number of animals en that grid
        
        #Daily Infected Count
        Infcount <- PopMatWB %>%
          group_by(Current_pixel, Infect_status) %>%
          filter(Infect_status == 2) %>%
          dplyr::summarise(n())
        Infcount <- t(Infcount)
        if(gTime==TimeSeedInf){
          print(which(Infcount>0))
          
          #Final Infected count (number of infected animals) each day in each cell
          Finfcount[gTime, Infcount[1,]] = Infcount[3,]
          
          #Carcass Daily Count
          CarcassCount <- PopMatWB %>%
            group_by(Current_pixel, Infect_status) %>%
            filter(Infect_status == 3) %>%
            dplyr::summarise(n())
          CarcassCount <- t(CarcassCount)
          Fcarcount[gTime, CarcassCount[1,]] = CarcassCount[3,]
          
          #Immune Count
          Immunecount <- PopMatWB %>%
            group_by(Current_pixel, Infect_status) %>%
            filter(Infect_status == 4) %>%
            dplyr::summarise(n())
          Immunecount <- t(Immunecount)
          Fimmunecount[gTime, Immunecount[1,]] == Immunecount[3,]
          
          NAMEINFC <- paste(runID, 'InfectedCount.txt', sep = '-')
          NAMECARC <- paste(runID, 'CarcassCount.txt', sep = '-')
          NAMEIMC <- paste(runID, 'ImmuneCount.txt', sep = '-')
          NAMEPOPC <- paste(runID, 'PopCount.txt', sep ='-')
          
          write.table(Finfcount, NAMEINFC, sep = ' ', col.names = T, row.names = T)
          write.table(Fcarcount, NAMECARC, sep = ' ', col.names = T, row.names = T)
          write.table(Fimmunecount, NAMEIMC, sep = ' ', col.names = T, row.names = T)
          write.table(Fpopcount, NAMEPOPC, sep = ' ', col.names = T, row.names = T)
          
        }
      }#while(gTime
      
      ## Make the summaries per iteration
      EpDuration[Iter]            <- gTime - TimeSeedInf
      cumDeath[Iter]              <- cumDeathPar + sum(PopMatWB[,2]>0)
      #print(Iter)
      
    }#Close for(Iter in...)
    NAMENA    <- paste(runID,"FNewInfAnimals.txt",sep="-")
    NAMENI    <- paste(runID,"FNewInfGroups.txt",sep="-")
    NAMEIG    <- paste(runID,"FDayOutInfMat.txt",sep="-")
    NAMETG    <- paste(runID,"FDayOutPopMat.txt",sep="-")
    NAMETA    <- paste(runID,"FDayOutAniMat.txt",sep="-")
    
    #(NewInfAnimals, 'NAMENA', append = T, sep=" ", col.names=F, row.names=F)
    write.table(NewInfAnimals, NAMENA, sep=" ", col.names=T, row.names=T)
    write.table(NewInfGroups, NAMENI, sep=" ", col.names=T, row.names=T)
    write.table(DayOutInfMat, NAMEIG, sep=" ", col.names=T, row.names=T)
    write.table(DayOutPopMat, NAMETG, sep=" ", col.names=T, row.names=T)
    write.table(DayOutAniMat, NAMETA, sep=" ", col.names=T, row.names=T)
    
    OutPutList <- cbind(EpDuration, FreqRelapse, cumDeath)
    NAMEOL     <- paste(runID, "FOutPutList.txt", sep="-")
    write.table(OutPutList, NAMEOL, sep=" ", col.names=T, row.names=T)
    
  }
}

WBModel()
