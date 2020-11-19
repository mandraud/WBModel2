Ageing<-function(gTime,PopMatWB){
  if(gTime > 365) Year <- ceiling(gTime/365)  ## Count the year 
  PopMatWB[PopMatWB[ ,5] > 0, 5] <- PopMatWB[PopMatWB[,5] > 0, 5] + 1 #sum & day in the age in days column
  
  PopMatWB[PopMatWB[ ,5] >= 365, 4]     <- 2 #convert piglets to subadults if they reach age
  PopMatWB[PopMatWB[ ,5] >= (365*2), 4] <- 3 #convert subadults to piglets
  PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
  
return(PopMatWB)
}


Mortality<-function(PopMatWB,ProbMortPigF,ProbMortPigM,ProbMortSAdF,ProbMortSAdM,ProbMortAdF,ProbMortAdM){
  
  ## We exclude infectious and death animals from disease in order to allow the carcass to persist
  AdultsFToDie    <- which(PopMatWB[ ,4] == 3 & PopMatWB[ ,3] == 1 &! PopMatWB[ ,9] %in% 2:3)
  AdultsFToDie    <- AdultsFToDie[rbinom(length(AdultsFToDie), 1, ProbMortAdF) == 1]
  AdultsMToDie    <- which(PopMatWB[ ,4] == 3 & PopMatWB[ ,3] == 0 &! PopMatWB[ ,9] %in% 2:3)
  AdultsMToDie    <- AdultsMToDie[rbinom(length(AdultsMToDie), 1, ProbMortAdM) == 1]
  
  SubAdultsFToDie <- which(PopMatWB[ ,4] == 2 & PopMatWB[ ,3] == 0 &! PopMatWB[ ,9] %in% 2:3)
  SubAdultsFToDie <- SubAdultsFToDie[rbinom(length(SubAdultsFToDie), 1, ProbMortSAdF) == 1]
  SubAdultsMToDie <- which(PopMatWB[ ,4] == 2 & PopMatWB[ ,3] == 1 &! PopMatWB[ ,9] %in% 2:3)
  SubAdultsMToDie <- SubAdultsMToDie[rbinom(length(SubAdultsMToDie), 1, ProbMortSAdM) == 1]
  
  PigsFToDie      <- which(PopMatWB[ ,4] == 1 & PopMatWB[ ,3] == 0 &! PopMatWB[ ,9] %in% 2:3)
  PigsFToDie      <- PigsFToDie[rbinom(length(PigsFToDie),1,ProbMortPigF)==1]
  PigsMToDie      <- which(PopMatWB[ ,4] == 1 & PopMatWB[ ,3] == 1 &! PopMatWB[ ,9] %in% 2:3)
  PigsMToDie      <- PigsMToDie[rbinom(length(PigsMToDie), 1, ProbMortPigM) == 1]
  TooOldAni       <- which(PopMatWB[ ,5] == (11*365))
  
  # Piglets that are less than 8 weeks and do not have a mother will die.
  # This is based on Petersen (1999) that showed that the majority of piglets would graze on their own when they reach 8 weeks 
  TooYoungNoMoth  <- which(PopMatWB[ ,5] < (8*7) &! (PopMatWB[ ,11] %in% PopMatWB[ ,1]) &! PopMatWB[ ,9] %in% 2:3)
  ToDieNormal <- c(AdultsFToDie, AdultsMToDie, SubAdultsFToDie, SubAdultsMToDie, PigsFToDie, PigsMToDie, TooYoungNoMoth, TooOldAni)
  
  
  return(ToDieNormal)
}

FemaleSplitting<-function(gTime,PopMatWB,WBMat,Year,Distance,DistThreshold){
  ## Occurs only once in week 28 Kramer-Schadt et al. (2009)
  if(gTime %in% ((28*7 + ((Year - 1)*365)) + 0:6)){
    
    ## Calculate the day in the week
    #DayInWeekSSplit <- (1-(ceiling((gTime-(365*(Year-1)))/7) - ((gTime-(365*(Year-1)))/7)))*7
    PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
    AnimalsWInG     <- unlist(sapply(unique(PopMatWB[ ,2]), function(x) 1:sum(PopMatWB[ ,2] == x)))
    BreedCapCells   <- WBMat[match(PopMatWB[ ,7],WBMat$ID) ,15]
    ## Identify the animals within the groups that will split
    ## Females sub-adults > breading capacity, did not split before and are not sick
    GroupSplitNum   <- cbind(sort(unique(PopMatWB[ ,2])), (tapply((AnimalsWInG > BreedCapCells &
                                                                     PopMatWB[ ,4] == 2 &
                                                                     PopMatWB[ ,3] == 1 &
                                                                     PopMatWB[,10] == 0 &
                                                                     PopMatWB[,9] < 2),  PopMatWB[ ,2], sum)))
    GroupSplitNum   <- GroupSplitNum[GroupSplitNum[ ,2] >= 2, , drop = FALSE]
    
    ## Identify groups where there are only males
    ## To allow females to join groups where there are only males.
    PixelMalesOnly  <- cbind(sort(unique(PopMatWB[ ,2])), (tapply((PopMatWB[ ,3] == 0), PopMatWB[ ,2], all)))
    PixelMalesOnly  <- PixelMalesOnly[PixelMalesOnly[ ,2] == 1, , drop = FALSE]
    OriginPixel     <- PopMatWB[match(GroupSplitNum[ ,1], PopMatWB[ ,2]), 7]
    TargetPixel  <- sapply(OriginPixel, function(x) {
      tmp1 <- which(Distance[x,] <= DistThreshold)
      tmp2 <- tmp1[(WBMat[tmp1, 5] == 1 & !(tmp1 %in% PopMatWB[!PopMatWB[ ,2] %in% PixelMalesOnly[ ,1], 7]))]
      if(length(tmp2) > 1)  tmp3 <- sample(tmp2, 1)
      if(length(tmp2) == 1) tmp3 <- tmp2
      if(length(tmp2) == 0) tmp3 <- 0
      tmp3
    })
    GroupSplitNumSD <- data.frame(group=GroupSplitNum[,1],nMoved=GroupSplitNum[,2], OriginPixel=OriginPixel, TargetPixel=TargetPixel)
    GroupSplitNumSD <- GroupSplitNumSD[GroupSplitNumSD$TargetPixel > 0, , drop=FALSE]
    
    ## Groups to split, may split during any day of the week. No need for a probability here the group splits only once a year.
    #if(dim(GroupSplitNumSD)[1]>0) GroupSplitNumSD <- GroupSplitNumSD[rbinom(dim(GroupSplitNumSD)[1],1,prob=ProbSplitSDS)==1,,drop=FALSE]
    PopMatWB[PopMatWB[ ,2] %in% GroupSplitNumSD$group, 10] <- 1  #Make all group split status 1??
    
    ## Make short term split to happen
    if(nrow(GroupSplitNumSD) > 0){
      
      ## This part of the code to allow movement
      
      ## A list to keep track for the pixels where the pigs have been. 
      ## Make sure that this information is exported on daily basis, because the 
      ## list will be re-initiate every day splitting may happen.
      PixelsMoved <- vector(mode = "list", length =nrow(GroupSplitNumSD))
      
      
      CurrentPos <- GroupSplitNumSD$OriginPixel
      for(i in 1:nrow(GroupSplitNumSD)){
        Trail <- 0
        prevEdge <- CurrentPos[i]
        while(CurrentPos[i] != GroupSplitNumSD$TargetPixel[i] & Trail < 10){
          Trail      <- Trail + 1
          Edges      <- unlist(WBMat[CurrentPos[i], 6:13])
          Edges      <- Edges[Edges > 0 & Edges != prevEdge]
          if(length(Edges) > 0){
            prevEdge    <- CurrentPos[i]
            NewPosition <- Edges[which.min(Distance[GroupSplitNumSD$TargetPixel[i],Edges])]
            if(length(NewPosition) > 1) NewPosition <- sample(NewPosition, 1)
            CurrentPos[i] <- NewPosition 
            # Here we keep track of the pixels where the pigs moved to.
            PixelsMoved[[i]]<- c(PixelsMoved[[i]], NewPosition)
            #print(c(i,OriginPixel[i],prevEdge,Trail,CurrentPos[i],TargetPixel[i]))
          }
        }
        if (CurrentPos[i] == GroupSplitNumSD$TargetPixel[i]){
          unique(BreedCapCells[PopMatWB[,2] == GroupSplitNumSD$group[i]])->RemainBreeders
          movingFemales <- which(PopMatWB[,2] == GroupSplitNumSD$group[i] & PopMatWB[,9]<2)
          movingFemales <-movingFemales[RemainBreeders+(1:GroupSplitNumSD$nMoved[i])]
          if (CurrentPos[i] %in% PopMatWB[,7]){ 
            newGroup <- unique(PopMatWB[PopMatWB[,7]==CurrentPos[i],2])
            #There might have a bug somewhere in males movements, 
            #or maybe it's normal as they keep moving, but
            #multiple groups can have the same home pixel. 
            #for now I rearrange to keep only one group ID per habitat cell 
            #when females arrive in such a group. 
            if (length(newGroup)>1) {
              newGroup<-min(newGroup)
              PopMatWB[PopMatWB[,7]==CurrentPos[i],2]<-newGroup
            }
          }else{
            newGroup <- max(PopMatWB[,2])+1
          }
          PopMatWB[movingFemales,2] <- newGroup
          PopMatWB[movingFemales,7] <- CurrentPos[i]
        }
      }
    }
    
  }
  return(PopMatWB)
}


Reproduction<-function(gTime, PopMatWB, NumOfSprProbList){
  if(gTime %in% PopMatWB[ ,6]){
    DelIndex        <- which(PopMatWB[ ,6] %in% gTime & PopMatWB[ ,9] != 3)
    NumOfSpring     <- sample(NumOfSprProbList$Number, length(DelIndex), rep = T, prob = NumOfSprProbList$Prob)
    DelIndex        <- DelIndex[NumOfSpring > 0]
    NumOfSpring     <- NumOfSpring[NumOfSpring > 0]
    if(length(NumOfSpring) > 0){
      newBorns<-cbind((max(PopMatWB[ ,1]) + 1):(max(PopMatWB[ ,1]) + sum(NumOfSpring)),
                      rep(PopMatWB[DelIndex,2],NumOfSpring),
                      rbinom(sum(NumOfSpring),1,0.5),
                      1,
                      1,
                      0,
                      rep(PopMatWB[DelIndex,7],NumOfSpring),
                      rep(PopMatWB[DelIndex,7],NumOfSpring),
                      0,
                      0,
                      rep(PopMatWB[DelIndex,1],NumOfSpring),
                      0,
                      0,
                      0,
                      0,
                      0)
      colnames(newBorns)=colnames(PopMatWB)
      PopMatWB <- rbind(PopMatWB,newBorns)         
      PopMatWB <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
      
    }
  }
  return(PopMatWB)
}