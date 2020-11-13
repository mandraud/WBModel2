defineWbMat<-function(file){
  WBMat <- read.table(file, sep = ",", header = T)
  WBMat <- cbind(WBMat, 0)
  WBMat <- WBMat[ ,-1]
  colnames(WBMat) <- c("ID", "Lon", "Lat", "Forest_cover", "Habitat", "Nb1", "Nb2",
                       "Nb3", "Nb4", "Nb5", "Nb6", "Nb7", "Nb8", "Colors", "Breeding")
  WBMat[WBMat[ ,5] == 1, 15] <- ceiling(rpert(sum(WBMat[ ,5] == 1), MinBCap, ModeBCap, MaxBCap))
  return(WBMat)  
}


InitPop<-function(InitialOccCells, 
                  Mat=WBMat,
                  Fem=FemPropGroup,
                  Mal=MalPropGroup,
                  SubA=SubAPropGroup,
                  Piglet=PigletPropGroup,AgeProb=AgeProbab){
  TMPHP <- Mat[ ,5] == 1
  homePixelsAll   <- sample(Mat[TMPHP, 1], InitialOccCells)
  PopMatWB=NULL
  size=0
  for(i in 1:InitialOccCells){
    ## The initial distribution of the groups was based on Merta et al. (2015). They were calculated relative to the average
    ## Female percentage. for instance males were average % males divided by average % females. 
    females       <- round(runif(MinBCap, ModeBCap, MaxBCap))
    males         <- round(females*Mal/Fem)
    subAdults     <- round(females*SubA/Fem)
    Piglets       <- round(females*Piglet/Fem)
    GroupID       <- rep(i, sum(females, males, subAdults, Piglets))
    Sex           <- c(rep(1, females), rep(0, males), rbinom(subAdults, 1, prob = 0.5), rbinom(Piglets, 1, prob = 0.5))
    AgeCat        <- c(rep(3, females),rep(3, males), rep(2, subAdults), rep(1, Piglets))
    Dam           <- c(rep(0, sum(females, males)), sample(1:females, subAdults, rep = T), sample(1:females, Piglets, rep = T))
    # Females deliver between Jan and June, so the age of the piglets will be between 183 and 365 days.
    # The same thing for sub-adults but with 365 days extra.
    tmpAgeSubA    <- sample((183:365) + 365, females, rep = T)
    tmpAgePig     <- sample((183:365), females, rep = T)
    Age           <- c(sample(2:10, size = sum(females, males), T, prob = AgeProb)*365, tmpAgeSubA[Dam[(females + males + 1):(females + males + subAdults)]],
                       tmpAgePig[Dam[(females + males + subAdults + 1):(females + males + subAdults + Piglets)]])
    # Adjust the Dam number to fit the actual ID of the DAMs 
    Dam[(females + males + 1):(females + males + subAdults)] <- Dam[(females + males + 1):(females + males + subAdults)] + size# -1 because we have extra row at start
    Dam[(females + males + subAdults + 1):(females + males + subAdults + Piglets)] <- Dam[(females + males + subAdults + 1):(females + males + subAdults + Piglets)] + size
    Breed         <- rep(0, sum(females, males, subAdults, Piglets))
    HomePixel     <- rep(homePixelsAll[i], sum(females, males, subAdults, Piglets))
    CurrPixel     <- HomePixel
    infectStatus  <- rep(0, sum(females, males, subAdults, Piglets))
    SplitStatus   <- rep(0, sum(females, males, subAdults, Piglets))
    SplitMale     <- rep(0, sum(females, males, subAdults, Piglets))
    TimeDeath     <- rep(0, sum(females, males, subAdults, Piglets))
    TimeToInfect  <- rep(0, sum(females, males, subAdults, Piglets))
    TimeToDeath   <- rep(0, sum(females, males, subAdults, Piglets))
    GroupMove     <- rep(0, sum(females, males, subAdults, Piglets))
    IDs           <- (size+ 1):(size + sum(females, males, subAdults, Piglets))
    
    InitMatWBPop  <- cbind(IDs, GroupID, Sex, AgeCat, Age, Breed, HomePixel, CurrPixel, infectStatus, SplitStatus, Dam, SplitMale, TimeDeath, TimeToInfect, TimeToDeath, GroupMove)
    PopMatWB   <- rbind(PopMatWB,InitMatWBPop)
    # Mat[unique(HomePixel), 15] <- females
    size=nrow(PopMatWB)
    
  }
  
  colnames(PopMatWB) <- c("IDs", "Group_ID", "Sex", "Age_Cat", "Age_days", "Breed",
                          "Home_pixel", "Current_pixel", "Infect_status", "Split_status", 
                          "Dam", "Split_male", "Time_death","Time_to_infect", "Time_to_death", "Group_move")
  return(PopMatWB)
}
