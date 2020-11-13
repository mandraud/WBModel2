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