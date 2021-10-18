# Final Scenarios Table ---------------------------------------------------
setwd("X:/Luis/Model")
library(dplyr)
library(magrittr)

ScenarioTable <- NULL
InfCellsTable <- NULL

for (Scenario in c(30:33, 39:42)){
  
  #Study Area
  if(Scenario %in% c(41, 42)) Study_Area = "Department_64" else 
    Study_Area = "France-Belgium_Border"
  
  #Load Data on Study Area
  if(Scenario %in% c(41,42)) load("X:/Luis/Model/Inputs/DEPA64/Dep64.RData") else
     load("X:/Luis/Model/Inputs/EA_9km/Input_EA9.RData")
  
  #File names
  if(Scenario %in% c(41:42)) filename = "Outputs/NewOutputs/WB_Model_Dep64_Scenario_" else 
    filename <- "Outputs/NewOutputs/WB_MovH_ASF_Scenario_"
  
  #Details
 if(Scenario %in% 30) Details = "Pigglets excluded from Hunt, ASF in day 730"
 if(Scenario %in% 31) Details = "Surival probability set to 70%"
 if(Scenario %in% 32) Details = "Surival probability set to 80%"
 if(Scenario %in% 33) Details = "Surival probability set to 90%"
 if(Scenario %in% 39) Details = "Group dispersal caused by hunting" 
 if(Scenario %in% 40) Details = "Carcass Persistence and Seasonality"
 if(Scenario %in% 41) Details = "Department 64 with hunting" 
 if(Scenario %in% 42) Details = "Department 64 without hunting" 
  
  # Mean of total WB Individuals
  animals <- read.delim(paste0(filename, Scenario, "-FDayOutAniMat.txt"), header = T, sep = " ", row.names = NULL)
    
  #Wild Boar Matrix
  if(Scenario %in% c(41, 42)) FileWildBoarMat = "Inputs/DEPA64/Input_dep64.csv" else 
    FileWildBoarMat = "Inputs/EA_9km/Input_EA9.csv"
  
  WBMat <- read.table(FileWildBoarMat, sep = ",", header = T)
  Area <- ((sum(WBMat[,"Habitat"] == 1)*9)/0.01)
  
  animean <- animals %>%
    transmute(WB_Mean = rowMeans(dplyr::select(., V1:V5))) %>%
    mutate(WB_density = round((WB_Mean*100)/Area, digits = 5))
  
  wbmean <- summarise_all(animean, .funs = mean)
  
  
  # Mean Incidence of Infected WB
  Incidence <- read.delim(paste0(filename, Scenario, "-FNewInfAnimals.txt"), header = T, sep = ' ', row.names = NULL) 
  Incidence <- Incidence %>%
    group_by(Iter, gTime) %>% 
    dplyr::summarise(n())
  colnames(Incidence) <- c('Iteration', 'Day', 'Infected_No')
  
  by(Incidence$Infected_No, Incidence$Iteration, cumsum) -> tmp
  Incidence$Infected_Cum <- do.call(what = 'c', tmp)
  MeanIC <- round(quantile(Incidence$Infected_Cum, probs = c(0.5, 0.025, 0.975)), digits = 5)
  
  
  # Mean Infected Carcasses
  Carcasses <- read.delim(paste0(filename, Scenario, "-FNewInfCarcass.txt"), header = T, sep = ' ', row.names = NULL)
  Carcasses <- Carcasses %>%
    group_by(Iter, gTime) %>% 
    dplyr::summarise(n())
  colnames(Incidence) <- c('Iteration', 'Day', 'Infected_Carc')
  
  by(Incidence$Infected_Carc, Incidence$Iteration, cumsum) -> tmp2
  Incidence$Infected_Carc <- do.call(what = 'c', tmp2)
  MeanCC <- round(quantile(Incidence$Infected_Carc, probs = c(0.5, 0.025, 0.975)), digits = 5)
  
  # Mean Epidemic DUration
  OPList <- read.delim(paste0(filename, Scenario,  "-FOutPutList.txt"), header = T, sep = " ", row.names = NULL)
  MeanED <- round(quantile(OPList$EpDuration, probs = c(0.5, 0.025, 0.975)), digits = 5)
  
  # Mean infected Cells
  InfGroups <- read.delim(paste0(filename, Scenario, "-FNewInfGroups.txt"), header = T, sep = ' ', row.names = NULL) 
  InfGroups <- InfGroups %>%
    group_by(Iter, gTime) %>% 
    dplyr::summarise(n())
  colnames(Incidence) <- c('Iteration', 'Day', 'Infected_Groups')
  
  by(Incidence$Infected_Groups, Incidence$Iteration, cumsum) -> tmp3
  Incidence$Infected_Groups <- do.call(what = 'c', tmp3)
  MeanIC <- round(quantile(Incidence$Infected_Groups, probs = c(0.5, 0.025, 0.975)), digits = 5)
  
  #Write table
  ScenarioTable <- rbind(ScenarioTable, c(Scenario       = Scenario,
                                          Study_Area     = Study_Area,
                                          Details        = Details,
                                          WB_Mean        = wbmean$WB_Mean,
                                          WB_Density     = wbmean$WB_density,
                                          Infected_Cum   = MeanIC[1],
                                          Carcass_Cum    = MeanCC[1],
                                          Ep_Duration    = MeanED[1],
                                          Infected_Cells = MeanIC[1]))
  
  #Make the map risk
  
 # InfCells <- read.delim(paste0(filename, Scenario, "-FNewInfGroups.txt"), header = T, sep = ' ', row.names = NULL)
  
#  ICells <- InfCells %>%
#    group_by(InfCells) %>%
#    summarise(n())
#  colnames(ICells) <- c("ID", "Reinfections")
  
 
  
#  habitats@data <- left_join(habitats@data, ICells, by = "ID", copy = F, suffix = c("habitats@data", "ICells"))
  
     
#  tmap_mode('view')
#  map <- tm_shape(habitats) + tm_polygons("Reinfections") + tm_layout(title = Details)
}

ScenarioTable
#map
