##########################################################################################################
####################                        Outputs Analysis and Statistics             ##################
                         
Sys.setenv(LANGUAGE="en")
library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(tmap)
library(tmaptools)
library(dplyr)
library(ggplot2)

########            Iterations of NOT INFECTED Population

Pop <- read.delim(file.choose(), header = T, sep = ' ') #Choose the FDayOutPopMat.txt
plot(x = c(1,nrow(Pop)), y = c(0,1300), xlab = 'Time', ylab = 'Groups', type ='n')
poplines <- for (i in 1:nrow(Pop)) {
  lines(1:nrow(Pop), Pop[ ,i], col = i)
}

wildboars <- read.delim(file.choose(), header = T, sep = ' ') #Choose the FDayOutAniMat.txt
plot(x = c(1, nrow(wildboars)), y = c(0, 24000), xlab = 'Time', ylab = 'Wild boars', type = 'n')
poplines <- for (i in 1:nrow(wildboars)) {
  lines(1:nrow(wildboars), wildboars[ ,i], col = i)
}

# Last iteration statistics NO INFECTION

Popcount <- read.delim(file.choose(), header = T, sep = ' ') #Choose the PopCount.txt
TIME <- 1:nrow(Popcount)
which(apply(Popcount, 2, sum) == 0) -> Indexempty
Popcount2 <- Popcount[,-Indexempty] #To remove empty cells
mean <- apply(Popcount2, 1, mean)
plot(TIME, mean, typ = 'l')

e <- ggplot(Popcount2, aes(TIME, mean))
e + geom_line()
e + stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=mean), alpha=0.3)


########            Iterations of Infected Population

#Total Population with ASF per Iteration Graph
PopI <- read.delim(file.choose(), header = T, sep = ' ') #Choose the FDayOutPopMat.txt
plot(x = c(1,nrow(PopI)), y = c(0,1500), xlab = 'Time', ylab = 'Groups', type ='n')
poplines <- for (i in 1:nrow(PopI)) {
  lines(1:nrow(PopI), PopI[ ,i], col = i)
}

#Incidence of ASF infections per iteration
Incidence <- read.delim(file.choose(), header = T, sep = ' ', row.names = NULL) #Choose FNewInfAnimals.txt
Incidence <- Incidence %>%
  group_by(Iter, gTime) %>% 
  dplyr::summarise(n())
colnames(Incidence) <- c('Iteration', 'Day', 'Infected_No')

plot(x = Incidence$Day, y = Incidence$Infected_No, type = 'n')
survival <- for (i in Incidence$Iteration) {
  lines(Incidence$Day[Incidence$Iteration==i], Incidence$Infected_No[Incidence$Iteration==i], col = rainbow(50)[i], lwd = 0.5)
}

s <- ggplot(Incidence, aes(Day, Infected_No))
s + geom_smooth()

#Cumulative ASF infections per day, per iteration
by(Incidence$Infected_No, Incidence$Iteration, cumsum) -> tmp
Incidence$Infected_Cum <- do.call(what = 'c', tmp)

plot(x = Incidence$Day, y = Incidence$Infected_Cum, type = 'n')
cumulative <- for (i in Incidence$Iteration) {
  lines(Incidence$Day[Incidence$Iteration==i], Incidence$Infected_Cum[Incidence$Iteration==i], col = rainbow(50)[i], lwd = 0.5)
}

c <- ggplot(Incidence, aes(Day, Infected_Cum))
c + geom_smooth()
c + geom_quantile()

#Infected Animals Matrix
Infected <- read.delim(file.choose(), header = T, sep = ' ') #ChooseFdayOutInfMat
plot(x = c(1,nrow(Infected)), y = c(0,700), xlab = 'Time', ylab = 'Infected Animals', type ='n')
inflines <- for (i in 1:nrow(Infected)) {
  lines(1:nrow(Infected), Infected[ ,i], col = i)
}

I <- ggplot(Infected, aes(x = 1:nrow(Infected), y = 1:nrow(Infected)))
I + geom_qq_line()

#To compare more than one lines in one graphic:
outIncidence <- NULL

for (Scenario in c(41, 42)) {
  Incidence <- read.delim(paste0("WB_Model_Dep64_Scenario_", Scenario, "-FNewInfAnimals.txt"), sep = " ", row.names = NULL)
  Incidence <- Incidence %>%
    group_by(Iter, gTime) %>% 
    dplyr::summarise(n())
  colnames(Incidence) <- c('Iteration', 'Day', 'Infected_No')
  
  by(Incidence$Infected_No, Incidence$Iteration, cumsum) -> tmp
  Incidence$Infected_Cum <- do.call(what = 'c', tmp)
  NumbInf <- Incidence$Infected_Cum
  Day     <- Incidence$Day
  
  load(file = paste0("WB_Model_Dep64_Scenario_", Scenario, "-FOutPutList.Rda"))
  Param      <- OutPutParamsList[[3]]
  
  outIncidence <- rbind(outIncidence, data.frame(Scenario = Scenario,
                                                 NI       = NumbInf,
                                                 Day      = Day,
                                                 Hunt     = Param[,"Hunt"], 
                                                 MovHunt  = Param[,"MovHunt"], 
                                                 Cells    = Param[,"Cells"], 
                                                 Init     =  paste0(Param[,"InitGroup"], "Cells")))
}

head(outIncidence)

c <- ggplot(outIncidence, aes(x = Day, y = NI, group = factor(Scenario), color = factor(Scenario)))
c + geom_smooth()

#To compare infected carcasses

outCarcass <- NULL

for (Scenario in c(41, 42)) {
  InfCarcass <- read.delim(paste0("WB_Model_Dep64_Scenario_", Scenario, "-FNewInfCarcass.txt"), sep = " ", row.names = NULL)
  InfCarcass <- InfCarcass %>%
    group_by(Iter, gTime) %>% 
    dplyr::summarise(n())
  colnames(InfCarcass) <- c('Iteration', 'Day', 'Infected_Carcass')
  
  by(InfCarcass$Infected_Carcass, InfCarcass$Iteration, cumsum) -> tmp
  InfCarcass$Infected_Carcass <- do.call(what = 'c', tmp)
  NumbCarc <- InfCarcass$Infected_Carcass
  Day     <- InfCarcass$Day
  
  load(file = paste0("WB_Model_Dep64_Scenario_", Scenario, "-FOutPutList.Rda"))
  Param      <- OutPutParamsList[[3]]
  
  outCarcass <- rbind(outCarcass, data.frame(Scenario = Scenario,
                                             NIC      = NumbCarc,
                                             Day      = Day,
                                             Hunt     = Param[,"Hunt"], 
                                             MovHunt  = Param[,"MovHunt"], 
                                             Cells    = Param[,"Cells"], 
                                             Init     = paste0(Param[,"InitGroup"], "Cells")))
}

head(outCarcass)

c <- ggplot(outCarcass, aes(x = Day, y = NIC, group = factor(Scenario), color = factor(Scenario)))
c + geom_smooth()

#Last iteration ASF statistics
Popcount <- read.delim(file.choose(), header = T, sep = ' ') #Choose PopCount
TIME <- 1:nrow(Popcount)
which(apply(Popcount, 2, sum)==0) -> Indexempty
Popcount2 <- Popcount[,-Indexempty]
mean <- apply(Popcount2, 1, mean)
plot(TIME, mean, typ = 'l')

e <- ggplot(Popcount2, aes(TIME, mean))
e + geom_line(col = 'blue')
e + stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=mean), alpha=0.3)

########            Density Estimates

animals <- read.delim(file.choose(), header = T, sep = ' ') #Choose DayOutAniMat
FileWildBoarMat = "Inputs/dep64.csv"
WBMat <- read.table(FileWildBoarMat, sep=",", header = T)
Area <- ((sum(WBMat[,"Habitat"] == 1)*9)/0.01) #cells are 9 km sq / 1 Ha is 0.01 km sq. Result is in Ha

#Mean of all iteration individuals per day and density of animals for 100 Ha. 
animean <- animals %>%
  transmute(WB_Mean = rowMeans(dplyr::select(., V1:V10))) %>%
  mutate(WB_density = ((WB_Mean*100)/Area))

wbmean <- summarise_all(animean, .funs = mean)

#1 animal/100 Ha is considered low density
#10 animals/100 Ha is considered high density
#animean <- animean %>%
#  mutate(WB_density = ((WB_Mean*100)/Area))

#Density graph
library(prismatic)
TIME <- 1:nrow(animean)
#d <- ggplot(animean, aes(TIME, WB_Mean))
#d + geom_col()

dh <- ggplot(animean, aes(TIME, WB_density))
dh + geom_col(aes(color = WB_density,
                  fill = after_scale(clr_mix(color, mix_in = 'green'))))

habcat <- WBMat %>% group_by(Habitat) %>% summarise(n())

