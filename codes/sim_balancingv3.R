library(dplyr)
if(grepl('codes',getwd())){setwd('..')}
source('functions/IIT_functions.R')
semilla <- 1546
set.seed(semilla)
N <- 10000 #Number of simulations
steps <- 2*10^5 #Max Number of steps each chain will move 
S <- 500 #Length of state space
mode <- 200 #Max for pi
min_pit <- 1 #Min for pi
reg <- 50 # a mid level between max and min for pi
pi <- matrix(NA,nrow=4,ncol = S)

#4 alternatives to pi
# 1 is basically unimodal with a low probability pit surrounding the mode
pi[1,] <- c(rep(reg,90),rep(min_pit,200),mode,rep(min_pit,110),rep(reg,99))
plot(pi[1,], xlab='state')
# 2 is multimodal with low probabiliti areas among modes
pi[2,] <- c(rep(reg,99),rep(min_pit,100),mode*2,rep(reg,99),mode,rep(min_pit,100),rep(reg,99),mode)
plot(pi[2,], xlab='state')
# 3 similar to 2 but with 1 global optima
pi[3,] <- c(rep(reg,99),mode*2,rep(min_pit,99),mode,rep(reg,99),mode*10,rep(min_pit,99),mode*5,rep(reg,99),mode)
plot(pi[3,], xlab='state')
# 4 a random graph with 1 global mode, 6 local modes
pi[4,] <- sample(c(min_pit,reg,mode),prob = c(.6,.35,.05), replace = T,size=S)
pi[4,sample(1:S,1)] <- mode*10
plot(pi[4,], xlab='state')
plot(pi[4,-which(pi[4,]==max(pi[4,]))], xlab='state',ylab='pi w/o global mode')

write.csv(pi,paste('results/pi',S,semilla,'.csv',sep='_'),row.names = F)

