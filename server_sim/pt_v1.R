###### Simulations ########
rm(list=ls())
source('QUBO_PT_optimization.R')
setwd('..')
#Big simulation
N=100 #Dimension of the QUBO problem
n_steps= 100000#Number of steps in every simulation
n_swaps= 500#After how many steps we try a swap
NumRep <- 1000 #Number of times to repeat the simulation
temps <- c(1,1 *1.09^(1:10)) #temperatures
h_baseline <- rep('h_min',length(temps))
h_m <- rep('h_max',length(temps))
h_2 <- rep('h_sq',length(temps))
if(length(temps)%%2==0){
  h_balancing <- c(rep('h_sq',length(temps)/2),rep('h_max',length(temps)/2))
}else{
  h_balancing <- c(rep('h_sq',(length(temps)-1)/2),rep('h_max',(length(temps)+1)/2))
}

matrix_models <- rbind(h_baseline,h_m,h_2,h_balancing)
set.seed(27101)
Sim_MaxCut_RF_PT(dimension=N, n_steps=n_steps, n_swaps=n_swaps, NumRep=NumRep, betas=1/temps, matrix_h=matrix_models)
