####
# Code to optimize a QUBO problem with dimension N using parallel tempering
# We use uniform proposal distributions
# We can use different balancing functions
# Mostly using prvious code from Sigeng Chen
####


#rm(list=ls())
library(tidyverse)
library(igraph) #To create random adjacency matrix (for MaxCut)
#setwd(here())
###### Functions ########
### Create a random QUBO matrix
RandomQuboMatrix <- function(std)     
{
  L = matrix(rnorm(N * N, mean = 0, sd = std), nrow = N)
  L[lower.tri(L)] <- 0
  return(L)
}

RandomMaxCutMatrix <- function(p=0.5)     
{
  #https://r.igraph.org/index.html
  temp <- sample_gnp(N,p)
  adjmat <- temp[]
  diag(adjmat) <- as.matrix(-1* adjmat%*%rep(1,N)) #to obtain the values for the diagonals
  return(-1*adjmat)
}

### Initialize the first state
RandomStart <- function()
{
  X = rbinom(N, 1, 0.5)
  return(X)
}

### Given a vector and a matrix, return the energy
Energy <- function(X, GivenMatrix, beta)
{
  return(as.numeric(beta * (t(X) %*% GivenMatrix %*% X)))
}

### Return the exponential of the energy
#beta is the inverse temperature
# beta = 1/T
# ExpEnergy <- function(X, GivenMatrix, beta)
# {
#   return(exp(Energy(X, GivenMatrix, beta)))
# }

alpha_temp <- function(GivenMatrix, XNow, h, beta=1){
  EnergyNow =  Energy(XNow, GivenMatrix,beta) #Energy at current state
  logProb = rep(0, N)
  #Check energy of neighbors
  #Neighbor set is composed by all states that can be reached by flipping 1 coordinate
  for(r in 1:N)
  {
    XNew = XNow 
    XNew[r] = 1 - XNew[r]
    EnergyNew = Energy(XNew, GivenMatrix, beta)
    logProb[r] = h(EnergyNew, EnergyNow) #un-normalized probability of accepting the flipping of coordinate r
  }
  return(logProb)
}

### Given a state, provides the weight of the current state and the next accepted state
RF_Update_temp <- function(GivenMatrix, XNow, h, beta=1)
{
  logProb <- alpha_temp(GivenMatrix, XNow, h, beta) #Calculate the logProb vector (contains log probabilities of transitioning from x to all neighbors)
  U = runif(N)
  logD = log(- log(U)) - logProb #Selecting a state proportional to the Prob vector
  #We are simply using the log of the original formula for simplicity.
  index = which.min(logD) #Select index to flip coordinate
  XNew = XNow
  XNew[index] = 1 - XNew[index] #Select the new state
  #Energy of the selected state
  EnergyNew = Energy(XNew, GivenMatrix, beta)
  #The weight is 1/Z_h 
  #Z_h is the sum, for all neighbors, of q(y|x) h(pi(y)/pi(x)) since we're using a uniform q(y|x)=1/N we can use the mean 
  #weight = 1/mean(exp(logProb))
  return(list(XNew, EnergyNew))
}

ReplicaSwapP <- function(GivenMatrix, X1,X2,beta1,beta2,h){
  #We use mean since we're using a uniform proposal distribution over all neighbors
  #We use exp since we're obtaining log probabilities
  a11 <- mean(exp(alpha_temp(GivenMatrix, X1, h, beta1)))
  a12 <- mean(exp(alpha_temp(GivenMatrix, X1, h, beta2)))
  
  a21 <- mean(exp(alpha_temp(GivenMatrix, X2, h, beta1)))
  a22 <- mean(exp(alpha_temp(GivenMatrix, X2, h, beta2)))
  
  EnergyCurT =  Energy(X1, GivenMatrix,beta1) + Energy(X2, GivenMatrix,beta2)
  EnergyNewT = Energy(X1, GivenMatrix,beta2) + Energy(X2, GivenMatrix,beta1)
  accept_p <- exp(EnergyNewT - EnergyCurT)*a12*a21/(a11*a22)
  
  return(accept_p) #Returns probability of accepting the swap
  #The probability is corrected for the rejection-free chain
}

RT_analysis <- function(M){
  #Take a matrix containing the results of swap proposals and acceptance
  #The matrix contains 2 columns, 
  #1st column indicate the state from which the prosal was evaluated (t to t+1)
  #2nd column indicates if the proposal was accepted
  M <- M[-1,] #Delete the first row that has NAs
  total_swaps <- nrow(M) #Calculate the total number of tried swaps
## Compute acceptance rate for the swaps
  M <- as.data.frame(M) #Convert to data frame
  colnames(M) <- c('state','result')
  acceptance_rate <- M |> group_by(state) |> summarise(acc_rate = mean(result))
  n_temps <- max(acceptance_rate$state)+1 #Calculate the total number of temperatures
## Compute the path for the 1st replica
  path <- c(1)

  c_s <- c(0,path[1]) #Current replica is 1, previous replica is 0
  #We will move this according to what the matrix M is saying
  M <- M |> filter(result==1) |> pull(state) #Consider only the accepted swaps
  check <- rep(0,length(M))
  for(i in 1:length(M)){
    if(M[i]==c_s[2]){ #If a swap was accepted for replica T 
      #so T and T+1 swapped
      c_s[1] <- c_s[2] #Update previous state to T
      c_s[2] <- M[i] + 1 #Current state T + 1
      path <- c(path,c_s[2]) #Add the new state
      check[i] <- 1
    }else{
      if(M[i]==c_s[2]-1){#If a swap was accepted for the previous state
        #so T-1 and T swapped
        c_s[1] <- c_s[2] #Update previous state to T
        c_s[2] <- M[i] #Current state T - 1
        path <- c(path,c_s[2]) #Add the new state
        check[i] <- 1
      }
    }
  }
  
## Compute round trip rate using the path
  rt_path <- path[path %in% c(1,n_temps)]
  round_trips <- 0
  place <- 'start'
  for(i in 1:length(rt_path)){
    if(rt_path[i]==n_temps && place=='start'){ place <- 'end'} #the path reached the last temperature
    if(rt_path[i]==1 && place=='end'){place <- 'start'; round_trips <- round_trips + 1} #The path reached the first temperature after a roundtrip
  }
  
  return(list(acceptance_rate,path,round_trips,total_swaps))
}

## h is a vector, we define 1 balancing function for each replica (each temperature)
RejectionFreePT <- function(steps,temps,GivenMatrix,h, swap_n){
  #print('Starting RF-PT')
  # Create an array to store the states visited
  n_rep <- length(temps) #Number of replicas
  current_state <- matrix(NA,nrow=N,ncol=n_rep) #each column is a replica
  max_state <- current_state #array to store the state corresponding to the maximum
  max_energy <- rep(-Inf,n_rep) #vector to store the max value found
  iter_found <- 0 #To record the iteration in which the maximum value was found
  swap_log <- matrix(c(NA,NA),nrow=1,ncol=2)
  count_swaps = 0
  #Randomly and independently initialize all replicas
  for (i in 1:n_rep){
    current_state[,i] <- RandomStart()
  }
  for(i in 1:steps){ #Number of steps that each replica will take
                              #########################################################
     if(i %% floor(steps/4) ==0){
       print(paste0('iteration # ',i," = ",round(100*i/steps,0),"%")) 
     }

### First step: Update the replicas
    for(k in 1:n_rep){
      #We provide a vector with the name of the balancing function to use and then
      #We use eval(parse(text)) to pass the actual function (not the string) to the RF_Update_temp function
      Ans <- RF_Update_temp(GivenMatrix, current_state[,k], h=eval(parse(text=h[k])), beta=temps[k])
      current_state[,k] = Ans[[1]]
      Value = Ans[[2]]
      if(max_energy[k] < Value) #New energy is bigger
      {
        max_energy[k] <-  Value
        max_state[,k] <-  current_state[,k]
        iter_found <- i 
      }
    }
### Second step: Try a swap of replicas
    if(i %% swap_n ==0 && i<steps){ #Every swap_n we try a swap of replicas, except in the very last iteration
      count_swaps <- count_swaps+1 #So next time we swap the replicas of different parity
                                  #########################################################      
      print(paste0('Trying replica swap #',count_swaps,' after ',i,' iterations'))
      #print(dim(swap_log))
      #identify the replicas to swap (even or odd)
      swap_t <- 1:(n_rep-1)
      swap_t <- swap_t[swap_t%%2 == (count_swaps%%2)]
      #Try to swap all replicas (that matches the parity of i) with the replica above them
      for(t in swap_t){
        accept_p <- ReplicaSwapP(GivenMatrix,current_state[,t],current_state[,t+1],temps[t],temps[t+1],h=eval(parse(text=h[k])))
        if(runif(1)<accept_p){#accept the swap of replicas
          temp <- current_state[,t]
          current_state[,t] <- current_state[,t+1]
          current_state[,t+1] <- temp
          swap_log <- rbind(swap_log,c(t,1)) #Record an accepted swap from t to t+1
        }else{
          swap_log <- rbind(swap_log,c(t,0)) #Record a rejected swap from t to t+1
        }
      }
    }
  }
  
  return(list(max_state,max_energy,iter_found, swap_log))
}

### some balancing functions
### Balancing functions to be used in the log-probabilities 
h_sq <- function(new,current){
  return((new-current)/2)}
h_min <- function(new,current){
  return(min(0,new-current))}
h_max <- function(new,current){
  return(max(0,new-current))}


Sim_MaxCut_RF_PT <- function(dimension, n_steps, n_swaps, NumRep, betas, matrix_h){
  #Matrix_h contains as many columns as betas
  #as many rows as the number of balancing functions combinations we want to try
  #Rows of the matrix_h should have the name of the model
  N <- dimension
  num_models <- nrow(matrix_h)
  num_temps <- length(betas)
  seed_s <- floor(runif(1)*1000) #Seed so that we can start each model in the same spot
  if(length(betas)!= ncol(matrix_h)){
    print('error in the dimension of betas and matrix_h')
    return(NA)
  }
  #To store when the algorithm find the optimum
  iter_find <- matrix(NA,nrow=NumRep,ncol=num_models)
  colnames(iter_find) <- rownames(matrix_h)
  #To store the maximum energy found
  max_energy <- matrix(NA,nrow=NumRep,ncol=num_models)
  colnames(max_energy) <- rownames(matrix_h)
  #To store the replica swap rates
  swap_rates <- array(-1,dim=c(NumRep,num_models,num_temps-1))
  #To store the number of roundtrips
  round_trips <- matrix(NA,nrow=NumRep,ncol=num_models)
  colnames(round_trips) <- rownames(matrix_h)
  
  #Start the simulation
  cat('Starting RF PT: \n Dim=',N,'\n# temps =',num_temps,'\n# models =',num_models,'\n# reps =',NumRep,'\n')
  for(i in 1:NumRep){
    set.seed(i+seed_s)
    MaxCut <- RandomMaxCutMatrix()
    for(m in 1:num_models){ #Each iteration repeat the simulation with same starting point for each model
      print(paste0('repetition #',i,', model ',m))
      h_local <- matrix_h[m,] #Choose model m
      set.seed(i+seed_s) #Use the same seed so both algorithms start at the same spot
      simPT <- RejectionFreePT(steps=n_steps,temps=betas,GivenMatrix=MaxCut,h=h_local,swap_n=n_swaps)
      #Store results
      iter_find[i,m] <- simPT[[3]] #Extract the number of simulation in which the maximum was identified
      max_energy[i,m] <- simPT[[2]][1] #Extract the maximum energy identified
      check_rt <- RT_analysis(simPT[[4]]) #Analyze the replica swapping path
      swap_rates[i,m,] <- check_rt[[1]]$acc_rate #Extract the swap rates
      round_trips[i,m] <- check_rt[[3]] #Calculate the number of round trips
    }
  }
  file_id <- paste0(Sys.Date(),'_',seed_s)
  #Export the results
  write.csv(iter_find,file=paste0('./results/',file_id,'iteration_find.csv'), row.names = F)
  write.csv(max_energy,file=paste0('./results/',file_id,'max_energy_PT.csv'), row.names = F)
  write.csv(round_trips,file=paste0('./results/',file_id,'round_trips.csv'), row.names = F)
  write.csv(swap_rates,file=paste0('./results/',file_id,'swap_rates.csv'), row.names = F)
}





# for(i in 1:NumRep){
#   print(paste0('repetition #',i))
#   set.seed(i+321)
#   MaxCut <- RandomMaxCutMatrix()
#   #Baseline
#   set.seed(i+321) #Use the same seed so both algorithms start at the same spot
#   h_local <- h_baseline
#   simPT <- RejectionFreePT(steps=n_steps,temps=betas,GivenMatrix=MaxCut,h=h_local,swap_n=n_swaps)
#   iter_find[i,1] <- simPT[[3]] #Extract the number of simulation in which the maximum was identified
#   max_energy[i,1] <- simPT[[2]][1] #Extract the maximum energy identified
#   check_rt <- RT_analysis(simPT[[4]]) #Analyze the replica swapping path
#   swap_rates[i,1,] <- check_rt[[1]]
#   round_trips[i,1] <- check_rt[[3]]
#   #Max
#   set.seed(i+321) #Use the same seed so both algorithms start at the same spot
#   h_local <- h_m
#   simPT <- RejectionFreePT(steps=n_steps,temps=betas,GivenMatrix=MaxCut,h=h_local,swap_n=n_swaps)
#   iter_find[i,2] <- simPT[[3]]
#   max_energy[i,2] <- simPT[[2]][1]
#   check_rt <- RT_analysis(simPT[[4]])
#   swap_rates[i,2,] <- check_rt[[1]]
#   round_trips[i,2] <- check_rt[[3]]
#   #Squared root
#   set.seed(i+321) #Use the same seed so both algorithms start at the same spot
#   h_local <- h_2
#   simPT <- RejectionFreePT(steps=n_steps,temps=betas,GivenMatrix=MaxCut,h=h_local,swap_n=n_swaps)
#   iter_find[i,3] <- simPT[[3]]
#   max_energy[i,3] <- simPT[[2]][1]
#   check_rt <- RT_analysis(simPT[[4]])
#   swap_rates[i,3,] <- check_rt[[1]]
#   round_trips[i,3] <- check_rt[[3]]
#   #Balancing
#   set.seed(i+321) #Use the same seed so both algorithms start at the same spot
#   h_local <- h_balancing
#   simPT <- RejectionFreePT(steps=n_steps,temps=betas,GivenMatrix=MaxCut,h=h_local,swap_n=n_swaps)
#   iter_find[i,4] <- simPT[[3]]
#   max_energy[i,4] <- simPT[[2]][1]
#   check_rt <- RT_analysis(simPT[[4]])
#   swap_rates[i,4,] <- check_rt[[1]]
#   round_trips[i,4] <- check_rt[[3]]
# }
# 
# swap_r_summary <- rbind(swap_rates[,1,],swap_rates[,2,],swap_rates[,3,],swap_rates[,4,])


##### Exporting data #####
