####
# Code to optimize a QUBO problem with dimension N using parallel tempering
# We use uniform proposal distributions
# We can use different balancing functions
# Mostly using prvious code from Sigeng Chen
####


rm(list=ls())
library(tidyverse)
library(here)
library(igraph) #To create random adjacency matrix (for MaxCut)
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
  diag(adjmat) <- -1* adjmat%*%rep(1,N) #to obtain the values for the diagonals
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
ExpEnergy <- function(X, GivenMatrix)
{
  return(exp(Energy(X, GivenMatrix)))
}


### Given a state, provides the weight of the current state and the next accepted state
RF_Update_temp <- function(GivenMatrix, XNow, h=h_min, beta=1)
{
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

## h is a vector, we define 1 balancing function for each replica (each temperature)
RejectionFreePT <- function(steps,temps,GivenMatrix,h){
  print('Starting RF-PT')
  # Create an array to store the states visited
  n_rep <- length(temps) #Number of replicas
  current_state <- matrix(NA,nrow=N,ncol=n_rep) #each column is a replica
  max_state <- current_state #array to store the state corresponding to the maximum
  max_energy <- rep(-Inf,n_rep) #vector to store the max value found
  iter_found <- 0 #To record the iteration in which the maximum value was found
  #Randomly and independently initialize all replicas
  for (i in 1:n_rep){
    current_state[,i] <- RandomStart()
  }
  for(i in 1:steps){ #Number of steps that each replica will take
    if(i %% floor(steps/4) ==0){
      print(paste0('iteration # ',i," = ",round(100*i/steps,0),"%")) 
    }

### First step to update the replica
    for(k in 1:n_rep){
      #We provide a vector with the name of the balancing function to use and then
      #We use eval(parse(text)) to pass the actual function (not the string) to the RF_Update_temp function
      Ans <- RF_Update_temp(GivenMatrix, current_state[,k], h=eval(parse(text=h[k])), beta=temps[k])
      current_state[,k] = Ans[[1]]
      Value = Ans[[2]]
      #Print = c(Print, Value)
      if(max_energy[k] < Value) #New energy is bigger
      {
        max_energy[k] <-  Value
        max_state[,k] <-  current_state[,k]
        iter_found <- i 
      }
    }
### Then try a swap of replicas
    
    #identify the replicas to swap (even or odd)
    swap_t <- 1:(n_rep-1)
    swap_t <- swap_t[swap_t%%2 == (i%%2)]
    
    #Try to swap all replicas (that mathces the parity of i) with the replica above them
    for(t in swap_t){
      EnergyCurT =  Energy(current_state[,t], GivenMatrix,beta=temps[t]) + Energy(current_state[,t+1], GivenMatrix,beta=temps[t+1])
      EnergyNewT = Energy(current_state[,t], GivenMatrix,beta=temps[t +1]) + Energy(current_state[,t+1], GivenMatrix,beta=temps[t])
      accept_logp <- EnergyNewT - EnergyCurT
      if(log(runif(1))<accept_logp){#accept the swap of replicas
        temp <- current_state[,t]
        current_state[,t] <- current_state[,t+1]
        current_state[,t+1] <- temp
      }
    }
  }
  
  return(list(max_state,max_energy,iter_found))
}

### some balancing functions
### Balancing functions to be used in the log-probabilities 
h_sq <- function(new,current){
  return((new-current)/2)}
h_min <- function(new,current){
  return(min(0,new-current))}
h_max <- function(new,current){
  return(max(0,new-current))}

###### Simulations ########
N=50 #Dimension of the QUBO problem
temps <- c(1,1 *1.29^(1:10)) #temperatures
betas <- 1/temps #Inverse temperatures
h_baseline <- rep('h_min',length(betas))
h_balancing <- c(rep('h_sq',5),rep('h_max',6))

# set.seed(321)
# MaxCut <- RandomMaxCutMatrix()
# set.seed(321)
# baseline <- RejectionFreePT(1000,betas,MaxCut,h_baseline)
# set.seed(321)
# compare <- RejectionFreePT(1000,betas,MaxCut,h_balancing)

NumRep <- 100
iter_find <- matrix(NA,nrow=NumRep,ncol=2)
colnames(iter_find) <- c('baseline','compare')
max_energy <- matrix(NA,nrow=NumRep,ncol=2)
colnames(max_energy) <- c('baseline','compare')

for(i in 1:NumRep){
  print(paste0('repetition #',i))
  set.seed(i+3)
  MaxCut <- RandomMaxCutMatrix()
  set.seed(i+3) #Use the same seed so both algorithms start at the same spot
  baseline <- RejectionFreePT(300,betas,MaxCut,h_baseline)
  iter_find[i,1] <- baseline[[3]]
  max_energy[i,1] <- baseline[[2]][1] 
  set.seed(i+3) #Use the same seed so both algorithms start at the same spot
  compare <- RejectionFreePT(300,betas,MaxCut,h_balancing)
  iter_find[i,2] <- compare[[3]]
  max_energy[i,2] <- compare[[2]][1] 
}

###### Plots ########
data <- as.data.frame(cbind(iter_find,max_energy))
colnames(data) <- c('baseline_iter','compare_iter','baseline_max','compare_max')
data$run <- 1:nrow(data)
data$best_iter <- ifelse(data$baseline_iter>data$compare_iter,1,ifelse(data$baseline_iter<data$compare_iter,0,0.5))
data$best_energy <- ifelse(data$baseline_max<data$compare_max,1,ifelse(data$baseline_max>data$compare_max,0,0.5))
data <- data |> pivot_longer(-run,names_to='variable',values_to='value')
data <- data |> separate(variable,c('model','variable'))


###
data |> filter(model=='best') |> group_by(variable) |> summarise(porcentage=sum(value)/NumRep)

data |> filter(model=='best', variable=='iter') |> 
  group_by(value) |> 
  summarise(count=n()) |> ungroup() |> 
  ggplot(aes(x=as.factor(value), y=count))+
    geom_col(fill='dodgerblue3')+
  geom_text(aes(label = count), vjust = -0.5)+
  labs(x='0=baseline, 0.5=tie, 1=balancing',
       title='Which found the max solution faster?')

data |> filter(model=='best', variable=='energy') |> 
  group_by(value) |> 
  summarise(count=n()) |> ungroup() |> 
  ggplot(aes(x=as.factor(value), y=count))+
  geom_col(fill='dodgerblue3')+
  geom_text(aes(label = count), vjust = -0.5)+
  labs(x='0=baseline, 0.5=tie, 1=balancing',
       title='Which found the best solution?')

data |> filter(model=='best', variable=='iter') |>ggplot(aes(x=as.factor(value)))+
  geom_bar(fill='dodgerblue3')+
  labs(x='0=baseline, 0.5=tie, 1=balancing',
       title='Which found the max solution faster?')

data |> filter(model=='best', variable=='energy') |>ggplot(aes(x=as.factor(value)))+
  geom_col(fill='dodgerblue3')+
  labs(x='0=baseline, 0.5=tie, 1=balancing',
       title='Which found the best solution?')
#This doesn't look that gives any improvement
###

data |> filter(variable=='max') |> 
  ggplot(aes(x=model,y=value))+geom_boxplot()

data |> filter(variable=='max') |> 
  ggplot(aes(x=run,y=value,color=model))+geom_point()
###### Saving results ########
write.csv(data,here('simulation results',paste0(Sys.Date(),'_PT_compare.csv')), row.names = F)
