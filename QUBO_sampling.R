####
# Code to sample from a QUBO problem with dimension at most 16
# Mostly using prvious code from Sigeng Chen
####
rm(list=ls())
library(tidyverse)
library(here)
###### Functions
### Create a random QUBO matrix
RandomQuboMatrix <- function(std)     
{
  L = matrix(rnorm(N * N, mean = 0, sd = std), nrow = N)
  L[lower.tri(L)] <- 0
  return(L)
}

### Initialize the first state
RandomStart <- function()
{
  X = rbinom(N, 1, 0.5)
  return(X)
}

### Given a vector and a matrix, return the energy
Energy <- function(X, GivenMatrix)
{
  return(t(X) %*% GivenMatrix %*% X)
}

### Return the exponential of the energy
ExpEnergy <- function(X, GivenMatrix)
{
  return(exp(Energy(X, GivenMatrix)))
}

### Map the N dimensional binary vector to its corresponding number (convert from base 2 to base 10)
VectorToNumber <- function(X)
{
  return(strtoi(paste(X, collapse=''), base=2)+1)
  #We add 1 since the 0 (in binary) is mapped to the 1 (in decimal) due to R vectors starting at index 1
}

### Map the number to the N dimensional binary vector (convert from base 10 to base 2)
NumberToVector <- function(number){
  if(number <= 0 || number > 2^N){
    print("Error")
    return(0)
  }
  X <- paste(rev(as.integer(intToBits(number - 1))), collapse='')
  #We substract 1 since the 0 (in binary) is mapped to the 1 (in decimal) due to R vectors starting at index 1
  entries <- nchar(X) #should be 32 due to the OS
  #trim the string to only have N entries
  vector <- substr(X,entries-N+1,entries)
  # Validation below is not needed since we already have an if excluding numbers > 2^N
  # remnant <- substr(X,1,entries-N)
  # 
  # if(VectorToNumber(remnant)!=0){
  #   print("Error")
  #   return(0)
  # }
  return(as.numeric(strsplit(vector,'')[[1]]))
}

### Calculate the true distribution given a matrix
TrueDistribution <- function(TrueMatrix)
{
  if(N > 16)
  {
    print("N too large; infeasible")
    return(0)
  }
  Distribution <- rep(0, 2^N)
  for(i in 1 : 2^N)
  {
    X = NumberToVector(i)
    Distribution[i] = ExpEnergy(X, TrueMatrix)
  }
  DistributionSum = sum(Distribution)
  Distribution = Distribution / DistributionSum
  return(Distribution)
}

### Calculate the Total variation distance of 2 probability vectors.
TVD <- function(D1, D2)
{
  nn = length(D1)
  if(length(D2) != nn)
  {
    print("Not the Same Length")
    return(1)
  }
  else
  {
    ans = 0
    for(i in 1 : nn)
    {
      ans = ans + abs(D1[i] - D2[i]) / 2
    }
    return(ans)
  }
}


### Given a state, provides the weight of the current state and the next accepted state
RejectionFreeUpdate <- function(GivenMatrix, XNow, h=h_min)
{
  EnergyNow =  Energy(XNow, GivenMatrix) #Energy at current state
  logProb = rep(0, N)
  #Check energy of neighbors
  #Neighbor set is composed by fliping 1 coordinate
  for(r in 1:N)
  {
    XNew = XNow 
    XNew[r] = 1 - XNew[r]
    EnergyNew = Energy(XNew, GivenMatrix)
    logProb[r] = h(EnergyNew, EnergyNow) #probability of accepting the flipping of coordinate r
  }
  U = runif(N)
  logD = log(- log(U)) - logProb #Selecting a state proportional to the Prob vector
  #We are simply using the log of the original formula for simplicity.
  index = which.min(logD) #Select index to flip coordinate
  XNew = XNow
  XNew[index] = 1 - XNew[index] #Select the new state
  
  #The weight is 1/Z_h 
  #Z_h is the sum, for all neighbors, of q(y|x) h(pi(y)/pi(x)) since we're using a uniform q(y|x)=1/N we can use the mean 
  weight = 1/mean(exp(logProb))
  #Below is the code to simulate the number of rejections
  # NumM = 1 + rgeom(1, mean(exp(logProb)))
  # if(is.na(NumM))
  # {
  #   print("Error in rgeom RF part")
  #   print(mean(exp(logProb)))
  #   return(list(XNew, 1000))
  # }
  # return(list(XNew, NumM))
  return(list(XNew, weight))
}

### Function to run a rejection free simulation
RejectionFreeSampling<- function(GivenMatrix, NumSamples, h)
{
  X = RandomStart() #uniform random initial point
  Distribution <- rep(0, 2^N) #Vector to store the estimation
  for(i in 1 : NumSamples) #Obtain NumSamples weights
  {
    Ans = RejectionFreeUpdate(GivenMatrix, X, h) #Do 1 step of rejection free update
    w = Ans[[2]]
    index = VectorToNumber(X)
    Distribution[index] = Distribution[index] + w #Update the estimated weight
    X = Ans[[1]] #Update the state
  } 
  DistributionSum = sum(Distribution)
  Distribution = Distribution / DistributionSum #Calculated the estimated density (using the weights)
  return(Distribution)
}

## Function to run many rejection free simulations.
SimulationGivenNumSamples <- function(TrueMatrix, NumRep, NumSamples,h)
{
  TD <- TrueDistribution(TrueMatrix) #True distribution
  ptm <- proc.time()
  Temp <- rep(0,NumRep)
  for(k in 1: NumRep)
  {
    Distribution <- RejectionFreeSampling(TrueMatrix, NumSamples,eval(parse(text=h)))
    #Here the parameter h is just a character with the name of a function
    #We use it like this so it's easier to create a vector with different functions to use and also
    #to present the name of the function in the summary
    Temp[k] =  TVD(TD, Distribution)
  }

  Ans <- tibble(TVD = median(Temp), 
                Time=(proc.time() - ptm)[1] / NumRep,
                NumSamples=NumSamples,
                Bfunction =h)
  return(Ans)
}

Simulation <- function(TrueMatrix, NumRep, NumSamplesList,b_functions)
{
  if(length(NumSamplesList) != length(b_functions)){print('List of samples and list of functions have different size'); break;}
  Table <- NULL
  for(i in 1:length(NumSamplesList))
  {
    print(paste0("Now, It's ", NumSamplesList[i]," samples for the function ",b_functions[i]))
    Ans = SimulationGivenNumSamples(TrueMatrix, NumRep, NumSamplesList[i], b_functions[i])
    Table <- rbind(Table, Ans)
  }
  return(Table)
}

### Balancing functions to be used in the log-probabilities 
h_sq <- function(new,current){
  return((new-current)/2)}
h_min <- function(new,current){
  return(min(0,new-current))}
h_max <- function(new,current){
  return(max(0,new-current))}
h_plus <- function(new,current){
  return(max(0,new-current))}

##### Inputs
semilla <- 202309
set.seed(semilla)

N = 8 # Dimension
NumSamplesList <- c(rep((1:5)*1000,3))#c(c(2:15) * 500,(8:10)*1000)
NumRep = 1000 # Number of times the experiment is repeated to get estimates of TVD and time
b_functions <- c(rep('h_sq',5),rep('h_min',5),rep('h_max',5))
#b_functions <- list(h_sq,h_sq,h_sq,h_max,h_max,h_max,h_min,h_min,h_min)

##### Simulations

Q <- RandomQuboMatrix(3)
sim <- Simulation(Q, NumRep, NumSamplesList,b_functions)
write.csv(sim,here('simulation results',paste0(Sys.Date(),'_b_f_compare.csv')), row.names = F)
report <- c(paste0('dimensons = ',N),paste0('seed = ',semilla),paste0('repetitions = ',NumRep))
write.table(report,here('simulation results',paste0(Sys.Date(),'_b_f_report.csv')))
##### Results