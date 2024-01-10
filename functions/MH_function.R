#Function to simulate 1 chain with Metropolis Hastings algorithm


##### Parameters #####
#INPUT
# S : Finite state space
# initial : (element of S) Initial state of the chain
# B : (int) Burn-in to consider
# M : (int) Number of sample size to simulate (after burn-in)
# pi: (vector) target distribution
# Q: (function) proposal distribution
#OUTPUT
# vector of size M 

##### AID FUNCTION #####
source('functions/Q_unif.R')

##### Function #####
mh_sim <- function(S,initial,B,M,pi,Q){
  state <- c() #Initialize vector to store the sample
  state[1] <- initial #initial state
  for(i in 1:(M+B)){ #burn-in steps + sample size
    n_i <- Q(state[i],S) #Considering the proposal distribution Q
    j <- sample(n_i[,1], 1,prob=n_i[,2]) #sample from the neighbors
    n_j <- Q(j,S) #Neighbors of state j (the proposed state)
    index_i <- which(n_j[,1]==state[i])
    index_j <- which(n_i[,1]==j)
    if (runif(1)< (pi[j]*n_j[index_i,2]/(pi[state[i]]*n_i[index_j,2]))){ #Check probability to accept or reject
      state[i+1] <- j
    }else{state[i+1] <- state[i]}
  } #End burn-in steps
  return(state[-(1:(B+1))]) #Exclude the initial state and the B steps of burn-in
}

#### Example #####
{
# S <- 1:5 #State space
# pi <- exp(S/2) #Proportions of target distribution
# ex <- mh_sim(S,1,50,10000,pi,Q_unif)
# compare <- tibble(sim=table(ex)/length(ex), target=pi/sum(pi),diff=abs(target-sim))
# #TVD
# 0.5*sum(compare$diff)
}