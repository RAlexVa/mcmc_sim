#Function to simulate 1 chain with Metropolis Hastings algorithm


##### Parameters #####
#INPUT
# S : Finite state space
# initial : (element of S) Initial state of the chain
# B : (int) Burn-in to consider
# M : (int) Number of sample size to simulate (after burn-in)
# pi: (vector) target distribution
# Q: (Matrix) proposal distribution
#OUTPUT
# vector of size M 


##### Function #####
mh_sim_M <- function(S,initial,B,M,pi,Q){
  state <- c() #Initialize vector to store the sample
  state[1] <- initial #initial state
  for(i in 1:(M+B)){ #burn-in steps + sample size
    new_s <- sample(S, 1,prob=Q[state[i],]) #from the neighbors sample using Q
    if (runif(1)< (pi[new_s]*Q[new_s,state[i]])/(pi[state[i]]*Q[state[i],new_s])) { #Check probability to accept or reject
      state[i+1] <- new_s
    }else{state[i+1] <- state[i]}
  } #End burn-in steps
  return(state[-(1:(B+1))]) #Exclude the initial state and the B steps of burn-in
}

##### Example #####
# S <- 1:5 #State space
# pi <- exp(S/2) #Proportions of target distribution
# #Q(i,j) ∝ 1
# Q <- matrix(1/length(S),ncol=length(S),nrow=length(S))
# diag(Q) <- 0
# ex <- mh_sim_M(S,1,50,10000,pi,Q)
# compare <- tibble(sim=table(ex)/length(ex), target=pi/sum(pi),diff=abs(target-sim))
# #TVD
# 0.5*sum(compare$diff)