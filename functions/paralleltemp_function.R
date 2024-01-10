#Function to simulate a parallel tempered metropolis hastings algorithm

##### PARAMETERS #####
#INPUT
# S : State space
# initial : Initial state of the chain
# B : Burn-in to consider
# M : Number of sample size to simulate (after burn-in)
# pi: target distribution
# Q: proposal distribution
# temps: vector with the temperatures to use
#OUTPUT
# matrix
# state contais the sample generated, ncol=number of temperatures, nrow=M

parallel_temp <- function(S, initial, B, M, pi, Q, temps){
  state <- matrix(NA,nrow=B+M,ncol=length(temps))
  state[1,] <- initial
  
  for(i in 1:(B+M-1)){#Run on the rows of the matrix
    for(t in 1:length(temps)){ #Run on the columns (for each temperature)
      new_s <- sample(S, 1,prob=Q[state[i,t],]) #propose a new state using Q
      tau <- temps[t] #current temperature
      if (runif(1)< (pi[new_s]^(1/tau)*Q[new_s,state[i,t]])/(pi[state[i,t]]^(1/tau)*Q[state[i,t],new_s])) { 
        state[i+1,t] <- new_s #Accept
      }else{state[i+1,t] <- state[i,t]} #Reject
    }
    #once the update happened for each temperature
    #We propose a swap in temperatures
    new_t <- sample(1:length(temps),2) #Choose 2 temperatures to swap
    t1 <- new_t[1]
    t2 <- new_t[2]
    #Probability of accepting the swap of temperatures
    prob_t <- (pi[state[i+1,t2]]^(1/temps[t1]) * pi[state[i+1,t1]]^(1/temps[t2]))/(pi[state[i+1,t1]]^(1/temps[t1]) * pi[state[i+1,t2]]^(1/temps[t2]))
    if(runif(1)< prob_t){
      #accept the swap
      temp <- state[i+1,t1]
      state[i+1,t1] <- state[i+1,t2]
      state[i+1,t2] <- temp
    }
  }
  return(state[-(1:B),])
}

##### Example #####
# S <- 1:5 #State space
# pi <- exp(S/2) #Proportions of target distribution
# #Q(i,j) ∝ 1
# Q <- matrix(1/length(S),ncol=length(S),nrow=length(S))
# diag(Q) <- 0
# temps <- c(1,2.5,3.5,5)
# ex <- parallel_temp(S,1,50,10000,pi,Q,temps)
# 
# probs <- matrix(NA,nrow=length(S),ncol=length(temps))
# for(i in 1:length(temps)){
#   probs[,i] <- (table(c(ex[,i],S))-1)/length(ex[,i])
# }
# 
# compare <- tibble(sim=probs[,1], target=pi/sum(pi),diff=abs(target-sim))
# #TVD
# 0.5*sum(compare$diff)
