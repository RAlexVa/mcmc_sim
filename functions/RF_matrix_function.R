#Function to simulate a rejection free Metropolis Hastings algorithm

##### PARAMETERS #####
#INPUT
# S : Finite state space
# initial : (element of S) Initial state of the chain
# B : (int) Burn-in to consider
# M : (int) Number of sample size to simulate (after burn-in)
# pi: (vector) target distribution
# Q: (matrix) proposal distribution
#OUTPUT
# tibble
# mul is the multiplicity
# sample is the state corresponding to that multiplicity
##### AID FUNCTION #####
source('functions/alpha_p_M.R')
##### FUNCTION #####
mh_jump <- function(S,initial,B,M,pi,Q){
  state <- c() #Initialize vector to store the sample
  multi <- c() #Initialize multiplicity chain
  state[1] <- initial #initial state
  esc_p <- alpha_p_M(pi,Q) #Obtain the escape probability and the transition prob.
  i <- 1 #initialize the index for the state and multi vectors
  #burn-in steps  
  while(sum(multi)<B){
    s_rep <- 1 + rgeom(1,esc_p$alpha[state[i]]) #Get the multiplicity for the current state
    if((sum(multi) + s_rep)>B){ #If the chain would remain more than the burn-in
      multi[i] <- B - sum(multi) #Stay in that state the remaining number of steps
      break #End the burn-in process
    }else{
      multi[i] <- s_rep #assign the value to the multiplicity list
      state[i+1] <- sample(S,1,prob=esc_p$matrix[state[i],]) #Get the next state
    }
    i <- i+1
  } #End burn-in steps
  burnin <- length(multi) #this identifies at which index we got the burn-in sample
  #Sampling steps
  while(sum(multi)<(B+M)){
    s_rep <- 1 + rgeom(1,esc_p$alpha[state[i]]) #Get the multiplicity for the current state
    if((sum(multi) + s_rep)>=(B+M)){ #If the chain would remain more than the burn-in
      multi[i] <- (B+M) - sum(multi) #Stay in that state the remaining number of steps
      break #End the sampling process
    }else{
      multi[i] <- s_rep #assign the value to the multiplicity list
      state[i+1] <- sample(S,1,prob=esc_p$matrix[state[i],]) #Get the next state
    }
    i <- i+1
  }
  return(tibble(mul=multi[-(1:burnin)], sample=state[-(1:burnin)])) #Return ignoring burnin
}

##### Example #####
# S <- 1:5 #State space
# pi <- exp(S/2) #Proportions of target distribution
# #Q(i,j) ∝ 1
# Q <- matrix(1/length(S),ncol=length(S),nrow=length(S))
# diag(Q) <- 0
# ex <- mh_jump(S,1,50,10000,pi,Q)
# est_prob <- ex |> #estimate the probability
#   group_by(sample) |> 
#   summarize(fre = sum(mul)) |> 
#   mutate(prob = fre/sum(fre)) |> 
#   pull(prob) 
# 
# compare <- tibble(sim=est_prob, target=pi/sum(pi),diff=abs(target-sim))
# #TVD
# 0.5*sum(compare$diff)
