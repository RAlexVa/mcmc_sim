#Function to simulate a rejection free Metropolis Hastings algorithm

##### PARAMETERS #####
#INPUT
# S : State space
# initial : Initial state of the chain
# B : Burn-in to consider
# M : Number of sample size to simulate (after burn-in)
# pi: target distribution
# Q: proposal distribution
#OUTPUT
# tibble
# mul is the multiplicity
# sample is the state corresponding to that multiplicity

##### AID FUNCTION #####
#Function to calculate the probability of escaping a state
#INPUT
# pi: target distribution
# Q: proposal distribution
#OUTPUT
# list
# alpha is the probability of escaping each state
# matrix is the new transition probabilities (considering that the chain always escapes from the current state)
alpha_p <- function(pi, Q){
  #pi doesn't need to be normalized
  qQ <- t(Q)/Q #Quotient of Q(j,i)/Q(i,j)
  diag(qQ) <- 0 #Set the diagonal to 0
  qQ <- sweep(qQ,MARGIN=2,pi,`*`) #Multiply by pi(i)
  qQ <- sweep(qQ,MARGIN=1,pi,`/`) #Divide by pi(j)
  qQ <- apply(qQ,c(1,2),function(x){min(1,x)}) #Getting the minimum of 1 and the entries on qQ
  
  qQ <- Q*qQ #multiply by Q(i,j)
  # Return the sum of rows = probability of escape
  # Return also the probabilities of transitioning
  return(list("alpha" = rowSums(qQ),"matrix" = qQ)) 
}

##### FUNCTION #####
mh_jump <- function(S,initial,B,M,pi,Q){
  state <- c() #Initialize vector to store the sample
  multi <- c() #Initialize multiplicity chain
  state[1] <- initial #initial state
  esc_p <- alpha_p(pi,Q) #Obtain the escape probability and the transition prob.
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
