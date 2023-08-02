#Function to simulate a rejection free Metropolis Hastings algorithm

##### PARAMETERS #####
#INPUT
# S : Finite state space
# initial : (element of S) Initial state of the chain
# B : (int) Burn-in to consider
# M : (int) Number of sample size to simulate (after burn-in)
# pi: (vector) target distribution
# Q: (function) proposal distribution
#OUTPUT
# tibble
# mul is the multiplicity
# sample is the state corresponding to that multiplicity

##### AID FUNCTION #####
#Function to calculate the probability of escaping a state
#INPUT
# x: state to calculate escape probability
# S: state space
# pi: target distribution
# Q: (function) proposal distribution
# adj: paramater for function Q
#OUTPUT
# Matrix
# 1st column is the neighbor
# 2nd column is the ptobability of moving from state X to the state in that row
alpha_p <- function(x,S,pi,Q, adj=0){
  n_x <- Q(x,S,adj)
  p_esc <- matrix(NA,ncol=2)
  for(y in n_x[,1]){
    n_y <- Q(y,S,adj)
    index_x <- which(n_y[,1]==x)
    index_y <- which(n_x[,1]==y)
    prob <- n_x[index_y,2]*min(1,(pi[y]*n_y[index_x,2])/(pi[x]*n_x[index_y,2]))
    p_esc <- rbind(p_esc,c(y,prob))
  }
  p_esc <- p_esc[-1,]
  return(p_esc) 
}

##### FUNCTION #####
mh_jump <- function(S,initial,B,M,pi,Q,adj=0){
  state <- c() #Initialize vector to store the sample
  multi <- c() #Initialize multiplicity chain
  state[1] <- initial #initial state
  i <- 1 #initialize the index for the state and multi vectors
  #burn-in steps  
  while(sum(multi)<B){
    esc_p <- alpha_p(state[i],S,pi,Q,adj) #Obtain the escape probability and the transition prob.
    s_rep <- 1 + rgeom(1,sum(esc_p[,2])) #Get the multiplicity for the current state
    if((sum(multi) + s_rep)>B){ #If the chain would remain more than the burn-in
      multi[i] <- B - sum(multi) #Stay in that state the remaining number of steps
      state[i+1] <- state[i] #The chain "jump" to the same state
      i <- i+1 #update the index
      break #End the burn-in process
    }else{
      multi[i] <- s_rep #assign the value to the multiplicity list
      state[i+1] <- sample(esc_p[,1],1,prob=esc_p[,2]) #Get the next state
    }
    i <- i+1
  } #End burn-in steps
  burnin <- length(multi) #this identifies at which index we got the burn-in sample
  #Sampling steps
  while(sum(multi)<(B+M)){
    esc_p <- alpha_p(state[i],S,pi,Q,adj) #Obtain the escape probability and the transition prob.
    s_rep <- 1 + rgeom(1,sum(esc_p[,2])) #Get the multiplicity for the current state
    if((sum(multi) + s_rep)>=(B+M)){ #If the chain would remain more than the burn-in
      multi[i] <- (B+M) - sum(multi) #Stay in that state the remaining number of steps
      break #End the sampling process
    }else{
      multi[i] <- s_rep #assign the value to the multiplicity list
      state[i+1] <- sample(esc_p[,1],1,prob=esc_p[,2]) #Get the next state
    }
    i <- i+1
  }
  return(tibble(mul=multi[-(1:burnin)], sample=state[-(1:burnin)])) #Return ignoring burnin
}

# #### Example #####
# library(dplyr)
# set.seed(123)
# S <- 1:5 #State space
# pi <- exp(S/2) #Proportions of target distribution
# #Q(i,j) ∝ 1
# Q_unif <- function(i,S, adj=0){
#   if(!(i %in% S)){print(paste("State",i,"is not in S"));return(NA)}else{#Check that state i is in S
#     if(adj < 0){print("adj parameter must be non-negative"); return(NA)}else{
#       if(adj == 0){ #Consider all neighbors
#         neighbors <- S[S!=i]
#       }else{ #Consider the number of adjacent neighbors defined by adj
#         index <- which(i==S) - 1 #index of the current state
#         adj_index <- c(-adj:-1,1:adj)
#         neighbors <- S[(adj_index + index)%%(length(S)) + 1]
#         neighbors <- unique(neighbors[neighbors != i]) #In case the # of adjacent neighbors is too big
#       }
#       #uniform proposal distribution
#       n_size <- length(neighbors)
#       prob <- rep(1/n_size,n_size)
#       return(cbind(neighbors,prob))
#     }
#   }
# }
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
