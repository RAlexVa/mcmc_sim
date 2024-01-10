#Function to simulate unbiased Partial Neighbor Search algorithm

##### PARAMETERS #####
#INPUT
# S : Finite state space
# initial : (element of S) Initial state of the chain
# B : (int) Burn-in to consider
# M : (int) Number of sample size to simulate (after burn-in)
# pi: (vector) target distribution
# Q: (function) proposal distribution
# L: (int) Number of samples to obtain for each PNSet
# nbr_func: (function) Function to select the Partial Neighbors of a state
#OUTPUT
# tibble
# mul is the multiplicity
# sample is the state corresponding to that multiplicity

##### AID FUNCTIONS #####
source('functions/Q_unif.R')
source('functions/alpha_pns.R')

##### FUNCTION #####
PNS_unbiased <- function(S,initial,M,L,pi,Q,nbr_func){
  state <- c() #Initialize vector to store the sample
  multi <- c() #Initialize multiplicity chain
  state[1] <- initial #initial state
  i <- 0 #initialize the index for the state and multi vectors
  n_count <- 0 #change of neighbor count
  last <- FALSE #flag for the last iteration
  
  while(sum(multi)<M){ #loop for the sample
    l_count <- L #Reseting the L to change neighborhoods
    if(sum(multi) + L >=M){l_count <- M-sum(multi); last <- TRUE}
    while(l_count >0){ #loop for changing PNSets
      i <- i+1
      neighbors <- nbr_func(state[i]) #calculate neighborhoods of current state
      PNSet <- neighbors[[(n_count%%length(neighbors))+1]] #Select the neighborhood depending on L
      esc_p <- alpha_pns(state[i],S,pi,Q,nbr=PNSet) #Obtain the escape probability and the transition prob for the partial neighbors
      s_rep <- 1 + rgeom(1,sum(esc_p[,2])) #Get the multiplicity for the current state
      if(s_rep>=l_count){ #Stop when you get L samples from that neighborhood
        multi[i] <- l_count #Stay in that state the remaining number of steps
        if(!last){state[i+1] <- state[i]}#The chain "jump" to the same state
        break #End the sampling process
      }else{
        multi[i] <- s_rep #assign the value to the multiplicity list
        state[i+1] <- sample(esc_p[,1],1,prob=esc_p[,2]) #Get the next state
        l_count <- l_count - s_rep
      }
    }
    #After finishing this loop we'll reset L
    n_count <- n_count + 1 #Count how many times we've restarted the L
  }
  return(tibble(mul=multi, sample=state)) 
}

#### Example #####
{
  # library(dplyr)
  # set.seed(123)
  # S <- 1:5 #State space
  # pi <- exp((S-1)/2) #Proportions of target distribution
  # pns_func <- function(x){
  #   indx <- which(S==x)
  #   n1 <- c(-1,0,1)
  #   n2 <- c(-2,0,2)
  #   n1 <- unique(S[((indx-1)+n1)%%length(S)+1])
  #   n2 <- unique(S[((indx-1)+n2)%%length(S)+1])
  #   if(!all(S%in%unique(c(n1,n2)))){
  #     print("Neighborhoods not covering S")
  #     return(NA)
  #   }else{
  #     return(list(n1,n2))
  #   }
  # }
  # ex <- PNS_unbiased(S,1,10000,100,pi,Q_unif,pns_func)
  # est_prob <- ex |> #estimate the probability
  #   group_by(sample) |>
  #   summarize(fre = sum(mul)) |>
  #   mutate(prob = fre/sum(fre)) |>
  #   pull(prob)
  # 
  # compare <- tibble(sim=est_prob, target=pi/sum(pi),diff=abs(target-sim))
  # #TVD
  # 0.5*sum(compare$diff)
}

