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

##### AID FUNCTION #####
# return uniform distribution for the specified neighbors
#Q(i,j) ∝ 1
#Option 1: Use parameter nbr to define specific neighbors to consider 
#Option 2: Use adj to define the number of adjacent neighbors (to the current state i) to consider
Q_unif <- function(i,S, adj=0, nbr=NA){
  #Returns neighbors and uniform distribution over those neighbors
  if(!all(is.na(nbr))){
    if(!all(nbr %in% S)){
      print("some of the specified neighbors are not in S");return(NA)
    }else{
      neighbors <- nbr[nbr!=i]
      n_size <- length(neighbors)
      prob <- rep(1/n_size,n_size)
      return(cbind(neighbors,prob))
    }
    
  }else{
    if(!(i %in% S)){print(paste("State",i,"is not in S"));return(NA)}else{#Check that state i is in S
      if(adj < 0){print("adj parameter must be non-negative"); return(NA)}else{
        if(adj == 0){ #Consider all neighbors
          neighbors <- S[S!=i]
        }else{ #Consider the number of adjacent neighbors defined by adj
          index <- which(i==S) - 1 #index of the current state
          adj_index <- c(-adj:-1,1:adj)
          neighbors <- S[(adj_index + index)%%(length(S)) + 1]
          neighbors <- unique(neighbors[neighbors != i]) #using unique in case the # of adjacent neighbors is too big
        }
        #uniform proposal distribution
        n_size <- length(neighbors)
        prob <- rep(1/n_size,n_size)
        return(cbind(neighbors,prob))
      }
    }
  }
}


##### AID FUNCTION #####
#Function to calculate the probability of escaping a state in the PNS setting
#INPUT
# x: state to calculate escape probability
# S: state space
# pi: target distribution
# Q: (function) proposal distribution
# nbr: Neighbors to consider
#OUTPUT
# Matrix
# 1st column is the neighbor
# 2nd column is the ptobability of moving from state X to the state in that row
alpha_pns <- function(x,S,pi,Q,nbr=S){
  p_esc <- matrix(NA,ncol=2) #To store probability of escaping
  n_x <- Q(x,S,adj=0) #Obtain transition probabilities considering all neighbors
  n_x <- n_x[n_x[,1]%in%nbr,] #Select only the specified neighbors
  n_x[,2] <- n_x[,2]/sum(n_x[,2]) #Re-distribute the probability considering only the specified neighbors
  #Previous line is done to create the new Q for the PNSet
  for(y in intersect(n_x[,1],nbr)){
    n_y <- Q(y,S,adj=0)
    n_y <- n_y[n_y[,1]%in%c(nbr,x),] #Select only the specified neighbors and the current state X
    n_y[,2] <- n_y[,2]/sum(n_y[,2]) #Re-distribute the probability considering only the specified neighbors
    #Previous line is done to create the new Q for the PNSet
    index_x <- which(n_y[,1]==x)
    index_y <- which(n_x[,1]==y)
    prob <- n_x[index_y,2]*min(1,(pi[y]*n_y[index_x,2])/(pi[x]*n_x[index_y,2]))
    p_esc <- rbind(p_esc,c(y,prob))
  }
  p_esc <- p_esc[-1,]
  return(p_esc) 
}

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

