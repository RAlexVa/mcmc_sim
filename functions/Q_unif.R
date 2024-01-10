##### AID FUNCTION #####
# Return uniform distribution for the specified neighbors
#Q(i,j) ∝ 1
##### PARAMETERS #####
#INPUT
# S : Finite state space
# initial : (element of S) Initial state of the chain
# adj or nbr, depending on the option we want
#OUTPUT
# matrix
# column 1 is the neighbor
# column 2 is the uniform probability of each neighbor calculated as 1/neighbor set size


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