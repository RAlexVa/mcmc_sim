# Function to do 1 iteration of IIT
# We consider uniform proposal distribution over the specified neighbors
# Input:
# pi: Vector with target distribution
# nbr: Vector with index of specified neighbors
# x: current state's index
# h: balancing function
# Output: New state choosen proportionally and previous state's weight
IITupdate <- function(pi, nbr,x,h){
  N <- length(nbr)
  probs <- c() #Vector to store weights
  for(i in 1:N){
    probs[i] <- h(pi[nbr[i]]/pi[x])
  }
  
  U = runif(N)
  D = -log(U)/probs #Selecting a state proportional to the Prob vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- nbr[index]
  weight <- 1/mean(probs)
  
  return(list(Xnew,weight))
}


#Function to simulate a IIT
# Considering a uniform proposal distribution over d adjacent neighbors
##### PARAMETERS #####
#INPUT
# initial : (index of element of S) Initial state's index of the chain
# M : (int) Number of jumps the chain performs
# pi: (vector) target distribution
# d: number of adjacent neighbors on the left and right side
# h: balancing function to use
#OUTPUT
# tibble
# w is the weight
# x is the state


##### FUNCTION #####
IIT_opt <- function(initial,M,pi,d,h){
  #Identify optima to halt loop
  optima <- which(pi==max(pi))
  n_opt <- length(optima)
  #Common setup
  N <- length(pi)
  weights <- c()
  state <- c()
  state[1] <- initial
  #Split depending on the number of optima

    for(i in 2:M){
      Xnow <- state[i-1]
      neighbors <- c(Xnow+seq(-d,-1),Xnow+seq(1,d))%%N #Neighbors excluding current state
      neighbors <- ifelse(neighbors==0,length(pi),neighbors) #Replace 0 mod S with the last element of S
      temp <- IITupdate(pi,neighbors,Xnow,h)
      state[i] <- temp[[1]] #Index of the new state
      weights[i-1] <- temp[[2]] #Weight of the previous state
      if(n_opt==1){ #If there's only 1 optimum, algorithm halts when it finds it
      if(state[i]==optima){weights[i] <- 0;break;} #Finish algorithm when you find the optimum
      #If you start in the optimum you still do at least 1 step.
      }else{#If there's more than 1 optima
        if(state[i]%in%optima){
          ind <- which(optima==state[i]) #Identify which optimum was visited
          optima <- optima[-ind] #Update optima vector
          n_opt <- length(optima)
        }
        
      }
    }
    return(list(state,weights))
  }
