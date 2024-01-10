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
