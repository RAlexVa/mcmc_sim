##### AID FUNCTION #####
#Function to calculate the probability of escaping a state
#INPUT
# x: state to calculate escape probability
# S: state space
# pi: target distribution
# Q: (function) proposal distribution
#OUTPUT
# Matrix
# 1st column is the neighbor
# 2nd column is the ptobability of moving from state X to the state in that row
alpha_p <- function(x,S,pi,Q){
  n_x <- Q(x,S)
  p_esc <- matrix(NA,ncol=2)
  for(y in n_x[,1]){
    n_y <- Q(y,S)
    index_x <- which(n_y[,1]==x)
    index_y <- which(n_x[,1]==y)
    prob <- n_x[index_y,2]*min(1,(pi[y]*n_y[index_x,2])/(pi[x]*n_x[index_y,2]))
    p_esc <- rbind(p_esc,c(y,prob))
  }
  p_esc <- p_esc[-1,]
  return(p_esc) 
}