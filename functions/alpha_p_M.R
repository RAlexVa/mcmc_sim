
##### AID FUNCTION #####
#Function to calculate the probability of escaping a state
#INPUT
# pi: target distribution
# Q: proposal distribution (proposl transition matrix)
#OUTPUT
# list
# alpha is the probability of escaping each state
# matrix is the new transition probabilities (considering that the chain always escapes from the current state)
alpha_p_M <- function(pi, Q){
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
