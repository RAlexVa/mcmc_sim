rm(list = ls())
library(tidyverse)
###### Specifications of the model #####
S <- 1:5 #State space
pi <- exp(S/2) #Proportions of target distribution
true_pi <- pi/sum(pi) #actual target distribution, just adding normalization
# Markov chain on a complete graph with 5 nodes

set.seed(123)
N <- 500 # # of simluations to repeat the experiment
B <- 51 #burn-in. Need to add 1 due to make the indexes match (initial state is at position 1, not position 0)
size <- c(500,1000,2000,5000,10000,20000) #Sample sizes to compare TVD

#Using 3 different proposal distributions
#1 Q(i,j) ∝ 1     uniform
#2 Q(i,j) ∝ i/j   more likely to move to low numbers
#3 Q(i,j) ∝ j/i   more likely to move to high numbers

#Building matrix with proposal distributions
#Q1 Q(i,j) ∝ 1
Q1 <- matrix(1/length(S),ncol=length(S),nrow=length(S))
diag(Q1) <- 0
#2 Q(i,j) ∝ i/j
Q2 <- matrix(0,ncol=length(S),nrow=length(S))
for(i in 1:length(S)){
  Q2[i,] <- (i/S)/sum(i/S[-i])}
diag(Q2) <- 0
#3 Q(i,j) ∝ j/i
Q3 <- matrix(0,ncol=length(S),nrow=length(S))
for(i in 1:length(S)){
  Q3[i,] <- (S/i)/sum(S[-i]/i)}
diag(Q3) <- 0

#Function to simulate a simple MCMC
mcmc_sim <- function(S,initial,B,sample_size,pi,Q){
  state <- c() #Initialize vector to store the sample
  state[1] <- initial #initial state
  
  for(i in 1:(B-1)){ #burn-in steps
    new_s <- sample(S, 1,prob=Q[state[i],]) #from the 4 neighbors sample using Q
    if (runif(1)< (pi[new_s]*Q[new_s,state[i]])/(pi[state[i]]*Q[state[i],new_s])) { #Check probability to accept or reject
      state[i+1] <- new_s
    }else{state[i+1] <- state[i]}
  } #End burn-in steps
  
  for(i in B:(B+sample_size-1)){
    new_s <- sample(S, 1,prob=Q[state[i],]) #from the 4 neighbors sample using Q
    if (runif(1)< (pi[new_s]*Q[new_s,state[i]])/(pi[state[i]]*Q[state[i],new_s])) { #Check probability to accept or reject
      state[i+1] <- new_s
    }else{state[i+1] <- state[i]}
  }
  return(state[-(1:B)])
}
#Function to calculate the probability of escaping a state
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
#Function to simulate a rejection free MCMC
mcmc_jump <- function(S,initial,B,sample_size,pi,Q){
  state <- c() #Initialize vector to store the sample
  multi <- c() #Initialize multiplicity chain
  state[1] <- initial #initial state
  esc_p <- alpha_p(pi,Q) #Obtain the escape probability and the transition prob.
  i <- 1 #initialize the entries for the state and multi vectors
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
  while(sum(multi)<(B+sample_size)){
    s_rep <- 1 + rgeom(1,esc_p$alpha[state[i]]) #Get the multiplicity for the current state
    if((sum(multi) + s_rep)>=(B+sample_size)){ #If the chain would remain more than the burn-in
      multi[i] <- (B+sample_size) - sum(multi) #Stay in that state the remaining number of steps
      break #End the sampling process
    }else{
      multi[i] <- s_rep #assign the value to the multiplicity list
      state[i+1] <- sample(S,1,prob=esc_p$matrix[state[i],]) #Get the next state
    }
    i <- i+1
  }
  return(tibble(mul=multi[-(1:burnin)], sample=state[-(1:burnin)])) #Return ignoring burnin
}


parallel_temp <- function(S, initial, B, sample_size, pi, Q, temps){
  state <- matrix(NA,nrow=B+sample_size,ncol=length(temps))
  state[1,] <- rep(initial,length(temps))
  
  for(i in 1:(B+sample_size-1)){#Run on the rows of the matrix
    for(t in 1:length(temps)){ #Run on the columns (for each temperature)
      new_s <- sample(S, 1,prob=Q[state[i,t],]) #propose a new state
      tau <- temps[t]
      if (runif(1)< (pi[new_s]^(1/tau)*Q[new_s,state[i,t]])/(pi[state[i,t]]^(1/tau)*Q[state[i,t],new_s])) { 
        state[i+1,t] <- new_s #Accept
      }else{state[i+1,t] <- state[i,t]} #Reject
    }
    #once the update happened for each temperature
    #We propose a swap in temperatures
    new_t <- sample(1:length(temps),2) #Choose 2 temperatures to swap
    t1 <- new_t[1]
    t2 <- new_t[2]
    
    prob_t <- (pi[state[i+1,t2]]^(1/temps[t1]) * pi[state[i+1,t1]]^(1/temps[t2]))/(pi[state[i+1,t1]]^(1/temps[t1]) * pi[state[i+1,t2]]^(1/temps[t2]))
  
    if(runif(1)< prob_t){
      #accept the swap
      temp <- state[i+1,t1]
      state[i+1,t1] <- state[i+1,t2]
      state[i+1,t2] <- temp
    }
    }
    return(state[-(1:B),])
}

rev <- parallel_temp(S, initial, B, sample_size, pi, Q2, temps)
#test
#sample <- mcmc_jump(S,1,B,size[k],pi,distributions[[d]])
# ex1 <- mcmc_sim(S,1,B,size[5],pi,Q1)
# ex2 <- mcmc_sim(S,1,B,size[5],pi,Q1)
# ex3 <- mcmc_sim(S,1,B,size[5],pi,Q1)
# 
# table(ex1)/length(ex1)
# table(ex2)/length(ex2)
# table(ex3)/length(ex3)
# true_pi
# 
# sum(abs(true_pi - table(ex3)/length(ex3)))
##### Measuring performance of normal MH ##### 
time_mcmc <- array(0,dim=c(N,length(size),3))
tvd <- array(0,dim=c(N,length(size),3))
distributions <- list(Q1,Q2,Q3)

#Measuriing normal MH
for(r in 1:N){# Times to repeat each experiment to get an average of TVD and time
 for(k in 1:length(size)){
   for(d in 1:length(distributions)){
     start_time <- Sys.time()
     sample <- mcmc_sim(S,1,B,size[k],pi,distributions[[d]]) #Run the simulation
     time_mcmc[r,k,d] <- Sys.time() - start_time #record the time it took to run the simulation
     est_prob <- table(sample)/size[k] #estimate the probability
     tvd[r,k,d] <- sum(abs(true_pi - est_prob)) #Record the estimated TVD 
   }
 } 
}

tvd_sum_mh <- apply(tvd,c(2,3),sum)
time_sum_mh <- apply(time_mcmc,c(2,3),sum)

##### Measuring performance of rejection free MH #####
time_mcmc <- array(0,dim=c(N,length(size),3))
tvd <- array(0,dim=c(N,length(size),3))
distributions <- list(Q1,Q2,Q3)


for(r in 1:N){# Times to repeat each experiment to get an average of TVD and time
  for(k in 1:length(size)){
    for(d in 1:length(distributions)){
      start_time <- Sys.time()
      sample <- mcmc_jump(S,1,B,size[k],pi,distributions[[d]]) #Run the simulation
      time_mcmc[r,k,d] <- Sys.time() - start_time #record the time it took to run the simulation
      est_prob <- sample |> #estimate the probability
                  group_by(sample) |> 
                  summarize(fre = sum(mul)) |> 
                  mutate(prob = fre/sum(fre)) |> 
                  pull(prob) 
      tvd[r,k,d] <- sum(abs(true_pi - est_prob)) #Record the estimated TVD 
    }
  } 
}

tvd_sum_rf <- apply(tvd,c(2,3),sum)
time_sum_rf <- apply(time_mcmc,c(2,3),sum)
