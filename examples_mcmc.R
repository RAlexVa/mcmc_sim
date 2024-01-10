rm(list = ls())
library(tidyverse)
###### Specifications of the model #####
S <- 1:5 #State space
pi <- c(1,2,exp(3),4,exp(5)) #Proportions of target distribution
true_pi <- pi/sum(pi) #actual target distribution, just adding normalization
# Markov chain on a complete graph with 5 nodes

##### Proposal distributions ##### 
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


##### IMPORTING FUNCTIONS FOR SIMULATION #####
source('functions/MH_matrix_function.R')
source('functions/RF_matrix_function.R')
source('functions/paralleltemp_function.R')

##### Measuring performance ##### 
set.seed(123)
N <- 500 # # of simulations to repeat the experiment
B <- 30 #burn-in
size <- c(100,500,1000,5000,10000) #Sample sizes to compare TVD
distributions <- list(Q1,Q2,Q3) #List with the distributions to consider

##### Measuring performance of Metropolis Hastings #####
time_mcmc <- array(0,dim=c(N,length(size),3)) #array to store execution time
tvd <- array(0,dim=c(N,length(size),3)) #array to store values of TVD
print('Metropolis Hastings simulation')
#progress bar 
#https://r-coder.com/progress-bar-r/#:~:text=Insert%20progress%20bar,-The%20txtProgressBar%20function&text=The%20most%20common%20functions%20used,arguments%20that%20you%20can%20customize.
pb <- txtProgressBar(min = 0,max = N,style = 3,width = 50,char = "=") 
for(r in 1:N){# Times to repeat each experiment to get an average of TVD and time
 for(k in 1:length(size)){
   for(d in 1:length(distributions)){
     start_time <- Sys.time()
     sample <- mh_sim(S,1,B,size[k],pi,distributions[[d]]) #Run the simulation
     time_mcmc[r,k,d] <- Sys.time() - start_time #record the time it took to run the simulation
     est_prob <- (table(c(sample,S))-1)/size[k] #estimate the stationary distribution
     #Added 1 of each state and then substract it so we correctly get the 0's
     tvd[r,k,d] <- 0.5*sum(abs(true_pi - est_prob)) #Record the estimated TVD 
   }
 }
  setTxtProgressBar(pb, r)#update progress bar  
}
close(pb)
tvd_sum_mh <- apply(tvd,c(2,3),sum)/length(tvd[,1,1])
time_sum_mh <- apply(time_mcmc,c(2,3),sum)

result <- as.data.frame(cbind(tvd_sum_mh,time_sum_mh))
colnames(result) <- c('TVD Q1','TVD Q2','TVD Q3','Time Q1','Time Q2','Time Q3')
rownames(result) <- paste0("S=",size)
write.csv(result,paste0("MH_sim_",Sys.Date(),'.csv'))
##### Measuring performance of rejection free MH #####
time_mcmc <- array(0,dim=c(N,length(size),3))
tvd <- array(0,dim=c(N,length(size),3))
print('Rejection Free MH simulation')
pb <- txtProgressBar(min = 0,max = N,style = 3,width = 50,char = "=") #progress bar  
for(r in 1:N){# Times to repeat each experiment to get an average of TVD and time
  for(k in 1:length(size)){
    for(d in 1:length(distributions)){
      start_time <- Sys.time()
      sample <- mh_jump(S,1,B,size[k],pi,distributions[[d]]) #Run the simulation
      time_mcmc[r,k,d] <- Sys.time() - start_time #record the time it took to run the simulation
      est_prob <- sample |> #estimate the probability
                  group_by(sample) |> 
                  summarize(fre = sum(mul)) |> 
                  mutate(prob = fre/sum(fre)) |> 
                  pull(prob) 
      tvd[r,k,d] <- 0.5*sum(abs(true_pi - est_prob)) #Record the estimated TVD 
    }
  }
  setTxtProgressBar(pb, r)#update progress bar  
}
close(pb)
tvd_sum_rf <- apply(tvd,c(2,3),sum)/length(tvd[,1,1])
time_sum_rf <- apply(time_mcmc,c(2,3),sum)

result <- as.data.frame(cbind(tvd_sum_rf,time_sum_rf))
colnames(result) <- c('TVD Q1','TVD Q2','TVD Q3','Time Q1','Time Q2','Time Q3')
rownames(result) <- paste0("S=",size)
write.csv(result,paste0("RF_sim_",Sys.Date(),'.csv'))
##### Measuring performance of parallel tempering #####
time_mcmc <- array(0,dim=c(N,length(size),3))
tvd <- array(0,dim=c(N,length(size),3))
temps <- c(1,2,3.5,5)
print('Parallel tempering simulation')
pb <- txtProgressBar(min = 0,max = N,style = 3,width = 50,char = "=") #progress bar  
for(r in 1:N){# Times to repeat each experiment to get an average of TVD and time
  for(k in 1:length(size)){
    for(d in 1:length(distributions)){
      start_time <- Sys.time()
      sample <- parallel_temp(S,1,B,size[k],pi,distributions[[d]], temps) #Run the simulation
      time_mcmc[r,k,d] <- Sys.time() - start_time #record the time it took to run the simulation
      est_prob <- (table(c(sample[,1],S))-1)/size[k] #estimate the stationary distribution
      tvd[r,k,d] <- 0.5*sum(abs(true_pi - est_prob)) #Record the estimated TVD 
    }
  }
  setTxtProgressBar(pb, r)#update progress bar  
}
close(pb)
tvd_sum_pt <- apply(tvd,c(2,3),sum)/length(tvd[,1,1])
time_sum_pt <- apply(time_mcmc,c(2,3),sum)

result <- as.data.frame(cbind(tvd_sum_pt,time_sum_pt))
colnames(result) <- c('TVD Q1','TVD Q2','TVD Q3','Time Q1','Time Q2','Time Q3')
rownames(result) <- paste0("S=",size)
write.csv(result,paste0("PT_sim_",Sys.Date(),'.csv'))
