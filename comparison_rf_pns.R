rm(list=ls())
set.seed(123)
S <- 1:10 #State space
#bimodal target distribution with 4 states between each modes
pi <- c(1,1,6,rep(1,4),6,1,1) #bimodal target distribution
true_pi <- pi/sum(pi)

#source('MH_function.r')
source('RF_function.r')
source('PNS_unbiased_function.r')

#Parameters
B <- 50 #Burn-in
M <- c(500,1000,2000,5000,10000,20000) #Sample sizes
N <- 200 #number of times to repeat the experiment

##### Simulation Rejection Free #####

#Distribution functions to consider
source('functions/Q_unif.R')
Q1 <- function(i,S){Q_unif(i,S,1)} #Uniform considering only 1 adjacent neighbor on each side
Q2 <- function(i,S){Q_unif(i,S,2)} #Uniform considering x adjacent neighbors on each side
Q3 <- function(i,S){Q_unif(i,S,3)}
Q4 <- function(i,S){Q_unif(i,S,4)}
dist_q <- list(Q_unif,Q1,Q2,Q3,Q4)

time_mcmc <- array(0,dim=c(N,length(M),length(dist_q))) #array to store execution time
tvd <- array(0,dim=c(N,length(M),length(dist_q))) #array to store values of TVD
no_sample <- c()
for(n in 1:N){ #number of simulations
  for(m in 1:length(M)){ #considering different sample sizes
    for(q in 1:length(dist_q)){#distinct proposal distributions
      start_time <- Sys.time()
      simulation <- mh_jump(S,3,B,M[m],pi,dist_q[[q]])
      time_mcmc[n,m,q] <- Sys.time() - start_time #record the time it took to run the simulation
      states_w_sample <- length(unique(simulation$sample))
      if(states_w_sample<length(S)){ #In case there was a state with no sample
        no_sample <- c(no_sample,paste0(length(S)-states_w_sample,' states with no sample, N:',n,' M:',M[m],' Q:',q ))
        simulation <- rbind(simulation,tibble(sample=S,mul=rep(1,length(S)))) #Add 1 to each state (to avoid states with 0 sample)
        est_prob <- simulation |> #estimate the probability
          group_by(sample) |> 
          summarize(fre = sum(mul)) |> 
          mutate(fre = fre-1) |> #Delete the additional sample added to avoid states with no 
          mutate(prob = fre/sum(fre)) |> 
          pull(prob) 
      }else{
        est_prob <- simulation |> #estimate the probability
          group_by(sample) |> 
          summarize(fre = sum(mul)) |> 
          mutate(prob = fre/sum(fre)) |> 
          pull(prob) 
      }
      tvd[n,m,q] <- 0.5*sum(abs(true_pi - est_prob))
    }
  }
  if(n%%10==0){print(paste("first",n,"simulations"))}
}

tvd_sum_mh <- apply(tvd,c(2,3),sum)/length(tvd[,1,1])
time_sum_mh <- apply(time_mcmc,c(2,3),sum)
result <- as.data.frame(cbind(tvd_sum_mh,time_sum_mh))
colnames(result) <- c(paste0('TVD Q',1:length(dist_q)),paste0('Time Q',1:length(dist_q)))
rownames(result) <- paste0("S=",M)
result <- rbind(result,c('Q1 unif on all neighbors', paste0('Q',2:5,' unif on ',1:4,' adjacent neighbors on each side'),rep(NA,length(dist_q))))
write.csv(result,paste0("simulation results/RF_sim_",Sys.Date(),".csv"))
if(length(no_sample)>0){write.csv(no_sample,paste0("simulation results/RF_no sample_",Sys.Date(),".csv"),row.names = F)}

##### Simulation Partial Neighbor Search #####

#PNSets to consider
#PNSets depend on the specific state space S 
#so this functions only work with |S|=10

PNS_def <- function(x,S,n1,n2){
  indx <- which(S==x)
  n1 <- unique(S[((indx-1)+n1)%%length(S)+1])
  n2 <- unique(S[((indx-1)+n2)%%length(S)+1])
  if(!all(S%in%unique(c(n1,n2)))){
    print("Neighborhoods not covering S")
    return(NA)
  }else{
    return(list(n1,n2))
  }
}
PNS_3 <- function(x){
  return(PNS_def(x,S,c(0,-2,2,5),c(0,1,3,4,-1,-3,-4)))
}
PNS_4 <- function(x){
  return(PNS_def(x,S,c(0,-2,-4,2,4),c(0,1,3,5,-1,-3)))
}
PNS_near <- function(x){
  return(PNS_def(x,S,c(0,-1,-2,1,2),c(0,3,4,5,-3,-4)))
}
PNS_5 <- function(x){
  indx <- which(S==x)
  n1 <- c(0,-2,2,5)
  n2 <- c(0,-1,1)
  n3 <- c(0,-3,-4,3,4)
  n1 <- unique(S[((indx-1)+n1)%%length(S)+1])
  n2 <- unique(S[((indx-1)+n2)%%length(S)+1])
  n3 <- unique(S[((indx-1)+n3)%%length(S)+1])
  n_list <- list(n1,n2,n3)
  if(!all(S%in%unique(c(n1,n2,n3)))){
    print("Neighborhoods not covering S")
    return(NA)
  }else{
    return(n_list)
  }
  
}

PNSets <- list(PNS_3,PNS_4,PNS_near,PNS_5)

#parameters
L <- 50
PNS_time_mcmc <- array(0,dim=c(N,length(M),length(PNSets))) #array to store execution time
PNS_tvd <- array(0,dim=c(N,length(M),length(PNSets))) #array to store values of TVD
no_sample <- c()
for(n in 1:N){ #number of simulations
  for(m in 1:length(M)){ #considering different sample sizes
    for(q in 1:length(PNSets)){#distinct proposal distributions
      start_time <- Sys.time()
      simulation <- PNS_unbiased(S,3,M[m],L,pi,Q_unif,PNSets[[q]])
      PNS_time_mcmc[n,m,q] <- Sys.time() - start_time #record the time it took to run the simulation
      states_w_sample <- length(unique(simulation$sample))
      if(states_w_sample<length(S)){ #In case there was a state with no sample
        no_sample <- c(no_sample,paste0(length(S)-states_w_sample,' states with no sample, N:',n,' M:',M[m],' PNS:',q ))
        simulation <- rbind(simulation,tibble(sample=S,mul=rep(1,length(S)))) #Add 1 to each state (to avoid states with 0 sample)
        est_prob <- simulation |> #estimate the probability
          group_by(sample) |> 
          summarize(fre = sum(mul)) |> 
          mutate(fre = fre-1) |> #Delete the additional sample added to avoid states with no 
          mutate(prob = fre/sum(fre)) |> 
          pull(prob) 
      }else{
        est_prob <- simulation |> #estimate the probability
          group_by(sample) |> 
          summarize(fre = sum(mul)) |> 
          mutate(prob = fre/sum(fre)) |> 
          pull(prob) 
      }
      PNS_tvd[n,m,q] <- 0.5*sum(abs(true_pi - est_prob))
    }
  }
  if(n%%10==0){print(paste("first",n,"simulations"))}
}

tvd_sum_mh <- apply(PNS_tvd,c(2,3),sum)/length(PNS_tvd[,1,1])
time_sum_mh <- apply(PNS_time_mcmc,c(2,3),sum)
result <- as.data.frame(cbind(tvd_sum_mh,time_sum_mh))
colnames(result) <- c(paste0('TVD PNS',1:length(PNSets)),paste0('Time PNS',1:length(PNSets)))
rownames(result) <- paste0("S=",M)

write.csv(result,paste0("simulation results/PNS_sim_",Sys.Date(),".csv"))
if(length(no_sample)>0){write.csv(no_sample,paste0("simulation results/PNS_no sample_",Sys.Date(),".csv"),row.names = F)}
