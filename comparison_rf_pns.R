rm(list=ls())
set.seed(123)
S <- 1:10 #State space
#bimodal target distribution with 4 states between each modes
pi <- c(1,1,6,rep(1,4),6,1,1) #bimodal target distribution
true_pi <- pi/sum(pi)

source('MH_function.r')
source('RF_function.r')
source('PNS_unbiased_function.r')

#Parameters
B <- 50 #Burn-in
M <- c(500,1000,2000,5000,10000,20000) #Sample sizes
N <- 200 #number of times to repeat the experiment

##### Simulation Rejection Free #####

#Distribution functions to consider
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
          neighbors <- unique(neighbors[neighbors != i]) #In case the # of adjacent neighbors is too big
        }
        #uniform proposal distribution
        n_size <- length(neighbors)
        prob <- rep(1/n_size,n_size)
        return(cbind(neighbors,prob))
      }
    }
  }
} #uniform on all neighbors
Q1 <- function(i,S){Q_unif(i,S,1)} #Uniform considering only 1 adjacent neighbor
Q2 <- function(i,S){Q_unif(i,S,2)} #Uniform considering x adjacent neighbors
Q3 <- function(i,S){Q_unif(i,S,3)}
Q4 <- function(i,S){Q_unif(i,S,4)}
dist_q <- list(Q_unif,Q1,Q2,Q3,Q4)
#dist_q <- list(Q_unif,Q2,Q4)

time_mcmc <- array(0,dim=c(N,length(M),length(dist_q))) #array to store execution time
tvd <- array(0,dim=c(N,length(M),length(dist_q))) #array to store values of TVD

for(n in 1:N){ #number of simulations
  for(m in 1:length(M)){ #considering different sample sizes
    for(q in 1:length(dist_q)){#distinct proposal distributions
      start_time <- Sys.time()
      simulation <- mh_jump(S,3,B,M[m],pi,dist_q[[q]])
      time_mcmc[n,m,q] <- Sys.time() - start_time #record the time it took to run the simulation
      est_prob <- simulation |> #estimate the probability
        group_by(sample) |> 
        summarize(fre = sum(mul)) |> 
        mutate(prob = fre/sum(fre)) |> 
        pull(prob) 
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

write.csv(result,paste0("simulation results/RF_sim_",Sys.Date(),".csv"))

##### Simulation Partial Neighbor Search #####

#PNSets to consider
PNS1 <- function(x){NBR_adj(x,S,nbr=2)} #Uniform considering only 1 adjacent neighbor
PNS2 <- function(x){NBR_adj(x,S,nbr=4)} #Uniform considering x adjacent neighbors
PNS3 <- function(x){NBR_adj(x,S,nbr=6)}
PNS4 <- function(x){NBR_adj(x,S,nbr=8)}
PNSets <- list(PNS1,PNS2,PNS3,PNS4)

#parameters
L <- 50
PNS_time_mcmc <- array(0,dim=c(N,length(M),length(PNSets))) #array to store execution time
PNS_tvd <- array(0,dim=c(N,length(M),length(PNSets))) #array to store values of TVD

for(n in 1:N){ #number of simulations
  for(m in 1:length(M)){ #considering different sample sizes
    for(q in 1:length(PNSets)){#distinct proposal distributions
      start_time <- Sys.time()
      simulation <- PNS_unbiased(S,3,M[m],L,pi,Q_unif,PNSets[[q]])
      PNS_time_mcmc[n,m,q] <- Sys.time() - start_time #record the time it took to run the simulation
      est_prob <- simulation |> #estimate the probability
        group_by(sample) |> 
        summarize(fre = sum(mul)) |> 
        mutate(prob = fre/sum(fre)) |> 
        pull(prob) 
      PNS_tvd[n,m,q] <- 0.5*sum(abs(true_pi - est_prob))
    }
  }
  if(n%%10==0){print(paste("first",n,"simulations"))}
}

tvd_sum_mh <- apply(PNS_tvd,c(2,3),sum)/length(PNS_tvd[,1,1])
time_sum_mh <- apply(PNS_time_mcmc,c(2,3),sum)
result <- as.data.frame(cbind(tvd_sum_mh,time_sum_mh))
colnames(result) <- c(paste0('TVD Q',1:length(PNSets)),paste0('Time Q',1:length(PNSets)))
rownames(result) <- paste0("S=",M)

write.csv(result,paste0("simulation results/PNS_sim_",Sys.Date(),".csv"))

