library(dplyr)
if(grepl('codes',getwd())){setwd('..')}
source('functions/IIT_functions.R')
semilla <- 169
set.seed(semilla)
N <- 10000 #Number of simulations
steps <- 2*10^5 #Max Number of steps each chain will move 
S <- 200 #Length of state space
mode <- 200 #Max for pi
min_pit <- 1 #Min for pi
reg <- 50 # a mid level between max and min for pi
pi <- matrix(NA,nrow=4,ncol = S)

#4 alternatives to pi
# 1 is basically unimodal with a low probability pit surrounding the mode
pi[1,] <- c(rep(reg,10),rep(min_pit,90),mode,rep(min_pit,89),rep(reg,10))
plot(pi[1,], xlab='state')
# 2 is multimodal with low probabiliti areas among modes
pi[2,] <- c(rep(reg,25),rep(min_pit,50),mode*2,rep(reg,50),mode,rep(min_pit,50),rep(reg,22),mode)
plot(pi[2,], xlab='state')
# 3 similar to 2 but with 1 global optima
pi[3,] <- c(rep(reg,19),mode*2,rep(min_pit,19),mode,rep(reg,19),mode*10,rep(min_pit,19),mode*5,rep(reg,19),mode)
plot(pi[3,], xlab='state')
# 4 a random graph with 1 global mode, 6 local modes
pi[4,] <- sample(c(min_pit,reg,mode),prob = c(.6,.35,.05), replace = T,size=S)
pi[4,sample(1:S,1)] <- mode*10
plot(pi[4,], xlab='state')
plot(pi[4,-which(pi[4,]==max(pi[4,]))], xlab='state',ylab='pi w/o global mode')

# For S=1000
# # 3 similar to 2 but with 1 global optima
# pi[3,] <- c(rep(reg,99),mode*2,rep(min_pit,99),mode,rep(reg,99),mode*10,rep(min_pit,99),mode*5,rep(reg,99),mode)
# plot(pi[3,], xlab='state')
# # 4 a random graph with 1 global mode, 6 local modes
# pi[4,] <- sample(c(min_pit,reg,mode),prob = c(.6,.35,.05), replace = T,size=S)
# pi[4,sample(1:S,1)] <- mode*10
# plot(pi[4,], xlab='state')
# plot(pi[4,-which(pi[4,]==max(pi[4,]))], xlab='state',ylab='pi w/o global mode')
# 


######### Copied from sim_balancing.R
### Balancing functions
h_max <- function(a){return(max(1,a))}
h_min <- function(a){return(min(1,a))}
h_sq <- function(a){return(sqrt(a))}
h_ratio <- function(a){return(a/(1+a))}
h_plus <- function(a){return(1+a)}

balancing <- list(h_max,h_min,h_sq,h_ratio,h_plus)
#Measurements
# running time
# iteration to solution
# If they discover all the modes
# top 10 
# Starting point

run_time_m <- matrix(nrow=N,ncol=nrow(pi)*length(balancing)+1)
its <- matrix(nrow=N,ncol=nrow(pi)*length(balancing)+1)


for(i in 1:N){ #Loop for simulations
  print(paste0('Simulation #',i))
  for(p in 1:nrow(pi)){ #Loop for target distributions
    for(h in 1:length(balancing)){# Loop for balancing functions
      if(i%%S==0){x=S}else{x <- i-S*floor(i/S)} #deterministic starting point 
      start <- Sys.time()
      simulation <- IIT_opt(initial=x,steps,pi[p,],d=1,h=balancing[[h]])
      run_time <- Sys.time() - start
      
      optima <- which(pi[p,]==max(pi[p,]))
      n_opt <- length(optima)
      if(n_opt==1){ #First time the solution was reached
        visits <- which(simulation[[1]]==which(pi[p,]==max(pi[p,])))
        if(length(visits)==0){iter_sol <- -1}else{iter_sol <- min(visits)}
        
      }else{
        print("There was more than one optimium in:")
        print(paste0('i=',i,'p=',p,'h=',h))
        iter_sol <- -2
        next;
        # indx_op <- which(pi==max(pi))
        # visits <- which(simulation[[1]]%in%indx_op)
        # summ <- simulation[[1]][visits]
        # opt_matrix <- tibble(estado=summ,iteration=visits)
        # res <- opt_matrix |> group_by(estado) |> summarise(iter=min(iteration))
      }
      
      #Storing information
      
      run_time_m[i,1] <- x
      run_time_m[i,(1+(p-1)*length(balancing)+h)] <- run_time
      its[i,1] <- x
      its[i,(1+(p-1)*length(balancing)+h)] <- iter_sol
    }
  }
}


write.csv(run_time_m,paste0('results/run_times_',semilla,'_',S,'.csv'),row.names = F)
write.csv(its,paste0('results/its_',semilla,'_',S,'.csv'),row.names = F)

summatrix <- matrix(nrow=7,ncol=2,
                  data=c('seed',semilla,
                         'Number simulations',N,
                         'max steps',steps,
                         'Length S',S,
                         'mode',mode,
                         'min_pit',min_pit,
                         'reg',reg), byrow=T)
write.table(summatrix, paste0('results/summary_',semilla,'.txt'), append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
