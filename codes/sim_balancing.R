library(dplyr)
if(grepl('codes',getwd())){setwd('..')}
source('functions/IIT_functions.R')
semilla <- 1546
set.seed(semilla)
N <- 2000 #Number of simulations
steps <- 10^5 #Max Number of steps each chain will move 
S <- 100 #Length of state space
mode <- 100 #Max for pi
min_pit <- 1 #Min for pi
reg <- 50 # a mid level between max and min for pi
pi <- matrix(NA,nrow=4,ncol = S)
#4 alternatives to pi
# 1 is basically unimodal with a low probability pit surrounding the mode
pi[1,] <- c(rep(reg,29),rep(min_pit,20),mode,rep(min_pit,21),rep(reg,29))
# 2 is multimodal with low probabiliti areas among modes
pi[2,] <- c(rep(reg,20),rep(min_pit,19),mode*2,rep(reg,19),mode,rep(min_pit,20),rep(reg,19),mode)
# 3 similar to 2 but with 1 global optima
pi[3,] <- c(rep(reg,19),mode*2,rep(min_pit,19),mode,rep(reg,19),mode*10,rep(min_pit,19),mode*5,rep(reg,19),mode)
# 4 a random graph with 1 global mode, 6 local modes
pi[4,] <- sample(c(min_pit,reg,mode),prob = c(.6,.35,.05), replace = T,size=S)
pi[4,sample(1:S,1)] <- mode*10
write.csv(pi,paste('results/pi',S,semilla,'.csv',sep='_'),row.names = F)

### Balancing functions
h_max <- function(a){return(max(1,a))}
h_min <- function(a){return(min(1,a))}
h_sq <- function(a){return(sqrt(a))}
h_ratio <- function(a){return(a/(1+a))}
h_plus <- function(a){return(1+a)}
balancing_names <- c('max','min','sq','a/1+a','a+1')
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

##### Export results #####
write.csv(run_time_m,paste('results/run_times',S,semilla,'.csv',sep='_'),row.names = F)
write.csv(its,paste('results/its',S,semilla,'.csv',sep='_'),row.names = F)
summatrix <- matrix(nrow=7,ncol=2,
                    data=c('seed',semilla,
                           'Number simulations',N,
                           'max steps',steps,
                           'Length S',S,
                           'mode',mode,
                           'min_pit',min_pit,
                           'reg',reg), byrow=T)
write.table(summatrix, paste('results/summary_',S,semilla,'.txt',sep='_'), append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
