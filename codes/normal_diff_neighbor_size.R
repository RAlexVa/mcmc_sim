
if(grepl('codes',getwd())){setwd('..')}
source('functions/IIT_functions.R')
library(tidyverse)
semilla <- 1546
set.seed(semilla)

N <- 1 #Number of simulations
steps <- 10^5 #Max Number of steps each chain will move
states <- 1:201
S <- length(states) # States from 1 to 201
sigma2 <- 25^2
pi <- exp(-(states-101)^2/(2*sigma2)) #Normal centered at 101
#Later we shift x axis so it's centered at 0
#plot(pi/sum(pi))
est_pi <- numeric(S)
visit_state <- numeric(S)

### Balancing functions
{
h_max <- function(a){return(max(1,a))}
h_min <- function(a){return(min(1,a))}
h_sq <- function(a){return(sqrt(a))}
h_ratio <- function(a){return(a/(1+a))}
h_plus <- function(a){return(1+a)}
balancing_names <- c('max','min','sq','a/1+a','a+1')
balancing <- list(h_max,h_min,h_sq,h_ratio,h_plus)
}
#Define the size of the neighbor set based on the distance to the current state
d <- 16
### First algorithm
#Start at the mean, 101
Xnow <- 101
for(i in 1:steps){
  #Define the neighbors to consider based on distance
  neighbors <- c(Xnow+seq(-d,-1),Xnow+seq(1,d))
  #Exclude if exceeds proposed range
  neighbors <- neighbors[neighbors %in% states]
  temp <- IITupdate(pi=pi, nbr=neighbors,x=Xnow,h=h_max)
  
  est_pi[Xnow] <- est_pi[Xnow] + temp[[2]] #Update the weight for the current state
  Xnow <- temp[[1]] #Choose the next state
  visit_state[Xnow] <-visit_state[Xnow]+1 
}
### Second algorithm
{
# est_pi_v2 <- numeric(S)
# #Start at the mean, 101
# Xnow <- 101
# for(i in 1:steps){
#   #Define the neighbors to consider based on distance
#   neighbors <- c(Xnow+seq(-d,-1),Xnow+seq(1,d)) #Neighbors excluding current state
#   #Identify if exceeds proposed range
#   neighbors_check <- neighbors[neighbors %in% states]
#   neighbors <- neighbors%%S 
#   neighbors <- ifelse(neighbors==0,length(pi),neighbors) #Replace 0 mod S with the last element of S
# 
#   outrange <- TRUE
#   while(outrange){
#     temp <- IITupdate(pi=pi, nbr=neighbors,x=Xnow,h=h_sq)
#     if(temp[[1]] %in% neighbors_check){
#       outrange <- FALSE
#     }
#   }
#   est_pi_v2[Xnow] <- est_pi[Xnow] + temp[[2]] #Update the weight for the current state
#   Xnow <- temp[[1]] #Choose the next state
# }
}
### Plots
#,y2=est_pi_v2/sum(est_pi_v2)
tibble(x=-100:100, y=est_pi/sum(est_pi),yvis=visit_state/sum(visit_state)) |> 
ggplot(aes(x=x,y=y))+
geom_line()+
geom_line(aes(y=yvis),col='blue')+
labs(y='probability', title=paste('Distance to neighbor = ',d))+
scale_x_continuous(breaks=seq(-100,100,by=25))+
scale_y_continuous(breaks=seq(0,max(y),by=0.002))+
theme_minimal()






