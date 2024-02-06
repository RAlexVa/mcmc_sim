#rm(list=ls())
semilla <- 1546
S <- 100
color_lines <- 'gray35'
chart_type <- 'hist'
#chart_type <- 'line'
##### Summarize findings #####
library(tidyverse)
library(Hmisc)
#its <- read.csv(paste0('results/its_',semilla,'.csv'))
its <- read.csv(paste('results/its',S,semilla,'.csv',sep='_'))
pi <- read.csv(paste('results/pi',S,semilla,'.csv',sep='_'))
#The CSV files contain:
# column 1 -> the starting point
# columns 2 onwards -> Information arranged by target distribution pi and balancing function
colnames(its) <- c('x0',paste0(rep(1:nrow(pi), each = length(balancing)),"-",balancing_names))

its <- its |> pivot_longer(-x0,names_to='comb', values_to='its') |> 
  separate(comb,into=c('pi','h'), sep='-') |> 
  mutate(hit=ifelse(its==-1,0,1))

### Identifying local and global modes
its$mode <- 0 #Initialize the column
for(i in 1:nrow(pi)){
  #Find local modes
  glo_mode <- which(pi[i,]==max(pi[i,]))
  loc_mode <- which(pi[i,]%nin%c(min_pit,reg,pi[i,glo_mode])) #integer(0) if we defined only 1 mode
  its <- its |> mutate(mode=ifelse(pi!=i,mode,ifelse(x0%in%glo_mode,2,ifelse(x0 %in% loc_mode,1,0)))) |> 
    mutate(mod_id=ifelse(mode==0,0,1))
}
its$mode_global <- ifelse(its$mode==2,its$hit,0)
its$mode_local <- ifelse(its$mode==1,its$hit,0)

#### Target distributions PI #### 

for(i in 1:nrow(pi)){
  jpeg(paste("plots/pi",S,semilla,i,".jpg",sep='_'), width = 700, height = 400) 
  plot(t(pi[i,]), xlab='state', main=paste('Target distribution Pi',i),ylab='pi')
  dev.off()
}
jpeg(paste("plots/pi",S,semilla,"4_0.jpg",sep='_'), width = 700, height = 400) 
plot(t(pi[4,-which(pi[4,]==max(pi[4,]))]), xlab='state',main='Pi 4 w/o global mode',ylab='pi') #mode at entry 58
dev.off()


#### Hit rate from the starting point #### 
for(i in 1:nrow(pi)){
  dummy <- its |> 
    group_by(x0, pi, h) |> 
    summarise(hit=mean(hit),mode_global=mean(mode_global),mode_local=mean(mode_local)) |> 
    filter(pi==as.character(i)) 

jpeg(paste("plots/hr",S,semilla,"pi",i,"1.jpg",sep='_'), width = 1305, height = 326)
if(chart_type=='hist'){
# Histogram ##################################
plot_sum <-  dummy |> filter(h %in%c('min','max','sq')) |>
    ggplot(aes(x=x0,
               y=hit))+
    geom_bar(stat = 'identity')+
    geom_bar(#data = dummy,
             aes(x=x0,
                 y=mode_global,fill='red'),stat = 'identity',show.legend = FALSE)+
    geom_bar(#data = dummy,
             aes(x=x0,
                 y=mode_local,fill='blue'),stat = 'identity',show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='Hit rate', x='Initial state', title=paste0('Hit rate by initial state, ','Pi ',i))+
    theme_minimal()
}  
if(chart_type=='line'){
# Line plot ##################################
plot_sum <- dummy|> filter(h %in%c('min','max','sq')) |> 
    ggplot(aes(x=x0,ymin=0,ymax=hit))+
    geom_linerange(color=color_lines)+
    geom_linerange(
      aes(x=x0,
          ymin=0,ymax=mode_global,color='red'),show.legend = FALSE)+
    geom_linerange(
      aes(x=x0,
          ymin=0,ymax=mode_local,color='blue'),show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='ItS', x='Initial state', title=paste0('ItS by initial state, ','Pi ',i))+
    theme_minimal()
}
  show(plot_sum)
dev.off()
### Second part ################################## 

jpeg(paste("plots/hr",S,semilla,"pi",i,"2.jpg",sep='_'), width = 1305, height = 326) 
if(chart_type=='hist'){
# Histogram ##################################
plot_sum <-  dummy |> filter(h %in%c('a/1+a','a+1')) |>
    ggplot(aes(x=x0,
               y=hit))+
    geom_bar(stat = 'identity')+
    geom_bar(#data = dummy,
             aes(x=x0,
                 y=mode_global,fill='red'),stat = 'identity',show.legend = FALSE)+
    geom_bar(#data = dummy,
             aes(x=x0,
                 y=mode_local,fill='blue'),stat = 'identity',show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='Hit rate', x='Initial state', title=paste0('Hit rate by initial state, ','Pi ',i))+
    theme_minimal()
}
if(chart_type=='line'){
# Line plot ##################################
plot_sum <- dummy|> filter(h %in%c('a/1+a','a+1')) |> 
  ggplot(aes(x=x0,ymin=0,ymax=hit))+
  geom_linerange(color=color_lines)+
  geom_linerange(
    aes(x=x0,
        ymin=0,ymax=mode_global,color='red'),show.legend = FALSE)+
  geom_linerange(
    aes(x=x0,
        ymin=0,ymax=mode_local,color='blue'),show.legend = FALSE)+
  facet_wrap(~h)+ #facet_wrap(~pi+h)+
  labs(y='ItS', x='Initial state', title=paste0('ItS by initial state, ','Pi ',i))+
  theme_minimal()
}
  show(plot_sum)
dev.off()
}

####  Iteration to Solution #### 
for(i in 1:nrow(pi)){
  dummy <- its |> 
    filter(pi==as.character(i)) |> #One target distribution
    # filter(h %in%c('min','max','sq')) |> #Specific balancing functions
    filter(hit==1) |> #Only instances where the global optimum was reached
    group_by(x0, pi, h, mode_global,mode_local) |> 
    summarise(its=mean(its), g_mode=max(mode_global)*its, l_mode=max(mode_local)*its) 
jpeg(paste("plots/its",S,semilla,"pi",i,"1.jpg",sep='_'), width = 1305, height = 326)
if(chart_type=='hist'){
# Histogram ################################## 
  plot_sum <-  dummy|> filter(h %in%c('min','max','sq')) |>
    ggplot(aes(x=x0,
               y=its))+
    geom_bar(stat = 'identity')+
    geom_bar(#data = dummy,
      aes(x=x0,
          y=g_mode,fill='red'),stat = 'identity',show.legend = FALSE)+
    geom_bar(#data = dummy,
      aes(x=x0,
          y=l_mode,fill='blue'),stat = 'identity',show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='ItS', x='Initial state', title=paste0('ItS by initial state, ','Pi ',i))+
    theme_minimal()
}
if(chart_type=='line'){
# Line plot ################################## 
plot_sum <- dummy|> filter(h %in%c('min','max','sq')) |> 
    ggplot(aes(x=x0,ymin=0,ymax=its))+
    geom_linerange(color=color_lines)+
    geom_linerange(
      aes(x=x0,
          ymin=0,ymax=g_mode,color='red'),show.legend = FALSE)+
    geom_linerange(
      aes(x=x0,
          ymin=0,ymax=l_mode,color='blue'),show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='ItS', x='Initial state', title=paste0('ItS by initial state, ','Pi ',i))+
    theme_minimal()
}
show(plot_sum)
dev.off()

### Second part ################################## 
  
jpeg(paste("plots/its",S,semilla,"pi",i,"2.jpg",sep='_'), width = 1305, height = 326) 
if(chart_type=='hist'){
# Histogram ################################## 
plot_sum <-  dummy|> filter(h %in% c('a/1+a','a+1')) |>
    ggplot(aes(x=x0,
               y=its))+
    geom_bar(stat = 'identity')+
    geom_bar(#data = dummy,
      aes(x=x0,
          y=g_mode,fill='red'),stat = 'identity',show.legend = FALSE)+
    geom_bar(#data = dummy,
      aes(x=x0,
          y=l_mode,fill='blue'),stat = 'identity',show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='ItS', x='Initial state', title=paste0('ItS rate by initial state, ','Pi ',i))+
    theme_minimal()
}
if(chart_type=='line'){ 
# Line plot ################################## 
plot_sum <- dummy|> filter(h %in%c('a/1+a','a+1')) |> 
    ggplot(aes(x=x0,ymin=0,ymax=its))+
    geom_linerange(color=color_lines)+
    geom_linerange(
      aes(x=x0,
          ymin=0,ymax=g_mode,color='red'),show.legend = FALSE)+
    geom_linerange(
      aes(x=x0,
          ymin=0,ymax=l_mode,color='blue'),show.legend = FALSE)+
    facet_wrap(~h)+ #facet_wrap(~pi+h)+
    labs(y='ItS', x='Initial state', title=paste0('ItS by initial state, ','Pi ',i))+
    theme_minimal()
}
show(plot_sum)
dev.off()
}

### General Summary ###
{
#ItS
jpeg(paste("plots/its",S,semilla,"all pi.jpg",sep='_'), width = 1305, height = 326) 
plot_sum <- its |> 
  group_by(pi, h) |> 
  summarise(its=mean(its), hit=mean(hit)) |> 
  ggplot(aes(x=h, y=its, fill=h))+
  geom_bar(stat = 'identity')+
  theme_minimal()+
  labs(title='ItS for each target distribution')+
  facet_wrap(~pi)
show(plot_sum)
dev.off()

#Hit rate
jpeg(paste("plots/hr",S,semilla,"all pi.jpg",sep='_'), width = 1305, height = 326) 
plot_sum <-its |> 
  group_by(pi, h) |> 
  summarise(its=mean(its), hit=mean(hit)) |> 
  ggplot(aes(x=h, y=hit, fill=h))+ylim(0,1)+
  geom_text(aes(label=hit), vjust=-0.5)+
  geom_bar(stat = 'identity')+
  theme_minimal()+
  labs(title='Hit Rate for each target distribution')+
  facet_wrap(~pi)+
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5))
show(plot_sum)
dev.off()
}
  