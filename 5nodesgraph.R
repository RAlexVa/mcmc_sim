library(tidyverse)
library(visNetwork)
library(networkD3)
library(RColorBrewer) 
# brew10 <- brewer.pal(10, "Set3")
# 
# col.dist <- function(inp, comp) sum( abs(inp - col2rgb(comp) ) )
# colors()[ apply(col2rgb(brew10), 2, 
#                 function(z) which.min( sapply(colors(), 
#                                               function(x) col.dist(inp=z, comp=x) ) ) ) ]
# 
# 

nodes <- data.frame(id = 1:3, 
                    label = c("  A  ", "  B  ", "  C  "),                                 # add labels on nodes
                    x = c(-50, 50, 0),
                    y = c(-30, -30, 30),
                    group = c("pi(A) = 1/6", "pi(B) = 2/6", "pi(C) = 3/6"),            # add groups on nodes 
                    #value = rep(10^1000, 3),                                   # size adding value
                    shape = c("circle", "circle", "circle"),                  # control shape of nodes
                    color = c("moccasin", "lightsteelblue", "salmon"),# color
                    shadow = rep(FALSE, 3))                  # shadow

nodes
edges <- data.frame(from = c(1, 2, 3), 
                    to = c(2, 3, 1),
                    label = c(rep(" ", 3)),
                    dashes = c(FALSE, FALSE, FALSE), 
                    color = c(rep("black", 3)),
                    shadow = rep(FALSE, 3))                    
edges

visNetwork(nodes, edges, width = "100%") %>%
   visNodes(fixed = TRUE) %>%
   visGroups(groupname = "pi(A) = 1/6", color = "moccasin", shape = "circle", value = 1) %>%
   visGroups(groupname = "pi(B) = 2/6", color = "lightsteelblue", shape = "circle", value = 1) %>%
   visGroups(groupname = "pi(C) = 3/6", color = "salmon", shape = "circle", value = 1) %>%
   visLegend(useGroups = TRUE, position = "right", main = "Target Density Values") 