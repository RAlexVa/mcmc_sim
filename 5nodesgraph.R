library(tidyverse)
library(visNetwork)
library(networkD3)

nodes <- data.frame(id = 1:5, 
                    label = c("  A  ", "  B  ", "  C  ", "  D  ", "  E  "),                                 # add labels on nodes
                    x = c(0,50,20,-20,-50),
                    y = c(30,0,-30,-30,0),
                    group = c("pi(A) = 1", "pi(B) = 1", "pi(C) = e^3", "pi(D) = 1", "pi(E) = e^5"),            # add groups on nodes 
                    #value = rep(10^1000, 3),                                   # size adding value
                    shape = c("circle", "circle", "circle", "circle", "circle"),                  # control shape of nodes
                    color = c("moccasin", "lightsteelblue", "salmon", "darkseagreen", "darkorchid"),# color
                    shadow = rep(FALSE, 5))                  # shadow

# edges <- data.frame(from = c(1,2,3,4,5), 
#                     to = c(2,3,4,5,1),
#                     label = c(rep(" ", 5)),
#                     dashes = c(rep(FALSE, 5)), 
#                     color = c(rep("black", 5)),
#                     shadow = rep(FALSE, 5))                    

edges <- data.frame(from = combn(5,2)[1,], 
                    to = combn(5,2)[2,])
edges$label <- " "
edges$dashes <- FALSE
edges$color <- "black"
edges$shadow <- FALSE

  visNetwork(nodes, edges, width = "100%") %>%
   #visNodes(fixed = TRUE) %>%
   visGroups(groupname = "pi(A) = 1", color = "moccasin", shape = "circle", value = 1) %>%
   visGroups(groupname = "pi(B) = 1", color = "lightsteelblue", shape = "circle", value = 1) %>%
   visGroups(groupname = "pi(C) = e^3", color = "salmon", shape = "circle", value = 1) %>%
   visGroups(groupname = "pi(D) = 1", color = "darkseagreen", shape = "circle", value = 1) %>%
   visGroups(groupname = "pi(E) = e^5", color = "darkorchid", shape = "circle", value = 1) %>%
   visLegend(useGroups = TRUE, position = "right", main = "Target Density Values") 
  