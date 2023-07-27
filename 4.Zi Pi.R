# Clear the workspace and set the working directory
rm(list = ls())
setwd("D:/huperzia/Network/flashweave-25/")
# Load  packages
library(igraph)
library(ggplot2)

# Define the zi.pi() function.
# This code is referenced from:https://github.com/xyz1396/MOBGCB/blob/0aebf7440fd069223406aaf8b8fa169263bb9b21/script/network.Rmd#L196
zi.pi<-function(nodes_bulk, z.bulk, modularity_class, degree){
  z.bulk[abs(z.bulk)>0]<-1
  module<-which(colnames(nodes_bulk)==modularity_class)
  module.max<-max(nodes_bulk[,module])
  degree<-which(colnames(nodes_bulk)==degree)
  
  #Split the matrix based on modules
  bulk.module<-list(NA)
  length(bulk.module)<-module.max
  
  for(i in 1:max(nodes_bulk[,module])){
    bulk.module[[i]]<-z.bulk[which(nodes_bulk[,module]==i),which(nodes_bulk[,module]==i)]
    bulk.module[[i]]<-as.data.frame(bulk.module[[i]])
    rownames(bulk.module[[i]])<-rownames(z.bulk)[which(nodes_bulk[,module]==i)]
    colnames(bulk.module[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
  }
  
  #Calculate within-module degree z
  z_bulk<-list(NA)
  length(z_bulk)<-module.max
  
  for(i in 1:length(z_bulk)){
    z_bulk[[i]]<-bulk.module[[i]][,1]
    z_bulk[[i]]<-as.data.frame(z_bulk[[i]])
    colnames(z_bulk[[i]])<-"z"
    rownames(z_bulk[[i]])<-rownames(bulk.module[[i]])
  }
  
  #Calculate z value
  for(i in 1:max(nodes_bulk[,module])){
    if(length(bulk.module[[i]])==1){
      z_bulk[[i]][,1]<-0
    }else if(sum(bulk.module[[i]])==0){
      z_bulk[[i]][,1]<-0
    }else{
      k<-rowSums(bulk.module[[i]]) -1
      mean<-mean(k)
      sd<-sd(k)
      if (sd==0){
        z_bulk[[i]][,1]<-0
      }else{
        z_bulk[[i]][,1]<-(k-mean)/sd
      }
    }
  }
  
  #Merge z values
  for(i in 2:max(nodes_bulk[,module])) {
    z_bulk[[i]]<-rbind(z_bulk[[i-1]],z_bulk[[i]])
  }
  z_bulk<-z_bulk[[module.max]]
  
  #Split the matrix columns based on modules
  bulk.module1<-list(NA)
  length(bulk.module1)<-module.max
  
  for(i in 1:max(nodes_bulk[,module])){
    bulk.module1[[i]]<-z.bulk[,which(nodes_bulk[,module]==i)]
    bulk.module1[[i]]<-as.data.frame(bulk.module1[[i]])
    rownames(bulk.module1[[i]])<-rownames(z.bulk)
    colnames(bulk.module1[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
  }
  
  #Calculate among-module connectivity c
  c_bulk<-list(NA)
  length(c_bulk)<-module.max
  
  for(i in 1:length(c_bulk)){
    c_bulk[[i]]<-z.bulk[,1]
    c_bulk[[i]]<-as.matrix(c_bulk[[i]])
    colnames(c_bulk[[i]])<-"c"
    rownames(c_bulk[[i]])<-rownames(z.bulk)
    c_bulk[[i]][,1]<-NA
  }
  
  #Calculate each node's connection count squared for each module
  for(i in 1:max(nodes_bulk[,module])){
    c_bulk[[i]]<-rowSums(bulk.module1[[i]])
    c_bulk[[i]]<-as.matrix(c_bulk[[i]])
    c_bulk[[i]][which(nodes_bulk$modularity == i),] = c_bulk[[i]][which(nodes_bulk$modularity == i),] -1
    c_bulk[[i]]<-c_bulk[[i]]*c_bulk[[i]]
    colnames(c_bulk[[i]])<-"c"
    rownames(c_bulk[[i]])<-rownames(z.bulk)
  }
  
  #Sum the squared connection counts for each module
  for(i in 2:max(nodes_bulk[,module])){
    c_bulk[[i]]<-c_bulk[[i]]+c_bulk[[i-1]]
  }
  c_bulk<-c_bulk[[module.max]]
  
  # Calculate the final c values 
  for(i in 1: length(c_bulk)){
    if(nodes_bulk$degree[i]==1){
      c_bulk[i] <- 0
    }else{
      c_bulk[i] <- 1-(c_bulk[i]/(nodes_bulk$degree[i]*nodes_bulk$degree[i]))
    }
  }
  colnames(c_bulk)<-"c"
  
  # Integrate z and c values
  z_c_bulk<-c_bulk
  z_c_bulk<-as.data.frame(z_c_bulk)
  z_c_bulk$z<-z_bulk[match(rownames(c_bulk),rownames(z_bulk)),]
  z_c_bulk<-z_c_bulk[,c(2,1)]
  names(z_c_bulk)[1:2]<-c('within_module_connectivities','among_module_connectivities')
  
  # Add node IDs to the result
  z_c_bulk$nodes_id<-rownames(z_c_bulk)
  nodes_bulk$nodes_id<-rownames(nodes_bulk)
  z_c_bulk<-merge(z_c_bulk,nodes_bulk,by='nodes_id')
  z_c_bulk
  
}

# Create an igraph object from the edge list
edge_list <- read.csv("D:/huperzia/Network/flashweave-25/edgelilst_netA.csv", header = TRUE)
graph <- graph.data.frame(edge_list,directed=FALSE)
adjacency_unweight <- as.matrix(get.adjacency(graph, attr = NULL))

# Calculate node degree and module
V(graph)$degree <- degree(graph)
V(graph)$module <- membership(cluster_fast_greedy(graph)) 
node_list <- data.frame(
  nodes_id = V(graph)$name,
  degree = V(graph)$degree,
  modularity = V(graph)$module)

# Calculate module-within connectivity (Zi) and among-module connectivity (Pi)
rownames(node_list) <- node_list$nodes_id
zi_pi_result <- zi.pi(node_list, adjacency_unweight, degree ="degree", modularity_class = "modularity" )

# Remove NA values from the result
zi_pi_result <- na.omit(zi_pi_result)

# Categorize nodes into four types based on connectivity values
zi_pi_result[zi_pi_result$within_module_connectivities < 2.5 & zi_pi_result$among_module_connectivities < 0.62, 'type'] <- 'Peripherals'
zi_pi_result[zi_pi_result$within_module_connectivities < 2.5 & zi_pi_result$among_module_connectivities > 0.62, 'type'] <- 'Connectors'
zi_pi_result[zi_pi_result$within_module_connectivities > 2.5 & zi_pi_result$among_module_connectivities < 0.62, 'type'] <- 'Module hubs'
zi_pi_result[zi_pi_result$within_module_connectivities > 2.5 & zi_pi_result$among_module_connectivities > 0.62, 'type'] <- 'Network hubs'

# Data visualization
ggplot(zi_pi_result, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.8, size = 4, shape = 20) +
  scale_y_continuous(limits = c(-2, 5)) +
  scale_color_manual(values = c("#8491B4FF","#91D1C2FF","#F39B7FFF", "#4DBBD5FF"), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs')) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62, linetype = 2, size = 1) +
  geom_hline(yintercept = 2.5, linetype = 2, size = 1) +  
  theme_bw() +
  theme(axis.ticks.length = unit(-0.25, "cm"), 
        axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")), 
        axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) )

# Save the results to a CSV file
write.csv(zi_pi_result, "zi_pi_result.csv", row.names = FALSE)
