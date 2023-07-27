# Clear the workspace and set the working directory
rm(list = ls())
setwd("D:/huperzia/figure/upset")
# Load packages
library(UpSetR)
library(tidyverse)

# Data Preparation, including otu abundance table and group information table
otu <- read.delim("otu.txt", header = T,  sep = '\t', stringsAsFactors = F, check.names = F)
group <- read.delim("group.txt", header = F, sep = '\t', stringsAsFactors = F)
names(group) <- c('sample', 'group')

# Transform abundance table into binary matrix
binary_matrix <- as.data.frame(apply(otu[, -1], 2, function(x) ifelse(x > 0, 1, 0)))
binary_matrix <- cbind(otu[, 1], binary_matrix)
colnames(binary_matrix)[1] = 'sample'

# Merge samples based on grouping information
merged_matrix <- binary_matrix %>%
  left_join(group, by = "sample") %>%
  group_by(group) %>%
  summarize(across(-sample, sum))

# Remove the Group column
merged_matrix <- t(merged_matrix)
colnames(merged_matrix)=merged_matrix[1,]
merged_matrix <- merged_matrix[-1,]
merged_matrix[merged_matrix>0]=1
write.table(merged_matrix, "merged_matrix.txt", sep = '\t', col.names = T,row.names = T, na = '', quote = F)
mm <- read.delim("merged_matrix.txt", header = T,  sep = '\t', stringsAsFactors = F, check.names = F)

# Generate the UpSetR plot
upset(mm, sets = c("Dc_leaf","Dc_stem","Ha_bulbil","Ha_normal_leaf","Ha_shoot","Ha_sporangium","Ha_stem","Ha_young_leaf"), 
      order.by = "freq", empty.intersections = "on") #sets.bar.color = "#56B4E9",