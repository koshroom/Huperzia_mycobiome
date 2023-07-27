#Nestedness analysis 
library(vegan)
library(bipartite) #this package contains the "nested" stat which looks at all possible ways to analyze nestedness
library(ggplot2)
library(reshape)
#install.packages('ggrastr')
#library(gridGraphics)
#library(readr)

rm(list = ls())
setwd("D:/personal/TXY/huperzia/1205koko/nestedness/hs")
setwd("D:/personal/TXY/huperzia/1205koko/nestedness/lc")
# loading data
# The rows (x) corresponds to a particular tissue
# The columns (y) correspond to a particular ASV
A <- as.matrix(read.csv("hs-asv-25.csv", header = TRUE, row.names = 1))
A <- as.matrix(read.csv("lc-asv-25.csv", header = TRUE, row.names = 1))
# obtain binary interaction matrix
# "1" indicates a species lives in a particular habitat, and '0' means the species was not found living in this particular habitat
A[A > 0] <- 1 

#Calculating nestedness using "nested"-------
nested(A, method = "binmatnest", rescale = FALSE, normalised = TRUE) #0 = cold = highly nested; 100 = hot = not nested at all)
nested(A, method = "NODF", rescale = FALSE, normalised = TRUE) #0 indicate non-nestedness,100 indicate perfect nesting. 
?nested

#Calculating nestedness using "nestedtemp" and plot -----
nest <- nestedtemp(A)
nestedness(A)
p <- plot(nest, col=rev(c("black", "white")), 
          names = TRUE, kind = "incidence", cex.axis = .75)

# extract this information to plot it with `ggplot
# get the presence-absence matrix
melted_nest = melt(nest$comm)
colnames(melted_nest) = c("Tissues", "ASVs", "Presence")
melted_nest$Presence = as.factor(melted_nest$Presence)
#  highlight potential HupA ASVs
write.csv(melted_nest,file = "melted_nest.csv")
melted_nest2 = read.csv("melted_nest.csv")
str(melted_nest)
str(melted_nest2)

melted_nest2$Tissues <- as.factor(melted_nest2$Tissues)
melted_nest2$ASVs <- as.factor(melted_nest$ASVs)
melted_nest2$Presence <- as.factor(melted_nest2$Presence)

# nest$smooth gives the expected distribution for nested communities
# extract it and put it into relation with the length of the matrix
smooth_x = nest$smooth$x * length(colnames(nest$comm))
smooth_y = length(rownames(nest$comm)) - 
  ((nest$smooth$y) * length(rownames(nest$comm))) + 0.5
smooth_data = as.data.frame(cbind(smooth_x, smooth_y))

g = ggplot(melted_nest2, aes(x = ASVs, y = Tissues, fill = Presence)) + 
  xlab(NULL) + 
  ylab(NULL) +
  geom_raster (show.legend = F) + # geom_raster 是条块大小相同时 geom_tile 的快速版本，而且当输出为 PDF 时所占空间也更小
  scale_y_discrete(limits = rev(rownames(nest$comm))) +
  geom_line(data = smooth_data, aes(x = smooth_x, y = smooth_y), color = "gray", size = 1,
            inherit.aes = F) + 
  scale_fill_manual(values = c("0" = "white", "1" = "black","2" = "red")) + 
  theme_minimal() + 
  theme(axis.text.y=element_text(size=12, face = "bold"), 
        axis.text.x = element_blank(), 
        panel.grid = element_blank())
g


# Nestedness Null Models ------

null_nest <- oecosimu(A, nestfun = "nestedtemp", 
                      method = "r00", #
                      nsimul = 1000, alternative = "less") #alternative：a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". Please note that the pp-value of two-sided test is approximately two times higher than in the corresponding one-sided test ("greater" or "less" depending on the sign of the difference).
null_nest
# Temperature of observed community is significantly less than (more nested) than null community simulation. 


# Comparing our nestedness model to null models ------

#To further see if our nestedness model is accurate, we will compare our results to null models
#Null models are constructed to test the hypothesis that the existence of a pattern (in this case, a nested community) is the result of random processes or can be expected by chance alone

#To do this, we will use the 'oecosimu' function, which takes a statistic or a vector of statistics in community and evaluates its significance in a series of simulated random communities.

#First, we will maintain the recorded number of presences in the matrix, but completely randomized the species – habitat distributions.
oecosimu(A, nestedtemp, "r00") #‘**’ 0.01 original data is significantly nested compared to null communities

#Second, we will maintain the species (rows) presences within the matrix but randomized the habitat distributions
oecosimu(A, nestedtemp, "r0")#‘*’ 0.05 

#Third, we will maintain the habitat (column) presences in the matrix but randomized the species distributions.
oecosimu(A, nestedtemp, "c0") #‘**’ 0.01


