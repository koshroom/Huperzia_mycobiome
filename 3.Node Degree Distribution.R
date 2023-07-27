# Load packages
rm(list = ls())
setwd("D:/huperzia/Network/flashweave-25/")
library(igraph)
library(ggplot2)
library(ggrepel)

# Load dataset and create an igraph object from the edge list
edge_list <- read.csv("D:/huperzia/Network/flashweave-25/edgelilst_netA.csv", header = TRUE)
graph <- graph.data.frame(edge_list)

# Calculate node degree and extract degree distribution statistics
V(graph)$degree <- degree(graph)
degree_dist <- table(V(graph)$degree)
degree_num <- as.numeric(names(degree_dist))
degree_count <- as.numeric(degree_dist)
data <- data.frame(degree = degree_num, count = degree_count)

# View node degree distribution
par(mfrow = c(1, 3))
hist(V(graph)$degree, xlab = 'Degree', ylab = 'Frequency', main = 'Degree Distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', main = 'Degree Distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', ylab = 'Log-count', main = 'Log-log Degree Distribution')

# Fit power law distribution and extract key values
fit <- nls(count ~ a*degree^b, data = data, start = list(a = 2, b = 1.5)) #manually specify the initial values of a and b, and try more if the data is different
summary(fit)
a <- round(coef(fit)[1], 3)
b <- round(coef(fit)[2], 3)

# Calculate R-squared value
fitted_values <- fitted(fit)
SSre <- sum((data$count - fitted_values)^2) #SSre：the sum of the squares of the distances of the points from the fit
SStot <- sum((data$count - mean(data$count))^2) #SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
R2 <- round(1 - SSre/SStot, 3)

# Calculate p-value using permutation test
p_num <- 1
dat_rand <- data
for (i in 1:999) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count - fitted_values)^2)
    SStot_rand <- sum((dat_rand$count - mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / 1000  # 999 + 1

# Data visualization with ggplot2
p <- ggplot(data, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')

# Add annotations for formula fitting, R-squared, and p-value
label <- data.frame(
    formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
    R2 = sprintf('italic(R^2) == %.3f', R2),
    p_value = sprintf('italic(P) < %.3f', p_value)
)

p + geom_text(x = 13, y = 210, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 13, y = 190, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 13, y = 170, aes(label = p_value), data = label, parse = TRUE, hjust = 0)
