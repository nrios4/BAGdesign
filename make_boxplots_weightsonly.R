rm(list = ls())

library(ggplot2)
library(latex2exp)

# setwd("C:/Users/nrios4/OneDrive - George Mason University - O365 Production/Papers in Progress/BAGDesign/BAGweightsonly")
setwd("~/GMU/Papers in Progress/BAGDesign/weights_only")


# all_weights = read.csv("BAGweights_n50sigma2binomialp2.csv",header = T)[,-1]
all_weights = read.csv("BAGweights_n50sigma2poissonp2.csv",header = T)[,-1]


# generate candidate design points
p = 2
xlist <- list()
for(i in 1:p){
  xlist[[i]] <- seq(from = -1, to = 1, by = 0.1)
  
}

candidates <- expand.grid(xlist)
xnames <- paste("x",1:p, sep = "")
colnames(candidates) <- xnames
candidates <- as.matrix(candidates)
N = nrow(candidates)

### Aggregation = Mean
BAG_weights_final = colMeans(all_weights[,1:(6*N)])
df_BAG <- data.frame(candidates, weights = BAG_weights_final, 
                     gamma = as.factor(c(rep(0,N), rep(0.01,N), rep(0.02,N),rep(0,N), rep(0.01,N), rep(0.02,N))),
                     aggregation = as.factor(c(rep("Mean Aggregation",3*N), rep("Geometric Median Aggregation",3*N))))
df_BAG_subset <- subset(df_BAG, df_BAG$weights > 0)
levels(df_BAG_subset$gamma) <- c('0' = TeX(r'($\gamma = 0$)'),
                                     '0.01' = TeX(r'($\gamma = 0.01$)'),
                                     '0.02' = TeX(r'($\gamma = 0.02$)'))
ggplot(df_BAG_subset, aes(x = x1, y = x2, fill = weights)) + 
geom_point(shape = 21, size = 1.5) + scale_fill_continuous(low = "lightblue", high = "darkblue") + 
  facet_grid(cols = vars(gamma), rows = vars(aggregation), labeller = labeller(gamma = label_parsed, aggregation = label_value)) +
  xlab(TeX('$x_1$')) + ylab(TeX('$x_2$')) 




