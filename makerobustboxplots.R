# Change directory to the location of csv files that are located in SimResults.zip
setwd("C:/Users/nrios4/OneDrive - George Mason University - O365 Production/Papers in Progress/BAGDesign")

make_boxplot <- function(res1, mainstr, xlab = "design", ylab = ""){
  res1 = res1[,-1]
  res1 = res1[complete.cases(res1),]
  res2 = res1[is.finite(rowSums(res1)),]
  colnames(res2) = c("local","BAG","cluster","maximin")
  relative_to_BAG = c(exp(res2[,1] - res2[,2]), exp(res2[,3] - res2[,2]), exp(res2[,4] - res2[,2]))
  label_col = c(rep("local", nrow(res2)), rep("cluster", nrow(res2)), rep("maximin", nrow(res2)))
  relative_df = data.frame(deffs = relative_to_BAG, design = label_col)
  boxplot(deffs ~ design, data = relative_df, ylab = ylab, xlab = xlab, main = mainstr)
  abline(h = 1, col = "red", lwd = 2)
  
  library(dplyr)
  relative_df %>% group_by(design) %>% summarize(mean = mean(deffs))
  
}


make_boxplot_v2 <- function(res1, mainstr, xlab = "design", ylab = ""){
  res1 = res1[,-1]
  res1 = res1[complete.cases(res1),]
  for(i in 1:nrow(res1)){
    for(j in 1:ncol(res1)){
      if(res1[i,j] == "#NAME?"){
        res1[i,j] = NA
      }
    }
  }
  res1[,1] = as.numeric(res1[,1])
  res2 = res1[is.finite(rowSums(res1)),]
  colnames(res2) = c("BAG","BAG 0.05", "BAG 0.025", "cluster","maximin")
  relative_to_local = c(exp(res2[,2] - res2[,1]), exp(res2[,3] - res2[,1]), exp(res2[,4] - res2[,1]),
                        exp(res2[,5] - res2[,1]), exp(res2[,6] - res2[,1]))
  label_col = c(rep("BAG", nrow(res2)), rep("BAG 0.05", nrow(res2)), rep("BAG 0.025", nrow(res2)), rep("cluster", nrow(res2)), rep("maximin", nrow(res2)))
  relative_df = data.frame(deffs = relative_to_local, design = label_col)
  boxplot(deffs ~ design, data = relative_df, ylab = ylab, xlab = xlab, main = mainstr)
  abline(h = 1, col = "red", lwd = 2)
  
  library(dplyr)
  relative_df %>% group_by(design) %>% summarize(mean = mean(deffs))
  
}

bin1 = read.csv("RobustBAGsim_n50sigma0.5binomialp2v4.csv")
bin2 = read.csv("RobustBAGsim_n50sigma2binomialp2v4.csv")
bin3 = read.csv("RobustBAGsim_n200sigma0.5binomialp2v4.csv")
bin4 = read.csv("RobustBAGsim_n200sigma2binomialp2v4.csv")
poi1 = read.csv("RobustBAGsim_n50sigma0.5poissonp2v4.csv")
poi2 = read.csv("RobustBAGsim_n50sigma2poissonp2v4.csv")
poi3 = read.csv("RobustBAGsim_n200sigma0.5poissonp2v4.csv")
poi4 = read.csv("RobustBAGsim_n200sigma2poissonp2v4.csv")


par(mfrow = c(2,2))

make_boxplot_v2(bin1, mainstr = "sigma = 0.5", ylab = "n = 50")
make_boxplot_v2(bin2, mainstr = "sigma = 2", ylab = "n = 50")
make_boxplot_v2(bin3, mainstr = "sigma = 0.5", ylab = "n = 200")
make_boxplot_v2(bin4, mainstr = "sigma = 2", ylab = "n = 200")

par(mfrow = c(2,2))

make_boxplot_v2(poi1, mainstr = "sigma = 0.5", ylab = "n = 50")
make_boxplot_v2(poi2, mainstr = "sigma = 2", ylab = "n = 50")
make_boxplot_v2(poi3, mainstr = "sigma = 0.5", ylab = "n = 200")
make_boxplot_v2(poi4, mainstr = "sigma = 2", ylab = "n = 200")


colMeans(bin1[,-(1:7)], na.rm = TRUE)
colMeans(bin2[,-(1:7)], na.rm = TRUE)
colMeans(bin3[,-(1:7)], na.rm = TRUE)
colMeans(bin4[,-(1:7)], na.rm = TRUE)

colMeans(poi1[,-(1:7)], na.rm = TRUE)
colMeans(poi2[,-(1:7)], na.rm = TRUE)
colMeans(poi3[,-(1:7)], na.rm = TRUE)
colMeans(poi4[,-(1:7)], na.rm = TRUE)

