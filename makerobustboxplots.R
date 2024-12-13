# makerobustboxplots.R

library(latex2exp)

# setwd("C:/Users/riosn/OneDrive/Documents/GMU/Papers in Progress/BAGDesign")
# Change working directory to location of the simulation csv files

make_boxplot <- function(res1, mainstr, xlab = "", ylab = "", ylimmax = 8){
  
  if(ylimmax == 8){
    buffer1 = 1.3
  } else{
    # buffer1 = 0.3
    buffer1 = 0.5
  }
  
  res1 = res1[,-1]
  res1 = res1[complete.cases(res1),]
  for(i in 1:nrow(res1)){
    for(j in 1:ncol(res1)){
      if(res1[i,j] == "#NAME?"){
        res1[i,j] = NA
      }
    }
  }
  res2 = res1[is.finite(rowSums(res1)),1:8]
  
  colnames(res2) = c("BAG","BAG 0.02", "BAG 0.01", "BAG Median", "BAG Median 0.02", "BAG Median 0.01", "cluster","maximin")
  relative_dfs <- c()
  label_lengths <- numeric(8)
  for(i in 1:8){
    
    res0 <- res2[,i]
    res0[res0 == Inf]  <- 7.9
    res0 <- res0[res0 < 8]
    label_lengths[i] <- length(res0)
    relative_dfs <- c(relative_dfs, res0)
    
  }
  label_col = c(rep("BAG", label_lengths[1]), rep("BAG 0.02", label_lengths[2]), rep("BAG 0.01", label_lengths[3]),
                rep("BAG_med", label_lengths[4]), rep("BAG_med 0.02", label_lengths[5]), rep("BAG_med 0.01", label_lengths[6]),
                    rep("cluster", label_lengths[7]), rep("maximin", label_lengths[8]))
  relative_df = data.frame(deffs = relative_dfs, design = label_col)
  boxplot(deffs ~ design, data = relative_df, ylab = ylab, xlab = xlab, main = mainstr, xaxt = "n", ylim = c(0,ylimmax))
  axis(side = 1, labels = FALSE, at = 1:8)
  text(x = 1:8,y = par("usr")[3]-buffer1, labels=c("BAG", "BAG0.01","BAG0.02","BAGmed","BAGmed0.01","BAGmed0.02","cluster","maximin"),
        xpd = NA, srt = 50)
  abline(h = 1, col = "red", lwd = 2)
  
  # library(dplyr)
  # relative_df %>% group_by(design) %>% summarize(median = median(deffs))
  
}

# mu = mu1
# bin1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5binomialp2robust1.csv")
# bin2 = read.csv("RobustBAGsim_withmedian_n50sigma2binomialp2robust1.csv")
# bin3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5binomialp2robust1.csv")
# bin4 = read.csv("RobustBAGsim_withmedian_n200sigma2binomialp2robust1.csv")
# poi1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5poissonp2robust1.csv")
# poi2 = read.csv("RobustBAGsim_withmedian_n50sigma2poissonp2robust1.csv")
# poi3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5poissonp2robust1.csv")
# poi4 = read.csv("RobustBAGsim_withmedian_n200sigma2poissonp2robust1.csv")

# mu = mu2
# bin1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5binomialp2robust2.csv")
# bin2 = read.csv("RobustBAGsim_withmedian_n50sigma2binomialp2robust2.csv")
# bin3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5binomialp2robust2.csv")
# bin4 = read.csv("RobustBAGsim_withmedian_n200sigma2binomialp2robust2.csv")
# poi1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5poissonp2robust2.csv")
# poi2 = read.csv("RobustBAGsim_withmedian_n50sigma2poissonp2robust2.csv")
# poi3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5poissonp2robust2.csv")
# poi4 = read.csv("RobustBAGsim_withmedian_n200sigma2poissonp2robust2.csv")

# Uncomment these for Appendix B

# mu = mu1
# bin1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5binomialp2robust1sd_betagrid1.csv")
# bin2 = read.csv("RobustBAGsim_withmedian_n50sigma2binomialp2robust1sd_betagrid1.csv")
# bin3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5binomialp2robust1sd_betagrid1.csv")
# bin4 = read.csv("RobustBAGsim_withmedian_n200sigma2binomialp2robust1sd_betagrid1.csv")
# poi1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5poissonp2robust1sd_betagrid1.csv")
# poi2 = read.csv("RobustBAGsim_withmedian_n50sigma2poissonp2robust1sd_betagrid1.csv")
# poi3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5poissonp2robust1sd_betagrid1.csv")
# poi4 = read.csv("RobustBAGsim_withmedian_n200sigma2poissonp2robust1sd_betagrid1.csv")

# mu = mu2
bin1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5binomialp2robust2sd_betagrid1.csv")
bin2 = read.csv("RobustBAGsim_withmedian_n50sigma2binomialp2robust2sd_betagrid1.csv")
bin3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5binomialp2robust2sd_betagrid1.csv")
bin4 = read.csv("RobustBAGsim_withmedian_n200sigma2binomialp2robust2sd_betagrid1.csv")
poi1 = read.csv("RobustBAGsim_withmedian_n50sigma0.5poissonp2robust2sd_betagrid1.csv")
poi2 = read.csv("RobustBAGsim_withmedian_n50sigma2poissonp2robust2sd_betagrid1.csv")
poi3 = read.csv("RobustBAGsim_withmedian_n200sigma0.5poissonp2robust2sd_betagrid1.csv")
poi4 = read.csv("RobustBAGsim_withmedian_n200sigma2poissonp2robust2sd_betagrid1.csv")



par(mfrow = c(2,2), cex.axis = 0.6)

make_boxplot(bin1, mainstr = TeX(r"($\sigma = 0.5$)"), ylab = "n = 50", ylimmax = 2)
make_boxplot(bin2, mainstr = TeX(r"($\sigma = 2$)"), ylab = "n = 50", ylimmax = 2)
make_boxplot(bin3, mainstr = TeX(r"($\sigma = 0.5$)"), ylab = "n = 200", ylimmax = 2)
make_boxplot(bin4, mainstr = TeX(r"($\sigma = 2$)"), ylab = "n = 200", ylimmax = 2)


round(rbind(colMeans(bin1[,-(1:9)], na.rm = TRUE),colMeans(bin2[,-(1:9)], na.rm = TRUE),
      colMeans(bin3[,-(1:9)], na.rm = TRUE),colMeans(bin4[,-(1:9)], na.rm = TRUE)), 1)


par(mfrow = c(2,2), cex.axis = 0.6)

make_boxplot(poi1, mainstr = TeX(r"($\sigma = 0.5$)"), ylab = "n = 50", ylimmax = 2)
make_boxplot(poi2, mainstr = TeX(r"($\sigma = 2$)"), ylab = "n = 50", ylimmax = 2)
make_boxplot(poi3, mainstr = TeX(r"($\sigma = 0.5$)"), ylab = "n = 200", ylimmax = 2)
make_boxplot(poi4, mainstr = TeX(r"($\sigma = 2$)"), ylab = "n = 200", ylimmax = 2)



 
round(rbind(colMeans(poi1[,-(1:9)], na.rm = TRUE), colMeans(poi2[,-(1:9)], na.rm = TRUE),
      colMeans(poi3[,-(1:9)], na.rm = TRUE), colMeans(poi4[,-(1:9)], na.rm = TRUE)),1)

