setwd("C:/Users/Karsten/Documents/GitHub/dbaccess")
source("BinningLossFunctions.R")

#----------------------------------------------------
#### Function Tests
#----------------------------------------------------

## Try 1d binning functions
xs <- runif(100,0,100)
cs <- rexp(100,.1)

# standard binning 
standbinnedxs <- StandRectBin1d(xs, origin=0, width=10)
qplot(xs,  standbinnedxs)
# random binning
randbinnedxs <- RandRectBin1d(xs, origin=0, width=10)
qplot(xs,  randbinnedxs)
# quantile binning 
quantbinnedcs <- QuantBin1d(cs, 7)
qplot(cs,  quantbinnedcs)


### Try 2d rectangular binning with MVN data
n=50000
mu <- 50
sig <- 15
xnorm <- rnorm(n, mean=mu, sd=sig)
ynorm <- rnorm(n, mean=mu, sd=sig)
qplot(xnorm, ynorm)
# create binned data for standard and random 
binout1 <- RectBin2d(xnorm,ynorm,0,0,10,10,type="standard")
binout2 <- RectBin2d(xnorm,ynorm,0,0,10,10,type="random")
# take a peek at binned data format
head(binout1[[1]])
head(binout2[[2]])
# compare spatial losses
binout1[[2]]
binout2[[2]]
#compare visualizations
qplot(binxs, binys, geom="tile", fill=binfreq, data=binout1[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
qplot(binxs, binys, geom="tile", fill=binfreq, data=binout2[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
#compare with shade by log counts
qplot(binxs, binys, geom="tile", fill=log(binfreq+1), data=binout1[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
qplot(binxs, binys, geom="tile", fill=log(binfreq+1), data=binout2[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")



### Try out frequency binning on the same data using standard/quantile/log-count binning
binoutFreqGroups <- freqBin(binout1, binType="standard", ncolor=5, logCount=FALSE)
binoutQuantFreqGroups <- freqBin(binout1, binType="quantile", ncolor=5, logCount=FALSE)
binoutLogFreqGroups <- freqBin(binout1, binType="standard", ncolor=5, logCount=TRUE)
# Peek into how frequency bins are defined
head(binoutFreqGroups[[1]])
head(binoutQuantFreqGroups[[1]])
head(binoutLogFreqGroups[[1]])
# Compare Frequency Losses (Note I dont know if these are all directly comparable after logging counts)
binoutFreqGroups[[2]]
binoutQuantFreqGroups[[2]]
binoutLogFreqGroups[[2]] 
# Plot with Orginal Frequencies
qplot(binxs, binys, geom="tile", fill=binfreq, data=binoutFreqGroups[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
# Plot with Binned Frequencies
qplot(binxs, binys, geom="tile", fill=freqgroup, data=binoutFreqGroups[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend", breaks=unique(binoutFreqGroups[[1]]$freqgroup))
# Plot with Quantile Binned Frequencies
qplot(binxs, binys, geom="tile", fill=freqgroup, data=binoutQuantFreqGroups[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend", breaks=unique(binoutQuantFreqGroups[[1]]$freqgroup))
# Plot with LogCount Binned Frequencies
qplot(binxs, binys, geom="tile", fill=freqgroup, data=binoutLogFreqGroups[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend", breaks=unique(binoutLogFreqGroups[[1]]$freqgroup))
