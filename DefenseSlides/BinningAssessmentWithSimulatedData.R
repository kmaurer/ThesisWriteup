### Preliminaries 
library(ggplot2)
library(mvtnorm)
library(plyr)
library(Hmisc)
library(gridExtra)


### Goal:
# use simulated bivariate data to compare binning algorithms under various scenarios

### Comparisons:
# Use spatial loss to compare standard and rectangular binning algorithms
# Use frequency loss to compare standard, quantile and log count frequency binning algorithms. 

### Scenarios:
# Simulate data from Bivariate Uniform, Normal and Exponential
# Round simulated data at different resolutions to mimic measurement
# Use a variety of bin sizes

#------------------------------------------------------------

### Try 2d binning with MVN data
n=100000
set.seed(31415)
unifdat <- data.frame(x= runif(n, min=0, max=100),
                      y= runif(n, min=0, max=100),
                      type="Uniform")
normdat <- data.frame(x= rnorm(n, mean=50, sd=11),
                     y= rnorm(n, mean=50, sd=11),
                     type="Normal")
expdat <- data.frame(x= rexp(n, .11),
                     y= rexp(n, .11),
                     type="Exponential")
datlist <- list(unifdat,normdat,expdat)
### round all data to resolution alpha
alpha <- 2
for(i in 1:3){
  datlist[[i]]$roundx <- ((datlist[[i]]$x + 1)%/%alpha)*alpha
  datlist[[i]]$roundy <- ((datlist[[i]]$y + 1)%/%alpha)*alpha
}


# plot the original data
# make rounded data find overlapping points to make plot rendering faster
alldat <- rbind(datlist[[1]],datlist[[2]],datlist[[3]])
smallround <- ddply(alldat,.(type,roundx,roundy),summarise,
      overlap = length(roundx))
head(smallround)
qplot(x,y,geom="point",data=alldat,facets=.~type)
qplot(roundx,roundy,geom="point",data=smallround,facets=.~type)



# create binned data for standard and random using variety of bin sizes (from 2 to 20, by 2)
#save spatial loss information
spatlossdat <- data.frame(NULL)
index=0
distindex = 0
distrib <- c("Uniform", "Normal", "Exponential")
for (dat in 1:3){
  distindex <- distindex + 1
  for (width in seq(2,20,by=2)){
      index <- index+1
      binout <- RectBin2d(datlist[[dat]]$x,datlist[[dat]]$y,0,0,width,width,type="standard")
      binout[[2]]$bintype <- "standard"
      binout[[2]]$distrib <- distrib[distindex]
      timeStart <- Sys.time()
      temp <- RectBin2dNoLoss(datlist[[dat]]$x,datlist[[dat]]$y,0,0,width,width,type="standard")
      binout[[2]]$binTime <- Sys.time() - timeStart
      spatlossdat <- rbind(spatlossdat, binout[[2]])
      index <- index+1
      binout <- RectBin2d(datlist[[dat]]$x,datlist[[dat]]$y,0,0,width,width,type="random")
      binout[[2]]$bintype <- "random"
      binout[[2]]$distrib <- distrib[distindex]
      timeStart <- Sys.time()
      temp <- RectBin2dNoLoss(datlist[[dat]]$x,datlist[[dat]]$y,0,0,width,width,type="random")
      binout[[2]]$binTime <- Sys.time() - timeStart
      spatlossdat <- rbind(spatlossdat, binout[[2]])
      print(c(distrib[distindex],width))
  }
}
# write.csv(spatlossdat,"SpatialLossData.csv", row.names=FALSE)
spatlossdatSquareBins <- spatlossdat
spatlossdatSquareBins <- read.csv("SpatialLossData.csv",header=TRUE)

qplot(widthx,propSpatialLoss,geom="path",colour=bintype, linetype=bintype,
      data=spatlossdatSquareBins,group=bintype, facets=.~distrib, size=I(1)) +
  theme_bw()

qplot(widthx,totalSpatialLoss,geom="path",colour=bintype, linetype=bintype,
      data=spatlossdatSquareBins,group=bintype, facets=.~distrib, size=I(1)) +
  theme_bw()

### round all data to resolution alpha
alpha <- 2
for(i in 1:3){
datlist[[i]]$roundx <- ((datlist[[i]]$x + 1)%/%alpha)*alpha
datlist[[i]]$roundy <- ((datlist[[i]]$y + 1)%/%alpha)*alpha
}
head(datlist[[1]])


## Creat visualization for why we should select binwidths that are scalar multiples of alpha 
# create binned data for standard and random with data resolution and binwidth mismatch
binout1 <- RectBin2d(datlist[[2]]$roundx,datlist[[2]]$roundy,0,0,1.5*alpha,1.5*alpha,type="standard")
binout2 <- RectBin2d(datlist[[2]]$roundx,datlist[[2]]$roundy,0,0,1.5*alpha,1.5*alpha,type="random")
#compare visualizations for standard and random

qplot(binxs, binys, geom="tile", fill=binfreq, data=binout1[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
qplot(binxs, binys, geom="tile", fill=binfreq, data=binout2[[1]]) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")





##  Check out the spatial loss for the data that is recorded to resolution = alpha at various origin offsets
# standard and random binning with bin sizes that are scalar multiples of the data resolution
# save spatial loss information
spatlossdatRound <- data.frame(NULL)
index=0
distindex = 0
distrib <- c("unif", "norm", "exp")
widths <- seq(alpha,5*alpha,by=alpha)
for (dat in 1:3){
  distindex <- distindex + 1
  for (width in widths){
    for (origin in seq(-width,0, by=width/60)){
      index <- index+1
      binout <- RectBin2d(datlist[[dat]]$roundx,datlist[[dat]]$roundy,origin,origin,width,width,type="standard")
      binout[[2]]$bintype <- "standard"
      binout[[2]]$distrib <- distrib[distindex]
      spatlossdatRound <- rbind(spatlossdatRound, binout[[2]])
      index <- index+1
      binout <- RectBin2d(datlist[[dat]]$roundx,datlist[[dat]]$roundy,origin,origin,width,width,type="random")
      binout[[2]]$bintype <- "random"
      binout[[2]]$distrib <- distrib[distindex]
      spatlossdatRound <- rbind(spatlossdatRound, binout[[2]])
      print(c(distrib[distindex],width,origin))
    }
  }
}
# write.csv(spatlossdatRound,"LossesByOriginBinwidth.csv", row.names=FALSE)
spatlossdatRound2 <- read.csv("LossesByOriginBinwidth.csv",header=TRUE)

tail(spatlossdatRound2)

# Origin offset in number of binwidths the origin is shifted below data resolution 
spatlossdatRound$originoffset <- spatlossdatRound$originx / spatlossdatRound$widthx
# Calculate minimum proportion of maximum spatial loss at each binwidth
minSpatialLoss <- ddply(spatlossdatRound, .(bintype,distrib,widthx),summarise,
                        minSL = min(propSpatialLoss),
                        x = originx[which.min(propSpatialLoss)] )
# Check out proportion of maximum spatial loss for each bin size and origin location combination
qplot(originoffset, widthx, geom="point", data=spatlossdatRound,size=propSpatialLoss, facets=distrib~bintype)


# proportion of spatial loss for each 
# distribution/binning/binwidth/origin offset combination
# points overlay for origin shift with minimum loss for each binwidth 
qplot(originx,propSpatialLoss,geom="path",group=widthx,colour=as.factor(widthx),
      data=spatlossdatRound)+  geom_point(aes(x=x,y=minSL), data=minSpatialLoss) +
  facet_grid(distrib~bintype,scales="free_y")




lossOrigin <- read.csv("LossesByOriginBinwidth.csv",header=TRUE)

qplot(originx,totalSpatialLoss,geom="path",group=distrib, colour=distrib,
      data=subset(lossOrigin,bintype=="standard"), size=I(1))+ 
  facet_grid(distrib~widthx) + theme_bw() + 
  theme(panel.grid= element_blank())+
  geom_vline(xintercept = -1, colour=I("darkgray")) +
  ylab("Net Spatial Loss") + xlab("Origin Offset")
  
  



##  Check out the spatial loss for the data that is recorded to resolution = alpha
# standard and random binning with bin sizes that are NOT scalar multiples of the data resolution
losswidths <- data.frame(NULL)
origin <- -alpha/2
index=0
distindex = 0
distrib <- c("unif", "norm", "exp")
widths <- seq(alpha,4*alpha,by=alpha/5)
for (dat in 1:3){
  distindex <- distindex + 1
  for (width in widths){
    index <- index+1
    binout <- RectBin2d(datlist[[dat]]$roundx,datlist[[dat]]$roundy,origin,origin,width,width,type="standard")
    binout[[2]]$bintype <- "standard"
    binout[[2]]$distrib <- distrib[distindex]
    losswidths <- rbind(losswidths, binout[[2]])
    index <- index+1
    binout <- RectBin2d(datlist[[dat]]$roundx,datlist[[dat]]$roundy,origin,origin,width,width,type="random")
    binout[[2]]$bintype <- "random"
    binout[[2]]$distrib <- distrib[distindex]
    losswidths <- rbind(losswidths, binout[[2]])
    print(c(distrib[distindex],width,origin))
  }
}

# write.csv(losswidths,"losswidths.csv", row.names=FALSE)
# losswidths <- read.csv("losswidths.csv",header=TRUE)

# proportion of spatial loss for each binwidth combination
# facetted by distribution/bin type

losswidths$DistLabel <- rep(c("Uniform","Normal","Exponential"),each=nrow(losswidths)/3)
qplot(widthx,propSpatialLoss,geom="path",color=bintype,data=losswidths) +
  facet_grid(.~DistLabel,scales="free_y") + theme_bw()






### Binned Scatterplots for Uniform Data to Show stripes
# create binned data for standard and random with binwidth 5
binoutfine <- RectBin2d(datlist[[1]]$x,datlist[[1]]$y,0,0,5,5,type="standard")
binoutcoarse <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,5,5,type="standard")
binoutfineRand <- RectBin2d(datlist[[1]]$x,datlist[[1]]$y,0,0,5,5,type="random")
binoutcoarseRand <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,5,5,type="random")

allbindat <- rbind(binoutfine[[1]],binoutcoarse[[1]],binoutfineRand[[1]],binoutcoarseRand[[1]])
allbindat$SimDataType <- rep(c("Fine","Coarse","Fine","Coarse"), each=nrow(binoutfine[[1]]))
allbindat$BinningType <- rep(c("Standard","Random"), each=2*nrow(binoutfine[[1]]))

qplot(binxs, binys, geom="tile", fill=binfreq, data=allbindat) +
  facet_grid(SimDataType ~ BinningType)+
  ggtitle("Binned Scatterplots of Uniform Data with 5X5 Bins") +  
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits=c(0, 450),name="Bin Count") + 
  theme_bw() + theme(legend.position="bottom",aspect.ratio=1) +
  xlab("X") + ylab("Y")

#compare plots
p1 <- qplot(binxs, binys, geom="tile", fill=binfreq, data=binoutfine[[1]]) +
  ggtitle("Fine Data") +  
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits=c(0, 450)) + 
  theme_bw() + theme(legend.position="none",plot.margin = unit(c(1,.1,.1,1),"cm"),
                     axis.title.y = element_text(size = rel(1.2))) +
  xlab("") + ylab("Standard Binning")
p2 <- qplot(binxs, binys, geom="tile", fill=binfreq, data=binoutcoarse[[1]]) +  
  ggtitle("Coarse Data") + 
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits=c(0, 450)) + 
  theme_bw() + theme(legend.position="none",plot.margin = unit(c(1,1,.1,.1),"cm")) +
  xlab("") + ylab("")
p3 <- qplot(binxs, binys, geom="tile", fill=binfreq, data=binoutfineRand[[1]]) +
  ggtitle("") +  
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits=c(0, 450)) + 
  theme_bw() + theme(legend.position="none",plot.margin = unit(c(.1,.1,1,1),"cm"),
                     axis.title.y = element_text(size = rel(1.2))) +
  xlab("") + ylab("Random Binning")
p4 <- qplot(binxs, binys, geom="tile", fill=binfreq, data=binoutcoarseRand[[1]]) +
  ggtitle("") + 
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits=c(0, 450)) + 
  theme_bw() + theme(legend.position="none",plot.margin = unit(c(.1,1,1,.1),"cm")) +
  xlab("") + ylab("")
grid.arrange(p1,p2,p3,p4,nrow=2)



### Create data to demonstrate stripes from improper binwidths
# binoutcoarse1 <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,1,1,type="standard")[[1]]
# binoutcoarse4 <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,4,4,type="standard")[[1]]
# binoutcoarse5 <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,5,5,type="standard")[[1]]
# binoutcoarse1Rand <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,1,1,type="random")[[1]]
# binoutcoarse4Rand <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,4,4,type="random")[[1]]
# binoutcoarse5Rand <- RectBin2d(datlist[[1]]$roundx,datlist[[1]]$roundy,0,0,5,5,type="random")[[1]]
# 
# binoutcoarse1$BinSize <- "1 X 1 Unit" ; binoutcoarse1$BinningType <- "Standard" 
# binoutcoarse4$BinSize <- "4 X 4 Unit" ; binoutcoarse4$BinningType <- "Standard" 
# binoutcoarse5$BinSize <- "5 X 5 Unit" ; binoutcoarse5$BinningType <- "Standard" 
# binoutcoarse1Rand$BinSize <- "1 X 1 Unit" ; binoutcoarse1Rand$BinningType <- "Random" 
# binoutcoarse4Rand$BinSize <- "4 X 4 Unit" ; binoutcoarse4Rand$BinningType <- "Random" 
# binoutcoarse5Rand$BinSize <- "5 X 5 Unit" ; binoutcoarse5Rand$BinningType <- "Random" 
# 
# allbindat <- rbind(binoutcoarse1,binoutcoarse4,binoutcoarse5,
#                    binoutcoarse1Rand,binoutcoarse4Rand,binoutcoarse5Rand)
# write.csv(allbindat,"StripeBinData.csv", row.names=FALSE)
allbindat <- read.csv("StripeBinData.csv",header=TRUE)



# create binned data for fine exponential with 5X5 bins at good and bad origin
binoutExp5 <- RectBin2d(datlist[[3]]$x,datlist[[3]]$y,0,0,5,5,type="standard")
binoutExp5offset <- RectBin2d(datlist[[3]]$x,datlist[[3]]$y,-4,-4,5,5,type="standard")


p1 <- qplot(binxs, binys, geom="tile", fill=log(binfreq), data=binoutExp5[[1]]) +  
  ggtitle("") +
  scale_fill_gradient(low="#56B1F7", high="#132B43") + 
  theme_bw() + theme(legend.position="none",plot.margin = unit(c(1,1,.1,.1),"cm")) +
  xlab("") + ylab("")

p2 <- qplot(binxs, binys, geom="tile", fill=log(binfreq), data=binoutExp5offset[[1]]) +  
  ggtitle("") + 
  scale_fill_gradient(low="#56B1F7", high="#132B43") + 
  theme_bw() + theme(legend.position="none",plot.margin = unit(c(1,1,.1,.1),"cm")) +
  xlab("") + ylab("")

grid.arrange(p1,p2,nrow=1, sub="X")

binoutExp5offset[[2]][5] / binoutExp5[[2]][5] 



### Try out frequency binning on the same data using standard/quantile/log-count binning
head(datlist[[2]])
binout1 <- RectBin2d(datlist[[2]]$x,datlist[[2]]$y,0,0,5,5,type="standard")
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


### Loop over ncolors 1 to 7 and also over data distribution
# find freq loss from standard, quantile, logCount
freqLossSimDat <- NULL
distrib <- c("Uniform", "Normal", "Exponential")
for (distidx in 1:3){
  binout1 <- RectBin2d(datlist[[distidx]]$x,datlist[[distidx]]$y,0,0,5,5,type="standard")
  for (ncol in 1:7){
    binoutFreq <- rbind(freqBin(binout1, binType="standard", ncolor=ncol, logCount=FALSE)[[2]],
                        freqBin(binout1, binType="quantile", ncolor=ncol, logCount=FALSE)[[2]], 
                        freqBin(binout1, binType="standard", ncolor=ncol, logCount=TRUE)[[2]])
    binoutFreq$FreqBinType <- c("Standard", "Quantile", "Log Count")
    binoutFreq$DistType <- distrib[distidx]
    binoutFreq$ncolor <- ncol 
    freqLossSimDat <- rbind(freqLossSimDat,binoutFreq)
  }
}

freqLossSimDat$FreqBinType <- factor(freqLossSimDat$FreqBinType , levels=c("Standard", "Quantile", "Log Count"))
plotdat <- freqLossSimDat[(freqLossSimDat$FreqBinType!="Log Count" & freqLossSimDat$DistType=="Exponential"),]
p1 <- qplot(ncolor, totalFreqLoss, geom="path",linetype=FreqBinType,
      data=plotdat , color = FreqBinType, size=I(1)) +
  facet_grid(.~DistType, scales="free_y") +  theme_bw() +
  scale_x_continuous(breaks=c(1,3,5,7)) + ylim(0,max(plotdat$totalFreqLoss)) +
  xlab("") + ylab("Total Frequency Loss") + 
  theme(legend.position="bottom")
plotdat <- freqLossSimDat[(freqLossSimDat$FreqBinType!="Log Count" & freqLossSimDat$DistType=="Normal"),]
p2 <-  qplot(ncolor, totalFreqLoss, geom="path",linetype=FreqBinType,
             data=plotdat , color = FreqBinType, size=I(1)) +
  facet_grid(.~DistType, scales="free_y") +  theme_bw() +
  scale_x_continuous(breaks=c(1,3,5,7)) + ylim(0,max(plotdat$totalFreqLoss)) +
  xlab("") + ylab("Total Frequency Loss") + 
  theme(legend.position="bottom")
plotdat <- freqLossSimDat[(freqLossSimDat$FreqBinType!="Log Count" & freqLossSimDat$DistType=="Uniform"),]
p3 <- qplot(ncolor, totalFreqLoss, geom="path", linetype=FreqBinType,
             data=plotdat , color = FreqBinType, size=I(1)) +
  facet_grid(.~DistType, scales="free_y") +  theme_bw() +
  scale_x_continuous(breaks=c(1,3,5,7)) + ylim(0,max(plotdat$totalFreqLoss)) +
  xlab("") + ylab("Total Frequency Loss") + 
  theme(legend.position="bottom")

grid.arrange(arrangeGrob(p1,p2,p3, nrow=1), sub=textGrob("Number of Frequency Bins"))



binoutFreqGroups <- freqBin(binout1, binType="standard", ncolor=5, logCount=FALSE)
binoutQuantFreqGroups <- freqBin(binout1, binType="quantile", ncolor=5, logCount=FALSE)
binoutLogFreqGroups <- freqBin(binout1, binType="standard", ncolor=5, logCount=TRUE)