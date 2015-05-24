### Preliminaries 

library(ggplot2)
library(mvtnorm)
library(plyr)
library(Hmisc)

#----------------------------------------------------
### Binning Functions
#----------------------------------------------------

## Standard 1d Rectangular Binning
StandRectBin1d <- function(xs, origin, width){
  binxs <- origin + width*( ceiling((xs-origin)/width) - 0.5 )
  binxs[xs == origin] <- origin + width/2
  return(binxs)
}

## Random 1d Rectangular Binning
# Establish bounding bin centers for each observation
# Then use Unif(0,1) draws compared to assignment probs to allocate
# Reassignment for values below first center (bx1) or
#     above highest center (bxJ)
RandRectBin1d <- function(xs, origin, width){
  bx1 <- origin + width/2
  bxJ <- origin + width*(floor((max(xs)-bx1)/width) + .5)
  lbs <- bx1 + width*floor((xs-bx1)/width)
  ubs <- bx1 + width*ceiling((xs-bx1)/width)
  # initially assign all values to upper bound
  binxs <- ubs
  # then use runif to reassign based on distance from 
  plower <- (ubs - xs)/width
  lowerindex <- (plower > runif(length(xs), 0, 1))
  binxs[lowerindex] <- lbs[lowerindex]
  binxs[xs < bx1] <- bx1
  binxs[xs > bxJ] <- bxJ
  return(binxs)
}


## Quantile 1d Binning
# used for binning the counts values by quantile
# define vector of counts and number of bin
QuantBin1d <- function(cs, nbin){
  bincenters <- quantile(cs, seq(1/(2*nbin) , 1- 1/(2*nbin), by=1/nbin))
  binbounds <- quantile(cs, seq(0 , 1, by=1/nbin))
  bincs <- rep(bincenters[1],length(cs))
   for (i in 2:length(bincenters)){
     bincs[(binbounds[i] < cs & cs <= binbounds[i+1])] <- bincenters[i]
   }   
  return(bincs)
}


## 2d Rectangular Binning (for either Standard or Random)
# standard is straight forward extension of 1d binning
# random binning needs post processing to calculate minimum spatial loss
RectBin2d <- function(xs,ys, originx, originy, widthx, widthy, type="standard"){
  tempdat <- data.frame(xs = xs, ys=ys,
                        binxs = StandRectBin1d(xs,originx,widthx),
                        binys = StandRectBin1d(ys,originy,widthy))
  
  if(type=="random"){
    tempdat<- data.frame( tempdat,
                          randbinxs = RandRectBin1d(xs,originx,widthx),
                          randbinys = RandRectBin1d(ys,originy,widthy))
    tempdat$binname <- paste(tempdat$binxs,tempdat$binys)
    tempdat$randbinname <- paste(tempdat$randbinxs,tempdat$randbinys)
    tempdat$randdist <- with(tempdat,sqrt((xs-randbinxs)^2 + (ys-randbinys)^2))
    tempdat$index <- 1:length(xs)
    # points where mismatch between standard and random binning
    mismatchindex <- which(tempdat$binxs != tempdat$randbinxs | tempdat$binys != tempdat$randbinys )
    mmdat <- tempdat[mismatchindex,]
    # identify which need to be swapped to standard for net spatial loss
    
    #loop over all neighboring bin pairs (shift all xs over by 1 then all )
    xbins <- seq(min(c(tempdat$binxs,tempdat$randbinxs)) , max(c(tempdat$binxs,tempdat$randbinxs)), by=widthx)
    ybins <- seq(min(c(tempdat$binys,tempdat$randbinys)) , max(c(tempdat$binys,tempdat$randbinys)), by=widthy)
    nbrs <- data.frame(binxs = rep(rep(xbins,length(ybins)),2) ,
                       binys = rep(rep(ybins,each=length(xbins)),2),
                       nbrsxs = c(rep(xbins+widthx,length(ybins)),rep(xbins,length(ybins))),
                       nbrsys = c(rep(ybins,each=length(xbins)), rep(ybins+widthy,each=length(xbins))) )
    nbrs$binname <- paste(nbrs$binxs,nbrs$binys)
    nbrs$nbrsname <- paste(nbrs$nbrsxs,nbrs$nbrsys)   
    
    swapindex <- NULL    
    for (i in 1:nrow(nbrs)){
      #id points in standard bin i assigned to bin j
      itoj <- which(mmdat$binname == nbrs$binname[i] & mmdat$randbinname == nbrs$nbrsname[i])
      #id points in standard bin j assigned to bin i
      jtoi <- which(mmdat$binname == nbrs$nbrsname[i] & mmdat$randbinname == nbrs$binname[i])
      # number to swap in bins i and j is equal to minimum misplaced
      nswap <- min(length(itoj), length(jtoi))
      # if there are points to swap, then pick the ones with largest distance
      # from point to random bin center
      if(nswap > 0){ 
        swapindex <- c(swapindex,mmdat$index[itoj][order(mmdat$randdist[itoj],decreasing=TRUE)[nswap]])
        swapindex <- c(swapindex,mmdat$index[jtoi][order(mmdat$randdist[jtoi],decreasing=TRUE)[nswap]])
      }
    }
    swapindex <- unique(swapindex)
    tempdat$binxs[!(tempdat$index %in% swapindex)] <- tempdat$randbinxs[!(tempdat$index %in% swapindex)]
    tempdat$binys[!(tempdat$index %in% swapindex)] <- tempdat$randbinys[!(tempdat$index %in% swapindex)]
  }
  outdat <- ddply(tempdat, .(binxs,binys), summarise,
                  binfreq = length(xs),
                  binspatialloss = sum(sqrt((xs-binxs)^2+(ys-binys)^2)) )
  summarydata <- data.frame( originx = originx, originy = originy, 
                             widthx = widthx, widthy = widthy,
                             totalSpatialLoss = sum(outdat$binspatialloss))
  templist <- list(bindat = outdat[,1:3],
                   summarydata)
  return(templist)
}

## Frequency Binning  
# allows for standard or quantile binning of either  bin counts or log bin counts (4 combinations)
# input requires binned data output, number of freq breaks and type of freq binning
# output of frequency bin values, labels and loss are attached the original then returned
freqBin <- function(binout, binType="standard", ncolor, logCount=FALSE){ 
  cs <- binout[[1]]$binfreq
  # add frequency bin centers and bin labels of the form "(a,b]" to binout
  # binning depends on type / log counts
  if (logCount)  cs <- log(binout[[1]]$binfreq)
  if(binType=="standard"){
      #width <- ceiling((max(cs)-min(cs))/ncolor)
      width <- (max(cs)-min(cs))/ncolor 
      binout[[1]]$freqgroup <- round(StandRectBin1d(cs, min(cs) , width),5)
      binout[[1]]$freqlabel <- paste("(",round(binout[[1]]$freqgroup - width/2,1),
                                     ",",round(binout[[1]]$freqgroup + width/2,1),"]",sep="")
      #close interval for smallest counts
      closeidx <- binout[[1]]$freqlabel == min(binout[[1]]$freqlabel)
      binout[[1]]$freqlabel[closeidx] <- paste("[",round(min(binout[[1]]$freqgroup)- width/2,1),
                                     ",",round(min(binout[[1]]$freqgroup) + width/2,1),"]",sep="")
  } 
  if(binType=="quantile"){
    binout[[1]]$freqgroup <- as.numeric(round(QuantBin1d(cs, ncolor),5))
    quantbounds <- unique(quantile(cs, (0:ncolor)/ncolor))
    if (length(quantbounds)-1 < ncolor) warning("two or more quantiles of the data have same value due to many bins with equivalent frequencies, color will be rendered equivalently for bins with matching quantiles")
    binout[[1]]$freqlabel <- cut(cs, breaks=quantbounds, include.lowest=TRUE)
  }
  # Total Freq Loss
  binout[[2]]$totalFreqLoss <- sum((cs - binout[[1]]$freqgroup)^2)
  # Max Freq Loss based on shading at data average
#  binout[[2]]$maxFreqLoss <- sum((cs - mean(cs))^2)
  # Proportion of Max Freq Loss
#  binout[[2]]$propFreqLoss <- sum((cs - binout[[1]]$freqgroup)^2) /  sum((cs - mean(cs))^2)
  return(binout)
}

# #----------------------------------------------------
# #### Function Tests
# #----------------------------------------------------
# 
# ## Try 1d binning functions
# xs <- runif(100,0,100)
# cs <- rexp(100,.1)
# 
# # standard binning 
# standbinnedxs <- StandRectBin1d(xs, origin=0, width=10)
# qplot(xs,  standbinnedxs)
# # random binning
# randbinnedxs <- RandRectBin1d(xs, origin=0, width=10)
# qplot(xs,  randbinnedxs)
# # quantile binning 
# quantbinnedcs <- QuantBin1d(cs, 7)
# qplot(cs,  quantbinnedcs)
# 
# 
# ### Try 2d rectangular binning with MVN data
# n=50000
# mu <- 50
# sig <- 15
# xnorm <- rnorm(n, mean=mu, sd=sig)
# ynorm <- rnorm(n, mean=mu, sd=sig)
# qplot(xnorm, ynorm)
# # create binned data for standard and random 
# binout1 <- RectBin2d(xnorm,ynorm,0,0,10,10,type="standard")
# binout2 <- RectBin2d(xnorm,ynorm,0,0,10,10,type="random")
# # take a peek at binned data format
# head(binout1[[1]])
# head(binout2[[2]])
# # compare spatial losses
# binout1[[2]]
# binout2[[2]]
# #compare visualizations
# qplot(binxs, binys, geom="tile", fill=binfreq, data=binout1[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
# qplot(binxs, binys, geom="tile", fill=binfreq, data=binout2[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
# #compare with shade by log counts
# qplot(binxs, binys, geom="tile", fill=log(binfreq+1), data=binout1[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
# qplot(binxs, binys, geom="tile", fill=log(binfreq+1), data=binout2[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
# 
# 
# 
# ### Try out frequency binning on the same data using standard/quantile/log-count binning
# binoutFreqGroups <- freqBin(binout1, binType="standard", ncolor=5, logCount=FALSE)
# binoutQuantFreqGroups <- freqBin(binout1, binType="quantile", ncolor=5, logCount=FALSE)
# binoutLogFreqGroups <- freqBin(binout1, binType="standard", ncolor=5, logCount=TRUE)
# # Peek into how frequency bins are defined
# head(binoutFreqGroups[[1]])
# head(binoutQuantFreqGroups[[1]])
# head(binoutLogFreqGroups[[1]])
# # Compare Frequency Losses (Note I dont know if these are all directly comparable after logging counts)
# binoutFreqGroups[[2]]
# binoutQuantFreqGroups[[2]]
# binoutLogFreqGroups[[2]] 
# # Plot with Orginal Frequencies
# qplot(binxs, binys, geom="tile", fill=binfreq, data=binoutFreqGroups[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend")
# # Plot with Binned Frequencies
# qplot(binxs, binys, geom="tile", fill=freqgroup, data=binoutFreqGroups[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend", breaks=unique(binoutFreqGroups[[1]]$freqgroup))
# # Plot with Quantile Binned Frequencies
# qplot(binxs, binys, geom="tile", fill=freqgroup, data=binoutQuantFreqGroups[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend", breaks=unique(binoutQuantFreqGroups[[1]]$freqgroup))
# # Plot with LogCount Binned Frequencies
# qplot(binxs, binys, geom="tile", fill=freqgroup, data=binoutLogFreqGroups[[1]]) +
#   scale_fill_gradient(low="#56B1F7", high="#132B43", guide="legend", breaks=unique(binoutLogFreqGroups[[1]]$freqgroup))





### Create binning function that does not calculate loss
# to be used to track computation time for binning

## 2d Rectangular Binning (for either Standard or Random)
# standard is straight forward extension of 1d binning
# random binning needs post processing to calculate minimum spatial loss
RectBin2dNoLoss <- function(xs,ys, originx, originy, widthx, widthy, type="standard"){
  if(type=="standard"){
  tempdat <- data.frame(xs = xs, ys=ys,
                        binxs = StandRectBin1d(xs,originx,widthx),
                        binys = StandRectBin1d(ys,originy,widthy))
  }
  if(type=="random"){
    tempdat<- data.frame(xs = xs, ys=ys,
                          binxs =  RandRectBin1d(xs,originx,widthx),
                          binys =  RandRectBin1d(ys,originy,widthy))
  }
  outdat <- ddply(tempdat, .(binxs,binys), summarise,
                  binfreq = length(xs))
  return(outdat[,1:3])
}

