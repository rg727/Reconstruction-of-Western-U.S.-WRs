

rm(list=ls())
library(zoo)
library(foreign)
library(nnet)
library(DirichletReg)
library(ncdf4)
library(maps)
library(colorRamps)
library(glmnet)
library(Rsolnp)
library(cmaes)
library(adagio)
library(ctmcd)

args <- commandArgs()
s=as.numeric(args[6]) #Seeding for the Borg Evolutionary Algorithm (Remove if not using Borg)

my.dir <- "/home/rg727"  #Set your directory

setwd(my.dir)
source("DirichReg_optim.R")
source("DirichReg_predict_cv.R")
source("DirichReg_reconstruct.R")
source("borg.R")

#######################################################################################

num.states <- 5

synoptic.state.assignments <- read.table("../data/synoptic_state_assignments.txt")
yr <- unique(synoptic.state.assignments$V2)
n.year <- length(yr)
WR <- array(NA,c(n.year,num.states))
for (i in 1:n.year) {
  WR[i,] <- tabulate(synoptic.state.assignments$V3[synoptic.state.assignments$V2==yr[i]],nbins=num.states)
}

#define WR fractions
WR_frac <- WR / apply(WR,FUN=sum,1)
WR_frac.DR <- DR_data(WR_frac)

#PCA on WRs
mypc <- prcomp(WR,center=TRUE,scale=TRUE)
(mypc$sdev^2)/sum(mypc$sdev^2)
WR.PCs <- mypc$x
#######################################################################################

#get reconstructed SPI
my.nc <- nc_open("../data/SPI_recon_NDJFM_0.5deg_Master_Extended_v2.nc")
print(my.nc)
spi_recon_org <- ncvar_get(my.nc,"spi_recon")
yr_recon <- ncvar_get(my.nc,"yr_recon")
vrsq_org <- ncvar_get(my.nc,"vrsq")
lat <- ncvar_get(my.nc,"lat")
lon <- ncvar_get(my.nc,"lon")
nc_close(my.nc)
lon_lat_org <- expand.grid(lon,lat)


#process SPI data into a matrix
spi_recon <- aperm(spi_recon_org,c(3,1,2))
dim(spi_recon) <- c(dim(spi_recon)[1],prod(dim(spi_recon)[2:3]))
spi_recon <- apply(spi_recon,2,function(x){scale(x)[,1]})
vrsq <- aperm(vrsq_org,c(3,1,2))
dim(vrsq) <- c(dim(vrsq)[1],prod(dim(vrsq)[2:3]))

#drop any grid cells that don't go back to 1400
keep <- apply(spi_recon,c(2),function(x) {length(which(is.na(x)))==0})
spi_recon <- spi_recon[,keep]
lon_lat <- lon_lat_org[keep,]
vrsq <- vrsq[,keep]

#PCs of SPI data
spi_recon.pca <- prcomp(spi_recon,center=TRUE,scale=TRUE)
spi_recon.pcs <- spi_recon.pca$x
colnames(spi_recon.pcs) <- paste("PC",1:ncol(spi_recon.pcs),sep="")
#######################################################################################



##########################################Optimization########################################

leave.k.out <- 10
smooth.list <- c(1,3,5,10)

  #Optimization: change the number of RBFs and can add more NFE
  NFE <- 2000
  span.lower.limit <- 0.004    
  span.upper.limit <- 0.032   
  lowerb <- c(2,span.lower.limit)
  upperb <- c(30,span.upper.limit) 
  mystopeval=NFE
  
  start.time <- Sys.time()
  my.optim <- borg(nvars, nobjs, 0, DirichReg_optim, mystopeval, epsilons=c(0.002,0.001,0.001,0.0003,1), lowerBounds=lowerb,upperBounds=upperb,WR.PCs=WR.PCs,num.states=num.states,WR_frac=WR_frac,WR_frac.DR=WR_frac.DR,yr=yr,yr_recon=yr_recon,spi_recon=spi_recon,smooth.list=smooth.list,use.smooth=TRUE,leave.k.out=leave.k.out)
  end.time <- Sys.time()
  print(end.time - start.time)
  
  filename=paste("optimization_results.rda",sep = "")
  save(my.optim,file=filename)
  
  
  ########################################Reconstruction in Instrumental Period#######################################
  fixed.num.pcs <- 11         #for the baseline (no smoothing) case
  my.param <- c(fixed.num.pcs,9999)  #no smoothing, second parameter doesnt matter
  DirichReg_predict_cv_result <- DirichReg_predict_cv(param=my.param,WR.PCs=WR.PCs,num.states=num.states,WR_frac=WR_frac,WR_frac.DR=WR_frac.DR,
                                                      yr=yr,yr_recon=yr_recon,spi_recon.pcs=spi_recon.pcs,smooth.list=smooth.list,
                                                      use.smooth=FALSE,leave.k.out=leave.k.out)
  
  WR.pred.cv <- DirichReg_predict_cv_result[[1]]
  MSE <- DirichReg_predict_cv_result[[2]]
  log.lik <- DirichReg_predict_cv_result[[3]]
  
  #Example plotting to see Instrumental Period Results
  
  kk <- 5
  par(mfrow=c(2,3),mar=c(2,2,1,1))
  for (i in 1:num.states) {
    plot(yr,WR_frac.DR[,i],main=paste("WR",i,"no smooth"),type="l")
    lines(rollmean(yr,kk),rollmean(WR_frac.DR[,i],kk),lwd=2)
    lines(yr,WR.pred.cv[,i],type="l",col="red",lwd=1)
    lines(rollmean(yr,kk),rollmean(WR.pred.cv[,i],kk),type="l",col="red",lwd=2)
  }
  
  #################Develop Paleo Reconstruction##################################
  recon.out <- DirichReg_reconstruct(param=my.param,WR.PCs=WR.PCs,num.states=num.states,
                                     WR_frac=WR_frac,WR_frac.DR=WR_frac.DR,yr=yr,yr_recon=yr_recon,spi_recon.pcs=spi_recon.pcs,
                                     use.smooth=TRUE) 
  my.recon <- recon.out[[1]]
  recon.model <- recon.out[[2]]
  my.lm <- recon.out[[3]]
  
  my.param <- c(11,9999)
  recon.out <- DirichReg_reconstruct(param=my.param,WR.PCs=WR.PCs,num.states=num.states,
                                     WR_frac=WR_frac,WR_frac.DR=WR_frac.DR,yr=yr,yr_recon=yr_recon,
                                     spi_recon.pcs=spi_recon.pcs,use.smooth=FALSE) 
  
  my.recon.no.smooth <- recon.out[[1]]
  recon.model.no.smooth <- recon.out[[2]]
  
  #Example plotting to see Paleo Period Results
  
  kk <- smooth.list[obj.smooth]
  par(mfrow=c(2,3),mar=c(2,2,1,1))
  for (i in 1:num.states) {
    plot(rollmean(yr_recon,kk),rollmean(my.recon[,i],kk),lwd=2,ylim=c(0,.45),main=paste("WR",i,"  roll",smooth.list[obj.smooth]),type="l")
    lines(rollmean(yr,kk),rollmean(WR_frac[,i],kk),lwd=1,col="red")
    lines(rollmean(yr_recon,kk),rollmean(my.recon.no.smooth[,i],kk),lwd=1,col="orange")
  }
  
  




