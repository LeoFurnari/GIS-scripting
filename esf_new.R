tool_exec <- function(in_params, out_params)
{
##
library(spdep)
library(RColorBrewer)
library(classInt)
library(MASS)

## ------------------------------------------------------------------------
wdPath <- in_params[[1]]
setwd(wdPath)

## ------------------------------------------------------------------------
csvFile <- in_params[[2]]
pr.f <- read.csv(csvFile)
galFile <- in_params[[3]]
pr.nb <- read.gal(galFile)

pr.listw <- nb2listw(pr.nb, style="W")
pr.listb <- nb2listw(pr.nb, style="B")

## ------------------------------------------------------------------------
n <- length(pr.nb)
M <- diag(n) - matrix(1,n,n)/n
B <- listw2mat(pr.listb)
MBM <- M %*% B %*% M
eig <- eigen(MBM,symmetric=T)
EV <- as.data.frame( eig$vectors[,eig$values/eig$values[1]>0.25])
colnames(EV) <- paste("EV", 1:NCOL(EV), sep="")

## ------------------------------------------------------------------------
##send the list of variables in pr.f to the ArcGIS GUI
##Have user select rain_mean and nofarms_07 values manu

##-------------------------------------------------------------------------
dependentVariable <- in_params[[4]]
independentVariable <- in_params[[5]]
fileOpen <- arc.open(csvFile)
featureSelect <- arc.select(fileOpen, c(dependentVariable, independentVariable))
print(featureSelect)
rain <- pr.f[[independentVariable]]

farm.den07 <- pr.f[[dependentVariable]]/pr.f$area
y.fd <- (farm.den07 - 0.12)^0.38
lm.fd <- lm(y.fd ~ rain)
lm.fd.s <- summary(lm.fd)

s2 <- round(lm.fd.s$sigma^2,5)
c1 <- round(0.5 * (-0.25+(1/0.38-2+1.5)^2),5)

lm.full <- lm(y.fd ~ rain + ., data=EV)
lm.sf <- lm.sf <- step(lm.full, direction='forward', lm(y.fd ~rain, data=EV), 0.1, verbose=F)

pred.sf <- lm.sf$fitted

s2.sf <- round(summary(lm.sf)$sigma^2,5)
y.e.sf <- pred.sf^(1/0.38) + c1*s2.sf + 0.12
lm.sf.bt <- lm(farm.den07 ~ y.e.sf)
summary(lm.sf.bt)


## ------------------------------------------------------------------------
summary(lm.sf)$r.squared

plot(y.e.sf, farm.den07, pch=20)
abline(lm.sf.bt)

}


