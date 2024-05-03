

load(file = 'dataMerged6.Rdata')

dim(dataMerged6)
# Get the desired subset of data
BoneMerged7 <- subset(dataMerged6, HSSEX==2 & DMARACER==1 & HSAGEIR>=50)
summary(BoneMerged7)
cor(BoneMerged7)
Idarm <- which((BoneMerged7$BMPARMC==8888.00))
BoneMerged8 <- BoneMerged7[-1*Idarm,]
Idwaist <- which((BoneMerged8$BMPWAIST==88888.0))
BoneMerged9 <- BoneMerged8[-1*Idwaist,]
Idbutt <- which((BoneMerged9$BMPBUTTO==88888.0))
BoneMerged10 <- BoneMerged9[-1*Idbutt,]
summary(BoneMerged10)
dim(BoneMerged10)


# Get subset of younger demographic for BMD data of a healthly individual
BoneD1 <- subset(dataMerged6,  HSAGEIR>=20 & HSAGEIR<=30)
dim(BoneD1)

save(BoneMerged10,file="BoneMerged10.RData")
save(BoneD1, file="BoneD1.RData")

#### preparation data


## BMPARMC and BMPBUTTO with corr = 0.81
x11 <- BoneMerged10$BMPARMC
x12 <- BoneMerged10$BMPBUTTO

x11 <- (x11-mean(x11))/sd(x11)
x12 <- (x12-mean(x12))/sd(x12)





##  Response Variable

Y0 <- as.matrix(BoneMerged10$BDRTOBMD)

##  Response Variable standardization

Y0 <- (Y0 - mean(Y0))/sd(Y0)

## BMD of norm
BMDavg <- mean(BoneD1$BDRTOBMD)


## sd of norm
BMDstd <- sd(BoneD1$BDRTOBMD)

## Set mu0 to be BMDavg - BMDstd
mu0 <- BMDavg - BMDstd

BMPARMC_BMPBUTTO<- data.frame(X1=x11,X2=x12,Y=Y0)

save(BMPARMC_BMPBUTTO,file="BMPARMC_BMPBUTTO.RData")


