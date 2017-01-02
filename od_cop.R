
# Try to fit bivariate copula to 
# max data. choise of
# elements is arbitrary as
# the goal is to learn R routine.

library(compositions)
library(copula)
library(VineCopula)
library(dplyr)
library(VIM)
library(ggplot2)
library(reshape2)

# read data in and impute missing values
hem <- read.csv("hem_sorted.txt", header=TRUE, sep="\t")
#average <- sapply(hem[,8:53],mean,na.rm=TRUE)
hem.imp <- kNN(hem[,8:53], 
           variable = c('P31','As75','Sr88','Zr90','Sb121','Ba137','Hf178', 'Ta181'), 
           k = 1, numFun = median)
#average.imp <- sapply(hem.imp[,1:46],mean)
#dif<- average - average.imp

# drop columns and rename
hem.imp <- subset(hem.imp,select=-c(P31_imp,Ta181_imp,As75_imp,Sr88_imp, 
                                    Zr90_imp,Sb121_imp,Ba137_imp,Hf178_imp))

hem.log <- log1p(hem.imp[,1:46])
hem.log[,47:53] <- hem[,1:7]

# box-cox
df.temp <- hem.imp+1
df.bx = data.frame(apply(df.temp,2,transform_bc))
df.bx[,47:54]<- hem[,1:7]

# data by textures
hem.crs <- df.bx %>% filter(Texture == 'Coarse Hm with U-min inclusions')
hem.rep <- hem.log %>% filter(Texture == "Hm - replacement carbonate? ")
hem.osc <- hem.log %>% filter(Texture == "Oscillatory-zoned")

wilcox <-function(a,b){
      result <- wilcox.test(x=a,y=b, alternative=c("two.sided"), paired = FALSE)
      if (result$p.value>0.05){
            print("do not reject Ho")
            return(result$p.value)
      }
      else{
            print("reject Ho")
            return(result$p.value)
      }
      }


wilcox(hem.crs$Hf178,hem.osc$Hf178)
wilcox(hem.crs$Mg24,hem.osc$Mg24)
wilcox(hem.crs$Al27,hem.osc$Al27)
wilcox(hem.crs$Ti49,hem.osc$Ti49)
wilcox(hem.crs$V51,hem.osc$V51)
wilcox(hem.crs$Ni60,hem.osc$Ni60)
wilcox(hem.crs$Zr90,hem.osc$Zr90)
wilcox(hem.crs$Lu175,hem.osc$Lu175)
wilcox(hem.crs$Zr90,hem.rep$Zr90)



# plot several scatterplots to choose
# variables for expiremental moelling
plot(hem.log$Ta181,hem.log$Nb93)
plot(hem.log$Nb93,hem.log$W182)
plot(hem.log$Co59,hem.log$Pb208)
plot(hem.log$V51,hem.log$As75)
plot(hem.log$La139,hem.log$Ni60)
plot(hem.log$Mn55,hem.log$Ga69)

p <- hem.log %>% ggplot(aes(x=hem.log$W182, y=hem.log$Mo95))+geom_point(aes(color=factor(hem.log$Texture)))
p

# plot pairs
p <- hem.crs %>% ggplot(aes(x=W182, y=U238))+geom_point()+geom_smooth(method='lm')
p
cor(hem.crs$W182,hem.crs$U238,method='spearman')

# create subgroups of elements
grans <-(hem.crs$W182+hem.crs$U238+hem.crs$Sn118+hem.crs$Mo95)
rey <- log1p(hem$sumREY[hem$Texture == "Coarse Hm with U-min inclusions"])
hfse <-(hem.crs$Th232 + hem.crs$Nb93 + hem.crs$V51+hem.crs$Ti49)
mafic <-(hem.crs$Ni60 + hem.crs$Ta181 +hem.crs$Co59)
groups <- as.data.frame(cbind(rey,grans,hfse,mafic))

p <- groups %>% ggplot(aes(x=mafic, y=hfse))+geom_point()+geom_smooth(method='lm')
p

p <- groups %>% ggplot(aes(x=grans, y=hfse))+geom_point()+geom_smooth(method='lm')
p



###############
#             #  
#    mafic    #
#     vs      #
#    hfse     #
#             #
###############
# calculate mean and sd for edf
mafic_mu <- mean(mafic)
mafic_sd <- sd(mafic)
hfse_mu <- mean(hfse)
hfse_sd <- sd(hfse)


# hist 1
g <- mafic
hist(g, density=20, breaks=20, prob=TRUE, 
     xlab="mafic", ylim=c(0, 1), 
     main="normal curve over histogram")
curve(dnorm(x, mean=mafic_mu, sd=mafic_sd), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

# hist 2
g <- hfse
hist(g, density=20, breaks=20, prob=TRUE, 
     xlab="hfse", ylim=c(0, 1.8), 
     main="normal curve over histogram")
curve(dnorm(x, mean=hfse_mu, sd=hfse_sd), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

###############
#             #  
#   W vs Mo   #
#   coarse    #
#  hematite   #
#             #
###############

# plot pairs
p <- hem.crs %>% ggplot(aes(x=W182, y=Mo95))+geom_point()+geom_smooth(method='lm')
cor(hem.crs$W182,hem.crs$Mo95,method='spearman')

# calculate mean and sd for edf
W_mu <- mean(hem.crs$W182)
W_sd <- sd(hem.crs$W182)
Mo_mu <- mean(hem.crs$Mo95)
Mo_sd <- sd(hem.crs$Mo95)

# hist 1
g <- hem.crs$W182
hist(g, density=20, breaks=20, prob=TRUE, 
     xlab="W", ylim=c(0, 1), 
     main="normal curve over histogram")
curve(dnorm(x, mean=W_mu, sd=W_sd), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
 
# hist 2
g <- hem.crs$Mo95
hist(g, density=20, breaks=20, prob=TRUE, 
     xlab="Mo", ylim=c(0, 1.8), 
     main="normal curve over histogram")
curve(dnorm(x, mean=Mo_mu, sd=Mo_sd), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

# decide which copula to select
u <- pobs(as.matrix(hem.crs$W182))
v <- pobs(as.matrix(hem.crs$Mo95))
selectedCopula <- BiCopSelect(u,v,familyset=NA)
selectedCopula
plot(u,v)

# Bivariate copula: t (par = 0.7, par2 = 2.62, tau = 0.49)
# create a 2D tCopula
t.cop <- tCopula(dim=2)
set.seed(500)
m <- pobs(as.matrix(cbind(hem.crs$W182,hem.crs$Mo95)))

# fitting copula
fit <- fitCopula(t.cop,m,method='ml')

# Find and save the coefficients
coef(fit)
trho = coef(fit)[1]
tdf  = coef(fit)[2]
persp(tCopula(trho, dim = 2, df = tdf),dCopula)

# Create a contour plot of a copula
contour(tCopula(trho, dim = 2, df = tdf),dCopula, main = "Student t")

# build copula and sample from it 253 samples
C <- rCopula(253, tCopula(dim=2,trho,df=tdf))
plot(C[,1],C[,2],pch='.', col="darkblue")
cor(C, method='spearman')

# simulate data
copula_dist <- mvdc(copula=tCopula(dim=2,trho,df=tdf), margins=c("norm","norm"),
                    paramMargins=list(list(mean=W_mu, sd=W_sd),
                                      list(mean=Mo_mu, sd=Mo_mu)))


sim <- rMvdc(253,copula_dist)

# plotting results
plot(hem.crs$W182,hem.crs$Mo95,main='observed')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

# goodness of fit test
mydata <- cbind(hem.crs$W182,hem.crs$Mo95)
gf <- gofCopula(normalCopula(dim = 2), as.matrix(mydata), N = 50)
gf

###############
#             #  
# first pair  #
#             #
###############

# plot pairs
plot(hem.log$La139,hem.log$Ni60)
abline(lm(hem.log$La139~hem.log$Ni60),col='red')
cor(hem.log$La139,hem.log$Ni60,method='spearman')

# decide which copula to select
u <- pobs(as.matrix(hem.log$La139))
v <- pobs(as.matrix(hem.log$Ni60))
selectedCopula <- BiCopSelect(u,v,familyset=NA)
selectedCopula

# Bivariate copula: BB7 
# (par = 1.12, par2 = 0.91, tau = 0.34)

# Bivariate normal
# create a 2D tCopul. pobs() converts data to uniform d.
n.cop <- normalCopula(dim=2)
set.seed(500)
m <-pobs(as.matrix(cbind(hem.log$La139,hem.log$Ni60)))

# fitting copula
fit <- fitCopula(n.cop,m,method='ml')
rho = coef(fit)[1]
persp(normalCopula(dim=2,rho),dCopula)


# build copula and sample from it 454 samples
u<- rCopula(454, normalCopula(dim=2,rho))
plot(u[,1],u[,2],pch='.', col="darkblue")
cor(u, method='spearman')

# model the margnals
La_mu <- mean(hem.log$La139)
La_sd <- sd(hem.log$La139)
Ni_mu <- mean(hem.log$Ni60)
Ni_sd <- sd(hem.log$Ni60)

# hist 1
g <- hem.log$La139
hist(g, density=20, breaks=10, prob=TRUE, 
     xlab="La", ylim=c(0, 0.7), 
     main="normal curve over histogram")
curve(dnorm(x, mean=La_mu, sd=La_sd), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

# hist 2
g <- hem.log$Ni60
hist(g, density=20, breaks=10, prob=TRUE, 
     xlab="Ni", ylim=c(0, 0.8), 
     main="normal curve over histogram")
curve(dnorm(x, mean=Ni_mu, sd=Ni_sd), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")


# simulate data
copula_dist <- mvdc(copula=normalCopula(rho,dim=2), margins=c("norm","norm"),
                    paramMargins=list(list(mean=La_mu, sd=La_sd),
                                      list(mean=Ni_mu, sd=Ni_mu)))


sim <- rMvdc(454,copula_dist)

# plotting results
plot(hem.log$La139,hem.log$Ni60,main='observed')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)




###############
#             #  
# second pair #
#             #
###############

# plot pairs
plot(hem.log$Mn55,hem.log$Ga69)
abline(lm(hem.log$Mn55~hem.log$Ga69),col='red',lwd=1)
cor(hem.log$Mn55,hem.log$Ga69,method='pearson')

# decide which copula to select
u <- pobs(as.matrix(hem.log$Mn55))
v <- pobs(as.matrix(hem.log$Ga69))
selectedCopula <- BiCopSelect(u,v,familyset=NA)
selectedCopula

# Bivariate copula: Rotated Tawn type 2 180 degrees 
# (par = 1.64, par2 = 0.16, tau = 0.11) 

# Bivariate normal
# create a 2D tCopul. pobs() converts data to uniform d.
n.cop <- normalCopula(dim=2)
set.seed(500)
m <-pobs(as.matrix(cbind(hem.log$Mn55,hem.log$Ga69)))

# fitting copula
fit <- fitCopula(n.cop,m,method='ml')
rho = coef(fit)[1]
persp(normalCopula(dim=2,rho),dCopula)

# build copula and sample from it 3965 samples
u<- rCopula(454, normalCopula(dim=2,rho))
plot(u[,1],u[,2],pch='.', col="darkblue")
cor(u, method='spearman')


# model the margnals
Mn_mu <- mean(hem.log$Mn55)
Mn_sd <- sd(hem.log$Mn55)
Ga_mu <- mean(hem.log$Ga69)
Ga_sd <- sd(hem.log$Ga69)

# plot hists and normal dist with mu and sigma
# first element
hist(hem.log$Mn55,breaks=30,main='Mn',col='cyan')

hist(hem.log$Ga69,breaks=30,main='Ga',col='cyan')

# simulate data
copula_dist <- mvdc(copula=normalCopula(rho,dim=2), margins=c("norm","norm"),
                    paramMargins=list(list(mean=Mn_mu, sd=Mn_sd),
                                      list(mean=Ga_mu, sd=Ga_sd)))


sim <- rMvdc(454,copula_dist)

# compute density
pdf_mvd <- dMvdc(sim,copula_dist)

# compute cdf
cdf_mvd <- pMvdc(sim,copula_dist)

# # 3D plain scatterplot of the multivariate distribution
# par(mfrow = c(1, 1))
# scatterplot3d(sim[,1],sim[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
# scatterplot3d(sim[,1],sim[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
# persp(copula_dist, dMvdc, xlim = c(0, 8), ylim=c(0, 8), main = "Density")
# contour(copula_dist, dMvdc, xlim = c(0, 8), ylim=c(0, 8), main = "Contour plot")
# persp(copula_dist, pMvdc, xlim = c(0, 8), ylim=c(0, 8), main = "CDF")
# contour(copula_dist, pMvdc, xlim = c(0, 8), ylim=c(0, 8), main = "Contour plot")
# 

# plotting results
plot(hem.log$Mn55,hem.log$Ga69,main='observed')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

mydata <- cbind(hem.log$Mn55,hem.log$Ga69)
gf <- gofCopula(normalCopula(dim = 2), as.matrix(mydata), N = 50)
gf



