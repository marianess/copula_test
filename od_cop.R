
# Try to fit bivariate copula to 
# max data. choise of
# elements is arbitrary as
# the goal is to learn R routine.

library(compositions)
library(copula)
library(VineCopula)
library(dplyr)

hem <- read.csv("hem_sorted.txt", header=TRUE, sep="\t")
hem <- hem %>% filter(Drillhole != "RD1988")
hem.log <- log1p(hem[,8:53])


# plot several scatterplots to choose
# variables for expiremental moelling
plot(hem.log$Ta181,hem.log$Nb93)
plot(hem.log$Nb93,hem.log$W182)
plot(hem.log$Co59,hem.log$Pb208)
plot(hem.log$V51,hem.log$GaAs75)
plot(hem.log$La139,hem.log$Ni60)
plot(hem.log$Mn55,hem.log$Ga69)

###############
#             #  
# first pair  #
#             #
###############

# plot pairs
plot(hem.log$La139,hem.log$Ni60)
abline(lm(hem.log$La139~hem.log$Ni60),col='red',lwd=1)
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






