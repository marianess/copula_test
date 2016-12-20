# generating mv normal distribution of 3 random variables
#given tho cov matrix sifma

library(MASS)
set.seed(100)

m <- 3
n <- 2000
# cov matrix sigma: var = 1
sigma <- matrix(c(1, 0.4, 0.2,
                 0.4, 1, -0.8,
                 0.2, -0.8, 1),nrow=3)

z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma,empirical=T)

# use spearman cor coefficient
library(psych)
cor(z,method = 'spearman')
pairs.panels(z)

# make 3 distributions uniform
u <- pnorm(z)
pairs.panels(u)

library(rgl)
plot3d(u[,1],u[,2],u[,3],pch=20,col='navyblue')

# define new distribution functions
x1 <- qgamma(u[,1],shape=2,scale=1)
x2 <- qbeta(u[,2],2,2)
x3 <- qt(u[,3],df=5)
plot3d(x1,x2,x3,pch=20,col='blue')

# look at spearman correlation
df <- cbind(x1,x2,x3)
pairs.panels(df)
cor(df,meth='spearman')

#####################
#                   #
# real-life example #
#                   #
#####################

library(copula)
library(VineCopula)

cree <- read.csv('cree_r.csv',header=F)$V2
yahoo <- read.csv('yahoo_r.csv',header=F)$V2

# plot regressionline
plot(cree,yahoo,pch='.')
abline(lm(yahoo~cree),col='red',lwd=1)
cor(cree,yahoo,method='spearman')

# decide which copula to select
u <- pobs(as.matrix(cbind(cree,yahoo)))[,1]
v <- pobs(as.matrix(cbind(cree,yahoo)))[,2]
selectedCopula <- BiCopSelect(u,v,familyset=NA)
selectedCopula

# create a 2D tCopul. pobs() converts data to uniform d.
t.cop <- tCopula(dim=2)
set.seed(500)
m <-pobs(as.matrix(cbind(cree,yahoo)))

# fitting copula
fit <- fitCopula(t.cop,m,method='ml')
coef(fit)

# save coefficients of fitted copula
rho <- coef(fit)[1]
df <- coef(fit)[2]
persp(tCopula(dim=2,rho,df=df),dCopula)

# build copula and sample from it 3965 samples
u<- rCopula(3965, tCopula(dim=2,rho,df=df))
plot(u[,1],u[,2],pch='.', col="darkblue")
cor(u, method='spearman')

# model the margnals
cree_mu <- mean(cree)
cree_sd <- sd(cree)
yahoo_mu <- mean(yahoo)
yahoo_sd <- sd(yahoo)

# plot hists and normal dist with mu and sigma
# cree
hist(cree,breaks=80,main='Cree returns',freq=F,density=30,col='cyan',ylim=c(0,20),xlim=c(-0.2,0.3))
lines(seq(-0.5,0.5,0.01),dnorm(seq(-0.5,0.5,0.01),cree_mu,cree_sd),col='red',lwd=2)
legend('topright',c('Fitted normal'),col=c('red'),lwd=2)

#yahoo
hist(yahoo,breaks=80,main='Yahoo returns',density=30,col='cyan',freq=F,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(seq(-0.5,0.5,0.01),dnorm(seq(-0.5,0.5,0.01),yahoo_mu,yahoo_sd),col='red',lwd=2)
legend('topright',c('Fitted normal'),col=c('red'),lwd=2)

# simulate data
copula_dist <- mvdc(copula=tCopula(rho,dim=2,df=df), margins=c("norm","norm"),
                    paramMargins=list(list(mean=cree_mu, sd=cree_sd),
                                      list(mean=yahoo_mu, sd=yahoo_sd)))


sim <- rMvdc(3965,copula_dist)


# plot
plot(cree,yahoo,main='Returns')
points(sim[,1],sim[,2],col='red')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

#####################
#                   #
#      example      #
#                   #
#####################

library(scatterplot3d)
library(grid)
library(ggplot2)
set.seed(235)

# If I know parameters already:
# in normal case param=rho
# in t case param = rho, df = degrees of freadom
normal <- normalCopula(param = 0.7, dim = 2)
stc <- tCopula(param = 0.8, dim = 2, df = 2)
frank <- frankCopula(dim = 2, param = 8)
gumbel <- gumbelCopula(dim = 3, param = 5.6)
clayton <- claytonCopula(dim = 4, param = 19)

# Select the copula
cp <- claytonCopula(param = c(3.4), dim = 2)

# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
multivariate_dist <- mvdc(copula = cp,
                          margins = c("norm", "t"),
                          paramMargins = list(list(mean = 2, sd=3),list(df = 2)))

print(multivariate_dist)

# Generate random samples from Copula
# rCopula() for asingl copula
# rMvdc() for multivariate distribution

# Generate random samples
fr <- rCopula(1000, frank)
gu <- rCopula(1000, gumbel)
cl <- rCopula(1000, clayton)

# plot random samples
p1 <- qplot(fr[,1], fr[,2], colour = fr[,1], 
            main="Frank copula random samples theta = 8", 
            xlab = "u", ylab = "v")

p2 <- qplot(gu[,1], gu[,2], colour = gu[,1], 
            main="Gumbel copula random samples theta = 5.6", 
            xlab = "u", ylab = "v") 

p3 <- qplot(cl[,1], cl[,2], colour = cl[,1], 
            main="Clayton copula random samples theta = 19", 
            xlab = "u", ylab = "v")

samples <- rMvdc(2000, multivariate_dist)

scatterplot3d(samples[,1], samples[,2], color = "blue",pch = "+")
plot3d(samples[,1], samples[,2],pch=20,col='blue')

# ---- # CDF and PDF # ---- #

# generate normal copula with coef 0.8
# and sample some observations from it
coef_ <- 0.8
mycopula <- normalCopula(coef_,dim=3)
u <- rCopula(2000, mycopula)

# compute the density PDF: dCopula()
pdf_<-dCopula(u, mycopula)

# compute the cumulative CDF: pCopula()
cdf <- pCopula(u, mycopula)


# generate random sample obs-s 
# from multivariate distribution
cp <- claytonCopula(param = c(3.4), dim = 2)

# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
multivariate_dist <- mvdc(copula = cp,
                  margins = c("norm", "t"),
                  paramMargins = list(list(mean = 2, sd=3),list(df = 2)))

v <- rMvdc(2000, multivariate_dist)

# compute density
pdf_mvd <- dMvdc(v,multivariate_dist)

# compute cdf
cdf_mvd <- pMvdc(v,multivariate_dist)

# ----# Graphics # ---- #
# 3D plain scatterplot of the multivariate distribution
par(mfrow = c(1, 1))
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
persp(multivariate_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Density")
contour(multivariate_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")
persp(multivariate_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "CDF")
contour(multivariate_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")

plot3d(v[,1],v[,2], pdf_mvd,pch=20,col='blue',main="Density", xlab = "u1", ylab="u2", zlab="pMvdc")
plot3d(v[,1],v[,2], cdf_mvd,pch=20,col='blue',main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc")

#####################
#                   #
#      example      #
#                   #
#####################


mydata <- read.csv('mydata.csv',header=T)

plot(mydata$x,mydata$y)
hist(mydata$x)
hist(mydata$y)

# Estimate x  gamma distribution parameters and 
# visually compare simulated vs observed data

x_mean <- mean(mydata$x)
x_var <- var(mydata$x)
x_rate <- x_mean / x_var
x_shape <- ( (x_mean)^2 ) / x_var

# hist of mydatax
hist(mydata$x, breaks = 20, col = "green", density = 20)
# hist of gamma with rate and shape of x
hist(rgamma(nrow(mydata), rate = x_rate, shape = x_shape), 
     breaks = 20,col = "blue", add = T, density = 20, angle = -45)

# Estimate y gamma distribution parameters and visually compare simulated vs observed data
y_mean <- mean(mydata$y)
y_var <- var(mydata$y)
y_rate <- y_mean / y_var
y_shape <- ( (y_mean)^2 ) / y_var


hist(mydata$y, breaks = 20, col = "green", density = 20)
hist(rgamma(nrow(mydata), rate = y_rate, shape = y_shape), 
     breaks = 20, col = "blue", add = T, density = 20, angle = -45)

# measure of association
cor(mydata$x,mydata$y, method='kendall')

# BiCopSelect() to select copula
# acepts psedo-observations as parameter.
# pobs() returns a matrix
var_a <- pobs(mydata$x)
var_b <- pobs(mydata$y)
selectedCopula <- BiCopSelect(var_a,var_b,familyset=NA)
selectedCopula

# Fitting process
# Estimate copula parameters
cop_model <- claytonCopula(dim = 2)
m <- pobs(as.matrix(mydata[,2:3]))
fit <- fitCopula(cop_model, m, method = 'ml')
coef(fit)

# Check Kendall's tau value for the Clayton copula with theta = 1.48
tau(claytonCopula(param = 1.48))

# Goodness of fit test GOF

# clayton copula
gfc <- gofCopula(claytonCopula(dim = 2), as.matrix(mydata[,2:3]), N = 50)
gfc

# bivariate normal copula
gf <- gofCopula(normalCopula(dim = 2), as.matrix(mydata[,2:3]), N = 50)
gf


# Building bivariate distribution
# using the copula

my_dist <- mvdc(claytonCopula(param = 1.48, dim = 2), 
                margins = c("gamma","gamma"), 
                paramMargins = list(list(shape = x_shape, rate = x_rate), 
                                    list(shape = y_shape, rate = y_rate)))

# Generate random sample observations 
# from the multivariate distribution
v <- rMvdc(5000, my_dist)

# Compute the density
pdf_mvd <- dMvdc(v, my_dist)
# Compute the CDF
cdf_mvd <- pMvdc(v, my_dist)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 1))
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
persp(my_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Density")
contour(my_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")
persp(my_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "CDF")
contour(my_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")

# Now we can sample some observations from the estimated 
# joint distribution and check how they compare to the one observed.

sim <- rMvdc(1000, my_dist)

# Plot the data for a visual comparison
plot(mydata$x, mydata$y, main = 'Test dataset x and y', col = "blue")
points(sim[,1], sim[,2], col = 'red')
legend('bottomright', c('Observed', 'Simulated'), col = c('blue', 'red'), pch=21)

cor(mydata[,2:3], method = "kendall")
cor(sim, method = "kendall")

























