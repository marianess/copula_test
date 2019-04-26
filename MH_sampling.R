library(copula)
## set seed so results reproducible
set.seed(42)
##pdfcalc(x) calculates the pdf of x for bivariate copula rho=p
p <- 0.9
mycopula <- normalCopula(p , dim = 2)

pdfcalc=function(x){
#pdfc = exp(-0.5*(x[1]^2-2*p*x[1]*x[2]+x[2]^2)/(1-p^2))
pdfc = dCopula(x, copula = mycopula)   
return(pdfc)
}

## setting the parameters
K=1000
x <- matrix(rep(0.5,2*K),nrow=K,ncol=2)
x[1,1] <- runif(1)#rnorm(1)
x[1,2] <- runif(1)#rnorm(1)

## sampling
for (i in 2:K){
    currentx = x[i-1,]
    proposedx = currentx + rnorm(2, mean=0,sd=0.20)
    ratio= pdfcalc(proposedx)/pdfcalc(currentx)
    if(runif(1) < ratio){
        x[i,] = proposedx # accept move with probabily min(1,A)
    }else{
        x[i,] = currentx  # otherwise "reject" move, and stay where we are
    }
    
}

plot(x)
plot(as.ts(x[,1]))
hist(x[,1])
print(cov(x))
print(c("means",mean(x[,1]),mean(x[,2])))


