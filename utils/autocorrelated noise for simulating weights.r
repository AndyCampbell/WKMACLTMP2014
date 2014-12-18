


cv <- 0.1
rho <- 0.6

x=rnorm(100) * cv
plot(x,type="l")
for (i in 2:length(x)) x[i] <- rho * x[i-1] + sqrt(1-rho^2) * x[i]
lines(x,col="red")
acf(x)


cW <- rep(5,100)
plot(cW,type="l")
cW <- cW * exp(x)
lines(cW,col="red")
