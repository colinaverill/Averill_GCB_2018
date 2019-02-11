#multiple stable combinations of parameters for growth.gaus.
#subset to some chains that work, validate.
#move forward with these chains for downstream predictions and plotting.

rm(list=ls())
source('required_products_utilities/master_path_list.r')
library(runjags)

#Update: putting intercept back in. Constraining gaussian parameters to be (+) because they don't make sense otherwise.
mod <- readRDS(local_growth.gaus_model_path.east)
mod <- readRDS(local_growth.gaus_model_path)
summary(mod)

#convert it to a list you can trim
z <- coda::as.mcmc.list(mod)
z <- window(z, start = 207001, end = 211000, thin = 1)
test <- z
for(i in 1:length(z)){
  test[[i]] <- z[[i]][6001:nrow(z[[i]]),]
}
test <- coda::as.mcmc.list(test)

#pick the happy chains.
y <- coda::mcmc.list(z[[12]],z[[10]],z[[7]],z[[6]],z[[4]],z[[3]],z[[2]])
y <- coda::mcmc.list(z[[10]],z[[3]],z[[2]])
plot(y)

#trim off the lead a row that is singular (intercept set to zero)
#y <- lapply(y, function(x)x[,c(2:ncol(x))])
#y <- coda::as.mcmc.list(y)
#get gelman diagnostics and plot.
coda::gelman.diag(y)
plot(y)

gaus.fun <- function(x, a, b, c, d) {a * exp(-1*((x-b)^2)/2*c^2) + d*x}
x1 <- seq(-4,4, by = 0.1)
x2 <- gaus.fun(x1, 5.6,7,7.9)
x2 <- gaus.fun(x1, 1,1,1,.2)
plot(x2 ~ x1)


#try with a negative exponential
x1 <- seq(-4,4, by = 0.1)
x2 <- x1*1 + -0.2 * x1^2
plot(x2 ~ x1)

1/(1 + exp(-x1))
