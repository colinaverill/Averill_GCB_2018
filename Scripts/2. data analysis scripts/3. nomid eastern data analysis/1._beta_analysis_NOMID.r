#modeling the relative abundance of EM vs. AM trees. 
#using beta regression, correct 0/1 values
#clear envrionment, load packages
rm(list=ls())
library(data.table)
library(runjags)
library(here)
#source('required_products_utilities/master_path_list.r')
master_path <- paste0(here(),'/master_path_list.r') #gets path to whatever directory you are already running in.
source(master_path)

#load data.
#d <- readRDS(Product_1.path)
d <- data.table(readRDS(remote_Product_1.E_nomid.path))

#timer functions
tic = function() assign("timer", Sys.time(), envir=.GlobalEnv)
toc = function() print(Sys.time()-timer)

#these are the remote file path destinations
model.path      = remote_beta_model_path.nomid
prediction.path = remote_beta_prediction_path.nomid
summary.path    = remote_beta_summary_path.nomid
filtered.path   = remote_beta_filtered_path.nomid
range.path      = remote_beta_range_path.nomid

#grab complete cases - 2,384 sites. 
dat <-   d[, .(mat, map, tot.15, cn, pH_H2O, cn, clay, relEM, BASAL,PLT_CN,latitude,longitude)]
dat <- dat[complete.cases(dat),]

#transform relEM from [0,1] to (0,1)
dat[,relEM.beta := (relEM * (nrow(dat) - 1) + 0.5) / nrow(dat)]

#save filtered data product
saveRDS(dat, filtered.path)

#name variables
y  <- dat$relEM.beta
x1 <- dat$mat
x2 <- dat$map
x3 <- dat$cn
x4 <- dat$pH_H2O
x5 <- dat$tot.15
x6 <- dat$clay

#get predictor means and ranges for predictions from model 
m.x1        <- mean(x1)
m.x2        <- mean(x2)
m.x3        <- mean(x3)
m.x4        <- mean(x4)
m.x5        <- mean(x5)
m.x6        <- mean(x6)

#ranges for parameters of interest.
range.x1 <- seq(min(x1), max(x1), by=(max(x1) - min(x1))/100)
range.x2 <- seq(min(x2), max(x2), by=(max(x2) - min(x2))/100)
range.x3 <- seq(min(x3), max(x3), by=(max(x3) - min(x3))/100)
range.x4 <- seq(min(x4), max(x4), by=(max(x4) - min(x4))/100)
range.x5 <- seq(min(x5), max(x5), by=(max(x5) - min(x5))/100)
range.x6 <- seq(min(x6), max(x6), by=(max(x6) - min(x6))/100)

#save ranges for downstream plotting. 
range.out <- data.frame(range.x1,range.x2,range.x3,range.x4,range.x5,range.x6)
colnames(range.out) <- c('mat.range','map.range','cn.range','ph.range','ndep.range','clay.range')
saveRDS(range.out, range.path)

jags.model = "
model{
# priors
a0 ~ dnorm(0, .001)
a1 ~ dnorm(0, .001)
a2 ~ dnorm(0, .001)
a3 ~ dnorm(0, .001)
a4 ~ dnorm(0, .001)
a5 ~ dnorm(0, .001)
a6 ~ dnorm(0, .001)
#t0 ~ dnorm(0, .01) 
#tau <- exp(t0)
tau ~ dgamma(.1,.1)


# likelihood for mu and tau- This is a beta regression on the continuous values between 0 and 1. mu2 is p, and tau is q. Both are shape parameters. 
for (i in 1:N){
y[i] ~ dbeta(p[i], q[i])
p[i] <- mu[i] * tau
q[i] <- (1 - mu[i]) * tau
logit(mu[i]) <- a0 + a1 * x1[i] + a2 * x2[i] + a3 * x3[i] + a4 * x4[i] + a5 * x5[i] + a6 * x6[i]
}

#predictions!
for (i in 1:N.pred){
z.x1[i]    <- a0 + a1 * range.x1[i] + a2 * m.x2 + a3 * m.x3 + a4 * m.x4 + a5 * m.x5 + a6 * m.x6
z.x2[i]    <- a0 + a1 * m.x1 + a2 * range.x2[i] + a3 * m.x3 + a4 * m.x4 + a5 * m.x5 + a6 * m.x6
z.x3[i]    <- a0 + a1 * m.x1 + a2 * m.x2 + a3 * range.x3[i] + a4 * m.x4 + a5 * m.x5 + a6 * m.x6
z.x4[i]    <- a0 + a1 * m.x1 + a2 * m.x2 + a3 * m.x3 + a4 * range.x4[i] + a5 * m.x5 + a6 * m.x6
z.x5[i]    <- a0 + a1 * m.x1 + a2 * m.x2 + a3 * m.x3 + a4 * m.x4 + a5 * range.x5[i] + a6 * m.x6
z.x6[i]    <- a0 + a1 * m.x1 + a2 * m.x2 + a3 * m.x3 + a4 * m.x4 + a5 * m.x5 + a6 * range.x6[i]
}
}  
"

jd <- list(x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, x6=x6,
           y=y,
           N=length(y),N.pred = length(range.x1),
           m.x1=m.x1,m.x2=m.x2,m.x3=m.x3,m.x4=m.x4,m.x5=m.x5, m.x6=m.x6,
           range.x1=range.x1,range.x2=range.x2,range.x3=range.x3,range.x4=range.x4,range.x5=range.x5, range.x6=range.x6)

#inits. Chosen based on output of betareg function from the betareg package. 
inits <- list()
inits[[1]] <- list(a0=2.89, a1 = -0.083, a2 = 0.0009, a3 = 0.019, a4 = -0.31, a5 = -0.33, a6= 0.01)
inits[[2]] <- list(a0=2.90, a1 = -0.084, a2 = 0.0010, a3 = 0.020, a4 = -0.32, a5 = -0.34, a6= 0.00)
inits[[3]] <- list(a0=2.88, a1 = -0.082, a2 = 0.0008, a3 = 0.018, a4 = -0.30, a5 = -0.32, a6=-0.01)


tic()
jags.out <- run.jags(jags.model,
                     data=jd,
                     n.chains=3,
                     inits=inits,
                     monitor=c('a0','a1','a2','a3','a4','a5','a6','deviance'),
                     method = 'rjparallel')
toc()
#save model output
saveRDS(jags.out, model.path)
jags.out <- readRDS(model.path)

#use extend.jags as in soil script to get prediction interval estimates.
tic()
prediction.out <- extend.jags(jags.out, add.monitor=c('z.x1','z.x2','z.x3','z.x4','z.x5','z.x6'))
toc()
saveRDS(prediction.out,prediction.path)

#summarize predictions and save- much smaller file size for plotting
#note inverse logit transform still has to happen. 
prediction.out <- readRDS(prediction.path)
tic()
out <- summary(prediction.out)
toc()
saveRDS(out,summary.path)