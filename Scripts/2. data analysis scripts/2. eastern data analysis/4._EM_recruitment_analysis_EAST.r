#This script runs a zero-inflated poisson model to fit adult tree recruitment from FIA scripts (no saplings).
#This script fits EM recruitment for all plots.
#predictions monitor both the probability of recruitment happening at all, and if it does how many, as separate things.
#clear environment, load packages.
rm(list=ls())
library(data.table)
library(runjags)
library(here)
#source('required_products_utilities/master_path_list.r')
master_path <- paste0(here(),'/master_path_list.r') #gets path to whatever directory you are already running in.
source(master_path)

#timer functions
tic = function() assign("timer", Sys.time(), envir=.GlobalEnv)
toc = function() print(Sys.time()-timer)

#Save paths.
model.path      = remote_EM_recruitment_model_path.east
prediction.path = remote_EM_recruitment_prediction_path.east
summary.path    = remote_EM_recruitment_summary_path.east
filtered.path   = remote_EM_recruitment_filtered_path.east
range.path      = remote_EM_recruitment_range_path.east

#load data for analysis
#d <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/Product_3.rds')
d <- readRDS(remote_Product_3.E.path)
d <- data.table(d)

#grab relEM and relAM
d[, relEM := BASAL.ECM/BASAL]
d[, relAM := BASAL.AM /BASAL]

#calculate soil CN
d[, cn := C.storage / N.storage]

#grab relevant data. 1878 unique sites. 
dat <- d[,.(recruit.adult.em, pH_H2O, cn, mat, map, relEM, relAM, tot.15, clay, REMPER,PLT_CN, PREV_PLT_CN, BASAL)]
dat <- dat[complete.cases(dat),]
#save for later
saveRDS(dat, filtered.path)

#grab variables for model
y  <- dat$recruit.adult.em
t  <- dat$REMPER
x1 <- dat$pH_H2O
x2 <- dat$cn
x3 <- dat$mat
x4 <- dat$map
x5 <- dat$tot.15
x6 <- dat$relEM
x7 <- dat$clay
x8 <- log(dat$BASAL)

#grab means of predictors for prediction.
x1.m <- mean(x1)
x2.m <- mean(x2)
x3.m <- mean(x3)
x4.m <- mean(x4)
x5.m <- mean(x5)
x6.m <- mean(x6)
x7.m <- mean(x7)
x8.m <- mean(x8)
t.m  <- mean(t)

#get sd of y and each predictor - necesssary for calculating beta factors.
 y.sd <- sd(y )
x1.sd <- sd(x1)
x2.sd <- sd(x2)
x3.sd <- sd(x3)
x4.sd <- sd(x4)
x5.sd <- sd(x5)
x7.sd <- sd(x7)

#ranges of predictors for prediction. 
x1.r <- seq(min(x1), max(x1), by = ((max(x1) - min(x1))/ 100))
x2.r <- seq(min(x2), max(x2), by = ((max(x2) - min(x2))/ 100))
x3.r <- seq(min(x3), max(x3), by = ((max(x3) - min(x3))/ 100))
x4.r <- seq(min(x4), max(x4), by = ((max(x4) - min(x4))/ 100))
x5.r <- seq(min(x5), max(x5), by = ((max(x5) - min(x5))/ 100))
x6.r <- seq(min(x6), max(x6), by = ((max(x6) - min(x6))/ 100))
x7.r <- seq(min(x7), max(x7), by = ((max(x7) - min(x7))/ 100))

#save ranges for downstream plotting.
range <- data.frame(x1.r, x2.r, x3.r, x4.r, x5.r, x6.r, x7.r)
saveRDS(range, range.path)


j.model ="
model{
for (i in 1:N){
#this fits the blended model to your data. 
y[i] ~ dpois(m[i]*t[i])

#This blends the poisson and zero inflation models
m[i] <- mu[i]*x[i] + 0.00001

#this is the bernoulli outcome of the zero inflation
x[i] ~ dbern(pro[i])

#this logit transforms the theta model for 0-1 probability of zero inflation
logit(pro[i]) <- theta[i]

#mu[i] is predictors for poisson model- log link function.
log(mu[i]) <- int.p + x1[i]*b1.p + x2[i]*b2.p + x3[i]*b3.p + x4[i]*b4.p + x5[i]*b5.p + x6[i]*b6.p + x7[i]*b7.p# + b8.p*exp(-1*((pow((x8[i] - b9.p),2)) / (2*pow(b10.p,2))))

#theta[i] is predictors for zero inflation
theta[i]   <- int.z + x1[i]*b1.z + x2[i]*b2.z + x3[i]*b3.z + x4[i]*b4.z + x5[i]*b5.z + x6[i]*b6.z + x7[i]*b7.z# + b8.z*exp(-1*((pow((x8[i] - b9.z),2)) / (2*pow(b10.z,2))))
}

#predictions
#This gets the range of pred.m blended model values, before they get sent through the poisson distribution, and corrected for time. 
#This also uses the probability of recruitment happening at all, pre-bernoulli transform, to avoid non-linear breaks in 95% CI.
for (i in 1:N.pred){
  ph.pred[i] <- pred1.mu[i]*pred1.pro[i] + 0.00001
  cn.pred[i] <- pred2.mu[i]*pred2.pro[i] + 0.00001
 mat.pred[i] <- pred3.mu[i]*pred3.pro[i] + 0.00001
 map.pred[i] <- pred4.mu[i]*pred4.pro[i] + 0.00001
ndep.pred[i] <- pred5.mu[i]*pred5.pro[i] + 0.00001
clay.pred[i] <- pred7.mu[i]*pred7.pro[i] + 0.00001

#probabilty of recruiting at all.
  ph.pro.pred[i] <- pred1.pro[i] + 0.00001
  cn.pro.pred[i] <- pred2.pro[i] + 0.00001
 mat.pro.pred[i] <- pred3.pro[i] + 0.00001
 map.pro.pred[i] <- pred4.pro[i] + 0.00001
ndep.pro.pred[i] <- pred5.pro[i] + 0.00001
clay.pro.pred[i] <- pred7.pro[i] + 0.00001

#number of recruits if you do recruit. 
  ph.mu.pred[i] <- pred1.mu[i]
  cn.mu.pred[i] <- pred2.mu[i]
 mat.mu.pred[i] <- pred3.mu[i]
 map.mu.pred[i] <- pred4.mu[i]
ndep.mu.pred[i] <- pred5.mu[i]
clay.mu.pred[i] <- pred7.mu[i]

pred1.x[i] ~ dbern(pred1.pro[i])
pred2.x[i] ~ dbern(pred2.pro[i])
pred3.x[i] ~ dbern(pred3.pro[i])
pred4.x[i] ~ dbern(pred4.pro[i])
pred5.x[i] ~ dbern(pred5.pro[i])
pred7.x[i] ~ dbern(pred7.pro[i])

logit(pred1.pro[i]) <- pred1.theta[i]
logit(pred2.pro[i]) <- pred2.theta[i]
logit(pred3.pro[i]) <- pred3.theta[i]
logit(pred4.pro[i]) <- pred4.theta[i]
logit(pred5.pro[i]) <- pred5.theta[i]
logit(pred7.pro[i]) <- pred7.theta[i]

log(pred1.mu[i]) <- int.p + x1.r[i]*b1.p + x2.m*b2.p + x3.m*b3.p + x4.m*b4.p + x5.m*b5.p + x6.m*b6.p + x7.m*b7.p# + b8.p*exp(-1*((pow((x8.m - b9.p),2)) / (2*pow(b10.p,2))))
log(pred2.mu[i]) <- int.p + x1.m*b1.p + x2.r[i]*b2.p + x3.m*b3.p + x4.m*b4.p + x5.m*b5.p + x6.m*b6.p + x7.m*b7.p# + b8.p*exp(-1*((pow((x8.m - b9.p),2)) / (2*pow(b10.p,2))))
log(pred3.mu[i]) <- int.p + x1.m*b1.p + x2.m*b2.p + x3.r[i]*b3.p + x4.m*b4.p + x5.m*b5.p + x6.m*b6.p + x7.m*b7.p# + b8.p*exp(-1*((pow((x8.m - b9.p),2)) / (2*pow(b10.p,2))))
log(pred4.mu[i]) <- int.p + x1.m*b1.p + x2.m*b2.p + x3.m*b3.p + x4.r[i]*b4.p + x5.m*b5.p + x6.m*b6.p + x7.m*b7.p# + b8.p*exp(-1*((pow((x8.m - b9.p),2)) / (2*pow(b10.p,2))))
log(pred5.mu[i]) <- int.p + x1.m*b1.p + x2.m*b2.p + x3.m*b3.p + x4.m*b4.p + x5.r[i]*b5.p + x6.m*b6.p + x7.m*b7.p# + b8.p*exp(-1*((pow((x8.m - b9.p),2)) / (2*pow(b10.p,2))))
log(pred7.mu[i]) <- int.p + x1.m*b1.p + x2.m*b2.p + x3.m*b3.p + x4.m*b4.p + x5.m*b5.p + x6.m*b6.p + x7.r[i]*b7.p# + b8.p*exp(-1*((pow((x8.m - b9.p),2)) / (2*pow(b10.p,2))))

pred1.theta[i]   <- int.z + x1.r[i]*b1.z + x2.m*b2.z + x3.m*b3.z + x4.m*b4.z + x5.m*b5.z + x6.m*b6.z + x7.m*b7.z# + b8.z*exp(-1*((pow((x8.m - b9.z),2)) / (2*pow(b10.z,2))))
pred2.theta[i]   <- int.z + x1.m*b1.z + x2.r[i]*b2.z + x3.m*b3.z + x4.m*b4.z + x5.m*b5.z + x6.m*b6.z + x7.m*b7.z# + b8.z*exp(-1*((pow((x8.m - b9.z),2)) / (2*pow(b10.z,2))))
pred3.theta[i]   <- int.z + x1.m*b1.z + x2.m*b2.z + x3.r[i]*b3.z + x4.m*b4.z + x5.m*b5.z + x6.m*b6.z + x7.m*b7.z# + b8.z*exp(-1*((pow((x8.m - b9.z),2)) / (2*pow(b10.z,2))))
pred4.theta[i]   <- int.z + x1.m*b1.z + x2.m*b2.z + x3.m*b3.z + x4.r[i]*b4.z + x5.m*b5.z + x6.m*b6.z + x7.m*b7.z# + b8.z*exp(-1*((pow((x8.m - b9.z),2)) / (2*pow(b10.z,2))))
pred5.theta[i]   <- int.z + x1.m*b1.z + x2.m*b2.z + x3.m*b3.z + x4.m*b4.z + x5.r[i]*b5.z + x6.m*b6.z + x7.m*b7.z# + b8.z*exp(-1*((pow((x8.m - b9.z),2)) / (2*pow(b10.z,2))))
pred7.theta[i]   <- int.z + x1.m*b1.z + x2.m*b2.z + x3.m*b3.z + x4.m*b4.z + x5.m*b5.z + x6.m*b6.z + x7.r[i]*b7.z# + b8.z*exp(-1*((pow((x8.m - b9.z),2)) / (2*pow(b10.z,2))))
}

#priors
b1.p ~ dnorm(0, .0001)
b2.p ~ dnorm(0, .0001)
b3.p ~ dnorm(0, .0001)
b4.p ~ dnorm(0, .0001)
b5.p ~ dnorm(0, .0001)
b6.p ~ dnorm(0, .0001)
b7.p ~ dnorm(0, .0001)
#b8.p ~ dnorm(0, .0001)
#b9.p ~ dnorm(0, .0001)
#b10.p ~ dnorm(0, .0001)
b1.z ~ dnorm(0, .0001)
b2.z ~ dnorm(0, .0001)
b3.z ~ dnorm(0, .0001)
b4.z ~ dnorm(0, .0001)
b5.z ~ dnorm(0, .0001)
b6.z ~ dnorm(0, .0001)
b7.z ~ dnorm(0, .0001)
#b8.z ~ dnorm(0, .0001)
#b9.z ~ dnorm(0, .0001)
#b10.z ~ dnorm(0, .0001)
int.p ~ dnorm(0, .0001)
int.z ~ dnorm(0, .0001)
}
"

#data object for jags.model.
data=list(y=y,t=t,
          x1=x1, x2=x2, x3=x3, x4=x4, x5=x5,x6=x6, x7=x7,# x8=x8,
          x1.m=x1.m, x2.m=x2.m, x3.m=x3.m, x4.m=x4.m, x5.m=x5.m, x6.m=x6.m, x7.m=x7.m,# x8.m=x8.m,
          x1.r=x1.r, x2.r=x2.r, x3.r=x3.r, x4.r=x4.r, x5.r=x5.r, x6.r=x6.r, x7.r=x7.r,
          N = length(y), N.pred = length(x5.r))

#run jags model
tic()
jags.out <- autorun.jags(j.model,
                     data=data,
                     n.chains=3,
                     adapt = 4000,
#                     inits=inits,
                     monitor=c('int.p','b1.p','b2.p','b3.p','b4.p','b5.p','b6.p','b7.p',#'b8.p','b9.p','b10.p',
                               'int.z','b1.z','b2.z','b3.z','b4.z','b5.z','b6.z','b7.z',#'b8.z','b9.z','b10.z',
                               'deviance'),
                     method = 'rjparallel')
toc()
#save model output.
saveRDS(jags.out, model.path)

#use extend.jags to get prediction interval estimates. 
tic()
prediction.out <- extend.jags(jags.out, add.monitor=c('ph.pred','cn.pred','mat.pred','map.pred','ndep.pred','clay.pred',
                                                      'ph.pro.pred','cn.pro.pred','mat.pro.pred','map.pro.pred','ndep.pro.pred','clay.pro.pred',
                                                      'ph.mu.pred','cn.mu.pred','mat.mu.pred','map.mu.pred','ndep.mu.pred','clay.mu.pred'))
toc()
saveRDS(prediction.out,prediction.path)

#summarize prediction output
prediction.out <- readRDS(prediction.path)
tic()
out <- summary(prediction.out)
toc()
saveRDS(out,summary.path)
