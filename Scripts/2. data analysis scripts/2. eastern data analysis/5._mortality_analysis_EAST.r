#!/usr/bin/env Rscript
#writing mortality model in rjags because runjags wont install with 32bit JAGS on scc1. 
#clear environment, load packages.
rm(list=ls())
library(data.table)
library(runjags)
library(here)
#source('required_products_utilities/master_path_list.r')
master_path <- paste0(here(),'/master_path_list.r') #gets path to whatever directory you are already running in.
source(master_path)

#load clock functions.
tic = function() assign("timer", Sys.time(), envir=.GlobalEnv)
toc = function() print(Sys.time()-timer)

#specify filepaths
model.path      = remote_mortality_model_path.east
prediction.path = remote_mortality_prediction_path.east
summary.path    = remote_mortality_summary_path.east
filtered.path   = remote_mortality_filtered_path.east
range.path      = remote_mortality_range_path.east


#load data - subset 
#d <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/Product_2.rds')
d <- readRDS(remote_Product_2.E.path)

#throw an error if data doesn't exist. 
if(exists('d') == FALSE) stop('The data did not load. Does not run w/o data. duh.')

#do some quick data processing
#only considering AM and EM trees. Add this information
d$em    <- ifelse(d$MYCO_ASSO %in% c('ECM'),1,
                  ifelse(d$MYCO_ASSO %in% c('AM'),0,NA))

d$cn <- d$C.storage / d$N.storage

#grab complete.cases. 1911 observations. 
dat <- d[,.(mat, map, tot.15, REMPER, em, pH_H2O, cn, clay, death, DIA, LAT, LON, PLT_CN,PREV_PLT_CN, TRE_CN)]
dat <- dat[complete.cases(dat),]


#remove places where more than 50% of trees died- Removes 202 sites, 1709 remain. 
dat <- dat[,m.rate := sum(death, na.rm=T) / (length(death)),by=PLT_CN]
dat <- dat[m.rate < 0.5,]

#remove tree death caused by insects and disease. 
#Insect mortality represents 1793 of 4430 tree deaths observed in this data set. 
to.remove<- rbind(d[AGENTCD==10,], d[AGENTCD==20,])
dat<- dat[!(TRE_CN %in% to.remove$TRE_CN),]

#save filtered output for mapping sites later
saveRDS(dat,filtered.path)

#assign predictors to codes for model
y  <- dat$death
t  <- as.numeric(dat$REMPER) #census interval, t
x1 <- dat$DIA
x2 <- dat$mat
x3 <- dat$map
x4 <- dat$pH_H2O
x5 <- dat$cn
x6 <- dat$tot.15
x7 <- dat$em
x8 <- dat$clay

#get predictor means and ranges for predictions from model 
m.t         <- mean(t )
m.x1        <- mean(x1)
m.x2        <- mean(x2)
m.x3        <- mean(x3)
m.x4        <- mean(x4)
m.x5        <- mean(x5)
m.x6        <- mean(x6)
m.x7        <- mean(x7)
m.x8        <- mean(x8)

#ranges for parameters of interest other than em, which will be run as 0 or 1 for AM or EM. 
range.x1 <- seq(min(x1), max(x1), by=(max(x1) - min(x1))/100)
range.x2 <- seq(min(x2), max(x2), by=(max(x2) - min(x2))/100)
range.x3 <- seq(min(x3), max(x3), by=(max(x3) - min(x3))/100)
range.x4 <- seq(min(x4), max(x4), by=(max(x4) - min(x4))/100)
range.x5 <- seq(min(x5), max(x5), by=(max(x5) - min(x5))/100)
range.x6 <- seq(min(x6), max(x6), by=(max(x6) - min(x6))/100)
range.x7 <- seq(0,1,by=0.01)
range.x8 <- seq(min(x8), max(x8), by=(max(x8) - min(x8))/100)

#save data ranges for downstream plotting later
range <- data.frame(range.x1, range.x2, range.x3, range.x4, range.x5, range.x6, range.x7, range.x8)
saveRDS(range,range.path)

jags.model = "
model {
#model!
for (i in 1:N){
y[i] ~ dbern(m[i])
m[i] <- (1 - (1 - p[i])^t[i]) #iterate annual mort probability to whatever interval tree was remeasured on. 
p[i] <-  1 / (1 + exp(-z[i])) #inverse logit transform.
z[i] <- a +   exp(b * x1[i]) + exp(-c * x1[i]) + d*x2[i] + e*x3[i] + f*x4[i] + g*x5[i] + h*x6[i] + j*x7[i] + k*x6[i]*x7[i] + l*x8[i]
#y[i] ~ dbern(p[i] ^ t[i])                                #y[i] is the mortality probabilty, drawn with prob p that varies between 0 and 1. t[i] is the recensus interval.
#p[i] <- 1 / (1 + exp(-z[i]))                             #p[i] is the inverse  logit transform of z[i], converting our model of z to 0-1 probability interval. 
#z[i] is a continuous variable. More negative numbers translate to lower mortality probabilities.
#z[i] <- a +   exp(b * x1[i]) + exp(-c * x1[i]) + d*x2[i] + e*x3[i] + f*x4[i] + g*x5[i] + h*x6[i] + j*x7[i] + k*x6[i]*x7[i]
}

#predictions!
for (i in 1:N.pred){
z.x1[i]    <- a +   exp(b * range.x1[i]) + exp(-c * range.x1[i]) + d*m.x2 + e*m.x3 + f*m.x4 + g*m.x5 + h*m.x6 + j*m.x7 + k*m.x6*m.x7 + l*m.x8
z.x2[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*range.x2[i] + e*m.x3 + f*m.x4 + g*m.x5 + h*m.x6 + j*m.x7 + k*m.x6*m.x7 + l*m.x8
z.x3[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*range.x3[i] + f*m.x4 + g*m.x5 + h*m.x6 + j*m.x7 + k*m.x6*m.x7 + l*m.x8
z.x4[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*range.x4[i] + g*m.x5 + h*m.x6 + j*m.x7 + k*m.x6*m.x7 + l*m.x8
z.x5[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*m.x4 + g*range.x5[i] + h*m.x6 + j*m.x7 + k*m.x6*m.x7 + l*m.x8
z.x6[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*m.x4 + g*m.x5 + h*range.x6[i] + j*m.x7 + k*range.x6[i]*m.x7 + l*m.x8
z.x6.am[i] <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*m.x4 + g*m.x5 + h*range.x6[i] + j*0    + k*range.x6[i]*0 + l*m.x8
z.x6.em[i] <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*m.x4 + g*m.x5 + h*range.x6[i] + j*1    + k*range.x6[i]*1 + l*m.x8
z.x7[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*m.x4 + g*m.x5 + h*m.x6 + j*range.x7[i] + k*m.x6*range.x7[i] + l*m.x8
z.x8[i]    <- a +   exp(b * m.x1) + exp(-c * m.x1) + d*m.x2 + e*m.x3 + f*m.x4 + g*m.x5 + h*m.x6 + j*m.x7 + k*m.x6*m.x7 + l*range.x8[i]
}

#priors
a ~ dnorm(0, .0001)
b ~ dnorm(0, .0001) I(0, )
c ~ dnorm(0, .0001) I(0, )
d ~ dnorm(0, .0001)
e ~ dnorm(0, .0001)
f ~ dnorm(0, .0001)
g ~ dnorm(0, .0001)
h ~ dnorm(0, .0001)
j ~ dnorm(0, .0001)
k ~ dnorm(0, .0001)
l ~ dnorm(0, .0001)
}
"


data = list( y = y,    
             x1=x1, x2=x2,x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, x8=x8,
             m.x1=m.x1,m.x2=m.x2,m.x3=m.x3,m.x4=m.x4,m.x5=m.x5,m.x6=m.x6,m.x7=m.x7,m.x8=m.x8,
             range.x1=range.x1,range.x2=range.x2,range.x3=range.x3,range.x4=range.x4,range.x5=range.x5,range.x6=range.x6,range.x7=range.x7,range.x8=range.x8,
             N.pred=length(range.x1),
             N = length(y),
             t=t
) 


#a, b and c inits Chosen based on an initial run modeling mortality only as a function of the diamter predictors, w/o inits and 1000 coda samples.
#remaining init values based on a glm binomial fit, ignoring REMPER period. 
inits <- list()
inits[[1]] <- list(a=1.00, b = 0.010, c = 0.010, d = -0.41, e = 0.0004, f = -0.09, g = -0.006, h = -0.004, j = -0.25, k = 0.04, l = 0.002)
inits[[2]] <- list(a=1.01, b = 0.011, c = 0.011, d = -0.42, e = 0.0005, f = -0.10, g = -0.007, h = -0.005, j = -0.26, k = 0.05, l = 0.003)
inits[[3]] <- list(a=1.02, b = 0.012, c = 0.012, d = -0.43, e = 0.0006, f = -0.11, g = -0.008, h = -0.006, j = -0.27, k = 0.06, l = 0.004)



#begin runjags model
cat("Intializing, burning in and sampling runjags model! \n")
tic()
jags.out <- run.jags(jags.model,
                     data=data,
                     inits=inits,
                     n.chains=3,
                     monitor=c('a','b','c','d','e','f','g','h','j','k','l'),
                     method = 'rjparallel')
toc()

#save model output
saveRDS(jags.out, model.path)
jags.out <- readRDS(model.path)

#use extend.jags as in soil script to get prediction interval estimates.
cat("Getting prediction estimates!! \n")
tic()
prediction.out <- extend.jags(jags.out, add.monitor=c('z.x1','z.x2','z.x3','z.x4','z.x5','z.x6','z.x7','z.x6.am','z.x6.em','z.x8'))
toc()
saveRDS(prediction.out,prediction.path)

#summarize predictions and save- much smaller file size for plotting
#note inverse logit transform still has to happen. 
prediction.out <- readRDS(prediction.path)
cat("SUMMARIZING AND SAVING PREDICTION OUTPUT!!! \n")
tic()
out <- summary(prediction.out)
toc()
saveRDS(out,summary.path)
