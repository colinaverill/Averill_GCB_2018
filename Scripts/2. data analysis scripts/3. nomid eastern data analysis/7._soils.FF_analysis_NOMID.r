#building bayesian analysis of forest floor soil C. Linear model accounting for increasing variance in C as a function of its observed value. 

#######################################################
##########     Load and filter data     ###############
#######################################################
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

#load data - local
#d<- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/Product_1.rds')
#load data - remote
d<- readRDS(remote_Product_1.E_nomid.path)


#get FF C as total - mineral.
d$ff.C.storage <- d$C.storage - d$m.C.storage
d$ff.N.storage <- d$N.storage - d$m.N.storage

hist(d$ff.C.storage)

#bump up values with smallest C and N storage numbers for downstream log transform.
d$ff.C.storage <- d$ff.C.storage + min(d[d$ff.C.storage > 0]$C.storage, na.rm = T)
d$ff.N.storage <- d$ff.N.storage + min(d[d$ff.C.storage > 0]$N.storage, na.rm = T)

#mod <- lm(log(m.C.storage) ~ mat + map + pH_H2O + clay + tot.15 + relEM, data = d)
#z <- d[,.(mat, map, pH_H2O, clay, tot.15, ff.N.storage)]
#z <- colMeans(z, na.rm = T)
#newdat <- data.frame(seq(0,1, by = 0.01)); colnames(newdat) <- 'relEM'
#for(i in 1:length(z)){
#  newdat[,i + 1] <- z[i]
#  colnames(newdat)[i + 1] <- names(z[i])
#}
#test <- predict(mod, newdata = newdat)
#plot(exp(test) ~ newdat$relEM)

#specify file paths
model.path      = remote_soils.ff_model_path.nomid
prediction.path = remote_soils.ff_prediction_path.nomid
summary.path    = remote_soils.ff_summary_path.nomid
filtered.path   = remote_soils.ff_filtered_path.nomid
range.path      = remote_soils.ff_range_path.nomid


#Convert C and N storage from kg / m2 to g /m2 for easier interpretation of log transformation.
d[,ff.C.storage := ff.C.storage * 1000]
d[,ff.N.storage := ff.N.storage * 1000]

#gather complete cases- everything that is used to model C.storage.
#complete cases drops us from 3621 observations to 2880. 
dat <-   d[,.(mat, map, pH_H2O, tot.15, relEM, clay, ff.N.storage, ff.C.storage)]
dat <- dat[complete.cases(dat),]
saveRDS(dat,filtered.path)


#######################################################
##########           Build model        ###############
#######################################################

#assign predictors to codes for model
y  <- dat$ff.C.storage
x1 <- dat$mat
x2 <- dat$map
x3 <- dat$pH_H2O
x4 <- dat$tot.15
x5 <- dat$relEM
x6 <- log(dat$ff.N.storage)
x7 <- dat$clay

#generate ranges and means for predictions
x1.range <- seq(min(x1), max(x1), by = ((max(x1) - min(x1)) / 100))
x2.range <- seq(min(x2), max(x2), by = ((max(x2) - min(x2)) / 100))
x3.range <- seq(min(x3), max(x3), by = ((max(x3) - min(x3)) / 100))
x4.range <- seq(min(x4), max(x4), by = ((max(x4) - min(x4)) / 100))
x5.range <- seq(0,1, by=0.01)
x7.range <- seq(min(x7), max(x7), by = ((max(x7) - min(x7)) / 100))

#save ranges for downstream plotting
range <- data.frame(x1.range, x2.range,x3.range,x4.range,x5.range,x7.range)
saveRDS(range, range.path)

m.x1 <- mean(x1)
m.x2 <- mean(x2)
m.x3 <- mean(x3)
m.x4 <- mean(x4)
m.x5 <- mean(x5)
m.x6 <- mean(x6)
m.x7 <- mean(x7)

#add in predictions for specific levels of EM and N.
em.range <- seq(0.7928644,0.1849294, -(0.7928644 - 0.1849294)/100)
n.range <- x4.range

#add in predictions for effect of EM at 3 levels of N.
n.low  <-  min(dat$tot.15)
n.mean <- mean(dat$tot.15)
n.high <-  18

#write JAGS model.
jags.model = "
model {
#model
for (i in 1:N){
#y[i] is the soil C stock, with a variance scaled to the observed value of tau. Bigger value, less precision. 
y[i] ~ dlnorm(y.hat[i], tau)

#y.hat contains the model.
y.hat[i] <- a + b*x1[i] + c*x2[i] + d*x3[i] + e*x4[i] + f*x5[i] + g*x6[i] + h*x5[i]*x6[i] + j*x4[i]*x5[i] + k*x7[i]
}

#predictions
for (i in 1:N.pred){
y.x1[i]    <- a + b*x1.range[i] + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x5*m.x6 + j*m.x4*m.x5 + k*m.x7
y.x2[i]    <- a + b*m.x1 + c*x2.range[i] + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x5*m.x6 + j*m.x4*m.x5 + k*m.x7
y.x3[i]    <- a + b*m.x1 + c*m.x2 + d*x3.range[i] + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x5*m.x6 + j*m.x4*m.x5 + k*m.x7
y.x4[i]    <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*x4.range[i] + f*m.x5 + g*m.x6 + h*m.x5*m.x6 + j*x4.range[i]*m.x5 + k*m.x7
y.x5[i]    <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*x5.range[i] + g*m.x6 + h*x5.range[i]*m.x6 + j*m.x4*x5.range[i] + k*m.x7
y.x7[i]    <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x5*m.x6 + j*m.x4*m.x5 + k*x7.range[i]
y.x4.em[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*x4.range[i] + f*1 + g*m.x6 + h*1*m.x6 + j*x4.range[i]*1 + k*m.x7
y.x4.am[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*x4.range[i] + f*0 + g*m.x6 + h*0*m.x6 + j*x4.range[i]*0 + k*m.x7
y.npred[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*n.range[i] + f*em.range[i] + g*m.x6 + h*em.range[i]*m.x6 + j*n.range[i]*em.range[i] + k*m.x7
y.nlow[i]  <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*n.low  + f*x5.range[i] + g*m.x6 + h*x5.range[i]*m.x6 + j*n.low*x5.range[i] + k*m.x7
y.nmean[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*n.mean + f*x5.range[i] + g*m.x6 + h*x5.range[i]*m.x6 + j*n.mean*x5.range[i] + k*m.x7
y.nhigh[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*n.high + f*x5.range[i] + g*m.x6 + h*x5.range[i]*m.x6 + j*n.high*x5.range[i] + k*m.x7
}


#priors - flat, uninformative
a ~ dnorm(0, .0001)
b ~ dnorm(0, .0001)
c ~ dnorm(0, .0001)
d ~ dnorm(0, .0001)
e ~ dnorm(0, .0001)
f ~ dnorm(0, .0001)
g ~ dnorm(0, .0001)
h ~ dnorm(0, .0001)
j ~ dnorm(0, .0001)
k ~ dnorm(0, .0001)
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)
}
"

data = list( y = y,    
             x1=x1, x2=x2,x3=x3, x4=x4, x5=x5, x6=x6, x7=x7,
             x1.range=x1.range, x2.range=x2.range, x3.range=x3.range, x4.range=x4.range, x5.range=x5.range, x7.range=x7.range, n.range=n.range, em.range=em.range,
             m.x1=m.x1, m.x2=m.x2, m.x3=m.x3, m.x4=m.x4, m.x5=m.x5, m.x6=m.x6, m.x7=m.x7, n.low=n.low, n.mean=n.mean, n.high=n.high,
             N = length(y), N.pred = length(x1.range)
) 

#######################################################
##########           Run a model        ###############
#######################################################
tic()
jags.out <- run.jags(jags.model,
                     data=data,
                     n.chains=3,
                     monitor=c('a','b','c','d','e','f','g','h','j','k','tau'),
                     method = 'rjparallel',
                     sample = 50000)
#jags.out <- extend.jags(jags.out, sample = 20000)
toc()
#save output
saveRDS(jags.out, model.path)

#use extend.jags to get samples from converged model for each prediction range.
tic()
prediction.out <- extend.jags(jags.out, add.monitor=c('y.x1','y.x2','y.x3','y.x4','y.x5','y.x7','y.x4.em','y.x4.am','y.npred','y.nlow','y.nmean','y.nhigh'))
toc()
#save output
saveRDS(prediction.out, prediction.path)

#summarize predictions, and save. Much smaller file to deal with for generating figures. 
prediction.out <- readRDS(prediction.path)
tic()
out <- summary(prediction.out)
toc()
saveRDS(out, summary.path)