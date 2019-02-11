#modeling growth of surviving trees.
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
model.path      = remote_growth.gaus_model_path
prediction.path = remote_growth.gaus_prediction_path
summary.path    = remote_growth.gaus_summary_path
filtered.path   = remote_growth.gaus_filtered_path
range.path      = remote_growth.gaus_range_path


#load data - local path
#d<- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/Product_3.rds')
#load data - remote path
d <- readRDS(remote_Product_3.path)

#calculate growth, basal.m2, C:N, exclude sites w/ C:N > 90.
d[,growth   := GROWTH.total / (REMPER * area.m2)]
d[,basal.m2 := BASAL / area.m2]
d[,cn       := C.storage / N.storage]
d = d[!cn>90,]

#remove 80 observations of shrinking trees at the plot level
#0s in basal.m2 reflect a few plots where all the trees died. 
d= d[growth   > 0,]
d= d[basal.m2 > 0,]

#calculate the relative abundance EM trees among trees that survived the remeasurement interval
d[,relEM:= BASAL.ECM / BASAL]

#losing 467 sites on complete cases. Is this just poor raster extraction? 
dat <- d[,.(growth, basal.m2, mat, map, tot.15, cn, pH_H2O, relEM, clay, latitude, longitude, PLT_CN,PREV_PLT_CN)]
dat <- dat[complete.cases(dat),]

#save filtered data for downstream mapping. 
saveRDS(dat, filtered.path)


#######################################################
##########        Build JAGS model      ###############
#######################################################

#assign predictors to codes for model
y  <- dat$growth
x1 <- log(dat$basal.m2)
x2 <- dat$mat
x3 <- dat$map
x4 <- dat$cn
x5 <- dat$pH_H2O
x6 <- dat$relEM
x7 <- dat$tot.15
x8 <- dat$clay

#generate ranges and means for predictions
x1.range <- seq(min(x1), max(x1), by = ((max(x1) - min(x1)) / 100))
x2.range <- seq(min(x2), max(x2), by = ((max(x2) - min(x2)) / 100))
x3.range <- seq(min(x3), max(x3), by = ((max(x3) - min(x3)) / 100))
x4.range <- seq(min(x4), max(x4), by = ((max(x4) - min(x4)) / 100))
x5.range <- seq(min(x5), max(x5), by = ((max(x5) - min(x5)) / 100))
x6.range <- seq(0, 1, by= 0.01)
x7.range <- seq(min(x7), max(x7), by = ((max(x7) - min(x7)) / 100))
x8.range <- seq(min(x8), max(x8), by = ((max(x8) - min(x8)) / 100))

#save ranges for downstream plotting.
range <- data.frame(x1.range,x2.range,x3.range,x4.range,x5.range,x6.range,x7.range,x8.range)
saveRDS(range, range.path)

m.x1 <- mean(x1)
m.x2 <- mean(x2)
m.x3 <- mean(x3)
m.x4 <- mean(x4)
m.x5 <- mean(x5)
m.x6 <- mean(x6)
m.x7 <- mean(x7)
m.x8 <- mean(x8)

#write JAGS model.
jags.model = "
model {
#model
for (i in 1:N){
y[i] ~ dlnorm(y.hat[i], tau) #y[i] is the growth value prediction, with precision tau 
#y.hat contains the model of the prediction.
y.hat[i] <- a + b*x1[i] + c*x2[i] + d*x3[i] + e*x4[i] + f*x5[i] + g*x6[i] + h*x7[i] + j*x6[i]*x7[i] + k*x8[i]  + g1 * exp(-1 * (pow((x7[i] - g2),2) / (2 * pow(g3, 2))))
}

#predictions
for (i in 1:N.pred){
y.x1[i] <- a + b*x1.range[i] + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x7 + j*m.x6*m.x7 + k*m.x8 + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
y.x2[i] <- a + b*m.x1 + c*x2.range[i] + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x7 + j*m.x6*m.x7 + k*m.x8 + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
y.x3[i] <- a + b*m.x1 + c*m.x2 + d*x3.range[i] + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x7 + j*m.x6*m.x7 + k*m.x8 + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
y.x4[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*x4.range[i] + f*m.x5 + g*m.x6 + h*m.x7 + j*m.x6*m.x7 + k*m.x8 + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
y.x5[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*x5.range[i] + g*m.x6 + h*m.x7 + j*m.x6*m.x7 + k*m.x8 + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
y.x6[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*x6.range[i] + h*m.x7 + j*x6.range[i]*m.x7 + k*m.x8 + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
y.x7[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*x7.range[i] + j*m.x6*x7.range[i] + k*m.x8 + g1 * exp(-1 * (pow((x7.range[i] - g2),2) / (2 * pow(g3, 2))))
y.x7.am[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*0 + h*x7.range[i] + j*0*x7.range[i] + k*m.x8 + g1 * exp(-1 * (pow((x7.range[i] - g2),2) / (2 * pow(g3, 2))))
y.x7.em[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*1 + h*x7.range[i] + j*1*x7.range[i] + k*m.x8 + g1 * exp(-1 * (pow((x7.range[i] - g2),2) / (2 * pow(g3, 2))))
y.x8[i] <- a + b*m.x1 + c*m.x2 + d*m.x3 + e*m.x4 + f*m.x5 + g*m.x6 + h*m.x7 + j*m.x6*m.x7 + k*x8.range[i] + g1 * exp(-1 * (pow((m.x7 - g2),2) / (2 * pow(g3, 2))))
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
g1 ~ dnorm(0, .0001) I(0, )
g2 ~ dnorm(0, .0001) I(0, )
g3 ~ dnorm(0, .0001) I(0, )
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)
}
"

data = list( y = y,    
             x1=x1, x2=x2,x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, x8=x8,
             x1.range=x1.range, x2.range=x2.range, x3.range=x3.range, x4.range=x4.range, x5.range=x5.range, x6.range=x6.range, x7.range=x7.range, x8.range=x8.range,
             m.x1=m.x1, m.x2=m.x2, m.x3=m.x3, m.x4=m.x4, m.x5=m.x5, m.x6=m.x6, m.x7=m.x7, m.x8=m.x8,
             N = length(y), N.pred = length(x1.range)
) 

#choose starting values
inits <-list()
inits[[1]] <- list(a = -5.8, b = 0, c = 0, d = 0, e = 0, f = 0, g=0, h=0, j=0, k=0,
                   g1=7, g2=16.0, g3=5.6)
inits[[2]] <- lapply(inits[[1]],"*",1.01)
inits[[3]] <- lapply(inits[[1]],"*",0.99)

#######################################################
##########          Run JAGS model      ###############
#######################################################
tic()
jags.out <- run.jags(jags.model,
                     data=data,
                     burnin = 200000,
                     sample = 10000,
                     inits = inits,
                     n.chains=12,
                     monitor=c('a','b','c','d','e','f','g','h','j','k','g1','g2','g3','tau'),
                     method='rjparallel'
)
toc()

#save if converged!
saveRDS(jags.out, model.path)

#now get predictions.
tic()
prediction.out <- extend.jags(jags.out, add.monitor = c('y.x1','y.x2','y.x3','y.x4','y.x5','y.x6','y.x7','y.x7.am','y.x7.em','y.x8'))
toc()

#save predictions!
saveRDS(prediction.out, prediction.path)
#save summary of predictions! (much smaller file to work with for plotting)
#inverse log transform so everything is in the right units on the other end. 
tic()
out <- summary(prediction.out)
out[,1:5] <- exp(out[,1:5])
toc()
saveRDS(out, summary.path)
