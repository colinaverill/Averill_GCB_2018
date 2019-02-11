#Plotting N-dep fits to detrended, binned means.
rm(list=ls())
library(betareg)
library(data.table)
library(runjags)
library(wesanderson) #wes anderson of course. 
library(boot) #for inv.logit command
source('required_products_utilities/master_path_list.r')

#load response curves to N-deposition with 95% credible intervals of the mean.
m.out    <- readRDS(local_mortality_summary_path)
r.out.am <- readRDS(local_AM_recruitment_summary_path)
r.out.em <- readRDS(local_EM_recruitment_summary_path)
b.out    <- readRDS(local_beta_summary_path)

#load rsq table
rsq.table <- readRDS(all_binned_rsq_summary_path)

#g.out is special due to unmixing chain.
g.out <- readRDS(local_growth.gaus_prediction_path)
pred.trim <- coda::as.mcmc.list(g.out)
pred.trim <- mcmc.list(pred.trim[[1]],pred.trim[[3]])
out <- summary(pred.trim)
mean.out <- out$statistics
cred.out <- out$quantiles
pred.out <- data.frame(cred.out[,1],cred.out[,3],cred.out[,5],mean.out[,1])
colnames(pred.out) <- c('2.5%','50%','97.5%','mean')
#grab the mean responses to N in AM and EM conditions
ndep <- pred.out[grep('y.x7',rownames(pred.out)),]
g.ndep.am <- exp(ndep[102:202,])
g.ndep.em <- exp(ndep[203:303,])
g.ndep    <- exp(ndep[1:101,])

#A few quick transformations.
#g.out <- exp(g.out[,1:5]) #no longer needed with the two unmixed chains. 
m.out <- (inv.logit(m.out[,1:5]))

#subsetresponses to N. 
g.ndep.am <-    g.out[grep('y.x7.am', rownames(g.out)),]
g.ndep.em <-    g.out[grep('y.x7.em', rownames(g.out)),]
m.ndep.am <-    m.out[grep('z.x6.am', rownames(m.out)),]
m.ndep.em <-    m.out[grep('z.x6.em', rownames(m.out)),]
m.ndep    <- m.out[grep('x6',    rownames(m.out)),]
m.ndep    <- m.ndep[1:101,]
r.ndep.am <- r.out.am[grep('ndep.pred', rownames(r.out.am)),]
r.ndep.em <- r.out.em[grep('ndep.pred', rownames(r.out.em)),]
r.ndep.am.pois <- r.out.am[grep('ndep.mu.pred', rownames(r.out.am)),]
r.ndep.em.pois <- r.out.em[grep('ndep.mu.pred', rownames(r.out.em)),]
r.ndep.am.bin  <- r.out.am[grep('ndep.pro.pred', rownames(r.out.am)),]
r.ndep.em.bin  <- r.out.em[grep('ndep.pro.pred', rownames(r.out.em)),]
b.ndep    <-    b.out[grep('x5'  , rownames(b.out)),]
b.mean <- inv.logit(b.ndep[,4])
b.upr  <- inv.logit(b.ndep[,3])
b.lwr  <- inv.logit(b.ndep[,1])

#we need the Ndep ranges of the original JAGS predictors predictor_ranges
range <- readRDS(local_growth.gaus_range_path)
g.range <- range$x7.range
range <- readRDS(local_mortality_range_path)
m.range <- range$range.x6
r.range <- readRDS(local_EM_recruitment_range_path)
r.range <- r.range$x5.r
b.range <- readRDS(local_beta_range_path)
b.range <- b.range$ndep.range

#load up the binned mean response data, after it has been detrended for other predictor effects.
b.bins <- readRDS(bins_beta_all_path)
g.bins.am <- readRDS(bins_growth.gaus_all_AM_path)
g.bins.em <- readRDS(bins_growth.gaus_all_EM_path)
r.bins <- readRDS(bins_recruitment_all_path)
m.bins <- readRDS(bins_mortality_all_path)
#g.bins$detrend.am <- exp(g.bins$detrend.am)
#g.bins$detrend.em <- exp(g.bins$detrend.em)
#g.bins[,11:14] <- exp(g.bins[,11:14])

#start plots.
#save dimensions, destination
png(filename=Figure2_all_path,width=8,height=8,units='in',res=300)

cols <- wes_palette("Zissou", 5)
trans <- 0.3
par(mfrow=c(3,3),
    mar=c(2,2,2,2),
    oma=c(5,3,1,0.5))
limx<- c(0,20)
#inner/outer label cex (size)
cex.ilab <- 1.0
cex.olab <- 1.2

#relative abundance plot
limy<- c(0,1)
plot(b.bins$relEM.detrend ~ b.bins$tot.15, ylim= limy, xlim=limx, cex.axis = cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(b.bins$tot.15, b.bins$relEM.detrend.lwr, b.bins$tot.15, b.bins$relEM.detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(b.mean ~ b.range) , lwd=2, col=cols[3])
polygon(c(b.range, rev(b.range)),c(b.upr, rev(b.lwr)), col=adjustcolor(cols[3], trans), lty=0)
lines(smooth.spline(b.upr ~ b.range), lty=2, lwd=0.5)
lines(smooth.spline(b.lwr ~ b.range), lty=2, lwd=0.5)
#label
r.sq <- rsq.table[analyses == 'beta',]$r.sq
mylabel = bquote(italic(R)^2 == .(format(r.sq, digits = 2)))
mtext(side=3, mylabel, adj = 0.05, line = -1.75,cex = cex.ilab)
mtext(side=2,expression(paste('EM basal area *')), las=0, line = 3, cex = 0.75)
mtext(side=2,expression(paste('(total basal area)'^'-1')), las=0, line = 2, cex = 0.75)
mtext('Relative Abundance EM', line = .75, cex=cex.ilab)
mtext('a.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)

#drop a mortality plot
limy <- c(0,0.025)
plot(m.bins$detrend ~ m.bins$tot.15, ylim=limy, xlim=limx, cex.axis = cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(m.bins$tot.15, m.bins$detrend.lwr, m.bins$tot.15, m.bins$detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline((m.ndep[,4]) ~ m.range) , lwd=2, col=cols[3])
polygon(c(m.range, rev(m.range)),c((m.ndep[,3]), rev((m.ndep[,1]))), col=adjustcolor(cols[3], trans), lty=0)
lines(smooth.spline((m.ndep[,3]) ~ m.range), lty=2, lwd=0.5)
lines(smooth.spline((m.ndep[,1]) ~ m.range), lty=2, lwd=0.5)
#r2 label
r.sq <- rsq.table[analyses == 'mortality',]$r.sq
mylabel = bquote(italic(R)^2 == .(format(r.sq, digits = 2)))
mtext(side=3, mylabel, adj = 0.05, line = -1.75,cex = cex.ilab)
#labels
mtext(side=2,expression(paste('year'^'-1')), las=0, line = 2, cex=0.75)
mtext('Mortality', line = .5, cex=cex.ilab)
mtext('b.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)

#drop an empty plot
plot.new()

#growth plots
limy <- c(0,0.25)
#AM
plot(g.bins.am$detrend ~ g.bins.am$tot.15, ylim=limy, xlim=limx, cex.axis = cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(g.bins.am$tot.15, g.bins.am$detrend.lwr, g.bins.am$tot.15, g.bins.am$detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(g.ndep.am[,4] ~ g.range) , lwd=2, col=cols[2])
polygon(c(g.range, rev(g.range)),c(g.ndep.am[,3], rev(g.ndep.am[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(g.ndep.am[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.am[,3] ~ g.range) , lty=2, lwd=0.5)
#r2 label
r.sq <- rsq.table[analyses == 'growth.gaus.AM',]$r.sq
mylabel = bquote(italic(R)^2 == .(format(r.sq, digits = 2)))
mtext(side=3, mylabel, adj = 0.05, line = -1.75,cex = cex.ilab)
#labels
mtext(side=2,expression(paste('cm'^'2','m'^'-2',' yr'^'-1')), las=0, line = 2, cex = 0.75)
mtext('Growth AM', line = .5, cex=cex.ilab)
mtext('c.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)

#EM
plot(g.bins.em$detrend ~ g.bins.em$tot.15, ylim=limy, xlim=limx, cex.axis = cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(g.bins.em$tot.15, g.bins.em$detrend.lwr, g.bins.em$tot.15, g.bins.em$detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(g.ndep.em[,4] ~ g.range) , lwd=2, col=cols[5])
polygon(c(g.range, rev(g.range)),c(g.ndep.em[,3], rev(g.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(g.ndep.em[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,3] ~ g.range) , lty=2, lwd=0.5)
#r2 label
r.sq <- rsq.table[analyses == 'growth.gaus.EM',]$r.sq
mylabel = bquote(italic(R)^2 == .(format(r.sq, digits = 2)))
mtext(side=3, mylabel, adj = 0.05, line = -1.75,cex = cex.ilab)
#labels
mtext('Growth EM', line = .5, cex=cex.ilab)
mtext('d.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)

#drop the curves crossing.
plot(g.bins.am$detrend ~ g.bins.am$tot.15, ylim=limy, xlim=limx, cex.axis = cex.olab, pch=16, cex =0, xlab=NA, ylab=NA)
lines(smooth.spline(g.ndep.am[,4] ~ g.range) , lwd=2, col=cols[2])
polygon(c(g.range, rev(g.range)),c(g.ndep.am[,3], rev(g.ndep.am[,1])), col=adjustcolor(cols[2], trans), lty=0)
lines(smooth.spline(g.ndep.am[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.am[,3] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,4] ~ g.range) , lwd=2, col=cols[5])
polygon(c(g.range, rev(g.range)),c(g.ndep.em[,3], rev(g.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(g.ndep.em[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,3] ~ g.range) , lty=2, lwd=0.5)
#labels
mtext('Growth', line = .5, cex=cex.ilab)
mtext('e.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)
#legend
legend(0, limy[2]*1.05,c('AM','EM'), lwd=2, col=c(cols[1],cols[5]), bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.5, seg.len=1.3)


#recruitment plots
limy <- c(0,0.65)
#AM
plot(r.bins$am.detrend ~ r.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(r.bins$tot.15, r.bins$am.detrend.lwr, r.bins$tot.15, r.bins$am.detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(r.ndep.am[,4] ~ r.range) , lwd=2, col=cols[1])
polygon(c(r.range, rev(r.range)),c(r.ndep.am[,3] , rev(r.ndep.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
lines(smooth.spline(r.ndep.am[,3]  ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.am[,1]  ~ r.range), lty=2, lwd=0.5)
#r2 label
r.sq <- rsq.table[analyses == 'recruitment.AM',]$r.sq
mylabel = bquote(italic(R)^2 == .(format(r.sq, digits = 2)))
mtext(side=3, mylabel, adj = 0.05, line = -1.75,cex = cex.ilab)
#labels
mtext('Recruitment AM', line = .5, cex=cex.ilab)
mtext(side=2,expression(paste('stems yr'^'-1')), las=0, line = 2, cex=0.75)
mtext('f.', side = 1, line = -1, adj = 0.95, cex = cex.ilab)

#EM
plot(r.bins$em.detrend ~ r.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(r.bins$tot.15, r.bins$em.detrend.lwr, r.bins$tot.15, r.bins$em.detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(r.ndep.em[,4] ~ r.range) , lwd=2, col=cols[5])
polygon(c(r.range, rev(r.range)),c(r.ndep.em[,3], rev(r.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(r.ndep.em[,1] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.em[,3] ~ r.range), lty=2, lwd=0.5)
#r2 label
r.sq <- rsq.table[analyses == 'recruitment.EM',]$r.sq
mylabel = bquote(italic(R)^2 == .(format(r.sq, digits = 2)))
mtext(side=3, mylabel, adj = 0.05, line = -1.75,cex = cex.ilab)
#labels
mtext('Recruitment EM', line = .5, cex=cex.ilab)
mtext('g.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)


#drop curves
plot(r.bins$am.detrend ~ r.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex =0, xlab=NA, ylab=NA)
lines(smooth.spline(r.ndep.am[,4] ~ r.range) , lwd=2, col=cols[1])
polygon(c(r.range, rev(r.range)),c(r.ndep.am[,3] , rev(r.ndep.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
lines(smooth.spline(r.ndep.am[,3]  ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.am[,1]  ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.em[,4] ~ r.range) , lwd=2, col=cols[5])
polygon(c(r.range, rev(r.range)),c(r.ndep.em[,3], rev(r.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(r.ndep.em[,1] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.em[,3] ~ r.range), lty=2, lwd=0.5)
#labels
mtext('Recruitment', line = .5, cex=cex.ilab)
mtext('h.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)
#legend
legend(0, limy[2]*1.05,c('AM','EM'), lwd=2, col=c(cols[1],cols[5]), bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.5, seg.len=1.3)


#outer labels
mtext('Nitrogen Deposition', side=1, out=T, cex=1.5, line = 1.2)
mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, out=T, line = 3.5, cex = 1)
dev.off()
