#Plotting N-dep fits to detrended, binned means.
rm(list=ls())
library(betareg)
library(data.table)
library(runjags)
library(wesanderson) #wes anderson of course. 
library(boot) #for inv.logit command

#output file path
file.save.path <- 'Figures/east_data_figures/Supp Fig X. binned mean fits.png'

#load response curves to N-deposition with 95% credible intervals of the mean.
g.out    <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_051517_growth.summary.Rdata')
m.out    <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_runjags_EAST_mortality.summary.Rdata')
r.out.am <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_AM.adult_EAST_linear_recruit.summary.Rdata')
r.out.em <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_EM.adult_EAST_linear_recruit.summary.Rdata')
b.out    <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_east_053017_rel.abundance.summary.rds')


#A few quick transformations.
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
range <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_051517_growth.predictor_ranges.rds')
g.range <- range$x7.range
range <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_runjags_EAST_mortality.predictor_ranges.rds')
m.range <- range$range.x6
r.range <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_EAST_052617_AM.adult_EAST_linear_recruit.predictor_ranges.rds')
r.range <- r.range$x5.r
b.range <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/FIA7_east_053017_rel.abundance.ranges.rds')
b.range <- b.range$ndep.range

#load up the binned mean response data, after it has been detrended for other predictor effects.
b.bins    <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_beta_bins.rds')
g.bins.am <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_growth.am_bins.rds')
g.bins.em <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_growth.em_bins.rds')
r.bins    <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_recruit_bins.rds')
m.bins    <- readRDS('/fs/data3/caverill/FIA7_Ndep_data.processed/east.data_analysis_output/east_mort_bins.rds')

#start plots.
#save dimensions, destination
png(filename=file.save.path,width=8,height=8,units='in',res=300)

cols <- wes_palette("Zissou", 5)
trans <- 0.3
par(mfrow=c(3,3),
    mar=c(2,2,2,2),
    oma=c(5,3,1,0.5))
limx<- c(2.8,20)
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
mtext(side=3,expression(paste('R'^'2','=0.86')), adj = 0.05, line = -1.75,cex = cex.ilab)
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
#labels
mtext(side=2,expression(paste('unitless')), las=0, line = 2, cex=0.75)
mtext('Mortality', line = .5, cex=cex.ilab)
mtext(side=3,expression(paste('R'^'2','= NA')), adj = 0.05, line = -1.75,cex = cex.ilab)
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
#labels
mtext(side=2,expression(paste('cm'^'2','m'^'-2',' yr'^'-1')), las=0, line = 2, cex = 0.75)
mtext('Growth AM', line = .5, cex=cex.ilab)
mtext(side=3,expression(paste('R'^'2','=0.54')), adj = 0.05, line = -1.75,cex = cex.ilab)
mtext('c.', side = 1, line = -1.25, adj = 0.975, cex = cex.ilab)

#EM
plot(g.bins.em$detrend ~ g.bins.em$tot.15, ylim=limy, xlim=limx, cex.axis = cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(g.bins.em$tot.15, g.bins.em$detrend.lwr, g.bins.em$tot.15, g.bins.em$detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(g.ndep.em[,4] ~ g.range) , lwd=2, col=cols[5])
polygon(c(g.range, rev(g.range)),c(g.ndep.em[,3], rev(g.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(g.ndep.em[,1] ~ g.range) , lty=2, lwd=0.5)
lines(smooth.spline(g.ndep.em[,3] ~ g.range) , lty=2, lwd=0.5)
#labels
mtext('Growth EM', line = .5, cex=cex.ilab)
mtext(side=3,expression(paste('R'^'2','=0.77')), adj = 0.05, line = -1.75,cex = cex.ilab)
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
legend(0, 0.26,c('AM','EM'), lwd=2, col=c(cols[1],cols[5]), bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.5, seg.len=1.3)


#recruitment plots
limy <- c(0,1.1)
#AM
plot(r.bins$am.detrend ~ r.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(r.bins$tot.15, r.bins$am.detrend.lwr, r.bins$tot.15, r.bins$am.detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(r.ndep.am[,4] ~ r.range) , lwd=2, col=cols[1])
polygon(c(r.range, rev(r.range)),c(r.ndep.am[,3] , rev(r.ndep.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
lines(smooth.spline(r.ndep.am[,3]  ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.am[,1]  ~ r.range), lty=2, lwd=0.5)
#labels
mtext('Recruitment AM', line = .5, cex=cex.ilab)
mtext(side=2,expression(paste('stems yr'^'-1')), las=0, line = 2, cex=0.75)
mtext(side=3,expression(paste('R'^'2','= NA')), adj = 0.05, line = -1.75,cex = cex.ilab)
mtext('f.', side = 1, line = -1, adj = 0.95, cex = cex.ilab)

#EM
plot(r.bins$em.detrend ~ r.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
arrows(r.bins$tot.15, r.bins$em.detrend.lwr, r.bins$tot.15, r.bins$em.detrend.upr, length=0.05, angle=90, code=3)
lines(smooth.spline(r.ndep.em[,4] ~ r.range) , lwd=2, col=cols[5])
polygon(c(r.range, rev(r.range)),c(r.ndep.em[,3], rev(r.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
lines(smooth.spline(r.ndep.em[,1] ~ r.range), lty=2, lwd=0.5)
lines(smooth.spline(r.ndep.em[,3] ~ r.range), lty=2, lwd=0.5)
#labels
mtext('Recruitment EM', line = .5, cex=cex.ilab)
mtext(side=3,expression(paste('R'^'2','=0.78')), adj = 0.05, line = -1.75,cex = cex.ilab)
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
legend(0, 0.625,c('AM','EM'), lwd=2, col=c(cols[1],cols[5]), bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.5, seg.len=1.3)


#mortality plots - when I plotted them separately by AM vs. EM
#limy <- c(0,0.035)
#AM
#plot(m.bins$detrend.am~m.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex =1.5, xlab=NA, ylab=NA)
#lines(smooth.spline(m.ndep.am[,4] ~ m.range) , lwd=2, col=cols[1], lty=1)
#polygon(c(m.range, rev(m.range)),c(m.ndep.am[,3], rev(m.ndep.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
#lines(smooth.spline(m.ndep.am[,1] ~ m.range), lty=2, lwd=0.5)
#lines(smooth.spline(m.ndep.am[,3] ~ m.range), lty=2, lwd=0.5)
#labels
#mtext(side=2,expression(paste('unitless')), las=0, line = 2, cex=0.75)
#mtext('Mortality AM', line = .5, cex=cex.ilab)
#EM
#plot(m.bins$detrend.em~m.bins$tot.15, ylim=limy, xlim=limx, cex.axis=cex.olab, pch=16, cex=1.5, xlab=NA, ylab=NA)
#lines(smooth.spline(m.ndep.em[,4] ~ m.range) , lwd=2, col=cols[5], lty=1)
#polygon(c(m.range, rev(m.range)),c(m.ndep.em[,3], rev(m.ndep.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
#lines(smooth.spline(m.ndep.em[,1] ~ m.range), lty=2, lwd=0.5)
#lines(smooth.spline(m.ndep.em[,3] ~ m.range), lty=2, lwd=0.5)
#labels
#mtext('Mortality EM', line = .5, cex=cex.ilab)

#outer labels
mtext('Nitrogen Deposition', side=1, out=T, cex=1.5, line = 1.2)
mtext(expression(paste('kg N ha'^'-2',' yr'^'-1')), side=1, out=T, line = 3.5, cex = 1)
dev.off()
