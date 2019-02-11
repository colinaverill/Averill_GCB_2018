#12-panel MAT, MAP, C:N and pH effects on growth, mortality and relative abundance.
rm(list=ls())
source('required_products_utilities/master_path_list.r')
#load output data for growth and mortality. 
#all outputs have been back transformed as required. 
g.out    <- readRDS(local_growth.gaus_summary_path)
m.out    <- readRDS(local_mortality_summary_path)
r.out.am <- readRDS(local_AM_recruitment_summary_path)
r.out.em <- readRDS(local_EM_recruitment_summary_path)
b.out    <- readRDS(local_beta_summary_path)


####GET DATA AND PREDICTOR RANGES####
#inverse logit transform m.out to get annual mortality probability.
require(boot)
m.out <- m.out[,1:5]
m.out <- (inv.logit(m.out))

#subset output by predictor for MAT, MAP, C:N and pH
g.mat     <- g.out[grep('y.x2', rownames(g.out)),]
g.map     <- g.out[grep('y.x3', rownames(g.out)),]
g.cn      <- g.out[grep('y.x4', rownames(g.out)),]
g.pH      <- g.out[grep('y.x5', rownames(g.out)),]
g.clay    <- g.out[grep('y.x8', rownames(g.out)),]
m.mat     <- m.out[grep('z.x2', rownames(m.out)),]
m.map     <- m.out[grep('z.x3', rownames(m.out)),]
m.cn      <- m.out[grep('z.x5', rownames(m.out)),]
m.pH      <- m.out[grep('z.x4', rownames(m.out)),]
m.clay    <- m.out[grep('z.x8', rownames(m.out)),]
r.mat.em  <- r.out.em[grep('mat.pred' , rownames(r.out.em)),]
r.map.em  <- r.out.em[grep('map.pred' , rownames(r.out.em)),]
r.ph.em   <- r.out.em[grep('ph.pred'  , rownames(r.out.em)),]
r.cn.em   <- r.out.em[grep('cn.pred'  , rownames(r.out.em)),]
r.clay.em <- r.out.em[grep('clay.pred', rownames(r.out.em)),]
r.mat.am  <- r.out.am[grep('mat.pred' , rownames(r.out.am)),]
r.map.am  <- r.out.am[grep('map.pred' , rownames(r.out.am)),]
r.ph.am   <- r.out.am[grep('ph.pred'  , rownames(r.out.am)),]
r.cn.am   <- r.out.am[grep('cn.pred'  , rownames(r.out.am)),]
r.clay.am <- r.out.am[grep('clay.pred', rownames(r.out.am)),]
b.mat     <- b.out[grep('z.x1',rownames(b.out)),]
b.map     <- b.out[grep('z.x2',rownames(b.out)),]
b.cn      <- b.out[grep('z.x3',rownames(b.out)),]
b.ph      <- b.out[grep('z.x4',rownames(b.out)),]
b.clay    <- b.out[grep('z.x6',rownames(b.out)),]

#we need the original ranges of each predictor.
range <- readRDS(local_growth.gaus_range_path)
g.mat.r  <- range$x2.range
g.map.r  <- range$x3.range
g.cn.r  <- range$x4.range
g.pH.r  <- range$x5.range
g.clay.r <- range$x8.range
range <- readRDS(local_mortality_range_path)
m.mat.r  <- range$range.x2
m.map.r  <- range$range.x3
m.cn.r   <- range$range.x5
m.pH.r   <- range$range.x4
m.clay.r <- range$range.x8
range <- readRDS(local_AM_recruitment_range_path)
r.mat.r <- range$x3.r
r.map.r <- range$x4.r
r.ph.r  <- range$x1.r
r.cn.r  <- range$x2.r
r.clay.r<- range$x7.r
range <- readRDS(local_beta_range_path)
b.mat.r <- range$mat.range
b.map.r <- range$map.range
b.cn.r  <- range$cn.range
b.ph.r  <- range$ph.range
b.clay.r<- range$clay.range

####START SETTING UP PLOTS####
#save dimensions, destination
png(filename=Supp_Figure6_all_path,width=8.1,height=7.2,units='in',res=300)

#begin the plots!
require(wesanderson) #wes anderson of course. 
cols <- wes_palette("Moonrise3", 5)
trans <- 0.3
par(mfrow=c(5,5), oma=c(5,12,2,1), mar=c(0,0,0,0),las=1)

#inner/outer label cex (size)
cex.ilab <- 1.0
cex.olab <- 1

####RELATIVE ABUNDANCE####
#y limits for relative abundance
limy=c(0,1)

#mat effect
plot(inv.logit(b.mat[,4]) ~ b.mat.r, lwd=0, ylim= limy, cex.axis = cex.olab, xaxt='n')
lines(smooth.spline(inv.logit(b.mat[,4]) ~ b.mat.r) , lwd=2, col=cols[1])
polygon(c(b.mat.r, rev(b.mat.r)),c(inv.logit(b.mat[,3]), rev(inv.logit(b.mat[,1]))), col=adjustcolor(cols[1], trans), lty=0)
#mtext('MAT', line = 0, cex=cex.ilab)
#mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('a.', adj=0.05,line = -1.7, side=1)
mtext(side=2,'Relative \n Abundance EM',cex.lab = cex.olab, las=0, line = 6)
mtext('MAT', line = 0, cex=cex.ilab)

#map effect
plot(inv.logit(b.map[,4]) ~ b.map.r, lwd=0, ylim= limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(inv.logit(b.map[,4]) ~ b.map.r) , lwd=2, col=cols[3])
polygon(c(b.map.r, rev(b.map.r)),c(inv.logit(b.map[,3]), rev(inv.logit(b.map[,1]))), col=adjustcolor(cols[3], trans), lty=0)
#mtext('MAP', line = 0, cex=cex.ilab)
#mtext(expression(paste('mm yr'^'-1')), side=1, line = 3)
mtext('MAP', line = 0, cex=cex.ilab)
mtext('b.', adj=0.05,line = -1.7, side=1)

#cn effect
plot(inv.logit(b.cn[,4]) ~ b.cn.r, lwd=0, ylim= limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(inv.logit(b.cn[,4]) ~ b.cn.r) , lwd=2, col=cols[2])
polygon(c(b.cn.r, rev(b.cn.r)),c(inv.logit(b.cn[,3]), rev(inv.logit(b.cn[,1]))), col=adjustcolor(cols[2], trans), lty=0)
#mtext(expression(paste('g C g'^'-1',' N')), side=1, line = 3)
mtext('C:N', line = 0, cex=cex.ilab)
mtext('c.', adj=0.05,line = -1.7, side=1)

#pH effect
plot(inv.logit(b.ph[,4]) ~ b.ph.r, lwd=0, ylim= limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(inv.logit(b.ph[,4]) ~ b.ph.r) , lwd=2, col=cols[4])
polygon(c(b.ph.r, rev(b.ph.r)),c(inv.logit(b.ph[,3]), rev(inv.logit(b.ph[,1]))), col=adjustcolor(cols[4], trans), lty=0)
mtext('pH', line = 0, cex=cex.ilab)
#mtext(expression(paste('-log'[10],' (mols H'^'+',')')), side=1, line = 3)
mtext('d.', adj=0.05,line = -1.7, side=1)

#clay effect - not significant!
plot(inv.logit(b.clay[,4]) ~ b.clay.r, lwd=0, ylim= limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(inv.logit(b.clay[,4]) ~ b.clay.r) , lwd=2, col=cols[5], lty=2)
polygon(c(b.clay.r, rev(b.clay.r)),c(inv.logit(b.clay[,3]), rev(inv.logit(b.clay[,1]))), col=adjustcolor(cols[5], trans), lty=0)
mtext('clay', line = 0, cex=cex.ilab)
#mtext(expression(paste('-log'[10],' (mols H'^'+',')')), side=1, line = 3)
mtext('e.', adj=0.05,line = -1.7, side=1)

####GROWTH####
#ylim values for growth plots
limy <- c(0,0.2)
#limy<- c(1,1.25)

#mat effect
plot(g.mat[,4] ~ g.mat.r, lwd=0, ylim= limy, cex.axis = cex.olab, xaxt='n')
lines(smooth.spline(g.mat[,4] ~ g.mat.r) , lwd=2, lty=1, col=cols[1])
polygon(c(g.mat.r, rev(g.mat.r)),c(g.mat[,3], rev(g.mat[,1])), col=adjustcolor(cols[1], trans), lty=0)
#mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('f.', adj=0.05,line = -1.7, side=1)
mtext(side=2,'Growth of \n Surviving Trees',cex.lab = cex.olab, las=0, line = 6)
mtext(side=2,expression(paste('cm'^'2',' m'^'-2',' yr'^'-1')), las=0, line = 4)

#map effect
plot(g.map[,4] ~ g.map.r, lwd=0, ylim= limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(g.map[,4] ~ g.map.r) , lwd=2, col=cols[3])
polygon(c(g.map.r, rev(g.map.r)),c(g.map[,3], rev(g.map[,1])), col=adjustcolor(cols[3], trans), lty=0)
#mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('g.', adj=0.05,line = -1.7, side=1)

#cn effect- not significant!
plot(g.cn[,4] ~ g.cn.r, lwd=0, ylim= limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(g.cn[,4] ~ g.cn.r) , lwd=2, lty=2, col=cols[2])
polygon(c(g.cn.r, rev(g.cn.r)),c(g.cn[,3], rev(g.cn[,1])), col=adjustcolor(cols[2], trans), lty=0)
#mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('h.', adj=0.05,line = -1.7, side=1)


#pH effect
plot(g.pH[,4] ~ g.pH.r, lwd=0, ylim= limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(g.pH[,4] ~ g.pH.r) , lwd=2, col=cols[4])
polygon(c(g.pH.r, rev(g.pH.r)),c(g.pH[,3], rev(g.pH[,1])), col=adjustcolor(cols[4], trans), lty=0)
#mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('i.', adj=0.05,line = -1.7, side=1)

#clay effect
plot(g.clay[,4] ~ g.clay.r, lwd=0, ylim= limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(g.clay[,4] ~ g.clay.r) , lwd=2, col=cols[4])
polygon(c(g.clay.r, rev(g.clay.r)),c(g.clay[,3], rev(g.clay[,1])), col=adjustcolor(cols[5], trans), lty=0)
#mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('j.', adj=0.05,line = -1.7, side=1)

####EM RECRUIT####
#limits for EM recruitment
limy=c(0,0.88)

#mat effect + +/- 1 sd here. 
plot(r.mat.em[,4] ~ r.mat.r, lwd=0, ylim=limy, cex.axis = cex.olab, xaxt='n')
lines(smooth.spline(r.mat.em[,4] ~ r.mat.r), lwd=2, col=cols[1])
polygon(c(r.mat.r, rev(r.mat.r)),c(r.mat.em[,3], rev(r.mat.em[,1])), col=adjustcolor(cols[1], trans), lty=0)
mtext('k.', adj = 0.05, line=-1.7, side=1)
mtext(side=2,'Recruitment \n EM Trees',cex.lab = cex.olab, las=0, line = 6)
mtext(side=2,expression(paste('stems plot'^'-1','yr'^'-1')), las=0, line = 4)

#map effect 
plot(r.map.em[,4] ~ r.map.r, lwd=0, ylim=limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(r.map.em[,4] ~ r.map.r), lwd=2, col=cols[3])
polygon(c(r.map.r, rev(r.map.r)),c(r.map.em[,3], rev(r.map.em[,1])), col=adjustcolor(cols[3], trans), lty=0)
mtext('l.', adj = 0.05, line=-1.7, side=1)

#cn effect
plot(r.cn.em[,4] ~ r.cn.r, lwd=0, ylim=limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(r.cn.em[,4] ~ r.cn.r), lwd=2, col=cols[2])
polygon(c(r.cn.r, rev(r.cn.r)),c(r.cn.em[,3], rev(r.cn.em[,1])), col=adjustcolor(cols[2], trans), lty=0)
mtext('m.', adj = 0.05, line=-1.7, side=1)

#pH effect
plot(r.ph.em[,4] ~ r.ph.r, lwd=0, ylim=limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(r.ph.em[,4] ~ r.ph.r), lwd=2, col=cols[4])
polygon(c(r.ph.r, rev(r.ph.r)),c(r.ph.em[,3], rev(r.ph.em[,1])), col=adjustcolor(cols[4], trans), lty=0)
mtext('n.', adj = 0.05, line=-1.7, side=1)

#clay effect
plot(r.clay.em[,4] ~ r.clay.r, lwd=0, ylim=limy, cex.axis = 1.5, xaxt='n', yaxt='n')
lines(smooth.spline(r.clay.em[,4] ~ r.clay.r), lwd=2, lty=2, col=cols[5])
polygon(c(r.clay.r, rev(r.clay.r)),c(r.clay.em[,3], rev(r.clay.em[,1])), col=adjustcolor(cols[5], trans), lty=0)
mtext('o.', adj = 0.05, line=-1.7, side=1)

####AM RECRUIT####
#limits for AM recruitment
limy=c(0,0.88)

#mat effect
plot(r.mat.am[,4] ~ r.mat.r, lwd=0, ylim=limy, cex.axis = cex.olab, xaxt='n')
lines(smooth.spline(r.mat.am[,4] ~ r.mat.r), lwd=2, lty=1, col=cols[1])
polygon(c(r.mat.r, rev(r.mat.r)),c(r.mat.am[,3], rev(r.mat.am[,1])), col=adjustcolor(cols[1], trans), lty=0)
mtext('p.', adj = 0.05, line=-1.7, side=1)
mtext(side=2,'Recruitment \n AM Trees',cex.lab = cex.olab, las=0, line = 6)
mtext(side=2,expression(paste('stems plot'^'-1','yr'^'-1')), las=0, line = 4)

#map effect 
plot(r.map.am[,4] ~ r.map.r, lwd=0, ylim=limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(r.map.am[,4] ~ r.map.r), lwd=2, col=cols[3])
polygon(c(r.map.r, rev(r.map.r)),c(r.map.am[,3], rev(r.map.am[,1])), col=adjustcolor(cols[3], trans), lty=0)
mtext('q.', adj = 0.05, line=-1.7, side=1)

#cn effect 
plot(r.cn.am[,4] ~ r.cn.r, lwd=0, ylim=limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(r.cn.am[,4] ~ r.cn.r), lwd=2, lty=2, col=cols[2])
polygon(c(r.cn.r, rev(r.cn.r)),c(r.cn.am[,3], rev(r.cn.am[,1])), col=adjustcolor(cols[2], trans), lty=0)
mtext('r.', adj = 0.05, line=-1.7, side=1)

#pH effect
plot(r.ph.am[,4] ~ r.ph.r, lwd=0, ylim=limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(r.ph.am[,4] ~ r.ph.r), lwd=2, col=cols[4])
polygon(c(r.ph.r, rev(r.ph.r)),c(r.ph.am[,3], rev(r.ph.am[,1])), col=adjustcolor(cols[4], trans), lty=0)
mtext('s.', adj = 0.05, line=-1.7, side=1)

#clay effect - not significant!
plot(r.clay.am[,4] ~ r.clay.r, lwd=0, ylim=limy, cex.axis = cex.olab, xaxt='n', yaxt='n')
lines(smooth.spline(r.clay.am[,4] ~ r.clay.r), lwd=2, col=cols[5], lty = 2)
polygon(c(r.clay.r, rev(r.clay.r)),c(r.clay.am[,3], rev(r.clay.am[,1])), col=adjustcolor(cols[5], trans), lty=0)
mtext('t.', adj = 0.05, line=-1.7, side=1)

####MORTALITY PLOT####
#set mortality y limits
limy = c(0,0.030)

#mat effect - not significant!
plot(m.mat[,4] ~ m.mat.r, lwd=0, ylim= limy, cex.axis = cex.olab)
lines(smooth.spline(m.mat[,4] ~ m.mat.r) , lwd=2, lty=2, col=cols[1])
polygon(c(m.mat.r, rev(m.mat.r)),c(m.mat[,3], rev(m.mat[,1])), col=adjustcolor(cols[1], trans), lty=0)
#mtext('MAT', line = 0, cex=cex.ilab)
mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('u.', adj=0.05,line = -1.7, side=1)
mtext(side=2,'Annual Mortality \n Probability',cex.lab = cex.olab, las=0, line = 6)


#map effect
plot(m.map[,4] ~ m.map.r, lwd=0, ylim= limy, cex.axis = cex.olab,yaxt='n')
lines(smooth.spline(m.map[,4] ~ m.map.r) , lwd=2, col=cols[3])
polygon(c(m.map.r, rev(m.map.r)),c(m.map[,3], rev(m.map[,1])), col=adjustcolor(cols[3], trans), lty=0)
#mtext('MAP', line = 0, cex=cex.ilab)
mtext(expression(paste('mm yr'^'-1')), side=1, line = 3)
mtext('v.', adj=0.05,line = -1.7, side=1)

#cn effect
plot(m.cn[,4] ~ m.cn.r, lwd=0, ylim= limy, cex.axis = cex.olab, yaxt='n')
lines(smooth.spline(m.cn[,4] ~ m.cn.r) , lwd=2, col=cols[2])
polygon(c(m.cn.r, rev(m.cn.r)),c(m.cn[,3], rev(m.cn[,1])), col=adjustcolor(cols[2], trans), lty=0)
#mtext('C:N', line = 0, cex=cex.ilab)
mtext(expression(paste('g C g'^'-1',' N')), side=1, line = 3)
mtext('w.', adj=0.05,line = -1.7, side=1)


#pH effect
plot(m.pH[,4] ~ m.pH.r, lwd=0, ylim= limy, cex.axis = cex.olab,yaxt='n')
lines(smooth.spline(m.pH[,4] ~ m.pH.r) , lwd=2, col=cols[4])
polygon(c(m.pH.r, rev(m.pH.r)),c(m.pH[,3], rev(m.pH[,1])), col=adjustcolor(cols[4], trans), lty=0)
#mtext('pH', line = 0, cex=cex.ilab)
mtext(expression(paste('-log'[10],' (mols H'^'+',')')), side=1, line = 3)
mtext('x.', adj=0.05,line = -1.7, side=1)

#clay effect - not significant!
plot(m.clay[,4] ~ m.clay.r, lwd=0, ylim= limy, cex.axis = cex.olab,yaxt='n')
lines(smooth.spline(m.clay[,4] ~ m.clay.r) , lwd=2, col=cols[5], lty=2)
polygon(c(m.clay.r, rev(m.clay.r)),c(m.clay[,3], rev(m.clay[,1])), col=adjustcolor(cols[5], trans), lty=0)
#mtext('pH', line = 0, cex=cex.ilab)
mtext('percent', side=1, line = 3)
mtext('y.', adj=0.05,line = -1.7, side=1)

####end plot####
dev.off()
