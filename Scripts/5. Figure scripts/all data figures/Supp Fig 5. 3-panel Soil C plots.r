#plotting soil effects
rm(list=ls())
source('required_products_utilities/master_path_list.r')
d     <- readRDS(local_soils_summary_path)
range <- readRDS(local_soils_range_path)


d <- d[,1:5]

#subset output by predictor for MAT, MAP, C:N and pH
mat     <- d[grep('y.x1', rownames(d)),];  mat <- exp(mat )
map     <- d[grep('y.x2', rownames(d)),];  map <- exp(map )
pH      <- d[grep('y.x3', rownames(d)),];   pH <- exp(pH  )
ndep    <- d[grep('y.x4', rownames(d)),]; ndep <- exp(ndep)
em      <- d[grep('y.x5', rownames(d)),];   em <- exp(em  )
clay    <- d[grep('y.x7', rownames(d)),]; clay <- exp(clay)

#grab ranges
 mat.r  <- range$x1.range
 map.r  <- range$x2.range
  pH.r  <- range$x3.range
ndep.r  <- range$x4.range
  em.r  <- range$x5.range * 100
clay.r  <- range$x7.range

#save dimensions, destination
png(filename=Supp_Figure5_all_path,width=7,height=7,units='in',res=300)

#begin the plots!
require(wesanderson) #wes anderson of course. 
cols <- wes_palette("Moonrise3", 5)
trans <- 0.3
par(mfrow=c(2,2), oma=c(1,8,0,1), mar=c(3,0,3,0),las=1) #number of panels and set margins

#inner/outer label cex (size)
cex.ilab <- 1.0
cex.olab <- 1.0

#ylim values for growth plots
limy <- c(0,10000)

#mat effect
plot(mat[,4] ~ mat.r, lwd=0, ylim= limy, cex.axis = 1.2)
lines(smooth.spline(mat[,4] ~ mat.r) , lwd=2, col=cols[4])
polygon(c(mat.r, rev(mat.r)),c(mat[,3], rev(mat[,1])), col=adjustcolor(cols[4], trans), lty=0)
mtext('MAT', line = -2, cex=cex.ilab)
mtext(expression(paste('C'^'o')), side=1, line = 3)
mtext('a.', adj=0.05,line = -1.7, side=1)

#map effect
plot(map[,4] ~ map.r, lwd=0, ylim= limy, cex.axis = 1.2, yaxt='n')
lines(smooth.spline(map[,4] ~ map.r) , lwd=2, col=cols[5])
polygon(c(map.r, rev(map.r)),c(map[,3], rev(map[,1])), col=adjustcolor(cols[5], trans), lty=0)
mtext('MAP', line = -2, cex=cex.ilab)
mtext(expression(paste('mm yr'^'-1')), side=1, line = 3)
mtext('b.', adj=0.05,line = -1.7, side=1)

#pH effect
plot(pH[,4] ~ pH.r, lwd=0, ylim= limy, cex.axis = 1.2)
lines(smooth.spline(pH[,4] ~ pH.r) , lwd=2, col=cols[1])
polygon(c(pH.r, rev(pH.r)),c(pH[,3], rev(pH[,1])), col=adjustcolor(cols[1], trans), lty=0)
mtext('pH', line = -2, cex=cex.ilab)
mtext(expression(paste('-log'[10],' (mols H'^'+',')')), side=1, line = 3)
mtext('c.', adj=0.05,line = -1.7, side=1)

#clay effect
plot(clay[,4] ~ clay.r, lwd=0, ylim= limy, cex.axis = 1.2, yaxt='n')
lines(smooth.spline(clay[,4] ~ clay.r) , lwd=2, col=cols[2])
polygon(c(clay.r, rev(clay.r)),c(clay[,3], rev(clay[,1])), col=adjustcolor(cols[2], trans), lty=0)
mtext('clay', line = -2, cex=cex.ilab)
mtext('percent', side=1, line = 3)
mtext('d.', adj=0.05,line = -1.7, side=1)

mtext(side=2,'Soil Carbon',cex = 1.5, las=0, line = 6, out = T)
mtext(side=2,expression(paste('g C m'^'-2')), las=0, line = 4, out = T)

#end figure
dev.off()