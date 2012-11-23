#!/usr/bin/env Rscript --vanilla --slave


invisible(library( ggplot2 ))
options(warn=-1)###Suppress wraning messages (there are some in printing the graphs due to data points off the axes).

tmap <- read.table( "trim.gtf.tmap", header = T )
dists <- read.table( "distances.txt", header=T )

dist.cov <- merge( dists, tmap[,c(2,10)], by.x="Gene", by.y="ref_id" )

###Bin genes by expression (20 bins) and calculate MADs for bins. plot MADS for bins, add a line showing MAD for entire dataset (Or MAD for bin - MAD for dataset)
mad.all <- mad( dist.cov$Distance)
dist.cov <- dist.cov[ order( dist.cov$gene.cov ), ]

bin <- cut( c(1:length(dist.cov$gene.cov)), 20 )
bin <- as.numeric( bin )
dist.cov <- cbind( dist.cov, bin )

mads <- tapply( dist.cov$Distance, dist.cov$bin, mad )
expr <- tapply( dist.cov$gene.cov, dist.cov$bin, median )
mads <- as.data.frame(cbind( mads, expr ))

###Calculate expression level cutoff (this is the median expression of genes in the first bin below the population MAD).
exp.cutoff <- mads$expr[ mads$mads < mad.all ][1]


print( c("The predicted minimum gene expression (per-base coverage) for accurate predictions is:-", exp.cutoff,"This is an estimate, and may not be correct. Check ExprCutoffPlot.jpg to decide how believable this is.") )

###Produce MAD-expression plots
jpeg("MadExpressionPlot.jpg")
ggplot( mads, aes( log10(expr), mads)) + geom_point( ) + theme_bw() + geom_vline( aes( xintercept=1 ) ) + geom_hline( aes( yintercept=mad.all), linetype="dashed") + coord_trans( y="log10" )  + scale_y_continuous(  name = "Bin MAD of predicted TSS to nearest annotated TSS" ) + scale_x_continuous( name="Bin median per-base expression (log10)" )
invisible(dev.off())

jpeg("ExprCutoffPlot.jpg")
ggplot( dist.cov, aes(Distance, gene.cov) ) + theme_bw() + scale_x_continuous( limits=c(-mads$mads[1],mads$mads[1])) + scale_y_continuous( limits=c(0,mads$expr[19])) + geom_point(alpha=0.4, size =1.5 ) + geom_hline( aes(yintercept=exp.cutoff), linetype="dashed", colour ="red", size=1 ) + labs( x="Distance from predicted TSS to nearest Ensembl TSS", y="Per-base gene coverage" )
invisible(dev.off())

jpeg("DistanceCutoffPlot.jpg")
ggplot( dist.cov, aes(Distance, gene.cov) ) + theme_bw() + scale_x_continuous( limits=c(-mads$mads[1],mads$mads[1])) + scale_y_continuous( limits=c(0,mads$expr[19])) + geom_point(alpha=0.4, size =1.5 ) + geom_vline( aes(xintercept=-500), linetype="dashed", colour ="red", size=1 ) + geom_vline( aes(xintercept=500), linetype="dashed", colour ="red", size = 1 ) + labs( x="Distance from predicted TSS to nearest Ensembl TSS", y="Per-base gene coverage" )
invisible(dev.off())
