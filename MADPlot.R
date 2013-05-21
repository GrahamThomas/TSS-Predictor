#!/usr/bin/env Rscript

library( ggplot2 )

results <- read.table("CutoffDistances.txt", header=T )
colnames(results) <- c(colnames(results)[1], sub( 'X', '', colnames(results)[-1] ))

mad.dist <- apply( results[,-1], 2, function(x)( mad(x, na.rm=T) ) )
mad.dist <- as.data.frame( as.matrix( mad.dist ) )

min.mad <- min( mad.dist$V1 ) 
bar.fill <- ifelse(( mad.dist$V1 == min.mad),"yes","no")
mad.dist <- cbind( mad.dist, bar.fill )

#jpeg("Optimizer-MADPlot.jpg")
pdf("Optimizer-MADPlot.pdf")
ggplot( mad.dist, aes(rownames(mad.dist), V1, fill=bar.fill ) ) + geom_bar( alpha = 0.85, stat="identity") + scale_fill_manual( values = c("darkgrey","red"), name="Minimum\nMAD" ) + coord_flip() + theme_bw() + labs( x="Cutoff Fraction", y="MAD between predicted \nand nearest annotated TSS for gene (bases)" )
dev.off()
