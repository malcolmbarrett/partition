##############################################################################################
# Program Name: Partition_plotFunc_StackedArea.R
# Purpose: Plot clustering as a function of minimum information captured.
# Programmer: Joshua Millstein
# Date: 12/01/17

plot_dr = function( dat, method, dist.type="p", var.vec = seq(.1, .5, length.out=10), out.file=NA, tol.missing=.2, impt=TRUE ){
	# pct.var = minimum variance captured in cluster
	# clst.num = number of clusters within count rng
	# clst.cnt.rng = range of number of elements within cluster, e.g., ">5"
	
	# remove rows with more than (tol.missing*100)% missing
	rind = rep(TRUE, nrow(dat))
	for(i in 1:nrow(dat)) if( sum(is.na(dat[i,])) > (tol.missing*ncol(dat)) ) rind[ i ] = FALSE
	dat = dat[ rind, ]
	
	# remove columns with more than (tol.missing*100)% missing
	cind = rep(TRUE, ncol(dat))
	for(j in 1:ncol(dat)) if( sum(is.na(dat[ , j])) > (tol.missing*nrow(dat)) ) cind[ i ] = FALSE
	dat = dat[ , cind]
	
	# impute missing with pamr
	if(impt){
		if( sum(is.na(as.numeric(as.matrix(dat)))) > .1 ){
			tmp = list(x = t(dat), y = rep(1, nrow(dat)) )
			dat = as.data.frame( t( pamr::pamr.knnimpute( tmp )$x ) )
		}
	} # end impt check
	
	c.rng = c("1", "2-4", "5-8", "9-16", "17-32", ">32")
	breaks = c(0,1,4,8,16,32,10^10)
	nms = c("pct.var.target", "pct.var.observed", "clst.num", "clst.cnt.rng" )
	dat.plot = as.data.frame( matrix(NA, nrow=length(var.vec)*length(c.rng), ncol=length(nms)) )
	names(dat.plot) = nms
	dat.plot[, "clst.cnt.rng"] = factor( rep(c.rng, length(var.vec)), levels=c.rng )
	dat.plot[, "pct.var.target"] = rep( var.vec, each=length(c.rng) )

	for( pct.var in var.vec ){
		print( paste( "observed",  Sys.time(), " pct.var = ", pct.var ) )
		tmp = partition( dat, pct.var,  method, dist.type )
		clstr = tmp[[2]]
		dat.plot[ is.element(dat.plot$pct.var.target, pct.var), "clst.num"] = hist(table(clstr[,"cluster"]),breaks=breaks,plot=FALSE)$counts
		dat.plot[ is.element(dat.plot$pct.var.target, pct.var), "pct.var.observed"] = mean( clstr[,"pct.var"], na.rm=TRUE )
		print(paste("Percent variance explained:",pct.var))
	}

	dat.plot[, "pct.var.observed1"] = round(dat.plot[, "pct.var.observed"] * 100)
	dat.plot[, "pct.var.target1"] = round(dat.plot[, "pct.var.target"] * 100)
	
	lab.dat = dat.plot[ !duplicated( dat.plot$pct.var.target1 ), ]
	lab.dat = lab.dat[ order( lab.dat$pct.var.target1 ), c("pct.var.target1", "pct.var.observed1") ]
	mypos = c( min(dat.plot[, "pct.var.target1"]), max(dat.plot[, "pct.var.target1"]), 
					mean( c(min(dat.plot[, "pct.var.target1"]), max(dat.plot[, "pct.var.target1"])) ) )
	
	pct.var.target1 <- clst.num <- clst.cnt.rng <- NULL
	p1 = ggplot2::ggplot(dat.plot, aes(x=pct.var.target1, y=clst.num, fill=clst.cnt.rng)) + theme_bw() +
		geom_area(colour="black", size=.2, alpha=.4) +
		scale_fill_brewer(palette="Set1") + scale_x_continuous(limits = mypos[1:2], 
		breaks=lab.dat$pct.var.target1,
		sec.axis = sec_axis(~ ., breaks=lab.dat$pct.var.target1, labels=lab.dat$pct.var.observed1)) +
		xlab("minimum variance explained (%)") + ylab("k") + labs(title="Observed") +
		theme(legend.position = c(0.1, .65)) +
		guides(fill=guide_legend(title="Partition Size")) +
		annotate("text", x=mypos[3], y=1.05*nrow(clstr), label="variance captured (%)")
    
	dat.plot.obs = dat.plot    
    
	# Recreate graph after randomizing data
	dat1 = dat
	dat.plot1 = dat.plot
	for(j in 1:ncol(dat1)) dat1[,j] = dat1[ sample(1:nrow(dat1), nrow(dat1)), j]

	for( pct.var in var.vec ){
		print( paste( "permuted",  Sys.time(), " pct.var = ", pct.var ) )
		tmp = partition( dat, pct.var,  method, dist.type )
		clstr = tmp[[2]]
		dat.plot[ is.element(dat.plot$pct.var.target, pct.var), "clst.num"] = hist(table(clstr[,"cluster"]),breaks=breaks,plot=FALSE)$counts
		dat.plot[ is.element(dat.plot$pct.var.target, pct.var), "pct.var.observed"] = mean( clstr[,"pct.var"], na.rm=TRUE )
		print(paste("Percent variance explained:",pct.var))
	}

	dat.plot[, "pct.var.observed1"] = round(dat.plot[, "pct.var.observed"] * 100)
	dat.plot[, "pct.var.target1"] = round(dat.plot[, "pct.var.target"] * 100)
	
	lab.dat = dat.plot[ !duplicated( dat.plot$pct.var.target1 ), ]
	lab.dat = lab.dat[ order( lab.dat$pct.var.target1 ), c("pct.var.target1", "pct.var.observed1") ]
	mypos = c( min(dat.plot[, "pct.var.target1"]), max(dat.plot[, "pct.var.target1"]), 
					mean( c(min(dat.plot[, "pct.var.target1"]), max(dat.plot[, "pct.var.target1"])) ) )

	p2 = ggplot2::ggplot(dat.plot, aes(x=pct.var.target1, y=clst.num, fill=clst.cnt.rng)) + theme_bw() +
		geom_area(colour="black", size=.2, alpha=.4) +
		scale_fill_brewer(palette="Set1") + scale_x_continuous(limits = mypos[1:2], 
		breaks=lab.dat$pct.var.target1,
		sec.axis = sec_axis(~ ., breaks=lab.dat$pct.var.target1, labels=lab.dat$pct.var.observed1)) +
		xlab("minimum variance explained (%)") + ylab("k") + labs(title="Randomized") +
		theme(legend.position = c(0.1, .65)) +
		guides(fill=guide_legend(title="Cluster Size")) +
		annotate("text", x=mypos[3], y=1.05*nrow(clstr), label="variance captured (%)")
	
	dat.plot.perm = dat.plot 

	if( is.na(out.file) ){
		gridExtra::grid.arrange( p1, p2, ncol=1)
	} else {
		pdf( out.file )
			gridExtra::grid.arrange( p1, p2, ncol=1)
		dev.off()
	}
	
	outp = vector('list', 2)
	outp[[ 1 ]] = dat.plot.obs
	outp[[ 2 ]] = dat.plot.perm
	return( outp )
	
} # End plot_dr



