############################################################################################
# Program Name: Partition_Agglomerative.R
# Purpose: unsupervised dimension-reduction with agglomerative approach to partition algorithm.
# Programmer: Joshua Millstein
# Date: 11/28/17

# Update of Partition_Agglomerative_funcs_v6.R to change mutual information method from mm to shrink

smry.mthd = function( tmpmat, method ){
	switch( method,
		ICC.mn = ICC.mn( tmpmat ),
		ICC.sm = ICC.log.sm( tmpmat ),
		MI.mn = MI.mn( tmpmat ),
		MI.sm = MI.sm( tmpmat ),
		minR2.mn = minR2.mn( tmpmat ),
		minR2.sm = minR2.sm( tmpmat ),
		PC1 = pc1( tmpmat ),
		PC1.log = pc1.log( tmpmat ) )
} # End function smry.mthd

# function to compute ICC after log( +1) transformation
# return ICC and sum of raw variables
ICC.log.sm = function( tmpmat ){
	dat.vec = log(as.numeric( as.matrix(tmpmat) )+1)
	nms = rep( paste( "V", 1:nrow(tmpmat), sep="" ), ncol(tmpmat) )
	outp = vector('list', 2 )
	outp[[1]] = ICC::ICCbareF( nms, dat.vec ) 
	outp[[2]] = apply(tmpmat, 1, sum, na.rm=TRUE)
	return( outp )
} # End ICC.log.sm

# return ICC and mean of raw variables, center and scale
ICC.mn = function( tmpmat ){
	dat.vec = as.numeric( as.matrix(tmpmat) )
	nms = rep( paste( "V", 1:nrow(tmpmat), sep="" ), ncol(tmpmat) )
	outp = vector('list', 2 )
	outp[[1]] = ICC::ICCbareF( nms, dat.vec ) 
	outp[[2]] = apply(tmpmat, 1, mean, na.rm=TRUE)
	return( outp )
} # End ICC.mn

# function to compute standardized mutual information and return mean
MI.mn = function( tmpmat ){
	tmpmat1 = infotheo::discretize( tmpmat, disc="equalfreq" )
	myy = apply(tmpmat, 1, mean, na.rm=TRUE)
	myy.d = infotheo::discretize( myy, disc="equalfreq" )
	tmp = infotheo::mutinformation( tmpmat1, myy.d, method="shrink" ) # methods: emp mm shrink sg
	etpy = infotheo::entropy(cbind(tmpmat1, myy.d), method="shrink" )
	std.mi = tmp / etpy
	outp = vector('list', 2 )
	outp[[1]] = std.mi   # standardized mutual information
	outp[[2]] = apply(tmpmat, 1, mean, na.rm=TRUE)
	return( outp )
} # End MI.sm

# function to compute standardized mutual information and return sum
MI.sm = function( tmpmat ){
	tmpmat1 = infotheo::discretize( tmpmat, disc="equalfreq" )
	myy = apply(tmpmat, 1, sum, na.rm=TRUE)
	myy.d = infotheo::discretize( myy, disc="equalfreq" )
	tmp = infotheo::mutinformation( tmpmat1, myy.d, method="shrink" ) # methods: emp mm shrink sg
	etpy = infotheo::entropy(cbind(tmpmat1, myy.d), method="shrink" )
	std.mi = tmp / etpy
	outp = vector('list', 2 )
	outp[[1]] = std.mi   # standardized mutual information
	outp[[2]] = apply(tmpmat, 1, sum, na.rm=TRUE)
	return( outp )
} # End MI.sm

# return min R^2 and mean of raw variables, center and scale
minR2.mn = function( tmpmat ){
	myy = apply(tmpmat, 1, mean, na.rm=TRUE)
	myyn = (myy - mean(myy, na.rm=TRUE)) / sd(myy, na.rm=TRUE) # center and scale myy
	cors = cor( cbind(myyn, tmpmat), method="pearson", use="pairwise.complete.obs")
	cors2 = cors^2
	minR2 = min( cors2[ -1, "myyn" ] )
	nms = rep( paste( "V", 1:nrow(tmpmat), sep="" ), ncol(tmpmat) )
	outp = vector('list', 2 )
	outp[[1]] = minR2
	outp[[2]] = myyn
	return( outp )
} # End minR2.mn

# return min Spearman R^2 and sum of raw variables
minR2.sm = function( tmpmat ){
	myy = apply(tmpmat, 1, sum, na.rm=TRUE)
	cors = cor( cbind(myy, tmpmat), method="spearman", use="pairwise.complete.obs")
	cors2 = cors^2
	minR2 = min( cors2[ -1, "myy" ] )
	nms = rep( paste( "V", 1:nrow(tmpmat), sep="" ), ncol(tmpmat) )
	outp = vector('list', 2 )
	outp[[1]] = minR2
	outp[[2]] = myy
	return( outp )
} # End minR2.sm

# function to compute PC1 and return percent variance explained
pc1 = function( mymat ){
	rot = prcomp( mymat, retx = TRUE, center=TRUE, scale=TRUE, tol=sqrt(.Machine$double.eps) )
	pct.var = rot$sdev[1]^2 / sum(rot$sdev^2)
	mypc1 = rot$x[, 1 ]
	mypc1 = (mypc1 - mean(mypc1, na.rm=TRUE)) / sd(mypc1, na.rm=TRUE) # center and scale PC1
	outp = vector('list', 2 )
	outp[[1]] = pct.var
	outp[[2]] = mypc1
	return( outp )
}

# function to compute PC1 and return percent variance explained
pc1.log = function( mymat ){
	mymat.log = log(as.matrix(mymat) + 1)
	rot = prcomp( mymat.log, retx = TRUE, center=TRUE, scale=TRUE, tol=sqrt(.Machine$double.eps) )
	pct.var = rot$sdev[1]^2 / sum(rot$sdev^2)
	mypc1 = rot$x[, 1 ]
	mypc1 = (mypc1 - mean(mypc1, na.rm=TRUE)) / sd(mypc1, na.rm=TRUE) # center and scale PC1
	outp = vector('list', 2 )
	outp[[1]] = pct.var
	outp[[2]] = mypc1
	return( outp )
}

# correlation distance
r.dist.s = function(x1, x2) 1 - cor(x1, x2, method="spearman", use="pairwise.complete.obs")
r.dist.p = function(x1, x2) 1 - cor(x1, x2, method="pearson", use="pairwise.complete.obs")

# function to update distance matrix with cluster variable replacing raw variables
updt.dist = function( dist.r, cluster.nm, clust.var.nms, dat.r, dist.type ){
	myvar = dat.r[, cluster.nm ]  
	nms1 = colnames(dat.r)[ !is.element( colnames(dat.r), cluster.nm ) ]
	dat.r.tmp = as.data.frame( dat.r[, nms1 ] )
	names( dat.r.tmp ) = nms1
	# compute distances between myvar and other variables
	if( !is.na(dim(dat.r.tmp)[1]) ){
		if(dist.type == "s"){
			tmp.dist = apply( dat.r.tmp, 2, r.dist.s, x2=myvar )
		} else {
			tmp.dist = apply( dat.r.tmp, 2, r.dist.p, x2=myvar )
		}
		names(tmp.dist) = colnames(dat.r.tmp)
	}
	
	# remove cluster component variables from distance matrix
	dist.r = dist.r[ !is.element( rownames(dist.r), clust.var.nms ), !is.element( colnames(dist.r), clust.var.nms ) ]
	if( !is.null(dim(dist.r)[1])){
		dist.r = cbind( dist.r[ names(tmp.dist), ], tmp.dist )
		colnames( dist.r )[ncol(dist.r)] = cluster.nm
		dist.r = rbind( dist.r, rep(NA, ncol(dist.r)) )
		rownames( dist.r )[ nrow(dist.r) ] = cluster.nm
	}
	return( dist.r )
} # End function updt.dist

assn.clustr = function( clust.vec, dist.r, dat.r, dat, pct.var, clusters, cluster.ind, method, dist.type ){
	success = FALSE
	clust.var.nms = colnames(dat.r[ , clust.vec ])
	clust.var.nms.raw = rownames(clusters)[ is.element( clusters[,"cluster"], colnames(dat.r[ , clust.vec ]) ) ]
	cluster.nm = paste( "C", cluster.ind, sep="" )
	
	# 2. compute summary and test against pct.var
	var.ind = which( is.element( colnames(dat), clust.var.nms.raw ) )
	tmp.svec = smry.mthd( dat[, var.ind], method )

	# 3. update, dist.r, dat.r, clusters, cluster.new.names
	if( tmp.svec[[ 1 ]] >= pct.var ){ 
			success = TRUE
			clusters[ clust.var.nms.raw, "cluster" ] = cluster.nm
			clusters[ clust.var.nms.raw, "pct.var" ] = tmp.svec[[ 1 ]]
			cluster.ind = cluster.ind + 1                                     # increment to form unique name for next formed cluster
			dat.r = as.data.frame( cbind( dat.r, tmp.svec[[ 2 ]] ) ) # add summary variable to dat.r
			names(dat.r)[ ncol(dat.r) ] = cluster.nm
			dat.r = dat.r[ , -clust.vec]
			aa = !is.null(dim(dat.r)[1])
			if( aa ) dist.r = updt.dist( dist.r, cluster.nm, clust.var.nms, dat.r, dist.type )
	} else {   # set distance to NA to avoid testing again
		dist.r[ clust.vec[ 1 ], clust.vec[ 2 ] ] = NA
	}
	
	outp = vector('list', 5 )
	outp[[1]] = clusters
	outp[[2]] = cluster.ind
	outp[[3]] = dist.r
	outp[[4]] = dat.r
	outp[[5]] = success
	return( outp )
} # End assn.clustr

## Main function
partition = function( mymat, pct.var, method, dist.type="p" ){
	mymat = as.data.frame( mymat )
	m = ncol(mymat)
	
	mymat.r = mymat 
	# reduced variable set mymat.r, sequentially replace variables with clusters if allowable
	
	# r-squared distance
	if(dist.type == "s"){
		mydist = 1 - cor( mymat, use="pairwise.complete.obs", method="spearman" ) 
	} else {
		mydist = 1 - cor( mymat, use="pairwise.complete.obs", method="pearson" ) 
	}
	mydist[ lower.tri(mydist, diag=TRUE) ] = NA
	mydist.r = mydist                                            # distance matrix to be updated and reduced as clusters are formed
	maxdist = max( mydist, na.rm=TRUE )
	
	# update cluster membership to this dataframe
	clusters = as.data.frame( matrix( NA, nrow=m, ncol=2 ) )
	rownames(clusters) = colnames(mymat)
	names(clusters) = c("cluster", "pct.var")
	clusters[, "cluster"] = colnames(mymat)
	clusters[, "pct.var" ] = 1
	cluster.ind = 1 # increment indicator for cluster name every time an allowable cluster is formed
	
	min0 = 0 # minimum distance with no qualifying pairs
	while( min0 < 1000 ){
		# identify the closest pairs
		#print( paste( iter, Sys.time(), "Min0:", min0) )
		
		if( sum(mydist.r, na.rm=TRUE) > 0 ) mymin = min(mydist.r, na.rm=TRUE) else break
		tmp.min = which( mydist.r <= mymin, arr.ind=TRUE )
		tmp.min = tmp.min[ 1, ]
		if( mymin >= maxdist ) break
		
		if(dim(mydist.r)[1] < 2) break
		# assign clusters
		tmpout =  assn.clustr( tmp.min, mydist.r, mymat.r, mymat, pct.var, clusters, cluster.ind, method, dist.type )
		
		clusters = tmpout[[ 1 ]]
		cluster.ind = tmpout[[ 2 ]]
		mydist.r = tmpout[[ 3 ]]
		mymat.r = tmpout[[ 4 ]]
		if( tmpout[[ 5 ]] ) min0 = 0 else min0 = min0 + 1
	} # End while tmp.min
	
	outp = vector('list', 2 )
	outp[[ 1 ]] = mymat.r
	outp[[ 2 ]] = clusters 
	return( outp )
} # End function partition

## Function to simulate block diagonal correlated data based on blocks of mvn data
sim_blk_diag_mvn = function( blk.vec, c.lb, c.ub, n  ){
	blk.num = length(blk.vec)
	ydat = as.data.frame( matrix( NA, nrow=n, ncol=blk.num ) )
	for( blk in 1:blk.num ){
		blk.size = blk.vec[ blk ]
		sigma = matrix( runif(blk.size*blk.size, c.lb, c.ub), ncol=blk.size, nrow=blk.size )
		sigma[lower.tri(sigma)] = t(sigma)[lower.tri(sigma)]
		diag(sigma) = 1
		tmp = MASS::mvrnorm( n, mu=rep(0,blk.size), Sigma=sigma )
		if( blk == 1 ) xdat = tmp else xdat = cbind( xdat, tmp )
	} # End blk loop
	
	xdat = as.data.frame( xdat )
	names( xdat ) = paste( "X.mvn", 1:ncol(xdat), sep="." )
	return( xdat )
} # End function sim_blk_diag_mvn
