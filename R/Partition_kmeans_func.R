##############################################################################################
# Program Name: Partition_kmeans_func.R
# Purpose: Partition-conforming function using k-means, by contraining all clusters to have ICC > nu
# Programmer: Joshua Millstein
# Date: 10/25/17

# Function to compute ICC values from k-means results
km_icc = function( mymat, k ){
	nms = c("variable", "cluster", "icc")
	cluster.info = as.data.frame( matrix(NA, ncol=length(nms), nrow=ncol(mymat)) )
	names(cluster.info) = nms
	
	if( k < ncol(mymat) ){
		clstrs = kmeans( t(mymat), centers=k, nstart=25 )$cluster
	} else {
		dat.r = mymat
		cluster.info[, "variable"] = colnames( mymat )
		cluster.info[, "cluster"] = colnames( mymat )
		cluster.info[, "icc"] = 1
		outp = vector('list', 2 )
		outp[[1]] = dat.r
		outp[[2]] = cluster.info
		return( outp )
	}
	
	cluster.info[, "variable"] = names( clstrs )
	cluster.info[, "cluster"] = paste( "C", clstrs, sep="" )
	dat.r = as.data.frame( matrix(NA, ncol=k, nrow=nrow(mymat)) )
	names(dat.r) = paste( "C", 1:k, sep="" )
	
	for( c in 1:k ){
		c.nm = paste( "C", c, sep="" )
		tmpdat = as.matrix( mymat[ , is.element( clstrs, c ) ] )
		if( ncol(tmpdat) > 1 ){
			tmp.icc = ICC.mn( tmpdat )
			cluster.info[ is.element(cluster.info[,"cluster"], c.nm) , "icc"] = tmp.icc[[ 1 ]]
			dat.r[, c.nm ] = tmp.icc[[ 2 ]]
		} else {
			cluster.info[ is.element(cluster.info[,"cluster"], c.nm) , "icc"] = 1
			dat.r[, c.nm ] = tmpdat[, 1 ]
		}
	} # end c loop

	outp = vector('list', 2 )
	outp[[1]] = dat.r
	outp[[2]] = cluster.info
	return( outp )
} # end km_icc

# Function to compute K means for ICC threshold
#   Return k, summary variables, and respective ICC values
kmeans_icc = function( mymat, threshold.icc ){
	m = ncol(mymat)
	min.icc.m = rep(0, m)
	n.low = 0
	for( k in (m-1):1 ){
		tmp.km = km_icc( mymat, k )
		min.icc = min( tmp.km[[ 2 ]][ , "icc" ] )
		min.icc.m[k] = min.icc
		if( min.icc < threshold.icc ) n.low = n.low + 1
		if( n.low > 3 ){ 
			k.final = which( min.icc.m > threshold.icc )[1]
			if( is.na(k.final) ){
				kmdat = km_icc( mymat, m )
			} else {
				kmdat = km_icc( mymat, k.final )
			}
			break
		}
	} # end k loop
	
	return( kmdat )
} # End kmeans.vars







