##############################################################################################
# Program Name: Partition_kmeans_func.R
# Purpose: Partition-conforming function using k-means, by contraining all clusters to have ICC > nu
# Programmer: Joshua Millstein
# Date: 10/25/17

# Function to compute ICC values from k-means results
km_icc = function( mymat, k){
	clusters = as.data.frame( matrix( NA, nrow=ncol(mymat), ncol=2 ) )
	rownames(clusters) = colnames(mymat)
	names(clusters) = c("cluster", "pct.var")
	
	if( k < ncol(mymat) ){
		clstrs = quiet_kmeans(mymat, clusters = k)$cluster
	} else {
		dat.r = mymat
		clusters[, "cluster"] = colnames( mymat )
		clusters[, "pct.var"] = 1
		outp = vector('list', 2 )
		outp[[1]] = dat.r
		outp[[2]] = clusters
		return( outp )
	}
	
	rownames(clusters) = names( clstrs )
	clusters[, "cluster"] = paste( "ReducedNewVar", clstrs, sep="" )
	dat.r = as.data.frame( matrix(NA, ncol=k, nrow=nrow(mymat)) )
	names(dat.r) = paste( "ReducedNewVar", 1:k, sep="" )
	
	for( c in 1:k ){
	  
		c.nm = paste( "ReducedNewVar", c, sep="" )
		tmpdat = as.matrix( mymat[ , is.element( clstrs, c ) ] )
		if (ncol(tmpdat) == 0) browser()
		if( ncol(tmpdat) > 1 ){
			tmp.icc = ICC( tmpdat )
			clusters[ is.element(clusters[,"cluster"], c.nm) , "pct.var"] = tmp.icc[[ 1 ]]
			dat.r[, c.nm ] = tmp.icc[[ 2 ]]
		} else {
			clusters[ is.element(clusters[,"cluster"], c.nm) , "pct.var"] = 1
			dat.r[, c.nm ] = tmpdat[, 1 ]
		}
	} # end c loop

	outp = vector('list', 2 )
	outp[[1]] = dat.r
	outp[[2]] = clusters
	return( outp )
} # end km_icc

# Function to compute K means for ICC threshold
#   Return k, summary variables, and respective ICC values


#' K-means dimension reduction
#' 
#' Dimension reduction function using k-means to partition variables into
#' subsets, summarizing each by their means. The reduction is subject to an
#' information loss constraint specified by the user and measured by the
#' intraclass correlation (ICC), thus k is computed rather than provided as an
#' input.
#' 
#' The algorithm works by reducing K to its smallest value that satisfies
#' threshold.icc. That is, a value of K - 1, would yield at least one cluster
#' where the ICC comparing cluster variables to the mean would be smaller than
#' threshold.icc.
#' 
#' @param mymat Full dataset, NxP matrix or dataframe with N samples in rows
#' and P variables in columns.
#' @param threshold.icc This value (0 < threshold.icc < 1) specifies the maximum
#'   information loss for any cluster as measured by the ICC. Using k-means, the
#'   input dataset is reduced to an NxK dataset, where K is as small as possible
#'   subject the constraint. Reducing the threshold.icc yields a more permissive
#'   constraint, which tends to reduce the dimension K in the output dataset.
#' @return A list which includes the following columns: \item{ dat.r }{A dataset
#'   (dataframe) with reduced number of variables, dimensions N x R.}
#'   \item{cluster.info }{A dataset (dataframe) specifying the mapping between
#'   the original variables and the reduced variables. Also included is a
#'   column, icc, with ICC estimates for proportion of variance captured by the
#'   summary variable, the mean.}
#' @author Joshua Millstein
#' @references Millstein J, et al.
#' @keywords first principal component, intraclass correlation coefficient,
#' mutual information
#' @examples
#' 
#' 
#' blk.vec = 2:20
#' c.lb = .2
#' c.ub = .4
#' n = 200
#' 
#' dat = sim_blk_diag_mvn( blk.vec, c.lb, c.ub, n  )
#' 
#' rslts = kmeans_icc( dat, threshold.icc=.4 )
#' 
#' @export
kmeans_icc = function( mymat, threshold.icc ){

	m = ncol(mymat)
	min.icc.m = rep(0, m)
	n.low = 0
	for( k in rev(seq_along(mymat))[-1] ){
		tmp.km = km_icc( mymat, k )
		min.icc = min( tmp.km[[ 2 ]][ , "pct.var" ] )
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







