############################################################################################
# Program Name: Partition_Agglomerative.R
# Purpose: unsupervised dimension-reduction with agglomerative approach to partition algorithm.
# Programmer: Joshua Millstein
# Date: 11/28/17

# Update of Partition_Agglomerative_funcs_v6.R to change mutual information method from mm to shrink

smry.mthd = function( tmpmat, method ){
	switch( method,
		ICC = ICC( tmpmat ),
		MI = MI( tmpmat ),
		minR2 = minR2( tmpmat ),
		PC1 = pc1( tmpmat ) )
} # End function smry.mthd

# return ICC and mean of raw variables, center and scale
ICC = function( tmpmat ){
  ICC_c(as.matrix(tmpmat))
	# outp = vector('list', 2 )
	# outp[[1]] = icc_c(tmpmat)
	# outp[[2]] = scale(rowMeans(tmpmat, na.rm = TRUE))
	# return( outp )
} # End ICC

# function to compute standardized mutual information and return mean
MI = function( tmpmat ){
	tmpmat1 = infotheo::discretize( tmpmat, disc="equalfreq" )
	myy = apply(tmpmat, 1, mean, na.rm=TRUE)
	myy.d = infotheo::discretize( myy, disc="equalfreq" )
	tmp = infotheo::mutinformation( tmpmat1, myy.d, method="shrink" ) # methods: emp mm shrink sg
	etpy = infotheo::entropy(cbind(tmpmat1, myy.d), method="shrink" )
	std.mi = tmp / etpy
	outp = vector('list', 2)
	outp[[1]] = std.mi   # standardized mutual information
	outp[[2]] = scale( myy )
	return( outp )
} # End MI

# return min R^2 and mean of raw variables, center and scale
minR2 = function( tmpmat ){
  minR2_c(as.matrix(tmpmat))
	# myyn = scale( apply(tmpmat, 1, mean, na.rm=TRUE) )
	# cors = cor( cbind(myyn, tmpmat), method="pearson", use="pairwise.complete.obs")
	# cors2 = cors^2
	# minR2 = min( cors2[ -1, "myyn" ] )
	# nms = rep( paste( "V", 1:nrow(tmpmat), sep="" ), ncol(tmpmat) )
	# outp = vector('list', 2 )
	# outp[[1]] = minR2
	# outp[[2]] = myyn
	# return( outp )
} # End minR2

# function to compute PC1 and return percent variance explained
pc1 = function( mymat ){
  pca_c(as.matrix(mymat))
	# rot = prcomp( mymat, retx = TRUE, center=TRUE, scale=TRUE, tol=sqrt(.Machine$double.eps) )
	# pct.var = rot$sdev[1]^2 / sum(rot$sdev^2)
	# outp = vector('list', 2 )
	# outp[[1]] = pct.var
	# outp[[2]] = scale( rot$x[, 1 ] ) 
	# return( outp )
} # End pc1

# correlation distance
r.dist.s = function(x1, x2) 1 - corr(x1, x2, spearman = TRUE)
r.dist.p = function(x1, x2) 1 - corr(x1, x2)

# function to update distance matrix with cluster variable replacing raw variables
updt.dist = function( dist.r, cluster.nm, clust.var.nms, dat.r, dist.type ){
 
	myvar = dat.r[, cluster.nm ]  
	nms = colnames(dat.r)[colnames(dat.r) %nin% cluster.nm]
	dat.r.tmp = as.data.frame( dat.r[, nms ] )
	names( dat.r.tmp ) = nms
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
	dist.r = dist.r[ rownames(dist.r) %nin% clust.var.nms, 
	                 colnames(dist.r) %nin% clust.var.nms ]
	if( !is.null(dim(dist.r)[1])){
	  
		dist.r = cbind( dist.r[ names(tmp.dist), ], tmp.dist )
		colnames( dist.r )[ncol(dist.r)] = cluster.nm
		dist.r = rbind( dist.r, rep(NA, ncol(dist.r)) )
		rownames( dist.r )[ nrow(dist.r) ] = cluster.nm
	}
	
	dist.r
} # End function updt.dist

assn.clustr = function( clust.vec, dist.r, dat.r, dat, pct.var, clusters, cluster.ind, method, dist.type ){
	success = FALSE
	clust.var.nms = colnames(dat.r[ , clust.vec ])
	mapping.index = clusters[,"cluster"] %in% colnames(dat.r[ , clust.vec ])
	clust.var.nms.raw = rownames(clusters)[mapping.index]
	cluster.nm = paste( "ReducedNewVar", cluster.ind, sep="" )
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


#' Agglomerative partitioning dimension reduction
#' 
#' An agglomerative partitioning strategy is used to reduce the dimensionality
#' of a dataset that has substantial dependencies among (at least some)
#' variables. An N x P dataset, with N samples and P variables is reduced to
#' and N x R dataset, where R < P. The amount of reduction, i.e., the dimension
#' of the reduced dataset, depends on a user specified information loss
#' constraint.
#' 
#' Clusters are formed and grown using an agglomerative approach, continuing
#' until the proportion of information explained by the summary variable for
#' each cluster drops below pct.var. The algorithm proceeds as follows,
#' 
#' \enumerate{
#'   \item compute a distance matrix, describing dissimilarity between all pairs
#'   of features
#'   \item join nearest two features into a candidate cluster/subset
#'   \item compute a summary variable for the candidate cluster
#'   \item check whether pct.var of information is captured in the new summary
#'   variable
#'   \item if pct.var of information is captured, then replace cluster features
#'   with the summary variable
#'   \item update distance matrix
#'   \item repeat steps 2-6 until the reduced dataset converges to a stable state
#' }
#' 
#' In step 4, the amount of information in the new summary variable is
#' compared against the set of original variable(s), not a summary variable.
#' The above general algorithm is implemented using a number of different
#' specific methods for evaluating proportion of information captured and
#' computing the summary variable,
#' 
#' \itemize{
#'   \item ICC.mn -- Intraclass correlation coefficient (ICC) is used to
#'   estimate the proportion of variance explained by the summary variable which
#'   is generated by taking the mean of values within a cluster.
#'   \item ICC.sm -- Same as ICC.mn, but data are transformed before the ICC
#'   computation by log(x + 1) and the summary variable is generated by summing
#'   non-transformed values within a cluster. This version is designed for count
#'   data.
#'   \item MI.mn -- Information
#' captured is evaluated by the standardized mutual information (MI), which is
#' computed using the infotheo package following discretization of original
#' variables into approximately nrow(mymat)^(1/3) bins. The summary variable is
#' the mean.  
#'   \item MI.sm -- Same as MI.mn, but the summary variable is generated by
#'   summing values within a cluster, designed for count data.
#'   \item minR2.mn -- Information captured is evaluated by the minimum r-squared
#' (Pearson correlation) between the summary variable and original variables.
#' Summary variable is the mean.  
#'   \item minR2.sm -- Same as minR2.mn but the summary variable is generated
#'   by summing values within a cluster and Spearman correlation is used to
#'   compute r-squared instead of Pearson. This version is designed for count
#'   data.
#'   \item PC1 -- The summary variable is the first principal component computed
#'   from variables within the cluster only. Information captured is evaluated
#'   by the percent variation of variables within the cluster explained by this
#'   summary variable.  
#'   \item PC1.log -- Same as PC1, but data are first log transformed after
#'   adding 1, log(x + 1), designed for count data.
#' }
#' 
#' @param mat Full dataset, NxP matrix or dataframe with N samples in rows
#' and P variables in columns.
#' @param pct.var Information constraint for reducing each candidate subset of
#' variables into a single new summary variable. pct.var must be a value in the
#' interval (0, 1]. For each candidate partition (subset of the variables),
#' pct.var specifies the minimum proportion of retained information allowable
#' in the summary variable from the total information contained in all
#' variables within that subset. Generally, the smaller this value, the fewer
#' the number of variables in the reduced set. It is possible to reduce the
#' number of variables without discarding any information at all (pct.var = 1)
#' if there are multiple variables that have the same values for each
#' observation.
#' @param method This argument determines the method used to compute the new
#' summary variable and evaluate the proportion of information captured. See
#' Details for descriptions of each method.
#' @param dist.type The agglomerative approach requires evaluating similarity
#' between variables. There are two methods currently implemented. "p"
#' specifies that similarity will be computed as, 1 - cor( mymat,
#' use="pairwise.complete.obs", method="pearson" ), whereas, "s" will be, 1 -
#' @param niter Number of iterations
#' cor( mymat, use="pairwise.complete.obs", method="spearman" ).
#' @return A list which includes the following columns:
#' \item{ mymat.r }{A dataset (dataframe) with reduced number of variables,
#' dimensions N x R.}
#' \item{clusters }{A dataset (dataframe) specifying the mapping between the
#' original variables and the reduced variables. Also included is a column,
#' pct.var, with estimates of the proportion of information of variables in a
#' cluster captured by the summary variable.}
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
#' rslts = partition( dat, pct.var=.8, method="PC1", dist.type="p" )
#' 
#' @export
partition = function( mat, pct.var, method, dist.type="p", niter = 1000){
	mymat = as.data.frame( mat )
	m = ncol(mymat)
	
	mymat.r = mymat 
	# reduced variable set mymat.r, sequentially replace variables with clusters if allowable
	
	# r-squared distance
	if(dist.type == "s"){
		mydist = 1 - corr(mymat, spearman = TRUE)
	} else {
		mydist = 1 - corr(mymat) 
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
	while ( min0 < niter ){
		# identify the closest pairs
		#print( paste( iter, Sys.time(), "Min0:", min0) )
		
		if ( sum(mydist.r, na.rm=TRUE) > 0 ) mymin = min(mydist.r, na.rm=TRUE) else break
		tmp.min = which( mydist.r <= mymin, arr.ind=TRUE )
		tmp.min = tmp.min[ 1, ]
		if ( mymin >= maxdist ) break
		
		if (dim(mydist.r)[1] < 2) break
		# assign clusters
		
		tmpout =  assign_clusters(tmp.min, as.matrix(mydist.r), as.matrix(mymat.r),
		                          as.matrix(mymat), pct.var, clusters, cluster.ind,
		                          method, dist.type, "ReducedNewVar")
	
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


#' Simulate block diagonal multivariate normal data
#' 
#' The MASS R package with function mvrnorm is used to simulate block diagonal
#' multivariate normal data, with correlations uniformly distributed within
#' blocks.
#' 
#' The MASS R package with function mvrnorm is used to simulate block diagonal
#' multivariate normal data, with correlations uniformly distributed within
#' blocks. Correlations within each block are randomly sampled from a U(c.lb,
#' c.ub) distribution prior to data generation.
#' 
#' @param blk.vec Vector of block sizes that determines the total number of
#' variables as well as number of variables within each block.
#' @param c.lb Within each block, variables are correlated according to values
#' sampled from a uniform distribution with lower bound c.lb.
#' @param c.ub Within each block, variables are correlated according to values
#' sampled from a uniform distribution with upper bound c.ub.
#' @param blk.names Character vector for column name prefix. Default is "X.mvn".
#' @param sep Character separating prefix and number in column names. Default is
#'  ".".
#' @param n Number of observations in the dataset.
#'
#' @return A dataframe with sum(blk.vec) columns and n rows.
#' @author Joshua Millstein
#' @references Millstein J, et al.
#' @examples
#' 
#' 
#' block_sizes <- 2:20
#' lower_bound <- .2
#' upper_bound <- .4
#' n <- 200
#' 
#' dat <- sim_blk_diag_mvn(block_sizes, lower_bound, upper_bound, n)
#' 
#' @export
sim_blk_diag_mvn <- function(blk.vec, c.lb, c.ub, n, blk.names = "X.mvn", 
                             sep = ".") {
	
	blks <- lapply(blk.vec, function(blk.size) {
    sigma <- matrix(runif(blk.size * blk.size, c.lb, c.ub), 
                    ncol = blk.size, 
                    nrow = blk.size)
    sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
    diag(sigma) <- 1
    MASS::mvrnorm(n, mu = rep(0, blk.size), Sigma = sigma)
	})
	
	sim_data  <- do.call(cbind.data.frame, blks)
	names(sim_data) <- paste(blk.names, seq_along(sim_data), sep = sep)
	
	sim_data
} # End function sim_blk_diag_mvn
