##########################################################################
# Program Name: Count_data_sumulate.R
# Purpose: Simulate microbiome data according to sample data, modeling negative binomial 
#	zero inflation count data
# Programer: Joshua Millstein
# Date: 3/19/18

# parameters computed in Microbiome_sim_v6.R based relative abundance microbiome data from 
#  "/Users/iTeams/Google\ Drive/Cozen/colonPolyps/data/Guoqin/RA_Genus.txt"
# excluded taxa with 0 reads in > 80% of samples

logit = function( myxb ){
	prob = 1 / (1 + exp( -myxb ))
	return( prob )
} # End function logit

#' Simulate NBZIN microbiome count data dependent on continuous covariate.
#' 
#' Data parameters were generated from genus summarized 16s rRNA count data 
#' adjusted for read depth, excluding taxa with > 80 percent zeros.
#' Samples included > 20,000 reads
#' NBZIN parameters in dataset params estimated using the zeroinfl() function of the pscl R package.
#' 
#' @param ntaxa Number of variables in the simulated microbiome dataset. Distribution
#' of these variables is negative binomial w/ zero inflation where the mean depends on 
#' xvec through beta.
#' @param xvec A vector of values that influence the mean of the simulated data through beta
#'  via a log-linear combination, that is, log(adjusted(mu)(i)) = mu.estimated + beta*xvec(i)*IQR(xvec),
#' for individual i, where IQR scales beta to the interquartile range of xvec.
#' @param beta Within each block, variables are correlated according to values
#' sampled from a uniform distribution with upper bound c.ub.
#' @param alt.ind Number of observations in the dataset.
#' 
#' @examples
#' data(params)
#' ntaxa = 5
#' xvec = rnorm( 200 )
#' beta = .1
#' sim.nb( ntaxa, xvec, beta, alt.ind=1:ntaxa ) # all alternative hypotheses
#' 
#' for( beta in seq(.5,2, length.out=10) ){
#' 	tmp = sim.nb( ntaxa, xvec, beta, alt.ind=1:ntaxa )
#' 	tmpc = cor( tmp, method="spearman" )[1,2:ntaxa]
#' 	print( paste("beta =", beta, "cor =", tmpc) )
#' 	for( j in 1:ntaxa ) print( paste( cor(xvec,tmp[,j],method="spearman") ) ) 
#' }
#' 
#' @export
sim.nb = function (ntaxa, xvec, beta=0, alt.ind=NULL) {
	ss = length(xvec)
	xvec1 = scale(xvec)
	# compute means conditional on covariates, simulate negbiom data
	mdat = as.data.frame(matrix( NA, nrow=0, ncol=ntaxa ) )
	
	for(tx in 1:ntaxa){
		# randomly sample parameter set for taxon from params
		ind = sample( 1:length(params), 1 )
		tmp = params[[ ind ]]
		mu.tx = tmp[ 1 ]
		theta.tx = exp( tmp[ 2 ] )
		xvec2 = xvec1 * mu.tx
		for( id in 1:ss){
			# if tx in fnullvec then adjust mean according to beta and xvec
			if( is.element( tx, alt.ind ) ){
				mu.adj = exp( mu.tx + beta*xvec1[ id ] )
				cnt = rnbinom( 1, size=theta.tx, mu=mu.adj )
			} else {          # End if fnullvec
				cnt = rnbinom( 1, size=theta.tx, mu=exp(mu.tx) )
			} # End else
			mdat[ id, tx ] = cnt
		} # End id loop
		
		## model zero inflation
		beta0z = tmp[ 3 ] - beta*xvec1 # the larger the xvec the smaller the prob of zero, assuming beta > 0
		probs = logit( beta0z )
		zeros = ifelse( probs > runif(ss), 0, 1 )
		mdat[ , tx ] = mdat[ , tx ] * zeros
	} # End tx loop
	
	return( mdat )
	
} # End function sim.nb

