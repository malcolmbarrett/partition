\name{sim_blk_diag_mvn}
\alias{sim_blk_diag_mvn}
\title{
  Simulate block diagonal multivariate normal data
}
\description{
   The MASS R package with function mvrnorm is used to simulate block diagonal multivariate normal data, with correlations uniformly distributed within blocks. 
}
\usage{
   sim_blk_diag_mvn( blk.vec, c.lb, c.ub, n  )
}
\arguments{
  \item{blk.vec}{
     Vector of block sizes that determines the total number of variables as well as number of variables within each block. 
}
  \item{c.lb}{
     Within each block, variables are correlated according to values sampled from a uniform distribution with lower bound c.lb.
}
  \item{c.ub}{
     Within each block, variables are correlated according to values sampled from a uniform distribution with upper bound c.ub.
}
  \item{n}{
     Number of observations in the dataset.
}
}
\details{
  The MASS R package with function mvrnorm is used to simulate block diagonal multivariate normal data, with correlations uniformly distributed within blocks. Correlations within each block are randomly sampled from a U(c.lb, c.ub) distribution prior to data generation.
}
\value{
  A dataframe with sum(blk.vec) columns and n rows.
}
\references{
 Millstein J, et al.
}
\author{
  Joshua Millstein
}

\examples{

blk.vec = 2:20
c.lb = .2
c.ub = .4
n = 200

dat = sim_blk_diag_mvn( blk.vec, c.lb, c.ub, n  )

}


