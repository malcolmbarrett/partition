##########################################################################
# Program Name: params-data.R
# Purpose: Dataset of estimated negative binomial zero inflation parameters from count data
#   from Dr. Wendy Cozen led study on microbiome and colon polyps
# Programer: Joshua Millstein
# Date: 3/19/18
#
#' Dataset of estimated negative binomial zero inflation parameters from count data
#'   from Dr. Wendy Cozen led study on microbiome and colon polyps
#'
#' Data parameters were generated from genus summarized 16s rRNA count data 
#' adjusted for read depth, excluding taxa with > 80 percent zeros.
#' Samples included > 20,000 reads
#' NBZIN parameters estimated using the zeroinfl() function of the pscl R package.
#'
#' @name params
#'
#' @docType data
#'
#' @usage data(params)
#'
#' @format A list object, with each element a vector of parameters including,
#'   NB intercept, NB log(theta), ZIN intercept, ZIN intercept p-value
#'
#' @keywords datasets
#'
#' @examples
#' data(params)
"params"