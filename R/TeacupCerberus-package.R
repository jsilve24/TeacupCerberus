#' TeacupCerberus
#'
#' Estimates covariation within and between two count datasets where the counts
#' contain multinomial variation (e.g., sequence count data like microbiome 16S or bulk/single-cell RNA-seq).
#' The model outputs Bayesian posterior samples over covariance matricies.
#' The entire posterior reflects uncertainty in the true covariation due to multinomial counting.
#'
#' @docType package
#' @author Justin D Silverman
#' @useDynLib TeacupCerberus
#' @import Rcpp
#' @import RcppEigen
#' @import RcppCoDA
#' @name RcppCoDA
NULL
