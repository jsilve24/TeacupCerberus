#' Multinomial Dirichlet Sampler
#' @param Y count data (D x N) (e.g., taxa x samples)
#' @param alpha D-vector of Dirichlet Prior parameters
#' @param n_samples number of posterior samples to collect
#' @return Posterior samples as array of dimensions: D x N x n_samples
MultDirichletResample <- function(Y, alpha=rep(0.65, nrow(Y)), n_samples=2000){
  X <- array(0, dim=c(nrow(Y), ncol(Y), n_samples))
  for (i in 1:ncol(Y)){
    X[,i,] <- rDirichlet(Y[,i]+alpha, n_samples)
  }
  return(X)
}

#' Dirichlet Sampler
#' @param alpha D-vector of Dirichlet parameters
#' @param n_samples number of samples to produce
#' @return array of samples of dimensions: D x n+samples
rDirichlet <- function(alpha, n_samples=2000){
  X <- array(0, dim=c(length(alpha), n_samples))
  for (i in 1:length(alpha)){
    X[i,] <- rgamma(n_samples, alpha[i], scale=1)
  }
  n <- colSums(X)
  return(sweep(X, 2, n, FUN=`/`))
}


#' Convert Sigma from teacup_cerberus to Correlations
#'
#' Simple wrapper for cov2cor function that operates
#' on all posterior samples of Sigma
#'
#' @param Sigma Array of Covariance matricies from output of teacup_cerberus
#' @return Array of Correlation matricies of same dimensions
#' @export
tc2cor <- function(Sigma){
  if (length(dim(Sigma)) != 3) stop("Sigma must be an array of dimension 3")
  for (i in 1:dim(Sigma)[3]){
    Sigma[,,i] <- cov2cor(Sigma[,,i])
  }
  return(Sigma)
}


#' Filter Sigma "Hits" Based on ECDF
#'
#' Set a value of lambda  and this function will filter
#' to only show elements of Sigma(i,j) where i!=j that satisfy the following inequality:
#' ECDF(Sigma(i,j))(0) > lambda or  ECDF(Sigma(i,j))(0) < 1-lambda.
#' @param Sigma output of teacup_cerberus
#' @param lambda scalar between 0 and 1 (default = 0.95)
#' @param use_names if TRUE tries to name the hits rather than just giving indexes i,j
#' @return data.frame of "hits" ranked by abs(ECDF(Sigma(i,j))(0)-0.5), first row is then
#'  the "hit" we are most confident is non-zero and large in magnitude.
#' @export
#' @importFrom driver gather_array
#' @import dplyr
filter_sigma_ecdf <- function(Sigma, lambda=0.95, use_names=TRUE){
  Sigma.tidy <- driver::gather_array(Sigma, val, Coord1, Coord2, iter)

  if (use_names){
    if (!is.null(rownames(Sigma))) {
      Sigma.tidy <- dplyr::mutate(Sigma.tidy, Coord1 = rownames(Sigma)[Coord1])
    }
    if (!is.null(colnames(Sigma))) {
      Sigma.tidy <- dplyr::mutate(Sigma.tidy, Coord2 = colnames(Sigma)[Coord2])
    }
  }

  Sigma.tidy %>%
    group_by(Coord1, Coord2) %>%
    summarise(ECDF0 = ecdf(val)(0)) %>%
    filter((ECDF0 >= lambda) | (ECDF0 <= 1-lambda)) %>%
    filter(Coord1 != Coord2) %>%
    arrange(-abs(ECDF0-0.5))
}


#' Transform output of TeacupCerberus to iqlr
#' @param Sigma output of TeacupCerberus (in CLR basis)
#' @param D1 (integer) dimension of the first dataset
#' @param D2 (integer) dimension of the second dataset
#' @param qLow1 lower quantile defining iqlr for first dataset
#' @param qHigh1 upper quantile defining iqlr for first dataset
#' @param qLow2 lower quantile defining iqlr for second dataset
#' @param qHigh2 upper quantile defining iqlr for second
#' @return array of same dimensions as Sigma
#' @export
tc2iqlrvar <- function(Sigma, D1, D2,
                        qLow1=0.25, qHigh1 = 0.75,
                        qLow2=0.25, qHigh2 = 0.75){
  b <- 1:2
  Sigma <- RcppCoDA::vec_to_array(Sigma)
  s <- dim(Sigma)
  Sigma <- RcppCoDA:::array_pre(Sigma, b)
  Sigma <- tc2iqlrvar_internal(Sigma, as.integer(D1), as.integer(D2),
                                qLow1, qHigh1, qLow2, qHigh2)
  Sigma <- RcppCoDA:::array_post(Sigma, b, s)
  return(Sigma)
}

#' Transform output of Teacup Cerberus to phi Statistics
#'
#' @param Sigma output of teacup_cerberus in (clr coordinates)
#' @return array of the same dimensions as Sigma
#' @details this is just a wrapper around RcppCoDA::clrvar2phi
#' @export
tc2phi <- function(Sigma){
  return(RcppCoDA::clrvar2phi(Sigma))
}
