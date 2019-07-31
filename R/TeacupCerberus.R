#' Teacup Cerberus Model
#'
#' Estimates Covariance betweeen and withing two datasets accounting for
#' multinomial count uncertainty/error. Returns covariance with respect to
#' CLR coordinates, this can easily be converted to alternative representations.
#' See vignette for details.
#'
#' @param Y1 count data (D1 x N) (e.g., taxa x samples)
#' @param Y2 count data (D2 x N) (e.g., food x samples)
#' @param alpha1 D1-vector prior for Dirichlet for Y1 (think of it as a "pseudo-count" like thing,
#'   must be greater than zero)
#'   default: rep(1, D1)
#' @param alpha2 D2-vector prior for Dirichlet for Y2 (think of it as a "pseudo-count" like thing,
#'   must be greater than zero)
#'   default: rep(1, D2)
#' @details This fits the following model
#' \deqn{Y_1 ~ Multinomial(\pi_1)}
#' \deqn{Y_2 ~ Multinomial(\pi_1)}
#' \deqn{\pi_1 ~ Dirichlet(\alpha_1)}
#' \deqn{\pi_2 ~ Dirichlet(\alpha_2)}
#' and then transforming posterior samples of that model via
#' \deqn{\eta_1 = CLR^{-1}(\eta_1)}
#' \deqn{\eta_2 = CLR^{-1}(\eta_2)}
#' and then the s-th sample of Sigma (as a correlation matrix is given by)
#' \deqn{\Sigma^s = corr(cbind(\eta_1^s, \eta_2^s))}
#' @return Array Sigma of dimension (D1+D2) x (D1+D2) x n_samples
#'   (Sample of Covariance Matricies)
#' @import RcppCoDA
#' @export
teacup_cerberus <- function(Y1, Y2, alpha1 = NULL,
                            alpha2 = NULL,
                            n_samples=2000){
  D1 <- nrow(Y1);  D2 <- nrow(Y2); N <- ncol(Y1)
  if (is.null(alpha1)) alpha1 <- rep(1, D1)
  if (is.null(alpha2)) alpha2 <- rep(1, D2)
  idx1 <- 1:(D1-1)
  idx2 <- D1:(D1-1+D2-1)
  if (N != ncol(Y1)) stop("Y1 and Y2 must have the same number of columns")
  Eta <- array(0, dim=c(D1-1+D2-1, N, n_samples))
  Eta[idx1, ,] <-  alr(MultDirichletResample(Y1, alpha1, n_samples))
  Eta[idx2, ,] <-  alr(MultDirichletResample(Y2, alpha2, n_samples))

  Sigma <- array(0, dim=c(D1-1+D2-1,D1-1+D2-1, n_samples))
  for (i in 1:n_samples){
    Sigma[,,i] <- cov(t(Eta[,,i]))
  }

  # convert to CLR
  V1 <- driver::bdiag(alrContrast(D1, D1, inv=TRUE),
                      alrContrast(D2, D2, inv=TRUE))
  V2 <- driver::bdiag(clrContrast(D1, inv=FALSE),
                      clrContrast(D2, inv=FALSE))
  Sigma <- glrvar2glrvar(Sigma, V1, V2)

  if (is.null(rownames(Y1))) {
    n1 <- paste0("a", 1:D1)
  } else {
    n1 <- rownames(Y1)
  }
  if (is.null(rownames(Y2))) {
    n2 <- paste0("b", 1:D2)
  } else {
    n2 <- rownames(Y2)
  }
  n <- c(n1,n2)
  rownames(Sigma) <- colnames(Sigma) <- n

  return(Sigma)
}
