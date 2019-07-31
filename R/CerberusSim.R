#' Simulate from Cerberus Model
#'
#' @param D1 Number of multinomial categories in first dataset
#' @param D2 Number of multinomial categories in second dataset
#' @param N number of samples
#' @param Q number of covariates
#' @param true_priors should Xi and upsilon be chosen to have mean at true
#'   simulated value
#' @param use_names should outputs have named dimensions?
#' @return list
#'
#' @export
#' @importFrom RcppCoDA alrInv
#' @importFrom stats rnorm rmultinom
#'
#' @examples
#' sim <- cerberus_sim()
cerberus_sim <- function(D1=5, D2=5, N=10, Q=2, true_priors=TRUE, use_names=TRUE){

  # Simulate Data
  Sigma <- diag(sample(1:8, D1-1+D2-1, replace=TRUE))
  Sigma[2, 3] <- Sigma[3,2] <- -1
  Gamma <- diag(sqrt(rnorm(Q)^2))
  Theta <- matrix(0, D1-1+D2-1, Q)
  Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D1-1+D2-1)), nrow=D1-1+D2-1)%*%chol(Gamma)
  X <- matrix(rnorm(N*(Q-1)), Q-1, N)
  X <- rbind(1, X)
  Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D1-1+D2-1)), nrow=D1-1+D2-1)
  Eta2 <- Eta[D1:(D1-1+D2-1),]
  Eta1 <- Eta[1:(D1-1),]
  Pi1 <- alrInv(Eta1)
  Pi2 <- alrInv(Eta2)
  Y1 <- matrix(0, D1, N)
  Y2 <- matrix(0, D2, N)
  for (i in 1:N) {
    Y1[,i] <- rmultinom(1, sample(5000:10000), prob = Pi1[,i])
    Y2[,i] <- rmultinom(1, sample(5000:10000), prob = Pi2[,i])
  }
  if (use_names){
    colnames(X) <- colnames(Y1) <- colnames(Y2) <- paste0("s", 1:N)
    rownames(Y1) <- paste0("c", 1:D1)
    rownames(Y2) <- paste0("d", 1:D2)
    rownames(X) <- paste0("x", 1:Q)
  }

  # Priors
  if (true_priors){
    upsilon <- D1+D2+10
    Xi <- Sigma*(upsilon-(D1-1+D1-1)-1)
  } else {
    upsilon <- D
    Xi <- diag(D1-1+D2-1)
  }

  # Precompute
  KInv <- solve(Xi)
  AInv <- solve(diag(N)+ t(X)%*%Gamma%*%X)

  return(list(Sigma=Sigma, Gamma=Gamma, D1=D1, D2=D2, N=N, Q=Q, Theta=Theta, Phi=Phi,
              X=X, Y1=Y1, Y2=Y2, Eta=Eta, upsilon=upsilon, Xi=Xi, KInv=KInv, AInv=AInv))
}
