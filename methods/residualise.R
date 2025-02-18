## Load libraries and functions ##
source("methods/libs_fns.R")

# function to calculate residualised distance matrices
residualise <- function(var, dist, I, H) {
  d <- dist[paste0(var), paste0(var)]
  A <- -0.5 * d^2
  B <- dbl_center(A)  #B is same as Gowers matrix G
  B_res <- ((I - H) %*% B %*% (I - H))
  colnames(B_res) <- colnames(B)
  rownames(B_res) <- rownames(B)
  B_j <- matrix(diag(B_res), nrow=ncol(H), ncol=ncol(H))
  d_res <- sqrt(abs(t(B_j) - 2*B_res + B_j))
  d_res <- d_res[levels(var), levels(var)]
  d_res
}
