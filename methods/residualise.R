# function to calculate residualised distance matrices

# load libraries
library(fastDummies)
source("methods/recipe_cap.R")                     # CAP method

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

## Function to calculate matrices ######
d_resid <- function(dat_resid, d_Hamming, residualised){
  
  ## model, hat, and identity matrices
  # calculate relative frequencies
  p <- dat_resid |> group_by(across(all_of(residualised))) |> summarise(p=1/n())
  
  # convert residualised variable to dummy
  resid_dummy <- dummy_cols(dat_resid |> select(id, class, all_of({{residualised}})), select_columns={{residualised}}) |> arrange(across(all_of(residualised)))
  
  # substitute relative frequency values with dummy 1s
  contrast_fn <- function(dat, i){
    data.frame(dat[,i]-dat[,i+1]) |> set_names(colnames(dat)[i])}
  resid_dat <- resid_dummy |> select(-c(id, class, all_of({{residualised}}))) * rep(p$p, each=nrow(dat_resid))
  resid_contrasts <- cbind(
    id=resid_dummy |> pull(all_of(id)),
    map_dfc(1:(nrow(p)-1), ~contrast_fn(resid_dat, i=.))
  )
  
  # align rows of X with rows of B
  X <- resid_contrasts |> arrange(id) |> select(-id) |> as.matrix() 
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  I <- H |> ncol() |> diag()
  
  # residualise each distance matrix
  var_cols <- dat_resid |> select(starts_with("Var"))
  resid_distance_matrices <- map2(var_cols, d_Hamming, residualise, I, H)
  resid_distance_matrices <- map(resid_distance_matrices, as.matrix)
  
  # output
  return(resid_distance_matrices)
}




