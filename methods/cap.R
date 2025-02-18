#### Prepare test data and training data for CAP method ####

# Load functions
source("methods/helpers.R")       # Functions used internally in the methods

epsilon <- sqrt(.Machine$double.eps)

# Function to map from variable levels to scores.

# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)

# ck is the number of k axes for cap
# cm is the number of m axes for cap
factor_to_CAP_score <- function(var, dist, class, ck, cm, cmp, c) {
  var_levels <- droplevels(var)
  d.train <- dist[levels(var_levels),levels(var_levels),drop=FALSE]
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)  #B is same as Gowers matrix G
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  # use only non-zero PCO eigenvalues
  nlambdas <- sum(eigen_B$values > epsilon) # this is not affected by axes (c) for CAP
  if (nlambdas == 0) {
    # No non-zero eigenvectors
    return(NULL)
  }
  ev <- eigen_B$values[seq_len(nlambdas)]
  lambdas_B <- filter_eigenvalues(ev, m=cm, mp=cmp) # subset PCO axes by cm or cmp (cm takes priority)
  Qo <- eigen_B$vectors[, seq_along(lambdas_B), drop=FALSE]  # note that this is different to the Q score in PCO method which is scaled by the sqrt(abs(lambdas_B))
  ct <- table(Var_level=var_levels, Class=class)
  H <- hat_fn(ct, k=ck) # restrict ct to ck axes, if ck is null then ck=ncol(ct)-1
  QHQ <- t(Qo) %*% H %*% Qo   # Combine Qo and H to get C_score
  eigen_QHQ <- eigen_decomp(QHQ, symmetric=TRUE)
  # select number of axes to retain - default is number of pco axes (NOTE CHANGE, used to be number of k axes)
  if(is.null(c)){c <- length(lambdas_B)}
  lambda_QHQ <- filter_eigenvalues(eigen_QHQ$values, m=c) # subset CAP axes for final RF variables
  U <- eigen_QHQ$vectors[,seq_along(lambda_QHQ),drop=FALSE]
  C_score <- Qo %*% U   # no need to scale by eigenvalue as only taking score, not doing ordination
  # Fill C_scores to individual isolates.  
  score <- left_join(data.frame(Var_level = var_levels), data.frame(C_score) |> rownames_to_column("Var_level"), by = "Var_level")
  output <- score |> dplyr::select(-Var_level)
  list(output = output,
       extra = list(d=dist, var_levels=levels(var_levels), diag.B.train=diag(B.train), lambdas_B=lambdas_B,  
                    U=U, Qo=Qo, C_score=C_score, dim = ncol(C_score), suffix = colnames(C_score)))
  }

prepare_training_cap <- function(data, vars, class, d, ck=NULL, cm=NULL, cmp=100, c=NULL, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data,all_of(vars))
  classes   <- data |> pull(class)
  # iterate over the var columns and distance matrices, and convert
  prepped <- map2(var_cols, d, factor_to_CAP_score, class = classes, ck, cm, cmp, c) |> compact() # removes empties
  output <- map(prepped,"output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select(all_of(class)) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  if(!is.null(residualised)) {
    prepped_data <- prepped_data |> bind_cols(data |> select(all_of(residualised)))
  }
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

predict_cap <- function(new.var_level, d, diag.B.train, Qo, lambdas_B, U) {
  # new.var_level is length one
  d.new <- d[new.var_level, names(diag.B.train)]
  d.gower <- diag.B.train - (d.new ^ 2)
  newQo <- d.gower %*% Qo / (2 * lambdas_B)
  colnames(newQo) = colnames(Qo)
  #for plotting, the canonical variable scores are standardized by the square root of their corresponding eigenvalue (lambda_QHQ). 
  #this is not necessary for the CAP method per se
  #newCscore <- newQo %*% U * sqrt(abs(lambda_QHQ))
  newCscore <- newQo %*% U
  new_var_score <- data.frame(Var_level = new.var_level, newCscore)
  new_var_score
}

# given variable information and a level map, do the mapping. 
impute_score_cap <- function(var, extra) {
  var_levels <- pluck(extra, "var_levels")
  var <- droplevels(var)
  new.var_levels <- setdiff(levels(var), var_levels)
  d <- pluck(extra, "d")
  diag.B.train <- pluck(extra, "diag.B.train")
  lambdas_B <- pluck(extra, "lambdas_B")
  Qo <- pluck(extra, "Qo")
  C_score <- pluck(extra, "C_score")
  U <- pluck(extra, "U")
  new_scores <- map_df(new.var_levels, ~predict_cap(new.var_level = {.}, d, diag.B.train, Qo, lambdas_B, U))
  var_level_score <- bind_rows(data.frame(C_score) |> rownames_to_column("Var_level"), new_scores)
  test_score = data.frame(Var_level = var) |> left_join(var_level_score, by = "Var_level") |> select(-Var_level)
  list(test_score = test_score)
}

prepare_test_cap <- function(data, extra, id, residualised=NULL) {
  id <- data |> select(all_of(id))
  var_cols <- data |> select(any_of(names(extra)))
  newdata_score <- map2(var_cols, extra, impute_score_cap)
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- if(!is.null(residualised)) {
    bind_cols(id, data |> select(all_of(residualised)), newdata_pred)
  } else {
      bind_cols(id, newdata_pred)
  }
  newdata_pred
}

