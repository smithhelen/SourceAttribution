#### Prepare test data and training data for pco method ####

# Load functions
source("methods/helpers.R")       # Functions used internally in the methods

# Function to map from variable levels to scores.

# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)
factor_to_pco_score <- function(var, dist, m, mp) {
  var_levels <- droplevels(var)
  d.train <- dist[levels(var_levels),levels(var_levels),drop=FALSE]
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  # use only non-zero eigenvalues
  nlambdas <- sum(eigen_B$values > epsilon)
  if (nlambdas == 0) {
    # No non-zero eigenvectors
    return(NULL)
  }
  # Restrict to a maximum of eigenvectors set by "m" or "mp (propG)" (default is mp=100% variation)
  ev <- eigen_B$values[seq_len(nlambdas)]
  lambdas_B <- filter_eigenvalues(ev, m=m, mp=mp) # restrict by axes or mp
  # Scale eigenvectors
  Qo <- eigen_B$vectors
  Q <- sweep(Qo[, seq_along(lambdas_B), drop=FALSE], 2, sqrt(abs(lambdas_B)), "*")
  score <- left_join(data.frame(Var_Level = var_levels), data.frame(Q) |> rownames_to_column("Var_Level"), by = "Var_Level") 
  output <- score |> select(-Var_Level)
  list(output = output,
       extra = list(d=dist, var_levels=levels(var_levels), diag.B.train=diag(B.train), Q=Q, lambdas_B=lambdas_B, 
                    dim = ncol(Q), suffix = colnames(Q), propG = cumsum(eigen_B$values)/sum(eigen_B$values)*100))
}

prepare_training_pco <- function(data, var_cols, class, d, m=NULL, mp=100, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data, all_of(var_cols))
  classes   <- data |> pull(class)
  # iterate over the var columns and distance matrices, and convert
  prepped <- map2(var_cols, d, factor_to_pco_score, m, mp) |> compact() # removes empties
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

predict_pco <- function(new.var_level, d, diag.B.train, Q, lambdas_B) {
  # new.var_level is length one
  d.new <- d[new.var_level, names(diag.B.train)]
  d.gower <- diag.B.train - (d.new ^ 2)
  newQ <- d.gower %*% Q / (2 * lambdas_B)
  colnames(newQ) = colnames(Q)
  new_var_level_score <- data.frame(Var_Level = new.var_level, newQ)
  new_var_level_score
}

# given var information and extra information, do the mapping. 
impute_score_pco <- function(var, extra) {
  var <- droplevels(var)
  var_levels <- pluck(extra, "var_levels")
  new.var_levels <- setdiff(levels(var), var_levels)
  d <- pluck(extra, "d")
  diag.B.train <- pluck(extra, "diag.B.train")
  Q <- pluck(extra, "Q")
  lambdas_B <- pluck(extra, "lambdas_B")
  new_scores <- map_df(new.var_levels, ~predict_pco(new.var_level = {.}, d, diag.B.train, Q, lambdas_B))
  var_level_score <- bind_rows(data.frame(Q) |> rownames_to_column("Var_Level"), new_scores)
  test_score <- data.frame(Var_Level = var) |> left_join(var_level_score, by = "Var_Level") |> select(-Var_Level) 
  list(test_score = test_score)
}

prepare_test_pco <- function(data, list_of_extras, id, residualised=NULL) {
  test_data <- data |> select(any_of(names(list_of_extras)))
  newdata_score <- map2(test_data, list_of_extras, impute_score_pco)
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- if(!is.null(residualised)) {
    bind_cols(data |> select(all_of(id), all_of(residualised)), newdata_pred)
  } else {
    bind_cols(data |> select(all_of(id)), newdata_pred)
  }
  newdata_pred
}
