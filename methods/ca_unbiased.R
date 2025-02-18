#### Prepare test data and training data for ca_unbiased method ####

# Function to map from variable levels to scores
epsilon <- sqrt(.Machine$double.eps)

# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)
factor_to_ca0_score <- function(var, class, k) {
  var_levels <- droplevels(var)
  if(nlevels(var_levels) < 2){
    # if we only have one var_level we can't do anything
    return(NULL)      }
  ct <- table(Var_Level=var_levels, Class=class)
  if(!length(k)){k <- ncol(ct)-1}
  P <- ct/rowSums(ct)
  S <- cov.wt(P, wt = rowSums(ct))$cov
  
  eigen_S <- eigen_decomp(S, symmetric=TRUE) ## PCA of weighted covariance matrix of class probabilites
  # Restrict to a maximum of eigenvectors set by "k" (the number of axes) (default is NULL = ncol(ct)-1)
  nlambdas <- min(sum(eigen_S$values > epsilon), k)
  # principal components
  pc <- eigen_S$vectors
  X <- P %*% pc[, seq_len(nlambdas), drop=FALSE] |> as.data.frame() 
  score <- left_join(data.frame(Var_Level = var_levels), X |> rownames_to_column("Var_Level"), by = "Var_Level") 
  output <- score |> select(-Var_Level)
  list(output = output,
       extra = list(X=X, k=k, var_levels=levels(var_levels), dim = ncol(X), suffix = colnames(X)))
}

# iterate over the variable columns and grab the output as a new data.frame to send into ca, and store the absent level stuff for later
prepare_training_ca0 <- function(data, var_cols, class, k=NULL, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- data |> select(all_of(var_cols))
  classes   <- data |> pull(all_of(class))
  # iterate over the var columns, and convert
  prepped <- map(var_cols, factor_to_ca0_score, class = classes, k = k) |> compact()
  output <- map(prepped, "output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select(all_of(class)) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  if(!is.null(residualised)) {
    prepped_data <- prepped_data |> bind_cols(data |> select(all_of(residualised)))
  }
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

# score absent levels as zero
predict_ca0 <- function(new.var_level, X) {
  newX <- matrix(NA, nrow = 1, ncol = ncol(X)) |> as.data.frame()
  colnames(newX) <- colnames(X)
  new_var_level_score <- data.frame(Var_Level = new.var_level, newX)
  new_var_level_score
}

# given variable information and a level map, do the mapping. Absent levels are mapped to 'new'
impute_score_ca0 <- function(var, extra) {
  var <- droplevels(var)
  var_levels <- pluck(extra, "var_levels")
  new.var_levels <- setdiff(levels(var), var_levels)
  X <- pluck(extra, "X")
  
  new_var_level_score <- map_df(new.var_levels, ~predict_ca0(new.var_level = {.}, X))
  var_level_score <- bind_rows(data.frame(X) |> rownames_to_column("Var_Level"), new_var_level_score)
  test_score <- data.frame(Var_Level = var) |> left_join(var_level_score, by = "Var_Level") |> select(-Var_Level) 
  
  # add noise to imputed zero scores
  #set.seed(123) # only for confirming that this method and the recipe method are the same
  r <- min(abs(test_score[test_score != 0]), na.rm=TRUE)/100 
  test_score <- test_score |> rowwise() |> mutate(across(everything(), \(x) replace_na(x, runif(n=1, min=-r, max=r)))) |> ungroup()
  
  list(test_score = test_score)
}

# iterate over the variable columns and use the extra info from before to remap the levels in the test data. 
prepare_test_ca0 <- function(data, list_of_extras, id, residualised=NULL) {
  # first remap the variable levels to the appropriate ordinal level
  test_data <- data |> select(any_of(names(list_of_extras)))
  newdata_score <- map2(test_data, list_of_extras, impute_score_ca0)
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- if(!is.null(residualised)) {
    bind_cols(data |> select(all_of(id), all_of(residualised)), newdata_pred)
  } else {
    bind_cols(data |> select(all_of(id)), newdata_pred)
  }
  newdata_pred
}
