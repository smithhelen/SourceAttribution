# Recipe for CAP method

source('methods/helpers.R')

# constructor function for our recipe step
step_cap_new <- 
  function(terms, role, trained, distances, ck, cm, cmp, c, objects, options, skip, id) {
    step(
      subclass = "cap", 
      terms = terms,
      role = role,
      trained = trained,
      distances = distances,
      ck = ck,
      cm = cm,
      cmp = cmp,
      c = c,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  }

# user facing function for our recipe
step_cap <- function(
    recipe, 
    ..., 
    role = "predictor", 
    trained = FALSE,
    distances,
    ck = NULL,
    cm = NULL,
    cmp = 100,
    c = NULL,
    objects = NULL, # info from our training data
    skip = FALSE,
    options = list(),
    id = rand_id("cap")
) {
  
  if (TRUE) { #!recipes::is_tune(axes)) { # TODO: Do we need to check for tuning here???
    cm <- as.integer(cm)
  }
  add_step(
    recipe, 
    step_cap_new(
      terms = enquos(...),
      trained = trained,
      role = role,
      distances = distances,
      ck = ck,
      cm = cm,
      cmp = cmp,
      c = c,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  )
}

# Prep step

# CAP scoring function
encode_cap <- function(var, outcome, distance, ck, cm, cmp, c) {
  var <- droplevels(var)
  d.train <- distance[levels(var),levels(var),drop=FALSE]
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  
  # use only non-zero eigenvalues
  nlambdas <- sum(eigen_B$values > epsilon)
  if (nlambdas == 0) {
    # No non-zero eigenvectors
    return(NULL)
  }

  ct <- table(level=var, outcome=outcome)
  H <- hat_fn(ct, k=ck) # restrict ct to ck axes, if ck is null then ck=ncol(ct)-1
  
  # Restrict to a maximum of eigenvectors set by "cm" or "cmp (propG)" (default is cmp=100% variation)
  lambdas_B <- filter_eigenvalues(eigen_B$values[seq_len(nlambdas)], m=cm, mp=cmp) # restrict by axes or cmp
  # Scale eigenvectors
  Qo <- eigen_B$vectors[, seq_along(lambdas_B), drop=FALSE]  # note that this is different to the Q score in PCO method which is scaled by the sqrt(abs(lambdas_B))
  QHQ <- t(Qo) %*% H %*% Qo   # Combine Qo and H to get C_score
  eigen_QHQ <- eigen_decomp(QHQ, symmetric=TRUE)
  
  # select number of axes to retain - default is number of pco axes (NOTE CHANGE, used to be number of k axes)
  if(is.null(c)){c <- length(lambda_B)}
  #if(is.null(c)){c <- ncol(ct)-1}
  lambda_QHQ <- filter_eigenvalues(eigen_QHQ$values, m=c)
  U <- eigen_QHQ$vectors[,seq_along(lambda_QHQ),drop=FALSE]
  C_score <- Qo %*% U
  
  objects <- list(levels = levels(var),
                  d=distance,
                  ct=ct,
                  H=H,
                  lambdas_B=lambdas_B,
                  eigen_B=eigen_B,
                  nlambdas=nlambdas,
                  U=U,
                  Qo=Qo,
                  QHQ=QHQ,
                  eigenQHQ=eigen_QHQ,
                  lambda_QHQ=lambda_QHQ,
                  C_score=C_score,
                  diag.B.train = diag(B.train),
                  propG = cumsum(eigen_B$values)/sum(eigen_B$values)*100)
  objects
  #print(objects)
}

prep.step_cap <- function(x, training, info = NULL, ...) {
  # x is the object from the step_pco function, 
  # training is the training set data (tibble),
  # and info is a tibble that has information on the current set of data eg variable name, type, and role
  
  # grab the columns we're going to prep
  col_names <- recipes_eval_select(x$terms, training, info)
  
  # grab the outcome column
  outcome_name <- info |> filter(role == "outcome") |> pull(variable)
  if (length(outcome_name) != 1) {
    rlang::abort("One variable with role 'outcome' is required")
  }
  
  # check the column types are what we want: we want factor variables
  # (or categorical variables?)
  check_type(training[, col_names], types = c("character", "factor"))
  check_type(training[, outcome_name], types = c("character", "factor"))
  
  # check we have the information we need for distances
  if (!is.list(x$distances)) {
    rlang::abort("distances should be a list object of matrices")
  }
  if (!all(col_names %in% names(x$distances))) {
    rlang::abort("distances should have entries for all columns, named after the column")
  }
  if (!all(map_lgl(x$distances, ~ any(class(.) %in% c('dist', 'matrix'))))) {
    rlang::abort("distances should be of class 'dist' or 'matrix'")
  }
  # convert distances to matrices
  x$distances <- map(x$distances, as.matrix)
  
  # OK, now do the actual CAP step on each column
  # This computes the m etc using the distance matrices.
  # The actual data transformation
  # of current variables is done in 'prep' or 'bake' not here
  objects <- purrr::map2(
    training[, col_names],
    x$distances[col_names],
    \(var, dist) encode_cap(var=var, 
                            outcome = training |> pull(outcome_name), 
                            distance = dist, 
                            ck=x$ck, 
                            cm = x$cm, 
                            cmp = x$cmp, 
                            c=x$c)
  )
  
  ## Use the constructor function to return the updated object. 
  ## Note that `trained` is now set to TRUE
  
  step_cap_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    distances = x$distances,
    ck = x$ck,
    cm = x$cm,
    cmp = x$cmp,
    c = x$c,
    objects = objects,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}

# Bake step: take our scores and apply them as needed to the columns
apply_cap_to_column <- function(var, encoding) {
  
  # Gower's transformation of new observations into PCO space
  new_level_to_cap <- function(new.var_level, d, diag.B.train, Qo, lambdas_B, U) {
    d.new <- d[new.var_level, names(diag.B.train)]
    d.gower <- diag.B.train - (d.new ^ 2)
    newQo <- d.gower %*% Qo / (2 * lambdas_B)
    #for plotting, the canonical variable scores are standardized by the square root of their corresponding eigenvalue (lambda_QHQ). 
    #this is not necessary for the CAP method per se
    #newCscore <- newQo %*% U * sqrt(abs(lambda_QHQ))
    newCscore <- newQo %*% U
    #colnames(newCscore) = colnames(Qo)
    new_var_level_score <- data.frame(level = new.var_level, newCscore)
    new_var_level_score
  }
  
  var <- droplevels(var) # ignore levels we don't have in these data
  
  # Now we figure out which levels are new, and which are the usual
  new_levels <- setdiff(levels(var), encoding$levels)
  new_scores <- map(new_levels, ~new_level_to_cap(new.var_level = {.},
                                                  d = encoding$d,
                                                  diag.B.train = encoding$diag.B.train,
                                                  Qo = encoding$Qo,
                                                  lambdas_B = encoding$lambdas_B,
                                                  U = encoding$U)) |>
    list_rbind()
  
  var_level_score <- bind_rows(data.frame(encoding$C_score) %>% rownames_to_column("level"), new_scores)
  new_cols <- data.frame(level = var) |> left_join(var_level_score, by="level") |> select(-level)
  new_cols
}

bake.step_cap <- function(object, new_data, ...) {
  # object is the updated step function that has been prepped, 
  # new_data is a tibble of data to be processed.
  
  col_names <- names(object$objects)
  check_new_data(col_names, object, new_data)

  # generate some new names
  new_names <- imap(object$objects, \(x, nm) { paste(nm, "cap", seq_along(colnames(x$C_score)), sep="_") })
  new_tbl <- tibble::new_tibble(x = list(), nrow=nrow(new_data))

  # iterate over and generate our new columns
  for (col_name in col_names) {
    i_col <- new_data[[col_name]]
    i_obj <- object$objects[[col_name]]
    if (!is.null(i_obj)) { # only if we need to include this column...
      new_cols <- apply_cap_to_column(var = i_col, encoding = i_obj)
      new_col_names <- new_names[[col_name]]
      colnames(new_cols) <- new_col_names
      new_tbl[new_col_names] <- new_cols
    }
  }
  
  # check the new names and produce our final dataset
  new_tbl <- check_name(new_tbl, new_data, object, names(new_tbl))
  new_data <- bind_cols(new_data, new_tbl)
  new_data <- dplyr::select(new_data, -dplyr::all_of(col_names))
  new_data
}

print.step_cap <-
  function(x, width = max(20, options()$width - 35), ...) {
    title <- "CAP transformation on "
    
    print_step(
      # Names after prep:
      tr_obj = names(x$objects),
      # Names before prep (could be selectors)
      untr_obj = x$terms,
      # Has it been prepped? 
      trained = x$trained,
      # What does this step do?
      title = title,
      # An estimate of how many characters to print on a line: 
      width = width
    )
    invisible(x)
  }

format_cap <- function(x) {
  if (is.null(x)) {
    # consistency
    tibble(
      value = character(),
      axis = integer(),
      cap = numeric(),
      prop_var = numeric()
    )
  } else {
    prop_var <- data.frame(axis = seq_along(x$propG), prop_var = x$propG)
    x$Q |> as.data.frame() |>
      tibble::rownames_to_column("value") |>
      pivot_longer(-value, names_to="axis", values_to="cap", names_prefix="V", names_transform = as.integer) |>
      left_join(prop_var, by='axis')
  }
}

tidy.step_cap <- function(x, ...) {
  if (is_trained(x)) {
    if (length(x$objects) == 0) {
      # We need to create consistent output when no variables were selected
      res <- tibble(
        term = character(),
        value = character(),
        axis = integer(),
        prop_var = numeric(),
        cap = numeric()
      )
    } else {
      res <- map(x$objects, format_cap) |> list_rbind(names_to = "term")
    }
  } else {
    term_names <- sel2char(x$terms)
    res <-
      tibble(
        term = term_names,
        value = rlang::na_chr,
        axis = rlang::na_int,
        prop_var = rlang::na_dbl,
        cap = rlang::na_dbl
      )
  }
  # Always return the step id: 
  res$id <- x$id
  res
}
