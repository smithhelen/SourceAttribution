# variable importance using independent hold-out

# load libraries
#library(tidyverse)
#source("methods/libs_fns.R")


permute <- function(var, data) {
  data |> mutate({{var}} := sample(.data[[{{var}}]]))
}

mc <- function(mod, dat, class){
  preds <- data.frame(prediction = predict(mod, data=dat, predict.all = FALSE)$predictions, truth = dat |> pull({{class}})) |> table()
  misclass <- 1-sum(diag(preds))/sum(preds)
  misclass
}

new_vi <- function(dat.one, dat.two, class, id, ntrees=500){
  
  # set seed
  set.seed(1234)
  
  ## 1st fold
  # run default ranger with ordered method
  r_mod_1 <- ranger(dependent.variable.name = class, data=(dat.one$train |> select(!any_of(id))), 
                    num.trees=ntrees, respect.unordered.factors = TRUE)

  # calculate misclassification rates
  misclass_1 <- mc(r_mod_1, dat.one$test, class=class)
  
  # permute variable in dat.one test
  vars_1 <- dat.one$test |> select(!all_of(c(id, class))) |> colnames()
  list.perm.dat_1 <- map(as.list(vars_1), permute, data=dat.one$test)
  misclass.perm_1 <- list.perm.dat_1 |> map(mc, mod=r_mod_1, class=class)
  names(misclass.perm_1) <- vars_1
  
  # calculate variable importance
  vi_1 <- misclass.perm_1 |> map(\(x) {x - misclass_1}) |> unlist()
  vi_1 <- vi_1 |> as.data.frame() |> rownames_to_column("var") 
  
  ## 2nd fold
  # run default ranger with ordered method
  r_mod_2 <- ranger(dependent.variable.name = class, data=dat.two$train |> select(!any_of(id)), 
                    num.trees=ntrees, respect.unordered.factors = TRUE)
  
  # calculate misclassification rates
  misclass_2 <- mc(r_mod_2, dat.two$test, class=class)
  
  # permute variable in dat.two test
  vars_2 <- dat.two$test |> select(!all_of(c(id, class))) |> colnames()
  list.perm.dat_2 <- map(as.list(vars_2), permute, data=dat.two$test)
  misclass.perm_2 <- list.perm.dat_2 |> map(mc, mod=r_mod_2, class=class)
  names(misclass.perm_2) <- vars_2
  
  # calculate variable importance
  vi_2 <- misclass.perm_2 |> map(\(x) {x - misclass_2}) |> unlist()
  vi_2 <- vi_2 |> as.data.frame() |> rownames_to_column("var") 
  
  # output results
  new_vi=left_join(vi_1, vi_2, by="var") |> mutate(vi = round(rowMeans(pick(where(is.numeric))), 6), .keep = "unused")
  new_vi
  #out <- list(new_vi = new_vi, vi1 = vi_1, vi2= vi_2, mc1 = misclass_1, mc2 = misclass_2)
  #out
}


