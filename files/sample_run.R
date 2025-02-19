### Source attribution using whole genome sequencing data and Hamming distances ###

## Load libraries
library(tidyverse)
library(ranger)
library(tidymodels)


## Load functions
source("methods/hamming.R")                 # calculate Hamming distance between pairs of alleles
source("methods/residualise.R")             # calculate Residualised distance between pairs of alleles
source("methods/ca_unbiased.R")             # CA method with new levels scored as zero
source("methods/pco.R")                     # PCO method
source("methods/cap.R")                     # CAP method
source("methods/vi_independent_holdout.R")  # Independent Holdout Variable Importance method


## Load sample data 
# data is genetic data with alleles (levels) from 10 genes (variables) across 3 sources (classes), n=10 for each source
load("data/SourceAttribution.RData") # sample data
load("data/SeqDat.RData") # sequence information for sample data

# training data is the source data
Dat.train <- SourceAttribution |> filter(class != "Human") |> select(where(~n_distinct(.) > 1))  |> droplevels() |> mutate(across(everything(), factor)) 

# test data is the human data to be attributed to a source
Dat.test <- SourceAttribution |> filter(class == "Human") |> droplevels() |> select(colnames(Dat.train)) |> mutate(across(everything(), factor)) 


## Prepare data
# calculate hamming distances between levels of predictor variable
genes <- colnames(Dat.train |> select(starts_with("Var")))
list_of_distance_matrices <- map(genes, ~ dfun(gene=., dat=SourceAttribution, seqdat=SeqDat))
names(list_of_distance_matrices) <- genes

# residualise option
d <- d_resid(SourceAttribution, d_Hamming = list_of_distance_matrices, residualised = "CC")


# set seed
set.seed(123)

# encode data using method of choice:
#' reminder of terminology: k = number ca axes (default is num.classes-1), 
#'                          m = number pco axes (default is num.classes-1), 
#'                          mp = propG (default is 100%), 
#'                          c = number cap axes (default is num.classes-1)
#'                          ck = number ca axes within cap (default is num.classes-1)
#'                          cm = number pco axes within cap (default is num.classes-1)
#'                          cmp = propG within cap (default is 100%)

# 1. ca unbiased method
train <- prepare_training_ca0(Dat.train, starts_with("Var"), "class")
test <- prepare_test_ca0(Dat.test, train$extra, "id")

# 2. pco method
train <- prepare_training_pco(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, mp=95)
test <- prepare_test_pco(Dat.test, train$extra, "id")

# 3. cap method
train <- prepare_training_cap(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, cmp=95)
test <- prepare_test_cap(Dat.test, train$extra, "id")


## generate random forest models
rf_mod <- ranger(class ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)


## make predictions for Human data
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() |> mutate(p = round(Freq/sum(Freq)*100, round(2))) # counts of predictions for each class

# Prediction intervals
# 1. probability forest
set.seed(234)
rf_mod2 <- ranger(class ~ ., data = train$training, num.trees=500, respect.unordered.factors = TRUE, probability = TRUE)

# probability for each class for each tree
probs <- predict(rf_mod2, data=test, predict.all = TRUE)$predictions  # probs[individual, class, tree]

# mean probability for each class for each tree (500 means)
class_probs <- map_df(1:500, \(tree) {map(1:3, \(class) probs[,class,tree] |> mean()) |> unlist() |> t() |> as.data.frame()}) |> 
  magrittr::set_colnames(c("Cattle", "Chicken", "Sheep")) |> as_tibble()

# quantiles
prob_quantiles <- class_probs |> map_df(quantile, prob=c(0.05, 0.95, 0.5, 0.025, 0.975), .id="Class")
prob_quantiles

# convert back to counts
prob_quantiles |> mutate(across(where(is.numeric), \(x) x * nrow(Dat.test)))

# plot
plot_dat <- class_probs |> rownames_to_column("tree") |> 
  pivot_longer(cols = Cattle:Sheep,names_to = "Class", values_to = "Prob")

plot_dat |> ggplot(aes(x=tree,y=Prob,colour = Class)) + geom_point() + theme_bw() +
  labs(title = "Probability of each source for each tree", x="Tree (x500)",y="Probability") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(data=prob_quantiles, aes(yintercept = `5%`, colour=Class), linetype=2) +
  geom_hline(data=prob_quantiles, aes(yintercept = `95%`, colour=Class), linetype=2) +
  geom_hline(data=prob_quantiles, aes(yintercept = `2.5%`, colour=Class), linetype=3) +
  geom_hline(data=prob_quantiles, aes(yintercept = `97.5%`, colour=Class), linetype=3) +
  scale_colour_brewer(palette = "Dark2")

plot_dat |> ggplot(aes(x=Class,y=Prob, group=Class, fill=Class)) + geom_boxplot() + 
  scale_fill_manual(values=c("#404f89","#44781E","#B8321A"), guide = "none") +
  scale_y_continuous("Percentage of human cases",expand=c(0.03,0.03),limits=c(0,1),labels=scales::percent) +
  theme_bw(base_size = 12) 


## Variable Importance ###############
# Split data into training and test sets using folds
set.seed(345)

flds <- vfold_cv(SourceAttribution, v = 2)$splits
train.one.dat <- flds[[1]] |> training()
test.one.dat <- flds[[1]] |> testing()
train.two.dat <- flds[[2]] |> training()
test.two.dat <- flds[[2]] |> testing()

train.one <- prepare_training_cap(train.one.dat, vars=starts_with("Var"), class="class", d=list_of_distance_matrices, cmp=95)
train.two <- prepare_training_cap(train.two.dat, vars=starts_with("Var"), class="class", d=list_of_distance_matrices, cmp=95)
test.one <- prepare_test_cap(test.one.dat, train.one$extra, "id")
test.two <- prepare_test_cap(test.two.dat, train.two$extra, "id")

dat.one <- list(train = train.one$training |> bind_cols(train.one.dat |> select(id)), 
                test = test.one |> left_join(test.one.dat |> select(c(id, class))))
dat.two <- list(train = train.two$training |> bind_cols(train.two.dat |> select(id)), 
                test = test.two |> left_join(test.two.dat |> select(c(id, class))))

vi_sacnz <- new_vi(dat.one, dat.two, class="class", id="id", ntrees=500)
# note variables may be in this list more than once as each variable is represented by multiple principal components/coordinates
vi_sacnz |> arrange(desc(vi)) |> slice_head(n=5) 

