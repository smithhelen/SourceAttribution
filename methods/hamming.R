#### Calculate Hamming distances between levels of predictor variable ####

# Load libraries
library(seqinr) # for Hamming distances
library(stringdist) # for Hamming distances

# function to create distance (Hamming) matrices between alleles for each gene ie seq differences on alleles not individuals
dfun <- function(gene, dat, seqdat){  #gene is column name
  alleles <- dat %>% pull({gene}) %>% unique
  x <- seqdat %>% select({gene}) %>% distinct %>% pull({gene})
  d <- stringdistmatrix(x, x, method="hamming")
  rownames(d) <- alleles
  colnames(d) <- alleles
  return(d)
}

