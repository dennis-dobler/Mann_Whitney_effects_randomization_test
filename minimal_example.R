# Copyright (c) 2026 Dennis Dobler, Paavo Sattler, Nils Hichert, 
#                    Jörg-Tobias Kuhn, Lubna Amro
# Licensed under the MIT License. See LICENSE file for details.

################################################################################
#------------------------------------------------------------------------------#
#----------------------------- Minimal examples -------------------------------#
#------------------------------------------------------------------------------#
################################################################################



################################################################################
#------------------------------------------------------------------------------#
#-------------------------- One-factor example --------------------------------#
#------------------------------------------------------------------------------#
################################################################################

## CODY data from the real-data example in the paper.
## The three intervention conditions are analyzed separately.
## Within each condition, equality of the 11 marginal Mann-Whitney effects
## is tested.


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load("cody.rda")

## The functions of the following R file are required for running the rest below.
source("Mann_Whitney_effects_functions.R")


## Convert the data to the format required by Mann_Whitney_effects_test().
prepare_data <- function(data) {
  X <- t(as.matrix(data))
  L <- !is.na(X)
  
  ## Missing values are represented through L.
  ## Their entries in X are therefore replaced by arbitrary finite values.
  X[!L] <- 0
  
  list(X = X, L = L)
}


conditions <- c("Control", "CODY", "CODY+NIP")

## Note: Setting a seed makes the output of the R code reproducible.
##  The R code in the present version is not parallelized.
##  Due to the nature of randomization tests, the resulting p-values of the
##  randomization tests depend on the chosen seed. With a larger number R
##  of randomization iterations, the variability of the p-values decrease.
set.seed(31415)

cody_results <- lapply(
  conditions,
  function(condition) {
    
    data_condition <- cody[
      cody$condition == condition,
      paste0("test", 1:11)
    ]
    
    data_condition <- prepare_data(data_condition)
    
    Mann_Whitney_effects_test(
      X = data_condition$X,
      L = data_condition$L,
      d_A = 11,
      d_B = 1,
      hypothesis_type = "equal",
      randomization_type = "all",
      test_type = "randomization",
      R = 2000
    )
  }
)

names(cody_results) <- conditions


## P-values
sapply(
  cody_results,
  function(x) x$p_value
)


################################################################################
#------------------------------------------------------------------------------#
#-------------------------- Two-factor example --------------------------------#
#------------------------------------------------------------------------------#
################################################################################

## Simulated 3 x 2 factorial setting corresponding to one of the simulation
## settings considered in the paper: Gumbel copula with Kendall's tau = 0.2,
## nonexchangeable components, and MCAR missingness generated using an
## exchangeable Gaussian copula with Kendall's tau = 0.2.

load("simulated_twofactor.rda")

twofactor_data <- prepare_data(simulated_twofactor)

## First approach: Randomize everything

## Note: Setting a seed makes the output of the R code reproducible.
##  The R code in the present version is not parallelized.
##  Due to the nature of randomization tests, the resulting p-values of the
##  randomization tests depend on the chosen seed. With a larger number R
##  of randomization iterations, the variability of the p-values decrease.
set.seed(31415)


## Equality of all six marginal Mann-Whitney effects.
result_equal <- Mann_Whitney_effects_test(
  X = twofactor_data$X,
  L = twofactor_data$L,
  d_A = 3,
  d_B = 2,
  hypothesis_type = "equal",
  randomization_type = "all",
  test_type = "randomization",
  R = 2000
)


## Main effect of factor A.
result_A <- Mann_Whitney_effects_test(
  X = twofactor_data$X,
  L = twofactor_data$L,
  d_A = 3,
  d_B = 2,
  hypothesis_type = "A",
  randomization_type = "all",
  test_type = "randomization",
  R = 2000
)


## Interaction effect A x B.
result_AB <- Mann_Whitney_effects_test(
  X = twofactor_data$X,
  L = twofactor_data$L,
  d_A = 3,
  d_B = 2,
  hypothesis_type = "AB",
  randomization_type = "all",
  test_type = "randomization",
  R = 2000
)


## P-values
c(
  equal = result_equal$p_value,
  A = result_A$p_value,
  AB = result_AB$p_value
)



## The same as above, but with "minimal" randomization

## Note: Setting a seed makes the output of the R code reproducible.
##  The R code in the present version is not parallelized.
##  Due to the nature of randomization tests, the resulting p-values of the
##  randomization tests depend on the chosen seed. With a larger number R
##  of randomization iterations, the variability of the p-values decrease.
set.seed(31415)


## Main effect of factor A.
result_A <- Mann_Whitney_effects_test(
  X = twofactor_data$X,
  L = twofactor_data$L,
  d_A = 3,
  d_B = 2,
  hypothesis_type = "A",
  randomization_type = "minimal",
  test_type = "randomization",
  R = 2000
)


## Interaction effect A x B.
result_AB <- Mann_Whitney_effects_test(
  X = twofactor_data$X,
  L = twofactor_data$L,
  d_A = 3,
  d_B = 2,
  hypothesis_type = "AB",
  randomization_type = "minimal",
  test_type = "randomization",
  R = 2000
)


## P-values
c(
  A = result_A$p_value,
  AB = result_AB$p_value
)
