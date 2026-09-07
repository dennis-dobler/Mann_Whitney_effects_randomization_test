# Copyright (c) 2026 Dennis Dobler, Paavo Sattler, Nils Hichert, 
#                    Jörg-Tobias Kuhn, Lubna Amro
# Licensed under the MIT License. See LICENSE file for details.

################################################################################
#------------------------------------------------------------------------------#
#----------------- This File contains all R helper functions ------------------#
#------------------------------------------------------------------------------#
################################################################################


## Note:
## The functions perform basic input checks but do not validate all assumptions
## required by the underlying statistical methods. Users are responsible for
## ensuring that the supplied data, hypotheses, and design specifications satisfy
## the respective methodological requirements. Although the code has been tested,
## minor implementation errors may still occur, particularly in special cases or
## for unusual combinations of input arguments.



if(!require("Rcpp"))
  install.packages("Rcpp")
## source helper functions implemented in C++
sourceCpp("helpers.cpp")


## H0 - Constructs the hypothesis matrix and hypothesis vector
##
## Input:
##   type:
##     Hypothesis type: "equal", "A", "B", "AB", or "custom".
##
##   d_A, d_B:
##     Numbers of levels of factors A and B.
##     For a one-factor design, set d_B = 1.
##     For two-factor designs, the components are assumed to be ordered as
##     (A1,B1), ..., (A1,Bd_B), ..., (Ad_A,B1), ..., (Ad_A,Bd_B).
##
##   hypothesis_matrix:
##     User-defined hypothesis matrix for type = "custom".
##     It must have d_A * d_B columns.
##
##   hypothesis_vector:
##     User-defined hypothesis vector for type = "custom".
##     Its length must equal the number of rows of hypothesis_matrix.
##
## Output:
##   A named list containing:
##     hypothesis_matrix - numeric hypothesis matrix
##     hypothesis_vector - numeric hypothesis vector
H0 <- function(
    type = c("equal", "A", "B", "AB", "custom"),
    d_A,
    d_B = 1,
    hypothesis_matrix = NULL,
    hypothesis_vector = NULL) {
  
  allowed_types <- c(
    "equal",
    "A",
    "B",
    "AB",
    "custom"
  )
  
  if (length(type) != 1L) {
    stop("The type of hypothesis must be specified.")
  }
  
  if (!(type %in% allowed_types)) {
    stop(
      paste0(
        "Unknown hypothesis type '",
        type,
        "'."
      )
    )
  }
  
  if (
    length(d_A) != 1L ||
    !is.numeric(d_A) ||
    !is.finite(d_A) ||
    d_A < 1 ||
    d_A != floor(d_A)
  ) {
    stop("'d_A' must be a positive integer.")
  }
  
  if (
    length(d_B) != 1L ||
    !is.numeric(d_B) ||
    !is.finite(d_B) ||
    d_B < 1 ||
    d_B != floor(d_B)
  ) {
    stop("'d_B' must be a positive integer.")
  }
  
  d_A <- as.integer(d_A)
  d_B <- as.integer(d_B)
  
  d <- d_A * d_B
  
  
  ## Constructs the centering matrix of dimension n.
  P_mat <- function(n) {
    diag(n) -
      matrix(
        1 / n,
        nrow = n,
        ncol = n
      )
  }
  
  
  ###########################################################################
  # User-defined hypothesis
  ###########################################################################
  
  if (type == "custom") {
    
    if (is.null(hypothesis_matrix)) {
      stop(
        paste(
          "'hypothesis_matrix' must be supplied",
          "for type = 'custom'."
        )
      )
    }
    
    if (
      !is.matrix(hypothesis_matrix) ||
      !is.numeric(hypothesis_matrix)
    ) {
      stop(
        "'hypothesis_matrix' must be a numeric matrix."
      )
    }
    
    if (
      anyNA(hypothesis_matrix) ||
      any(!is.finite(hypothesis_matrix))
    ) {
      stop(
        "'hypothesis_matrix' must contain only finite values."
      )
    }
    
    if (nrow(hypothesis_matrix) < 1L) {
      stop(
        "'hypothesis_matrix' must contain at least one row."
      )
    }
    
    if (ncol(hypothesis_matrix) != d) {
      stop(
        paste0(
          "The number of columns of 'hypothesis_matrix' ",
          "must equal d_A * d_B = ",
          d,
          "."
        )
      )
    }
    
    if (is.null(hypothesis_vector)) {
      stop(
        paste(
          "'hypothesis_vector' must be supplied",
          "for type = 'custom'."
        )
      )
    }
    
    if (!is.numeric(hypothesis_vector)) {
      stop(
        "'hypothesis_vector' must be numeric."
      )
    }
    
    if (
      anyNA(hypothesis_vector) ||
      any(!is.finite(hypothesis_vector))
    ) {
      stop(
        "'hypothesis_vector' must contain only finite values."
      )
    }
    
    if (
      length(hypothesis_vector) !=
      nrow(hypothesis_matrix)
    ) {
      stop(
        paste(
          "The length of 'hypothesis_vector' must equal",
          "the number of rows of 'hypothesis_matrix'."
        )
      )
    }
    
    C <- hypothesis_matrix
    c0 <- as.numeric(hypothesis_vector)
  }
  
  
  ###########################################################################
  # Predefined hypotheses
  ###########################################################################
  
  if (type != "custom") {
    
    if (type == "equal") {
      C <- P_mat(d)
    }
    
    if (d_B == 1L) {
      
      if (type %in% c("B", "AB")) {
        stop(
          "Hypotheses for factor B or the interaction AB are not defined when d_B = 1."
        )
      }
      
      if (type == "A") {
        C <- P_mat(d_A)
      }
      
    } else {
      
      if (type == "A") {
        C <- kronecker(
          P_mat(d_A),
          t(rep(1, d_B))
        )
      }
      
      if (type == "B") {
        C <- kronecker(
          t(rep(1, d_A)),
          P_mat(d_B)
        )
      }
      
      if (type == "AB") {
        C <- kronecker(
          P_mat(d_A),
          P_mat(d_B)
        )
      }
    }
    
    ## All predefined hypotheses use a zero hypothesis vector.
    c0 <- numeric(nrow(C))
  }
  
  
  ###########################################################################
  # Return hypothesis matrix and hypothesis vector
  ###########################################################################
  
  list(
    hypothesis_matrix = C,
    hypothesis_vector = c0
  )
}

## myPermutation - Generates permutation indices for one- or two-factor designs
##
## Input:
##   n:
##     [integer] Number of observational units.
##
##   d_A, d_B:
##     [integer] Numbers of levels of factors A and B.
##     For a one-factor design, set d_B = 1.
##     For two-factor designs, the components are assumed to be ordered as
##     (A1,B1), ..., (A1,Bd_B), ..., (Ad_A,B1), ..., (Ad_A,Bd_B).
##
##   hypothesis_type:
##     [character] Hypothesis type: "equal", "A", "B", or "AB".
##
##   randomization_type:
##     [character] Randomization scheme. "all" permutes all components,
##     whereas "minimal" permutes only the components required for the
##     hypotheses "A", "B", or "AB".
##     The "minimal" scheme is not defined for hypothesis_type = "equal".
##
##   n_perm:
##     [integer] Number of permutations.
##
## Output:
##   [integer array] An array of dimension d x n x n_perm, where
##   d = d_A * d_B. Each array slice [, , i] contains the linear indices
##   corresponding to one permutation of the data matrix.
myPermutation <- function(
  n, d_A, d_B = 1, 
  hypothesis_type = c("equal", "A", "B", "AB"),
  randomization_type = "all",
  n_perm = 2000){
  
  if (!(hypothesis_type[1] %in% c("equal", "A", "B", "AB")))
    stop("Unknown hypothesis_type")
  
  if (!(randomization_type[1] %in% c("all", "minimal")))
    stop("Unknown randomization_type")
  
  d <- d_A * d_B
  
  if (randomization_type[1] == "minimal" &&
      hypothesis_type[1] == "equal") {
    stop("randomization_type = 'minimal' is not defined for hypothesis_type = 'equal'.")
  }
  
  ## Case 1: randomization_type is 'all'
  if(randomization_type == "all"){
    ind <- array(replicate(n * n_perm, sample.int(d)), dim = c(d, n, n_perm))
  }
  if(randomization_type == "minimal"){
    ## Case 2.1: hypothesis is "A"
    if(hypothesis_type[1] == "A"){
      ind <- array(
        replicate(
          n*n_perm, 
          rep(sample(d_A), each = d_B) * d_B - ((d_B-1):0)), 
        dim = c(d, n, n_perm))
    }
    
    ## Case 2.2: hypothesis is "B"
    if(hypothesis_type[1] == "B"){
      ind <- array(
        replicate(
          n*n_perm, 
          rep(sample(d_B), times = d_A) + rep((0:(d_A-1))*d_B, each = d_B)), 
        dim = c(d, n, n_perm))
    }
    
    ## Case 2.3: hypothesis is "AB"
    if(hypothesis_type[1] == "AB"){
      ind1 <- array(
        replicate(
          n*n_perm, 
          rep(sample(d_A), each = d_B) * d_B - ((d_B-1):0)), 
        dim = c(d, n, n_perm))
      ind2 <- array(
        replicate(
          n*n_perm, 
          rep(sample(d_B), times = d_A) + rep((0:(d_A-1))*d_B, each = d_B)), 
        dim = c(d, n, n_perm))
      ind <- array(ind1[ind2], dim = c(d, n, n_perm))
    }
  }
  
  
  
  sweep(ind, 2, d*(1:n - 1), "+")
}


## Mann_Whitney_effects_test - Tests hypotheses on marginal Mann-Whitney effects
##    for one- or two-factor designs
##
## Input:
##   X:
##     [numeric matrix] Data matrix of dimension d x n, where rows correspond
##     to factor-level combinations and columns to observational units.
##     For two-factor designs, the rows are assumed to be ordered as
##     (A1,B1), ..., (A1,Bd_B), ..., (Ad_A,B1), ..., (Ad_A,Bd_B).
##
##   L:
##     [logical matrix] Missingness indicator matrix of dimension d x n.
##     TRUE (or 1) indicates an observed value and FALSE (or 0) a missing value.
##
##   d_A, d_B:
##     [integer] Numbers of levels of factors A and B.
##     For a one-factor design, set d_B = 1.
##
##   hypothesis_type:
##     [character] Hypothesis type: "equal", "A", "B", "AB", or "custom".
##
##   hypothesis_matrix:
##     [numeric matrix] User-defined hypothesis matrix for
##     hypothesis_type = "custom". It must have d_A * d_B columns.
##
##   hypothesis_vector:
##     [numeric vector] User-defined hypothesis vector for
##     hypothesis_type = "custom". Its length must equal the number of rows
##     of hypothesis_matrix.
##
##   randomization_type:
##     [character] Randomization scheme used for test_type = "randomization".
##     "all" permutes all components, whereas "minimal" permutes only the
##     components required for the selected predefined hypothesis.
##
##   test_type:
##     [character] Method used to determine the reference distribution:
##     "randomization", "bootstrap", or "asymptotic".
##
##   R:
##     [integer] Number of randomization or bootstrap samples.
##     Ignored for test_type = "asymptotic".
##
##   alpha:
##     [numeric] Significance level. Must be in [0, 1].
##
##   conf_int:
##     [logical] If TRUE, computes a confidence interval for the linear
##     contrast defined by the first row of the hypothesis matrix.
##     Not available for test_type = "asymptotic".
##
## Output:
##   A named list containing:
##     p_value            - p-value of the test
##     decision           - logical test decision at significance level alpha
##     only_zeros_counter - number of resampling iterations in which at least
##                          one row of the missingness matrix contains only zeros
##     phat               - estimated value of C %*% p, where p denotes the
##                          vector of marginal Mann-Whitney effects
##     CI_estimate        - estimated contrast for the first row of C,
##                          if a confidence interval is requested
##     CI_lower, CI_upper - lower and upper confidence limits, if requested
##     test_statistic     - value of the test statistic
Mann_Whitney_effects_test <- function(
    X,
    L,
    d_A,
    d_B = 1,
    hypothesis_type = c("equal", "A", "B", "AB", "custom"),
    hypothesis_matrix = NULL,
    hypothesis_vector = NULL,
    randomization_type = c("all", "minimal"),
    test_type = c("randomization", "bootstrap", "asymptotic"),
    R = 2000,
    alpha = 0.05,
    conf_int = FALSE) {
  
  if (!(hypothesis_type[1] %in% c("equal", "A", "B", "AB", "custom"))) {
    stop('hypothesis_type muss "equal", "A", "B", "AB" oder "custom" sein')
  }
  
  if (!(test_type[1] %in% c("randomization", "bootstrap", "asymptotic"))) {
    stop('test_type muss "randomization", "bootstrap" oder "asymptotic" sein')
  }
  
  if (conf_int && test_type[1] == "asymptotic") {
    stop("Confidence intervals are not implemented for test_type = 'asymptotic'.")
  }
  
  if (hypothesis_type[1] == "custom") {
    
    if (test_type[1] == "asymptotic") {
      stop("The asymptotic test is not available for custom hypotheses.")
    }
    
    if (is.null(hypothesis_matrix)) {
      stop("'hypothesis_matrix' must be supplied for custom hypotheses.")
    }
    
    if (is.null(hypothesis_vector)) {
      stop("'hypothesis_vector' must be supplied for custom hypotheses.")
    }
    
    if (ncol(hypothesis_matrix) != nrow(X)) {
      stop(
        "The number of columns of 'hypothesis_matrix' must equal nrow(X)."
      )
    }
    
    if (length(hypothesis_vector) != nrow(hypothesis_matrix)) {
      stop(
        "The length of 'hypothesis_vector' must equal nrow(hypothesis_matrix)."
      )
    }
  }
  
  if(alpha < 0 | alpha > 1) stop("alpha must be in [0, 1].")
  
  if (nrow(X) != d_A * d_B) {
    stop("nrow(X) must equal d_A * d_B.")
  }
  
  ## Construct the hypothesis matrix C.
  hypothesis <- H0(
    type = hypothesis_type[1],
    d_A = d_A,
    d_B = d_B,
    hypothesis_matrix = hypothesis_matrix,
    hypothesis_vector = hypothesis_vector
  )
  
  C <- hypothesis$hypothesis_matrix
  c0 <- hypothesis$hypothesis_vector
  
  ## Extract the dimensions from the data matrix.
  n <- ncol(X)
  d <- nrow(X)
  
  
  #### Estimate p and V from the original data.
  p <- est_p_cpp(X, L)
  V <- est_V_cpp(X, L)
  phat <- drop(C %*% p)
  difference_orig <- drop(C %*% p - c0)
  
  T_orig <- drop(
    t(difference_orig) %*%
      MASS::ginv(C %*% V %*% t(C)) %*%
      difference_orig
  ) * n
  ## If the test statistic is negative, return NA as the test decision.
  if(T_orig < 0) return(list(decision =  NA, only_zeros_counter = 0)) ## This biases the only_zeros_counter!
  if (test_type[1] == "asymptotic") {
    
    df <- if (hypothesis_type[1] == "A") {
      d_A - 1
    } else if (hypothesis_type[1] == "B") {
      d_B - 1
    } else if (hypothesis_type[1] == "AB") {
      (d_A - 1) * (d_B - 1)
    } else {
      d - 1
    }
    
    p_value <- 1 - pchisq(T_orig, df = df)
    
    return(list(
      p_value = p_value,
      decision = p_value <= alpha,
      only_zeros_counter = 0,
      phat = phat,
      test_statistic = T_orig
    ))
  }
  
  ## Generate all randomization or bootstrap samples at once.
  if(test_type[1] == "randomization"){
    ## Randomize the indices R times.
    if (hypothesis_type[1] == "custom") {
      permutation_hypothesis <- "equal"
      randomization_type_used <- "all"
    } else {
      permutation_hypothesis <- hypothesis_type[1]
      randomization_type_used <- randomization_type[1]
    }
    
    ind <- myPermutation(
      n = n,
      d_A = d_A,
      d_B = d_B,
      hypothesis_type = permutation_hypothesis,
      randomization_type = randomization_type_used,
      n_perm = R
    )
  }else if(test_type[1] == "bootstrap"){
    ## Draw R bootstrap samples.
    ind <- array(replicate(n*R, sample(1:d, replace=TRUE))
                 + rep(0:(n-1), each = d)*d, dim =c( c(d,n,R)))
  }
  
  if (conf_int) {
    ## Use the first row of C for the confidence interval.
    C_CI <- C[1, , drop = FALSE]
    T_pi_CI <- numeric(R)
  }
  
  T_pi <- numeric(R)
  
  ## Count how often a row of L consists only of zeros.
  only_zeros_counter <- 0
  
  for (r in 1:R) {
    ## Compute the test statistic for each resampling iteration.
    i <- ind[, , r]
    L_pi <- matrix(L[i], ncol = n)
    
    if (any(rowSums(L_pi) == 0)) {
      T_pi[r] <- Inf
      
      if (conf_int) {
        T_pi_CI[r] <- Inf
      }
      
      only_zeros_counter <- only_zeros_counter + 1
      
    } else {
      X_pi <- matrix(X[i], ncol = n)
      
      p_pi <- est_p_cpp(X_pi, L_pi)
      V_pi <- est_V_cpp(X_pi, L_pi)
      
      ## Center the randomized Mann-Whitney effects at their randomization
      ## reference value 1/2.
      centered_p_pi <- p_pi - 0.5
      
      difference_pi <- drop(
        C %*% centered_p_pi
      )
      
      T_pi[r] <- drop(
        t(difference_pi) %*%
          MASS::ginv(
            C %*% V_pi %*% t(C)
          ) %*%
          difference_pi
      ) * n
      
      if (conf_int) {
        difference_pi_CI <- drop(
          C_CI %*% centered_p_pi
        )
        
        Sigma_pi_CI <- drop(
          C_CI %*% V_pi %*% t(C_CI)
        )
        
        T_pi_CI[r] <- n *
          difference_pi_CI^2 /
          Sigma_pi_CI
      }
    }
  }
  
  
  T_pi <- T_pi[!is.na(T_pi)]
  
  if (length(T_pi) < 1) {
    return(
      list(
        decision = NA,
        only_zeros_counter = only_zeros_counter
      )
    )
  }
  
  p_value <- mean(
    T_pi >= T_orig
  )
  
  critical_value <- quantile(
    T_pi,
    1 - alpha,
    names = FALSE
  )
  
  decision <- T_orig > critical_value
  
  
  if (conf_int) {
    T_pi_CI <- T_pi_CI[!is.na(T_pi_CI)]
    
    if (length(T_pi_CI) < 1) {
      return(
        list(
          p_value = p_value,
          decision = decision,
          only_zeros_counter = only_zeros_counter,
          phat = phat,
          CI_estimate = NA,
          CI_lower = NA,
          CI_upper = NA,
          test_statistic = T_orig
        )
      )
    }
    
    q_star_CI <- quantile(
      T_pi_CI,
      1 - alpha,
      names = FALSE
    )
    
    phat_CI <- drop(
      C_CI %*% p
    )
    
    Sigma_CI <- drop(
      C_CI %*% V %*% t(C_CI)
    )
    
    CI_lower <- phat_CI -
      sqrt(q_star_CI) *
      sqrt(Sigma_CI) /
      sqrt(n)
    
    CI_upper <- phat_CI +
      sqrt(q_star_CI) *
      sqrt(Sigma_CI) /
      sqrt(n)
    
    return(
      list(
        p_value = p_value,
        decision = decision,
        only_zeros_counter = only_zeros_counter,
        phat = phat,
        CI_estimate = phat_CI,
        CI_lower = CI_lower,
        CI_upper = CI_upper,
        test_statistic = T_orig
      )
    )
  }
  
  
  list(
    p_value = p_value,
    decision = decision,
    only_zeros_counter = only_zeros_counter,
    phat = phat,
    test_statistic = T_orig
  )
}
