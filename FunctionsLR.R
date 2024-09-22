# Function that implements multi-class logistic regression.
#############################################################
# Description of supplied parameters:
# X - n x p training data, 1st column should be 1s to account for intercept
# y - a vector of size n of class labels, from 0 to K-1
# Xt - ntest x p testing data, 1st column should be 1s to account for intercept
# yt - a vector of size ntest of test class labels, from 0 to K-1
# numIter - number of FIXED iterations of the algorithm, default value is 50
# eta - learning rate, default value is 0.1
# lambda - ridge parameter, default value is 1
# beta_init - (optional) initial starting values of beta for the algorithm, should be p x K matrix

function_beta <- function(i, Pk_mat, y){
  indicator_sum <- sum(log(Pk_mat[y == (i-1), i]))
}

## Return output
##########################################################################
# beta - p x K matrix of estimated beta values after numIter iterations
# error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
# error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
# objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
LRMultiClass <- function(X, y, Xt, yt, numIter = 50, eta = 0.1, lambda = .5, beta_init = NULL){
  ## Check the supplied parameters as described. You can assume that X, Xt are matrices; y, yt are vectors; and numIter, eta, lambda are scalars. You can assume that beta_init is either NULL (default) or a matrix.
  ###################################
  n1 <- dim(X)[1]
  n2 <- dim(Xt)[1]
  K <- dim(X)[2]
  
  objective <- rep(0, numIter+1)
  error_train <- rep(0, numIter+1)
  error_test <- rep(0, numIter+1)
  
  # Check that the first column of X and Xt are 1s, if not - display appropriate message and stop execution.
  
  # Check for compatibility of dimensions between X and Y
  
  # Check for compatibility of dimensions between Xt and Yt
  
  # Check for compatibility of dimensions between X and Xt
  
  # Check eta is positive
  
  # Check lambda is non-negative
  
  # Check whether beta_init is NULL. If NULL, initialize beta with p x K matrix of zeroes. If not NULL, check for compatibility of dimensions with what has been already supplied.
  if(is.null(beta_init)){
    beta <- matrix(0, nrow = dim(X)[2], ncol = K)
  }
  else{
    beta <- beta_init
  }
  
  ## Calculate corresponding pk, objective value f(beta_init), training error and testing error given the starting point beta_init
  ##########################################################################
  
  X_beta <- X %*% beta
  exp_X_beta <- exp(X_beta)
  
  Pk_mat <- (1/rowSums(exp_X_beta)) * exp_X_beta
  
  sum_log <- sapply(1:K, \(i){function_beta(i, Pk_mat, y)})
  objective[1] <- sum(beta * beta) * (lambda / 2) - sum(sum_log)
  
  y_guess <- sapply(1:n1, \(i){which.max(Pk_mat[i, ])-1})
  error_train[1] <- (1 - sum(as.numeric(y == y_guess))/n1) * 100
  
  Xt_beta <- Xt %*% beta
  exp_Xt_beta <- exp(Xt_beta)
  
  Pk_test <- (1/rowSums(exp_Xt_beta)) * exp_Xt_beta
  yt_guess <- sapply(1:n2, \(i){which.max(Pk_test[i, ])-1})
  error_test[1] <- (1 - sum(as.numeric(yt == yt_guess))/n2) * 100
  
  ## Newton's method cycle - implement the update EXACTLY numIter iterations
  ##########################################################################
  
  # Within one iteration: perform the update, calculate updated objective function and training/testing errors in %
  
  for(i in 1:numIter){
    
    for(k in 1:K){
      
      w <- Pk_mat[ , k] * (1-Pk_mat[ , k])
      hessian <- crossprod(X, w*X)
      diag(hessian) <- diag(hessian) + lambda
      hessian_inv <- chol2inv(chol(hessian))
      
      prob_vec <- Pk_mat[ , k] - as.numeric(y == (k-1))
      jacobian <- crossprod(X, prob_vec) + lambda * beta[ , k]
      
      beta[ , k] <- beta[ , k] - eta * hessian_inv %*% jacobian
      
    }
    
    X_beta <- X %*% beta
    exp_X_beta <- exp(X_beta)
    Pk_mat <- (1/rowSums(exp_X_beta)) * exp_X_beta
    
    y_guess <- sapply(1:n1, \(i){which.max(Pk_mat[i, ])-1})
    error_train[i+1] <- (1 - sum(as.numeric(y == y_guess))/n1) * 100
    
    sum_log <- sapply(1:K, \(i){function_beta(i, Pk_mat, y)})
    objective[i+1] <- sum(beta * beta) * (lambda / 2) - sum(sum_log)
    
    Xt_beta <- Xt %*% beta
    exp_Xt_beta <- exp(Xt_beta)
    Pk_test <- (1/rowSums(exp_Xt_beta)) * exp_Xt_beta
    
    yt_guess <- sapply(1:n2, \(i){which.max(Pk_test[i, ])-1})
    error_test[i+1] <- (1 - sum(as.numeric(yt == yt_guess))/n2) * 100
  }
  
  ## Return output
  ##########################################################################
  # beta - p x K matrix of estimated beta values after numIter iterations
  # error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
  # error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
  # objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
  return(list(beta = beta, error_train = error_train, error_test = error_test, objective =  objective))
  #return(beta)
}
