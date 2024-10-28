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
  p <- dim(X)[2]
  K <- length(unique(y))
  
  X_trans <- t(X)
  
  objective <- rep(0, numIter+1)
  error_train <- rep(0, numIter+1)
  error_test <- rep(0, numIter+1)
  
  # Check that the first column of X and Xt are 1s, if not - display appropriate message and stop execution.
  
  if(sum(as.numeric(X[ ,1] == 1)) == n1 && sum(as.numeric(Xt[ ,1] == 1)) == n2){}
  else{
    stop("Error: First column of X and Xt must be 1's.")
  }
  
  # Check for compatibility of dimensions between X and Y
  
  if(length(y) == n1){}
  else{
    stop("Error: incompatible dimensions with X and Y.")
  }
  
  # Check for compatibility of dimensions between Xt and Yt
  
  if(length(yt) == n2){}
  else{
    stop("Error: incompatible dimensions with Xt and Yt.")
  }
  
  # Check for compatibility of dimensions between X and Xt
  
  if(dim(X)[2] == dim(Xt)[2]){}
  else{
    stop("Error: incompatible dimensions with X and Xt")
  }
  
  # Check eta is positive
  
  if(eta <= 0){
    stop("Error: eta must be positive.")
  }
  
  # Check lambda is non-negative
  
  if(lambda < 0){
    stop("Error: lambda must be non-negative")
  }
  
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
  
  Pk_mat <- exp_X_beta/rowSums(exp_X_beta)
  
  sum_log <- sapply(1:K, \(i){function_beta(i, Pk_mat, y)})
  objective[1] <- sum(beta * beta) * (lambda / 2) - sum(sum_log)
  
  y_guess <- max.col(Pk_mat)-1
  error_train[1] <- (1 - sum(as.numeric(y == y_guess))/n1) * 100
  
  Xt_beta <- Xt %*% beta
  exp_Xt_beta <- exp(Xt_beta)
  
  Pk_test <- exp_Xt_beta/rowSums(exp_Xt_beta)
  yt_guess <- max.col(Pk_test)-1
  error_test[1] <- (1 - sum(as.numeric(yt == yt_guess))/n2) * 100
  
  ## Newton's method cycle - implement the update EXACTLY numIter iterations
  ##########################################################################
  
  # Within one iteration: perform the update, calculate updated objective function and training/testing errors in %
  
  for(i in 1:numIter){
    
    for(k in 1:K){
      
      w <- Pk_mat[ , k] * (1-Pk_mat[ , k])
      hessian <- X_trans %*% (w * X)
      diag(hessian) <- diag(hessian) + lambda
      #hessian_inv <- chol2inv(chol(hessian))
      hessian_inv <- solve(hessian)
      
      prob_vec <- Pk_mat[ , k] - as.numeric(y == (k-1))
      jacobian <- X_trans %*% prob_vec + lambda * beta[ , k]
      
      beta[ , k] <- beta[ , k] - eta * hessian_inv %*% jacobian
      
    }
    
    X_beta <- X %*% beta
    exp_X_beta <- exp(X_beta)
    Pk_mat <- exp_X_beta/rowSums(exp_X_beta)
    
    y_guess <- max.col(Pk_mat)-1
    error_train[i+1] <- (1 - sum(as.numeric(y == y_guess))/n1) * 100
    
    sum_log <- sapply(1:K, \(i){function_beta(i, Pk_mat, y)})
    objective[i+1] <- sum(beta * beta) * (lambda / 2) - sum(sum_log)
    
    Xt_beta <- Xt %*% beta
    exp_Xt_beta <- exp(Xt_beta)
    Pk_test <- exp_Xt_beta/rowSums(exp_Xt_beta)
    
    yt_guess <- max.col(Pk_test)-1
    error_test[i+1] <- (1 - sum(as.numeric(yt == yt_guess))/n2) * 100
  }
  
  ## Return output
  ##########################################################################
  # beta - p x K matrix of estimated beta values after numIter iterations
  # error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
  # error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
  # objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
  return(list(beta = beta, error_train = error_train, error_test = error_test, objective =  objective))
}
# Application of multi-class logistic to letters data
# Load the letter data
#########################
# Training data
letter_train <- read.table("Data/letter-train.txt", header = F, colClasses = "numeric")
Y <- letter_train[, 1]
X <- as.matrix(letter_train[, -1])

# Testing data
letter_test <- read.table("Data/letter-test.txt", header = F, colClasses = "numeric")
Yt <- letter_test[, 1]
Xt <- as.matrix(letter_test[, -1])

# [ToDo] Make sure to add column for an intercept to X and Xt

vec_1 <- rep(1, dim(X)[1])
vec_2 <- rep(1, dim(Xt)[1])

X <- cbind(vec_1, X)
Xt <- cbind(vec_2, Xt)

# Source the LR function
source("FunctionsLR.R")

# [ToDo] Try the algorithm LRMultiClass with lambda = 1 and 50 iterations. Call the resulting object out, i.e. out <- LRMultiClass(...)

out <- LRMultiClass(Xt, Yt, X, Y, numIter = 50, eta = 0.1, lambda = 1, beta_init = NULL)
out

# The code below will draw pictures of objective function, as well as train/test error over the iterations
plot(out$objective, type = 'o')
plot(out$error_train, type = 'o')
plot(out$error_test, type = 'o')

# Feel free to modify the code above for different lambda/eta/numIter values to see how it affects the convergence as well as train/test errors

# [ToDo] Use microbenchmark to time your code with lambda=1 and 50 iterations. To save time, only apply microbenchmark 5 times.

library("microbenchmark")
#microbenchmark(LRMultiClass(X, Y, Xt, Yt, numIter = 50, eta = 0.1, lambda = 1, beta_init = NULL), times = 5)

# [ToDo] Report the median time of your code from microbenchmark above in the comments below

# Median time: 2.313614 seconds