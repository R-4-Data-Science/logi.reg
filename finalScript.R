
#outline for package code from class

set.seed(123)
x<-matrix(rnorm(600), nrow = 200)
epsilon <- rnorm(200, 0, sd = 0.25)
x <- cbind(rep(1, 200),x)
beta <- c(2, -3, 4, 3) 
y <- x%*%beta + epsilon

beta_ls <- function(response, predictors, beta) {
  
  t(response - predictors%*%beta)%*%(response - predictors%*%beta)
  
}

beta_ls(response = y, predictors = x, beta = beta)

optim(rep(0, 4), fn = beta_ls, response = y, predictors = x)

solve(t(x)%*%x)%*%t(x)%*%y