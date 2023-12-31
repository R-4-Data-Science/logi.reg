#' @title Logistic Regression Using Numerical Optimization
#'
#' @description Estimates a beta coefficient from a set of explanatory and response variables
#' with basic logistic regression using numerical optimization. Also displays initial coefficients used for
#' the regression from the least squares formula.
#' @param x a numeric vector or matrix containing explanatory variables
#' @param y a numeric vector containing response variable
#' @return a \code{vector} that contains estimated coefficients of the logistic
#' regression model
#' @importFrom stats runif
#' @export
#' @examples
#' data <- read.csv("/Users/ebyrne1/Desktop/iris_csv.csv")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
# 'data$class=ifelse(data$class=="Iris-setosa",1,0)
# 'x <- as.matrix(data[,1:4])
# 'y <- data[,5]
#' logreg_optim(x, y)
logreg_optim <- function(x, y, alpha, n = 20) {

  # Calculate initial coefficients from least-squares formula
  init_weights <- function(x, y){
    weight <- solve(t(x)%*%x)%*%t(x)%*%y
    weight <- as.vector(weight)
    return(weight)
  }
  # Add intercept term to x
  x <- cbind(1, x)

  # Begin objective function
  logi_obj <- function(init_weights, x, y) {

    # Probability function
    probs <- 1 / (1 + exp(-x%*%init_weights))

    # Logical regression function
    cost <- (1/nrow(x))*sum((-y*log(probs)) - ((1-y)*log(1-probs)))

    return(cost)
  }
  # Optimization using the logistic regression objective function
  optimal <- optim(par = init_weights(x, y), fn = logi_obj, x=x,y=y)

  # Report estimated optimal coefficients
  est_weights <- optimal$par
  return(est_weights)
}


#' Bootstrap Confidence Interval
#' @description Calculates Bootstrap confidence intervals for the beta from \code{logreg_optim}
#' @param alpha the designated significance level for the intervals
#' @param n the number of bootstraps (default is 20)
#' @param x a numeric vector or matrix containing explanatory variables
#' @param y a numeric vector containing response variable
#' @return a numeric matrix containing the bootstrap confidence interval
#' @export
#' @examples
#' data <- read.csv("~/Desktop/Classes/STAT 6120/Final/iris_csv.numbers")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#' x <- as.matrix(data[,1:4])
#' y <- data[,5]
#' bootstrap_CI(x = x, y, alpha = 0.05,n = 20)
bootstrap_CI = function(x, y, alpha, n=20){

  numvar <- ncol(x)+1
  #Create and empty matrix for the bootstrapped coefficients
  boot_beta <- matrix(NA, nrow = n, ncol = numvar)

  for(i in 1:n){

    #Create a matrix with x new data
    boot_newx <- x[sample(1:n, nrow(x), replace = TRUE),]

    #Create a vector with Y new data
    boot_newy <- y[sample(1:n, length(y), replace = TRUE)]

    #Fill the matrix for the new betas from the regression
    boot_beta <- logreg_optim(boot_newx, boot_newy)
  }

  CI_log = matrix(NA, nrow = n, ncol = 2)
  lower<-(alpha/2)*100
  upper<-(1 - alpha/2)*100
  CI_log[1:n,] <- quantile(boot_beta, c(alpha/2, 1 - alpha/2))
  return(t(CI_log))
}

#' Logistic Regression Curve Plot
#' @description shows a visual of logistic regression values against a variable.
#' @param y the numeric y values of the test dataset
#' @param x the numeric x values of the test dataset
#' @param i the column of x values to be graphed
#' @return a line plot of the y values and one independent variable (indicated by `i`)
#' @export
#' @examples
#' data <- read.csv("~/Desktop/Classes/STAT 6120/Final/iris_csv.numbers")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#' x <- as.matrix(data[,1:4])
#' y <- as.vector(data[,5])
#' glm<- glm_plot(x,y,2)
glm_plot<- function(x, y, i){

  x1=x[,i]
  #fit logistic regression model
  coeff <- logreg_optim(x1,y)

  #define new data frame that contains predictor variable
  newdata <- data.frame(x=seq(min(x1), max(x1),len=500))
  #use fitted model to predict values
  weights=coeff[2]
  intercept=coeff[1]
  pred=as.vector(weights*newdata$x+intercept)
  newdata$y = 1/(1+exp(-pred))

  #plot logistic regression curve
  plot(y ~ x1, col="steelblue", main = "Logistic Regression Curve", xlab= "x", ylab = "p")
  lines(y ~ x,newdata)
}

#' Confusion Matrix
#' @description confusion matrix and relevent metrics with cutoff feature.
#' @param y the predicted y values of the test dataset
#' @param x the numeric y values of the test dataset
#' @param cutoff the cutoff value which by default is 0.5
#' @return a list of confusion matrix and metrics such a accuracy, specificity etc.
#' @export
#' @examples
#' data=read.csv("~/Desktop/Classes/STAT 6120/Final/iris_csv.numbers")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#' n_features=ncol(data)
#' dt = sort(sample(nrow(data), nrow(data)*.7))
#' train<-data[dt,]
#' test<-data[-dt,]
#' x=as.matrix(train[,1:n_features-1])
#' y=train[,n_features]
#' xtest=test[,1:n_features-1]
#' ytest=test[,n_features]
#' conf_matrix(x,y)
conf_matrix=function(x, y, cutoff = 0.5){
  coeff=logreg_optim(x,y)
  weights=coeff[2:length(coeff)]
  intercept=coeff[1]
  pred=as.vector(weights%*%t(x)+intercept)
  pred_prob = 1/(1+exp(-pred))
  pred_class=ifelse(pred>cutoff,1,0)
  confusion_matrix=confusionMatrix(as.factor(pred_class), as.factor(y))
  metric_data=c()
  metric_data[1]=as.vector(confusion_matrix$byClass[8])
  metric_data[2]=as.vector(confusion_matrix$overall[1])
  metric_data[3]=as.vector(confusion_matrix$byClass[1])
  metric_data[4]=as.vector(confusion_matrix$byClass[2])
  metric_data[5]=1-as.vector(confusion_matrix$byClass[3])
  metric_data[6]=metric_data[2]*metric_data[3]/((1-metric_data[2])*(1-metric_data[3]))
  metric_data[7]=cutoff
  names(metric_data)=c("Prevalence","Accuracy","Sensitivity","Specificity","False Discovery Rate","Diagnostic Odds Ratio","cutoff")
  return(list(confusion_matrix$table,metric_data))
}

#' Metric Plot
#' @description provides any required metric plot against a range of cutoff values.
#' @param y the predicted y values of the test dataset
#' @param x the numeric y values of the test dataset
#' @param metric is the name of mr=etric for which the plot is required, like Accuracy.
#' @return any required metric plot against a range of cutoff values from 0.1-0.9
#' @export
#' @examples
#' data <- read.csv("~/Desktop/Classes/STAT 6120/Final/iris_csv.numbers")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#' x <- as.matrix(data[,1:4])
#' y <- as.vector(data[,5])
#' metricplot(x,y,metric="Accuracy")
metricplot=function(x,y,metric){
  metric_data=matrix(nrow=9,ncol=7)
  for(i in 1:9){
    metric_data[i,]=unlist(conf_matrix(x,y,cutoff=i/10)[2])
  }
  colnames(metric_data)=c("Prevalence","Accuracy","Sensitivity","Specificity","False Discovery Rate","Diagnostic Odds Ratio","cutoff")
  return(plot(y=metric_data[,metric],x=metric_data[,7],type="l",xlab="cutoff",ylab=metric,main=paste("Cutoff vs",metric)))
}

n_features=ncol(data)
dt = sort(sample(nrow(data), nrow(data)*.7))
train<-data[dt,]
test<-data[-dt,]
x=as.matrix(train[,1:n_features-1])
y=train[,n_features]
xtest=test[,1:n_features-1]
ytest=test[,n_features]
metricplot(x,y,metric="Accuracy")


