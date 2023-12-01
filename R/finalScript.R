#' @title Logistic Regression Using Numerical Optimization
#'
#' @description Estimates a beta coefficient from a set of explanatory and response variables
#' with basic logistic regression using numerical optimization. Also displays initial coefficients used for
#' the regression from the least squares formula.
#' @param X a numeric dataset containing explanatory variables
#' @param y a vector containing response variable
#' @return A \code{list} containing the following attributes:
#' \describe{
#'      \item{Initial}{Initial coefficients from least squares}
#'      \item{Estimated}{Estimated coefficients from logistic regression}
#' }
#' @importFrom stats
#' @export
#' @examples
#' data=read.csv("iris_csv.csv")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#' dt = sort(sample(nrow(data), nrow(data)*.7))
#' train<-data[dt,]
#' X=as.matrix(train[,1:n_features-1])
#' y=train[,n_features]
#' logreg_optim <- function(X, y)
logreg_optim <- function(X, y) {

  # Calculate initial coefficients from least-squares formula
  init_weights <- solve(t(X)%*%X)%*%t(X)%*%y

  # Add intercept term to X
  X <- cbind(1, X)

  # Begin objective function
  logi_obj <- function(init_weights) {

    # Probability function
    probs <- 1 / (1 + exp(-X%*%init_weights))

    # Logical regression function
    cost <- (1/nrow(X))*sum((-y*log(probs)) - ((1-y)*log(1-probs)))

    return(cost)
  }

  # Optimization using the logistic regression objective function
  optimal <- optim(par = init_weights, fn = logi_obj, X=X,y=y)

  # Report estimated optimal coefficients
  est_weights <- optimal$par

  out = list(Initial = init_weights, Estimated = est_weights)
  out
}



#' Bootstrap Confidence Interval
#' @description calculates Bootstrap confidence intervals
#' @param x a numeric dataset containing explanatory variables
#' @param y a vector containing response variable
#' @param alpha the designated significance level for the intervals
#' @param n the number of bootstraps (default is 20)
#'
#' @return a numeric matrix containing the bootstrap confidence interval
#' \item{beta}{estimates of the logistic regression model with different samples}
#' \item{CI}{confidence intervals for each variable and intercept}
#' @export
#'
#' @examples
#'
#' data=read.csv("iris_csv.csv")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#' dt = sort(sample(nrow(data), nrow(data)*.7))
#' train<-data[dt,]
#' x=as.matrix(train[,1:n_features-1])
#' y=train[,n_features]
#'
#' bootstrapCI(x,y,alpha=0.05,n=20)
#'

bootstrap = function(x,y,alpha,n=20){
  data=data.frame("y"=y,"x"=x)
  n_features=ncol(data)
  beta=matrix(nrow=n,ncol=n_features)
  for(i in 1:n){
    dt = data[sample(nrow(data), size=nrow(data),replace=TRUE),]
    x_boot=dt[,2]
    y_boot=dt[,1]
    beta[i,]=logreg_optim(x_boot,y_boot)
  }
  CI=matrix(nrow=n_features,ncol=3)
  i=0
  for(i in 1:n_features){
    b1<-boot(beta[i,],function(u,i) mean(u[i]),R=n)
    g=boot.ci(b1,type=c("norm","basic","perc"),conf=1-alpha)$norm
    CI[i,]=as.vector(g)
  }
  rownames(CI)=c("intercept",colnames(x))
  colnames(CI)=c("Level","Lower","Upper")
  return(CI)
}
# Number of boostrap replications
B <- 500

# Compute the length of vector
n <- length(student_work)

# Confidence level
alpha <- 0.05

# Initialisation of
boot_mean <- rep(NA, B)

# Step 1
for (i in 1:B){
  # Step 2
  student_work_star <- student_work[sample(1:n, replace = TRUE)]

  # Step 3
  boot_mean[i] <- mean(student_work_star)
}

# Step 4
quantile(boot_mean, c(alpha/2, 1 - alpha/2))

#' Logistic curve plot
#' @description show a visual of logistic regression values against a variable.
#' @param ytest the numeric y values of the test dataset
#' @param xtest the numeric x values of the test dataset
#' @param i the column of x values to be graphed
#'
#' @return a line plot of the y values and one independent variable (indicated by `i`)
#' @export
#'
#' @examples
#' data=read.csv("iris_csv.csv")
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
#'
#' glm_plot(ytest,xtest,2)
#'
glm_plot<- function(x, y,i)
{
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
  plot(y ~ x1, col="steelblue", main = "Logistic Regression Curve", xlab= "X", ylab = "p")
  lines(y ~ x,newdata)
}
#' Confusion Matrix
#' @description confusion matrix and relevent metrics with cutoff feature.
#' @param y the predicted y values of the test dataset
#' @param x the numeric y values of the test dataset
#' @param cutoff the cutoff value which by default is 0.5
#' @return a list of confusion matrix and metrics such a accuracy, specificity etc.
#' @export
#'
#' @examples
#' data=read.csv("iris_csv.csv")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#'
#' n_features=ncol(data)
#' dt = sort(sample(nrow(data), nrow(data)*.7))
#' train<-data[dt,]
#' test<-data[-dt,]
#' x=as.matrix(train[,1:n_features-1])
#' y=train[,n_features]
#' xtest=test[,1:n_features-1]
#' ytest=test[,n_features]
#' conf_matrix(x,y)

conf_matrix=function(x,y,cutoff=0.5){
  coeff=logreg_optim(x,y)
  weights=coeff[2:length(coeff)]
  intercept=coeff[1]
  pred=as.vector(weights%*%t(x)+intercept)
  pred_prob = 1/(1+exp(-pred))
  pred_class=ifelse(pred_prob>cutoff,1,0)
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
#' metric Plot
#' @description provides any required metric plot against a range of cutoff values.
#' @param y the predicted y values of the test dataset
#' @param x the numeric y values of the test dataset
#' @param metric is the name of mr=etric for which the plot is required, like Accuracy.
#' @return any required metric plot against a range of cutoff values from 0.1-0.9
#' @export
#'
#' @examples
#' data=read.csv("iris_csv.csv")
#' data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
#' data$class=ifelse(data$class=="Iris-setosa",1,0)
#'
#' n_features=ncol(data)
#' dt = sort(sample(nrow(data), nrow(data)*.7))
#' train<-data[dt,]
#' test<-data[-dt,]
#' x=as.matrix(train[,1:n_features-1])
#' y=train[,n_features]
#' xtest=test[,1:n_features-1]
#' ytest=test[,n_features]
#' metricplot(x,y,metric="Accuracy")

metricplot=function(x,y,metric){
  metric_data=matrix(nrow=9,ncol=7)
  for(i in 1:9){
    metric_data[i,]=unlist(conf_matrix(x,y,cutoff=i/10)[2])
  }
  colnames(metric_data)=c("Prevalence","Accuracy","Sensitivity","Specificity","False Discovery Rate","Diagnostic Odds Ratio","cutoff")
  return(plot(y=metric_data[,metric],x=metric_data[,7],type="l",xlab="cutoff",ylab=metric,main=paste("Cutoff vs",metric)))
}


data=read.csv("~/Desktop/iris_csv.csv")
data=data[data$class=="Iris-setosa" | data$class=="Iris-virginica",]
data$class=ifelse(data$class=="Iris-setosa",1,0)

n_features=ncol(data)
dt = sort(sample(nrow(data), nrow(data)*.7))
train<-data[dt,]
test<-data[-dt,]
x=as.matrix(train[,1:n_features-1])
y=train[,n_features]
xtest=test[,1:n_features-1]
ytest=test[,n_features]
metricplot(x,y,metric="Accuracy")


