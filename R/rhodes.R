#' Logistic Curve Fitting by Rhodes Method
#' @export
#' @param p a numeric vector
#' @author Arnab Roy, Debarghya Baul.
#' @importFrom graphics legend matplot
#' @description This function  fits the Logistic Curve in population Data by Rhodes Method along with estimates of the parameters and predicted value.
#' @details Suppose we have n observations from population size corresponding to n equivalent time points say, at t=0,1,...,n-1. Here we assume the Logistic law of population growth , p=L/(1+exp(r*(beta-t))).
#' @return r.hat , L.hat , beta.hat   :    the estimated values of the parameters r, L and beta.
#' @return  predicted.values  :   the predicted values of p
#' @examples  u=c(12,15,16,18,16,21,25,27,29,30,35,36)
#' @examples  rhodes.curve(u)






rhodes.curve=function(p){
  n=length(p)
  t=0:(n-1)
  x=(1/p)[-n]
  y=(1/p)[-1]

  B_hat=sqrt(sum((y-mean(y))^2)/sum((x-mean(x))^2))
  A_hat=mean(y)-B_hat*mean(x)

  r_hat=-log(B_hat)
  L_hat=(1-B_hat)/A_hat

  beta_hat=(1/(n*r_hat))*sum(log((L_hat/p)-1))+(n-1)/2

  p_hat=L_hat/(1+exp(r_hat*(beta_hat-t)))

  matplot(t,cbind(p,p_hat),xlab = "YEAR" , ylab="y",main="LOGISTIC CURVE (RHODES METHOD)"

          , type="l",lwd=c(2,3),lty=c(1,2),col=c(1,2))
  legend("topleft",legend=c("Actual Plot", "Logistic Curve"),
         lty=c(1,2),
         col=c(1,2))

  list(r.hat=r_hat,L.hat=L_hat,beta.hat=beta_hat,predicted.values=p_hat)


}
