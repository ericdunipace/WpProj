#include <RcppEigen.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List timesTwo(Eigen::MatrixXd x, Eigen::MatrixXd theta, Eigen::MatrixXd Y) {
  Eigen::MatrixXd mu = x*theta;
  for(int i = 0; i < mu.size(); i ++){
    if(mu(i) != Y(i)) Rcpp::Rcout << i << ", ";
  }
  return Rcpp::List::create(mu, Y);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
n <- 100
p <- 10
s <- 1000

x <- matrix(rnorm(p*n), nrow=n, ncol=p)
beta <- (1:10)/10
y <- x %*% beta + rnorm(n)

#posterior
prec <- crossprod(x) + diag(1,p,p)*1
mu_post <- solve(prec, crossprod(x,y))
alpha <- 1 + n/2
beta <- 1 + 0.5 * (crossprod(y) + t(mu_post) %*% prec %*% mu_post )
sigma_post <- 1/rgamma(s, alpha, 1/beta)
theta <- sapply(sigma_post, function(ss) mu_post + t(chol(ss * solve(prec))) %*% matrix(rnorm(p, 0, 1),p,1))
mu <- x %*% theta
mu_e <- timesTwo(x, theta, mu)
all.equal(mu_e[[1]],mu)
all.equal(mu_e[[2]],mu)
all.equal(mu_e[[1]], mu_e[[2]])
*/
