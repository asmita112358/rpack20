#' d_covR
#'
#' @param X - vector of random samples fron the distribution of X
#' @param Y - vector of random samples from the distribution of Y
#' @param N - scalar positive integer valued parameter to be supplied by the user, must be less than sqrt(length(x))
#'
#' @return distance covariance between X and Y via a fast algorithm
#' @export
#'
#' @examples
d_covR = function(X,Y,N)
{
  if(N > sqrt(length(X)))
    stop(print("N is supposed to be less than square root of length of X"))
  return(d_cov(X, Y, N))
}