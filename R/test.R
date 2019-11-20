#' test
#'
#' @param X scaler to be squared
#'
#' @return the square of X
#' @export
#'
#' @examples
#' temp = test(4)
test = function(X)
{
  temp = timesTwo(X)
  return(temp)
}

#' test1
#'
#' @param z a thing
#' @param y  another thing
#'
#' @return returns the sum of things
#' @export
#'
#' @examples
#' temp = test1(1,2)
test1 = function(z, y)
{
  return(z + y)
}


#' Title
#'
#' @param x 
#' @param y 
#'
#' @return the product of 
#' @export
#'
#' @examples
test2 =function(x,y)
{
  return(x*y)
}