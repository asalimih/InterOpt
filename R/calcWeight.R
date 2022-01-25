#' Calculate weights for a set of internal controls
#'
#' @param x A matrix with each row corresponding to an internal control
#'
#' @return a vector containing the weights of each internal control
#' @export
#'
#' @examples
#' x <- matrix(1:6, 2, 3)
#' w <- calcWeight(x)
calcWeight = function(x){
	# some comments
	d = nrow(x)
	return(rep(1/d, d))
}
