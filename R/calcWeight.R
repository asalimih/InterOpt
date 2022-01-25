#' Calculate weights for a set of internal controls
#'
#' @param x a matrix with each row corresponding to an internal control
#' @param ctVal a logical specifining the input matrix being in qPCR CT values (T) or normalized high-throughput data (F)
#' @param weight_method a string defining the optimization method to calculate weights.
#' choose one of 'geom_sd', 'geom_cv', 'arith_sd', 'arith_cv', 'arith', 'geom', 'random'
#'
#' @return a vector containing the weights of each internal control
#' @export
#'
#' @examples
#' x <- matrix(1:6, 2, 3)
#' w <- calcWeight(x, TRUE, 'geom')
calcWeight = function(x, ctVal=T, weight_method='geom_sd'){
	data = x
	k = nrow(data)

	if(ctVal) {
		data_ct = data
		data = 2^(mean(data)-data) # the mean(data) is important because it prevents very small numbers in clean_g equations
	}else{
		data_ct = log2(1+data)
	}

	if(weight_method=="arith_cv")
		w = arith_cv(data)
	else if(weight_method=='arith_sd')
		w = arith_sd(data)
	else if(weight_method=="arith" | weight_method=="geom")
		w = rep(1/k, k)
	else if(weight_method=="geom_cv")
		w = geom_cv(log(data))
	else if(weight_method=="random"){
		rr = stats::runif(k)
		w = rr/sum(rr)
	}else if(weight_method=="geom_sd")
		w = geom_sd(data_ct)
	else if(weight_method=="sd_simple")
		w = sd_simple(data_ct)

	return(w)
}
