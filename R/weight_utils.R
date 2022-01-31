geom_sd = function(mat){
	mat = mat-rowMeans(mat)
	mat2 = t(unlist(mat[1,])-t(mat))
	mat2[1,] = mat[1,]
	rownames(mat2) = as.character(1:nrow(mat2))
	rownames(mat2)[1]='y'
	fit = stats::lm(y~.+0,data=data.frame(t(mat2)))
	w = as.vector(c(1-sum(fit$coefficients, na.rm=T),fit$coefficients))
	w[is.na(w)] = 0
	w = w/sum(w)
	return(w)
}

arith_cv = function(M){
	n = ncol(M)
	A = rep(1, n)
	w_greg = MASS::ginv( (M-rowMeans(M))%*%t(M) )%*%(as.matrix(rowSums(M)))
	w_greg = c(w_greg/sum(w_greg))

	negs = which(w_greg<0)
	poses = which(w_greg>=0)
	if(length(negs)>=1){# the problem with negative w is it can cause negative expression on different datas.
		w_greg[negs] = 0
		if(length(poses)>1)
			w_greg[poses] = arith_cv(M[poses,])
		w_greg = c(w_greg/sum(w_greg))
	}
	return(w_greg)
}

cv_grad = function(x,M){
	n = ncol(M)
	d = nrow(M)
	alpha = sum(x)
	w = x/alpha
	# next line produces NaN when sum(x) is a small number which causes w to have big numbers like 30
	y = exp(colSums(w*M)) # y is column matrix but here is a vector
	beta = sum(y)
	meany = mean(y)
	z = y-meany
	phi = stats::sd(y)/meany
	#Q = eachrow(M, y, "*")
	Q = M * outer(rep.int(1L, nrow(M)), y)
	#Q = M%*%diag(y) # MY
	#g = (n/(phi*alpha*(beta^3))) * (diag(d)-matrix(rep(w,d),nrow=d,byrow=T)) %*% (beta*(Q-rowMeans(Q))-rowSums(Q)%*%t(z)) %*% z
	g = (n/(phi*alpha*(beta^3))) * (  t(diag(d)-w) %*% (  (beta*(Q-rowMeans(Q))-tcrossprod(rowSums(Q),z)) %*% z  )  )
	return(list(grad=as.vector(g), cv=phi))
}

geom_cv = function(M, resistance = 1e-07, tuner = 0.001, x0=NULL, times_alg_ran=0) {

	if(times_alg_ran>=100){
		warning('wgmean not converged in a combination. gmean method used instead')
		return(rep(1/nrow(M),nrow(M)))
	}
	if(is.null(x0))
		x0 = rep(1/nrow(M),nrow(M))
	gradres0 = cv_grad(x0, M)
	g0 = gradres0[['grad']] # maybe its needed to check if g0 is small enough we can just report x0
	cv0 = gradres0[['cv']]
	w0 = x0/sum(x0)
	if(is.na(g0[1]))
		return(run_again_geom_cv(M, resistance, tuner, times_alg_ran+1))

	g0_norm = (0.1/norm(g0,'2'))
	# if cv0 of the combination becomes too small this `step` becomes small then
	# x wont change that much gradient wont change jimb = 0 the next `step` becomes NaN
	step = (tuner*cv0)/sum(g0^2)
	if(abs(step)>g0_norm)
		step = g0_norm
	x_cur = x0-step*g0
	x_cur = x_cur/sum(x_cur)
	x_pre = x0
	g_pre = g0
	cv_cur = cv0 # changing -1 to cv0 solved the NaN problem above
	for(i in seq(500)){
		gradres = cv_grad(x_cur, M)
		g_cur = gradres[['grad']]

		if(i==500 | is.na(cv_cur) | is.na(g_cur[1]))
			return(run_again_geom_cv(M, resistance, tuner, times_alg_ran+1))

		if(abs(cv_cur-gradres[['cv']])<=resistance)
			break

		cv_cur = gradres[['cv']]
		jimb = g_cur-g_pre # if an error occured it may be because jimb became 0
		A = as.vector(crossprod(x_cur-x_pre,jimb))
		g_cur_norm = 0.1/norm(g_cur,'2')

		step = A/sum(jimb^2)
		if(abs(step)>g_cur_norm)
			step = g_cur_norm

		x_next = x_cur-step*g_cur
		x_next = x_next/sum(x_next)
		x_pre = x_cur
		x_cur = x_next
		g_pre = g_cur
	}
	return(x_cur/sum(x_cur))
}

gen_new_x_pos = function(k){
	new_x0 = stats::runif(k)
	new_x0 = new_x0/sum(new_x0)
	if(max(abs(new_x0))>5)
		return(gen_new_x_pos(k))
	return(new_x0)
}

gen_new_x_neg = function(k){
	new_x0 = (stats::runif(k)*2)-1
	new_x0 = new_x0/sum(new_x0)
	if(max(abs(new_x0))>5)
		return(gen_new_x_neg(k))
	return(new_x0)
}

f_cv = function(x){
	return(stats::sd(x)/mean(x))
}

# g stands for gmean
dirty_g = function(mat){
	vec_cv = c()
	min_cv = 99
	w_opt = NULL
	for(i in -500:500) {
		ww = c(0.01*i,1-0.01*i)
		if(any(mat<=0))
			stop('zero or negative expression not allowed')
		newRef <- exp(colSums(ww*log(mat)))
		cur_cv = f_cv(newRef)
		if(min_cv>cur_cv){
			min_cv = cur_cv
			w_opt = ww
		}
		vec_cv = c(vec_cv, cur_cv)
	}
	return(w_opt)
}

run_again_geom_cv = function(M, resistance = 1e-07, tuner = 0.001, times_alg_ran=0){
	k = nrow(M)
	if(k==2 & times_alg_ran>20)
		return(dirty_g(exp(M)))

	if(times_alg_ran>15)
		new_x0 = gen_new_x_neg(k)
	else
		new_x0 = gen_new_x_pos(k)

	if(times_alg_ran>10)
		tuner = ((stats::runif(1)*9)+1)/100

	return(geom_cv(M, resistance, tuner, new_x0, times_alg_ran))
}

#imput should be in linear form (not ct or log2)
arith_sd = function(mat){
	if(nrow(mat)>2)
		stop('arith_sd only works for a combination of two internal controls')
	vec_sd = c()
	min_sd = 99
	w_opt = NULL
	for(i in 0:100) {
		ww = c(0.01*i,1-0.01*i)
		if(any(mat<=0))
			stop('zero or negative expression not allowed')

		newRef <- log2(colSums(ww*mat))
		cur_sd = stats::sd(newRef)
		if(min_sd>cur_sd){
			min_sd = cur_sd
			w_opt = ww
		}
		vec_sd = c(vec_sd, cur_sd)
	}
	return(w_opt)
}

sd_simple = function(mat){
	sds = matrixStats::rowSds(mat)
	w = 1/sds
	w = w/sum(w)
	return(w)
}

