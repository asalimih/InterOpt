
# vectorized calculation of sd and cv for all combinations
calc_cv_sd2 = function(wmat, combs, data_norm, ctVal, k, weight_method){
	combs = t(combs) # because the bug that happens in data_norm[combs[,rg]]
	# the data_norm should be already normalized
	wmat = as.matrix(wmat) #this is so crucial

	if(ctVal)
		data_norm = 2^(mean(data_norm)-data_norm)

	# the cv of each batch is calculated vectorizedly :)
	bsize = 128
	nn = nrow(wmat)
	bNum = ceiling(nn/bsize)
	startInd = 1
	endInd = min(nn,bsize)

	cvs = c()
	sds = c()
	for(bInd in 1:bNum) {
		rg = c(startInd:endInd)
		len = length(rg)
		if(weight_method=="arith_cv" | weight_method=="arith" | weight_method=="arith_sd"){# arithmatic
			A = as.vector(t(wmat[rg,])) * data_norm[combs[,rg],]
			A = rowsum(A, rep(1:len,each=k))
		}else{											#geometric
			A = as.vector(t(wmat[rg,])) * log(data_norm[combs[,rg],])
			A = exp(rowsum(A, rep(1:len,each=k)))
		}
		# each row in A is a reference made of a combination
		sds = c(sds,matrixStats::rowSds(log2(A)))
		cvs = c(cvs,matrixStats::rowSds(A)/rowMeans(A))

		startInd = startInd+bsize
		endInd = min(startInd+bsize-1,nn)
		if(startInd>endInd) break;
	}
	return(cbind(CV=cvs, SD=sds))
}


comb_weights2 = function(data, ctVal, k, weight_method="arith", sub_ind=NULL, lim = NULL, combs = NULL, weights_from_raw = F, mc.cores=32){

	if(is.null(combs))
		stop('combs cant be null!')

	n = nrow(data)

	if(ctVal) {
		data_ct = data
		data = 2^(mean(data)-data) # the mean(data) is important because it prevents very small numbers in geom_cv equations
	}else{
		data_ct = log2(1+data)
	}
	# data is the actual expression (cpm or 2^(40-ct))
	# data_ct is log2 scale form (log2(cpm+1) or ct)
	weights = c()
	lim = ifelse(is.null(lim),nrow(combs),lim)
	if(lim>nrow(combs))
		lim = nrow(combs)

	if(k>=2 & weight_method=="geom_cv"){
		data_log = log(data)
	}

	calc_single_comb_weights = function(i){
		if(weight_method=="arith_cv")
			w = arith_cv(data[combs[i,],])
		else if(weight_method=='arith_sd')
			w = arith_sd(data[combs[i,],])
		else if(weight_method=="arith" | weight_method=="geom")
			w = rep(1/k, k)
		else if(weight_method=="geom_cv_exh")
			w = dirty_g(data[combs[i,],])
		else if(weight_method=="geom_cv"){
			w = geom_cv(data_log[combs[i,],])
		}else if(weight_method=="random"){
			rr = stats::runif(k)
			w = rr/sum(rr)
		}else if(weight_method=="geom_sd"){
			w = geom_sd(data_ct[combs[i,],])
		}else if(weight_method=="geom_sd_soft"){
			w = geom_sd.soft(data_ct[combs[i,],])
		}else if(weight_method=="geom_sd_hybrid"){
			w = geom_sd.hybrid(data_ct[combs[i,],])
		}else if(weight_method=="sd_simple"){
			w = sd_simple(data_ct[combs[i,],])
		}
		if(i%%500==0)
			cat(i,'/', lim, '              \r')
		return(w)
	}


	if(mc.cores<=1){ # Serial mode
		FinalW <<- lapply(seq(1,lim,1), calc_single_comb_weights)
	}else if(Sys.info()['sysname']=='Windows'){ # Parallel on windows
		tryCatch( {
			this.env <- environment()
			cl = parallel::makeCluster(mc.cores, type='PSOCK')
			parallel::clusterExport(cl,  varlist=c('arith_cv','geom_cv', 'dirty_g','geom_sd',
												   'geom_sd.soft','geom_sd.hybrid','arith_sd',
												   'sd_simple', 'f_cv','run_again_geom_cv',
												   'gen_new_x_neg','gen_new_x_pos','cv_grad'),
									envir = this.env)
			ivec = seq(1,lim,1)
			FinalW = parallel::parLapply(cl, X=ivec, calc_single_comb_weights)
		}, error = function(e){
			message('Parallel computation failed, switching to serial mode.')
			FinalW <<- lapply(seq(1,lim,1), calc_single_comb_weights)
		}, finally = {
			parallel::stopCluster(cl)
		})

	}else{ # Parallel on linux
		FinalW <- parallel::mclapply(seq(1,lim,1), calc_single_comb_weights, mc.cores=mc.cores)
	}

	weights = do.call(rbind, FinalW)
	return(weights)
}

write.gct = function(df, destFile){
	writeLines(c("#1.2",paste0(nrow(df), '\t', ncol(df))), destFile)
	suppressWarnings(utils::write.table(cbind(gene = rownames(df), gene_d = rownames(df), df)
										, destFile, sep = '\t', quote = F, row.names = F, append = T))
}

checkInputType = function(data, ctVal){
	if(!is.numeric(data)){
		stop('The input expression data must be a numeric matrix.')
	}
	if(!ctVal & max(data)<30){
		stop('none-CT input data should not be in log2')
	}
	if(ctVal & max(data)>40){
		stop('CT input data has values larger than 40')
	}
	minimum = min(data)
	if(ctVal & minimum<1){
		warning(paste0('the none-CT input data has minimum of ',minimum,'. make the minimum zero by subtracting it from data'))
	}
}

# this is to perform the exact filtering in cuda code filterInputData() so that no filtering is done there
filterInputData = function(data, ctVal, logarithm=F){
	if(ctVal){
		exc = (rowSums(data<0)>0)
	}else{
		if(logarithm){
			exc = (rowSums(data<0)>0)
		} else{
			exc = (rowSums(data<1)>0)
		}
	}
	if(any(exc)){
		message(paste0(sum(exc),' rows from data were removed.'))
		return(data[-which(exc),])
	}else{
		return(data)
	}
}

prep_normalized = function(data, ctVal, data_norm=NULL, norm_method='median_sd', high_exp_thr=35){
	if(!is.null(data_norm))
		data_norm = data_norm[rownames(data),]
	else if(ctVal){
		if(norm_method=='median_sd'){
			sds = matrixStats::rowSds(data)
			thr = stats::quantile(sds, 0.5)
			ref = colMeans(data[sds<thr,])
		}else if(norm_method=='high_exp'){
			ref = apply(data,2, function(x) mean(x[x<high_exp_thr]))
		}else{
			ref = colMeans(data)
		}
		data_norm = t(t(data) - ref)


	}else
		data_norm = data # data is not ct so it must have been normalized before
	return(data_norm)
}


saveFlat = function(df, destFile){
	flat = as.vector(t(df))
	utils::write.table(flat, destFile, quote = F, row.names = F, col.names = F)
}

generate_combs_inds = function(mirs_for_comb, all_mirs, combs_name_mat, k){
	if(!is.null(combs_name_mat)){
		combs_idx_mat = sapply(combs_name_mat,function(x) match(x, all_mirs))
	}else{
		# generate all combinations of mirs_for_comb based on their index in all_mirs
		mirs_for_comb_ind = match(mirs_for_comb, all_mirs)
		combs_idx_mat = t(utils::combn(mirs_for_comb_ind, k))
	}
	#each row is a combination. indices are based on all_mirs
	return(combs_idx_mat)
}

saveGroupVec = function(gr, destFile){
	utils::write.table(paste0(as.numeric(factor(gr))-1,collapse = ''),
					   destFile, quote = F, row.names = F, col.names = F)
}

saveUnstablesIdx = function(data_norm, destFile, genorm_k_stables = 10) {
	inds = which(matrixStats::rowSds(data_norm)>1.2)
	the_sds = matrixStats::rowSds(data_norm)
	inds = setdiff(seq(nrow(data_norm)), order(the_sds)[1:genorm_k_stables])
	#inds = which(matrixStats::rowSds(data_norm)>1.2)
	# one based
	utils::write.table(paste0(inds,collapse = ' '),
					   destFile, quote = F, row.names = F, col.names = F)
}

readValidationResult = function(fileRes){
	df = utils::read.table(fileRes, sep = '\t', header = T)
	return(df$Stability)
}


run_algor_cuda = function(data, gr, ctVal, k, alg, wmethod, comb_num, tmpFolder,
						  tag="source", verbose=T, remove_left_over=T, cuda_kernel='SOURCE/./InterOptCuda')
{
	ExecFile           = cuda_kernel
	COMBINATION_LENGTH = k
	ERROR_FILE_NAME    = file.path(tmpFolder, 'errors.err')
	META_FILE_NAME     = file.path(tmpFolder, 'metas.meta')
	MIRS               = nrow(data)
	FILE_NAME          = file.path(tmpFolder, paste0(tag,'_data_processed.gct'))
	SAMPLES            = ncol(data)
	METHOD             = ifelse(ctVal, 1, 2)
	GROUP_FILE         = file.path(tmpFolder, paste0(tag,'_groupvec.txt'))

	ALGORITHM =  which(c('Genorm', 'NormFinder', 'BestKeeper')==alg)
	COMBINATION_NUMBER = comb_num
	OUTPUT_FILE_NAME   = file.path(tmpFolder, paste0(tag,'_',wmethod,'_',alg,'_result.out'))
	GEOMETRIC          = as.numeric(wmethod=='geom' | wmethod=='geom_cv' | wmethod=='geom_cv_exh' |
										wmethod=='geom_sd' | wmethod=='sd_simple' |
										wmethod=='geom_sd_soft' | wmethod=='geom_sd_hybrid' | wmethod=='random')
	WEIGHT_FILE_NAME   = file.path(tmpFolder, paste0(wmethod,'_flatweights.txt'))
	COMBS_FILE_NAME    = file.path(tmpFolder, paste0(tag, '_flatcombs.txt')) # combs_idx is saved in there beforehand
	UNSTABLES_FILE_NAME= file.path(tmpFolder, paste0(tag, '_unstables.txt'))
	unstable_ids = utils::read.table(UNSTABLES_FILE_NAME, sep='\t', stringsAsFactors=F)$V1
	UNSTABLES_NUM      = length(strsplit(unstable_ids, ' ')[[1]])

	commandStr <- paste(ExecFile, ALGORITHM, COMBINATION_LENGTH, COMBINATION_NUMBER, FILE_NAME,
						OUTPUT_FILE_NAME, ERROR_FILE_NAME, META_FILE_NAME,
						MIRS, SAMPLES, METHOD, WEIGHT_FILE_NAME, COMBS_FILE_NAME, GEOMETRIC,
						UNSTABLES_NUM, UNSTABLES_FILE_NAME,
						paste0('$(cat ',GROUP_FILE,')'), paste0('> ',tmpFolder,'/cuda_log.txt') )

	#cat("Command: ", commandStr,"\n")
	system(commandStr) #run cuda code
	df = utils::read.table(OUTPUT_FILE_NAME, sep = '\t', header = T)
	if(k==1){ # the order of genes changes in the Genorm method
		rownames(df) = df$Name
		df = df[rownames(data),]
	}
	stability =  df$Stability

	if(remove_left_over){
		system(paste("rm", META_FILE_NAME))
		if(file.exists(ERROR_FILE_NAME))
			system(paste("rm", ERROR_FILE_NAME))
		system(paste("rm", file.path(tmpFolder, '*.out')))
	}

	return(stability)
}

# used in iterative mode
next_combMat = function(kBestMat, genes_idx, keep){
	keep = min(nrow(kBestMat),keep)

	combMat = c()
	for(i in 1:keep) {
		exc <- kBestMat[i,] #indexes to exclude from all genes list
		if(i>1) {
			#checking previous candidates to know which index to exclude
			for(pre in 1:i-1) {
				dif <- setdiff(kBestMat[pre,],kBestMat[i,])
				if(length(dif)==1)
					exc = c(exc,dif)
			}
		}
		idxCheck =  setdiff(genes_idx,exc) #indexes to add to this row of kBestMat
		if(length(idxCheck)>0)
		{
			combMat=rbind(combMat,cbind(t(matrix(rep(kBestMat[i,],length(idxCheck)),
												 ncol =length(idxCheck),nrow = ncol(kBestMat))),
										ind=idxCheck))
		}
	}
	return(combMat)
}

#' get weights for all combinations
#'
#' @param data_source a matrix of genes expression. rows are genes and columns are the samples. its \code{rownames} should be the names of the genes
#' @param gr_source a vector of characters showing the groups of the samples
#' @param ctVal_source logical. if \code{TRUE}, the elements in data_source are considered as qPCR CT values. if \code{FALSE}, they are considered as normalized expression of an RNA-seq experiment (count per million)
#' @param tmpFolder a temporary directory to store intermediate files while calculating weights. If not specified an automatic temporary directory is built and used
#' @param sub_names a character vector of gene names to consider in combinations. if \code{NULL}, all genes are considered
#' @param combs_name_mat a matrix of characters. each row shows a combination of genes. if \code{NULL}, all combinations of all genes or a subset of them specified by \code{sub_names} is used.
#' @param sub_samples_for_weights a vector of sample indices to consider for calculating aggregation weights. this argument is used only for benchmarking purposes
#' @param data_target a matrix of external data genes expression. This parameter is used when you want to examine the calculated weights on a separate external data. rows are genes and columns are the samples.
#' @param gr_target a vector of characters showing the groups of the samples of the external data.
#' @param ctVal_target logical. if \code{TRUE}, the elements in data_target are considered as qPCR CT values. if \code{FALSE}, they are considered as normalized expression of an RNA-seq experiment (count per million)
#' @param k integer, the number of genes in each combination.
#' @param iter logical, if \code{TRUE}, an iterative approach is utilized. instead of calculating all combinations for k>2, in each iteration the top most stable combinations (defined by the `keep` parameter) are crossed with other genes to make the new combinations. useful for k>3 when the number of combinations is very high
#' @param keep integer, the number of genes to keep in each iteration of iterative mode
#' @param retain_iters if TRUE, all intermediate iterations results are also reported in the output. (only for iterative mode)
#' @param retain_thr integer, the number of iterations to retain in the output result. (only when \code{retain_iters=TRUE})
#' @param genorm_k_stables integer, number of top stable genes (in terms of standard deviation) to consider in the modified Genorm stability measure. default is 10
#' @param weight_methods a character vector of methods to calculate aggregation weights. available methods are 'arith', 'random','arith_cv','geom','geom_cv', 'geom_cv_exh','geom_sd','geom_sd_soft','geom_sd_hybrid','arith_sd','sd_simple'
#' @param algors a character vector of stability measures. available measures are 'SD', 'CV', 'Genorm', 'NormFinder'
#' @param data_source_norm a matrix of normalized genes expression. if \code{NULL}, the default automatic normalization is used. aggregation weights are calculated based on the normalized data by default unless \code{weight_from_raw=TRUE}. Moreover the stability measures are also calculated based on the normalized data
#' @param data_target_norm a matrix of normalized external genes expression
#' @param norm_method a single character string. 'median_sd' or 'high_exp'. In 'median_sd' method, the average of half of the genes with lower standard deviation is used as a reference gene to normalize the data. In 'high_exp' method only in each sample only the ct values larger than \code{norm_method_exp_thr} are considered.
#' @param norm_method_exp_thr integer, the CT threshold to use in high_exp normalization method.
#' @param weights_from_raw if \code{TRUE}, aggregation weights are calculated based on the raw CT values instead of the normalized data.
#' @param val_on_source if \code{TRUE}, the stability measures are calculated on data_source.
#' @param val_on_target if \code{TRUE}, the stability measures are calculated on data_target.
#' @param verbose logical, print the calculation process in console
#' @param remove_left_over logical, if \code{FALSE}, the intermediate files used for CUDA calculation are not removed in the tmpFolder. (used for debugging purposes)
#' @param saveRDS logical, if \code{TRUE}, the final output is also saved as an RDS file in the tmpFolder
#' @param mc.cores number of the cpu cores to use for calculation SD and CV stability measures.
#' @param cuda_kernel a single character string for the InterOpt cuda kernel executable. defauly is 'InterOptCuda'
#'
#' @return
#' @include weight_utils.R
#' @export
#'
#' @examples
run_experiment = function(data_source, gr_source, ctVal_source, tmpFolder=NULL,
						  sub_names=NULL, combs_name_mat=NULL, sub_samples_for_weights=NULL,
						  data_target=NULL, gr_target=NULL, ctVal_target=NULL,
						  k=2, iter=F, keep=50, retain_iters = F, retain_thr = 10, genorm_k_stables = 10,
						  weight_methods = c('geom_sd_hybrid'),
						  algors = c('SDCV'),
						  data_source_norm=NULL, data_target_norm=NULL, norm_method='high_exp', norm_method_exp_thr=35,
						  weights_from_raw=F, val_on_source=T, val_on_target=T,
						  verbose=T, remove_left_over=T, saveRDS=F, mc.cores=10, cuda_kernel='InterOptCuda')
{
	if(any(!weight_methods %in% c('arith', 'random','arith_cv','geom','geom_cv', 'geom_cv_exh','geom_sd','geom_sd_soft','geom_sd_hybrid','arith_sd','sd_simple')))
		stop('wrong weight_methods element!')
	cuda_algor_flag = any(c('Genorm', 'NormFinder', 'BestKeeper')%in%algors)
	if(is.null(data_target))
		val_on_target=F
	if(k==1){
		weight_methods=c('arith')
		message('k:1, No Optimization, results are in arith.')
	}
	if(k>2 & ("geom_sd_hybrid" %in% weight_methods))
		stop("geom_sd_hybrid is not available for k>2 please set weight_methods parameter to 'geom_sd'")
	if(iter){
		if(k<=2)
			stop('k must be larger than 2 in iterative mode.')
		if(!is.null(combs_name_mat))
			warning('combs_name_mat is ignored in iterative mode')
	}
	if(is.null(sub_samples_for_weights)){
		sub_samples_for_weights = 1:ncol(data_source)
	}
	if(ncol(data_source)!=length(gr_source))
		stop('The number of columns in the data must be equal to the number of elements in gr')
	if(is.null(tmpFolder)){
		tmpFolder = tempdir()
	}
	suppressWarnings(dir.create(tmpFolder, recursive=TRUE))
	## Preprocess
	# 	If normalzied data is not prepared in the input and its type is ct then the average is considered
	# 	as normalized data.

	# Weights are calculated based on data_source_norm unless weights_from_raw
	# CV and SD are calculated based on data_source_norm and data_target_norm
	# (nomfinder and genorm get the data and do the normalization themselves)

	# if sub_names==NULL mirs_for_comb is all the genes which are not filtered
	# if sub_names!=NULL mirs_for_comb is those sub_names which are still available in data_source after intial filter
	# if combs_name_mat!=NULL mirs_for_comb is ignored - all the names in combs_name_mat must be available in filtered data_source


	## Preprocess Source data
	checkInputType(data_source, ctVal_source)
	data_source = filterInputData(data_source, ctVal_source)
	genes_source = rownames(data_source)
	data_source_norm = prep_normalized(data_source, ctVal_source, data_source_norm, norm_method, norm_method_exp_thr)

	## Preprocess Target data
	if(val_on_target){
		checkInputType(data_target, ctVal_target)
		data_target = filterInputData(data_target, ctVal_target)
		genes_target = rownames(data_target)
		data_target_norm = prep_normalized(data_target, ctVal_target, data_target_norm, norm_method, norm_method_exp_thr)
	}

	## Define combination mirs
	if(is.null(sub_names))
		mirs_for_comb = genes_source
	else
		mirs_for_comb = intersect(sub_names, genes_source) # sub_name must come first so the results would have sub_name order
	if(val_on_target)
		mirs_for_comb = intersect(mirs_for_comb, genes_target)

	## Save required files
	if(cuda_algor_flag){
		if(val_on_source){
			FILE_NAME  = file.path(tmpFolder, 'source_data_processed.gct')
			GROUP_FILE = file.path(tmpFolder, 'source_groupvec.txt')
			write.gct(data_source, FILE_NAME)
			saveGroupVec(gr_source, GROUP_FILE)
			COMBS_FILE_NAME_S = file.path(tmpFolder, 'source_flatcombs.txt')
			UNSTABLES_FILE_NAME = file.path(tmpFolder, 'source_unstables.txt')
			saveUnstablesIdx(data_source_norm, UNSTABLES_FILE_NAME, genorm_k_stables)
		}
		if(val_on_target){
			FILE_NAME  = file.path(tmpFolder, 'target_data_processed.gct')
			GROUP_FILE = file.path(tmpFolder, 'target_groupvec.txt')
			write.gct(data_target, FILE_NAME)
			saveGroupVec(gr_target, GROUP_FILE)
			COMBS_FILE_NAME_T = file.path(tmpFolder, 'target_flatcombs.txt')
			UNSTABLES_FILE_NAME = file.path(tmpFolder, 'target_unstables.txt')
			saveUnstablesIdx(data_source_norm, UNSTABLES_FILE_NAME, genorm_k_stables)
		}
	}

	proc_combs = function(gene_names, k, destFile) {
		#(must be calculated for source and target seperately but the order is the same as weights)
		combs_idx_mat = generate_combs_inds(mirs_for_comb, gene_names, combs_name_mat, k)
		if(cuda_algor_flag){
			saveFlat(combs_idx_mat-1, destFile)
		}
		return(combs_idx_mat)
	}

	proc_combs_iter = function(kBestMat, genes_idx, keep, destFile) {
		combs_idx_mat = next_combMat(kBestMat, genes_idx, keep)
		if(cuda_algor_flag){
			saveFlat(combs_idx_mat-1, destFile)
		}
		return(combs_idx_mat)
	}

	proc_weights = function(data, ctVal, combs_idx_mat, k, data_norm=NULL) {
		if(verbose)
			cat('Calculating weights(',wmethod,')[k',k,']...                         \r', sep='')
		if(is.null(data_norm) | weights_from_raw){ # data_source_norm can happen to be null just when source data is not CT
			weights_k = comb_weights2(data[,sub_samples_for_weights], ctVal, k, weight_method = wmethod, combs=combs_idx_mat, weights_from_raw=T, mc.cores=mc.cores)
		}else{
			weights_k = comb_weights2(data_norm, ctVal, k, weight_method = wmethod, combs=combs_idx_mat, mc.cores=mc.cores)
		}
		if(verbose)
			cat('(',wmethod,')[k',k,'] weights Calculated!                         \n', sep='')
		if(cuda_algor_flag){
			WEIGHT_FILE_NAME   = file.path(tmpFolder, paste0(wmethod,'_flatweights.txt'))
			saveFlat(weights_k, WEIGHT_FILE_NAME)
		}
		return(weights_k)
	}

	proc_algor = function(data, ctVal, gr, k, data_norm, combs_mat, weights_mat, alg, tag){
		n_comb = nrow(combs_mat)
		if(verbose)
			cat('Calculating measure(',alg,')[k',k,'][',alg,'][',n_comb,'] ',tag,' data ...                         \r', sep='')
		if(alg=='SDCV' | alg=='SD' | alg=='CV'){
			measures = calc_cv_sd2(weights_mat, combs_mat, data_norm, ctVal, k, wmethod)
			if(verbose)
				cat('(',wmethod,')[k',k,'][',alg,'][',n_comb,'] ',tag,' data Done!                         \n', sep='')
			if(iter){ # in this case alg is either SD or CV (because of the separation)
				measure = matrix(measures[,alg])
				colnames(measure) = alg
				return(measure)
			}else
				return(as.matrix(as.data.frame(measures)))
		}else{
			stab = run_algor_cuda(data, gr, ctVal, k, alg, wmethod, n_comb, tmpFolder, tag, verbose, remove_left_over, cuda_kernel)
			if(verbose)
				cat('(',wmethod,')[k',k,'][',alg,'][',n_comb,'] ',tag,' data Done!                         \n', sep='')
			stab = matrix(stab) # just to be able to give a name to the vector so when binded to other vectors it would keep it
			colnames(stab) = alg
			return(stab)
		}

	}

	proc_resdf = function(gene_names, combs_mat, weights_mat, stabs_mat, k){
		resdf = data.frame(apply(combs_mat, 2, function(x) gene_names[x]),
						   weights_mat, stabs_mat, stringsAsFactors = F)
		colnames(resdf) = c(paste0("Gene",1:k), paste0("w",1:k), colnames(stabs_mat))
		rownames(resdf) = c(1:nrow(resdf))
		return(resdf)
	}


	res_source = list()
	res_target = list()

	#--------- Main --------
	if(!iter){
		## Combinations
		combs_mat_source = proc_combs(genes_source, k, COMBS_FILE_NAME_S)
		if(val_on_target)
			combs_mat_target = proc_combs(genes_target, k, COMBS_FILE_NAME_T)

		for(wmethod in weight_methods){
			## Weights
			weights_mat = proc_weights(data_source, ctVal_source, combs_mat_source, k, data_source_norm) # weights_mat is accessed in functions

			## Algors
			stabs_mat_source = c()
			stabs_mat_target = c()
			for(alg in algors){
				if(val_on_source){
					stab = proc_algor(data_source, ctVal_source, gr_source, k, data_source_norm, combs_mat_source, weights_mat, alg, 'source')
					stabs_mat_source = cbind(stabs_mat_source, stab)
				}
				if(val_on_target){
					stab = proc_algor(data_target, ctVal_target, gr_target, k, data_target_norm, combs_mat_target, weights_mat, alg, 'target')
					stabs_mat_target = cbind(stabs_mat_target, stab)
				}
			}

			## Results DF
			if(val_on_source){
				res_source[[wmethod]] = proc_resdf(genes_source, combs_mat_source, weights_mat, stabs_mat_source, k)
			}
			if(val_on_target){
				res_target[[wmethod]] = proc_resdf(genes_target, combs_mat_target, weights_mat, stabs_mat_target, k)
			}
		}
	}

	#---------Iterative---------
	if(iter){
		if('SDCV'%in%algors){
			i = which(algors=='SDCV')
			algors = algors[-i]
			algors = append(algors, c('SD','CV'), after=(i-1))

		}
		genes_source_idx = match(mirs_for_comb, genes_source)
		if(val_on_target)
			genes_target_idx = match(mirs_for_comb, genes_target)
		for(wmethod in weight_methods){

			if(val_on_source)
				res_source[[wmethod]] = list()
			if(val_on_target)
				res_target[[wmethod]] = list()

			for(alg in algors){
				if(val_on_source){
					res_source[[wmethod]][[alg]] = list()
					res_source[[wmethod]][[alg]][['stats']] = c()
				}
				if(val_on_target){
					res_target[[wmethod]][[alg]] = list()
					res_target[[wmethod]][[alg]][['stats']] = c()
				}

				## Weights 2
				# although Weights 2 is calculated here it should be outside the algors loop. but because the k_iter loop
				# destroys the weights that are saved for now this can't be done. one solutions is to save weights tagged
				# by the k. but the problem then would be that in k_iter loop every steps weights would be saved separately.

				for(k_i in seq(2, k)){
					## Combination
					if(k_i==2)
						combs_mat_source = proc_combs(genes_source, 2, COMBS_FILE_NAME_S) #full 2 combination
					else
						combs_mat_source = proc_combs_iter(kBestMat, genes_source_idx, keep, COMBS_FILE_NAME_S)
					## Weights
					weights_iter = proc_weights(data_source, ctVal_source, combs_mat_source, k_i, data_source_norm)
					## Score
					stab_source_iter = proc_algor(data_source, ctVal_source, gr_source, k_i, data_source_norm, combs_mat_source, weights_iter, alg, 'source')
					if(val_on_target){
						combs_mat_target = apply(combs_mat_source, 2, function(x) genes_target_idx[match(x, genes_source_idx)])
						if(cuda_algor_flag)
							saveFlat(combs_mat_target-1, COMBS_FILE_NAME_T)
						stab_target_iter = proc_algor(data_target, ctVal_target, gr_target, k_i, data_target_norm, combs_mat_target, weights_iter, alg, 'target')
					}

					# Save Stats
					if(val_on_source)
						res_source[[wmethod]][[alg]][['stats']] = rbind(res_source[[wmethod]][[alg]][['stats']],
																		c(Min=min(stab_source_iter), Mean=mean(stab_source_iter)))
					if(val_on_target)
						res_target[[wmethod]][[alg]][['stats']] = rbind(res_target[[wmethod]][[alg]][['stats']],
																		c(Min=min(stab_target_iter), Mean=mean(stab_target_iter)))

					## New Combinations
					sorted_order = order(stab_source_iter)
					kBestMat = combs_mat_source[sorted_order,]

					## Results
					if(k>2)
						sorted_order = sorted_order[1:keep]
					if((retain_iters & k_i<=retain_thr) | k_i==k){
						if(val_on_source){
							res_source[[wmethod]][[alg]][[paste0('k',k_i)]] = proc_resdf(genes_source,combs_mat_source[sorted_order,],
																						 weights_iter[sorted_order,], stab_source_iter[sorted_order], k_i)
						}
						if(val_on_target){
							res_target[[wmethod]][[alg]][[paste0('k',k_i)]] = proc_resdf(genes_target, combs_mat_target[sorted_order,],
																						 weights_iter[sorted_order,], stab_target_iter[sorted_order], k_i)
						}
					}
				}
			}
		}
	}

	## Integrate results
	res_exper = list(res_source=res_source, res_target=res_target, genes=mirs_for_comb)
	if(saveRDS)
		saveRDS(res_exper, file.path(tmpFolder, 'res.rds'))

	if(remove_left_over & cuda_algor_flag) {
		system(paste("rm", file.path(tmpFolder, '*.gct')))
		system(paste("rm", file.path(tmpFolder, '*.txt')))
		if(verbose)
			message('Left over files removed.')
	}
	return(res_exper)
}

