#set of functions to standardise the model outputs and sensibly name things

standardise_jag_out <- function(jags_obj){
	#data frame to hold everything in
	df <- data.frame(stand_name = c('int_pool', 'int_herb', 'int_tree', 'mat_dim', 'life_expect', 'PC_temp', 'log(AI)', 'p_sea', 'life_expect|Tree',
		'PC_temp|Tree', 'log(AI)|Tree', 'p_sea|Tree', 'life_expect x PC_temp' ,'life_expect x log(AI)', 'life_expect x p_sea',
		'PC_temp x log(AI)', 'PC_temp x p_sea', 'log(AI) x p_sea', 'spp_sd', 'global_sd', 'alpha', 'rho', 'realtives', 'neighbours'), 
		used_name = NA, used_struct = NA, struct_ind = NA)
	#list of predictor names used in the models
	pred_names = c('mat_dim[n]', 'life_expect[n]', 'PC_temp[n]', 'AI[n]', 'p_sea[n]', 'life_expect[n]*PC_temp[n]', 'life_expect[n]*AI[n]', 
		'life_expect[n]*p_sea[n]', 'PC_temp[n]*AI[n]', 'PC_temp[n]*p_sea[n]', 'AI[n]*p_sea[n]', 'phy_pred[n]', 'geo_pred[n]', 'ac_pred[n]')
	#read in the model as raw text
	mod = readLines(jags_obj$model.file)
	#find the line in the model that defines the mean of the likliehood and extract the name for the intercept
	eq  = mod[grep('+', mod, fixed = TRUE)]
	if(length(eq) == 2) eq = paste0(eq[1], ' ', eq[2])
	used_pars = strsplit(strsplit(eq, split = '<-')[[1]][2], split = ' + ', fixed = TRUE)[[1]]
	if(used_pars[1] == ' a') df[df$stand_name == 'int_pool', 2:4] <- c('a', 'a', 1)
	if(used_pars[1] == ' a[GF[n]]'){
		df[df$stand_name == 'int_herb', 2:4] <- c('a[1]', 'a', 1)
		df[df$stand_name == 'int_tree', 2:4] <- c('a[2]', 'a', 2)
	}
	#fill in the rest of 
	row_ind = NA
	for(i in pred_names){
		all = strsplit(used_pars[grep(i, used_pars, fixed = TRUE)], split = '*', fixed = TRUE)
		if(length(all) > 0){
			for(j in 1:length(all)){
				#test if this interaction or main effect 
				if(length(all[[j]]) == 2 | length(all[[j]]) == 1){#main effect or pooled growth form
					if(length(grep('b[', all[[j]], fixed = TRUE)) > 0 | length(grep('b', all[[j]], fixed = TRUE)) > 0){# simple main effect
						all_non = strsplit(all[[j]], split = '[n]', fixed = TRUE)
						#specieal cases
						if(all_non[[2]] == 'AI') all_non[[2]] = 'log(AI)'
						if(all_non[[2]] == 'phy_pred') all_non[[2]] = 'realtives'
						if(all_non[[2]] == 'geo_pred') all_non[[2]] = 'neighbours'
						if(all_non[[2]] == 'ac_pred') all_non[[2]] = 'neighbours'
												
						row_ind = grep(all_non[[2]], df$stand_name, fixed = TRUE)[1] #take the first 
						if(is.na(df$used_name[row_ind])) df[row_ind, 2:4] <- c(all_non[[1]], substr(x = all_non[[1]], start = 1, stop = 1), substr(x = all_non[[1]], start = 3, stop = 3))
					}
					if(length(grep('B[', all[[j]], fixed = TRUE)) > 0){#split growth type
						all_non = strsplit(all[[j]], split = '[n]', fixed = TRUE)
						#specieal cases
						if(all_non[[2]] == 'AI') all_non[[2]] = 'log(AI)'
						if(all_non[[2]] == 'phy_pred') all_non[[2]] = 'realtives'
						if(all_non[[2]] == 'geo_pred') all_non[[2]] = 'neighbours'
						if(all_non[[2]] == 'ac_pred') all_non[[2]] = 'neighbours'
					
						row_ind = grep(all_non[[2]], df$stand_name, fixed = TRUE)[1:2] #take the first 2 
						if(is.na(df$used_name[row_ind[1]])) df[row_ind[1], 2:4] <- c(paste0(all_non[[1]], '[1]]')[1], substr(x = all_non[[1]][1], start = 1, stop = 1), substr(x = all_non[[1]][1], start = 3, stop = 3))
						if(is.na(df$used_name[row_ind[2]])) df[row_ind[2], 2:4] <- c(paste0(all_non[[1]], '[2]]')[1], substr(x = all_non[[1]][1], start = 1, stop = 1), substr(x = all_non[[1]][1], start = 3, stop = 3))
					}
				}
				if(length(all[[j]]) == 3){#interaction
					all_non = strsplit(all[[j]], split = '[n]', fixed = TRUE)
					#specieal cases
					if(all_non[[2]] == 'AI') all_non[[2]] = 'log(AI)'
					if(all_non[[3]] == 'AI') all_non[[3]] = 'log(AI)'
					if(all_non[[2]] == 'phy_pred') all_non[[2]] = 'realtives'
					if(all_non[[2]] == 'geo_pred') all_non[[2]] = 'neighbours'
					if(all_non[[2]] == 'ac_pred') all_non[[2]] = 'neighbours'
					
					#paste the names together
					int_name = paste0(all_non[[2]], ' x ', all_non[[3]]) 
					if(all_non[[2]] == 'AI') int_name = paste0('log(', all_non[[2]], ') x ', all_non[[3]])
					if(all_non[[3]] == 'AI') int_name = paste0(all_non[[2]], 'x log(', all_non[[3]], ')')
					row_ind = grep(int_name, df$stand_name, fixed = TRUE)[1] #take the first 
					if(is.na(df$used_name[row_ind])) df[row_ind, 2:4] <- c(all_non[[1]], substr(x = all_non[[1]], start = 1, stop = 1), substr(x = all_non[[1]], start = 3, stop = 3))
				}
			}
		}
	}
	if(length(mod[grep('sd.spp', mod)]) > 0) df[grep('spp_sd', df$stand_name), 2:4] = c('sd.spp', 'sd.spp', 1)
	if(length(mod[grep('alpha', mod)]) > 0) df[grep('alpha', df$stand_name), 2:4] = c('alpha', 'alpha', 1)
	if(length(mod[grep('rho', mod)]) > 0) df[grep('rho', df$stand_name), 2:4] = c('rho', 'rho', 1)
	if(length(mod[grep('sd.global', mod)]) > 0) df[grep('global_sd', df$stand_name), 2:4] = c('sd.global', 'sd.global', 1)

	return(list(df = df, jags_obj = jags_obj))
}

#another function that takes the output of standardise_jag_out() and the model object used to build it and adds crebible intervals for each predictor
CI_std_out <- function(std_jags){
	#seperates out the sub objects
	df = std_jags$df
	jags_obj = std_jags$jags_obj

	#add new coloumns to df
	df$lower_01 = NA
	df$lower_025 = NA
	df$lower_05 = NA
	df$median = NA
	df$mean = NA
	df$upper_05 = NA
	df$upper_975 = NA
	df$upper_99 = NA

	#go through and get all the mean and median for each predictor, can use sinple names since thedata are stored in structures with the same 
	#dimetionality as the predictor names
	split_names = strsplit(df$used_name, split = ',', fixed = TRUE)
	simp_var_names = sapply(split_names, FUN = function(x){
			if(length(x) == 1) return(x[1])
			if(length(x) == 2) return(paste0(x[1], ' ,', gsub(']]', '', gsub(' GF[', '', x[2], fixed = T), fixed = T), ']'))
		})

	for(i in 1:length(simp_var_names)){
		if(!is.na(simp_var_names[i])){
			#use some string parsing to locate the right variables in the jags object
			df$mean[i] = eval(parse(text = paste0('jags_obj$BUGSoutput$mean$', simp_var_names[i]))) 
			df$median[i] = eval(parse(text = paste0('jags_obj$BUGSoutput$median$', simp_var_names[i])))
		} 
	}
	#get quantiles, need to add a dimention to variable names at the data structures have an extra dimention to hold the simulations
	comp_var_names = gsub('[', '[,', simp_var_names, fixed = TRUE)
	for(i in 1:length(comp_var_names)){
		if(!is.na(comp_var_names[i])){
			#use some string parsing to locate the right variables in the jags object
			sims = eval(parse(text = paste0('jags_obj$BUGSoutput$sims.list$', comp_var_names[i])))
			df[i, c('lower_01', 'lower_025', 'lower_05', 'upper_05', 'upper_975', 'upper_99')] <-  quantile(sims, probs = c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99))
		} 
	}
	
	return(df)	
}

custom_jags_object <- function(jags_obj, obs_y){
	quant_df = CI_std_out(standardise_jag_out(jags_obj)) #standard names and quanttiles

	#reconstrcut the simulated and mean predictions from the residual
	mus = matrix(nrow = dim(jags_obj$BUGSoutput$sims.list$residual)[1], ncol = dim(jags_obj$BUGSoutput$sims.list$residual)[2])
	for(i in seq_along(obs_y)) mus[, i] = obs_y[i] - jags_obj$BUGSoutput$sims.list$residual[,i]
	mean_mu = apply(mus, MARGIN = 2, FUN = mean)

	#construct a vector of R^2 values
	SS_tot = sum((obs_y - mean(obs_y))^2)
	SS_res = apply(jags_obj$BUGSoutput$sims.list$residual, MARGIN = 1, FUN = function(x) sum(x^2))
	R_sq = 1 - (SS_res / SS_tot)

	bp_value = sum(jags_obj$BUGSoutput$sims.list$bpvalue)/length(jags_obj$BUGSoutput$sims.list$bpvalue)

	#return a named list with the summaries plus put afew things buried in the jags_obj up on the surface layer
	return(list(quant_df = quant_df, pred_matrix = mus, pred_means = mean_mu, R_sq = R_sq, model_file = jags_obj$model.file, 
		DIC = jags_obj$BUGSoutput$DIC, pD = jags_obj$BUGSoutput$pD, bp_value = bp_value, residual = jags_obj$BUGSoutput$sims.list$residual, 
		obs_y = obs_y, full_obj = jags_obj))
}                                                             

