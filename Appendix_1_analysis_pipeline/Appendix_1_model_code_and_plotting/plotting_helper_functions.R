#Set of helper functions for plotting multipule models

diag_plots <- function(jags_ob_names, output_file = 'diagnositc_summary_plots.pdf', output_loc, working_files_loc, Y){
  setwd(output_loc)
  lmat<-rbind(c(0, 1, 0.9, 1), c(0, 0.5, 0, 0.45), c(0.5, 1, 0, 0.45), c(0, 0.5, 0.45, 0.9), c(0.5, 1, 0.45, 0.9))
  text_offset = 0.25
  pdf(file = output_file, paper = 'a4', height = 14, width = 12)
    for(this_model in jags_ob_names){	
      setwd(working_files_loc)
      print(this_model)
      jag_obj <- custom_jags_object(eval(parse(text = this_model)), Y)
      setwd(output_loc)
      split.screen(lmat)
      screen(1)#a bit of text to say what the model is and some stats on how it performs
	par(mar = c(0.5, 0.5, 0.5, 0.5))
	plot(x = 1:10, y = 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
	text(x = 1, y = 8, labels = paste0('model name: ', this_model), pos = 4)
	text(x = 7, y = 8, labels = paste0('DIC: ', round(jag_obj$DIC, 3)), pos = 4)
	text(x = 1, y = 2, labels = paste0('Bp_value: ', round(jag_obj$bp_value, 3)), pos = 4)
	text(x = 5, y = 2, labels = paste0('pD: ', round(jag_obj$pD, 3)), pos = 4)

      screen(2) #R^2
	plot(density(jag_obj$R_sq), type = 'n', main = 'R squared\nn = 1500', xlab = expression(R^2), bty = 'n')
	polygon(density(jag_obj$R_sq), col = grey(0.6), border = FALSE)

      screen(3) #credible intervals on parameters
	used_inds = which(!is.na(jag_obj$quant_df$used_name))
	pred_names = jag_obj$quant_df$stand_name[used_inds]
	dat_x_min = min(jag_obj$quant_df$lower_01[used_inds])
	dat_x_max = max(jag_obj$quant_df$upper_99[used_inds])
	plot_x_min = dat_x_min - ((dat_x_max - dat_x_min) * text_offset) 
	plot_y_max = length(used_inds)
			  
	#set up blank plotting area
	par(mar = c(4,0,0,0), oma = c(0,0,0,0))
	plot(x = seq(plot_x_min, dat_x_max, length = plot_y_max + 2), y = 0:(plot_y_max + 1), type = "n", bty = "n", xlab = "coefficients", ylab = "", xaxt = "n", yaxt = "n")
	#put alternate grey bars on the plot
	altRow<-seq(plot_y_max, 0, -2)
	rect(xleft = plot_x_min, xright = dat_x_max*1.05, ybottom = altRow - 1, ytop = altRow, border = NA, col = grey(0.95)) 
	text(x = plot_x_min, y = plot_y_max + 0.5, labels = 'Predictor', pos = 4, cex = 1)
	lines(x = c(0, 0), y = c(0, plot_y_max), col = 'red')
	#do names and quantiles 
	for(j in 1:length(used_inds)){
	  text(x = plot_x_min, y = j - 0.5, labels = pred_names[j], pos = 4, cex = 0.7)
	  lines(x = c(jag_obj$quant_df$lower_01[used_inds[j]], jag_obj$quant_df$upper_99[used_inds[j]]), y = c(j - 0.5, j - 0.5), lend = 2, col = grey(0.4))
	  lines(x = c(jag_obj$quant_df$lower_025[used_inds[j]], jag_obj$quant_df$upper_975[used_inds[j]]), y = c(j - 0.5, j - 0.5), lend = 2, lwd = 1.5)
	  points(x = jag_obj$quant_df$mean[used_inds[j]], y = j - 0.5, pch = 19)
	}
	axis(side = 1, at = round(c(seq(dat_x_min, 0, length = 4), seq(0, dat_x_max, length = 4)), 2), cex = 1, tck = 0.02, line = 1)

      screen(4) #residual plot
	plot(jag_obj$pred_means, apply(jag_obj$residual, MARGIN = 2, FUN = mean), main = "Mean residual vs mean fitted value", 
	  xlab = "predicted value", ylab = "residual", frame.plot = FALSE)
	abline(0, 0, col = "black")

      screen(5)#fitted vs predcited, with uncertianty in prediction
	mu_ci <- t(apply(jag_obj$pred_matrix, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))#get uncertianty in the estimated value using credible intervals
	max_y <- max(mu_ci) #gets max and min values for the plot
	min_y <- min(mu_ci)
	plot(jag_obj$obs_y, jag_obj$pred_means, main = 'Observed vs predicted (95%CI)', ylim = c(min_y, max_y), bty = 'n', 
	  type = 'n', xlab = 'observed', ylab = 'predicted')
	x_couples <- cbind(jag_obj$obs_y, jag_obj$obs_y)
	for(j in seq_along(mu_ci[,1])) lines(x = x_couples[j, ], y = mu_ci[j, ], lwd = 0.5, col = grey(0.6))
	points(jag_obj$obs_y, jag_obj$pred_means, cex = 0.2, pch = 19)
	abline(0, 1, col = "red", lwd = 2)
      close.screen(all.screens = TRUE)
    }
  dev.off()
}

#Function to make R_sq plot for each response type shows if we can predict the same for each measure
R_sq_plotter <- function(model_spec, color_palett_fill, color_palett_border, list_jags_obj, file_name = 'R_sq_all_responses.pdf', output_loc){
  number_of_measures = length(list_jags_obj)
  names_of_measures = names(list_jags_obj)
  setwd(output_loc)
  pdf(file = file_name, paper = 'a4', height = 14, width = 12)
    par(mfrow = c(number_of_measures, 1))
    for(i in 1:number_of_measures){
      list_R2_density <- lapply(list_jags_obj[[i]], FUN = function(x) density = density(x$R_sq))
      plot_max = max(sapply(list_R2_density, FUN = function(x) max(x$y)))
      plot(x = seq(0, 1, length = 10), y = seq(0, plot_max, length = 10), type = 'n', bty = 'n', xlab = expression(R^2), ylab = '', main = '', tck = 0.015, cex.axis = 2, cex.lab = 2.5, 
	ylim = c(0, plot_max), xlim = c(0, 0.85), yaxt = 'n')
      for(j in 1:length(list_R2_density)){
	polygon(list_R2_density[[j]], col = color_palett_fill[j], border = color_palett_border[j])
      }
      text(x = 0, y = plot_max, labels = names_of_measures[i], pos = 4)
      if(i == 1) legend(x = 'topright', legend = model_spec, fill = color_palett_fill, border = color_palett_border, bty = 'n', cex = 1)
    }
  dev.off()
}

#alternate more space effeicent R^2 plot
R_sq_quant_plot <- function(list_jags_obj, color_palett_fill, color_palett_border, model_spec, include_models = seq_along(model_spec), file_name = 'R_sq_quant_plot_all.pdf', output_loc){
  number_of_measures = length(list_jags_obj)
  names_of_measures = names(list_jags_obj)
  names_of_measures = c('SPG', expression(lambda), expression(rho), expression('CV('~lambda~')'))
  number_of_models = length(include_models)
  y_max = number_of_measures
  x_max = 0.7
  x_offset = 0.4
  x_min = 0 - (x_max * x_offset)
  within_pred_offset = seq(0.1, 0.9, length = number_of_models)
  measure_y_pos = seq(0.5, y_max, 1)
  r2_qants = list()
  for(i in 1:number_of_measures) r2_qants[[i]] <- sapply(list_jags_obj[[i]], FUN = function(x) density = quantile(x$R_sq, probs = c(0.01, 0.025, 0.5, 0.975, 0.99)))
  #make the blank plot
  setwd(output_loc)
  pdf(file = file_name, height = 8, width = 10)
    par(mar = c(4, 0, 1.5, 1))
    plot(x = seq(x_min, x_max, length = 10), y = seq(0, y_max, length = 10), type = 'n', bty = 'n', xlab = expression(~~~~~~~~~~~~'                    '~R^2), ylab = '', main = '', tck = 0.015, cex.axis = 2, cex.lab = 2, 
	  ylim = c(0, y_max), xlim = c(x_min, x_max), yaxt = 'n', xaxt = 'n')
    axis(side = 1, at = seq(0, x_max, 0.1), labels = paste0(seq(0, x_max, 0.1)*100,'%'), tck = 0.015, cex.axis = 1.5)
    #make light grey section dividers
    altRow <- seq(y_max - 1, 1, -2)
    rect(xleft = x_min, xright = x_max*1.01, ybottom = altRow - 1, ytop = altRow, border = NA, col = grey(0.95)) 
    #put text labels
    text(x = x_min, y = y_max + 0.05, labels = 'Measure', pos = 4, cex = 1.7)
    text(x = -0.01, y = y_max + 0.05, labels = 'Model', pos = 2, cex = 1.7)
    
    #put the lines and names on the plot
    for(i in 1:number_of_measures){
      text(x = x_min, y = measure_y_pos[i], labels = names_of_measures[i], pos = 4, cex = 1.5)
      for(j in seq_along(include_models)){
	text(x = -0.01, y = i - within_pred_offset[j], labels = model_spec[include_models[j]], pos = 2, col = color_palett_border[include_models[j]])
	lines(x = c(r2_qants[[i]][c('1%', '99%'), include_models[j]]), y = c(i - within_pred_offset[j], i - within_pred_offset[j]), lend = 2, col = color_palett_fill[include_models[j]], lwd = 2)
	lines(x = c(r2_qants[[i]][c('2.5%', '97.5%'), include_models[j]]), y = c(i - within_pred_offset[j], i - within_pred_offset[j]), lend = 2, col = color_palett_border[include_models[j]], lwd = 2)
	points(x = r2_qants[[i]][c('50%'), include_models[j]], y = i - within_pred_offset[j], pch = 19, col = color_palett_border[include_models[j]], cex = 1)
      }
    }
  dev.off()
}      

#Function to make a quantile plot of all the measures and  
quant_plots <- function(list_jags_obj, color_palett_fill, color_palett_border, model_spec, file_name = 'quant_coef_plot_all.pdf', output_loc){
  number_of_measures = length(list_jags_obj)
  names_of_measures = names(list_jags_obj)
  number_of_models = length(model_spec)
  y_max = max(length(list_jags_obj[[1]][[1]]$quant_df$stand_name)) + 1
  x_offset = 0.25
  pred_names = rev(as.character(list_jags_obj[[1]][[1]]$quant_df$stand_name))
  within_pred_offset = seq(0.1, 0.9, length = number_of_models)
  
  setwd(output_loc)
  pdf(file = file_name,  paper = 'a4', height = 14, width = 12)
    plot(x = 1:10, y = 1:10, bty = 'n', type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
    legend(x = 'center', legend = model_spec, fill = color_palett_fill, border = color_palett_border, bty = 'n', cex = 2, title = 'model specification')
    for(i in 1:number_of_measures){
      x_max =  max(sapply(list_jags_obj[[i]], FUN = function(x) max(x$quant_df$upper_99, na.rm = TRUE)))
      x_min =  min(sapply(list_jags_obj[[i]], FUN = function(x) min(x$quant_df$lower_01, na.rm = TRUE)))
      plot_x_min = x_min - ((x_max - x_min) * x_offset)
      #set up blank plotting area 
      par(mar = c(4,0,1,0), oma = c(0,0,0,0))
      plot(x = seq(plot_x_min, x_max, length = 10), y = seq(0, y_max, length = 10), type = 'n', bty = 'n', xlab = 'coefficients', ylab = '', xaxt = 'n', yaxt = 'n')
      #put alternate grey bars on the plot
      altRow <- seq(y_max - 1, 1, -2)
      rect(xleft = plot_x_min, xright = x_max*1.05, ybottom = altRow - 1, ytop = altRow, border = NA, col = grey(0.95)) 
      text(x = plot_x_min, y = y_max + 0.5, labels = 'Predictor', pos = 4, cex = 1)
      text(x = 0, y = y_max + 0.7, labels = paste0('measure: ', names_of_measures[i]), cex = 1.2)
      lines(x = c(0, 0), y = c(0, y_max - 1), col = adjustcolor('red', alpha.f = 0.2))
      text(x = x_min, y = 0.5:length(pred_names), labels = pred_names, pos = 2, cex = 0.8)
      axis(side = 1, at = round(c(seq(x_min, 0, length = 5), seq(0, x_max, length = 5)), 3), tck = 0.02)
      #plot the lines
      for(k in 1:number_of_models){
	#pull the quantile dataframe out of the object and rename for convenience
	quant_df <- list_jags_obj[[i]][[k]]$quant_df
	used_names = as.character(quant_df$stand_name[!is.na(quant_df$lower_01)])
	plot_names = used_names
	#sort out the specicial case where in geo only model the auto_corr is put in the rho value by swapping rho for alpha and neighbours for reltives for phy models
	if(grepl('phy', model_spec[k])){
	  plot_names[plot_names == 'rho'] <- 'alpha'
	  plot_names[plot_names == 'neighbours'] <- 'realtives'
	}
	
	for(j in seq_along(used_names)){
	  lines(x = quant_df[which(quant_df$stand_name == used_names[j]), c('lower_01', 'upper_99')], y = c(which(pred_names == plot_names[j]), which(pred_names == plot_names[j])) - within_pred_offset[k], lend = 2, col = color_palett_fill[k])
	  lines(x = quant_df[which(quant_df$stand_name == used_names[j]), c('lower_025', 'upper_975')], y = c(which(pred_names == plot_names[j]), which(pred_names == plot_names[j])) - within_pred_offset[k], lend = 2, col = color_palett_border[k])
	  points(x = quant_df[which(quant_df$stand_name == used_names[j]), c('mean')], y = which(pred_names == plot_names[j]) - within_pred_offset[k], pch = 20, cex = 0.5, col = color_palett_border[k])
	}
	
      }
    }
  dev.off()
}


#function to plot the decay parameters of the different models
decay_plots <- function(list_jags_obj, color_palett_samp_lines, model_spec, names_of_models = c('main_int', 'phygeo-all_pops', 'phygeo-no_self', 'main', 'geo-all_pops', 'geo-no_self', 'phy-all_pops', 'phy-no_self'),
  names_of_measures = names(list_jags_obj), file_name = 'decay_models_plotted_all.pdf', output_loc){
  setwd(output_loc)
  pdf(file = file_name, height = 14, width = 14)  
    number_of_measures = length(list_jags_obj)
    number_of_models = length(model_spec)
    layout_mat = rbind(c(0, 1, 0, 1), c(0.007, 0.239, 0.6, 0.9), c(0.239, 0.487, 0.6, 0.9), #specifies the split-plot layout
      c(0.007, 0.239, 0.3, 0.6), c(0.239, 0.487, 0.3, 0.6), c(0.007, 0.239, 0, 0.3), c(0.239, 0.487, 0, 0.3), 
      c(0.52, 0.76, 0.6, 0.9), c(0.76, 1, 0.6, 0.9), c(0.52, 0.76, 0.3, 0.6), c(0.76, 1, 0.3, 0.6), 
      c(0.52, 0.76, 0, 0.3), c(0.76, 1, 0, 0.3)) 
    geo_x_seq = seq(0, 100, length = 100) 
    phy_x_seq = seq(0, 150, length = 100)
    geo_inds = which(!grepl('phy', model_spec))
    phy_inds = which(!grepl('geo', model_spec))
    
    for(i in 1:number_of_measures){
      split.screen(layout_mat)
      screen(1)
	par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
	plot(x = 0:1, y = 0:1, bty = 'n', type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
	rect(xleft = c(-0.1, 0.51), ybottom = c(-0.1, -0.1), xright = c(0.49, 1.1), ytop = c(1.1, 1.1), col = grey(0.95), border = grey(0.95))
	text(x = 0.5, y = 0.95, labels = names_of_measures[i], cex = 1.5)
      screen_count = 2
      for(k in geo_inds){
	w <- sapply(geo_x_seq, FUN = function(dist) exp(-list_jags_obj[[i]][[k]]$full_obj$BUGSoutput$sims.list$rho * dist)) #creates the matrix of weights for each sanple by applying every element of the dist to all samples of rho 
	screen(screen_count)
	  par(mar = c(4, 4, 2, 1), oma = c(0, 0 ,0, 0))
	  plot(x = geo_x_seq, y = seq(0, 1, length = length(geo_x_seq)), type = 'n', bty = 'n', xlab = 'Geographic dist. (km)', ylab = 'weight', main = names_of_models[k], xlim = c(0, max(geo_x_seq)), ylim = c(0, 1), tck = 0.02)
	  for(j in 1:dim(w)[1]) lines(x = geo_x_seq, y = w[j, ], lwd = 0.5, col = color_palett_samp_lines[k])#sample lines to show uncertianity 
	  lines(x = geo_x_seq, y = exp(-list_jags_obj[[i]][[k]]$full_obj$BUGSoutput$mean$rho * geo_x_seq), lwd = 2, col = grey(0.3))#mean line
	screen_count = screen_count + 1
      }
      
      for(k in phy_inds){
	if(grepl('phy', model_spec[k])){#handels the case where phy uses rho instead of alpha as the decay parameter
	  w <- sapply(phy_x_seq, FUN = function(dist) exp(-list_jags_obj[[i]][[k]]$full_obj$BUGSoutput$sims.list$rho * dist))
	}else{
	  w <- sapply(phy_x_seq, FUN = function(dist) exp(-list_jags_obj[[i]][[k]]$full_obj$BUGSoutput$sims.list$alpha * dist)) #creates the matrix of weights for each sanple by applying every element of the dist to all samples of rho 
	}
	screen(screen_count)
	  par(mar = c(4, 4, 2, 1), oma = c(0, 0 ,0, 0))
	  plot(x = phy_x_seq, y = seq(0, 1, length = length(phy_x_seq)), type = 'n', bty = 'n', xlab = 'Phylogenetic dist. (million yrs)', ylab = 'weight', main = names_of_models[k], xlim = c(0, max(phy_x_seq)), ylim = c(0, 1), tck = 0.02)
	  for(j in 1:dim(w)[1]) lines(x = phy_x_seq, y = w[j, ], lwd = 0.5, col = color_palett_samp_lines[k])#sample lines to show uncertianity 
	  #mean line
	  if(grepl('phy', model_spec[k])){#handels the case where phy uses rho instead of alpha as the decay parameter
	    lines(x = phy_x_seq, y = exp(-list_jags_obj[[i]][[k]]$full_obj$BUGSoutput$mean$rho * phy_x_seq), lwd = 2, col = grey(0.3))
	}else{
	  lines(x = phy_x_seq, y = exp(-list_jags_obj[[i]][[k]]$full_obj$BUGSoutput$mean$alpha * phy_x_seq), lwd = 2, col = grey(0.3))
	}
	screen_count = screen_count + 1
      }
      close.screen(all.screens = TRUE)
    }
  dev.off()
}

#Function to plot the decay plots for a specified model for each measure. Indented to be used for publication 
decay_plots_single_model <- function(list_jags_obj, color, model_spec, plotted_model, names_of_measures = names(list_jags_obj), file_name = 'decay_models_single_model.pdf', output_loc){
  number_of_measures = length(list_jags_obj)
  letter_labels = paste0(letters[1:8], ')') 
  metric_order = c(2, 4, 1, 3) #sets the plot order so it goes lambda, CV(lambda), SGG, damping
  layout_mat = rbind(c(0, 1, 0, 1),#big overall plot area to print some rectangles and text to demark different areas 
    c(0.007, 0.247, 0.48, 0.96), c(0.247, 0.487, 0.48, 0.96), c(0.007, 0.247, 0, 0.48), c(0.247, 0.487, 0, 0.48), #geo plot positions
    c(0.52, 0.76, 0.48, 0.96), c(0.76, 1, 0.48, 0.96), c(0.52, 0.76, 0, 0.48), c(0.76, 1, 0, 0.48))#phy plot positions
  
  geo_x_seq = seq(0, 100, length = 100) 
  phy_x_seq = seq(0, 150, length = 100)
  geo_inds = which(!grepl('phy', model_spec[plotted_model]))
  phy_inds = which(!grepl('geo', model_spec[plotted_model]))
  setwd(output_loc)
  pdf(file = file_name, height = 7, width = 14) 
    #make a base to plot put the rest on top of 
    split.screen(layout_mat)
    screen(1)
      par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
      plot(x = 0:1, y = 0:1, type = 'n', bty = 'n', axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), main = '', ylab = '', xlab = '')
      rect(xleft = c(-0.1, 0.51), ybottom = c(-0.1, -0.1), xright = c(0.49, 1.1), ytop = c(1, 1), col = grey(0.95), border = grey(0.95))
      #text(x = c(0.25, 0.8), y = c(1, 1), labels = c('Geographic distance', 'Phylogenetic distance'), cex = 1.5) #TAKEN OUT AS IS A BIT REDUNDANT
    if(grepl('geo', model_spec[plotted_model])){ #test that there are some geo terms in the model being requested
      screen_count = 2
      for(i in metric_order){
	w <- sapply(geo_x_seq, FUN = function(dist) exp(-list_jags_obj[[i]][[plotted_model]]$full_obj$BUGSoutput$sims.list$rho * dist)) #creates the matrix of weights for each sanple by applying every element of the dist to all samples of rho 
	screen(screen_count)
	  par(mar = c(4, 4, 2, 0.5), oma = c(0, 0 ,0, 0))
	  plot(x = geo_x_seq, y = seq(0, 1, length = length(geo_x_seq)), type = 'n', bty = 'n', xlab = 'Geographic dist. (km)', ylab = 'Weight', main = '', xlim = c(0, max(geo_x_seq)), ylim = c(0, 1), tck = 0.02, cex.lab = 1.2)
	  for(j in 1:dim(w)[1]) lines(x = geo_x_seq, y = w[j, ], lwd = 0.5, col = color)#sample lines to show uncertianity 
	  lines(x = geo_x_seq, y = exp(-list_jags_obj[[i]][[plotted_model]]$full_obj$BUGSoutput$mean$rho * geo_x_seq), lwd = 2, col = grey(0.3))#mean line
	  text(x = 50, y = 0.95, labels = names_of_measures[i], cex = 1.5)
	  mtext(text = letter_labels[screen_count - 1] , adj = 0, side = 3, line = 0, cex = 1.2)
	screen_count = screen_count + 1
      }
    }
    #now plot the phy modesl
    if(grepl('phy', model_spec[plotted_model])){ #test that there are some geo terms in the model being requested
      for(i in metric_order){
	if(!grepl('geo', model_spec[plotted_model])){#handels the case where phy uses rho instead of alpha as the decay parameter
	    w <- sapply(phy_x_seq, FUN = function(dist) exp(-list_jags_obj[[i]][[plotted_model]]$full_obj$BUGSoutput$sims.list$rho * dist))
	}else{
	  w <- sapply(phy_x_seq, FUN = function(dist) exp(-list_jags_obj[[i]][[plotted_model]]$full_obj$BUGSoutput$sims.list$alpha * dist)) #creates the matrix of weights for each sanple by applying every element of the dist to all samples of rho 
	}
	screen(screen_count)
	  par(mar = c(4, 4, 2, 0.5), oma = c(0, 0 ,0, 0))
	  plot(x = phy_x_seq, y = seq(0, 1, length = length(phy_x_seq)), type = 'n', bty = 'n', xlab = 'Phylogenetic dist. (million yrs)', ylab = 'Weight', main = '', xlim = c(0, max(phy_x_seq)), ylim = c(0, 1), tck = 0.02, cex.lab = 1.2)
	  for(j in 1:dim(w)[1]) lines(x = phy_x_seq, y = w[j, ], lwd = 0.5, col = color)#sample lines to show uncertianity 
	  #mean line
	  if(!grepl('geo', model_spec[plotted_model])){#handels the case where phy uses rho instead of alpha as the decay parameter
	    lines(x = phy_x_seq, y = exp(-list_jags_obj[[i]][[plotted_model]]$full_obj$BUGSoutput$mean$rho * phy_x_seq), lwd = 2, col = grey(0.3))
	  }else{
	    lines(x = phy_x_seq, y = exp(-list_jags_obj[[i]][[plotted_model]]$full_obj$BUGSoutput$mean$alpha * phy_x_seq), lwd = 2, col = grey(0.3))
	  }
	  text(x = 75, y = 0.95, labels = names_of_measures[i], cex = 1.5)
	  mtext(text = letter_labels[screen_count - 1] , adj = 0, side = 3, line = 0, cex = 1.2)
	screen_count = screen_count + 1
      }
    }
    close.screen(all.screens = TRUE)
  dev.off()
 }
  
  
  
## PLOT FOR PAPER, TEMPATURE SPACE AND PC1 SOCRES MAPPED WITH STUDY LOCATIONS
#biplot function customised using ggplot2
library(ggplot2)
library(grid)

PCbiplot <- function(PC, x = "PC1", y = "PC2", colors = c('black', 'black', 'red', 'red')) {
    # PC being a prcomp object
    data <- data.frame(obsnames=row.names(PC$x), PC$x)
    plot <- ggplot(data, aes_string(x = x, y = y)) + geom_point(shape = 19, size = 3, aes(label = obsnames), color = colors[1])
    plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2, color=colors[2])
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
    datapc <- transform(datapc,
            v1 = .7 * mult * (get(x)),
            v2 = .7 * mult * (get(y))
            )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color=colors[3])
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color=colors[4])
    plot
}

#PCbiplot(temp_pca, color = c('black', 'black', 'blue', 'blue'))

temp_space_PCA_plotter <- function(env_pred, temp_pca, precip_pca, data_source, darkness = 0, filename = 'temp_season_PCA_mapped_550.pdf', output_loc){  
  #make a plot of the map and the mean_temp and temp_seasonality environmetal space
  #rescale a few things to colour rnages and numbers are easier to work with between plots and between variables
  pc1_scaled <- ((1 - darkness)*(data_source$pc1_temp - min(data_source$pc1_temp)) / (max(data_source$pc1_temp) - min(data_source$pc1_temp))) + darkness
  temp_range <- env_pred$bio_1@data@max - env_pred$bio_1@data@min
  temp_min <- env_pred$bio_1@data@min
  print('rescaling some layers this can take some time')
  temp_scaled <- (env_pred$bio_1 - temp_min) / temp_range 
  print('Finished 1 of 3')
  var_range <- env_pred$bio_4@data@max - env_pred$bio_4@data@min
  var_min <- env_pred$bio_4@data@min
  var_scaled <- (env_pred$bio_4 - var_min) / var_range
  obs_temp_scaled <- (data_source$bio_1 - temp_min) / temp_range
  obs_var_scaled <- (data_source$bio_4 - var_min) / var_range
  print('Finished 2 of 3')
  layer_of_0 <- env_pred$bio_1*0
  temp_var_stack <- stack(temp_scaled, var_scaled, layer_of_0)
  print('Finsihed 3 of 3')

  #plot the study points and PC1 scores and bi-plot of the PCA
  layout_mat <- rbind(c(0, 0.85, 0.45, 1), c(0.85, 1, 0.45, 1), c(0, 0.333, 0, 0.45), c(0.333, 0.66, 0, 0.45), c(0.666, 1, 0, 0.45))
  
  setwd(output_loc)
  pdf(file = filename, width = 22, height = 15)
  #try and work out a sensible color map for two varibles mean_temp (bio_1) and temp_varibility (bio_4)
  #make temp represented by red, higher temps = redder.
  #make varibility represented by blue with more varible places = bluer
    split.screen(layout_mat)
    screen(1)#first screen is the map
      plotRGB(temp_var_stack, r = 1, g = 3, b = 2, scale = 1)
      points(data_source[, 'lon_dd'], data_source[,'lat_dd'], pch = 20, col = grey(c(1 - pc1_scaled)))
      text(x = temp_var_stack$bio_1.1@extent@xmin+10, y = temp_var_stack$bio_1.1@extent@ymax - 3, labels = 'a)', cex = 2)
    
    screen(2)#plot the PC1 legend
      dummy_pc1 <- seq(1, darkness, length = 500)
      par(mar = c(1, 10, 2, 1))
      image(t(as.matrix(dummy_pc1)), axes = TRUE, bty = 'n', xaxt = 'n', yaxt = 'n', ylab = '', xlab = '', col = grey(1-dummy_pc1))
      pc1_temp_range <- max(data_source$pc1_temp) - min(data_source$pc1_temp)
      axis(side = 2, at = seq(0, 1, length = 5), labels = round(seq(min(data_source$pc1_temp), max(data_source$pc1_temp), length = 5), 1), lwd = 0, tck = 0, , line = 1.5)
      mtext(c('colder\nmore seasonal', 'hotter\nless seasonal'), side = 2, line = 7, at = c(0, 1), adj = 0.5, padj = 0.5, las = 2, cex = 1.2)
      mtext('PC1 score', side = 2, line = 4, adj = 0.5, padj = 0.5, cex = 1.3)
      #give the density of points along the axis
      axis(side = 2, at = 1 - pc1_scaled, labels = rep('', length(pc1_scaled)), tck = -0.1)
      
    screen(3)# plot the environmetnal space
      #make 2D legend for colours
      dummy_temp <- raster(matrix(rep(seq(0, 1, length = 500), 500), ncol = 500, byrow = TRUE))
      dummy_var <- raster(matrix(rep(seq(1, 0, length = 500), 500), ncol = 500, byrow = FALSE))
      dummy_0 <- raster(matrix(rep(rep(0, 500), 500), ncol = 500)) 
      legend_stack <- stack(dummy_temp, dummy_var, dummy_0)
      par(mar = c(5, 6, 3, 1))
      plot(seq(0, 1, 0.2), seq(0, 1, 0.2), type = 'n', bty = 'n', axes = FALSE, xlab = 'Mean annual temperature (C)', ylab = 'Temperature seasonality (sd of temperature)', cex.lab = 1.5, cex.axis = 1.2,
	main = 'Temperature environmental space ', cex.main = 1.5)
      mtext('b)', side = 3, adj = -0.2, cex = 2)
      
      plotRGB(legend_stack, r = 1, g = 3, b = 2, scale = 1, add = TRUE)#, axes = TRUE, xlab = 'mean annual tempature', ylab = 'tempature seasonality')
      axis(side = 1, at = seq(0, 1, length = 5), labels = round(seq(temp_min, env_pred$bio_1@data@max, length = 5)/10, 1), line = 0.7, tck = 0.02) 
      axis(side = 2, at = seq(0, 1, length = 5), labels = round(seq(var_min, env_pred$bio_4@data@max, length = 5)/100, 1), line = 0.7, tck = 0.02) 
      #plot the points of where my study points lie in that environmental space
      #color code gray by PCA1 score 
      points(obs_temp_scaled, obs_var_scaled, pch = 20, col = grey(c(1 - pc1_scaled)))

    screen(4)
      par(mar = c(4, 3, 4, 0.5))
      biplot(temp_pca, col = c('black', 'blue'), cex = c(0.4, 1), main = 'Temperature PCA', cex.main = 1.5, cex.lab = 1.5) #is a bit hard to interprit but seems like variability (bio_4 and bio_7) and mean temp (bio_1)strech out along the PC1
      mtext('c)', side = 3, adj = -0.17, cex = 2)
    
    screen(5)
      par(mar = c(4, 3, 4, 1))
      biplot(precip_pca, col = c('black', 'blue'), cex = c(0.4, 1), ylim = c(-0.15, 0.2), main = 'Precipitation PCA', cex.main = 1.5, cex.lab = 1.5) #is a bit hard to interprit but seems like variability (bio_4 and bio_7) and mean temp (bio_1)strech out along the PC1
      mtext('d)', side = 3, adj = -0.17, cex = 2)
    
    close.screen(all.screens = TRUE)    # exit split-screen mode
  dev.off()
}

