#set of helper functions in processing model output from jags

#model evaluation done in a couple of ways
jags_model_eval = function(jags_obj, output_file_name = 'default_name_jags_model_output.pdf', obs_Y){
  pdf(file = output_file_name, paper = 'a4')
    plot(jags_obj$BUGSoutput$sims.list$fit, jags_obj$BUGSoutput$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, 
      main = paste0('R_hats min - max: ', paste(summary(jags_obj$BUGSoutput$summary[,8]), collapse = ', '), '\nBp-value: ', sum(jags_obj$BUGSoutput$sims.list$bpvalue)/length(jags_obj$BUGSoutput$sims.list$bpvalue)))
    abline(0, 1, col = "black")

    mus = matrix(nrow = dim(jags_obj$BUGSoutput$sims.list$residual)[1], ncol = dim(jags_obj$BUGSoutput$sims.list$residual)[2])
    for(i in seq_along(obs_Y)) mus[, i] = obs_Y[i] - jags_obj$BUGSoutput$sims.list$residual[,i]
    mean_mu <- apply(mus, MARGIN = 2, FUN = mean)

    #residual plot 
    plot(mean_mu, jags_obj$BUGSoutput$mean$residual, main = "", xlab = 
      "predicted value", ylab = "residual", frame.plot = FALSE)
    abline(0, 0, col = "black")

    #do a fitted versus predicted plot
    mu_ci <- t(apply(mus, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))#get uncertianty in the estimated value using credible intervals
    x_couples <- cbind(obs_Y, obs_Y)
    plot(obs_Y, mean_mu, main = 'Observed vs Predicted (95%CI)', ylim = c(min(mu_ci), max(mu_ci)), bty = 'n', 
      type = 'n', xlab = 'observed Y', ylab = 'predicted Y', tck = 0.015, cex.lab = 2.5, cex.axis = 2, cex.main = 2.5)
    for(j in seq_along(mu_ci[,1])) lines(x = x_couples[j, ], y = mu_ci[j, ], lwd = 1, col = grey(0.6))
    points(obs_Y, mean_mu, cex = 1, pch = 19)
    abline(0, 1, col = "black")

    #look at the chains
    print('printing chains, this takes a while')
    jags_obj_mc <- as.mcmc(jags_obj)
    plot(jags_obj_mc, ask = FALSE)
    print('printing autocorrelation plots, these are faster but still take a while')
    autocorr.plot(jags_obj_mc, ask = FALSE, lag.max = 10)
  dev.off()
  dev.off()
  dev.off()
}

#set up the distance matricies
dist_matrix_setup = function(phylo_tree, long_lat, species_name){
  #caclulate distance to every other population, both geographic and phylogenetic. 
  geo_dist <- matrix(nrow = dim(long_lat)[1], ncol = dim(long_lat)[1] - 1) #minus 1 becasue distance of each pop to itself is 0
  phy_dist <- matrix(nrow = dim(long_lat)[1], ncol = dim(long_lat)[1] - 1) #minus 1 becasue distance of each pop to itself is 0_dist <- matrix(nrow = length(sp_data$vr_elast_PC1), ncol = length(sp_data$vr_elast_PC1) - 1) #minus 1 becasue distance of each pop to itself is 0
  TLCA_all <- cophenetic(phylo_tree) #Time to last comon ancestor for all species, need to pick out each pair of species in the actual data set 
  
  print('Finding geogrpahic and phylogentic distances')
  #first pass to make the distance and neighbour value matricies
  for(y in seq_along(long_lat[,1])){ 
    count = 1
    for(x in seq_along(long_lat[,1])[-y]){
      geo_dist[y, count] <- Ellipsoidal.Distance(olat = long_lat[y, 2], olon = long_lat[y, 1], tlat = long_lat[x, 2], tlon = long_lat[x, 1])$dist #geographic distance
      phy_dist[y, count] <- as.numeric(TLCA_all[species_name[y], species_name[x]])#gets phylogenetic distance 
      count = count + 1
    }
  } 
  
  print('removing the 0 distances from both distance matrices independently')
  #second pass to take out all the zero distances 
  geo_dist_noz <- matrix(data = 0, nrow = dim(long_lat)[1], ncol = dim(long_lat)[1]) #force last coloumn to be 0 so that can find the index of the last non-zero input
  phy_dist_noz <- matrix(data = 0, nrow = dim(long_lat)[1], ncol = dim(long_lat)[1]) #force last coloumn to be 0 so that can find the index of the last non-zero input0
 
  for(y in seq_along(long_lat[,1])){ 
    count_geo = 1
    count_phy = 1
    for(x in 1:(dim(long_lat)[1] - 1)){
      if(geo_dist[y, x] > 0.001){
	geo_dist_noz[y, count_geo] <- geo_dist[y, x]
	count_geo = count_geo + 1
      }
      if(phy_dist[y, x] > 0.01){
	phy_dist_noz[y, count_phy] <- phy_dist[y, x]
	count_phy = count_phy + 1
      } 
    }
  } 
  #find the maximum index that each prediction goes to, that is how many other popualtions is each populations prediction based on 
  geo_max_ind <- apply(geo_dist_noz, MARGIN = 1, FUN = function(x) which(x == 0)[1] - 1)
  phy_max_ind <- apply(phy_dist_noz, MARGIN = 1, FUN = function(x) which(x == 0)[1] - 1)
  
  return(list(geo_dist = geo_dist, phy_dist = phy_dist, geo_dist_noz = geo_dist_noz, phy_dist_noz = phy_dist_noz, 
    geo_max_ind = geo_max_ind, phy_max_ind = phy_max_ind))
}

get_neigh_values = function(geo_dist, phy_dist, Y){
  Y_neigh <- matrix(nrow = length(Y), ncol = length(Y) - 1) #says what the value of Y elast is for each neighbour or relative  
  print('Finding Y values for neighbour')
    for(y in seq_along(Y)){ 
    count = 1
    for(x in seq_along(Y)[-y]){
      Y_neigh[y, count] <- Y[x]
      count = count + 1
    }
  }
  
  Y_neigh_phy <- matrix(data = 0, nrow = length(Y), ncol = length(Y) - 1) #says what the value of Y elast is for each neighbour or relative  
  Y_neigh_geo <- matrix(data = 0, nrow = length(Y), ncol = length(Y) - 1) #says what the value of Y elast is for each neighbour or relative  
  print('Finding Y values for neighbour with zero distances ignored')
  for(y in seq_along(Y)){ 
    count_geo = 1
    count_phy = 1
    for(x in 1:(length(Y) - 1)){
      if(geo_dist[y, x] > 0.001){
	Y_neigh_geo[y, count_geo] <- Y_neigh[y, x]
	count_geo = count_geo + 1
      }
      if(phy_dist[y, x] > 0.01){
	Y_neigh_phy[y, count_phy] <- Y_neigh[y, x]
	count_phy = count_phy + 1
      } 
    }
  }  
  return(list(Y_neigh = Y_neigh, Y_neigh_geo = Y_neigh_geo, Y_neigh_phy = Y_neigh_phy))
}
