#Script to find how several indicies of population dynamics can be preidcited
library(GEOmap)
library(gtools)
library(ape)
library(caper)
library(phytools)
library(sp)
library(maptools)
library(dismo)
library(raster)
library(rgdal)
library(rgl)
data(wrld_simpl)
library(lattice)
library(coda)
library(R2WinBUGS)
#install.packages(pkgs = 'R2WinBUGS', repos = 'http://cran.us.r-project.org')
library(R2jags)
#install.packages(pkgs = 'R2jags', repos = 'http://cran.us.r-project.org')
library(MASS)
#set the working directory to hold the various working objects on disk

# This is an attempt to automaticaly find the file location of the extracted folder, this requires root access, it works on Linux (and so probably on mac as well), 
# not sure how it will go on window
print('this is an automated attemp to find and set the model folder location')
print('It requires root access. The file path can also be set manually at line')
print('32 in file pop_metric_prediction.R')
file_loc = paste0(strsplit(system('sudo updatedb;locate combined_bc_ai_all_responses.Rdata', intern = TRUE), '/')[[1]], '/')
model_file_loc = ''
for(i  in 1:(length(file_loc) - 1)) model_file_loc = paste0(model_file_loc, file_loc[i])
#MANUAL SETTING
#model_file_loc = '/your/path/to/folder/here'
print('this is an automated attemp to find and set the data folder location')
print('It requires root access. The file path can also be set manually at line')
print('107 in file data_extract_clean_up_first_pass.R')
file_loc = paste0(strsplit(system('sudo updatedb;locate PMDbl260314A1.nex', intern = TRUE), '/')[[1]], '/')
data_file_loc = ''
for(i  in 1:(length(file_loc) - 1)) data_file_loc = paste0(data_file_loc, file_loc[i])
#
#first get the data I need
phyloIN <- data_file_loc
compadreIN <- model_file_loc
analysis_ind <- model_file_loc

setwd(analysis_ind)
source('non_elast_predict_models.R')
setwd(analysis_ind)
source('setup_processing_helper_functions.R')
model_names <- get_model_names()
#get the data_sets I need
#ge the tree and 
setwd(phyloIN)
Tree <- read.nexus("PMDbl260314A1.nex")

setwd(compadreIN)
obj_name <- load('combined_bc_ai_all_responses.Rdata')

dist_mats <- dist_matrix_setup(phylo_tree = Tree, long_lat = cr_pops[, c('lon_dd', 'lat_dd')], species_name = cr_pops$SpeciesAuthor)
neigh_values_lambda <- get_neigh_values(geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y = cr_pops$lambda)
neigh_values_dr <- get_neigh_values(geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y = cr_pops$damping_ratio)
neigh_values_PC1 <- get_neigh_values(dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y = cr_pops$elast_PC1)

#create the reduced data for cv_lambda
cv_pops <- cr_pops[which(!is.na(cr_pops$cv_lambda)),]
dist_mats_cv <- dist_matrix_setup(phylo_tree = Tree, long_lat = cv_pops[, c('lon_dd', 'lat_dd')], species_name = cv_pops$SpeciesAuthor)
neigh_values_cv <- get_neigh_values(geo_dist = dist_mats_cv$geo_dist, phy_dist = dist_mats_cv$phy_dist, Y = cv_pops$cv_lambda)

################################################################################################################################################################################################################
#### FIT MODELS WITH ENVIRONMENTAL PREDICTORS FOR LAMBDA #######################################################################################################################################################
#full interactions model, all pops ### ### ####### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = log(neigh_values_lambda$Y_neigh),
  Y_neigh_phy = log(neigh_values_lambda$Y_neigh), num_data = length(cr_pops$lambda), geo_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), GF = cr_pops$herb_fact, num_GF = max(cr_pops$herb_fact), mat_dim = cr_pops$MatrixDimension, 
  life_expect = cr_pops$life_expect, PC_temp = cr_pops$pc1_temp, AI = log(cr_pops$AI), p_sea = cr_pops$bio_15, min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(max(cr_pops$herb_fact), -20, 20), b = runif(9, -20, 20), B = cbind(runif(4, -20, 20), runif(4, -20, 20)), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", 'B', "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_full_allpops <- jags(jags.data, inital, monitored, "interactions_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_full_allpops, file = 'lambda_full_allpops.Rdata')
jags_model_eval(jags_obj = lambda_full_allpops, output_file_name = 'lambda_full_allpops_eval_output.pdf', obs_Y = log(cr_pops$lambda))

#ALL CONVERGED

#main effects only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = log(neigh_values_lambda$Y_neigh),
  Y_neigh_phy = log(neigh_values_lambda$Y_neigh), num_data = length(cr_pops$lambda), geo_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), GF = cr_pops$herb_fact, num_GF = max(cr_pops$herb_fact), 
  mat_dim = cr_pops$MatrixDimension, life_expect = cr_pops$life_expect, PC_temp = cr_pops$pc1_temp, AI = log(cr_pops$AI), p_sea = cr_pops$bio_15, min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(max(cr_pops$herb_fact), -20, 20), b = runif(7, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_noint_allpops <- jags(jags.data, inital, monitored, "noint_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_noint_allpops, file = 'lambda_noint_allpops.Rdata')
jags_model_eval(jags_obj = lambda_noint_allpops, output_file_name = 'lambda_noint_allpops_eval_output.pdf', obs_Y = log(cr_pops$lambda))

#ALL CONVERGED

#### NO FIXED EFFECTS LAMBDA
#no fixed effects, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = log(neigh_values_lambda$Y_neigh),
  Y_neigh_phy = log(neigh_values_lambda$Y_neigh), num_data = length(cr_pops$lambda), geo_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), min_rho = 0.001, max_rho = 1, min_alpha = 0.001, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.001, 1), rho = runif(1, 0.001, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_nofixed_allpops <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_nofixed_allpops, file = 'lambda_nofixed_allpops.Rdata')
jags_model_eval(jags_obj = lambda_nofixed_allpops, output_file_name = 'lambda_nofixed_allpops_eval_output.pdf', obs_Y = log(cr_pops$lambda))

#ALL CONVERGED

#no fixed effects, noself ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), geo_dist = dist_mats$geo_dist_noz, phy_dist = dist_mats$phy_dist_noz, Y_neigh_geo = log(neigh_values_lambda$Y_neigh_geo),
  Y_neigh_phy = log(neigh_values_lambda$Y_neigh_phy), num_data = length(cr_pops$lambda), geo_max_ind = dist_mats$geo_max_ind, phy_max_ind = dist_mats$phy_max_ind, 
  min_rho = 0.01, max_rho = 1, min_alpha = 0.03, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.03, 1), rho = runif(1, 0.01, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 400000
nThin <- 500
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_nofixed_noself <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

#lambda_nofixed_noself <- update(lambda_nofixed_noself, n.iter = 100000, n.thin = 500)

save(lambda_nofixed_noself, file = 'lambda_nofixed_noself.Rdata')
jags_model_eval(jags_obj = lambda_nofixed_noself, output_file_name = 'lambda_nofixed_noself_eval_output.pdf', obs_Y = log(cr_pops$lambda))

##ALL CONVERGED

#no fixed effects geo only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), ac_dist = dist_mats$geo_dist, Y_neigh = log(neigh_values_lambda$Y_neigh),
  num_data = length(cr_pops$lambda), ac_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_onlygeo_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_onlygeo_allpops, file = 'lambda_onlygeo_allpops.Rdata')
jags_model_eval(jags_obj = lambda_onlygeo_allpops, output_file_name = 'lambda_onlygeo_allpops_eval_output.pdf', obs_Y = log(cr_pops$lambda))

#ALL CONVERGED

#no fixed effects phy only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), ac_dist = dist_mats$phy_dist, Y_neigh = log(neigh_values_lambda$Y_neigh),
  num_data = length(cr_pops$lambda), ac_max_ind = rep(dim(neigh_values_lambda$Y_neigh)[2], dim(neigh_values_lambda$Y_neigh)[1]), min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_onlyphy_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_onlyphy_allpops, file = 'lambda_onlyphy_allpops.Rdata')
jags_model_eval(jags_obj = lambda_onlyphy_allpops, output_file_name = 'lambda_onlyphy_allpops_eval_output.pdf', obs_Y = log(cr_pops$lambda))

#ALL CONVERGED

#no fixed effects geo only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), ac_dist = dist_mats$geo_dist_noz, Y_neigh = log(neigh_values_lambda$Y_neigh_geo),
  num_data = length(cr_pops$lambda), ac_max_ind = dist_mats$geo_max_ind, min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_onlygeo_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_onlygeo_noself, file = 'lambda_onlygeo_noself.Rdata')
jags_model_eval(jags_obj = lambda_onlygeo_noself, output_file_name = 'lambda_onlygeo_noself_eval_output.pdf', obs_Y = log(cr_pops$lambda))

#ALL CONVERGED

#no fixed effects phy only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$lambda), ac_dist = dist_mats$phy_dist_noz, Y_neigh = log(neigh_values_lambda$Y_neigh_phy),
  num_data = length(cr_pops$lambda), ac_max_ind = dist_mats$phy_max_ind, min_rho = 0.01, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0.01, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
lambda_onlyphy_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(lambda_onlyphy_noself, file = 'lambda_onlyphy_noself.Rdata')
jags_model_eval(jags_obj = lambda_onlyphy_noself, output_file_name = 'lambda_onlyphy_noself_eval_output.pdf', obs_Y = log(cr_pops$lambda))

##ALL CHAINS CONVERGED

#######################################################################################################################################################################################################################################################
### FIT MODELS WITH ENVIRONMENTAL PREDICTORS FOR damping ratio  #######################################################################################################################################################################################
#full interactions model, all pops ### ### ####### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = log(neigh_values_dr$Y_neigh),
  Y_neigh_phy = log(neigh_values_dr$Y_neigh), num_data = length(cr_pops$damping_ratio), geo_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), GF = cr_pops$herb_fact, num_GF = max(cr_pops$herb_fact), mat_dim = cr_pops$MatrixDimension, 
  life_expect = cr_pops$life_expect, PC_temp = cr_pops$pc1_temp, AI = log(cr_pops$AI), p_sea = cr_pops$bio_15, min_alpha = 0.001, max_alpha = 1, min_rho = 0.0001, max_rho =1)

#intial values
inital <- function() list(a = runif(max(cr_pops$herb_fact), 1, 1), b = runif(13, 0, 1), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.3, 0.3), rho = runif(1, 0.005, 0.005)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_full_allpops <- jags(jags.data, inital, monitored, "redu_int_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_full_allpops, file = 'damp_ratio_full_allpops.Rdata')
jags_model_eval(jags_obj = damp_ratio_full_allpops, output_file_name = 'damp_ratio_full_allpops_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED

#main effects only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = log(neigh_values_dr$Y_neigh),
  Y_neigh_phy = log(neigh_values_dr$Y_neigh), num_data = length(cr_pops$damping_ratio), geo_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), GF = cr_pops$herb_fact, num_GF = max(cr_pops$herb_fact), 
  mat_dim = cr_pops$MatrixDimension, life_expect = cr_pops$life_expect, PC_temp = cr_pops$pc1_temp, AI = log(cr_pops$AI), p_sea = cr_pops$bio_15, min_alpha = 0.01, 
  max_alpha = 1, min_rho = 0.01, max_rho = 1)

#intial values
inital <- function() list(a = runif(max(cr_pops$herb_fact), 0.5, 0.5), b = runif(7, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.15, 0.15), rho = runif(1, 0.15, 0.15)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_noint_allpops <- jags(jags.data, inital, monitored, "noint_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_noint_allpops, file = 'damp_ratio_noint_allpops.Rdata')
jags_model_eval(jags_obj = damp_ratio_noint_allpops, output_file_name = 'damp_ratio_noint_allpops_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED 

#### NO FIXED EFFECTS damp_ratio
#no fixed effects, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = log(neigh_values_dr$Y_neigh),
  Y_neigh_phy = log(neigh_values_dr$Y_neigh), num_data = length(cr_pops$damping_ratio), geo_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), min_alpha = 0, max_alpha = 1, min_rho = 0, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_nofixed_allpops <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_nofixed_allpops, file = 'damp_ratio_nofixed_allpops.Rdata')
jags_model_eval(jags_obj = damp_ratio_nofixed_allpops, output_file_name = 'damp_ratio_nofixed_allpops_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED

#no fixed effects geo only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), ac_dist = dist_mats$geo_dist, Y_neigh = log(neigh_values_dr$Y_neigh),
  num_data = length(cr_pops$damping_ratio), ac_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), min_alpha = 0, 
  max_alpha = 1, min_rho = 0, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_onlygeo_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_onlygeo_allpops, file = 'damp_ratio_onlygeo_allpops.Rdata')
jags_model_eval(jags_obj = damp_ratio_onlygeo_allpops, output_file_name = 'damp_ratio_onlygeo_allpops_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED

#no fixed effects, noself ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), geo_dist = dist_mats$geo_dist_noz, phy_dist = dist_mats$phy_dist_noz, Y_neigh_geo = log(neigh_values_dr$Y_neigh_geo),
  Y_neigh_phy = log(neigh_values_dr$Y_neigh_phy), num_data = length(cr_pops$damping_ratio), geo_max_ind = dist_mats$geo_max_ind, phy_max_ind = dist_mats$phy_max_ind, 
  min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.001, 1), rho = runif(1, 0.001, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 450000
nThin <- 800
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damping_ratio_nofixed_noself <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damping_ratio_nofixed_noself, file = 'damping_ratio_nofixed_noself.Rdata')
jags_model_eval(jags_obj = damping_ratio_nofixed_noself, output_file_name = 'damping_ratio_nofixed_noself_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#CONVERGED BUT AUTO_CORRELATION PROBLEMS IN SOME OF THE CHAINS, WILL HAVE TO SIMULATE LONGER TO GET ENOUGH SAMPLES, need to thin at every 800 
#WORKED AT 800 THINING

#no fixed effects phy only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), ac_dist = dist_mats$phy_dist, Y_neigh = log(neigh_values_dr$Y_neigh),
  num_data = length(cr_pops$damping_ratio), ac_max_ind = rep(dim(neigh_values_dr$Y_neigh)[2], dim(neigh_values_dr$Y_neigh)[1]), min_rho = 0, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_onlyphy_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_onlyphy_allpops, file = 'damp_ratio_onlyphy_allpops.Rdata')
jags_model_eval(jags_obj = damp_ratio_onlyphy_allpops, output_file_name = 'damp_ratio_onlyphy_allpops_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED

#no fixed effects geo only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), ac_dist = dist_mats$geo_dist_noz, Y_neigh = log(neigh_values_dr$Y_neigh_geo),
  num_data = length(cr_pops$damping_ratio), ac_max_ind = dist_mats$geo_max_ind, min_rho = 0, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_onlygeo_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_onlygeo_noself, file = 'damp_ratio_onlygeo_noself.Rdata')
jags_model_eval(jags_obj = damp_ratio_onlygeo_noself, output_file_name = 'damp_ratio_onlygeo_noself_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED

#no fixed effects phy only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cr_pops$damping_ratio), ac_dist = dist_mats$phy_dist_noz, Y_neigh = log(neigh_values_dr$Y_neigh_phy),
  num_data = length(cr_pops$damping_ratio), ac_max_ind = dist_mats$phy_max_ind, min_rho = 0, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 0.1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
damp_ratio_onlyphy_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(damp_ratio_onlyphy_noself, file = 'damp_ratio_onlyphy_noself.Rdata')
jags_model_eval(jags_obj = damp_ratio_onlyphy_noself, output_file_name = 'damp_ratio_onlyphy_noself_eval_output.pdf', obs_Y = log(cr_pops$damping_ratio))

#ALL CONVERGED

################################################################################################################################################################################################################
#### FIT MODELS WITH ENVIRONMENTAL PREDICTORS FOR COEFFICENT OF VARIATION FOR LAMBDA #######################################################################################################################################################
#full interactions model, all pops ### ### ####### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda + 1), geo_dist = dist_mats_cv$geo_dist, phy_dist = dist_mats_cv$phy_dist, Y_neigh_geo = log(neigh_values_cv$Y_neigh + 1),
  Y_neigh_phy = log(neigh_values_cv$Y_neigh + 1), num_data = length(cv_pops$lambda), geo_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), GF = cv_pops$herb_fact, num_GF = max(cv_pops$herb_fact), mat_dim = cv_pops$MatrixDimension, 
  life_expect = cv_pops$life_expect, PC_temp = cv_pops$pc1_temp, AI = log(cv_pops$AI), p_sea = cv_pops$bio_15, min_rho = 0.0001, max_rho = 1, min_alpha = 0.0001, max_alpha = 1)

#intial values
inital <- function() list(a = runif(max(cv_pops$herb_fact), 0, 1), b = runif(9, 0, 1), B = cbind(runif(4, -20, 20), runif(4, -20, 20)), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.02, 0.1), rho = runif(1, 0.02, 0.1)) 

# Parameters monitored
monitored <- c('a', "b", 'B', "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 150000
nThin <- 200
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_full_allpops <- jags(jags.data, inital, monitored, "interactions_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_full_allpops, file = 'cv_full_allpops.Rdata')
jags_model_eval(jags_obj = cv_full_allpops, output_file_name = 'cv_full_allpops_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#main effects only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda + 1), geo_dist = dist_mats_cv$geo_dist, phy_dist =dist_mats_cv$phy_dist, Y_neigh_geo = log(neigh_values_cv$Y_neigh + 1),
  Y_neigh_phy = log(neigh_values_cv$Y_neigh + 1), num_data = length(cv_pops$lambda), geo_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), GF = cv_pops$herb_fact, num_GF = max(cv_pops$herb_fact), 
  mat_dim = cv_pops$MatrixDimension, life_expect = cv_pops$life_expect, PC_temp = cv_pops$pc1_temp, AI = log(cv_pops$AI), p_sea = cv_pops$bio_15, min_rho = 0.0001, 
  max_rho = 1, min_alpha = 0.0001, max_alpha = 1)

#intial values
inital <- function() list(a = runif(max(cv_pops$herb_fact), 0, 1), b = runif(7, 0, 1), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.0001, 1), rho = runif(1, 0.0001, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_noint_allpops <- jags(jags.data, inital, monitored, "noint_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_noint_allpops, file = 'cv_noint_allpops.Rdata')
jags_model_eval(jags_obj = cv_noint_allpops, output_file_name = 'cv_noint_allpops_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#### NO FIXED EFFECTS CV LAMBDA
#no fixed effects, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda + 1), geo_dist = dist_mats_cv$geo_dist, phy_dist =dist_mats_cv$phy_dist, Y_neigh_geo = log(neigh_values_cv$Y_neigh + 1),
  Y_neigh_phy = log(neigh_values_cv$Y_neigh + 1), num_data = length(cv_pops$lambda), geo_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), min_rho = 0.0001, max_rho = 1, min_alpha = 0.0001, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, 0, 0), b = runif(2, 0, 1), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.2, 0.2), rho = runif(1, 0.001, 0.5)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 150000
nThin <- 200
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_nofixed_allpops <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_nofixed_allpops, file = 'cv_nofixed_allpops.Rdata')
jags_model_eval(jags_obj = cv_nofixed_allpops, output_file_name = 'cv_nofixed_allpops_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#no fixed effects, noself ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda+1), geo_dist = dist_mats_cv$geo_dist_noz, phy_dist = dist_mats_cv$phy_dist_noz, Y_neigh_geo = log(neigh_values_cv$Y_neigh_geo + 1),
  Y_neigh_phy = log(neigh_values_cv$Y_neigh_phy + 1), num_data = length(cv_pops$lambda), geo_max_ind = dist_mats_cv$geo_max_ind, phy_max_ind =dist_mats_cv$phy_max_ind, 
  min_rho = 0, max_rho = 1, min_alpha = 0.01, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.01, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 550000
nThin <- 1000
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_nofixed_noself <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

cv_nofixed_noself <- update(cv_nofixed_noself, n.iter = 1000000, n.thin = 2000)

save(cv_nofixed_noself, file = 'cv_nofixed_noself.Rdata')
jags_model_eval(jags_obj = cv_nofixed_noself, output_file_name = 'cv_nofixed_noself_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

##ALL COVERGED

#no fixed effects geo only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda+1), ac_dist =dist_mats_cv$geo_dist, Y_neigh = log(neigh_values_cv$Y_neigh + 1),
  num_data = length(cv_pops$lambda), ac_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_onlygeo_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_onlygeo_allpops, file = 'cv_onlygeo_allpops.Rdata')
jags_model_eval(jags_obj = cv_onlygeo_allpops, output_file_name = 'cv_onlygeo_allpops_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#no fixed effects phy only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda+1), ac_dist =dist_mats_cv$phy_dist, Y_neigh = log(neigh_values_cv$Y_neigh + 1),
  num_data = length(cv_pops$lambda), ac_max_ind = rep(dim(neigh_values_cv$Y_neigh)[2], dim(neigh_values_cv$Y_neigh)[1]), min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_onlyphy_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_onlyphy_allpops, file = 'cv_onlyphy_allpops.Rdata')
jags_model_eval(jags_obj = cv_onlyphy_allpops, output_file_name = 'cv_onlyphy_allpops_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#no fixed effects geo only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda+1), ac_dist = dist_mats_cv$geo_dist_noz, Y_neigh = log(neigh_values_cv$Y_neigh_geo + 1),
  num_data = length(cv_pops$lambda), ac_max_ind = dist_mats_cv$geo_max_ind, min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_onlygeo_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_onlygeo_noself, file = 'cv_onlygeo_noself.Rdata')
jags_model_eval(jags_obj = cv_onlygeo_noself, output_file_name = 'cv_onlygeo_noself_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#no fixed effects phy only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = log(cv_pops$cv_lambda+1), ac_dist = dist_mats_cv$phy_dist_noz, Y_neigh = log(neigh_values_cv$Y_neigh_phy + 1),
  num_data = length(cv_pops$lambda), ac_max_ind = dist_mats_cv$phy_max_ind, min_rho = 0.001, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, 0.4, 0.4), b = runif(1, 0.75, 0.75), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0.002, 0.002)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
cv_onlyphy_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(cv_onlyphy_noself, file = 'cv_onlyphy_noself.Rdata')
jags_model_eval(jags_obj = cv_onlyphy_noself, output_file_name = 'cv_onlyphy_noself_eval_output.pdf', obs_Y = log(cv_pops$cv_lambda+1))

#ALL CONVERGED

#DO ELAST WITH THIS FRAMEWORK so know all models are the same and will make the methods easier to write
################################################################################################################################################################################################################
#### FIT MODELS WITH ENVIRONMENTAL PREDICTORS FOR PC1 #######################################################################################################################################################
#full interactions model, all pops ### ### ####### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = neigh_values_PC1$Y_neigh,
  Y_neigh_phy = neigh_values_PC1$Y_neigh, num_data = length(cr_pops$elast_PC1), geo_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), GF = cr_pops$herb_fact, num_GF = max(cr_pops$herb_fact), mat_dim = cr_pops$MatrixDimension, 
  life_expect = cr_pops$life_expect, PC_temp = cr_pops$pc1_temp, AI = log(cr_pops$AI), p_sea = cr_pops$bio_15, min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(max(cr_pops$herb_fact), -20, 20), b = runif(9, -20, 20), B = cbind(runif(4, -20, 20), runif(4, -20, 20)), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", 'B', "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_full_allpops <- jags(jags.data, inital, monitored, "interactions_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_full_allpops, file = 'PC1_full_allpops.Rdata')
jags_model_eval(jags_obj = PC1_full_allpops, output_file_name = 'PC1_full_allpops_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#main effects only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = neigh_values_PC1$Y_neigh,
  Y_neigh_phy = neigh_values_PC1$Y_neigh, num_data = length(cr_pops$elast_PC1), geo_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), GF = cr_pops$herb_fact, num_GF = max(cr_pops$herb_fact), 
  mat_dim = cr_pops$MatrixDimension, life_expect = cr_pops$life_expect, PC_temp = cr_pops$pc1_temp, AI = log(cr_pops$AI), p_sea = cr_pops$bio_15, min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(max(cr_pops$herb_fact), -20, 20), b = runif(7, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_noint_allpops <- jags(jags.data, inital, monitored, "noint_2ac_decay_models.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_noint_allpops, file = 'PC1_noint_allpops.Rdata')
jags_model_eval(jags_obj = PC1_noint_allpops, output_file_name = 'PC1_noint_allpops_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#### NO FIXED EFFECTS PC1
#no fixed effects, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, geo_dist = dist_mats$geo_dist, phy_dist = dist_mats$phy_dist, Y_neigh_geo = neigh_values_PC1$Y_neigh,
  Y_neigh_phy = neigh_values_PC1$Y_neigh, num_data = length(cr_pops$elast_PC1), geo_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), 
  phy_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0.001, 1), rho = runif(1, 0.001, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 150000
nThin <- 200
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_nofixed_allpops <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_nofixed_allpops, file = 'PC1_nofixed_allpops.Rdata')
jags_model_eval(jags_obj = PC1_nofixed_allpops, output_file_name = 'PC1_nofixed_allpops_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#no fixed effects, noself ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, geo_dist = dist_mats$geo_dist_noz, phy_dist = dist_mats$phy_dist_noz, Y_neigh_geo = neigh_values_PC1$Y_neigh_geo,
  Y_neigh_phy = neigh_values_PC1$Y_neigh_phy, num_data = length(cr_pops$elast_PC1), geo_max_ind = dist_mats$geo_max_ind, phy_max_ind = dist_mats$phy_max_ind, 
  min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_nofixed_noself <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_nofixed_noself, file = 'PC1_nofixed_noself.Rdata')
jags_model_eval(jags_obj = PC1_nofixed_noself, output_file_name = 'PC1_nofixed_noself_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#no fixed effects geo only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, ac_dist = dist_mats$geo_dist, Y_neigh = neigh_values_PC1$Y_neigh,
  num_data = length(cr_pops$elast_PC1), ac_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_onlygeo_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_onlygeo_allpops, file = 'PC1_onlygeo_allpops.Rdata')
jags_model_eval(jags_obj = PC1_onlygeo_allpops, output_file_name = 'PC1_onlygeo_allpops_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#no fixed effects phy only, all pops ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, ac_dist = dist_mats$phy_dist, Y_neigh = neigh_values_PC1$Y_neigh,
  num_data = length(cr_pops$elast_PC1), ac_max_ind = rep(dim(neigh_values_PC1$Y_neigh)[2], dim(neigh_values_PC1$Y_neigh)[1]), min_rho = 0, 
  max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_onlyphy_allpops <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_onlyphy_allpops, file = 'PC1_onlyphy_allpops.Rdata')
jags_model_eval(jags_obj = PC1_onlyphy_allpops, output_file_name = 'PC1_onlyphy_allpops_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#no fixed effects geo only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, ac_dist = dist_mats$geo_dist_noz, Y_neigh = neigh_values_PC1$Y_neigh_geo,
  num_data = length(cr_pops$elast_PC1), ac_max_ind = dist_mats$geo_max_ind, min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_onlygeo_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_onlygeo_noself, file = 'PC1_onlygeo_noself.Rdata')
jags_model_eval(jags_obj = PC1_onlygeo_noself, output_file_name = 'PC1_onlygeo_noself_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

#no fixed effects phy only, no-self ### ### ### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
jags.data <- list(Y = cr_pops$elast_PC1, ac_dist = dist_mats$phy_dist_noz, Y_neigh = neigh_values_PC1$Y_neigh_phy,
  num_data = length(cr_pops$elast_PC1), ac_max_ind = dist_mats$phy_max_ind, min_rho = 0, max_rho = 1, min_alpha = 0, max_alpha = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(1, -20, 20), tau.global = runif(1, 0.0001, 100), 
  rho = runif(1, 0, 0.1)) 

# Parameters monitored
monitored<- c('a', "b", "rho", "sd.global", "residual", "fit.new", "fit", "bpvalue")

#mcmc params
nIter <- 100000
nThin <- 100
nBurnIn <- 50000
nChains <- 3

#call to JAGS
set.seed(123)
setwd(analysis_ind)
PC1_onlyphy_noself <- jags(jags.data, inital, monitored, "nofixed_1ac_decay_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

save(PC1_onlyphy_noself, file = 'PC1_onlyphy_noself.Rdata')
jags_model_eval(jags_obj = PC1_onlyphy_noself, output_file_name = 'PC1_onlyphy_noself_eval_output.pdf', obs_Y = cr_pops$elast_PC1)

#ALL CONVERGED

