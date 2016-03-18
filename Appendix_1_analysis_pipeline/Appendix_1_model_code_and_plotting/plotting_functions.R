
#plot the credible intervals for each predictor with sensible names to help see what these results are really telling us.
#get all the model objects for those taht converged
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
library(colorspace)
#automated attempt to find the right file path, might fail and manual path will have to be entered 
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

pre_run_mods <- paste0(model_file_loc, 'pre-run_model_output') #if you want use models you have fit yourself the file path needs to be changed
print('You will need to download ESIR 30 arc-second BioClim data from http://worldclim.org/current')
bioclimIN = 'your/file/path/to/bioClim/data'

#get the data_sets I need
#ge the tree and 
setwd(data_file_loc)
Tree <- read.nexus("PMDbl260314.nex")

setwd(model_file_loc)
obj_name <- load('combined_bc_ai_all_responses.Rdata')

#===============================================================================
#LIST OF MODELS THAT HAVE CONVEREGED

#get the file names from the workingFiles Directory where all the models are kept
setwd(pre_run_mods)
PC1_elast_file_names <- list.files(pattern = '^PC1_.*\\.Rdata$')
lambda_file_names <- list.files(pattern = '^lambda_.*\\.Rdata$')
damp_ratio_file_names <- list.files(pattern = '^damp.*\\.Rdata$')
cv_file_names <- list.files(pattern = '^cv.*\\.Rdata$')

PC1_jags_ob_names <- character()
lambda_jags_ob_names <- character()
damp_jags_ob_names <- character()
cv_jags_ob_names <- character()

for(i in 1:length(PC1_elast_file_names)) PC1_jags_ob_names[i] <- load(PC1_elast_file_names[i])
for(i in 1:length(lambda_file_names)) lambda_jags_ob_names[i] <- load(lambda_file_names[i])
for(i in 1:length(damp_ratio_file_names)) damp_jags_ob_names[i] <- load(damp_ratio_file_names[i])
for(i in 1:length(cv_file_names)) cv_jags_ob_names[i] <- load(cv_file_names[i])

setwd(model_file_loc)
source('custom_jags_object_function.r')#function to make custom objects that the rest of the plotting requires
setwd(model_file_loc)
source('plotting_helper_functions.R')
#DIAGNOSTIC PLOTS FOR EACH MODEL, ONE SET OF 4 ON EACH PAGE WITH A PAGE TITLE THAT CONTAINS THE MODE SPECIFICIATION, THE DIC AND THE BP-VALUE 
#FOR THE PLOTS HAVE THE FITTED VS RESIDUAL, FITTED VS, PREDICTED (QUANTILES), R^2 HISTOGRAM, LINE PLOT OF COEFFCIENT QUANTILES
#make a set of diagnostic plots for each mode
diag_plots(jags_ob_names = PC1_jags_ob_names, output_file = 'diagnositc_summary_plots_elast_PC1.pdf', output_loc = model_file_loc, working_files_loc = pre_run_mods, Y = cr_pops$elast_PC1)
diag_plots(jags_ob_names = lambda_jags_ob_names, output_file = 'diagnositc_summary_plots_lambda.pdf', output_loc = model_file_loc, working_files_loc = pre_run_mods, Y = log(cr_pops$lambda))
diag_plots(jags_ob_names = damp_jags_ob_names, output_file = 'diagnositc_summary_plots_damp_ratio.pdf', output_loc = model_file_loc, working_files_loc = pre_run_mods, Y = log(cr_pops$damping_ratio))
diag_plots(jags_ob_names = cv_jags_ob_names, output_file = 'diagnositc_summary_plots_cv_lambda.pdf', output_loc = model_file_loc, working_files_loc = pre_run_mods, Y = log(cr_pops$cv_lambda[!is.na(cr_pops$cv_lambda)]+1))

#make a list of jags objects to pass to plotting plotting_helper_functions
setwd(pre_run_mods)
list_jags_ob <- list()
list_jags_ob$PC1_elast <- list()
list_jags_ob$lambda <- list()
list_jags_ob$damp_ratio <- list()
list_jags_ob$cv_lambda <- list()

for(i in seq_along(PC1_jags_ob_names)) list_jags_ob$PC1_elast[[i]] <- custom_jags_object(eval(parse(text = as.character(PC1_jags_ob_names[i]))), cr_pops$elast_PC1)  
for(i in seq_along(lambda_jags_ob_names)) list_jags_ob$lambda[[i]] <- custom_jags_object(eval(parse(text = as.character(lambda_jags_ob_names[i]))), log(cr_pops$lambda))  
for(i in seq_along(damp_jags_ob_names)) list_jags_ob$damp_ratio[[i]] <- custom_jags_object(eval(parse(text = as.character(damp_jags_ob_names[i]))), log(cr_pops$damping_ratio))  
for(i in seq_along(cv_jags_ob_names)) list_jags_ob$cv_lambda[[i]] <- custom_jags_object(eval(parse(text = as.character(cv_jags_ob_names[i]))), log(cr_pops$cv_lambda[!is.na(cr_pops$cv_lambda)] + 1))  

model_spec <- c('main_int', 'phygeo-all_pops', 'phygeo-no_self', 'main', 'geo-all_pops', 'geo-no_self', 'phy-all_pops', 'phy-no_self')
color_palett_fill = rainbow_hcl(n = length(model_spec), start = 30, end = 290, alpha = 0.5)
color_palett_border = rainbow_hcl(n = length(model_spec), start = 30, end = 290, alpha = 1)
color_palett_samp_lines = rainbow_hcl(n = length(model_spec), start = 30, end = 290, alpha = 0.2)
#make a plot with the density of R-sq for each measure 
R_sq_plotter(model_spec = model_spec, color_palett_fill = color_palett_fill, color_palett_border = color_palett_border, list_jags_obj = list_jags_ob, 
  file_name = 'R_sq_all_responses.pdf', output_loc = model_file_loc)
  
#alternate R^2 plot using quantiles
R_sq_quant_plot(list_jags_obj = list_jags_ob, color_palett_fill = color_palett_fill, color_palett_border = color_palett_border, model_spec = model_spec, include_models = c(1,2,3,5,6,7,8),
  file_name = 'R_sq_quant_plot_all.pdf', output_loc = model_file_loc)

#make a plot of parameter estimates 
quant_plots(list_jags_obj = list_jags_ob, color_palett_fill = color_palett_fill, color_palett_border = color_palett_border, model_spec = model_spec,
  file_name = 'quant_coef_plot_all.pdf', output_loc = model_file_loc)

#plot the decay curves for each model for each measure
decay_plots(list_jags_obj = list_jags_ob, color_palett_samp_lines = color_palett_samp_lines, model_spec = model_spec, names_of_measures = c('SPG', expression('ln('~lambda~')'), expression('ln('~rho~')'), expression('CV('~lambda~')')), 
  file_name = 'decay_models_plotted_all.pdf', output_loc = model_file_loc)  

#make a plot of just the nofixed_noself models, put the rest in the appendix
decay_plots_single_model(list_jags_obj = list_jags_ob, color = color_palett_samp_lines[3], model_spec = model_spec, plotted_model = 3, 
  names_of_measures = c('SPG', expression(lambda), expression(rho), expression('CV('~lambda~')')), file_name = 'decay_models_nofixed_noself.pdf', output_loc = model_file_loc)

#make the plot of the environemntal variables 
setwd(bioclimIN)
files <- list.files(bioclimIN)
files <- files[!files == "info"]
files <- files[!files == "pairsBioClimPlot.pdf"]
env_pred_vars <- stack(files) #makes the raster stack of bioclim files

#creat the PCA of tempature to make the byplot. Each site should only be used once in each 
cr_pops$unique_locations <- sapply(seq_along(cr_pops$lon_dd), FUN = function(x) paste0(cr_pops$lon_dd[x], ':', cr_pops$lat_dd[x]))
unique_loc <- unique(cr_pops$unique_locations)
#finds the inds of those unique_AI's and takes only the first ind from each one
site_inds <- as.numeric(sapply(unique_loc, FUN = function(x) which(cr_pops$unique_locations == x)[1]))
#get unique site data for pca
pca_data <- cr_pops[site_inds,]
#look at the tempature populations
temp_pca <- prcomp(formula = ~ bio_1 + bio_3 + bio_4 + bio_5 + bio_6 + bio_7 + bio_8 + bio_9, data = pca_data, center = TRUE, scale = TRUE)
#make the precip pca
precip_pca = prcomp(formula =  ~ bio_12 + bio_13 + bio_14 + bio_15 + bio_16 + bio_17 + bio_18 + bio_19 + AI, data = pca_data, center = TRUE, scale = TRUE)

temp_space_PCA_plotter(env_pred = env_pred_vars, temp_pca = temp_pca, precip_pca = precip_pca, data_source = cr_pops, darkness = 0, filename = 'temp_season_PCA_mapped_550.pdf', 
  output_loc = model_file_loc)
