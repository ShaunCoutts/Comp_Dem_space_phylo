#Data cleaning script for the full dataset of trees, shrubs, herbs and ferns, succulents and palms to try and automate at least part of the clean up
#get the libraries needed, mostly to test phylogeny 
library(ape)
library(caper)
library(phytools)
library(gtools)
library(GEOmap)
library(sp)
library(maptools)
library(dismo)
library(raster)
library(rgdal)
library(rgl)
data(wrld_simpl)
library(rworldmap)
library(abind)
library(Matrix)
library(plyr)

#two functions, the first to calculate elasticites and sensitivies, second to calculate life expectancy amoung other things
matrixElementPerturbation <- function(matU, matF, matC=NULL,pert=0.001){
	#Function to calculate matrix element level sensitivities and elasticities
	matA=matU+matF+matC
	aDim=dim(matA)[1]
	fakeA=matA
	sensA=elasA=matrix(NA,aDim,aDim)
	lambda=Re(eigen(matA)$values[1])
	propU=matU/matA
	propU[is.nan(propU)]=NA
	propProg=propRetrog=propU
	propProg[upper.tri(propU,diag=T)]=NA
	propRetrog[lower.tri(propU,diag=T)]=NA
	propStasis=matrix(Diagonal(aDim)*diag(propU),aDim,aDim)
	propF=matF/matA
	propF[is.nan(propF)]=NA
	propC=matC/matA
	propC[is.nan(propC)]=NA
	for (i in 1:aDim){
		for (j in 1:aDim){
			fakeA=matA
			fakeA[i,j]=fakeA[i,j]+pert
			lambdaPert=eigen(fakeA)$values[1]
			sensA[i,j]=(lambda-lambdaPert)/(matA[i,j]-fakeA[i,j])
		}
	}
	sensA=Re(sensA)
	elasA=sensA*matA/lambda
	out = data.frame("SStasis"=NA,"SProgression"=NA,"SRetrogression"=NA,"SFecundity"=NA,"SClonality"=NA,
	"EStasis"=NA,"EProgression"=NA,"ERetrogression"=NA,"EFecundity"=NA,"EClonality"=NA)
	out$SStasis=sum(sensA*propStasis,na.rm=T)
	out$SRetrogression=sum(sensA*propRetrog,na.rm=T)
	out$SProgression=sum(sensA*propProg,na.rm=T)
	out$SFecundity=sum(sensA*propF,na.rm=T)
	out$SClonality=sum(sensA*propC,na.rm=T)
	out$EStasis=sum(elasA*propStasis,na.rm=T)
	out$EProgression=sum(elasA*propProg,na.rm=T)
	out$ERetrogression=sum(elasA*propRetrog,na.rm=T)
	out$EFecundity=sum(elasA*propF,na.rm=T)
	out$EClonality=sum(elasA*propC,na.rm=T)
	return(out)
}
lifeTimeRepEvents <- function(matU, matF, startLife = 1){
	#Function to determine probability of reaching reproduction, age at maturity and reproductive lifespan (Code adapted from H. Caswell's matlab code):
	uDim = dim(matU)[1]
	surv = colSums(matU)
	repLifeStages = colSums(matF)
	repLifeStages[which(repLifeStages>0)] = 1
	if(missing(matF) | missing(matU)){stop('matU or matF missing')}
	if(sum(matF,na.rm=T)==0){stop('matF contains only 0 values')}
	#Probability of survival to first reprod event
	Uprime = matU
	Uprime[,which(repLifeStages==1)] = 0
	Mprime = matrix(0,2,uDim)
	for (p in 1:uDim[1]) {
		if (repLifeStages[p]==1) Mprime[2,p] = 1 else
		Mprime[1,p] = 1-surv[p]
	}	
	Bprime = Mprime%*%(ginv(diag(uDim)-Uprime))
	pRep = Bprime[2,startLife]
	out = data.frame(pRep = pRep)
	#Age at first reproduction (La; Caswell 2001, p 124)
	  D = diag(c(Bprime[2,]))
	Uprimecond = D%*%Uprime%*%ginv(D)
	expTimeReprod = colSums(ginv(diag(uDim)-Uprimecond))
	La = expTimeReprod[startLife]
	out$La = La
	#Mean life expectancy conditional on entering the life cycle in the first reproductive stage
	firstRepLifeStage = min(which(repLifeStages==1))
	N = solve(diag(uDim[1])-matU)
	meanRepLifeExpectancy = colSums(N)[firstRepLifeStage]
	out$meanRepLifeExpectancy = meanRepLifeExpectancy
	#Life expectancy from mean maturity
	remainingMatureLifeExpectancy = colSums(N)[startLife]-La
	out$remainingMatureLifeExpectancy = remainingMatureLifeExpectancy
	return(out)
}

# This is an attempt to automaticaly find the file location of the extracted folder, this requires root access, it works on Linux (and so probably on mac as well), 
# not sure how it will go on window
print('this is an automated attemp to find and set the folder location')
print('It requires root access. The file path can also be set manually at line')
print('107 in file data_extract_clean_up_first_pass.R')
file_loc = paste0(strsplit(system('sudo updatedb;locate PMDbl260314A1.nex', intern = TRUE), '/')[[1]], '/')
data_file_loc = ''
for(i  in 1:(length(file_loc) - 1)) data_file_loc = paste0(data_file_loc, file_loc[i])
#MANUAL SETTING
#data_file_loc = '/your/path/to/folder/here'
# bioclimIN <- paste0(data_file_loc, 'bioClim_data')
# aiIN <- paste0(data_file_loc, 'AI_annual/ai_yr') 
print('You will need to download ESIR 30 arc-second BioClim data from http://worldclim.org/current')
print('You will also need Aridity Index from http://www.cgiar-csi.org/data/global-aridity-and-pet-database')
bioclimIN = 'your/file/path/to/bioClim/data'
aiIN = 'your/file/path/to/bioClim/data'

setwd(data_file_loc)
Tree <- read.nexus("PMDbl260314A1.nex")
obj_name <- load('COMPADRE_Oct_24_2014_Version_.RData')

#seems to be some weridness with the numbers since numbers don't match up so will look for the dimnetions to see where they stop matching up
#mat_dims <- compadre$metadata$MatrixDimension[-c(2792:2798)]
#mat_dims <- mat_dims[-c(5377:5379)]
#a <- cbind(sapply(1:length(compadre$mat), FUN = function(x) dim(compadre$mat[[x]]$matA)[1]), mat_dims)
#a <- cbind(a, a[, 1] == a[, 2], 1:length(a[,1])) 
#get full matchup taking out those two small sections, likey to be missed species in this early version of COMPADRE

used_dat <- compadre$metadata
#remove a few rows so I can get the rows to match up with the list of matricies, and have the right dimentions, so lieky to be match
used_dat <- used_dat[-c(2792:2798), ]
used_dat <- used_dat[-c(5377:5379), ]
#test we get match up
a <- cbind(sapply(1:length(compadre$mat), FUN = function(x) dim(compadre$mat[[x]]$matA)[1]), used_dat$MatrixDimension)
a <- cbind(a, a[, 1] == a[, 2])
prod(a[,3])
#start cleaning data
used_dat$row_id <- 1:dim(used_dat)[1]
used_dat$spp = as.character(used_dat$SpeciesAuthor)
used_dat$spp <- gsub('_1', '', x = used_dat$spp, fixed = TRUE)
used_dat$spp <- gsub('_2', '', x = used_dat$spp, fixed = TRUE)
used_dat$spp <- gsub('_3', '', x = used_dat$spp, fixed = TRUE)
used_dat$spp <- gsub('_4', '', x = used_dat$spp, fixed = TRUE)
used_dat$spp <- gsub('_5', '', x = used_dat$spp, fixed = TRUE)

#series of filters to get the species with cirteria I need
used_dat <- used_dat[which(used_dat$GrowthType == 'Herbaceous perennial' | used_dat$GrowthType == 'Tree' | used_dat$GrowthType == 'Shrub' | used_dat$GrowthType == 'Succulent' | used_dat$GrowthType == 'Palm'), ]
used_dat <- used_dat[!is.na(used_dat$LatMin), ]
used_dat <- used_dat[-which(used_dat$LonMin == 0 & used_dat$LonSec == 0 & used_dat$LatMin == 0 & used_dat$LatSec == 0), ]
used_dat <- used_dat[which(used_dat$StudyDuration > 2), ]
used_dat <- used_dat[which(used_dat$MatrixSplit == 'Divided'), ]
used_dat <- used_dat[which(used_dat$MatrixComposite == 'Mean'), ]
used_dat <- used_dat[which((used_dat$MatrixEndYear - used_dat$MatrixStartYear) > 1), ]
used_dat <- used_dat[-grep(';', used_dat$Population), ] #remove populaitons with ';' as they show compound populaitons 
used_dat <- used_dat[!grepl('Mowing', used_dat$MatrixTreatment),]
used_dat <- used_dat[!grepl('Seed', used_dat$MatrixTreatment),]
used_dat <- used_dat[!grepl('Herbicide', used_dat$MatrixTreatment),]
used_dat <- used_dat[!grepl('burn', used_dat$Observation, ignore.case = TRUE),]
used_dat <- used_dat[!grepl('burn', used_dat$MatrixTreatment, ignore.case = TRUE),]
used_dat <- used_dat[!grepl('fire', used_dat$MatrixTreatment, ignore.case = TRUE),]
used_dat <- used_dat[!grepl('NDY', used_dat$MatrixTreatment, ignore.case = TRUE),]
used_dat <- used_dat[!grepl('The GPS coordinates were', used_dat$Observation), ]
used_dat <- used_dat[!grepl('approximated', used_dat$Observation), ]
used_dat <- used_dat[!grepl('Reproduction rates not given', used_dat$Observation), ]
used_dat <- used_dat[!grepl('Survival higher than one because many transitions are asexual reproduction', used_dat$Observation), ]
used_dat <- used_dat[!grepl('Fecundity is missing', used_dat$Observation), ]
used_dat <- used_dat[!grepl('Geolocation', used_dat$Observation), ]
used_dat <- used_dat[!grepl('latitude', used_dat$Observation), ]
used_dat <- used_dat[!grepl('Chamaecrista_keyensis', used_dat$spp), ]#very confusing about what is going on with this one 
used_dat <- used_dat[which(used_dat$SurvivalIssue < 1.05),]
used_dat <- used_dat[which(used_dat$MatrixDimension > 2),]

#check that all the species in used_dat are in the phylogeny 
used_dat = used_dat[used_dat$SpeciesAuthor %in% Tree$tip.label, ] #just one drops out
used_data_ns = used_dat #used latter to get the matricies for demographic metric calculations

#extract matricies using the row_id
used_mats <- list()
for(i in seq_along(used_dat$row_id)){
	used_mats[[i]] <- compadre$mat[[used_dat$row_id[i]]]
	used_mats[[i]]$row_id <- used_dat$row_id[i]
}

pert_wrapper <- function(x){
	a <- matrixElementPerturbation(matU = used_mats[[x]]$matU, matF = used_mats[[x]]$matF, matC = used_mats[[x]]$matC)
	cbind(a, row_id = used_mats[[x]]$row_id)
}
perts <- sapply(1:length(used_mats), FUN = function(x) tryCatch(pert_wrapper(x), error = function(e) NA, warning = function(w) NA))
perts <- ldply(perts, data.frame)
used_dat <- merge(x = used_dat, y = perts[, 1:11], by = 'row_id')

life_wrapper <- function(x){
	a[c('life_expect')] <- lifeTimeRepEvents(matU = used_mats[[x]]$matU, matF = used_mats[[x]]$matF)$meanRepLifeExpectancy
	a[c('row_id')] <- used_mats[[x]]$row_id
	b <- a[c('life_expect', 'row_id')]
	b
}
life <- sapply(1:length(used_mats), FUN = function(x) tryCatch(life_wrapper(x), error = function(e) NA, warning = function(w) NA))
life_mat <- matrix(ncol = 2, nrow = length(life))
for(i in 1:length(life)){
	life_mat[i, ] <- life[[i]]
}
life_df <- data.frame(life_expect = life_mat[, 1], row_id = life_mat[, 2])
used_dat <- merge(x = used_dat, y = life_df, by = 'row_id')
#note merge() also pushes out the NA's 
used_dat$SumE <- apply(used_dat[, c('EStasis', 'EProgression', 'ERetrogression', 'EFecundity', 'EClonality')], MARGIN = 1, FUN = sum)
#There are a few degenerate matricies that don't work too well
used_dat <- used_dat[which(used_dat$EStasis >= 0 & used_dat$EStasis <= 1 & used_dat$EProgression >= 0 & used_dat$EProgression <= 1 & used_dat$ERetrogression >= 0 & used_dat$ERetrogression <= 1 & used_dat$EFecundity >= 0 & used_dat$EFecundity <= 1 & used_dat$EClonality >= 0 & used_dat$EClonality <= 1), ]
used_dat <- used_dat[-which(used_dat$SumE < 0.99 | used_dat$SumE > 1.05), ]
used_dat <- used_dat[which(used_dat$life_expect > 0), ]

geo_converter <- function(DMS){
	latSec <- ifelse(is.na(DMS$LatSec), 0, DMS$LatSec)
	lonSec <- ifelse(is.na(DMS$LonSec), 0, DMS$LonSec)
	lat_dd <- ifelse(DMS$LatNS == 'N', DMS$LatDeg + (DMS$LatMin/60) + (latSec/3600), -(DMS$LatDeg + (DMS$LatMin/60) + (latSec/3600))) 
	lon_dd <- ifelse(DMS$LonWE == 'E', DMS$LonDeg + (DMS$LonMin/60) + (lonSec/3600), -(DMS$LonDeg + (DMS$LonMin/60) + (lonSec/3600)))
	c(lat_dd = lat_dd, lon_dd = lon_dd, row_id = DMS$row_id)
}
geo_dd <- t(sapply(seq_along(used_dat$row_id), FUN = function(x) geo_converter(used_dat[x, c('LatDeg', 'LatMin', 'LatSec', 'LatNS', 'LonDeg', 'LonMin', 'LonSec', 'LonWE', 'row_id')])))
used_dat <- merge(x = used_dat, y = geo_dd, by = 'row_id')

#take out a few of the coloumns I am not going to use
used_dat <- used_dat[, c('row_id', 'SpeciesAuthor', 'SpeciesAccepted', 'GenusAccepted', 'Genus', 'Family', 'Order', 'Class', 
	'Phylum', 'GrowthType', 'StudyDuration', 'StudyStart', 'StudyEnd', 'MatrixStartYear', 'MatrixEndYear', 'Population', 'Altitude', 'MatrixDimension', 'spp', 
	'SStasis', 'SProgression', 'SRetrogression', 'SFecundity', 'SClonality', 'EStasis', 'EProgression', 'ERetrogression', 'EFecundity', 
	'EClonality', 'life_expect', 'SumE', 'lat_dd', 'lon_dd')]


#turn the used_pops into a spatial object
coordinates(used_dat) <- ~lon_dd + lat_dd #turn sppLoc into a SpatialPointsDataFrame
proj4string(used_dat) <- '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs' #set the CRS to world wide lat long system, because that is what bio_clim data is in


#get bio_clim variables
setwd(bioclimIN)
#BIO1 = Annual Mean Temperature (C*10)
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) (C*10)
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month (C*10)
#BIO6 = Min Temperature of Coldest Month (C*10)
#BIO7 = Temperature Annual Range (BIO5-BIO6) (C*10)
#BIO8 = Mean Temperature of Wettest Quarter (C*10)
#BIO9 = Mean Temperature of Driest Quarter (C*10)
#BIO10 = Mean Temperature of Warmest Quarter (C*10)
#BIO11 = Mean Temperature of Coldest Quarter (C*10)
#BIO12 = Annual Precipitation (mm)
#BIO13 = Precipitation of Wettest Month (mm)
#BIO14 = Precipitation of Driest Month (mm)
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter (mm)
#BIO17 = Precipitation of Driest Quarter(mm)
#BIO18 = Precipitation of Warmest Quarter (mm)
#BIO19 = Precipitation of Coldest Quarter (mm)

#get files, This part is long so I save the two spatial points objects produced
files <- list.files(bioclimIN)
files <- files[!files == "info"]
files <- files[!files == "pairsBioClimPlot.pdf"]
env_pred <- stack(files) #makes the raster stack of bioclim files
bc_2km_mean <- extract(x = env_pred, y = used_dat, method = 'simple', buffer = 2000, fun = mean, sp = TRUE) #takes values from surrounding cells to help offset location uncertianty, and stop NA vlaues for populations on the coast
setwd(data_file_loc)
save(bc_2km_mean, file = 'spatialPoints_bc.Rdata')

setwd(aiIN)
files <- list.files(aiIN)
ai_image <- raster('w001001.adf')
#plot(ai_image*0.0001)
ai_2km_mean <- extract(x = ai_image, y = used_dat, method = 'simple', buffer = 2000, fun = mean, sp = TRUE) #takes values from surrounding cells to help offset location uncertianty, and stop NA vlaues for populations on the coast
setwd(data_file_loc) 
save(ai_2km_mean, file = 'spatialPoints_ai.Rdata')

#################################################################################################################################################
#################################################################################################################################################
#get extracted objects so I can just start here in future
#setwd('data_file_loc') 
#loadedObj <- load('spatialPoints_ai.Rdata')
#loadedObj <- load('spatialPoints_bc.Rdata')

#do a bit of object suffeling to ensure the objects are preserved and I merge the right ones
cr_pops <- bc_2km_mean
cr_pops <- merge(x = bc_2km_mean@data, y = ai_2km_mean@data[, c('w001001', 'row_id')], by = 'row_id', all.x = TRUE)
cr_pops[, c('lon_dd', 'lat_dd')] <- bc_2km_mean@coords
cr_pops <- cr_pops[!is.na(cr_pops$bio_1), ]
cr_pops$AI <- cr_pops$w001001
cr_pops$EFecClo <- cr_pops$EFecundity + cr_pops$EClonality
cr_pops$SFecClo <- cr_pops$SFecundity + cr_pops$SClonality
#get the PCA scores
#look for a axis in the elasticities that could summarise demography in one number that could be used as a response 
elast_pca <- prcomp(formula = ~EProgression + ERetrogression + EStasis + EFecClo, data = cr_pops, center = TRUE, scale = TRUE)
print(elast_pca)
summary(elast_pca) 
#Rotations
#Rotation:
#                      PC1         PC2         PC3       PC4
#EProgression   -0.5609341  0.08210908  0.67625458 0.4704156
#ERetrogression -0.1774500  0.87753796 -0.39585345 0.2043004
#EStasis         0.6495852 -0.03844182  0.01539378 0.7591603
#EFecClo        -0.4815636 -0.47085822 -0.62108198 0.4008070

#creat the PCA of tempature to make the byplot. Each site should only be used once in each 
cr_pops$unique_locations <- sapply(seq_along(cr_pops$lon_dd), FUN = function(x) paste0(cr_pops$lon_dd[x], ':', cr_pops$lat_dd[x]))
unique_loc <- unique(cr_pops$unique_locations)
#finds the inds of those unique_AI's and takes only the first ind from each one
site_inds <- as.numeric(sapply(unique_loc, FUN = function(x) which(cr_pops$unique_locations == x)[1]))
#get unique site data for pca
pca_data <- cr_pops[site_inds,]
#
temp_pca <- prcomp(formula = ~bio_1 + bio_3 + bio_4 + bio_5 + bio_6 + bio_7 + bio_8 + bio_9, data = pca_data, center = TRUE, scale = TRUE)
print(temp_pca)
summary(temp_pca)
#Importance of components:
#                           PC1    PC2     PC3     PC4     PC5     PC6     PC7 PC8
# Standard deviation     2.3793 1.1810 0.79878 0.49288 0.20570 0.12285 0.07745 3.447e-16
# Proportion of Variance 0.7076 0.1743 0.07976 0.03037 0.00529 0.00189 0.00075 0.000e+00
# Cumulative Proportion  0.7076 0.8820 0.96171 0.99207 0.99736 0.99925 1.00000 1.000e+00
#                              
#                             
biplot(temp_pca, col = c("gray", "black")) 

#PC1, which explains 71% of the variation, seems like variability (bio_4) and temp range (bio_7) are stronly positivly correlated, while mean temp is correlated 
#with the other tempature measures, all of which are negativly correlated with seasonality. Represents a troplical vs temperate gradient in the data 

#calculate the PCA1 and 2 score for each population
cr_pops$pc1_temp <- as.numeric(sapply(cr_pops$unique_locations, FUN = function(x) temp_pca$x[which(pca_data$unique_locations == x),'PC1']))
cr_pops$pc2_temp <- as.numeric(sapply(cr_pops$unique_locations, FUN = function(x) temp_pca$x[which(pca_data$unique_locations == x),'PC2'])) 

precip_pca <- prcomp(formula = ~ bio_12 + bio_13 + bio_14 + bio_15 + bio_16 + bio_17 + bio_18 + bio_19 + AI, data = pca_data, center = TRUE, scale = TRUE)
summary(precip_pca)
#Importance of components:
#                          PC1    PC2     PC3    PC4     PC5    PC6     PC7     PC8     PC9
#Standard deviation     2.4006 1.3908 0.81358 0.6264 0.40079 0.2700 0.09190 0.06170 0.05355
#Proportion of Variance 0.6403 0.2149 0.07355 0.0436 0.01785 0.0081 0.00094 0.06170 0.05355
#Cumulative Proportion  0.6403 0.8552 0.92877 0.9724 0.99022 0.9983 0.99926 0.99968 1.00000
biplot(precip_pca, col = c("gray", "black")) #looks like Aridity index is corelated with all the rainfall bio-clim vars, except for precip seasonality which is orthogonal, no need to use PC_socres

#make the herb/not herb factor a [1 2] numerical variable 
cr_pops$herb_fact <- ifelse(cr_pops$GrowthType == 'Herbaceous perennial', 1, 2)

#calculate other metric of population performance
pop_growth = numeric(length = dim(cr_pops)[1])
damping_ratio = numeric(length = dim(cr_pops)[1])
cv_pop_growth = numeric(length = dim(cr_pops)[1])

k = 1
while(k <= dim(cr_pops)[1]){
	matched_dat_all = used_dat_ns[which(used_dat_ns$SpeciesAuthor == cr_pops$SpeciesAuthor[k]), ] 
	matched_spp_pop = ifelse(is.na(matched_dat_all$Population), 0, matched_dat_all$Population)
	cr_pops_pop = ifelse(is.na(cr_pops$Population[k]), 0, cr_pops$Population[k])
	matched_dat_all = matched_dat_all[which(matched_spp_pop == cr_pops_pop), ]
	#split matched data becasue there are a few populatins were there is more than 1 population with the same name
	#find the chuncks of populations defined by rows with 'Mean' matrix composition
	mean_inds = c(which(matched_dat_all$MatrixComposite == 'Mean'), dim(matched_dat_all)[1] + 1)
	matched_dat = list()
	for(i in 1:(length(mean_inds) - 1)) matched_dat[[i]] = matched_dat_all[mean_inds[i]:(mean_inds[i + 1] - 1), ]

	for(i in 1:length(matched_dat)){ 
		mean_mat = compadre$mat[[matched_dat[[i]]$mat_ind[which(matched_dat[[i]]$MatrixComposite == 'Mean')]]]
		indiv_inds = matched_dat[[i]]$mat_ind[which(matched_dat[[i]]$MatrixComposite == 'Individual')]
		ind_mat  = lapply(indiv_inds, FUN = function(x) compadre$mat[[x]])
		pop_growth[k] = lambda(mean_mat$matA) 
		damping_ratio[k] = damping.ratio(mean_mat$matA)
		#need to replace NA with 0 for a few wierd pops
		if(length(ind_mat) > 0){
			for(j in 1:length(ind_mat)) ind_mat[[j]]$matA[is.na(ind_mat[[j]]$matA)] = 0
			cv_pop_growth[k] = cv(sapply(ind_mat, FUN = function(x) lambda(x$matA)))
		}else cv_pop_growth[k] = NA
		k = k + 1
	}
}

par(mfrow = c(1, 3))
hist(log(pop_growth))
hist(log(damping_ratio))
hist(log(cv_pop_growth+1))

cr_pops$lambda <- pop_growth
cr_pops$damping_ratio <- damping_ratio
cr_pops$cv_lambda <- cv_pop_growth

setwd(data_file_loc)
save(cr_pops, file = 'combined_bc_ai_all_responses.Rdata')

#Find the numner of authors for each species
auth_spp <- lapply(unique(cr_pops$SpeciesAccepted), FUN = function(x) list(spp = x, auth = auth_data$Authors[auth_data$SpeciesAccepted == x]))
auth_per_spp <- sapply(auth_spp, FUN = function(x) sum(!is.na(unique(x$auth))))
pro_one = sum(auth_per_spp == 1) / length(auth_per_spp) #0.919 spp have one author


