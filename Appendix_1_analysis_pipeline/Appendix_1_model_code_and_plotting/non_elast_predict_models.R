#models predicting generic response from two distance matrices (phylo and geo) pluse matirx of neighbours response values
sink("interactions_2ac_decay_models.jags")
	cat('
		model{
			#PRIORS
			tau.global ~ dgamma(0.0001, 0.0001)#overall percison
			sd.global <- 1/sqrt(tau.global)#overall sd

			for(i in 1:9){
				b[i] ~ dnorm(0, 0.0001)#coefficents for 1 intercept, and 6 slopes, two of which are related to autocorrelation terms
			}
			#the effect of being a herb vs non-herb
			for(i in 1:num_GF){
				a[i] ~ dnorm(0, 0.0001)#coefficents for intercept of herbs and non-herbs
				for(j in 1:4){
					B[j, i] ~ dnorm(0, 0.0001) #coeefcients for herb and non-herbs for 4 predictors (life_expect, PC_temp, AI, p_sea) 
				}
			}
			#auto correlation parameters
			alpha ~ dunif(min_alpha, max_alpha) #exponent on the phlyogenetic prediction model which is a weighted average of other species, weighted by distance to most recent common ancestor, alpha controls the decay of the weight with distance 
			rho ~ dunif(min_rho, max_rho) #exponent on the geographic prediction model which is a weighted average of other population, weighted by distance to between study sites, rho controls the decay of the weight with distance  
			
			#LIKLIEHOOD 
			for(n in 1:num_data){
				#fit the linear model
				Y[n] ~ dnorm(mu[n], tau.global) 
				mu[n] <- a[GF[n]] + b[1]*phy_pred[n] + b[2]*mat_dim[n] + B[1, GF[n]]*life_expect[n] + B[2, GF[n]]*PC_temp[n] + B[3, GF[n]]*AI[n] + B[4, GF[n]]*p_sea[n] +
					b[3]*life_expect[n]*PC_temp[n] + b[4]*life_expect[n]*AI[n] + b[5]*life_expect[n]*p_sea[n] + b[6]*PC_temp[n]*AI[n] + b[7]*PC_temp[n]*p_sea[n] + b[8]*AI[n]*p_sea[n] + b[9]*geo_pred[n] 

				#fit the phylo and geo models
				for(k in 1:(num_data - 1)){
					phy_weight[n, k] <- exp(-alpha*phy_dist[n, k])
					geo_weight[n, k] <- exp(-rho*geo_dist[n, k]) 
				}
				#make the sum of weights in denominator > 0
				phy_pred[n] <- sum(Y_neigh_phy[n, 1:phy_max_ind[n]]*phy_weight[n, 1:phy_max_ind[n]]) / max(0.00000000000001, sum(phy_weight[n, 1:phy_max_ind[n]]))
				geo_pred[n] <- sum(Y_neigh_geo[n, 1:geo_max_ind[n]]*geo_weight[n, 1:geo_max_ind[n]]) / max(0.00000000000001, sum(geo_weight[n, 1:geo_max_ind[n]]))
			}

			#MODEL PERFORMANCE	
			for(n in 1:num_data){
				#model performance predictors to make sure the model is not compleatly rubbish
				residual[n] <- Y[n] - mu[n]#takes the residual out of the eq below so can be tracked
				sq[n] <- pow(residual[n], 2)#residual squared 
				Y.new[n] ~ dnorm(mu[n], tau.global)#new data
				sq.new[n] <- pow(Y.new[n] - mu[n], 2)#residual squared for the new data
			}
			fit <- sum(sq[])#sum of squares for data
			fit.new <- sum(sq.new[])#sum of squares for simulated data
			test <- step(fit.new - fit)
			bpvalue <- mean(test)
		}',fill = TRUE)
sink()

#REDUCED INTERACTIONS BECASUE SOME MODELS HAVE A HARD TIME CONVERGING WITH THE FULL SET 
sink("redu_int_2ac_decay_models.jags")
	cat('
		model{
			#PRIORS
			tau.global ~ dgamma(0.0001, 0.0001)#overall percison
			sd.global <- 1/sqrt(tau.global)#overall sd

			for(i in 1:13){
				b[i] ~ dnorm(0, 0.0001)#coefficents for 1 intercept, and 6 slopes, two of which are related to autocorrelation terms
			}
			#the effect of being a herb vs non-herb
			for(i in 1:num_GF){
				a[i] ~ dnorm(0, 0.0001)#coefficents for intercept of herbs and non-herbs
			}
			#auto correlation parameters
			alpha ~ dunif(min_alpha, max_alpha) #exponent on the phlyogenetic prediction model which is a weighted average of other species, weighted by distance to most recent common ancestor, alpha controls the decay of the weight with distance 
			rho ~ dunif(min_rho, max_rho) #exponent on the geographic prediction model which is a weighted average of other population, weighted by distance to between study sites, rho controls the decay of the weight with distance  
			
			#LIKLIEHOOD 
			for(n in 1:num_data){
				#fit the linear model
				Y[n] ~ dnorm(mu[n], tau.global) 
				mu[n] <- a[GF[n]] + b[1]*phy_pred[n] + b[2]*mat_dim[n] + b[3]*life_expect[n] + b[4]*PC_temp[n] + b[5]*AI[n] + b[6]*p_sea[n] +
					b[7]*life_expect[n]*PC_temp[n] + b[8]*life_expect[n]*AI[n] + b[9]*life_expect[n]*p_sea[n] + b[10]*PC_temp[n]*AI[n] + b[11]*PC_temp[n]*p_sea[n] + b[12]*AI[n]*p_sea[n] + b[13]*geo_pred[n] 

				#fit the phylo and geo models
				for(k in 1:(num_data - 1)){
					phy_weight[n, k] <- exp(-alpha*phy_dist[n, k])
					geo_weight[n, k] <- exp(-rho*geo_dist[n, k]) 
				}
				#make the sum of weights in denominator > 0
				phy_pred[n] <- sum(Y_neigh_phy[n, 1:phy_max_ind[n]]*phy_weight[n, 1:phy_max_ind[n]]) / max(0.00000000000001, sum(phy_weight[n, 1:phy_max_ind[n]]))
				geo_pred[n] <- sum(Y_neigh_geo[n, 1:geo_max_ind[n]]*geo_weight[n, 1:geo_max_ind[n]]) / max(0.00000000000001, sum(geo_weight[n, 1:geo_max_ind[n]]))
			}

			#MODEL PERFORMANCE	
			for(n in 1:num_data){
				#model performance predictors to make sure the model is not compleatly rubbish
				residual[n] <- Y[n] - mu[n]#takes the residual out of the eq below so can be tracked
				sq[n] <- pow(residual[n], 2)#residual squared 
				Y.new[n] ~ dnorm(mu[n], tau.global)#new data
				sq.new[n] <- pow(Y.new[n] - mu[n], 2)#residual squared for the new data
			}
			fit <- sum(sq[])#sum of squares for data
			fit.new <- sum(sq.new[])#sum of squares for simulated data
			test <- step(fit.new - fit)
			bpvalue <- mean(test)
		}',fill = TRUE)
sink()

#NO INTERACTION NO SPECIES 
sink("noint_2ac_decay_models.jags")
	cat('
		model{
			#PRIORS
			tau.global ~ dgamma(0.0001, 0.0001)#overall percison
			sd.global <- 1/sqrt(tau.global)#overall sd

			for(i in 1:7){
				b[i] ~ dnorm(0, 0.0001)#coefficents for 1 intercept, and 6 slopes, two of which are related to autocorrelation terms
			}
			#the effect of being a herb vs non-herb
			for(i in 1:num_GF){
				a[i] ~ dnorm(0, 0.0001)#coefficents for intercept of herbs and non-herbs
			}
			#auto correlation parameters
			alpha ~ dunif(min_alpha, max_alpha) #exponent on the phlyogenetic prediction model which is a weighted average of other species, weighted by distance to most recent common ancestor, alpha controls the decay of the weight with distance 
			rho ~ dunif(min_rho, max_rho) #exponent on the geographic prediction model which is a weighted average of other population, weighted by distance to between study sites, rho controls the decay of the weight with distance  
			
			#LIKLIEHOOD 
			for(n in 1:num_data){
				#fit the linear model
				Y[n] ~ dnorm(mu[n], tau.global) 
				mu[n] <- a[GF[n]] + b[1]*phy_pred[n] + b[2]*mat_dim[n] + b[3]*life_expect[n] + b[4]*PC_temp[n] + b[5]*AI[n] + b[6]*p_sea[n] + b[7]*geo_pred[n] 

				#fit the phylo and geo models
				for(k in 1:(num_data - 1)){
					phy_weight[n, k] <- exp(-alpha*phy_dist[n, k])
					geo_weight[n, k] <- exp(-rho*geo_dist[n, k]) 
				}
				#make the sum of weights in denominator > 0
				phy_pred[n] <- sum(Y_neigh_phy[n, 1:phy_max_ind[n]]*phy_weight[n, 1:phy_max_ind[n]]) / max(0.00000000000001, sum(phy_weight[n, 1:phy_max_ind[n]]))
				geo_pred[n] <- sum(Y_neigh_geo[n, 1:geo_max_ind[n]]*geo_weight[n, 1:geo_max_ind[n]]) / max(0.00000000000001, sum(geo_weight[n, 1:geo_max_ind[n]]))
		
			}

			#MODEL PERFORMANCE	
			for(n in 1:num_data){
				#model performance predictors to make sure the model is not compleatly rubbish
				residual[n] <- Y[n] - mu[n]#takes the residual out of the eq below so can be tracked
				sq[n] <- pow(residual[n], 2)#residual squared 
				Y.new[n] ~ dnorm(mu[n], tau.global)#new data
				sq.new[n] <- pow(Y.new[n] - mu[n], 2)#residual squared for the new data
			}
			fit <- sum(sq[])#sum of squares for data
			fit.new <- sum(sq.new[])#sum of squares for simulated data
			test <- step(fit.new - fit)
			bpvalue <- mean(test)
		}',fill = TRUE)
sink()


#NO fixed effects POOLED GROWTH TYPE NO SPECIES
sink("nofixed_2ac_decay_model.jags")
	cat('
		model{
			#PRIORS
			tau.global ~ dgamma(0.0001, 0.0001)#overall percison
			sd.global <- 1/sqrt(tau.global)#overall sd

			a ~ dnorm(0, 0.0001)
			for(i in 1:2){
				b[i] ~ dnorm(0, 0.0001)
			}
			
			#auto correlation parameters
			alpha ~ dunif(min_alpha, max_alpha) #exponent on the phlyogenetic prediction model which is a weighted average of other species, weighted by distance to most recent common ancestor, alpha controls the decay of the weight with distance 
			rho ~ dunif(min_rho, max_rho) #exponent on the geographic prediction model which is a weighted average of other population, weighted by distance to between study sites, rho controls the decay of the weight with distance  
			
			#LIKLIEHOOD 
			for(n in 1:num_data){
				#fit the linear model
				Y[n] ~ dnorm(mu[n], tau.global) 
				mu[n] <- a + b[1]*phy_pred[n] + b[2]*geo_pred[n] 

				#fit the phylo and geo models
				for(k in 1:(num_data - 1)){
					phy_weight[n, k] <- exp(-alpha*phy_dist[n, k])
					geo_weight[n, k] <- exp(-rho*geo_dist[n, k]) 
				}
				#make the sum of weights in denominator > 0
				phy_pred[n] <- sum(Y_neigh_phy[n, 1:phy_max_ind[n]]*phy_weight[n, 1:phy_max_ind[n]]) / max(0.00000000000001, sum(phy_weight[n, 1:phy_max_ind[n]]))
				geo_pred[n] <- sum(Y_neigh_geo[n, 1:geo_max_ind[n]]*geo_weight[n, 1:geo_max_ind[n]]) / max(0.00000000000001, sum(geo_weight[n, 1:geo_max_ind[n]]))
		
			}

			#MODEL PERFORMANCE	
			for(n in 1:num_data){
				#model performance predictors to make sure the model is not compleatly rubbish
				residual[n] <- Y[n] - mu[n]#takes the residual out of the eq below so can be tracked
				sq[n] <- pow(residual[n], 2)#residual squared 
				Y.new[n] ~ dnorm(mu[n], tau.global)#new data
				sq.new[n] <- pow(Y.new[n] - mu[n], 2)#residual squared for the new data
			}
			fit <- sum(sq[])#sum of squares for data
			fit.new <- sum(sq.new[])#sum of squares for simulated data
			test <- step(fit.new - fit)
			bpvalue <- mean(test)
		}',fill = TRUE)
sink()

#NO Phy prediction POOLED GROWTH TYPE NO SPECIES
sink("nofixed_1ac_decay_model.jags")
	cat('
		model{
			#PRIORS
			tau.global ~ dgamma(0.0001, 0.0001)#overall percison
			sd.global <- 1/sqrt(tau.global)#overall sd

			a ~ dnorm(0, 0.0001)
			b ~ dnorm(0, 0.0001)
						
			#auto correlation parameters
			rho ~ dunif(min_rho, max_rho) #exponent on the geographic prediction model which is a weighted average of other population, weighted by distance to between study sites, rho controls the decay of the weight with distance  
			
			#LIKLIEHOOD 
			for(n in 1:num_data){
				#fit the linear model
				Y[n] ~ dnorm(mu[n], tau.global) 
				mu[n] <- a + b*ac_pred[n] 

				#fit the phylo and geo models
				for(k in 1:(num_data - 1)){
					ac_weight[n, k] <- exp(-rho*ac_dist[n, k]) 
				}
				#make the sum of weights in denominator > 0
				ac_pred[n] <- sum(Y_neigh[n, 1:ac_max_ind[n]]*ac_weight[n, 1:ac_max_ind[n]]) / max(0.00000000000001, sum(ac_weight[n, 1:ac_max_ind[n]]))
			}

			#MODEL PERFORMANCE	
			for(n in 1:num_data){
				#model performance predictors to make sure the model is not compleatly rubbish
				residual[n] <- Y[n] - mu[n]#takes the residual out of the eq below so can be tracked
				sq[n] <- pow(residual[n], 2)#residual squared 
				Y.new[n] ~ dnorm(mu[n], tau.global)#new data
				sq.new[n] <- pow(Y.new[n] - mu[n], 2)#residual squared for the new data
			}
			fit <- sum(sq[])#sum of squares for data
			fit.new <- sum(sq.new[])#sum of squares for simulated data
			test <- step(fit.new - fit)
			bpvalue <- mean(test)
		}',fill = TRUE)
sink()

#functions to get important information needed to run the model so I can rember what they are with out having to open up this file
get_model_names = function(){
  return(c('interactions_2ac_decay_models.jags', 'noint_2ac_decay_models.jags', 'nofixed_2ac_decay_model.jags', 'nofixed_1ac_decay_model.jags', 'redu_int_2ac_decay_models.jags'))
}
get_model_input = function(model){
  if(grepl('nofixed_2ac', model)) return(c('Y', 'geo_dist', 'phy_dist', 'Y_neigh_geo', 'Y_neigh_phy', 'num_data', 'geo_max_ind', 'phy_max_ind', 'min_rho', 'max_rho', 'min_alpha', 'max_alpha'))
  if(grepl('nofixed_1ac', model)) return(c('Y', 'ac_dist', 'Y_neigh', 'num_data', 'ac_max_ind', 'min_rho', 'max_rho'))
  if(grepl('int', model)) return(c('Y', 'geo_dist', 'phy_dist', 'Y_neigh_geo', 'Y_neigh_phy', 'num_data', 'geo_max_ind', 'phy_max_ind', 'GF', 'num_GF',
    'mat_dim[n] + life_expect[n] + PC_temp[n] + AI[n] + p_sea[n]', 'min_rho', 'max_rho', 'min_alpha', 'max_alpha'))
}
get_model_parameters = function(model){
  if(grepl('nofixed_2ac', model)) return(c('a = [1]:dnorm', 'b = [1:2]:dnorm', 'rho = [1]:dunif(0,1)', 'alpha = [1]:dunif(0,1)', 'tau.global = [1]:dgamma(0,inf)'))
  if(grepl('redu_int_2ac', model)) return(c('a = [1:num_GF]:dnorm', 'b = [1:13]:dnorm', 'rho = [1]:dunif(0,1)', 'alpha = [1]:dunif(0,1)', 'tau.global = [1]:dgamma(0,inf)'))
  if(grepl('nofixed_1ac', model)) return(c('a = [1]:dnorm', 'b = [1]:dnorm', 'rho = [1]:dunif(0,1)', 'tau.global = [1]:dgamma(0,inf)'))
  if(grepl('noint', model)) return(c('a = [1:num_GF]:dnorm', 'b = [1:7]:dnorm', 'rho = [1]:dunif(0,1)', 'alpha = [1]:dunif(0,1)', 'tau.global = [1]:dgamma(0,inf)'))
  if(grepl('interactions', model)) return(c('a = [1:num_GF]:dnorm', 'b = [1:9]:dnorm', 'B = [1:4, num_GF]:dnorm',  'rho = [1]:dunif(0,1)', 'alpha = [1]:dunif(0,1)', 'tau.global = [1]:dgamma(0,inf)'))
}


print('use get_model_names() to see names of models')
print('use get_model_input(model) to see input data for the model named')
print('use get_model_parameters(model) to see parameter format for the named model, required for initalisation and monitoring')
