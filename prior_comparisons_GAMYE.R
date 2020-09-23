###################################################
# Adam C. Smith & Brandon P.M. Edwards
# GAM Paper Script
# This is the main computational script that fits each of the four models to each of the 10 species
# then it runs the 15-fold cross validation for each species and model
## to run as coded here requires a computer with at least 15 cores and at least 64 GB of RAM
# Created July 2019
# Last Updated April 2020
###################################################

###################################################
# Initial Setup + Setting Constants
###################################################

remove(list = ls())
K <- 15 #k = 15 fold X-valid
# These are just the defaults for bbsBayes run_model. Modify as needed
n_iter = 10000
n_thin = 10
n_burnin = 10000
n_chains = 3
n_adapt = 1000

dir.create("output", showWarnings = F)

# Install bbsBayes from Zenodo 
#Edwards, B.P.M. and A.C. Smith (2020). bbsBayes v2.1.0 (Version 2.1.0). Zenodo. 485 doi:10.5281/zenodo.3727279

library(bbsBayes)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(patchwork)
library(tidybayes)

stratified_data <- stratify(by = "bbs_usgs")


species_to_run <- c("Chestnut-collared Longspur")

models <- c("gamye","gamyealtprior")


###################################################
# Analysis by Species X Model Combination
###################################################

for (species in species_to_run)
{
  sp.dir = paste0("output_prior_comp/", species)
  dir.create(sp.dir)
  
  
  #### identifying the K folds for cross-validation
    ## selecting stratified samples that remove 10% of data within each stratum
  jags_data <- prepare_jags_data(strat_data = stratified_data,
                                 species_to_run = species,
                                 min_max_route_years = 3,
                                 model = "slope")
  
  sp.k = paste0(sp.dir, "/sp_k.RData")
  
  if(file.exists(sp.k) == FALSE){
  
    kk = vector(mode = "integer",length = jags_data$ncounts)
    for(i in 1:jags_data$nstrata){
      set.seed(2019)
      wstrat = which(jags_data$strat == i)
      
      kk[wstrat] = as.integer(ceiling(runif(length(wstrat),0,K))) 
    }
    
    ### saving the k-fold identifiers so that they're the same across all models
    save(kk,file = sp.k)
  }
    

    jags_mod_out <- vector(mode = "list",length = 2)
for(m in c(2,1)){
      
      model = models[m]
      model_dir <- paste0(sp.dir,
                        "/",
                        model)
    dir.create(model_dir)
    
    jags_data <- prepare_jags_data(strat_data = stratified_data,
                                   species_to_run = species,
                                   min_max_route_years = 3,
                                   model = "gamye",
                                   parallel = TRUE)
    
    
    
    load(sp.k)
    jags_data$ki <- kk
    
    save(jags_data,file = paste0(model_dir,"/jags_data.RData"))
    ##################### FULL MODEL RUN ##########################
  

        #inits = function()
    model_file <- paste0("loo-models/", model, ".jags") #using hte normal-tailed error distribution
    
    jags_mod_full <- run_model(jags_data = jags_data,
                               model_file_path = model_file,
                               parameters_to_save = c("B.X","beta.x","sdX","n3","n3a","YE"),
                               n_iter = n_iter,
                               #n_adapt = n_adapt,
                               n_burnin = n_burnin,
                               n_chains = n_chains,
                               n_thin = n_thin,
                               parallel = TRUE,
                               modules = NULL)
     save(list = c("jags_mod_full",
                   "jags_data"), file = paste0(model_dir, "/jags_mod_full.RData"))
    
     jags_mod_out[[m]] <- jags_mod_full
    ##################### TRENDS AND TRAJECTORIES #################
     ### this trend calculation and plotting is not necessary, only for model-checking 

    
    }#end of full model parallel loop
    


    indsa <- generate_indices(jags_mod = jags_mod_out[[2]],jags_data = jags_data,alternate_n = "n3a")
    indso <- generate_indices(jags_mod = jags_mod_out[[1]],jags_data = jags_data,alternate_n = "n3a")

    tra <- generate_trends(indices = indsa)
    tro <- generate_trends(indices = indso)
    
    plta <- plot_indices(indices_list = indsa,add_observed_means = T,add_number_routes = T,species = "Alternate Prior")
    plto <- plot_indices(indices_list = indso,add_observed_means = T,add_number_routes = T,species = "Original Prior")
    
    
  pdf(paste0(sp.dir,species,"comparative index plots.pdf"),
      width = 11,height = 8)
  for(j in 1:length(plta)){
    print(plta[[j]] + plto[[j]])
      }
dev.off()

jags_mod_out[[1]]$summary["sdX",]
jags_mod_out[[2]]$summary["sdX",]

YE1 = jags_mod_out[[1]]$summary[paste0("YE[",jags_data$ymin:jags_data$ymax,"]"),]
YE2 = jags_mod_out[[2]]$summary[paste0("YE[",jags_data$ymin:jags_data$ymax,"]"),]

plot(YE1[,1],YE2[,1])
abline(0,1)  
 


BX1 = jags_mod_out[[1]]$summary[paste0("B.X[",1:jags_data$nknots,"]"),]
BX2 = jags_mod_out[[2]]$summary[paste0("B.X[",1:jags_data$nknots,"]"),]

plot(BX1[,1],BX2[,1])
abline(0,1)  



bx1 <- gather_draws(jags_mod_out[[1]]$samples,B.X[k])

bx2 <- gather_draws(jags_mod_out[[2]]$samples,B.X[k])

sbx1 <- bx1 %>% group_by(.draw) %>% 
  summarise(sbx = sum(.value)) %>% 
  ungroup() %>% 
  summarise(mean = mean(sbx),
            lci = quantile(sbx,0.025),
            uci = quantile(sbx,0.975))

sbx2 <- bx2 %>% group_by(.draw) %>% 
  summarise(sbx = sum(.value)) %>% 
  ungroup() %>% 
  summarise(mean = mean(sbx),
            lci = quantile(sbx,0.025),
            uci = quantile(sbx,0.975))

    }#end species loop










