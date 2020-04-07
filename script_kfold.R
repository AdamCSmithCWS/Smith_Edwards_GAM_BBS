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

# Only need to run this if you don't have BBS data saved
# in directory on your computer
# yes immediately following to agree to terms and conditions
####################################
# fetch_bbs_data()
# yes

stratified_data <- stratify(by = "bbs_usgs")


species_to_run <- c("Barn Swallow",
                    "Wood Thrush",
                    "American Kestrel",
                    "Chimney Swift",
                    "Ruby-throated Hummingbird",
                    "Chestnut-collared Longspur",
                    "Cooper's Hawk",
                    "Canada Warbler",
                    "Carolina Wren",
                    "Pine Siskin")

models <- c("firstdiff", "gam", "gamye","slope")


###################################################
# Analysis by Species X Model Combination
###################################################

for (species in species_to_run)
{
  sp.dir = paste0("output/", species)
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
    
    # Set up parallel stuff
  # requires only 4 cores to run each model in parallel for species
    n_cores <- length(models)
    cluster <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cluster)
    
    
foreach(m = 1:4,
        .packages = 'bbsBayes',
        .inorder = FALSE,
        .errorhandling = "pass") %dopar%
    {
      
      model = models[m]
      model_dir <- paste0(sp.dir,
                        "/",
                        model)
    dir.create(model_dir)
    
    jags_data <- prepare_jags_data(strat_data = stratified_data,
                                   species_to_run = species,
                                   min_max_route_years = 3,
                                   model = model)
    
    
    
    load(sp.k)
    jags_data$ki <- kk
    
    save(jags_data,file = paste0(model_dir,"/jags_data.RData"))
    ##################### FULL MODEL RUN ##########################
  

        #inits = function()
    model_file <- paste0("loo-models/", model, "-t.jags") #using hte heavy-tailed error distribution
    
    jags_mod_full <- run_model(jags_data = jags_data,
                               model_file_path = model_file,
                               n_iter = n_iter,
                               #n_adapt = n_adapt,
                               n_burnin = n_burnin,
                               n_chains = n_chains,
                               n_thin = n_thin,
                               parallel = FALSE,
                               modules = NULL)
     save(jags_mod_full, file = paste0(model_dir, "/jags_mod_full.RData"))
    
    ##################### TRENDS AND TRAJECTORIES #################
     ### this trend calculation and plotting is not necessary, only for model-checking 
    dir.create(paste0(model_dir, "/plots"))
    
    # Stratum level
    strat_indices <- generate_strata_indices(jags_mod = jags_mod_full)
    strat_trends_full <- generate_strata_trends(indices = strat_indices)
    strat_trends_10yr <- generate_strata_trends(indices = strat_indices,
                                                min_year = strat_indices$y_max-10)
    s_plots <- plot_strata_indices(strat_indices)
    
    for (i in 1:length(s_plots))
    {
      png(filename = paste0(model_dir, "/plots/", names(s_plots[i]), ".png"))
      print(s_plots[[i]])
      dev.off() 
    }
    
    # National level functionality doesn't exist in bbsBayes (yet)
    
    # Continental level
    cont_indices <- generate_cont_indices(jags_mod = jags_mod_full)
    cont_trend_full <- generate_cont_trend(indices = cont_indices)
    cont_trend_10yr <- generate_cont_trend(indices = cont_indices,
                                           min_year = strat_indices$y_max-10)
    c_plot <- plot_cont_indices(indices = cont_indices)
    
    png(filename = paste0(model_dir, "/plots/1continental.png"))
    print(c_plot)
    dev.off()
    
    outtrends = rbind(cont_trend_full,cont_trend_10yr,strat_trends_full,strat_trends_10yr)
    outtrends$species = species
    write.csv(outtrends,paste0(model_dir,"/all trends ",species,".csv"))
    
    
    
    }#end of full model parallel loop
    
stopCluster(cl = cluster)














for(model in models){
    
    ##################### CROSS VALIDATION ########################
  model_dir <- paste0(sp.dir,
                      "/",
                      model)
  dir.create(paste0(model_dir, "/cv"))
    
  load(paste0(model_dir,"/jags_data.RData"))
  load(paste0(model_dir, "/jags_mod_full.RData"))
  
    inits <- get_final_values(model = jags_mod_full)
    model_file <- paste0("loo-models/", model, "-loo.jags")
    
    # Set up parallel stuff
    # requires 15 cores to run each cross-validation fold in parallel
    n_cores <- K
    cluster <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cluster)
    
    posterior <- foreach(kk = 1:K,
                         .packages = 'bbsBayes',
                         .inorder = FALSE,
                         .errorhandling = "pass",
                         .combine = 'cbind') %dopar%
      {
        
        indices_to_remove <- which(jags_data$ki == kk)
        true_count_k <- jags_data$count[indices_to_remove]
        
        n_remove <- as.integer(length(indices_to_remove))
        
        jags_data_loo <- jags_data
        jags_data_loo$count[indices_to_remove] <- NA
        jags_data_loo$I <- indices_to_remove
        jags_data_loo$Y <- true_count_k
        jags_data_loo$nRemove <- n_remove
        
        # Run model and track some LOOCV variables
        # Save the original models to disk using model_to_file() then
        # modify the models to add in the LOOCV variable tracking. 
        # Then give run_model() the path to that new model.
        jags_mod_loo <- run_model(jags_data = jags_data_loo,
                                  model_file_path = model_file,
                                  parameters_to_save = c("LambdaSubset"),
                                  track_n = FALSE,
                                  inits = inits,
                                  n_iter = n_iter,
                                  #n_adapt = n_adapt,
                                  n_burnin = 1000,
                                  n_chains = n_chains,
                                  n_thin = n_thin,
                                  parallel = FALSE,
                                  modules = NULL)
        
        # Just comment this line out if you don't want to save individual
        # model runs for each year left out
        save(jags_mod_loo, 
             file = paste0(model_dir, "/cv/k_", kk, " removed.RData"))
        
        jags_mod_loo$sims.list$LambdaSubset
      }
    
    stopCluster(cl = cluster)
    

  }

}#end species loop
