###################################################
# Adam C. Smith & Brandon P.M. Edwards
# GAM Paper Script
# gam-kfold-bbsBayes.R
# Created July 2019
# Last Updated July 2019
###################################################

###################################################
# Initial Setup + Setting Constants
###################################################

remove(list = ls())

# These are just the defaults for bbsBayes run_model. Modify as needed
n_iter = 10000
n_thin = 10
n_burnin = 20000
n_chains = 3
n_adapt = 500

dir.create("output", showWarnings = F)

# Install v1.1.2 from Github (comment these three lines out if already installed:
install.packages("devtools")
library(devtools)
devtools::install_github("BrandonEdwards/bbsBayes")#, ref = "v1.1.2")

library(bbsBayes)
library(foreach)
library(doParallel)

# Only need to run this if you don't have BBS data saved
# in directory on your computer
# yes immediately following to agree to terms and conditions
fetch_bbs_data()
yes

stratified_data <- stratify(by = "bbs_cws")

species_to_run <- c("Barn Swallow",
                    "Tree Swallow",
                    "Cliff Swallow",
                    "Northern Rough-winged Swallow",
                    "Bank Swallow",
                    "Purple Martin",
                    "Chimney Swift",
                    "Eastern Whip-poor-will")

models <- c("slope", "firstdiff", "gam", "gamye")

###################################################
# Analysis by Species X Model Combination
###################################################

for (species in species_to_run)
{
  dir.create(paste0("output/", species))
  for (model in models)
  {
    model_dir <- paste0("output/",
                        species,
                        "/",
                        model)
    dir.create(model_dir)
    
    ##################### FULL MODEL RUN ##########################
    jags_data <- prepare_jags_data(strat_data = stratified_data,
                                   species_to_run = species,
                                   model = model)
    
    jags_mod_full <- run_model(jags_data = jags_data,
                               n_iter = n_iter,
                               n_adapt = n_adapt,
                               n_burnin = n_burnin,
                               n_chains = n_chains,
                               n_thin = n_thin,
                               parallel = TRUE)
    save(jags_mod_full, file = paste0(model_dir, "/jags_mod_full.RData"))
    
    ##################### TRENDS AND TRAJECTORIES #################
    dir.create(paste0(model_dir, "/plots"))
    
    # Stratum level
    strat_indices <- generate_strata_indices(jags_mod = jags_mod_full)
    strat_trends_full <- generate_strata_trends(indices = strat_indices)
    strat_trends_10yr <- generate_strata_trends(indices = strat_indices,
                                                min_year = 2008)
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
                                           min_year = 2008)
    c_plot <- plot_cont_indices(indices = cont_indices)
    
    png(filename = paste0(model_dir, "/plots/1continental.png"))
    print(c_plot)
    dev.off()
    
    ##################### CROSS VALIDATION ########################
    dir.create(paste0(model_dir, "/cv"))
    
    inits <- get_final_values(model = jags_mod_full)
    model_file <- paste0("loo-models/", model, "-loo.jags")
    
    # Set up parallel stuff
    n_cores <- 3
    cluster <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cluster)

    posterior <- foreach(year = jags_data$ymin:jags_data$ymax,
                         .packages = 'bbsBayes',
                         .combine = 'cbind') %dopar%
    {
      indices_to_remove <- which(jags_data$year == year)
      true_count_yr <- jags_data$count[indices_to_remove]
      
      n_remove <- as.integer(length(indices_to_remove))
      
      jags_data_loo <- jags_data
      jags_data_loo$count[indices_to_remove] <- NA
      jags_data_loo$I <- indices_to_remove
      jags_data_loo$Y <- true_count_yr
      jags_data_loo$nRemove <- n_remove
      
      # Run model and track some LOOCV variables
      # Save the original models to disk using model_to_file() then
      # modify the models to add in the LOOCV variable tracking. 
      # Then give run_model() the path to that new model.
      jags_mod_loo <- run_model(jags_data = jags_data_loo,
                                model_file_path = model_file,
                                parameters_to_save = c("logprob",
                                                       "LambdaSubset"),
                                track_n = FALSE,
                                inits = inits,
                                n_iter = n_iter,
                                n_adapt = n_adapt,
                                n_burnin = 0,
                                n_chains = n_chains,
                                n_thin = n_thin,
                                parallel = FALSE)
      
      # Just comment this line out if you don't want to save individual
      # model runs for each year left out
      save(jags_mod_loo, 
           file = paste0(model_dir, "/cv/year", year, "removed.RData"))
      
      jags_mod_loo$sims.list$LambdaSubset
    }
    
    stopCluster(cl = cluster)
    
    true_count <- jags_data$count
    true_index <- NULL
    for (year in jags_data$ymin:jags_data$ymax)
    {
      true_index <- c(true_index, which(jags_data$year == year))
    }
    
    loocv <- numeric(length = length(true_count))
    for (i in 1:length(true_count))
    {
      loocv[i] <- log(mean(dpois(true_count[true_index[i]], posterior[,i])))
    }
    
    save(loocv, file = paste0(model_dir, "/cv/loocv.RData"))
  }
}