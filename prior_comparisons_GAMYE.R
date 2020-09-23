###################################################
# Adam C. Smith & Brandon P.M. Edwards
# GAM Paper Script
# This is short demonstration of the insensitivity of the GAM results to the prior on sigma_B

n_iter = 10000
n_thin = 10
n_burnin = 10000
n_chains = 3
n_adapt = 1000

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


models <- c("gam","gam_alt_prior")

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

n_cores <- 5
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


foreach(species = species_to_run,
        .packages = c("bbsBayes","tidybayes"),
        .inorder = FALSE,
        .errorhandling = "pass") %dopar%
  {
    
###################################################
# Analysis by Species X Model Combination
###################################################

# for (species in species_to_run)
# {
  sp.dir = paste0("output_prior_comp/", species)
  dir.create(sp.dir)
  
  
  #### identifying the K folds for cross-validation
    ## selecting stratified samples that remove 10% of data within each stratum
  jags_data <- prepare_jags_data(strat_data = stratified_data,
                                 species_to_run = species,
                                 min_max_route_years = 3,
                                 model = "slope")
  
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
                                   model = "gam",
                                   parallel = TRUE)
    
    
    
    ##################### FULL MODEL RUN ##########################
  

        #inits = function()
    model_file <- paste0("loo-models/", model, ".jags") #using hte normal-tailed error distribution
    
    jags_mod_full <- run_model(jags_data = jags_data,
                               model_file_path = model_file,
                               parameters_to_save = c("B.X","beta.x","sdX","n","na"),
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
    


#     indsa <- generate_indices(jags_mod = jags_mod_out[[2]],jags_data = jags_data,alternate_n = "n3a")
#     indso <- generate_indices(jags_mod = jags_mod_out[[1]],jags_data = jags_data,alternate_n = "n3a")
# 
#     tra <- generate_trends(indices = indsa)
#     tro <- generate_trends(indices = indso)
#     
#     plta <- plot_indices(indices_list = indsa,add_observed_means = T,add_number_routes = T,species = "Alternate Prior")
#     plto <- plot_indices(indices_list = indso,add_observed_means = T,add_number_routes = T,species = "Original Prior")
#     
#     
#   pdf(paste0(sp.dir,species,"comparative index plots.pdf"),
#       width = 11,height = 8)
#   for(j in 1:length(plta)){
#     print(plta[[j]] + plto[[j]])
#       }
# dev.off()
# 
# jags_mod_out[[1]]$summary["sdX",]
# jags_mod_out[[2]]$summary["sdX",]
# 
# YE1 = jags_mod_out[[1]]$summary[paste0("YE[",jags_data$ymin:jags_data$ymax,"]"),]
# YE2 = jags_mod_out[[2]]$summary[paste0("YE[",jags_data$ymin:jags_data$ymax,"]"),]
# 
# plot(YE1[,1],YE2[,1])
# abline(0,1)  
#  
# 
# 
# BX1 = jags_mod_out[[1]]$summary[paste0("B.X[",1:jags_data$nknots,"]"),]
# BX2 = jags_mod_out[[2]]$summary[paste0("B.X[",1:jags_data$nknots,"]"),]
# 
# plot(BX1[,1],BX2[,1])
# abline(0,1)  
# 
# 
# 
# bx1 <- gather_draws(jags_mod_out[[1]]$samples,B.X[k])
# 
# bx2 <- gather_draws(jags_mod_out[[2]]$samples,B.X[k])
# 
# sbx1 <- bx1 %>% group_by(.draw) %>% 
#   summarise(sbx = sum(.value)) %>% 
#   ungroup() %>% 
#   summarise(mean = mean(sbx),
#             lci = quantile(sbx,0.025),
#             uci = quantile(sbx,0.975))
# 
# sbx2 <- bx2 %>% group_by(.draw) %>% 
#   summarise(sbx = sum(.value)) %>% 
#   ungroup() %>% 
#   summarise(mean = mean(sbx),
#             lci = quantile(sbx,0.025),
#             uci = quantile(sbx,0.975))

    }#end species loop










