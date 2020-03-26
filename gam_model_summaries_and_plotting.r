
# These are just the defaults for bbsBayes run_model. Modify as needed
n_iter = 10000
n_thin = 10
n_burnin = 0
n_chains = 3
n_adapt = NULL

library(bbsBayes)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)


stratified_data <- stratify(by = "bbs_usgs")

species_to_run <- c("Barn Swallow",
                    "Wood Thrush",
                    "American Kestrel",
                    "Chimney Swift",
                    "Ruby-throated Hummingbird",
                    "Chestnut-collared Longspur",
                    "Cooper's Hawk",
                    "Wild Turkey",
                    "Horned Lark")


models <- c("gam", "gamye")
#models <- c("gamye")

for(species in species_to_run[6:7]){

  sp.dir = paste0("output/", species)
  
  n_cores <- length(models)
  cluster <- makeCluster(n_cores, type = "PSOCK")
  registerDoParallel(cluster)
  
  
  foreach(m = 1:2,
          .packages = 'bbsBayes',
          .inorder = FALSE,
          .errorhandling = "pass") %dopar%
    {
      
     model = models[m]
     
model_dir <- paste0(sp.dir,
                    "/",
                    model)

if(model == "gam")
  {
  parms = c("B.X","beta.X","taubeta","tauX")
}else{
  parms = c("n3","B.X","beta.X","taubeta","tauX")
}


load(paste0(model_dir,"/jags_data.RData"))
load(paste0(model_dir, "/jags_mod_full.RData"))

inits <- get_final_values(model = jags_mod_full)


model_file <- paste0("loo-models/", model, "-t.jags") #using hte heavy-tailed error distribution

jags_mod_param <- run_model(jags_data = jags_data,
                           model_file_path = model_file,
                           n_iter = n_iter,
                           inits = inits,
                           #n_adapt = n_adapt,
                           n_burnin = n_burnin,
                           n_chains = n_chains,
                           n_thin = n_thin,
                           parallel = TRUE,
                           modules = NULL,
                           parameters_to_save = parms,
                           track_n = FALSE)

save(jags_mod_param, file = paste0(model_dir, "/parameter_model_run.RData"))


}#models
  
  
  stopCluster(cl = cluster)
  
  
}#species













