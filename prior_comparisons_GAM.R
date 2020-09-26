###################################################
# Adam C. Smith & Brandon P.M. Edwards
# GAM Paper Script
# This is simple demonstration of the insensitivity of the GAM results to the prior on sigma_B

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
    
  
    ##################### TRENDS AND TRAJECTORIES #################
     ### this trend calculation and plotting is not necessary, only for model-checking 

    
    }#end of full model parallel loop
    
  }








# plotting ----------------------------------------------------------------

library(ggplot2)
library(ggrepel)


prior_plots <- vector(mode = "list",length = length(species_to_run))
jj = 0

for(species in species_to_run[c(2,4,5,6)]){
jj = jj+1
    sp.dir = paste0("output_prior_comp/", species)
 
  jags_mod_out <- vector(mode = "list",length = 2)
  
  for(m in c(2,1)){
    
    model = models[m]
    model_dir <- paste0(sp.dir,
                        "/",
                        model)
    
    
    load(paste0(model_dir, "/jags_mod_full.RData"))
    jags_mod_out[[m]] <- jags_mod_full
    
  }
    indsa <- generate_indices(jags_mod = jags_mod_out[[2]],jags_data = jags_data,alternate_n = "na")
    indso <- generate_indices(jags_mod = jags_mod_out[[1]],jags_data = jags_data,alternate_n = "na")

    # tra <- generate_trends(indices = indsa)
    # tro <- generate_trends(indices = indso)

    # plta <- plot_indices(indices_list = indsa,add_observed_means = T,add_number_routes = T,species = "Alternate Prior")
    # plto <- plot_indices(indices_list = indso,add_observed_means = T,add_number_routes = T,species = "Original Prior")

    indsa_o <- generate_indices(jags_mod = jags_mod_out[[2]],jags_data = jags_data,alternate_n = "n")
    indso_o <- generate_indices(jags_mod = jags_mod_out[[1]],jags_data = jags_data,alternate_n = "n")
    
    int2 <- indsa_o$data_summary
    int1 <- indso_o$data_summary
    datc <- data.frame(rg = jags_data$strat_name,
                        Year = jags_data$r_year)
    plt <- vector(mode = "list",length = length(unique(indsa_o$data_summary$Region)))
    
    i = 0
    for(rg in unique(indsa_o$data_summary$Region)){
      i = i+1
      if(rg == "Continental"){
        datcc = datc
        
      }else{
        datcc = filter(datc,rg == rg)
        
      }
      int1t <- filter(int1,Region == rg)
      int1t$prior <- "gamma(0.01,0.0001)"
      int2t <- filter(int2,Region == rg)
      int2t$prior <- "gamma(2,0.2)"
      
      rtlab = "Number of Routes * 10"
      if(max(int1t$nrts) > 50){
        datcc <- slice_sample(datcc,prop = 0.1)
        rtlab = "Number of Routes * 10"
      
      if(max(int1t$nrts) > 300){
        datcc <- slice_sample(datcc,prop = 0.02)
        rtlab = "Number of Routes * 50"
      }
      }
      ip <- bind_rows(int1t,int2t)
      
      uylim <- max(ip$Index_q_0.975)
      
      cont_over = ggplot(data = ip,aes(x = Year,y = Index,group = prior,colour = prior,fill = prior))+
        theme_classic()+
        labs(title = paste(rg,species))+
        theme(legend.position = "top")+
        xlab(label = "")+
        ylab(label = "Predicted mean count")+
        geom_point(data = int1t, aes(x = Year,y = obs_mean),inherit.aes = FALSE,colour = grey(0.8),size = 0.8)+
        coord_cartesian(ylim = c(0,uylim))+
        scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
        scale_y_continuous(expand = c(0,0))+
        annotate(geom = "text",x = 1985,y = uylim*0.05,label = rtlab,colour = grey(0.6)) +
        #geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
        #geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = 5)+
        scale_colour_viridis_d(end = 0.8,begin = 0.2, aesthetics = c("colour","fill"),name = "Prior on 1/sigma2_B",option = "inferno")+
        
        geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = prior),alpha = 0.2,inherit.aes = FALSE)+
        geom_line(aes(colour = prior),size = 1.2)+
        geom_dotplot(data = datcc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)
      
      plt[[i]] <- cont_over
      
    }
    
  pdf(paste0(sp.dir,"/",species,"comparative index plots.pdf"),
      width = 11,height = 8)
  for(j in 1:length(plt)){
    print(plt[[j]])
      }
dev.off()

prior_plots[[jj]] <- plt[[length(plt)]]


# YE1 = jags_mod_out[[1]]$summary[paste0("YE[",jags_data$ymin:jags_data$ymax,"]"),]
# YE2 = jags_mod_out[[2]]$summary[paste0("YE[",jags_data$ymin:jags_data$ymax,"]"),]

# plot(YE1[,1],YE2[,1])
# abline(0,1)


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


pdf(paste0("Figures/supplement/comparative index plots.pdf"),
    width = 11,height = 8)
for(j in 1:length(prior_plots)){
  print(prior_plots[[j]])
}
dev.off()

save(list = "prior_plots",file = "figures/supplement/prior_plots.RData")







