### comparing the results of the k-fold cross validation
library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)


models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

species_to_run = c("Horned Lark","Wood Thrush", "American Kestrel","Barn Swallow","Chestnut-collared Longspur","Cooper's Hawk","Ruby-throated Hummingbird")

for(species in species_to_run[1]){

sp_dir = paste0("output/",species,"/")
#### calculate all annual indices (strata and continental)
#### and compile into a single data.frame
all_data = list()
length(all_data) = length(models)
names(all_data) = models


all_inds = list()
length(all_inds) = length(models)
names(all_inds) = models


### colour pallette

source("colourblind safe qualitative pallete.r")
model_pallete <- safe.pallet[[4]] 
model_pallete <- model_pallete[c(2,1,3,4)]

names(model_pallete) <- models

### ggplot2::scale_colour_manual(values = map_palette, aesthetics = c("colour","fill"))+

K = 15
n.iter = 3000 


for(m in models){
  m_dir = paste0(sp_dir,m,"/")
  m_dir_ext = paste0("D:/GAM Paper Script/",m_dir)
  
    load(paste0(m_dir,"jags_data.RData"))

    load(paste0(m_dir,"jags_mod_full.RData"))

  tinds = generate_regional_indices(jags_data = jags_data,
                                    jags_mod = jags_mod_full,
                                    max_backcast = NULL)
  
  ttrends = generate_regional_trends(indices = tinds,slope = T)
  
  tinds$data_summary$species = species
  tinds$data_summary$model = m
  
  
  ttrends$species = species
  ttrends$model = m
  
  
  if(m == models[1]){
    indsout = tinds$data_summary
    trendsout = ttrends
  }else{
    indsout = rbind(indsout,tinds$data_summary)
    trendsout = rbind(trendsout,ttrends)
    
  }
  
  all_data[[m]] = jags_data
  all_inds[[m]] = tinds
  
  rm(list = c("jags_mod_full"))
  
  
  
  ###### calculating the point-wise log probability for removed counts in each of k cross-validations
  
  t1 = Sys.time()
  true_count <- jags_data$count
  ki <- jags_data$ki
  
  
  loo <- matrix(NA,nrow = n.iter,ncol = length(true_count))
  
  
  for(kk in 1:K){
    
    
 true_index <- which(ki == kk)

   load(file = paste0(m_dir_ext, "cv/k_", kk, " removed.RData"))
    
    lambda.posterior = jags_mod_loo$sims.list$LambdaSubset
    
    for(i in 1:length(true_index)){
    loo[,true_index[i]] = dpois(true_count[true_index[i]], lambda.posterior[,i],log = F)
  }
    
  }
  
  #save(list = c("loo"),file = paste0(m_dir,"loo.RData"))
  t2 = Sys.time()
  t2-t1
  
  
  dat.df = get_prepared_data(jags_data = jags_data)
  dat.df$ki = jags_data$ki
  
  dat.df[,"mean.loo"] <- log(apply(loo,MARGIN = 2,FUN = mean)) #the point-wise mean of the posterior distributions of the log probability of the left-out counts given the model and the parameter estimates
  #dat.df[,"sd.loo"] <- apply(loo,MARGIN = 2,FUN = sd) #the point-wise sd of the posterior distributions of the log probability of the left-out counts given the model and the parameter estimates
  #dat.df[,"prec.loo"] <- 1/(dat.df[,"sd.loo"]^2) #the point-wise precision of the posterior distributions of the log probability of the left-out counts given the model and the parameter estimates
  
  for(q in c(0.5,0.025,0.975)){
    dat.df[,paste0("q",q,".loo")] <- log(apply(loo,MARGIN = 2,FUN = quantile,probs = q))
    
  } #calculates the quantiles of the posterior distributions of the log probability of the left-out counts given the model and the parameter estimates
  
  dat.df$model = m
  
  write.csv(dat.df,paste0(m_dir," log point-wise posterior prob.csv"))
  
  
  if(m == models[1]){
    alldat = dat.df
    datt = get_prepared_data(jags_data = jags_data)
    write.csv(datt,paste0(sp_dir,"original data file.csv"))
  }else{
    alldat = rbind(alldat,dat.df)
  }

  
  
} ### end indices and trends calculations

alldat$unit = factor(paste(alldat$Stratum,alldat$Route,alldat$Year,sep = "_"))

# x11()
# plot(alldat$prec.loo,1/((alldat$q0.975.loo - alldat$q0.025.loo)/(1.96*2))^2)
# plot(alldat$q0.5.loo,alldat$mean.loo)
# plot(alldat$q0.5.loo,((alldat$q0.975.loo - alldat$q0.025.loo)/(1.96*2))^2)

# abline(0,1)
#above suggests that the raw sd calculation is an overestimate of the error, and prone to some extreme values, likely because the point-wise loo stats have some very large tails
# replace with an alternative measure of precision that should be less sensitive to the tails

# alldat$prec.loo <- 1/((alldat$q0.975.loo - alldat$q0.025.loo)/(1.96*2))^2


write.csv(alldat,paste0(sp_dir," all models point-wise log prob.csv"))





tosave = list(model_pallete = model_pallete,
              all_inds = all_inds,
              all_data = all_data,
              species = species,
              models = models,
              indsout = indsout,
              trendsout = trendsout,
              alldat = alldat,
              datt = datt)


indcont = indsout[which(indsout$Region_type == "continental"),]
indcont2 = indcont[which(indcont$model == "slope"),]

uylim = max(c(indcont$Index_q_0.975,indcont$obs_mean))
indcont2$prts = (indcont2$nrts/indcont2$nrts_total)*uylim



labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
labl_obs$label = "Observed mean counts"
cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = model))+
  theme_classic()+
  labs(title = paste(species,"Continental"))+
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
  geom_line(aes(colour = model),size = 2)+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.7))+
  coord_cartesian(ylim = c(0,uylim))+
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim)+
  #annotate(geom = "text",x = labl_obs$Year,y = labl_obs$obs_mean,label = "Observed mean counts")+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
  geom_col(data = indcont2,aes(x = Year,y = prts),width = 0.2,inherit.aes = F,fill = "darkorange",alpha = 0.2)+
  #geom_dotplot(data = datt,mapping = aes(x = Year),binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = "darkorange",alpha = 0.2,dotsize = 0.4)+
  annotate(geom = "text",x = 1990,y = 0.02*uylim,label = paste("total of",max(indcont2$nrts_total),"routes"),colour = "darkorange",alpha = 0.4)

tosave = c(tosave,
           list(cont_over = cont_over))

#print(cont_over)

### facet plot of the trajectories overlaid

indstrata = indsout[which(indsout$Region_type != "continental"),]
indstrat1 = indsout[which(indsout$Region_type != "continental" & indsout$model == "gamye"),]

uylim = max(c(indstrata$Index_q_0.975,indstrata$obs_mean))
labl_obs = unique(indstrata[which(indstrata$Year == 1970),c("Year","obs_mean","Region")])
labl_obs$label = "Observed mean counts"

nreg = length((unique(indstrata$Region)))

pdf(paste0(sp_dir,"overplot facets.pdf"),
    height = 8.5,
    width = 11)
print(cont_over)
for(pp in 1:ceiling(nreg/12)){
strat_over = ggplot(data = indstrata,aes(x = Year,y = Index,group = model))+
  theme_classic()+
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.1)+
  geom_line(aes(colour = model),size = 1)+
  geom_point(data = indstrat1,aes(x = Year,y = obs_mean),colour = grey(0.7),size = 0.5)+
  #geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim)+
  #annotate(geom = "text",x = labl_obs$Year,y = labl_obs$obs_mean,label = "Observed mean counts")+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
  facet_wrap_paginate(facets = ~Region,nrow = 3,ncol = 4,page = pp,scales = "free")
print(strat_over)
}
dev.off()

pdf(paste0(sp_dir,"overplot by strat.pdf"),
    height = 8.5,
    width = 11)

print(cont_over)

for(pp in unique(indstrata$Region)){
  indstrat = indsout[which(indsout$Region == pp),]
  indstrat1 = indsout[which(indsout$Region == pp & indsout$model == "gamye"),]
  
  uylim = max(c(indstrat$Index_q_0.975,indstrat$obs_mean))
  labl_obs = unique(indstrat[which(indstrat$Year == 1970),c("Year","obs_mean","Region")])
  labl_obs$label = "Observed mean counts"
  
  uylim = max(c(indstrat$Index_q_0.975,indstrat$obs_mean))
  indstrat1$prts = (indstrat1$nrts/indstrat1$nrts_total)*uylim
  
  datt1 = datt[which(datt$Stratum == pp),] 
  
  labl_obs = unique(indstrat[which(indstrat$Year == 1970),c("Year","obs_mean")])
  labl_obs$label = "Observed mean counts"
  strat_over = ggplot(data = indstrat,aes(x = Year,y = Index,group = model))+
    theme_classic()+
    labs(title = paste(species,pp))+
    geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
    geom_line(aes(colour = model),size = 2)+
    geom_point(aes(x = Year,y = obs_mean),colour = grey(0.7))+
    coord_cartesian(ylim = c(0,uylim))+
    geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim)+
    #annotate(geom = "text",x = labl_obs$Year,y = labl_obs$obs_mean,label = "Observed mean counts")+
    scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
    geom_dotplot(data = datt1,mapping = aes(x = Year),binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = "darkorange",alpha = 0.2,dotsize = 0.4)+
    annotate(geom = "text",x = 1990,y = -0.02*uylim,label = paste("total of",max(indstrat1$nrts_total),"routes"),colour = "darkorange",alpha = 0.4)
 
   print(strat_over)
}
dev.off()


##### load the cross-validation estimates


#### overall fit comparisons.

#alldat$prec.loo = 1/((alldat$q0.0975.loo - alldat$q0.025.loo)/(1.96*2))^2

save(list = c("tosave"),file = paste0(sp_dir,"saved objects.RData"))
torm = c(names(tosave)[-which(names(tosave) %in% c("species","models"))],"jags_data","loo","lambda.posterior","dat.df","indcont","indcont2","indstrat")
rm(list = c("tosave",torm))

}#species loop
























library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)



models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

# species_to_run = c("Wood Thrush", "American Kestrel","Barn Swallow","Chestnut-collared Longspur","Cooper's Hawk","Ruby-throated Hummingbird")


for(species in species_to_run[1]){
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
############ Bayesian model estimating the difference in fit among models while accounting for the uncertainty in the point-wise loo
alldat = tosave$alldat



######## standard Z-score pairwise comparisons



sum.loo <- alldat %>% group_by(model) %>% summarise(sum = sum(mean.loo), mean = mean(mean.loo),sd = sd(mean.loo))

sum.loo.y <- alldat %>% group_by(model,Year) %>% summarise(sum = sum(mean.loo), mean = mean(mean.loo))


loo.point <- alldat %>% select(.,Year:ki,model,mean.loo) %>%
  pivot_wider(names_from = model,values_from = c(mean.loo),values_fn = list(mean.loo = mean))


# alldat$modelprec = paste(alldat$model,"prec",sep = "_")
# prec.loo.point <- alldat %>% select(.,Year:ki,modelprec,prec.loo) %>%
#   pivot_wider(names_from = modelprec,values_from = c(prec.loo),values_fn = list(prec.loo = mean))
# 
# 
# loo.point <- left_join(med.loo.point,prec.loo.point)

### the mean function is required because (apparrently) there are two route-year combinations that are repeated in the dataset

for(i in 1:nrow(loo.point)){
  loo.point[i,"best"] <- models[which.max(loo.point[i,models])]
}

# contr_names = c(paste(models[1],models[2],sep = "_"),
#                 paste(models[1],models[3],sep = "_"),
#                 paste(models[1],models[4],sep = "_"),
#                 paste(models[2],models[3],sep = "_"),
#                 paste(models[2],models[4],sep = "_"),
#                 paste(models[3],models[4],sep = "_")
# )
# 
# contrast_full_names = gsub(gsub(toupper(contr_names),pattern = "FIRSTDIFF",replacement = "DIFFERENCE",fixed = T),pattern = "_",replacement = " vs ",fixed = T)
# names(contrast_full_names) = contr_names


loo.point[,"gamye_gam"] <- loo.point[,"gamye"] -loo.point[,"gam"]

loo.point[,"gamye_firstdiff"] <- loo.point[,"gamye"] -loo.point[,"firstdiff"]

loo.point[,"gamye_slope"] <- loo.point[,"gamye"] -loo.point[,"slope"]

loo.point[,"gam_firstdiff"] <- loo.point[,"gam"] -loo.point[,"firstdiff"]

loo.point[,"gam_slope"] <- loo.point[,"gam"] -loo.point[,"slope"]

loo.point[,"firstdiff_slope"] <- loo.point[,"firstdiff"] -loo.point[,"slope"]

write.csv(loo.point,paste0(sp_dir,"wide form lppd.csv"))



qq = ggplot(data = loo.point,aes(sample = gamye_firstdiff))+
  geom_qq()+
  geom_qq_line()

pdf(file = paste0(sp_dir,"qq plot gamye_firstdiff.pdf"))
print(qq)
dev.off()

qq = ggplot(data = loo.point,aes(sample = gamye_slope))+
  geom_qq()+
  geom_qq_line()

pdf(file = paste0(sp_dir,"qq plot gamye_slope.pdf"))
print(qq)
dev.off()


qq = ggplot(data = loo.point,aes(sample = gamye_firstdiff))+
  geom_qq(distribution = stats::qt,
          dparams = list(df = 1.5))+
  geom_qq_line(distribution = stats::qt,
               dparams = list(df = 1.5))

pdf(file = paste0(sp_dir,"qq plot t-df-1.5 gamye_firstdiff.pdf"))
print(qq)
dev.off()

qq = ggplot(data = loo.point,aes(sample = gamye_slope))+
  geom_qq(distribution = stats::qt,
          dparams = list(df = 1.5))+
  geom_qq_line(distribution = stats::qt,
               dparams = list(df = 1.5))

pdf(file = paste0(sp_dir,"qq plot t-df-1.5 gamye_slope.pdf"))
print(qq)
dev.off()


}#species


# 
# hist(loo.point$gamye_firstdiff)
# #########
# 
# 
#  


# library(bbsBayes)
# library(ggplot2)
# library(ggrepel)
# library(ggforce)
# library(tidyverse)
library(foreach)
library(doParallel)



# models = c("gamye","gam","firstdiff","slope")
# heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance
# 
# species_to_run = c("Wood Thrush", "American Kestrel","Barn Swallow","Chestnut-collared Longspur","Cooper's Hawk","Ruby-throated Hummingbird")


# running comparison models in parallel -----------------------------------

# Set up parallel stuff
n_cores <- length(species_to_run)
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)
########################
##############
############## consider fixing nu at 3, sensu "robust regression" model in Gelman BDA pg 440.

foreach(m = 1:length(species_to_run),
        .packages = 'jagsUI',
        .inorder = FALSE,
        .errorhandling = "pass") %dopar%
  {
    

    species = species_to_run[m]

    sp_dir = paste0("output/",species,"/")
    
loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"),stringsAsFactors = F)

    
year = loo.point$Year-(min(loo.point$Year)-1)
nyears = max(year)

strat = loo.point$Stratum_Factored
nstrat = max(strat)

#fit2 = loo.point$firstdiff
 
 for(comp in c("gamye_firstdiff","gamye_slope")){

   dif = as.numeric(unlist(loo.point[,comp]))
   ncounts = length(dif)
   

 
   

jg.dat = list(
  ncounts = ncounts,
dif = dif,
group = year,
ngroups = nyears
)
# 
# 
# ############ Bayesian model estimating the difference in fit among models and years while accounting for the uncertainty in the point-wise loo
# 
m.year = jagsUI::jags(data = jg.dat,
                      model.file = "summary_models/jags.mod.loo.txt",
                      parameters.to.save = c("nu","difmod","difmod_group","tau"),
                      n.chains = 3,
                      n.burnin = 2000,
                      n.iter = 10000,
                      n.thin = 10,
                      parallel = F)
# 
# 
 tosave2 = c(list(m.year = m.year))
# 
# ############ Same as above but by strataum: Bayesian model estimating the difference in fit among models while accounting for the uncertainty in the point-wise loo
 jg.dat = list(
   ncounts = ncounts,
   dif = dif,
   group = strat,
   ngroups = nstrat
 )
 # 
 # 
 # ############ Bayesian model estimating the difference in fit among models and years while accounting for the uncertainty in the point-wise loo
 # 
 m.strat = jagsUI::jags(data = jg.dat,
                       model.file = "summary_models/jags.mod.loo.txt",
                       parameters.to.save = c("nu","difmod","difmod_group","tau"),
                       n.chains = 3,
                       n.burnin = 2000,
                       n.iter = 12000,
                       n.thin = 10,
                       parallel = F)
 # 
  tosave2 = c(tosave2,
            list(m.strat = m.strat))
# 
# 
# 
# # 
# 
# 
# ############ Same as above but overall: Bayesian model estimating the difference in fit among models
 
  jg.dat = list(
    ncounts = ncounts,
    dif = dif
  )
  # 
  # 
  # ############ Bayesian model estimating the difference in fit among models and years while accounting for the uncertainty in the point-wise loo
  # 
  m.overall = jagsUI::jags(data = jg.dat,
                         model.file = "summary_models/jags.mod.loo.overall.txt",
                         parameters.to.save = c("nu","difmod","tau"),
                         n.chains = 3,
                         n.burnin = 2000,
                         n.iter = 10000,
                         n.thin = 10,
                         parallel = T)
  
  
# 
# ############ Bayesian model estimating the difference in fit among models
# 
 tosave2 = c(tosave2,
             list(m.overall = m.overall))

 
 
 jg.dat = list(
   ncounts = ncounts,
   dif = dif,
   group1 = year,
   ngroups1 = nyears,
   group2 = strat,
   ngroups2 = nstrat,
   
 )
 # 
 # 
 # ############ Bayesian model estimating the difference in fit among models by year and stratum while accounting for the uncertainty in the point-wise loo
 # 
 m.both = jagsUI::jags(data = jg.dat,
                       model.file = "summary_models/jags.mod.loo.both.txt",
                       parameters.to.save = c("nu","difmod","difmod_group","tau","taugroup1","taugroup2"),
                       n.chains = 3,
                       n.burnin = 2000,
                       n.iter = 10000,
                       n.thin = 10,
                       parallel = F)
 # 
 # 
 tosave2 = c(list(m.both = m.both))
 
 
 
 if(comp == "gamye_firstdiff"){
   tosave2out = list(gamye_firstdiff = tosave2,
                   gamye_slope = list())
 }else{
                     tosave2out[["gamye_slope"]] = tosave2
                     }

 }

# 
# 
 save(list = c("tosave2out"),file = paste0(sp_dir,"saved objects2.RData"))
# 
  }


stopCluster(cl = cluster)


# 
# #mixed model examining the effect of strata and model on the fit
# # m_cont = lmer(data = alldat,formula = q0.5.loo ~ model + (1|unit) + (model|Stratum),weights = prec)
# # summary(m_cont)
# #  
# # 
# # 
# # 
# # 
# # 
# # plot(alldat$q0.5.loo,log(alldat$prec)) ### this plot demonstrates that all of the extreme logprob values are very low precision
# # 
# # plot(log(alldat$Count),alldat$q0.5.loo)
# 
# 
# write.csv(loo.point,paste0(sp_dir,"pointwise median  posterior log prob.csv"))


# 
# an.comp = ggplot(data = sum.loo.ry,aes(x = Year,))


#### outliers in fit statistics



##### yearly fit comparisons, overall and by strata


##### strata fit comparisons


#### fit vs count




############### pooling factors on the gam betas

# pooling factor (pg 478, Gelman and Hill) for each stratum and year = [in R code] (sd(yearcar[y,j])/sdyear[y])^2

#min(((beta.X[,k]-B.X[k])/posteriormean(sdbeta))^2,1)


#}




