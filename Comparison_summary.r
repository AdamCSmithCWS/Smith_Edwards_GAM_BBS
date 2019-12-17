#comparison_summary script



library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(lme4)


models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

for(species in c("American Kestrel","Barn Swallow")){
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
  
  load(paste0(sp_dir,"saved objects2.RData"))

  
  
  jags_data = tosave2$m.year$model$cluster1$data()



wparam = paste0("difmod_y[",rep(1:6,times = jags_data$nyears),",",rep(1:jags_data$nyears,each = 6),"]")
dif_mod_strat = as.data.frame(tosave2$m.year$summary[wparam,])
names(dif_mod_strat)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod_strat$Contrast = rep(1:6,times = jags_data$nyears)
dif_mod_strat$Year = rep(1:jags_data$nyears,each = 6)+1965

contr_names = c(paste(models[1],models[2],sep = "_"),
                paste(models[1],models[3],sep = "_"),
                paste(models[1],models[4],sep = "_"),
                paste(models[2],models[3],sep = "_"),
                paste(models[2],models[4],sep = "_"),
                paste(models[3],models[4],sep = "_")
)

dif_mod_strat$Contrast_name = rep(contr_names,times = jags_data$nyears)

an_contr = ggplot(data = dif_mod_strat[which(dif_mod_strat$Contrast < 4),],aes(x = Year,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 1))+
  geom_linerange(aes(x = Year,ymin = lci,ymax = uci),position = position_dodge(width = 1))

print(an_contr)




jags_data = tosave2$m.strat$model$cluster1$data()



wparam = paste0("difmod_s[",rep(1:6,times = jags_data$nstrat),",",rep(1:jags_data$nstrat,each = 6),"]")
dif_mod_strat = as.data.frame(tosave2$m.strat$summary[wparam,])
names(dif_mod_strat)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod_strat$Contrast = rep(1:6,times = jags_data$nstrat)
dif_mod_strat$Stratum = rep(1:jags_data$nstrat,each = 6)

contr_names = c(paste(models[1],models[2],sep = "_"),
                paste(models[1],models[3],sep = "_"),
                paste(models[1],models[4],sep = "_"),
                paste(models[2],models[3],sep = "_"),
                paste(models[2],models[4],sep = "_"),
                paste(models[3],models[4],sep = "_")
)

dif_mod_strat$Contrast_name = rep(contr_names,times = jags_data$nstrat)

ncounts = tapply(jags_data$unit,jags_data$strat,FUN = function(x){length(unique(x))})
for(i in 1:nrow(dif_mod_strat)){
  dif_mod_strat[i,"ncounts"] <- ncounts[dif_mod_strat[i,"Stratum"]]
}


an_contr = ggplot(data = dif_mod_strat[which(dif_mod_strat$Contrast < 4),],aes(x = ncounts,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 1))+
  #geom_linerange(aes(x = ncounts,ymin = lci,ymax = uci),position = position_dodge(width = 1))+
  coord_cartesian(ylim = c(-0.1,0.1))

print(an_contr)





}


