#comparison_summary script



library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(lme4)


models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

for(species in c("American Kestrel","Barn Swallow","Wood Thrush","Chestnut-collared Longspur")){
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
  
  load(paste0(sp_dir,"saved objects2.RData"))

  
  
  jags_data = tosave2$m.year$model$cluster1$data()



wparam = paste0("difmod_y[",rep(1:6,times = jags_data$nyears),",",rep(1:jags_data$nyears,each = 6),"]")
dif_mod_year = as.data.frame(tosave2$m.year$summary[wparam,])
names(dif_mod_year)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod_year$Contrast = rep(1:6,times = jags_data$nyears)
dif_mod_year$Year = rep(1:jags_data$nyears,each = 6)+1965

contr_names = c(paste(models[1],models[2],sep = "_"),
                paste(models[1],models[3],sep = "_"),
                paste(models[1],models[4],sep = "_"),
                paste(models[2],models[3],sep = "_"),
                paste(models[2],models[4],sep = "_"),
                paste(models[3],models[4],sep = "_")
)

dif_mod_year$Contrast_name = rep(contr_names,times = jags_data$nyears)

an_contr = ggplot(data = dif_mod_year[which(dif_mod_year$Contrast %in% c(2,3)),],aes(x = Year,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  geom_linerange(aes(x = Year,ymin = lci,ymax = uci),position = position_dodge(width = 0.75))+
  coord_cartesian(ylim = c(-0.03,0.03))+
  theme_minimal()+
  labs(title = paste(species,"annual cross validation results"))+
  xlab("")+
  ylab("")+
  scale_x_continuous(breaks = c(seq(1970,2010,by = 10),2018))+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  annotate(geom = "text",x = 1985,y = 0.02,label = "> 0 Favours GAMYE")+
  annotate(geom = "text",x = 1985,y = -0.02,label = "< 0 Favours Alternate")  
  
  
pdf(paste0(sp_dir,species,"annual cross validation.pdf"),
    width = 11,
    height = 8.5)
print(an_contr)
dev.off()



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



str_contr = ggplot(data = dif_mod_strat[which(dif_mod_strat$Contrast %in% c(2,3)),],aes(x = ncounts,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  #geom_linerange(aes(x = ncounts,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.2)+
  coord_cartesian(ylim = c(-0.05,0.05))+
  geom_smooth(stat = "smooth",se = F)+
  theme_minimal()+
  labs(title = paste(species,"geographic cross validation results by sample size"))+
  xlab("Number of route*year counts")+
  ylab("")+
  #scale_x_continuous(breaks = 6)+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  annotate(geom = "text",x = max(dif_mod_strat$ncounts,na.rm = T)/2,y = 0.02,label = "> 0 Favours GAMYE")+
  annotate(geom = "text",x = max(dif_mod_strat$ncounts,na.rm = T)/2,y = -0.02,label = "< 0 Favours Alternate")  


pdf(paste0(sp_dir,species,"geographic cross validation.pdf"),
    width = 11,
    height = 8.5)
print(str_contr)
dev.off()






}


