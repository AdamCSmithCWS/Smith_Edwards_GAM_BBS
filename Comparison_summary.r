#comparison_summary script



library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(lme4)




models = c("gamye","gam","firstdiff","slope")


contr_names = c(paste(models[1],models[2],sep = "_"),
                paste(models[1],models[3],sep = "_"),
                paste(models[1],models[4],sep = "_"),
                paste(models[2],models[3],sep = "_"),
                paste(models[2],models[4],sep = "_"),
                paste(models[3],models[4],sep = "_")
)

contrast_full_names = gsub(gsub(toupper(contr_names),pattern = "FIRSTDIFF",replacement = "DIFFERENCE",fixed = T),pattern = "_",replacement = " vs ",fixed = T)
names(contrast_full_names) = contr_names


heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance


demo_sp <- c("American Kestrel",
             "Barn Swallow",
             "Wood Thrush",
             "Chestnut-collared Longspur",
             "Cooper's Hawk",
             "Ruby-throated Hummingbird")

for(species in demo_sp){
  

  
  
  sp_dir = paste0("output/",species,"/")
   
  load(paste0(sp_dir,"saved objects.RData"))
  
  load(paste0(sp_dir,"saved objects2.RData"))


# extract the comparison model results and compile for graphing -----------

  
  for(comp in c("gamye_firstdiff","gamye_slope")){

  jags_data = tosave2out[[comp]]$m.year$model$cluster1$data()
 m.year = tosave2out[[comp]]$m.year


# distribution of the BPIC values -----------------------------------------

#  alldat = tosave$alldat
  
  
  
  
wparam = paste0("difmod_group[",1:jags_data$ngroups,"]")
dif_mod_year1 = as.data.frame(m.year$summary[wparam,])
names(dif_mod_year1)[3:7] <- c("lci","lqrt","med","uqrt","uci")


dif_mod_year1$Contrast_name = comp
dif_mod_year1$Year = (1:jags_data$ngroups)+(2018-jags_data$ngroups)
dif_mod_year1$species = species
dif_mod_year1$Contrast_full_name = contrast_full_names[dif_mod_year1$Contrast_name]

if(comp == "gamye_firstdiff"){
  dif_mod_year = dif_mod_year1
}else{
  dif_mod_year = rbind(dif_mod_year,dif_mod_year1)
}

}




lbl = dif_mod_year[which(dif_mod_year$Year == 1995),]
lbl[which(lbl$Contrast == 2),"Year"] = lbl[which(lbl$Contrast_name == "gamye_firstdiff"),"Year"]-0.25
lbl[which(lbl$Contrast == 3),"Year"] = lbl[which(lbl$Contrast_name == "gamye_slope"),"Year"]+0.25



an_contr = ggplot(data = dif_mod_year,aes(x = Year,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  geom_linerange(aes(x = Year,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.5)+
  coord_cartesian(ylim = c(-0.03,0.03))+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(title = paste(species,"GAMYE cross validation comparison, by year"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  scale_x_continuous(breaks = c(seq(1970,2010,by = 10),2018))+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_text_repel(data = lbl,aes(x = Year,y = lci,label = Contrast_full_name),nudge_y = -0.005)+
   annotate(geom = "text",x = 2017-jags_data$ngroups,y = 0.017,label = "Favours GAMYE",angle = 90)+
   annotate(geom = "text",x = 2017-jags_data$ngroups,y = -0.017,label = "Favours Alternate",angle = 90)  
  
  
pdf(paste0(sp_dir,species," annual cross validation.pdf"),
    width = 11,
    height = 8.5)
print(an_contr)
dev.off()



## plot the overall differences in model fit by pairwise comparisons


wparam = paste0("difmod[",6:1,"]")
dif_mod = as.data.frame(tosave2$m.year$summary[wparam,])
names(dif_mod)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod$Contrast = 6:1


dif_mod$Contrast_name = contr_names[dif_mod$Contrast]
dif_mod$Contrast_full_name = contrast_full_names[dif_mod$Contrast]

overcont.y = ggplot(data = dif_mod,aes(x = Contrast_full_name,y = med))+
  coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))+
  labs(title = paste(species,"Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point()+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),alpha = 0.5)+
annotate(geom = "text",x = 0.2,y = 0.017,label = "Favours first",angle = 90)+
  annotate(geom = "text",x = 0.2,y = -0.017,label = "Favours second",angle = 90)  



pdf(paste0(sp_dir,species," overall using annual model cross validation.pdf"),
    width = 7,
    height = 5)
print(overcont.y)
dev.off()





# Geographic Summaries ----------------------------------------------------






jags_data = tosave2$m.strat$model$cluster1$data()

strat.list = unique(tosave$alldat[,c("Stratum","Stratum_Factored")])
names(strat.list) = c("Stratum_name","Stratum")

wparam = paste0("difmod_s[",rep(1:6,times = jags_data$nstrat),",",rep(1:jags_data$nstrat,each = 6),"]")
dif_mod_strat = as.data.frame(tosave2$m.strat$summary[wparam,])
names(dif_mod_strat)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod_strat$Contrast = rep(1:6,times = jags_data$nstrat)
dif_mod_strat$Stratum = rep(1:jags_data$nstrat,each = 6)
dif_mod_strat$species = species

dif_mod_strat$Contrast_name = rep(contr_names,times = jags_data$nstrat)
dif_mod_strat$Contrast_full_name = rep(contrast_full_names,times = jags_data$nstrat)

ncounts = tapply(jags_data$unit,jags_data$strat,FUN = function(x){length(unique(x))})
for(i in 1:nrow(dif_mod_strat)){
  dif_mod_strat[i,"ncounts"] <- ncounts[dif_mod_strat[i,"Stratum"]]
}

#merge original strata names into output
dif_mod_strat = merge(dif_mod_strat,strat.list,by = "Stratum",sort = F)


lbl = dif_mod_strat[which(dif_mod_strat$ncounts == quantile(dif_mod_strat$ncounts,1) & dif_mod_strat$Contrast %in% c(2,3)),]
# lbl[which(lbl$Contrast == 2),"ncounts"] = lbl[which(lbl$Contrast == 2),"ncounts"]*1.05
# lbl[which(lbl$Contrast == 3),"ncounts"] = lbl[which(lbl$Contrast == 3),"ncounts"]*0.95


dif_mod_strat_p = dif_mod_strat[which(dif_mod_strat$Contrast %in% c(2,3)),]

str_contr = ggplot(data = dif_mod_strat_p,aes(x = ncounts,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  #geom_linerange(aes(x = ncounts,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.2)+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  geom_smooth(stat = "smooth",se = F)+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(title = paste(species,"geographic cross validation results by sample size"))+
  xlab("Number of route*year counts")+
  ylab("")+
  #scale_x_continuous(breaks = 6)+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_text_repel(data = lbl,aes(x = ncounts,y = med,label = Contrast_full_name),nudge_y = -0.01)+
  annotate(geom = "text",x = min(dif_mod_strat_p$ncounts,na.rm = T)*0.90,y = max(dif_mod_strat_p$med,na.rm = T)*0.90,label = "Favours GAMYE",angle = 90,hjust = 1)+
  annotate(geom = "text",x = min(dif_mod_strat_p$ncounts,na.rm = T)*0.90,y = min(dif_mod_strat_p$med,na.rm = T)*0.90,label = "Favours Alternate",angle = 90,hjust = 0)  


pdf(paste0(sp_dir,species," geographic cross validation.pdf"),
    width = 11,
    height = 8.5)
print(str_contr)
dev.off()









wparam = paste0("difmod[",6:1,"]")
dif_mod = as.data.frame(tosave2$m.strat$summary[wparam,])
names(dif_mod)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod$Contrast = 6:1


dif_mod$Contrast_name = contr_names[dif_mod$Contrast]
dif_mod$Contrast_full_name = contrast_full_names[dif_mod$Contrast]

overcont.y = ggplot(data = dif_mod,aes(x = Contrast_full_name,y = med))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))+
  labs(title = paste(species,"Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point()+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),alpha = 0.5)+
  annotate(geom = "text",x = 0.2,y = 0.017,label = "Favours first",angle = 90)+
  annotate(geom = "text",x = 0.2,y = -0.017,label = "Favours second",angle = 90)  



pdf(paste0(sp_dir,species,"  overall using strata model cross validation.pdf"),
    width = 7,
    height = 5)
print(overcont.y)
dev.off()



## overall comparison, no stratification


### load overal model


wparam = paste0("difmod[",6:1,"]")
dif_mod = as.data.frame(tosave2$m.overall$summary[wparam,])
names(dif_mod)[3:7] <- c("lci","lqrt","med","uqrt","uci")
dif_mod$Contrast = 6:1
dif_mod$species = species

dif_mod$Contrast_name = contr_names[dif_mod$Contrast]
dif_mod$Contrast_full_name = contrast_full_names[dif_mod$Contrast]

overall.y = ggplot(data = dif_mod,aes(x = Contrast_full_name,y = med))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))+
  labs(title = paste(species,"Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point()+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),alpha = 0.5)+
  annotate(geom = "text",x = 0.2,y = 0.017,label = "Favours first",angle = 90)+
  annotate(geom = "text",x = 0.2,y = -0.017,label = "Favours second",angle = 90)  



pdf(paste0(sp_dir,species,"  overall simple cross validation.pdf"),
    width = 7,
    height = 5)
print(overall.y)
dev.off()



if(species == demo_sp[1]){
  dif_mod_year_out = dif_mod_year
  dif_mod_strat_out = dif_mod_strat
  dif_mod_out = dif_mod
}else{
  dif_mod_year_out = rbind(dif_mod_year_out,dif_mod_year)
  dif_mod_strat_out = rbind(dif_mod_strat_out,dif_mod_strat)
  dif_mod_out = rbind(dif_mod_out,dif_mod)
}



}
save(list = c("dif_mod_year_out","dif_mod_strat_out","dif_mod_out","models","contrast_full_names","demo_sp"),
     file = "comparison_summary_output.RData")








