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

jj = 1
for(species in demo_sp){
  

  
  
  sp_dir = paste0("output/",species,"/")
   
  load(paste0(sp_dir,"saved objects.RData"))
  
  load(paste0(sp_dir,"saved objects2.RData"))
  loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))
  

# extract the comparison model results and compile for graphing -----------

# comparison accounting for annual variation ------------------------------


  
  for(comp in c("gamye_firstdiff","gamye_slope")){

  jags_data = tosave2out[[comp]]$m.year$model$data()
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
  
  # for(gp in 1:jags_data$ngroups){
  # df = data.frame(dif = jags_data$dif[which(jags_data$group == gp)])
  # dfm = m.year$mean$nu+0.5
  # ncp = m.year$mean$difmod_group[gp]
  # qq = ggplot(data = df,aes(sample = dif))+
  #   geom_qq(distribution = stats::qt,
  #           dparams = list(df = dfm,ncp = ncp))+
  #   geom_qq_line(distribution = stats::qt,
  #                dparams = list(df = dfm,ncp = ncp),
  #                line.p = c(0.25,0.75))
  # png(file = paste0(sp_dir,gp," qq plot t-df-modeled gamye_firstdiff.png"))
  # print(qq)
  # dev.off()
  # 
  # }
  # 
  
  qq = ggplot(data = df,aes(sample = dif))+
      geom_qq(distribution = stats::qt,
              dparams = list(df = m.year$mean$nu))+
      geom_qq_line(distribution = stats::qt,
                   dparams = list(df = m.year$mean$nu),
                   line.p = c(0.25,0.75))
  png(file = paste0(sp_dir,"qq plot t-df-modeled gamye_firstdiff.png"))
  print(qq)
  dev.off()
  
  qq = ggplot(data = loo.point,aes(sample = gamye_firstdiff))+
    geom_qq()+
    geom_qq_line()
  
  png(file = paste0(sp_dir,"qq plot gamye_firstdiff.png"))
  print(qq)
  dev.off()
  
  dif_mod_year = dif_mod_year1
}else{
  
  qq = ggplot(data = loo.point,aes(sample = gamye_slope))+
    geom_qq(distribution = stats::qt,
            dparams = list(df = 2))+
    geom_qq_line(distribution = stats::qt,
                 dparams = list(df = 2),
                 line.p = c(0.001,0.999))
  
  png(file = paste0(sp_dir,"qq plot t-df-modeled gamye_slope.png"))
  print(qq)
  dev.off()
  
  qq = ggplot(data = loo.point,aes(sample = gamye_slope))+
    geom_qq()+
    geom_qq_line()
  
  png(file = paste0(sp_dir,"qq plot gamye_slope.png"))
  print(qq)
  dev.off()
  
  
  dif_mod_year = rbind(dif_mod_year,dif_mod_year1)
}







}




lbl = dif_mod_year[which(dif_mod_year$Year == 1995),]
lbl[which(lbl$Contrast == 2),"Year"] = lbl[which(lbl$Contrast_name == "gamye_firstdiff"),"Year"]-0.25
lbl[which(lbl$Contrast == 3),"Year"] = lbl[which(lbl$Contrast_name == "gamye_slope"),"Year"]+0.25



an_contr = ggplot(data = dif_mod_year,aes(x = Year,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  geom_linerange(aes(x = Year,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.5)+
  #coord_cartesian(ylim = c(-0.03,0.03))+
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



## plot the overall differences in model fit accounting for annual variation


for(comp in c("gamye_firstdiff","gamye_slope")){
  


  jags_data = tosave2out[[comp]]$m.year$model$data()
  m.year = tosave2out[[comp]]$m.year
  

  
  wparam = "difmod"
  dif_mod_year_over1 = data.frame(matrix(m.year$summary[wparam,],nrow = 1))
  names(dif_mod_year_over1) <- names(m.year$summary[1,])
  names(dif_mod_year_over1)[3:7] <- c("lci","lqrt","med","uqrt","uci")
  
  
  dif_mod_year_over1$Contrast_name = comp
  dif_mod_year_over1$species = species
  dif_mod_year_over1$Contrast_full_name = contrast_full_names[dif_mod_year_over1$Contrast_name]
  
  if(comp == "gamye_firstdiff"){
    dif_mod_year_over = dif_mod_year_over1
  }else{
    dif_mod_year_over = rbind(dif_mod_year_over,dif_mod_year_over1)
  }
  
}




overcont.y = ggplot(data = dif_mod_year_over,aes(x = Contrast_full_name,y = med))+
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
annotate(geom = "text",x = 0.7,y = 0.5*max(dif_mod_year_over$uci),label = "Favours first",angle = 90)+
  annotate(geom = "text",x = 0.7,y = -0.5*max(dif_mod_year_over$uci),label = "Favours second",angle = 90)  



pdf(paste0(sp_dir,species," overall using annual model cross validation.pdf"),
    width = 7,
    height = 5)
print(overcont.y)
dev.off()






# Geographic Summaries ----------------------------------------------------


for(comp in c("gamye_firstdiff","gamye_slope")){
  
  jags_data = tosave2out[[comp]]$m.strat$model$data()
  m.strat = tosave2out[[comp]]$m.strat
  
  
  
  wparam = paste0("difmod_group[",1:jags_data$ngroups,"]")
  dif_mod_strat1 = as.data.frame(m.strat$summary[wparam,])
  names(dif_mod_strat1)[3:7] <- c("lci","lqrt","med","uqrt","uci")
  
 
  
  dif_mod_strat1$Contrast_name = comp
  dif_mod_strat1$Stratum_Factored = (1:jags_data$ngroups)
  dif_mod_strat1$species = species
  dif_mod_strat1$Contrast_full_name = contrast_full_names[dif_mod_strat1$Contrast_name]
  
  if(comp == "gamye_firstdiff"){
    dif_mod_strat = dif_mod_strat1
  }else{
    dif_mod_strat = rbind(dif_mod_strat,dif_mod_strat1)
  }
  
}




strat.list = unique(loo.point[,c("Stratum","Stratum_Factored")])

dif_mod_strat <- left_join(dif_mod_strat,strat.list,by = "Stratum_Factored")

nbystrat <- group_by(loo.point,Stratum) %>% summarise(.,ncounts = n())
dif_mod_strat <- left_join(dif_mod_strat,nbystrat,by = "Stratum")


# lbl = dif_mod_strat[which(dif_mod_strat$strat == 1995),]
# lbl[which(lbl$Contrast == 2),"strat"] = lbl[which(lbl$Contrast_name == "gamye_firstdiff"),"strat"]-0.25
# lbl[which(lbl$Contrast == 3),"strat"] = lbl[which(lbl$Contrast_name == "gamye_slope"),"strat"]+0.25
lbl = dif_mod_strat[which(dif_mod_strat$ncounts == quantile(dif_mod_strat$ncounts,1)),]


str_contr = ggplot(data = dif_mod_strat,aes(x = ncounts,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  geom_linerange(aes(x = ncounts,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.5)+
  coord_cartesian(ylim = c(-0.03,0.03))+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(title = paste(species,"GAMYE cross validation comparison, by strat"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  #scale_x_continuous(breaks = c(seq(1970,2010,by = 10),2018))+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
     geom_text_repel(data = lbl,aes(x = ncounts,y = med,label = Contrast_full_name),nudge_y = -0.01)+
     annotate(geom = "text",x = min(dif_mod_strat$ncounts,na.rm = T)*0.90,y = max(dif_mod_strat$med,na.rm = T)*0.90,label = "Favours GAMYE",angle = 90,hjust = 1)+
     annotate(geom = "text",x = min(dif_mod_strat$ncounts,na.rm = T)*0.90,y = min(dif_mod_strat$med,na.rm = T)*0.90,label = "Favours Alternate",angle = 90,hjust = 0)  
  

pdf(paste0(sp_dir,species," geographic cross validation.pdf"),
    width = 11,
    height = 8.5)
print(str_contr)
dev.off()


# mapping distribution of the BPIC values -----------------------------------------



# overall using geographic stratification  -----------------------------------------

for(comp in c("gamye_firstdiff","gamye_slope")){
  
  
  
  jags_data = tosave2out[[comp]]$m.strat$model$data()
  m.strat = tosave2out[[comp]]$m.strat
  
  wparam = "difmod"
  dif_mod_strat_over1 = data.frame(matrix(m.strat$summary[wparam,],nrow = 1))
  names(dif_mod_strat_over1) <- names(m.strat$summary[1,])
  names(dif_mod_strat_over1)[3:7] <- c("lci","lqrt","med","uqrt","uci")
  
  
  dif_mod_strat_over1$Contrast_name = comp
  dif_mod_strat_over1$species = species
  dif_mod_strat_over1$Contrast_full_name = contrast_full_names[dif_mod_strat_over1$Contrast_name]
  
  if(comp == "gamye_firstdiff"){
    dif_mod_strat_over = dif_mod_strat_over1
  }else{
    dif_mod_strat_over = rbind(dif_mod_strat_over,dif_mod_strat_over1)
  }
  
}





overcont.y = ggplot(data = dif_mod_strat_over,aes(x = Contrast_full_name,y = med))+
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
  annotate(geom = "text",x = 0.7,y = 0.5*max(dif_mod_strat_over$uci),label = "Favours first",angle = 90)+
  annotate(geom = "text",x = 0.7,y = -0.5*max(dif_mod_strat_over$uci),label = "Favours second",angle = 90)  



pdf(paste0(sp_dir,species,"  overall using strata model cross validation.pdf"),
    width = 7,
    height = 5)
print(overcont.y)
dev.off()



















# overall comparison, no stratification -----------------------------------



### load overal model
for(comp in c("gamye_firstdiff","gamye_slope")){
  
  jags_data = tosave2out[[comp]]$m.overall$model$cluster1$data()
  m.overall = tosave2out[[comp]]$m.overall
  
  
  # distribution of the BPIC values -----------------------------------------
  
  #  alldat = tosave$alldat
  
  
  
  
  wparam = "difmod"
  dif_mod_over1 = data.frame(matrix(m.overall$summary[wparam,],nrow = 1))
  names(dif_mod_over1) <- names(m.overall$summary[1,])
  names(dif_mod_over1)[3:7] <- c("lci","lqrt","med","uqrt","uci")
  
  
  dif_mod_over1$Contrast_name = comp
  dif_mod_over1$species = species
  dif_mod_over1$Contrast_full_name = contrast_full_names[dif_mod_over1$Contrast_name]
  
  if(comp == "gamye_firstdiff"){
    dif_mod = dif_mod_over1
  }else{
    dif_mod = rbind(dif_mod,dif_mod_over1)
  }
  
}



overall = ggplot(data = dif_mod,aes(x = Contrast_full_name,y = med))+
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
  annotate(geom = "text",x = 0.7,y = 0.5*max(dif_mod$uci),label = "Favours first",angle = 90)+
  annotate(geom = "text",x = 0.7,y = -0.5*max(dif_mod$uci),label = "Favours second",angle = 90)  



pdf(paste0(sp_dir,species,"  overall simple cross validation.pdf"),
    width = 7,
    height = 5)
print(overall)
dev.off()



if(jj == 1){
  dif_mod_year_out = dif_mod_year
  dif_mod_strat_out = dif_mod_strat
  dif_mod_year_over_out = dif_mod_year_over
  dif_mod_strat_over_out = dif_mod_strat_over
  dif_mod_out = dif_mod
}else{
  dif_mod_year_out = rbind(dif_mod_year_out,dif_mod_year)
  dif_mod_strat_out = rbind(dif_mod_strat_out,dif_mod_strat)
  dif_mod_year_over_out = rbind(dif_mod_year_over_out,dif_mod_year_over)
  dif_mod_strat_over_out = rbind(dif_mod_strat_over_out,dif_mod_strat_over)
  dif_mod_out = rbind(dif_mod_out,dif_mod)
}

jj = jj+1

}
save(list = c("dif_mod_year_out","dif_mod_year_over_out","dif_mod_strat_out","dif_mod_strat_over_out","dif_mod_out","models","contrast_full_names","demo_sp"),
     file = "comparison_summary_output.RData")








