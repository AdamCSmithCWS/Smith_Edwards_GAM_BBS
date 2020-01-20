#comparison_summary script



library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(lme4)




models = c("gamye","gam","firstdiff","slope")

z_score_sum = function(x){
  summdif = sum(x)
  sddif = sd(x)
  sesumdif = sddif*sqrt(length(x))
  z1 = summdif/sesumdif
  return(list(dif = summdif,
              se = sesumdif,
              z = z1,
              p = pnorm(z1)))
}

p_z_score_mean = function(x){
  mdif = mean(x)
  sddif = sd(x)
  semdif = sddif/sqrt(length(x))
  z1 = mdif/semdif
  return(p = pnorm(z1))
}


z_func = function(x){
  mdif = mean(x)
  sddif = sd(x)
  semdif = sddif/sqrt(length(x))
  z1 = mdif/semdif
  return(z1)
}

lci_func = function(x){
  mdif = mean(x)
  semdif = sd(x)/sqrt(length(x))
  lci = mdif-(1.96*semdif)
  return(lci)
}

uci_func = function(x){
  mdif = mean(x)
  semdif = sd(x)/sqrt(length(x))
  uci = mdif+(1.96*semdif)
  return(uci)
}

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
  
  loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))
  
  
  
nyears = tosave$all_data$gamye$ymax

  loo.point = mutate(loo.point,count_type = if_else(Count > 0,1,0))
  loo.point$all = "All"

# distribution of the BPIC values -----------------------------------------
  my_summarise2 <- function(df, 
                            grp = Year,
                            expr1,expr2,expr3,expr4,expr5,expr6) {
    grp <- enquo(grp)
    expr1 <- enquo(expr1)
    
     out1 = df %>% 
       group_by(!! grp) %>% 
      summarise(.,
                mean = mean(!! expr1),
              lci = lci_func(!! expr1),
              uci = uci_func(!! expr1),
              z = z_func(!! expr1),
              p = p_z_score_mean(!! expr1)
    )
     out1$Contrast = 1
     out1$Contrast_name = quo_name(expr1)
     
 
     expr2 <- enquo(expr2)
     
     out2 = df %>% 
       group_by(!! grp) %>% 
       summarise(.,
                 mean = mean(!! expr2),
                 lci = lci_func(!! expr2),
                 uci = uci_func(!! expr2),
                 z = z_func(!! expr2),
                 p = p_z_score_mean(!! expr2)
       )
     out2$Contrast = 2
     out2$Contrast_name = quo_name(expr2)
     
     expr3 <- enquo(expr3)
     
     out3 = df %>% 
       group_by(!! grp) %>% 
       summarise(.,
                 mean = mean(!! expr3),
                 lci = lci_func(!! expr3),
                 uci = uci_func(!! expr3),
                 z = z_func(!! expr3),
                 p = p_z_score_mean(!! expr3)
       )
     out3$Contrast = 3
     out3$Contrast_name = quo_name(expr3)
     
     expr4 <- enquo(expr4)
     
     out4 = df %>% 
       group_by(!! grp) %>% 
       summarise(.,
                 mean = mean(!! expr4),
                 lci = lci_func(!! expr4),
                 uci = uci_func(!! expr4),
                 z = z_func(!! expr4),
                 p = p_z_score_mean(!! expr4)
       )
     out4$Contrast = 4
     out4$Contrast_name = quo_name(expr4)
     
     expr5 <- enquo(expr5)
     
     out5 = df %>% 
       group_by(!! grp) %>% 
       summarise(.,
                 mean = mean(!! expr5),
                 lci = lci_func(!! expr5),
                 uci = uci_func(!! expr5),
                 z = z_func(!! expr5),
                 p = p_z_score_mean(!! expr5)
       )
     out5$Contrast = 5
     out5$Contrast_name = quo_name(expr5)
     
     expr6 <- enquo(expr6)
     
     out6 = df %>% 
       group_by(!! grp) %>% 
       summarise(.,
                 mean = mean(!! expr6),
                 lci = lci_func(!! expr6),
                 uci = uci_func(!! expr6),
                 z = z_func(!! expr6),
                 p = p_z_score_mean(!! expr6)
       )
     out6$Contrast = 6
     out6$Contrast_name = quo_name(expr6)
     
     
     out = bind_rows(out1,out2,out3,out4,out5,out6)
   return(out)
  }
  

  dif_mod_year <- my_summarise2(loo.point,
                                grp = Year,
                                expr1 = gamye_gam,
                                expr2 = gamye_firstdiff,
                                expr3 = gamye_slope,
                                expr4 = gam_firstdiff,
                                expr5 = gam_slope,
                                expr6 = firstdiff_slope)
  nbyyear <- group_by(loo.point,Year) %>% summarise(.,ncounts = n())
  dif_mod_year <- left_join(dif_mod_year,nbyyear,by = "Year")
  
  dif_mod_strat <- my_summarise2(loo.point,
                                grp = Stratum,
                                expr1 = gamye_gam,
                                expr2 = gamye_firstdiff,
                                expr3 = gamye_slope,
                                expr4 = gam_firstdiff,
                                expr5 = gam_slope,
                                expr6 = firstdiff_slope)
  nbystrat <- group_by(loo.point,Stratum) %>% summarise(.,ncounts = n())
  dif_mod_strat <- left_join(dif_mod_strat,nbystrat,by = "Stratum")
  
  dif_mod_count_type <- my_summarise2(loo.point,
                                 grp = count_type,
                                 expr1 = gamye_gam,
                                 expr2 = gamye_firstdiff,
                                 expr3 = gamye_slope,
                                 expr4 = gam_firstdiff,
                                 expr5 = gam_slope,
                                 expr6 = firstdiff_slope)
  
  
  
  dif_mod <- my_summarise2(loo.point,
                           grp = all,
                           expr1 = gamye_gam,
                           expr2 = gamye_firstdiff,
                           expr3 = gamye_slope,
                           expr4 = gam_firstdiff,
                           expr5 = gam_slope,
                           expr6 = firstdiff_slope
                           )
  
#dif_mod_year$Year = dif_mod_year$Year+(2018-nyears)
dif_mod_year$species = species
dif_mod_year$Contrast_full_name = contrast_full_names[dif_mod_year$Contrast_name]

lbl = dif_mod_year[which(dif_mod_year$Year == 1995 & dif_mod_year$Contrast %in% c(2,3)),]
lbl[which(lbl$Contrast == 2),"Year"] = lbl[which(lbl$Contrast == 2),"Year"]-0.25
lbl[which(lbl$Contrast == 3),"Year"] = lbl[which(lbl$Contrast == 3),"Year"]+0.25



an_contr = ggplot(data = dif_mod_year[which(dif_mod_year$Contrast %in% c(2,3)),],aes(x = Year,y = mean,group = Contrast_name,colour = Contrast_name))+
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
   annotate(geom = "text",x = 2017-nyears,y = 0.017,label = "Favours GAMYE",angle = 90)+
   annotate(geom = "text",x = 2017-nyears,y = -0.017,label = "Favours Alternate",angle = 90)  
  
  
pdf(paste0(sp_dir,species," simple annual cross validation.pdf"),
    width = 11,
    height = 8.5)
print(an_contr)
dev.off()







# Geographic Summaries ----------------------------------------------------





dif_mod_strat$species = species
dif_mod_strat$Contrast_full_name = contrast_full_names[dif_mod_strat$Contrast_name]



#merge original strata names into output


lbl = dif_mod_strat[which(dif_mod_strat$ncounts == quantile(dif_mod_strat$ncounts,1,na.rm= T) & dif_mod_strat$Contrast %in% c(2,3)),]
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
  geom_text_repel(data = lbl,aes(x = ncounts,y = mean,label = Contrast_full_name),nudge_y = -0.01)+
  annotate(geom = "text",x = min(dif_mod_strat_p$ncounts,na.rm = T)*0.90,y = max(dif_mod_strat_p$mean,na.rm = T)*0.90,label = "Favours GAMYE",angle = 90,hjust = 1)+
  annotate(geom = "text",x = min(dif_mod_strat_p$ncounts,na.rm = T)*0.90,y = min(dif_mod_strat_p$mean,na.rm = T)*0.90,label = "Favours Alternate",angle = 90,hjust = 0)  


pdf(paste0(sp_dir,species," simple geographic cross validation.pdf"),
    width = 11,
    height = 8.5)
print(str_contr)
dev.off()










## overall comparison, no stratification


dif_mod$species = species

dif_mod$Contrast_full_name = contrast_full_names[dif_mod$Contrast_name]

overall.y = ggplot(data = dif_mod,aes(x = Contrast_full_name,y = mean))+
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



pdf(paste0(sp_dir,species," simple overall simple cross validation.pdf"),
    width = 7,
    height = 5)
print(overall.y)
dev.off()



dif_mod_count_type$species = species
dif_mod_count_type$Contrast_full_name = contrast_full_names[dif_mod_count_type$Contrast_name]

# lbl = dif_mod_year[which(dif_mod_year$Year == 1995 & dif_mod_year$Contrast %in% c(2,3)),]
# lbl[which(lbl$Contrast == 2),"Year"] = lbl[which(lbl$Contrast == 2),"Year"]-0.25
# lbl[which(lbl$Contrast == 3),"Year"] = lbl[which(lbl$Contrast == 3),"Year"]+0.25



count_type_contr = ggplot(data = dif_mod_count_type[which(dif_mod_count_type$Contrast %in% c(2,3)),],aes(x = count_type,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  geom_linerange(aes(x = count_type,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.5)+
  #coord_cartesian(xlim = c(-0.5,1.5))+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(title = paste(species,"GAMYE cross validation comparison, by count-type"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  # geom_text_repel(data = lbl,aes(x = Year,y = lci,label = Contrast_full_name),nudge_y = -0.005)+
  annotate(geom = "text",x = 2,y = 0.017,label = "Favours GAMYE",angle = 90)+
  annotate(geom = "text",x = 2,y = -0.017,label = "Favours Alternate",angle = 90)  


pdf(paste0(sp_dir,species," simple count type cross validation.pdf"),
    width = 11,
    height = 8.5)
print(count_type_contr)
dev.off()





if(species == demo_sp[1]){
  dif_mod_year_out = dif_mod_year
  dif_mod_strat_out = dif_mod_strat
  dif_mod_out = dif_mod
  dif_mod_count_type_out = dif_mod_count_type
}else{
  dif_mod_year_out = rbind(dif_mod_year_out,dif_mod_year)
  dif_mod_strat_out = rbind(dif_mod_strat_out,dif_mod_strat)
  dif_mod_out = rbind(dif_mod_out,dif_mod)
  dif_mod_count_type_out = rbind(dif_mod_count_type_out,dif_mod_count_type)
}



}
save(list = c("dif_mod_year_out","dif_mod_strat_out","dif_mod_out","dif_mod_count_type_out","models","contrast_full_names","demo_sp"),
     file = "simple_comparison_summary_output.RData")








