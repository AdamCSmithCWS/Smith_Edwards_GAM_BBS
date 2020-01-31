

library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(lme4)
library(sf)
library(RColorBrewer)



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


demo_sp <- c("Horned Lark",
             "American Kestrel",
             "Barn Swallow",
             "Wood Thrush",
             "Chestnut-collared Longspur",
             #"Cooper's Hawk",
             "Ruby-throated Hummingbird")




load("comparison_summary_output.RData")



source("colourblind safe qualitative pallete.r")
species_pallete <- safe.pallet[[length(demo_sp)]] 
names(species_pallete) <- demo_sp




fmod <- function(x){
  str_sub(x,start = 1,end = str_locate(x,pattern = " vs ")[,1]-1)
}
smod <- function(x){
  str_sub(x,start = str_locate(x,pattern = " vs ")[,2]+1,end = str_length(x))
}

dif_mod_year_over_out$M1 = paste("Favours",fmod(dif_mod_out$Contrast_full_name))
dif_mod_year_over_out$M2 = paste("Favours",smod(dif_mod_out$Contrast_full_name))


dif_mod_labs = filter(dif_mod_year_over_out,species == demo_sp[1])
dif_mod_labs$M1loc = 0.0025
dif_mod_labs$M2loc = -0.0025


# overall summary graphs --------------------------------------------------

overall.comparison = ggplot(data = dif_mod_year_over_out,aes(x = Contrast_full_name,y = mean,group = species,colour = species))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
  # theme(legend.position = "none",
  #       axis.text.x = element_text(angle = 90))+
  #labs(title = paste("Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point(position = position_dodge(width = 0.4))+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),
                 alpha = 0.8,position = position_dodge(width = 0.4))+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M1loc,label = M1),show.legend = F,colour = grey(0.3),nudge_x = -0.3)+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M2loc,label = M2),show.legend = F,colour = grey(0.3),nudge_x = -0.3)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour"))+
  # annotate(geom = "text",x = 0.5,y = 0.005,label = "Favours first")+
  # annotate(geom = "text",x = 0.5,y = -0.005,label = "Favours second")+
  #scale_x_discrete(position = "top")+
  guides(colour = guide_legend(reverse = T))+
  coord_flip()


pdf(paste0("overall annual-strat simple cross validation.pdf"),
    width = 10,
    height = 5)
print(overall.comparison)
dev.off()






# Annual Comparison graphs ------------------------------------------------

dif_mod_year_out$M1 = paste("Favours",fmod(dif_mod_year_out$Contrast_full_name))
dif_mod_year_out$M2 = paste("Favours",smod(dif_mod_year_out$Contrast_full_name))

# selecting only two of the possible model comparisons
dif_y_p = filter(dif_mod_year_out,Contrast_name %in% c("gamye_firstdiff","gamye_slope"))

dif_mod_labs = filter(dif_y_p,species == demo_sp[1] & Year == 1994)
dif_mod_labs$M1loc = 0.015
dif_mod_labs$M2loc = -0.015



yearly.comparison = ggplot(data = dif_y_p,aes(x = Year,y = mean,group = species,colour = species))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
  # theme(legend.position = "none",
  #       axis.text.x = element_text(angle = 90))+
  #labs(title = paste("Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point(position = position_dodge(width = 0.4))+
  geom_linerange(aes(x = Year,ymin = lci,ymax = uci),
                 alpha = 0.3,position = position_dodge(width = 0.4))+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Year,y = M1loc,label = M1),show.legend = F,colour = grey(0.3),nudge_x = 0.5,angle = 0)+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Year,y = M2loc,label = M2),show.legend = F,colour = grey(0.3),nudge_x = 0.5,angle = 0)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour"))+
  # annotate(geom = "text",x = 0.5,y = 0.005,label = "Favours first")+
  # annotate(geom = "text",x = 0.5,y = -0.005,label = "Favours second")+
  #scale_x_discrete(position = "top")+
  guides(colour = guide_legend(reverse = F))+
  #coord_flip()+
  facet_wrap(facets = ~Contrast_full_name,nrow = 2,scales = "free",shrink = T)



pdf(paste0("annual cross validation.pdf"),
    width = 10,
    height = 8)
print(yearly.comparison)
dev.off()




# Strata sample size comparison -------------------------------------------


dif_mod_strat_out$M1 = paste("Favours",fmod(dif_mod_strat_out$Contrast_full_name))
dif_mod_strat_out$M2 = paste("Favours",smod(dif_mod_strat_out$Contrast_full_name))

# selecting only two of the possible model comparisons
dif_s_p = filter(dif_mod_strat_out,Contrast_name %in% c("gamye_firstdiff"))


mxcounts = max(dif_s_p$ncounts)


dif_mod_labs = dif_s_p[1,]
dif_mod_labs$M1loc = 0.015
dif_mod_labs$M2loc = -0.015
dif_mod_labs$ncounts = mxcounts*0.8


n.comparison = ggplot(data = dif_s_p,aes(x = ncounts,y = mean,group = species,colour = species))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
   theme(legend.position = "none")+
  #       axis.text.x = element_text(angle = 90))+
  #labs(title = paste("Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point()+
  #geom_smooth(method = "loess",formula = y ~ x,se = F)+
  # geom_linerange(aes(x = Year,ymin = lqrt,ymax = uqrt),
  #                alpha = 0.3,position = position_dodge(width = 0.4))+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = ncounts,y = M1loc,label = M1),show.legend = F,colour = grey(0.3),nudge_x = 0.5,angle = 0)+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = ncounts,y = M2loc,label = M2),show.legend = F,colour = grey(0.3),nudge_x = 0.5,angle = 0)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour"))+
  # annotate(geom = "text",x = 0.5,y = 0.005,label = "Favours first")+
  # annotate(geom = "text",x = 0.5,y = -0.005,label = "Favours second")+
  #scale_x_discrete(position = "top")+
  guides(colour = guide_legend(reverse = F))+
  #coord_flip()+
  facet_wrap(facets = ~species,ncol = 2,nrow = 3,scales = "free_y",shrink = T)



pdf(paste0("sample size by strata cross validation.pdf"),
    width = 7,
    height = 8)
print(n.comparison)
dev.off()






# Strata mapping comparison -----------------------------------------------

# load the strata map


mp = read_sf(dsn = "map/BBS_USGS_strata.shp")


for(comp in c("gamye_firstdiff","gamye_slope")){
  pdf(file = paste0(comp,"map of fit comparison.pdf"))
  for(spp in demo_sp){
  dif_s_p = filter(dif_mod_strat_out,Contrast_name == comp & species == spp)
  dif_s_p$ST_12 = as.character(dif_s_p$Stratum)
  dif_s_p$cl = cut(dif_s_p$mean,breaks = c(-Inf,seq(-0.01,0.01,by = 0.005),Inf))
  
  modcomp_pallete <- brewer.pal(length(levels(dif_s_p$cl)),"RdYlBu") 
  names(modcomp_pallete) <- levels(dif_s_p$cl)
  
  
  mp2 = left_join(x = mp,y = dif_s_p,by = "ST_12")
   
  mp.plot = ggplot()+
    geom_sf(data = mp2,aes(fill = cl),colour = grey(0.6))+
    theme_minimal()+
    labs(title = paste(spp,comp))+
    guides(fill = guide_legend(reverse = T))+
    scale_colour_manual(values = modcomp_pallete, aesthetics = c("fill"))
    
    
  print(mp.plot)
  
  }
  dev.off()
}










