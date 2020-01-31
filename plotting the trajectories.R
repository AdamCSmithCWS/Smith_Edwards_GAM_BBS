


# plotting trajectories for GAMYE ------------------------------------------

### comparing the results of the k-fold cross validation
library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)


models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

species_to_run = c("Horned Lark","Wood Thrush", "American Kestrel","Barn Swallow","Chestnut-collared Longspur","Cooper's Hawk","Ruby-throated Hummingbird")

### colour pallette

source("colourblind safe qualitative pallete.r")
model_pallete <- safe.pallet[[4]] 
model_pallete <- model_pallete[c(2,1,3,4)]

names(model_pallete) <- models

for(species in species_to_run[1:3]){
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
  
inds = tosave$indsout
  
# names(inds)[which(names(inds) %in% c("Index_q_0.025","Index_q_0.05","Index_q_0.25","Index_q_0.75",
#                                      "Index_q_0.95","Index_q_0.975"))] <- c("lci","lci5","lqrt",
#                                                                             "uqrt","uci5","uci")

datt = tosave$datt
ncby_y = table(datt$Year)
ncby_y = round(ncby_y/50)

dattc = data.frame(Year = rep(as.integer(names(ncby_y)),times = ncby_y))

  
indcont = inds[which(inds$Region == "Continental"),]  

indcont2 = indcont[which(indcont$model == "slope"),]

uylim = max(c(indcont$Index_q_0.975,indcont$obs_mean))

labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
labl_obs$label = "Observed mean counts"

mxy = tapply(indcont[which(indcont$Year %in% c(1970:2013)),"Index"],indcont[which(indcont$Year %in% c(1970:2013)),"model"],max)

labl_mods = indcont[which(indcont$Index %in% mxy),]
labl_mods$Model = toupper(labl_mods$model)
labl_mods$Model[which(labl_mods$Model == "FIRSTDIFF")] <- "DIFFERENCE"

cont_over = ggplot(data = indcont[-which(indcont$model %in% models[c(3,4)]),],aes(x = Year,y = Index,group = model))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(label = "")+
  ylab(label = "Predicted mean count")+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  coord_cartesian(ylim = c(0,uylim))+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_text_repel(data = labl_mods[4,],aes(x = Year,y = Index,label = Model),colour = grey(0.5), nudge_y = 0.075*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_mods[1:2,],aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = 5)+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
  geom_ribbon(data = indcont[which(indcont$model %in% models[c(4)]),],aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975),fill = grey(0.5),alpha = 0.2)+
  geom_line(data = indcont[which(indcont$model %in% models[c(4)]),],aes(x = Year,y = Index),colour = grey(0.7),size = 1.2)+
  
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
  geom_line(aes(colour = model),size = 1.2)+
  geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)

pdf(file = paste0(sp_dir,"Fig 1.pdf"),
    width = 5,
    height = 4)
print(cont_over)
 dev.off()
 
  
 pdf(file = paste0(sp_dir,"all Fig 1.pdf"),
     width = 5,
     height = 4)
 

for(pp in unique(inds$Region)){
  
  if(pp == "Continental"){dattt = dattc}else{dattt = datt[which(datt$Stratum == pp),] }
  indcont = inds[which(inds$Region == pp),]  
  
  indcont2 = indcont[which(indcont$model == "slope"),]
  
  uylim = max(c(indcont$Index_q_0.975,indcont$obs_mean))
  
  labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
  labl_obs$label = "Observed mean counts"
  
  mxy = tapply(indcont[which(indcont$Year %in% c(1970:2013)),"Index"],indcont[which(indcont$Year %in% c(1970:2013)),"model"],max)
  
  labl_mods = indcont[which(indcont$Index %in% mxy),]
  labl_mods$Model = toupper(labl_mods$model)
  labl_mods$Model[which(labl_mods$Model == "FIRSTDIFF")] <- "DIFFERENCE"
  
  cont_over = ggplot(data = indcont[-which(indcont$model %in% models[c(3,4)]),],aes(x = Year,y = Index,group = model))+
    theme_classic()+
    theme(legend.position = "none")+
    xlab(label = "")+
    ylab(label = "Predicted mean count")+
    labs(title = paste(pp))+
    geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
    coord_cartesian(ylim = c(0,uylim))+
    scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    geom_text_repel(data = labl_mods[4,],aes(x = Year,y = Index,label = Model),colour = grey(0.5), nudge_y = 0.075*uylim, nudge_x = 5)+
    geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
    geom_text_repel(data = labl_mods[1:2,],aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = 5)+
    scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
    geom_ribbon(data = indcont[which(indcont$model %in% models[c(4)]),],aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975),fill = grey(0.5),alpha = 0.2)+
    geom_line(data = indcont[which(indcont$model %in% models[c(4)]),],aes(x = Year,y = Index),colour = grey(0.7),size = 1.2)+
    
    geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
    geom_line(aes(colour = model),size = 1.2)+
    geom_dotplot(data = dattt,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)
  
  print(cont_over)

  
}
 
 dev.off()
 
 
 
 
 
 

# GAMYE components plot ---------------------------------------------------

load(paste0(sp_dir,"gamye/parameter_model_run.RData"))
 load(paste0(sp_dir,"gamye/jags_mod_full.RData"))
 
  load(paste0(sp_dir,"gamye/jags_data.RData"))
 
inds_gamye <- generate_regional_indices(jags_mod = jags_mod_full,
                                        jags_data = jags_data,
                                        regions = c("continental","stratum"),
                                        max_backcast = NULL,
                                        alternate_n = "n") 

inds_gamnoye <- generate_regional_indices(jags_mod = jags_mod_param,
                                        jags_data = jags_data,
                                        regions = c("continental","stratum"),
                                        max_backcast = NULL,
                                        alternate_n = "n3") 


 indge = inds_gamye$data_summary
indge$decomp = "Including Year Effects"

indgnoe = inds_gamnoye$data_summary
indgnoe$decomp = "Excluding Year Effects"

indsall = rbind(indge,indgnoe)


source("colourblind safe qualitative pallete.r")



pdf(file = paste0(sp_dir,"all Fig 4.pdf"),
    width = 5,
    height = 4)

colye = c(model_pallete["gamye"],safe.pallet[[6]][5] )
names(colye) <- unique(indsall$decomp)


for(pp in unique(indsall$Region)){
  
  indcont = indsall[which(indsall$Region == pp),]  
  if(pp == "Continental"){
    dattt <- data.frame(Year = rep(indcont$Year,times = round(indcont$nrts/50)))
  }else{
    dattt <- data.frame(Year = rep(indcont$Year,times = indcont$nrts))
  }
  
  
  uylim = max(c(indcont$Index_q_0.975,indcont$obs_mean))*1.2
  
  labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
  labl_obs$label = "Observed mean counts"
  
  mxy = tapply(indcont[which(indcont$Year %in% c(1970:2013)),"Index"],indcont[which(indcont$Year %in% c(1970:2013)),"decomp"],max)
  
  labl_mods = indcont[which(indcont$Index %in% mxy),]

  cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = decomp))+
    theme_classic()+
    theme(legend.position = "none")+
    xlab(label = "")+
    ylab(label = "Predicted mean count")+
    labs(title = paste(pp))+
    geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
    coord_cartesian(ylim = c(0,uylim))+
    scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
    geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = decomp,colour = decomp), nudge_y = 0.075*uylim, nudge_x = 5)+
    scale_colour_manual(values = colye, aesthetics = c("colour","fill"))+

    geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = decomp),alpha = 0.2)+
    geom_line(aes(colour = decomp),size = 1.2)+
    geom_dotplot(data = dattt,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)
  
  print(cont_over)
  
  
}

dev.off()




# Rolling trends ----------------------------------------------------------



  inds_gamye <- generate_regional_indices(jags_mod = jags_mod_full,
                                          jags_data = jags_data,
                                          regions = c("continental","national"),
                                          max_backcast = NULL,
                                          alternate_n = "n") 
  
  inds_gamnoye <- generate_regional_indices(jags_mod = jags_mod_param,
                                            jags_data = jags_data,
                                            regions = c("continental","national"),
                                            max_backcast = NULL,
                                            alternate_n = "n3") 
  
  
  load(paste0(sp_dir,"slope/jags_mod_full.RData"))
  
  inds_slope <- generate_regional_indices(jags_mod = jags_mod_full,
                                            jags_data = jags_data,
                                            regions = c("continental","national"),
                                            max_backcast = NULL,
                                            alternate_n = "n") 
  
  
  

  
  # indge = inds_gamye$data_summary
  # indge$decomp = "GAMYE with Year Effects"
  # 
  # indgnoe = inds_gamnoye$data_summary
  # indgnoe$decomp = "GAMYE no Year Effects"
  # 
  # indsl = inds_slope$data_summary
  # indsl$decomp = "SLOPE"
  # 
  
  fy = min(jags_data$r_year)
  short_time = 10
  YYYY = max(jags_data$r_year)
  rollTrend = "Trend"
  
  
  c_orng = RColorBrewer::brewer.pal(9,"Set1")[5]
  c_red = RColorBrewer::brewer.pal(9,"Set1")[1]
  c_blue = RColorBrewer::brewer.pal(9,"Set1")[2]
  c_purp = RColorBrewer::brewer.pal(9,"Set1")[4]
  c_green = RColorBrewer::brewer.pal(9,"Set1")[3]
  
for(ly2 in c((fy+short_time):YYYY)){
  trst_gamye = generate_regional_trends(indices = inds_gamye,
                                  Min_year = ly2-short_time,
                                  Max_year = ly2,
                                  #quantiles = qs,
                                  slope = F,
                                  prob_decrease = c(0,25,30,50),
                                  prob_increase = c(0,33,100))
  trst_gamye$decomp <- "Including Year Effects"
  
  trst_gamnoye = generate_regional_trends(indices = inds_gamnoye,
                                        Min_year = ly2-short_time,
                                        Max_year = ly2,
                                        #quantiles = qs,
                                        slope = F,
                                        prob_decrease = c(0,25,30,50),
                                        prob_increase = c(0,33,100))
  trst_gamnoye$decomp <- "Excluding Year Effects"
  
  
  trst_sl = generate_regional_trends(indices = inds_slope,
                                        Min_year = ly2-short_time,
                                        Max_year = ly2,
                                        #quantiles = qs,
                                        slope = F,
                                        prob_decrease = c(0,25,30,50),
                                        prob_increase = c(0,33,100))
  trst_sl$decomp <- "SLOPE"
  
  
  
  if(ly2 == fy+short_time){
    tcos = rbind(trst_gamye,trst_gamnoye,trst_sl)
    
  }else{
    tcos = rbind(tcos,trst_gamye,trst_gamnoye,trst_sl)
  }
  
}

tcos$rolt = tcos[,rollTrend]
tcos$roltlci = tcos[,paste0(rollTrend,"_Q0.025")]
tcos$roltlci2 = tcos[,paste0(rollTrend,"_Q0.25")]
tcos$roltuci = tcos[,paste0(rollTrend,"_Q0.975")]
tcos$roltuci2 = tcos[,paste0(rollTrend,"_Q0.75")]



thresh30 = (0.7^(1/short_time)-1)*100
thresh50 = (0.5^(1/short_time)-1)*100

threshs = data.frame(thresh = c(thresh30,thresh50),
                     p_thresh = c(paste("-30% over",short_time,"years"),
                                  paste("-50% over",short_time,"years")),
                     Year = rep(min(tcos$End_year),2))


colye = c(model_pallete["gamye"],safe.pallet[[6]][5],model_pallete["slope"] )
names(colye) <- unique(tcos$decomp)


pdf(paste0(paste0(sp_dir,"sub Rolling_Trends.pdf")),
    width = 8.5,
    height = 6)
for(rg in unique(tcos$Region_alt)){
  
  tmp = tcos[which(tcos$Region_alt == rg),]
  
  tmpend4 = tmp[nrow(tmp)-4,]
  tmpend4$lab50 = "50% CI"
  tmpend4$lab95 = "95% CI"
  
  st_exc <- ""
  
  tmpend = tmp[nrow(tmp),]
  
  pth_30_labs = paste0(signif(100*tmpend[,"prob_decrease_30_percent"],2),"% probability of 30% decrease") 
  pth_50_labs = paste0(signif(100*tmpend[,"prob_decrease_50_percent"],2),"% probability of 50% decrease") 
  tmpend$pdec = paste(signif(tmpend[,"Percent_Change"],2),"% Change over",short_time,"years") 
  
  
  cpt = ggplot(data = tmp,aes(x = End_year,y = rolt,group = decomp,colour = decomp))+
    theme_classic()+
    theme(legend.position = "none")+
    labs(title = paste(species,"rolling",short_time,"year trends",rg,st_exc),
         subtitle = paste("Based on",rollTrend,"in",YYYY,":",pth_30_labs,"and",pth_50_labs))+
    xlab(paste("Ending year of",short_time,"trend"))+
    ylab(paste(short_time,"year trends"))+
    geom_hline(yintercept = thresh30,colour = c_orng)+
    geom_hline(yintercept = thresh50,colour = c_red)+
    geom_hline(yintercept = 0,colour = grey(0.5))+
     geom_label_repel(data = threshs,aes(x = Year,y = thresh,label = p_thresh), inherit.aes = F,position = "nudge")+
     geom_linerange(aes(x = End_year,ymin = roltlci,ymax = roltuci),alpha = 0.4,size = 0.9,position = position_dodge(width = 0.4))+
     geom_point(aes(x = End_year,y = rolt),size = 1,position = position_dodge(width = 0.4))+
    scale_colour_manual(values = colye, aesthetics = c("colour","fill"))
  
  
  
  
  ## update the theme?
  print(cpt)
  
}
dev.off()





pdf(paste0(paste0(sp_dir,"Figure 5.pdf")),
    width = 8.5,
    height = 6)
rg = "Canada"
  tmp = tcos[which(tcos$Region_alt == rg & tcos$End_year > 2000),]
  

  st_exc <- ""
  
  tmpend = tmp[which(tmp$End_year == 2011),]
  tmpend$lably = tmpend$roltlci
  tmpend[which(tmpend$decomp == "Excluding Year Effects"),"lably"] <- tmpend[which(tmpend$decomp == "Excluding Year Effects"),"roltuci"]
  

  
  cpt = ggplot(data = tmp,aes(x = End_year,y = rolt,group = decomp,colour = decomp))+
    theme_classic()+
    theme(legend.position = "none")+
    xlab(paste0("Year of Status Assessment (End of ",short_time,"-year trend)"))+
    ylab(paste0(short_time,"-year trends"))+
    geom_hline(yintercept = thresh30,colour = c_orng)+
    geom_hline(yintercept = thresh50,colour = c_red)+
    geom_hline(yintercept = 0,colour = grey(0.7))+
    annotate(geom = "text",x = 2005, label = threshs[1,"p_thresh"],y = threshs[1,"thresh"]-0.5,colour = c_orng)+
    annotate(geom = "text",x = 2005, label = threshs[2,"p_thresh"], y = threshs[2,"thresh"]-0.5,colour = c_red)+
    geom_linerange(aes(x = End_year,ymin = roltlci,ymax = roltuci),alpha = 0.4,size = 0.9,position = position_dodge(width = 0.4))+
    geom_point(aes(x = End_year,y = rolt),size = 1,position = position_dodge(width = 0.4))+
    geom_text_repel(data = tmpend,inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_x = -75)+
    scale_colour_manual(values = colye, aesthetics = c("colour","fill"))
  
  
  
  ## update the theme?
  print(cpt)
dev.off()

}



