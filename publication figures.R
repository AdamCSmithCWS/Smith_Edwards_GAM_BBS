

library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(lme4)
library(sf)
library(RColorBrewer)
library(facetscales)

load("comparison_summary_output.RData")


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


demo_sp <- c("Pine Siskin",
             "Horned Lark",
             "Carolina Wren",
             "American Kestrel",
             "Barn Swallow",
             "Wood Thrush",
             "Chestnut-collared Longspur",
             "Cooper's Hawk",
             "Ruby-throated Hummingbird",
             "Chimney Swift")







# colours -----------------------------------------------------------------


source("colourblind safe qualitative pallete.r")
species_pallete <- safe.pallet[[length(demo_sp)]] 
names(species_pallete) <- demo_sp

model_pallete <- safe.pallet[[length(models)]] 
model_pallete <- model_pallete[c(2,1,3,4)]
names(model_pallete) <- models


colye = c(model_pallete["gamye"],safe.pallet[[6]][5] )
names(colye) <- c("Including Year Effects","Smooth Only")


c_orng = brewer.pal(9,"Set1")[5]
c_red = brewer.pal(9,"Set1")[1]
c_blue = brewer.pal(9,"Set1")[2]
c_purp = brewer.pal(9,"Set1")[4]
c_green = brewer.pal(9,"Set1")[3]





# Figure 1 ----------------------------------------------------------------


species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

inds = tosave$indsout


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

pdf(file = paste0("Figures/Fig 1.pdf"),
    width = 5,
    height = 4)
print(cont_over)
dev.off()



# END figure 1 ------------------------------------------------------------







# Figure 2 ----------------------------------------------------------------


species = "Barn Swallow"
model = "gamye"


sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,model,"/jags_data.RData"))  

load(paste0(sp_dir,model,"/parameter_model_run.RData"))  
raw.dat = jags_mod_param$model$cluster1$data()
nknots = raw.dat$nknots
nstrata = raw.dat$nstrata
ndraws = length(jags_mod_param$sims.list$taubeta)


sumr = as.data.frame(jags_mod_param$summary)

B.X = sumr[paste0("B.X[",1:nknots,"]"),]

beta.X = sumr[grep(row.names(sumr),pattern = "beta.X",fixed = T),]


BXmat = jags_mod_param$sims.list[["B.X"]]
basis = raw.dat$X.basis
nyears = nrow(basis)
cont.sm = matrix(NA,ncol = nyears,nrow = ndraws)
for(y in 1:nyears){
  
  tmp = matrix(NA,nrow = ndraws,ncol = nknots)
  for(k in 1:nknots){
    tmp[,k] <- (BXmat[,k]*basis[y,k])
  }
  cont.sm[,y] <- exp(apply(tmp,1,sum))
  
}

conti = data.frame(year = 1:ncol(cont.sm)+1965,
                   med = apply(cont.sm,2,median),
                   lci = apply(cont.sm,2,quantile,probs = 0.025),
                   uci = apply(cont.sm,2,quantile,probs = 0.975))

betaxmat = jags_mod_param$sims.list[["beta.X"]]


for(s in 1:dim(betaxmat)[2]){
  strat.sm = cont.sm
  
  for(y in 1:nrow(basis)){
    tmp = matrix(NA,nrow = ndraws,ncol = nknots)
    for(k in 1:nknots){
      tmp[,k]  <- (betaxmat[,s,k]*basis[y,k])
    }
    strat.sm[,y] <- exp(apply(tmp,1,sum))
  }#y
  
  tt = data.frame(year = 1:ncol(strat.sm)+1965,
                  med = apply(strat.sm,2,median),
                  lci = apply(strat.sm,2,quantile,probs = 0.025),
                  uci = apply(strat.sm,2,quantile,probs = 0.975),
                  strat = s)    
  if(s == 1){
    strati = tt
  }else{
    strati = rbind(strati,tt)
  }
  
}

strati$strat = factor(strati$strat)


labls = conti[which(conti$year == 1980),]
labls$labl = "Survey-wide mean trajectory"

betaplot = ggplot(data = conti,aes(x = year,y = med))+
  theme_classic()+
  labs(title = "")+#paste0(species," GAM components including hyperparameter"))+
  ylab("Population change based on GAM smooth (linear scale)")+
  theme(legend.position = "none")+
  geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.075)+
  geom_text_repel(data = labls,aes(x = year,y = med,label = labl),nudge_x = 5,nudge_y = -0.75)+
  coord_cartesian(ylim = c(0,max(conti$med)*2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  geom_line(colour = grey(0.2),size = 1.4)

pdf(file = paste0("Figures/Fig 2.pdf"),
    width = 5,
    height = 4)
print(betaplot)
dev.off()



# END figure 2 ------------------------------------------------------------










# Figure 3 ----------------------------------------------------------------

species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

inds = tosave$indsout





indsel = inds[grepl(inds$Region,pattern = "-23",fixed = T),]  

indsel2 = indsel[which(indsel$Region == "US-IL-23"),]




# 
# allysc = list()
# length(allysc) = length(unique(indsel$Region))
# names(allysc) = unique(indsel$Region)
# 
# for(s in unique(indsel$Region)){
#   
#   tty = max(c(indsel[which(indsel$Region == s),"Index_q_0.975"],indsel[which(indsel$Region == s),"obs_mean"]))
#   allysc[[s]] <- scale_y_continuous(limits = c(0,tty))
#   
# }
# 



labl_obs = unique(indsel2[which(indsel2$Year == 2000 & indsel2$model == "gamye"),c("Year","obs_mean","Region")])
labl_obs$label = "Observed mean counts"

mxy = tapply(indsel2[which(indsel2$Year %in% c(1970:2000)),"Index"],indsel2[which(indsel2$Year %in% c(1970:2000)),"model"],max)

labl_mods = indsel2[which(indsel2$Index %in% mxy),]
labl_mods$Model = toupper(labl_mods$model)
labl_mods$Model[which(labl_mods$Model == "FIRSTDIFF")] <- "DIFFERENCE"



datt = tosave$datt
datt = filter(datt,grepl(datt$Stratum,pattern = "-23",fixed = T))
datt$Region = datt$Stratum



cont_over = ggplot(data = indsel[-which(indsel$model %in% models[c(3,4)]),],aes(x = Year,y = Index,group = model))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(label = "")+
  ylab(label = "Predicted mean count")+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_text_repel(data = labl_mods[4,],aes(x = Year,y = Index,label = Model),colour = grey(0.5), nudge_y = -10, nudge_x = -5)+
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label,group = Region),colour = grey(0.5),inherit.aes = F, nudge_y = -5, nudge_x = 1)+
  geom_text_repel(data = labl_mods[1:2,],aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 15, nudge_x = 2)+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
  geom_ribbon(data = indsel[which(indsel$model %in% models[c(4)]),],aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975),fill = grey(0.5),alpha = 0.2)+
  geom_line(data = indsel[which(indsel$model %in% models[c(4)]),],aes(x = Year,y = Index),colour = grey(0.7),size = 1.2)+
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
  geom_line(aes(colour = model),size = 1.2)+
  geom_dotplot(data = datt,mapping = aes(x = Year),drop = T,binaxis = "x", 
               stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.5)+
  facet_wrap(facets = ~Region,nrow = 2,ncol = 3, scales = "free")
  
  #facet_grid_sc(rows = vars(Region),scales = list(y = allysc))


pdf(file = paste0("Figures/Fig 3.pdf"),
    width = 8.5,
    height = 6)
print(cont_over)
dev.off()





# END figure 3 ------------------------------------------------------------




# Figure 4 ----------------------------------------------------------------
load("comparison_summary_output.RData")







sp4 = data.frame(species = demo_sp,
                 sp = c("AMKE",
                        "BARS",
                        "WOTH",
                        "CCLO",
                        "CAWA",
                        "CAWR",
                        "PISI",
                        "RTHU",
                        "CHSW"),
                 stringsAsFactors = F)

source("colourblind safe qualitative pallete.r")
species_pallete <- safe.pallet[[length(demo_sp)]] 
names(species_pallete) <- sp4$sp

fmod <- function(x){
  str_sub(x,start = 1,end = str_locate(x,pattern = " vs ")[,1]-1)
}
smod <- function(x){
  str_sub(x,start = str_locate(x,pattern = " vs ")[,2]+1,end = str_length(x))
}

dif_mod_year_over_out$M1 = paste("Favours",fmod(dif_mod_year_over_out$Contrast_full_name))
dif_mod_year_over_out$M2 = paste("Favours",smod(dif_mod_year_over_out$Contrast_full_name))




dif_mod_year_over_out2 = filter(dif_mod_year_over_out,Contrast_full_name %in% c("GAMYE vs SLOPE","GAMYE vs GAM"))#,"GAMYE vs DIFFERENCE"))
dif_mod_year_over_out2$Contrast_full_name = factor(dif_mod_year_over_out2$Contrast_full_name)

dif_mod_year_over_out2 = left_join(dif_mod_year_over_out2, sp4,by = "species")

dif_mod_year_over_outsplab = filter(dif_mod_year_over_out2,Contrast_full_name %in% c("GAMYE vs SLOPE"))

dif_mod_labs = filter(dif_mod_year_over_out2,species == demo_sp[1])
dif_mod_labs$M1loc = 0.02
dif_mod_labs$M2loc = -0.015

# overall summary graphs --------------------------------------------------

overall.comparison = ggplot(data = dif_mod_year_over_out2,aes(x = Contrast_full_name,y = mean,group = sp,colour = sp))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
   theme(axis.text.y = element_blank(),
         legend.position = "right",
         legend.text = element_text(size = 6))+
  #labs(title = paste("Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point(position = position_dodge(width = 0.4))+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),
                 alpha = 0.8,position = position_dodge(width = 0.4))+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M1loc,label = M1),show.legend = F,colour = grey(0.3),nudge_x = -0.5,size = 3)+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M2loc,label = M2),show.legend = F,colour = grey(0.3),nudge_x = -0.5,size = 3)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour"))+
  scale_y_continuous(limits = c(-0.03,0.04))+ # annotate(geom = "text",x = 0.5,y = 0.005,label = "Favours first")+
  # geom_label_repel(data = dif_mod_year_over_outsplab,
  #                 aes(x = Contrast_full_name,y = mean,group = species,colour = species,
  #                     label = sp),position = position_dodge(width = 0.4),
  #                 #ylim = c(0.02,0.05),
  #                 ylim = c(-0.05,-0.01),
  #                 size = 3,
  #                 segment.alpha = 0.3)+

  # annotate(geom = "text",x = 0.5,y = -0.005,label = "Favours second")+
  #scale_x_discrete(position = "top")+
  guides(colour = guide_legend(reverse = T))+
  coord_flip()


pdf(paste0("figures/Fig 4.pdf"),
    width = 5,
    height = 3)
print(overall.comparison)
dev.off()






# END Figure 4 ------------------------------------------------------------








# Figure 5 ----------------------------------------------------------------
species = "Barn Swallow"
model = "gamye"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,model,"/jags_data.RData"))  

  load(paste0(sp_dir,model,"/jags_mod_full.RData")) 
  indx2 = generate_regional_indices(jags_mod = jags_mod_full,
                                    jags_data = jags_data,
                                    #quantiles = qs,
                                    regions = c("continental"),
                                    max_backcast = NULL,
                                    alternate_n = "n")
  
  
  load(paste0(sp_dir,model,"/parameter_model_run.RData"))  
  jags_mod_full = jags_mod_param
  indx1 = generate_regional_indices(jags_mod = jags_mod_full,
                                    jags_data = jags_data,
                                    #quantiles = qs,
                                    regions = c("continental"),
                                    max_backcast = NULL,
                                    alternate_n = "n3")
  


fy = min(jags_data$r_year)
short_time = 10
YYYY = max(jags_data$r_year)
rollTrend = "Trend"






# plot the Smooth and full plot from GAMYE --------------------------------

indcont = indx1$data_summary
indcont$version = "Smooth Only"
ttind = indx2$data_summary
ttind$version = "Including Year Effects"
indcont = rbind(indcont,ttind)

indcont2 = indcont[which(indcont$version == "Including Year Effects"),]

uylim = max(c(indcont$Index_q_0.975,indcont$obs_mean))

labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
labl_obs$label = "Observed mean counts"

mxy = tapply(indcont[which(indcont$Year %in% c(1970:2013)),"Index"],indcont[which(indcont$Year %in% c(1970:2013)),"version"],max)

labl_mods = indcont[which(indcont$Index %in% mxy),]


ncby_y = table(jags_data$r_year)
ncby_y = ceiling(ncby_y/50)

dattc = data.frame(Year = rep(as.integer(names(ncby_y)),times = ncby_y))



cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = version))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(label = "")+
  ylab(label = "Predicted mean count")+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  coord_cartesian(ylim = c(0,uylim*1.2))+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  #geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = version),colour = grey(0.5), nudge_y = 0.075*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = version,colour = version), nudge_y = 0.075*uylim, nudge_x = 5)+
  scale_colour_manual(values = colye, aesthetics = c("colour","fill"))+

  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version),size = 1.2)+
  geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)

pdf(file = paste0("Figures/Fig 5.pdf"),
    width = 5,
    height = 4)
print(cont_over)
dev.off()



# End figure 5 ------------------------------------------------------------























# Figure 6 ----------------------------------------------------------------

species = "Wood Thrush"
model = "gamye"
sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,model,"/jags_data.RData"))  

load(paste0(sp_dir,model,"/jags_mod_full.RData")) 

load(paste0(sp_dir,model,"/parameter_model_run.RData"))  


inds_gamye <- generate_regional_indices(jags_mod = jags_mod_full,
                                        jags_data = jags_data,
                                        regions = c("national"),
                                        max_backcast = NULL,
                                        alternate_n = "n") 

inds_gamnoye <- generate_regional_indices(jags_mod = jags_mod_param,
                                          jags_data = jags_data,
                                          regions = c("national"),
                                          max_backcast = NULL,
                                          alternate_n = "n3") 


load(paste0(sp_dir,"slope/jags_mod_full.RData"))

inds_slope <- generate_regional_indices(jags_mod = jags_mod_full,
                                        jags_data = jags_data,
                                        regions = c("national"),
                                        max_backcast = NULL,
                                        alternate_n = "n") 




fy = 1990
short_time = 10
YYYY = max(jags_data$r_year)
rollTrend = "Trend"


for(ly2 in c((fy+short_time):YYYY)){
  trst_gamye = generate_regional_trends(indices = inds_gamye,
                                        Min_year = ly2-short_time,
                                        Max_year = ly2,
                                        #quantiles = qs,
                                        slope = F,
                                        prob_decrease = c(0,25,30,50),
                                        prob_increase = c(0,33,100))
  trst_gamye$decomp <- "GAMYE - Including Year Effects"
  
  trst_gamnoye = generate_regional_trends(indices = inds_gamnoye,
                                          Min_year = ly2-short_time,
                                          Max_year = ly2,
                                          #quantiles = qs,
                                          slope = F,
                                          prob_decrease = c(0,25,30,50),
                                          prob_increase = c(0,33,100))
  trst_gamnoye$decomp <- "GAMYE - Smooth only"
  
  
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



rg = "Canada"
tmp = tcos[which(tcos$Region_alt == rg & tcos$End_year > 2000),]


st_exc <- ""

tmpend = tmp[which( (tmp$End_year == 2011 & tmp$decomp == "SLOPE") |
                      (tmp$End_year == 2014 & tmp$decomp == "GAMYE - Smooth only") |
                      (tmp$End_year == 2005 & tmp$decomp == "GAMYE - Including Year Effects")  ),]
tmpend$lably = tmpend$rolt
sm.v = which(tmpend$decomp == "GAMYE - Smooth only")

#tmpend[sm.v,"lably"] <- tmpend[sm.v,"roltuci"]
sl.v = which(tmpend$decomp == "SLOPE")

tmpend[sl.v,"End_year"] <- tmpend[sl.v,"End_year"] +0.1
sy.v <- which(tmpend$decomp == "GAMYE - Including Year Effects")
tmpend[sy.v,"End_year"] <- tmpend[sy.v,"End_year"] - 0.15

threshplot = data.frame(roltlci = c(rep(thresh30,2),rep(thresh50,2)),
                        roltuci = c(rep(thresh50,2),rep(-15,2)),
                        year = rep(c(1999,2019),2),
                        cat = rep(c("Threatened","Endangered"),each = 2))
lylim = min(tmp$roltlci)-0.75
uylim = max(tmp$roltuci)+0.75


cpt = ggplot(data = tmp,aes(x = End_year,y = rolt,group = decomp,colour = decomp))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(paste0("Year of Status Assessment (End of ",short_time,"-year trend)"))+
  ylab(paste0(short_time,"-year trends"))+
  geom_hline(yintercept = thresh30,colour = c_orng,size = 1)+
  geom_hline(yintercept = thresh50,colour = c_red,size = 1)+
  coord_cartesian(ylim = c(lylim,uylim),xlim = c(2000,2018))+
  geom_ribbon(data = threshplot[1:2,],aes(x = year,ymin = roltlci,ymax = roltuci),fill = c_orng ,alpha = 0.1,inherit.aes = F)+
  geom_ribbon(data = threshplot[3:4,],aes(x = year,ymin = roltlci,ymax = roltuci),fill = c_red ,alpha = 0.1,inherit.aes = F)+
  geom_hline(yintercept = 0,colour = grey(0.7))+
  annotate(geom = "text",x = 2005, label = "Threatened",y = threshs[1,"thresh"]-0.3,colour = c_orng,size = 4,alpha = 1)+
  annotate(geom = "text",x = 2005, label = "Endangered", y = threshs[2,"thresh"]-0.3,colour = c_red,size = 4,alpha = 1)+
  geom_linerange(aes(x = End_year,ymin = roltlci,ymax = roltuci),alpha = 0.4,size = 0.2,position = position_dodge(width = 0.4))+
  geom_point(aes(x = End_year,y = rolt),size = 0.8,position = position_dodge(width = 0.4))+
  geom_line(aes(x = End_year,y = rolt),size = 0.2,position = position_dodge(width = 0.4),alpha = 0.4)+
  geom_text_repel(data = tmpend[sy.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_x = -0.5,nudge_y = +2)+
  geom_text_repel(data = tmpend[sl.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_x = 2,nudge_y = -1)+
  geom_text_repel(data = tmpend[sm.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_y = 4,nudge_x = +1.5)+
  scale_colour_manual(values = colye, aesthetics = c("colour","fill"))



## update the theme?
pdf(paste0(paste0("Figures/Fig 6.pdf")),
    width = 5,
    height = 4)
print(cpt)
dev.off()


# END Figure 6 -----------------------------------------------------------






# Figure 7 ----------------------------------------------------------------

species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

inds = tosave$indsout


datt = tosave$datt
ncby_y = table(datt$Year)
ncby_y = round(ncby_y/50)

dattc = data.frame(Year = rep(as.integer(names(ncby_y)),times = ncby_y))

indcont = inds[which(inds$Region == "Continental"),]  

indcont2 = indcont[which(indcont$model == "slope"),]

uylim = max(c(indcont$Index_q_0.975,indcont$obs_mean))

labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
labl_obs$label = "Observed mean counts"

# mxy = tapply(indcont[which(indcont$Year %in% c(1970:2013)),"Index"],indcont[which(indcont$Year %in% c(1970:2013)),"model"],max)

ww = which((indcont$model == "slope" & indcont$Year == 1970) |
             (indcont$model == "gamye" & indcont$Year == 1980) |
             (indcont$model == "gam" & indcont$Year == 1989) |
             (indcont$model == "firstdiff" & indcont$Year == 1983))
  

labl_mods = indcont[ww,]
labl_mods$Model = toupper(labl_mods$model)
labl_mods$Model[which(labl_mods$Model == "FIRSTDIFF")] <- "DIFFERENCE"

cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = model))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(label = "")+
  ylab(label = "Predicted mean count")+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  coord_cartesian(ylim = c(0,uylim))+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
   geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = 5)+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+

  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
  geom_line(aes(colour = model),size = 1.2)+
  geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)

pdf(file = paste0("Figures/Fig 7.pdf"),
    width = 5,
    height = 4)
print(cont_over)
dev.off()




# END Figure 7 ------------------------------------------------------------




# Figure 8 ----------------------------------------------------------------

species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

load(paste0(sp_dir,"saved objects4.RData"))
loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))


# extract the comparison model results and compile for graphing -----------

# comparison accounting for annual variation ------------------------------



for(comp in names(tosave2out)[1:2]){

  
  jags_data = tosave2out[[comp]]$m.both$model$data()
  m.both = tosave2out[[comp]]$m.both
  
  wparam = paste0("difmod_g1_full[",1:jags_data$ngroups1,"]")
  dif_mod_year1 = as.data.frame(m.both$summary[wparam,])
  names(dif_mod_year1)[3:7] <- c("lci","lqrt","med","uqrt","uci")
  
  dif_mod_year1$Contrast_name = comp
  dif_mod_year1$Year = (1:jags_data$ngroups1)+(2018-jags_data$ngroups1)
  dif_mod_year1$species = species
  dif_mod_year1$Contrast_full_name = contrast_full_names[dif_mod_year1$Contrast_name]

  if(comp == "gamye_firstdiff"){
    dif_mod_year = dif_mod_year1
  }else{
    
    dif_mod_year = rbind(dif_mod_year,dif_mod_year1)
  }
  
  }#end compe


ww = which((dif_mod_year$Year == 1970 & dif_mod_year$Contrast_name == "gamye_firstdiff")|
             (dif_mod_year$Year == 1971 & dif_mod_year$Contrast_name == "gamye_slope"))

lbl = dif_mod_year[ww,]
lbl[which(lbl$Contrast_name == "gamye_firstdiff"),"Year"] = lbl[which(lbl$Contrast_name == "gamye_firstdiff"),"Year"]-0.25
lbl[which(lbl$Contrast_name == "gamye_slope"),"Year"] = lbl[which(lbl$Contrast_name == "gamye_slope"),"Year"]+0.25



an_contr = ggplot(data = dif_mod_year,aes(x = Year,y = mean,group = Contrast_name,colour = Contrast_name))+
  geom_point(position = position_dodge(width = 0.75))+
  geom_linerange(aes(x = Year,ymin = lci,ymax = uci),position = position_dodge(width = 0.75),alpha = 0.5)+
  #coord_cartesian(ylim = c(-0.03,0.03))+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(title = "")+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  scale_x_continuous(breaks = c(seq(1970,2010,by = 10),2018))+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_text_repel(data = lbl[1,],aes(x = Year,y = lci,label = Contrast_full_name),nudge_y = -0.005)+
  geom_text_repel(data = lbl[2,],aes(x = Year,y = uci,label = Contrast_full_name),nudge_y = +0.005,nudge_x = 8)+
  annotate(geom = "text",x = 2001,y = 0.04,label = "Positive favours GAMYE",size = 3.5,colour = grey(0.6))+
  annotate(geom = "text",x = 2001,y = -0.012,label = "Negative favours alternate",size = 3.5,colour = grey(0.6))  


pdf(paste0("Figures/Fig 8.pdf"),
    width = 5,
    height = 4)
print(an_contr)
dev.off()




# END Figure 8 ------------------------------------------------------------




# Figure 9 ----------------------------------------------------------------


species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

load(paste0(sp_dir,"saved objects4.RData"))
loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))


# extract the comparison model results and compile for graphing -----------

# comparison accounting for annual variation ------------------------------


for(comp in names(tosave2out)[1:2]){
  
  jags_data = tosave2out[[comp]]$m.both$model$data()
  m.both = tosave2out[[comp]]$m.both
  
  
  
  wparam = paste0("difmod_g2_full[",1:jags_data$ngroups2,"]")
  dif_mod_strat1 = as.data.frame(m.both$summary[wparam,])
  names(dif_mod_strat1)[3:7] <- c("lci","lqrt","med","uqrt","uci")
  
  
  
  dif_mod_strat1$Contrast_name = comp
  dif_mod_strat1$Stratum_Factored = (1:jags_data$ngroups2)
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



mp = read_sf(dsn = "map/BBS_USGS_strata.shp")


for(i in c("A","B")){
  if(i == "A"){
    dif_s_p = filter(dif_mod_strat,Contrast_name %in% c("gamye_firstdiff") & species == species)
    modcomp_pallete <- model_pallete[c("firstdiff","gamye")]
    
  }else{
    dif_s_p = filter(dif_mod_strat,Contrast_name %in% c("gamye_slope") & species == species)
    modcomp_pallete <- model_pallete[c("slope","gamye")]
    
      }

    dif_s_p$ST_12 = as.character(dif_s_p$Stratum)
    
    dif_s_p$cl[which(dif_s_p$mean > 0)] = "gamye"
    
    if(i == "A"){   
      dif_s_p$cl[which(dif_s_p$mean < 0)] = "firstdiff"
      
    dif_s_p$cl = factor(dif_s_p$cl ,levels = c("firstdiff","gamye"))
    }else{
      dif_s_p$cl[which(dif_s_p$mean < 0)] = "slope"
      
      dif_s_p$cl = factor(dif_s_p$cl ,levels = c("slope","gamye"))
    }
  
    
    names(modcomp_pallete) <- levels(dif_s_p$cl)
    
    
    mp2 = left_join(x = mp,y = dif_s_p,by = "ST_12")
    mp2gamye = filter(mp2, ST_12 %in% c("US-NC-27"))
    mp2gamye$lbl = "GAMYE"
    mp2diff = filter(mp2, ST_12 %in% c("US-CA-32"))
    mp2diff$lbl = "DIFFERENCE"
    
    mp.plot = ggplot()+
      geom_sf(data = mp2,aes(fill = cl),colour = grey(0.4),size = 0.1)+
      theme_minimal()+
      ylab("")+
      xlab("")+
      theme(legend.position = "none", line = element_line(size = 0.4), rect = element_rect(size = 0.1),
            axis.text = element_blank(),axis.line = element_blank())+
      geom_sf_text(data = mp2gamye,aes(label = lbl),colour = modcomp_pallete[2],hjust = -0.2)+
      geom_sf_text(data = mp2diff,aes(label = lbl),colour = modcomp_pallete[1],hjust = 1.1)+
      #geom_sf_text(data = mp2,aes(label = ST_12))+
      # annotate(geom = "text",y = 50,x = -155,label = "GAMYE",colour = model_pallete["gamye"])+
      # annotate(geom = "text",y = 50,x = -155,label = "DIFFERENCE",colour = model_pallete["firstdiff"])+
      # #labs(title = paste(spp,comp))+
      #guides(fill = guide_legend(reverse = T))+
      scale_colour_manual(values = modcomp_pallete, aesthetics = c("fill"))

      
      pdf(file = paste0("Figures/Fig 9",i,".pdf"),
          width = 4,
          height = 4)
    
    
    print(mp.plot)
    

  dev.off()

}



# END Figure 9 ------------------------------------------------------------












