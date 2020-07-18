# This script produces all of the figures in the manuscript, and provides the .RData file
# used to generate the supplemental material figures in Supplmental Material Figures.Rmd



library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(ggthemes)
library(tidyverse)
library(lme4)
library(sf)
library(RColorBrewer)
library(facetscales)
library(patchwork)

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


# demo_sp <- c("Pine Siskin",
#              "Canada Warbler",
#              "Carolina Wren",
#              "American Kestrel",
#              "Barn Swallow",
#              "Wood Thrush",
#              "Chestnut-collared Longspur",
#              "Cooper's Hawk",
#              "Ruby-throated Hummingbird",
#              "Chimney Swift")









# colours -----------------------------------------------------------------


#source("colourblind safe qualitative pallete.r")
# species_pallete <- ggthemes::colorblind_pal()
# species_pallete <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
#                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#species_pallete <- RColorBrewer::brewer.pal(length(demo_sp),"Paired")
#species_pallete <- safe.pallet[[length(demo_sp)]] 


species_pallete <- rep(viridis::viridis(5),each = 2)
scales::show_col(species_pallete)

names(species_pallete) <- demo_sp




model_pallete <- viridis::viridis(length(models))
model_pallete <- model_pallete[c(2,4,3,1)]
names(model_pallete) <- models
scales::show_col(model_pallete)


colye = c(model_pallete["gamye"],viridis::viridis(7)[6] )
names(colye) <- c("Including Year Effects","Smooth Only")
scales::show_col(colye)


c_orng = viridis::inferno(7)[6]
c_red = viridis::inferno(7)[5]
#scales::show_col( viridis::inferno(7))




# run time summary --------------------------------------------------------

runtime = data.frame(species = rep(demo_sp,each = length(models)),
                     model = rep(models,times = length(demo_sp)))
for(species in demo_sp){

sp_dir = paste0("output/",species,"/")

for(model in models){
  mod_dir = paste0(sp_dir,model,"/")
  
  load(paste0(mod_dir,"jags_mod_full.RData"))
  ww = which(runtime$species == species & runtime$model == model)
  
  runtime[ww,"minutes"] <- jags_mod_full$mcmc.info$elapsed.mins
  runtime[ww,"iterations"] <- jags_mod_full$mcmc.info$n.iter
  runtime[ww,"burnin"] <- jags_mod_full$mcmc.info$n.burnin
  runtime[ww,"samples"] <- jags_mod_full$mcmc.info$n.samples
  runtime[ww,"thin"] <- jags_mod_full$mcmc.info$n.thin
  
  
}


}#species

runtime$hours = round(runtime$minutes/60,1)
#labl_mods = runtime[which(runtime$species == "Ruby-throated Hummingbird"),]
runtime$MOD = toupper(runtime$model)

runtime$MOD[which(runtime$MOD == "FIRSTDIFF")] <- "DIFFERENCE"
model_pallete2 = model_pallete

names(model_pallete2) <- toupper(names(model_pallete))
names(model_pallete2)[which(names(model_pallete2) == "FIRSTDIFF")] <- "DIFFERENCE"


runp = ggplot(data = runtime,aes(x = species,y = hours,colour = MOD,fill = MOD,group = MOD))+
  geom_bar(position = "dodge",stat = "identity")+
  scale_colour_manual(name = "",
                      values = model_pallete2, aesthetics = c("colour","fill"),
                      guide = guide_legend(reverse = TRUE,nrow = 1,
                                           label.hjust = -1,
                                           label.theme = element_text(size = 8)))+
  xlab("")+
  ylab("Hours to run 20 000 iterations")+
  coord_flip()+
  theme_classic()+
  theme(legend.position = "bottom",legend.key.width = grid::unit(0.2,units = "inches"))
  
  
pdf(file = paste0("Figures/supplement/Fig 10alt.pdf"),
    width = 5,
    height = 7)
print(runp)
dev.off()

save(list = "runp",file = "Figures/supplement/Fig 10alt.RData")


p.min = function(x){round(x/x[4],2)}

pch <- runtime %>% 
  group_by(species) %>% 
  mutate(pmin = p.min(minutes))

tapply(pch$pmin,pch$model,mean)

write.csv(pch,"figures/supplement/runtime summary across models and species.csv")







# Figure 1 ----------------------------------------------------------------


svplots = list()
length(svplots) = 6
names(svplots) = c("Barn Swallow",
                   "Wood Thrush",
                   "Cooper's Hawk",
                   "Carolina Wren",
                   "Ruby-throated Hummingbird",
                   "Chimney Swift")

species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

inds = tosave$indsout


datt = tosave$datt
ncby_y = table(datt$Year)
nstrat_scl = max(datt$Stratum_Factored)
ncby_y = ceiling(ncby_y/50)

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
             (indcont$model == "firstdiff" & indcont$Year == 1994))


labl_mods = indcont[ww,]
labl_mods$Model = toupper(labl_mods$model)
labl_mods$Model[which(labl_mods$Model == "FIRSTDIFF")] <- "DIFFERENCE"

cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = model))+
  theme_classic()+
  labs(title = species)+
  theme(legend.position = "none")+
  xlab(label = "")+
  ylab(label = "Predicted mean count")+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  coord_cartesian(ylim = c(0,uylim))+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  annotate(geom = "text",x = 1985,y = uylim*0.05,label = "Routes/50",colour = grey(0.6)) +
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = 5)+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
  
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
  geom_line(aes(colour = model),size = 1.2)+
  geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)


svplots[[species]] <- cont_over

for(species in names(svplots)[-1]){
  
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
  
  inds = tosave$indsout
  
  
  datt = tosave$datt
  ncby_y = table(datt$Year)
  nstrat_scl = max(datt$Stratum_Factored)
  ncby_y = ceiling(ncby_y/50)
  
  dattc = data.frame(Year = rep(as.integer(names(ncby_y)),times = ncby_y))
  
  indcont = inds[which(inds$Region == "Continental"),]  
  
  indcont2 = indcont[which(indcont$model == "slope"),]
  
  uylim = min(c(max(c(indcont$Index_q_0.975,indcont$obs_mean)),max(indcont$Index * 3)))
  
  labl_obs = unique(indcont[which(indcont$Year == 1970),c("Year","obs_mean")])
  labl_obs$label = "Observed mean counts"
  
  # 
  
  cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = model))+
    theme_classic()+
    theme(legend.position = "none")+
    xlab(label = "")+
    ylab(label = "Predicted mean count")+
    labs(title = species)+
    geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
    coord_cartesian(ylim = c(0,uylim))+
    scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    #geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
    #geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = runif(1,-3,3))+
    scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
    
    geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.15)+
    geom_line(aes(colour = model),size = 1.2)+
    geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)
  
  svplots[[species]] <- cont_over
 
}



pdf(file = paste0("Figures/Fig 1.pdf"),
    width = 7,
    height = 9.5)
print(svplots[[1]] + svplots[[2]] + svplots[[3]] + svplots[[4]] + svplots[[5]] + svplots[[6]] + plot_layout(ncol = 2))
dev.off()


# species versions for supplement ------------------------------------------

svplots = list()
length(svplots) = length(demo_sp)
names(svplots) = demo_sp


pdf(file = paste0("Figures/supplement/Fig 1 all species.pdf"),
    width = 5,
    height = 4)

for(species in demo_sp){
  
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
  
  inds = tosave$indsout
  
  
  
  datt = tosave$datt
  ncby_y = table(datt$Year)
  nstrat_scl = max(datt$Stratum_Factored)
  ncby_y = ceiling(ncby_y/50)
  
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
               (indcont$model == "firstdiff" & indcont$Year == 1986))
  
  
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
    annotate(geom = "text",x = 1985,y = uylim*0.05,label = "Routes/50",colour = grey(0.6)) +
    geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
    geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 0.075*uylim, nudge_x = 5)+
    scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
    
    geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
    geom_line(aes(colour = model),size = 1.2)+
    geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.5)
  
  print(cont_over)
  
  svplots[[species]] <- cont_over
  
}

dev.off()

save(list = "svplots",file = "Figures/supplement/Fig 1 all species.RData")





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


labls = conti[which(conti$year == 1985),]
labls$labl = "Survey-wide mean trajectory"

betaplot = ggplot(data = conti,aes(x = year,y = med))+
  theme_classic()+
  labs(title = "")+#paste0(species," GAM components including hyperparameter"))+
  ylab("Trajectory from GAM smooth (linear scale)")+
  theme(legend.position = "none")+
  geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.075)+
  geom_text_repel(data = labls,aes(x = year,y = med,label = labl),nudge_x = 5,nudge_y = -0.85)+
  coord_cartesian(ylim = c(0,max(conti$med)*2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  geom_line(colour = grey(0.2),size = 1.4)

pdf(file = paste0("Figures/Fig 2.pdf"),
    width = 3.5,
    height = 3)
print(betaplot)
dev.off()




# species versions for supplement -----------------------------------------



svplots = list()
length(svplots) = length(demo_sp)
names(svplots) = demo_sp


pdf(file = paste0("Figures/supplement/Fig 2 all species.pdf"),
    width = 5,
    height = 4)
model = "gamye"

for(species in demo_sp){

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
  labs(title = species)+#paste0(species," GAM components including hyperparameter"))+
  ylab("Population change based on GAM smooth (linear scale)")+
  theme(legend.position = "none")+
  geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.075)+
  geom_text_repel(data = labls,aes(x = year,y = med,label = labl),nudge_x = 5,nudge_y = -0.75)+
  coord_cartesian(ylim = c(0,max(conti$med)*2))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  geom_line(colour = grey(0.2),size = 1.4)
print(betaplot)
svplots[[species]] <- betaplot

}

dev.off()

save(list = "svplots",file = "Figures/supplement/Fig 2 all species.RData")










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

mxy = tapply(indsel2[which(indsel2$Year %in% c(1970:2006)),"Index"],indsel2[which(indsel2$Year %in% c(1970:2006)),"model"],max)

wmxy = which(indsel2$Index %in% mxy)

names(wmxy) <- indsel2[wmxy,"model"]
wmxy["slope"] <- which(indsel2$Year == 1990 & indsel2$model == "slope")

labl_mods = indsel2[wmxy,]


labl_mods$Model = toupper(labl_mods$model)
labl_mods$Model[which(labl_mods$Model == "FIRSTDIFF")] <- "DIFFERENCE"



datt = tosave$datt
datt = filter(datt,grepl(datt$Stratum,pattern = "-23",fixed = T))
datt$Region = datt$Stratum



cont_over = ggplot(data = indsel[-which(indsel$model %in% models[c(2)]),],aes(x = Year,y = Index,group = model))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(label = "")+
  ylab(label = "Predicted mean count")+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = model),alpha = 0.2)+
  geom_line(aes(colour = model),size = 1.2)+
  geom_text_repel(data = labl_mods[c(3,4),],aes(x = Year,y = Index,label = Model,colour = model), nudge_y = -10, nudge_x = -5)+
 # geom_text_repel(data = labl_mods[3,],aes(x = Year,y = Index,label = Model),colour = grey(0.5), nudge_y = -10, nudge_x = -5)+
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label,group = Region),colour = grey(0.5),inherit.aes = F, nudge_y = -5, nudge_x = 1)+
  geom_text_repel(data = labl_mods[1,],aes(x = Year,y = Index,label = Model,colour = model), nudge_y = 15, nudge_x = 2)+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))+
  #geom_ribbon(data = indsel[which(indsel$model %in% models[c(4)]),],aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975),fill = grey(0.5),alpha = 0.2)+
  #geom_line(data = indsel[which(indsel$model %in% models[c(4)]),],aes(x = Year,y = Index),colour = grey(0.7),size = 1.2)+
 geom_dotplot(data = datt,mapping = aes(x = Year),drop = T,binaxis = "x", 
               stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.5)+
  facet_wrap(facets = ~Region,nrow = 2,ncol = 3, scales = "free")
  
  #facet_grid_sc(rows = vars(Region),scales = list(y = allysc))


pdf(file = paste0("Figures/Fig 3.pdf"),
    width = 7,
    height = 6)
print(cont_over)
dev.off()





# Figure 4 --------------------------------------------------------------

species = "Cooper's Hawk"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

inds = tosave$indsout

tr = tosave$trendsout

trS = filter(tr,Region_type == "stratum")
trS$Model = toupper(trS$model)
trS$Model[which(trS$Model == "FIRSTDIFF")] <- "DIFFERENCE"

trS$abs_Trend = abs(trS$Trend)
#trS$log_Number_Routes = log(trS$Mean_Number_of_Routes,base = 10)
#trS$log_Number_Routes = log(trS$Relative_Abundance*trS$Mean_Number_of_Routes,base = 10)

#trS$data_cat = cut(trS$log_Number_Routes,breaks = c(0,1,2,3,4,5))
trpl = ggplot(data = trS,aes(x = Mean_Number_of_Routes,y = abs_Trend,colour = model))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size = 6))+
  ylab("Absolute value of long-term trend (%/year)")+
  xlab("Mean number of routes in the stratum/year")+
  scale_x_continuous(trans = "log", breaks = (c(1,3,10,25,50)))+
  scale_y_continuous(limits = c(0,NA),expand = c(0,0))+
  facet_wrap(facets = ~Model,nrow = 2)+
  #geom_smooth(method = "lm")+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))
  

# trpl2 = ggplot(data = trS,aes(x = data_cat,y = abs_Trend,colour = model))+
#   geom_boxplot()+
#   facet_wrap(facets = ~model,nrow = 2)
pdf(paste0("figures/Fig 4.pdf"),
    width = 7,
    height = 5)
print(trpl)
dev.off()


# species versions of Figure 4 --------------------------------------------------------------
svplots = list()
length(svplots) = length(demo_sp)
names(svplots) = demo_sp


pdf(file = paste0("Figures/supplement/Fig 4 all species.pdf"),
    width = 5,
    height = 4)

for(species in demo_sp){

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

inds = tosave$indsout

tr = tosave$trendsout

trS = filter(tr,Region_type == "stratum")
trS$Model = toupper(trS$model)
trS$Model[which(trS$Model == "FIRSTDIFF")] <- "DIFFERENCE"

trS$abs_Trend = abs(trS$Trend)
trS$log_Number_Routes = log(trS$Mean_Number_of_Routes,base = 10)
#trS$log_Number_Routes = log(trS$Relative_Abundance*trS$Mean_Number_of_Routes,base = 10)

#trS$data_cat = cut(trS$log_Number_Routes,breaks = c(0,1,2,3,4,5))
trpl = ggplot(data = trS,aes(x = Mean_Number_of_Routes,y = abs_Trend,colour = model))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size = 6))+
  ylab("Absolute value of long-term trend (%/year)")+
  xlab("Mean number of routes in the stratum/year")+
  scale_x_continuous(trans = "log", breaks = (c(1,3,10,25,50)))+
  scale_y_continuous(limits = c(0,NA),expand = c(0,0))+
  facet_wrap(facets = ~Model,nrow = 2)+
  #geom_smooth(method = "lm")+
  scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))


print(trpl)
svplots[[species]] <- trpl

}

dev.off()

save(list = "svplots",file = "Figures/supplement/Fig 4 all species.RData")


# Figure 5 ----------------------------------------------------------------




demo_sp = sort(demo_sp)
sp = c("AMKE",
       "BARS",
       "CAWA",
       "CAWR",
       "CCLO",
       "CHSW",
       "COHA",
       "PISI",
       "RTHU",
       "WOTH")

sp4 = data.frame(species = demo_sp,
                 sp = sp,
                 stringsAsFactors = F)

pchs = rep(c(16,17),times = 5)
names(pchs) = sp4$sp

# source("colourblind safe qualitative pallete.r")
# species_pallete <- safe.pallet[[length(demo_sp)]] 
names(species_pallete) <- sp4$sp

fmod <- function(x){
  str_sub(x,start = 1,end = str_locate(x,pattern = " vs ")[,1]-1)
}
smod <- function(x){
  str_sub(x,start = str_locate(x,pattern = " vs ")[,2]+1,end = str_length(x))
}

dif_mod_year_over_out$M1 = paste("Favours",fmod(dif_mod_year_over_out$Contrast_full_name))
dif_mod_year_over_out$M2 = paste("Favours",smod(dif_mod_year_over_out$Contrast_full_name))




dif_mod_year_over_out2 = filter(dif_mod_year_over_out,Contrast_full_name %in% c("GAMYE vs SLOPE","GAMYE vs GAM","GAMYE vs DIFFERENCE"))#,"GAM vs DIFFERENCE"))#,"GAMYE vs DIFFERENCE"))


dif_mod_year_over_out2$Contrast_full_name = factor(dif_mod_year_over_out2$Contrast_full_name)

dif_mod_year_over_out2 = left_join(dif_mod_year_over_out2, sp4,by = "species")

dif_mod_year_over_outsplab = filter(dif_mod_year_over_out2,Contrast_full_name %in% c("GAMYE vs SLOPE"))

dif_mod_labs = filter(dif_mod_year_over_out2,species == demo_sp[1])
dif_mod_labs$M1loc = 0.02
dif_mod_labs$M2loc = -0.02

dif_mod_border = dif_mod_labs
dif_mod_border$Contrast_full_name = as.integer(dif_mod_border$Contrast_full_name)-0.4
dif_mod_border$M1loc = dif_mod_border$M1loc+0.0175
dif_mod_border$M2loc = dif_mod_border$M2loc-0.015

# overall summary graphs --------------------------------------------------

gl <- guide_legend(title = "Species",reverse = T,
                   nrow = 2,title.position = "top",
                   label.hjust = 0)

overall.comparison = ggplot(data = dif_mod_year_over_out2,aes(x = Contrast_full_name,y = mean,group = sp,colour = sp,fill = sp))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
   theme(axis.text.y = element_blank(),
         legend.position = "bottom",
         legend.text = element_text(size = 6))+
  #guide_legend(nrow = 2)+
  #labs(title = paste("Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point(aes(shape = sp),position = position_dodge(width = 0.4))+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),
                 alpha = 0.8,position = position_dodge(width = 0.4))+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M1loc,label = M1),show.legend = F,colour = grey(0.3),nudge_x = -0.3,size = 3)+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M2loc,label = M2),show.legend = F,colour = grey(0.3),nudge_x = -0.3,size = 3)+
  geom_linerange(inherit.aes = F,data = dif_mod_border,aes(x = Contrast_full_name,ymin = M2loc,ymax = M1loc),show.legend = F,colour = grey(0.5),size = 0.5)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour","fill"),name = "Species")+
  scale_shape_manual(values = pchs, name = "Species")+
  scale_y_continuous(limits = c(-0.04,0.04))+ # annotate(geom = "text",x = 0.5,y = 0.005,label = "Favours first")+
  # geom_label_repel(data = dif_mod_year_over_outsplab,
  #                 aes(x = Contrast_full_name,y = mean,group = species,colour = species,
  #                     label = sp),position = position_dodge(width = 0.4),
  #                 #ylim = c(0.02,0.05),
  #                 ylim = c(-0.05,-0.01),
  #                 size = 3,
  #                 segment.alpha = 0.3)+

  # annotate(geom = "text",x = 0.5,y = -0.005,label = "Favours second")+
  #scale_x_discrete(position = "top")+
  guides(colour = gl,fill = gl,shape = gl)+
  coord_flip()


pdf(paste0("figures/Fig 5.pdf"),
    width = 3.5,
    height = 5)
print(overall.comparison)
dev.off()








# model versions Figure 5 all models supplement ---------------------------------


dif_mod_year_over_out2 = dif_mod_year_over_out


dif_mod_year_over_out2$Contrast_full_name = factor(dif_mod_year_over_out2$Contrast_full_name)

dif_mod_year_over_out2 = left_join(dif_mod_year_over_out2, sp4,by = "species")

dif_mod_year_over_outsplab = filter(dif_mod_year_over_out2,Contrast_full_name %in% c("GAMYE vs SLOPE"))

dif_mod_labs = filter(dif_mod_year_over_out2,species == demo_sp[1])
dif_mod_labs$M1loc = 0.02
dif_mod_labs$M2loc = -0.02

dif_mod_border = dif_mod_labs
dif_mod_border$Contrast_full_name = as.integer(dif_mod_border$Contrast_full_name)-0.45
dif_mod_border$M1loc = dif_mod_border$M1loc+0.01
dif_mod_border$M2loc = dif_mod_border$M2loc-0.01

gl <- guide_legend("Species",reverse = T)

overall.comparison = ggplot(data = dif_mod_year_over_out2,aes(x = Contrast_full_name,y = mean,group = sp,colour = sp,fill = sp))+
  #coord_cartesian(ylim = c(-0.05,0.05))+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 6))+
  #labs(title = paste("Overall cross validation comparison"))+
  ylab("Mean difference in point-wise log-probability")+
  xlab("")+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_point(aes(shape = sp),position = position_dodge(width = 0.4))+
  geom_linerange(aes(x = Contrast_full_name,ymin = lci,ymax = uci),
                 alpha = 0.8,position = position_dodge(width = 0.4))+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M1loc,label = M1),show.legend = F,colour = grey(0.3),nudge_x = -0.35,size = 3)+
  geom_text(inherit.aes = F,data = dif_mod_labs,aes(x = Contrast_full_name,y = M2loc,label = M2),show.legend = F,colour = grey(0.3),nudge_x = -0.35,size = 3)+
  geom_linerange(inherit.aes = F,data = dif_mod_border,aes(x = Contrast_full_name,ymin = M2loc,ymax = M1loc),show.legend = F,colour = grey(0.5),size = 0.5)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour","fill"),name = "Species")+
  scale_shape_manual(values = pchs, name = "Species")+
  #scale_y_continuous(limits = c(-0.03,0.04))+ # annotate(geom = "text",x = 0.5,y = 0.005,label = "Favours first")+
  # geom_label_repel(data = dif_mod_year_over_outsplab,
  #                 aes(x = Contrast_full_name,y = mean,group = species,colour = species,
  #                     label = sp),position = position_dodge(width = 0.4),
  #                 #ylim = c(0.02,0.05),
  #                 ylim = c(-0.05,-0.01),
  #                 size = 3,
  #                 segment.alpha = 0.3)+
  
  # annotate(geom = "text",x = 0.5,y = -0.005,label = "Favours second")+
  #scale_x_discrete(position = "top")+
  guides(colour = gl,fill = gl,shape = gl)+
  coord_flip()


pdf(paste0("figures/supplement/Fig 5 all models.pdf"),
    width = 5,
    height = 7)
print(overall.comparison)
dev.off()


save(list = "overall.comparison",file = "Figures/supplement/Fig 5 all models.RData")





#

# Figure 6 ----------------------------------------------------------------
species = "Carolina Wren"
model = "gamye"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,model,"/jags_data.RData"))  

  load(paste0(sp_dir,model,"/jags_mod_full.RData")) 
  indx2 = generate_indices(jags_mod = jags_mod_full,
                                    jags_data = jags_data,
                                    #quantiles = qs,
                                    regions = c("continental"),
                                    max_backcast = NULL,
                                    alternate_n = "n")
  
  
  load(paste0(sp_dir,model,"/parameter_model_run.RData"))  
  jags_mod_full = jags_mod_param
  indx1 = generate_indices(jags_mod = jags_mod_full,
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

mxy = which((indcont$Year == c(1992) & indcont$version == "Smooth Only")| (indcont$Year == c(2005) & indcont$version == "Including Year Effects"))

labl_mods = indcont[mxy,]


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
  geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = version,colour = version), nudge_y = 0.15*uylim, nudge_x = -5)+
  scale_colour_manual(values = colye, aesthetics = c("colour","fill"))+

  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version),size = 1.2)+
  geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)

pdf(file = paste0("Figures/Fig 6.pdf"),
    width = 3.5,
    height = 3)
print(cont_over)
dev.off()






# supplement versions for all species -------------------------------------


svplots = list()
length(svplots) = length(demo_sp)
names(svplots) = demo_sp

pdf(file = paste0("Figures/supplement/Fig 6 all species.pdf"),
    width = 5,
    height = 4)

for(species in demo_sp){
model = "gamye"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,model,"/jags_data.RData"))  

load(paste0(sp_dir,model,"/jags_mod_full.RData")) 
indx2 = generate_indices(jags_mod = jags_mod_full,
                                  jags_data = jags_data,
                                  #quantiles = qs,
                                  regions = c("continental"),
                                  max_backcast = NULL,
                                  alternate_n = "n")


load(paste0(sp_dir,model,"/parameter_model_run.RData"))  
jags_mod_full = jags_mod_param
indx1 = generate_indices(jags_mod = jags_mod_full,
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
  labs(title = species)+
  geom_point(aes(x = Year,y = obs_mean),colour = grey(0.8),size = 0.8)+
  coord_cartesian(ylim = c(0,uylim*1.2))+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  #geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = version),colour = grey(0.5), nudge_y = 0.075*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_obs,aes(x = Year,y = obs_mean,label = label),colour = grey(0.5),inherit.aes = F, nudge_y = -0.1*uylim, nudge_x = 5)+
  geom_text_repel(data = labl_mods,aes(x = Year,y = Index,label = version,colour = version), nudge_y = 0.075*uylim, nudge_x = 5)+
  scale_colour_manual(values = colye, aesthetics = c("colour","fill"))+
  
  geom_ribbon(aes(x = Year,ymin = Index_q_0.025,ymax = Index_q_0.975,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version),size = 1)+
  geom_dotplot(data = dattc,mapping = aes(x = Year),drop = T,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = F,fill = grey(0.6),colour = grey(0.6),alpha = 0.2,dotsize = 0.3)
print(cont_over)

svplots[[species]] <- cont_over

}

dev.off()

save(list = "svplots",file = "Figures/supplement/Fig 6 all species.RData")























# Figure 7 ----------------------------------------------------------------



for(species in demo_sp){
  model = "gamye"
  
  
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,model,"/jags_data.RData"))  
  
  load(paste0(sp_dir,model,"/jags_mod_full.RData")) 
  indx = generate_indices(jags_mod = jags_mod_full,
                                   jags_data = jags_data,
                                   #quantiles = qs,
                                   regions = c("continental"),
                                   max_backcast = NULL,
                                   alternate_n = "n")
  fy = min(jags_data$r_year)
  short_time = 10
  YYYY = max(jags_data$r_year)
  rollTrend = "Trend"
  
  
  for(y in (fy+short_time):YYYY){
    if(y == fy+short_time){
      tmp = generate_trends(indx,
                                     Min_year = y-short_time,
                                     Max_year = y)
    }else{
      tmp2 = generate_trends(indx,
                                      Min_year = y-short_time,
                                      Max_year = y)
      tmp = rbind(tmp,tmp2)
    }
    
    
  }#y
  tmp$model = "gamye"
  tmp$species = species
  
  
  load(paste0(sp_dir,model,"/parameter_model_run.RData"))  
  jags_mod_full = jags_mod_param
  indx = generate_indices(jags_mod = jags_mod_full,
                                   jags_data = jags_data,
                                   #quantiles = qs,
                                   regions = c("continental"),
                                   max_backcast = NULL,
                                   alternate_n = "n3")
  
  
  
  for(y in (fy+short_time):YYYY){
    
    tmp2 = generate_trends(indx,
                                    Min_year = y-short_time,
                                    Max_year = y)
    tmp2$model = "gamye_alt"
    tmp2$species = species
    
    tmp = rbind(tmp,tmp2)
    
    
  }#y
  
  
  for(model in models[-1]){
    load(paste0(sp_dir,model,"/jags_data.RData"))  
    
    load(paste0(sp_dir,model,"/jags_mod_full.RData")) 
    indx = generate_indices(jags_mod = jags_mod_full,
                                     jags_data = jags_data,
                                     #quantiles = qs,
                                     regions = c("continental"),
                                     max_backcast = NULL,
                                     alternate_n = "n")
    
    for(y in (fy+short_time):YYYY){
      
      tmp2 = generate_trends(indx,
                                      Min_year = y-short_time,
                                      Max_year = y)
      tmp2$model = model
      tmp2$species = species
      
      tmp = rbind(tmp,tmp2)
      
      
    }#y
    
  }
  
  if(species == demo_sp[1]){
    tmpout = tmp
  }else{
    tmpout = rbind(tmpout,tmp)
  }
  
  
}#species

modnames = data.frame(model = c("gamye",
                                "gam",
                                "firstdiff",
                                "slope",
                                "gamye_alt"),
                      version = c("GAMYE - Including Year Effects",
                                  "GAM",
                                  "DIFFERENCE",
                                  "SLOPE",
                                  "GAMYE - Smooth Only"))

modnames$version = factor(modnames$version,ordered = T,levels = modnames$version[c(5,2,3,4,1)])
trends = merge(tmpout,modnames,by = "model")

trends = trends[order(trends$species,trends$version,trends$End_year),]


mean_abs_dif <- function(x){
  mean(abs(diff(x)))
}
sd_abs_dif <- function(x){
  sd(abs(diff(x)))
}


trend_difs <- trends %>% filter(species %in% c("Carolina Wren","Pine Siskin")) %>% group_by(version,species) %>% summarise(m_dif = mean_abs_dif(Trend),
                                                                                                                           sd_dif = sd_abs_dif(Trend))
 # trend_difs <- trends %>% group_by(species,version) %>% summarise(m_dif = mean_abs_dif(Trend),
#                                                                                                                            sd_dif = sd_abs_dif(Trend))
# trends2 <- trends %>% filter(species %in% c("Carolina Wren","Pine Siskin"))
# tt = ggplot(data = trends2,aes(x = End_year,y = Trend,colour = species))+
#   geom_point()+geom_line()+facet_wrap(facets = ~version)
# print(tt)

#trend_difs <- arrange(trend_difs,species,m_dif)

coltrends = c(model_pallete,viridis::viridis(7)[6] )
names(coltrends) <- modnames$version
coltrends <- coltrends[c(2,5,3,4,1)]

lbls = filter(trend_difs,species == "Carolina Wren")
lbls$xx = 0.65+(as.integer(lbls$version)-1)*0.18
# lbls[which(lbls$version == "GAMYE − Including Year Effects"),"xx"] = lbls[which(lbls$version == "GAMYE − Including Year Effects"),"xx"]+0.2

trenddif = ggplot(data = trend_difs,aes(x = species,y = m_dif,colour = version,fill = version,group = version))+
  geom_bar(position = "dodge",stat = "identity")+
  scale_colour_manual(name = "",
                      values = coltrends, aesthetics = c("colour","fill"))+
  scale_y_continuous(limits = c(0,7)) +
  xlab("")+
  ylab("Mean absolute change in annual 10-year trends")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90,hjust = 0.5),
        axis.title.x = element_text(size = 9))+
  coord_flip()+
  geom_text_repel(data = lbls,aes(x = xx,y = m_dif, label = version),nudge_y = 3,size = 3)



pdf(file = paste0("Figures/Fig 7.pdf"),
    width = 3.5,
    height = 2.5)
print(trenddif)
dev.off()








# Figure 8 ----------------------------------------------------------------



species = "Wood Thrush"
model = "gamye"
sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,model,"/jags_data.RData"))  

load(paste0(sp_dir,model,"/jags_mod_full.RData")) 

load(paste0(sp_dir,model,"/parameter_model_run.RData"))  


inds_gamye <- generate_indices(jags_mod = jags_mod_full,
                                        jags_data = jags_data,
                                        regions = c("national"),
                                        max_backcast = NULL,
                                        alternate_n = "n") 

inds_gamnoye <- generate_indices(jags_mod = jags_mod_param,
                                          jags_data = jags_data,
                                          regions = c("national"),
                                          max_backcast = NULL,
                                          alternate_n = "n3") 


load(paste0(sp_dir,"slope/jags_mod_full.RData"))

inds_slope <- generate_indices(jags_mod = jags_mod_full,
                                        jags_data = jags_data,
                                        regions = c("national"),
                                        max_backcast = NULL,
                                        alternate_n = "n") 



load(paste0(sp_dir,"firstdiff/jags_mod_full.RData"))

inds_firstdiff <- generate_indices(jags_mod = jags_mod_full,
                               jags_data = jags_data,
                               regions = c("national"),
                               max_backcast = NULL,
                               alternate_n = "n") 




fy = 1990
short_time = 10
YYYY = max(jags_data$r_year)
rollTrend = "Trend"


for(ly2 in c((fy+short_time):YYYY)){
  trst_gamye = generate_trends(indices = inds_gamye,
                                        Min_year = ly2-short_time,
                                        Max_year = ly2,
                                        #quantiles = qs,
                                        slope = F,
                                        prob_decrease = c(0,25,30,50),
                                        prob_increase = c(0,33,100))
  trst_gamye$decomp <- "GAMYE - Including Year Effects"
  
  trst_gamnoye = generate_trends(indices = inds_gamnoye,
                                          Min_year = ly2-short_time,
                                          Max_year = ly2,
                                          #quantiles = qs,
                                          slope = F,
                                          prob_decrease = c(0,25,30,50),
                                          prob_increase = c(0,33,100))
  trst_gamnoye$decomp <- "GAMYE - Smooth only"
  
  
  trst_sl = generate_trends(indices = inds_slope,
                                     Min_year = ly2-short_time,
                                     Max_year = ly2,
                                     #quantiles = qs,
                                     slope = F,
                                     prob_decrease = c(0,25,30,50),
                                     prob_increase = c(0,33,100))
  trst_sl$decomp <- "SLOPE"
  
  trst_firstdiff = generate_trends(indices = inds_firstdiff,
                                     Min_year = ly2-short_time,
                                     Max_year = ly2,
                                     #quantiles = qs,
                                     slope = F,
                                     prob_decrease = c(0,25,30,50),
                                     prob_increase = c(0,33,100))
  trst_firstdiff$decomp <- "DIFFERENCE"
  
  
  
  if(ly2 == fy+short_time){
    tcos = rbind(trst_gamye,trst_gamnoye,trst_sl,trst_firstdiff)
    
  }else{
    tcos = rbind(tcos,trst_gamye,trst_gamnoye,trst_sl,trst_firstdiff)
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


colye2 = c(colye,model_pallete[c("slope","firstdiff")] )
names(colye2) <- unique(tcos$decomp)



rg = "Canada"
tmp = tcos[which(tcos$Region_alt == rg & tcos$End_year > 2000),]


st_exc <- ""

tmpend = tmp[which( (tmp$End_year == 2011 & tmp$decomp == "DIFFERENCE") |
                      (tmp$End_year == 2011 & tmp$decomp == "SLOPE") |
                      (tmp$End_year == 2014 & tmp$decomp == "GAMYE - Smooth only") |
                      (tmp$End_year == 2005 & tmp$decomp == "GAMYE - Including Year Effects")  ),]
tmpend$lably = tmpend$rolt
sm.v = which(tmpend$decomp == "GAMYE - Smooth only")

#tmpend[sm.v,"lably"] <- tmpend[sm.v,"roltuci"]
sl.v = which(tmpend$decomp == "SLOPE")

tmpend[sl.v,"End_year"] <- tmpend[sl.v,"End_year"] +0.1
sy.v <- which(tmpend$decomp == "GAMYE - Including Year Effects")
tmpend[sy.v,"End_year"] <- tmpend[sy.v,"End_year"] - 0.15


sd.v = which(tmpend$decomp == "DIFFERENCE")
tmpend[sd.v,"End_year"] <- tmpend[sd.v,"End_year"]-0.2


threshplot = data.frame(roltlci = c(rep(thresh30,2),rep(thresh50,2)),
                        roltuci = c(rep(thresh50,2),rep(-15,2)),
                        year = rep(c(1999,2019),2),
                        cat = rep(c("Threatened","Endangered"),each = 2))
lylim = min(tmp$roltlci)-0.75
uylim = max(tmp$roltuci)+0.75


cpt = ggplot(data = tmp,aes(x = End_year,y = rolt,group = decomp,colour = decomp))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(paste0("Year of Assessment (End of ",short_time,"-year trend)"))+
  ylab(paste0(short_time,"-year trends"))+
  geom_hline(yintercept = thresh30,colour = c_orng,size = 1)+
  geom_hline(yintercept = thresh50,colour = c_red,size = 1)+
  coord_cartesian(ylim = c(lylim,uylim),xlim = c(2000,2018))+
  geom_ribbon(data = threshplot[1:2,],aes(x = year,ymin = roltlci,ymax = roltuci),fill = c_orng ,alpha = 0.1,inherit.aes = F)+
  geom_ribbon(data = threshplot[3:4,],aes(x = year,ymin = roltlci,ymax = roltuci),fill = c_red ,alpha = 0.1,inherit.aes = F)+
  geom_hline(yintercept = 0,colour = grey(0.7))+
  annotate(geom = "text",x = 2002, label = "Threatened",y = threshs[1,"thresh"]-0.275,colour = c_orng,size = 3,alpha = 1)+
  annotate(geom = "text",x = 2002, label = "Endangered", y = threshs[2,"thresh"]-0.28,colour = c_red,size = 3,alpha = 1)+
  geom_linerange(aes(x = End_year,ymin = roltlci,ymax = roltuci),alpha = 0.15,size = 0.2,position = position_dodge(width = 0.4))+
  geom_point(aes(x = End_year,y = rolt),size = 0.8,position = position_dodge(width = 0.4))+
  geom_line(aes(x = End_year,y = rolt),size = 0.2,position = position_dodge(width = 0.4),alpha = 0.4)+
  geom_text_repel(data = tmpend[sy.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_x = -0.5,nudge_y = +2,size = 4)+
  geom_text_repel(data = tmpend[sl.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_x = 2,nudge_y = -1,size = 4)+
  geom_text_repel(data = tmpend[sm.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_y = 3,nudge_x = +1.5,size = 4)+
  geom_text_repel(data = tmpend[sd.v,],inherit.aes = F,aes(x = End_year,y = lably, label = decomp,colour = decomp),nudge_x = -4,nudge_y = -1.7,size = 4)+
  scale_colour_manual(values = colye2, aesthetics = c("colour","fill"))



## update the theme?
pdf(paste0(paste0("Figures/Fig 8.pdf")),
    width = 3.5,
    height = 3)
print(cpt)
dev.off()







# Figure 9 ----------------------------------------------------------------

species = "Barn Swallow"

sp_dir = paste0("output/",species,"/")

load(paste0(sp_dir,"saved objects.RData"))

load(paste0(sp_dir,"saved objects4.RData"))
loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))



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
  coord_cartesian(ylim = c(-0.02,0.06))+
  scale_x_continuous(breaks = c(seq(1970,2010,by = 10),2018))+
  scale_colour_viridis_d(option = "C",end = 0.5)+
  geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
  geom_text_repel(data = lbl[1,],aes(x = Year,y = lci,label = Contrast_full_name),nudge_y = -0.005)+
  geom_text_repel(data = lbl[2,],aes(x = Year,y = uci,label = Contrast_full_name),nudge_y = +0.01,nudge_x = 13)+
  annotate(geom = "text",x = 2001,y = 0.03,label = "Positive favours GAMYE",size = 3.5,colour = grey(0.6))+
  annotate(geom = "text",x = 2001,y = -0.02,label = "Negative favours alternate",size = 3.5,colour = grey(0.6))  


pdf(paste0("Figures/Fig 9.pdf"),
    width = 3.5,
    height = 4)
print(an_contr)
dev.off()


# species versions in supplment -------------------------------------------


svplots = list()
length(svplots) = length(demo_sp)
names(svplots) = demo_sp

pdf(paste0("Figures/supplement/Fig 9 all species.pdf"),
    width = 5,
    height = 4)
for(species in demo_sp){
  
  sp_dir = paste0("output/",species,"/")
  
  load(paste0(sp_dir,"saved objects.RData"))
  
  load(paste0(sp_dir,"saved objects4.RData"))
  loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))
  
  
  
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
    labs(title = species)+
    ylab("Mean difference in point-wise log-probability")+
    xlab("")+
    scale_x_continuous(breaks = c(seq(1970,2010,by = 10),2018))+
    scale_colour_viridis_d(option = "C",end = 0.5)+
    geom_hline(yintercept = 0,colour = grey(0.2),alpha = 0.2)+
    geom_text_repel(data = lbl[1,],aes(x = Year,y = lci,label = Contrast_full_name),nudge_y = -0.005)+
    geom_text_repel(data = lbl[2,],aes(x = Year,y = uci,label = Contrast_full_name),nudge_y = +0.005,nudge_x = 8)+
    annotate(geom = "text",x = 2001,y = 0.04,label = "Positive favours GAMYE",size = 3.5,colour = grey(0.6))+
    annotate(geom = "text",x = 2001,y = -0.012,label = "Negative favours alternate",size = 3.5,colour = grey(0.6))  
  

  print(an_contr)
  svplots[[species]] <- an_contr
  
}
dev.off()
save(list = "svplots",file = paste0("Figures/supplement/Fig 9 all species.RData"))



# Figure 10 ----------------------------------------------------------------


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

      if(i == "B"){
      pdf(file = paste0("Figures/Fig 10",i,".pdf"),
          width = 3.5,
          height = 3.5)
      }else{
        pdf(file = paste0("Figures/Fig 10.pdf"),
            width = 3.5,
            height = 3.5)
      }
    
    
    print(mp.plot)
    

  dev.off()

}





# elpd Distribution demo for supplement --------------------------------------------------




species = "Barn Swallow"
  
  
  sp_dir = paste0("output/",species,"/")
  
  loo.point = read.csv(paste0(sp_dir,"wide form lppd.csv"))
  
  
  loo.point_stack <- pivot_longer(data = loo.point,cols = gamye_gam:firstdiff_slope,
                                  names_to = "model_comparison",
                                  values_to = "delta_elpd")
  
  contrast_full_names_df <- data.frame(model_comparison = names(contrast_full_names),
                                       comp = contrast_full_names)
  
  
  loo.point_stack <- left_join(loo.point_stack,contrast_full_names_df)
  
  qqout = ggplot(data = loo.point_stack,aes(sample = delta_elpd,group = comp))+
    geom_qq()+
    geom_qq_line()+
    facet_wrap(~comp,nrow = 2, ncol = 3)+    
   theme_minimal()+
    theme(line = element_line(size = 0.4), rect = element_rect(size = 0.1),
          axis.line = element_blank())
    
print(qqout)

  
  
  





# supplmental figure compile ----------------------------------------------



load(file = "c:/GAM_Paper_Script/figures/supplement/Fig 1 all species.RData")
svplots1 = svplots
rm(svplots)
load(file = "c:/GAM_Paper_Script/figures/supplement/Fig 2 all species.RData")
svplots2 = svplots
rm(svplots)
load(file = "c:/GAM_Paper_Script/figures/supplement/Fig 4 all species.RData")
svplots4 = svplots
rm(svplots)
load(file = "c:/GAM_Paper_Script/figures/supplement/Fig 5 all models.RData")

load(file = "c:/GAM_Paper_Script/figures/supplement/Fig 6 all species.RData")
svplots6 = svplots
rm(svplots)
load(file = "c:/GAM_Paper_Script/figures/supplement/Fig 9 all species.RData")
svplots9 = svplots
rm(svplots)


save(list = c("svplots1",
              "svplots2",
              "overall.comparison",
              "svplots4",
              "svplots6",
              "svplots9",
              "qqout"),
     file = "figures/supplement/all_suppl_figures.RData")








