library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)
library(RColorBrewer)
#devtools::install_github("zeehio/facetscales")
library(facetscales)

models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

species_to_run = c("Horned Lark","Cooper's Hawk","Wood Thrush", "American Kestrel","Barn Swallow","Chestnut-collared Longspur","Ruby-throated Hummingbird")


model = models[1]

pooling = list()
length(pooling) = length(species_to_run)
names(pooling) <- species_to_run

for(species in species_to_run){
  
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
  
#plot(conti$year,conti$med)
 

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
strati$prec = 1/((strati$uci-strati$lci)/(1.96*2))^2
scounts = table(raw.dat$strat)
for(s in 1:nstrata){
  strati[which(strati$strat == s),"alpha"] = (1/(sqrt(scounts[s])/sqrt(min(scounts))))*0.5
  strati[which(strati$strat == s),"ncounts"] = (scounts[s])
  strati[which(strati$strat == s),"strat_name"] <- unique(jags_data$strat_name[which(jags_data$strat == s)])
  
}

betaplot = ggplot(data = conti,aes(x = year,y = med))+
  theme_classic()+
  labs(title = "")+#paste0(species," GAM components including hyperparameter"))+
  ylab("Population change based on GAM smooth (linear scale)")+
  theme(legend.position = "none")+
  geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.1)+
  coord_cartesian(ylim = c(0,max(conti$med)*2))+
  geom_line(colour = grey(0.2),size = 1.4)

pdf(file = paste0(sp_dir,species," ",model," components.pdf"),
    width = 7,
    height = 5)
 print(betaplot)
 dev.off()
 
 betaplot = ggplot(data = conti,aes(x = year,y = med))+
   theme_minimal()+
   labs(title = paste0(species," GAM components including hyperparameter"))+
   ylab("Population change based on GAM smooth (linear scale)")+
   theme(legend.position = "none")+
   geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.1)+
   #coord_cartesian(ylim = c(0,max(conti$med)*2))+
   scale_y_log10()+
   geom_line(colour = grey(0.2),size = 1.4)
 
 pdf(file = paste0(sp_dir,species," ",model," logy components.pdf"),
     width = 7,
     height = 5)
 print(betaplot)
 dev.off()
 
 
 
 strati$species = species
 conti$species = species
 
 
 

# pooling factor calculation ----------------------------------------------

 sdbeta.mat = 1/sqrt(jags_mod_param$sims.list[["taubeta"]])
 BXmat = jags_mod_param$sims.list[["B.X"]]
 beta.x.mat = jags_mod_param$sims.list[["beta.X"]]
 
 poola = NA
 pool = matrix(NA,nrow = raw.dat$nstrata,ncol = raw.dat$nknots)
 for (k in 1:raw.dat$nknots){
   bdif = beta.x.mat[,,k]-BXmat[,k]
   pm_var_strat_bdif = mean (apply (bdif, 1, var))
   
   var_pmean_bdifa = var(apply(bdif,2,mean))
   
   #poola[k] <- min(var_pmean_bdifa/pm_var_strat_bdif  ,1)
   
   for(s in 1:raw.dat$nstrata){
     pvar_mean_bdif = var(bdif[,s])

 pool[s,k] <- min(pvar_mean_bdif/pm_var_strat_bdif  ,1)


   }
   
   
 }
 
 pool_by_strat = rowMeans(pool[,-c(4:10)]) #summarizing the pooling for the outer knots
 
 
 for(s in 1:nstrata){
   strati[which(strati$strat == s),"pool"] = (pool_by_strat[s])
   
 }
 
 
 strati$pool_q = NA
 strati[which(strati$pool < quantile(strati$pool,0.15)),"pool_q"] <- c("More counts")
 strati[which(strati$pool > quantile(strati$pool,0.85)),"pool_q"] <- c("Fewer counts")

 strati$data_q = NA
 strati[which(strati$ncounts < quantile(strati$ncounts,0.15)),"data_q"] <- c("Fewer counts")
 strati[which(strati$ncounts > quantile(strati$ncounts,0.85)),"data_q"] <- c("More counts")
 
 # strati$qual = paste(strati$data_q,strati$pool_q,sep = "_")
 conti2 = conti
 conti2$data_q <- "Fewer counts"
 conti3 = conti
 conti3$data_q <- "More counts"
 contiplot = rbind(conti2,conti3)
 betaplot = ggplot(data = strati,aes(x = year,y = med, group = data_q))+
   theme_minimal()+
   labs(title = paste0(species," GAM components including hyperparameter"))+
   ylab("Population change based on GAM smooth (linear scale)")+
   theme(legend.position = "none")+
   geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.1)+
   #coord_cartesian(ylim = c(0.1,max(conti$med)*2))+
   scale_y_log10()+
   geom_line(data = contiplot,aes(x = year,y = med),colour = grey(0.2),size = 1.4)+
   facet_col(facets = ~data_q,drop = T)
 
 pdf(file = paste0(sp_dir,species," ",model,"split components.pdf"),
     width = 5,
     height = 8)
 print(betaplot)
 dev.off()
 
 
 
 
 pooling[[species]] <- pool
 
 if(species == species_to_run[1]){
   stratiall = strati
   contiall = conti
 }else{
   stratiall = rbind(stratiall,strati)
   contiall = rbind(contiall,conti)
   
 }
 
} 
  #load(paste0(sp_dir,model,"/jags_mod_full.RData"))




source("colourblind safe qualitative pallete.r")
species_pallete <- safe.pallet[[length(species_to_run)]] 
names(species_pallete) <- species_to_run


allysc = list()
length(allysc) = length(species_to_run)
names(allysc) = species_to_run

for(s in species_to_run){
  tty = max(contiall[which(contiall$species == s),"med"])*1.75
  #allysc[[s]] <- coord_cartesian(ylim = c(0,tty))
  allysc[[s]] <- scale_y_continuous(limits = c(0,tty))
  
}


betaplot = ggplot(data = contiall,aes(x = year,y = med))+
  theme_classic()+
  theme(legend.position = "none",
        strip.text.y = element_text(size = 7),
        strip.background.y = element_rect(linetype = 0))+
  labs(title = "",
       x = "")+
  ylab("Population change based on GAM smooth (linear scale)")+
  theme(legend.position = "none")+
  geom_line(data = stratiall,aes(x = year, y = med,group = strat),alpha = 0.05)+
  geom_line(data = stratiall,aes(x = year, y = med,group = strat,colour = species),alpha = 0.13)+
  scale_colour_manual(values = species_pallete, aesthetics = c("colour"))+
  geom_line(aes(colour = species),size = 1.4)+
  geom_line(size = 1.4,alpha = 0.05)+
  facet_grid_sc(rows = vars(species),scales = list(y = allysc))#,switch = "y")
 
pdf(file = paste0("All species ",model,"-based change Figure 2.pdf"),
    height = 9,
    width = 4)
print(betaplot)
dev.off()













# rolling trends to show annual variability in trend estimates ------------


source("colourblind safe qualitative pallete.r")
model_pallete <- safe.pallet[[length(models)]] 
model_pallete <- model_palletec(2,1,3,4)
names(model_pallete) <- models


c_orng = brewer.pal(9,"Set1")[5]
c_red = brewer.pal(9,"Set1")[1]
c_blue = brewer.pal(9,"Set1")[2]
c_purp = brewer.pal(9,"Set1")[4]
c_green = brewer.pal(9,"Set1")[3]

######################## Figure 5
#add a script to plot the continental full and smooth-only trajectories for BARS, with all the regular additional info
########################
for(species in species_to_run){
  
  
  for (model in models){
    
    sp_dir = paste0("output/",species,"/")
    
    load(paste0(sp_dir,model,"/jags_data.RData"))  
    if(model == models[1]){
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
      
    }else{
    load(paste0(sp_dir,model,"/jags_mod_full.RData"))  
    }
    

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
    
    cont_over = ggplot(data = indcont,aes(x = Year,y = Index,group = version))+
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
    
    
    
    
fy2 = min(1995,fy)

if(model == models[1]){
  
indscos = generate_regional_indices(jags_mod = jags_mod_full,
                                    jags_data = jags_data,
                                    #quantiles = qs,
                                    regions = c("continental","national"),
                                    startyear = fy2,
                                    max_backcast = NULL,
                                    alternate_n = "n3")


}else{
  indscos = generate_regional_indices(jags_mod = jags_mod_full,
                                      jags_data = jags_data,
                                      #quantiles = qs,
                                      regions = c("continental","national"),
                                      startyear = fy2,
                                      max_backcast = NULL) 
}
for(ly2 in c((fy2+short_time):YYYY)){
  trst = generate_regional_trends(indices = indscos,
                                  Min_year = ly2-short_time,
                                  Max_year = ly2,
                                  #quantiles = qs,
                                  slope = F,
                                  prob_decrease = c(0,25,30,50),
                                  prob_increase = c(0,33,100))
  if(ly2 == fy2+short_time){
    tcos = trst
  }else{
    tcos = rbind(tcos,trst)
  }
  
}

tcos$rolt = tcos[,rollTrend]
tcos$roltlci = tcos[,paste0(rollTrend,"_Q0.025")]
tcos$roltlci2 = tcos[,paste0(rollTrend,"_Q0.25")]
tcos$roltuci = tcos[,paste0(rollTrend,"_Q0.975")]
tcos$roltuci2 = tcos[,paste0(rollTrend,"_Q0.75")]
tcos$model = model
if(model == models[1]){
  tcosplot = tcos
}else{
  tcosplot = rbind(tcosplot,tcos)
}

  }#models

write.csv(tcosplot,paste0(sp_dir,"Rolling Trends.csv"),row.names = F)

tcosplot$species = species
if(species == species_to_run[1]){
  roll_trends = tcosplot
}else{
  roll_trends = rbind(roll_trends,tcosplot)
}

thresh30 = (0.7^(1/short_time)-1)*100
thresh50 = (0.5^(1/short_time)-1)*100

threshs = data.frame(thresh = c(thresh30,thresh50),
                     p_thresh = c(paste("-30% over",short_time,"years"),
                                  paste("-50% over",short_time,"years")),
                     Year = rep(min(tcos$End_year),2))

pdf(paste0(paste0(sp_dir,"Rolling_Trends.pdf")),
    width = 8.5,
    height = 6)
for(rg in unique(tcosplot$Region_alt)){
  
  tmp = tcosplot[which(tcosplot$Region_alt == rg),]
  
  st_exc = unique(tmp$Strata_excluded)
  if(st_exc != ""){
    if(nchar(st_exc) > 20){
      stx2 = unlist(strsplit(st_exc,split = " ; "))
      st_exc <- paste("Excluding",length(stx2),"strata")
    }else{
      st_exc <- paste("Excluding",st_exc)
    }}
  
  
  tmpend4 = tmp[nrow(tmp)-4,]
  tmpend4$lab50 = "50% CI"
  tmpend4$lab95 = "95% CI"
  
  tmpend = tmp[nrow(tmp),]
  
  pth_30_labs = paste0(signif(100*tmpend[,"prob_decrease_30_percent"],2),"% probability of 30% decrease") 
  pth_50_labs = paste0(signif(100*tmpend[,"prob_decrease_50_percent"],2),"% probability of 50% decrease") 
  tmpend$pdec = paste(signif(tmpend[,"Percent_Change"],2),"% Change over",short_time,"years") 
  
  
  cpt = ggplot(data = tmp,aes(x = End_year,y = rolt,group = model,colour = model))+
    theme_minimal()+
    theme(legend.position = "none")+
    labs(title = paste(species,"rolling",short_time,"year trends",rg,st_exc),
         subtitle = paste("Based on",rollTrend,"in",YYYY,":",pth_30_labs,"and",pth_50_labs))+
    xlab(paste("Ending year of",short_time,"trend"))+
    ylab(paste(short_time,"year trends"))+
    geom_hline(yintercept = thresh30,colour = c_orng)+
    geom_hline(yintercept = thresh50,colour = c_red)+
    geom_hline(yintercept = 0,colour = grey(0.5))+
    geom_label_repel(data = threshs,aes(x = Year,y = thresh,label = p_thresh),position = "nudge")+
    geom_linerange(aes(x = End_year,ymin = roltlci,ymax = roltuci),alpha = 0.4,size = 0.9,position = position_dodge(width = 0.3))+
    geom_point(aes(x = End_year,y = rolt),size = 1,position = position_dodge(width = 0.3))+
    scale_colour_manual(values = model_pallete, aesthetics = c("colour","fill"))

  
  
  
  ## update the theme?
  print(cpt)
  
}
dev.off()



  }#species



roll_trends_sort <- arrange(roll_trends,species,model,Region_alt,End_year)

write.csv(roll_trends_sort,"Rolling trends for all species and models.csv",row.names = F)





# plotting rolling trend results ------------------------------------------


roll_trends_sort = read.csv("Rolling trends for all species and models.csv",stringsAsFactors = F)



my_acf <- function(x,lg = 1){
  ac = acf(x,lag.max = lg)
  ac = ac$acf[,,1][lg+1]
  return(ac)
}


my_diff <- function(x,lg = 1){
  ac = diff(x,lag = lg)
  mac = mean(abs(ac))
  return(mac)
}

# 
# my_diff <- function(x){
#   ac = diff(x)
#   mac = mean(abs(ac))
#   return(mac)
# }
# 
# my_diff2 <- function(x){
#   ac = diff(x,lag = 2)
#   mac = mean(abs(ac))
#   return(mac)
# }



acf_by_sp <- roll_trends_sort %>% group_by(species,model,Region_alt) %>% 
  summarise(.,acf = my_acf(Trend),acf10 = my_acf(Trend,10))

dif_by_sp <- roll_trends_sort %>% group_by(species,model,Region_alt) %>% 
  summarise(.,dif = my_diff(Trend), dif10 = my_diff(Trend,10))


  acfp = ggplot(data = acf_by_sp[which(acf_by_sp$Region_alt == "Continental"),],
                aes(x = model,y = acf,group = species,colour = species))+
    geom_point()+
    geom_line()
  x11()
  print(acfp)
  
  
  acfp10 = ggplot(data = acf_by_sp[which(acf_by_sp$Region_alt == "Continental"),],
                aes(x = model,y = acf10,group = species,colour = species))+
    geom_point()+
    geom_line()
  x11()
  print(acfp10)

  dffp = ggplot(data = dif_by_sp[which(dif_by_sp$Region_alt == "Continental"),],
                aes(x = model,y = dif,group = species,colour = species))+
    geom_point()+
    geom_line()
  x11()
  print(dffp)

  
  dffp10 = ggplot(data = dif_by_sp[which(dif_by_sp$Region_alt == "Continental"),],
                aes(x = model,y = dif10,group = species,colour = species))+
    geom_point()+
    geom_line()
  x11()
  print(dffp10)
  







