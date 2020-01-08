library(bbsBayes)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(tidyverse)


models = c("gamye","gam","firstdiff","slope")
heavy_tailed = TRUE #all models use the t-distribution to model extra-Poisson variance

species_to_run = c("Wood Thrush", "American Kestrel","Barn Swallow","Chestnut-collared Longspur","Ruby-throated Hummingbird")


model = models[1]



for(species in species_to_run){
  
  sp_dir = paste0("output/",species,"/")

  
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
}

betaplot = ggplot(data = conti,aes(x = year,y = med))+
  theme_minimal()+
  labs(title = paste0(species," GAM components including hyperparameter"))+
  ylab("Population change based on GAM smooth (linear scale)")+
  theme(legend.position = "none")+
  geom_line(data = strati,aes(x = year, y = med,group = strat),alpha = 0.1)+
  geom_line(colour = grey(0.2),size = 1.4)+
  coord_cartesian(ylim = c(0,max(conti$med)*2))

pdf(file = paste0(sp_dir,species," GAM components.pdf"),
    width = 7,
    height = 5)
 print(betaplot)
 dev.off()
 
} 
  #load(paste0(sp_dir,model,"/jags_mod_full.RData"))







