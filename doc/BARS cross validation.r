pkgs = c("sp","maptools","rgdal","maptools","RColorBrewer","GISTools","ggplot2",
"dplyr","tidyr","scales","geofacet","stringr","tools",
"ggrepel","e1071")
inst.pck = installed.packages()
  inst.pck2 = names(inst.pck[,2])

  if(length(pkgs[-which(pkgs %in% inst.pck2)])>0){
install.packages(pkgs[-which(pkgs %in% inst.pck2)])
  }
  
for(i in 1:length(pkgs)){
library(pkgs[i],character.only = T)
}

#for skewness function
source("c:/functions/transparency function.r")


  
  
  
  
  # sp = "BARS"
  # setwd(paste0("C:/BBSCV/",sp))
  # load("full GAM Jan 27.Rdata")
  # 
  ###plot the stratum-level smooths along with the hyper-smooth
  ## use bbsBayes to run the GAM
  ## extract n and add a variable "N" that tracks the range-wide smooth
  ## consider centering
  




for(sp in c("EWPW","CONI","AMKE","BARS")){
#setwd("M:/My Documents/Coop Staff/BARS")
setwd(paste0("C:/BBSCV/",sp))
source("c:/functions/multiplot.r")
year_list <- NULL
index <- NULL
index_neg <- NULL
index_pos <- NULL
model_used <- NULL

aw <- read.csv(paste0("strat.original",sp,".txt"))
if("stratum" %in% names(aw)){
  names(aw)[which(names(aw) == "stratum")] = "strat"
}
strats = aw[,c("strat.name","strat","stratcode")]
strats = strats[order(strats$strat),]
for(i in 1:51){
  if(i == 1){
  stratsy = strats}else{
    stratsy = bind_rows(stratsy,strats)
  }
}
stratsy$year = rep(c(1:51),each = nrow(strats))
##### compiling the indices for all models

### rerun the loop below to use the median instead of hte mean
### rerun the AMKE standard model using jags 4.2

models <- c("Standard", "FD", "GAM","GAMYE")
owd <- getwd()


setwd(paste0("F:/Summer 2018/Model Analysis/input/",sp,"/trend"))

for (model in models)
{
  load(paste(model, ".Rdata", sep = ""))

 n <- jagsModFull$sims.list$n
  n_samples <- dim(n)[1]
  n_strata <- dim(n)[2]
  n_weight <- n

  for (i in 1:n_samples)
  {
    for (j in 1:n_strata)
    {
      n_weight[i,j,] <- (n_weight[i,j,] * aw$Area[which(aw$strat == j)])/ sum(aw$Area)
    }
  }

  N = apply(n_weight, c(1,3),sum)

  n_mean <- apply(N, 2, median)
  n_25 <- apply(N, 2, quantile, probs = 0.025)
  n_975 <- apply(N, 2, quantile, probs = 0.975)

  year_list <- c(seq(1:51))
  indextc <- c((n_mean))
  index_negtc <- c((n_25))
  index_postc <- c((n_975))
  model_used <- c(rep(model, 51))


  index.cont.t <- data.frame(Year = year_list+1965,
                             Index = indextc,
                             Index_neg = index_negtc,
                             Index_pos = index_postc,
                             Model = model_used)

  index.cont.t$Model <- factor(index.cont.t$Model, levels = models)


  ################### end continental index calculation

  #################### beginning strata-level index calculation

  n_mean <- as.data.frame(jagsModFull$q50$n)
  n_25 <- as.data.frame(jagsModFull$q2.5$n)
  n_975 <- as.data.frame(jagsModFull$q97.5$n)



  indext = gather(n_mean[,paste0("V",1:51)],
                 key = "yearv",
                 value = "index",
                 factor_key = T)
  indext = cbind(stratsy,indext)

  lci = gather(n_25[,paste0("V",1:51)],
               key = "yearvl",
               value = "lci",
               factor_key = T)
  indext = cbind(indext,lci)
  uci = gather(n_975[,paste0("V",1:51)],
               key = "yearvu",
               value = "uci",
               factor_key = T)
  indext = cbind(indext,uci)

  if(any(indext$yearv != indext$yearvl | indext$yearvl != indext$yearvu)){
    print(paste(model,"years don't line up"))
    break}
  indext = indext[,-which(names(indext) %in% c("yearv","yearvl","yearvu"))]
  indext$model = model

  if(model == models[1]){
     index = indext
  index.cont = index.cont.t
  }else{
    index = rbind(index,indext)
    index.cont = rbind(index.cont,index.cont.t)
  }


}
index$rYear = index$year + 1965
index.cont$rYear = index.cont$Year

setwd(owd)

save(list = c("index",
              "index.cont",
              "models"),
     file = paste(sp,"all model indices.RData"))


}#temp sp loop end

  
load(paste(sp,"all model indices.RData"))


############# facet plot for strata
############# facet plot for strata
############# facet plot for strata
############# facet plot for strata
############# facet plot for strata
############# facet plot for strata

# load("full GAM Jan 24.RData")

stprovfacet = read.csv("c:/BBSCV/strata/BBB_StateProvCWS_facet_grid2.csv",stringsAsFactors = F)
#stprovfacet = stprovfacet[-which(stprovfacet$code == "HI"),]
### trajectory plots by prov state
stprovfacet[nrow(stprovfacet)+1,"row"] = 1
stprovfacet[nrow(stprovfacet),"col"] = 5
stprovfacet[nrow(stprovfacet),"code"] = "BCR7"
stprovfacet[nrow(stprovfacet),"name"] = "Bcr7" 

country = rep("USA",nrow(stprovfacet))
country[which(stprovfacet$code %in% c("BC",
                                          "YT",
                                          "NT",
                                          "NU",
                                          "BCR7",
                                          "NSPE",
                                          "PE",
                                          "NS",
                                          "ON",
                                          "AB",
                                          "SK",
                                          "MB",
                                          "QC",
                                          "NB",
                                          "NL"))] = "CAN"
strcol = rep("darkslategray1",length(country))
strcol[which(stprovfacet$country == "CAN")] = "lightpink"

# for(m in models){
#   indt = index[which(index$model == m),]
#   # indt$prst = gsub(str_extract(indt$stratcode,"-(.*)-"),pattern = "-",replacement = "",fixed = T)
#   # indt$prst = as.character(factor(indt$prst,levels = unique(stprovfacet$code)))
#   indt$prst = gsub(str_extract(indt$strat.name,".*-"),pattern = "-",replacement = "",fixed = T)
#   indt$prst = toTitleCase(tolower(indt$prst))
#   indt$prst = as.character(factor(indt$prst,levels = unique(stprovfacet$name)))
#   indt$bcr = gsub(str_extract(indt$strat.name,"-.*"),pattern = "-BCR",replacement = "",fixed = T)
#   labsbcr = indt[which(indt$rYear == max(indt$rYear,na.rm = T)),]
#    ptraj <- ggplot(data = indt) + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank(),
#           axis.text = element_text(colour = grey(0.2)),
#           strip.background = element_rect(fill = grey(0.97)),#strcol #, colour = grey(0.9), size = NULL, linetype = NULL, color = NULL, inherit.blank = FALSE
#           axis.line = element_line(colour = "black"),
#           legend.position = "none") +
#     labs(title = paste(sp,"Population trajectories by strata (separate BCRS) within Provinces and States", m, "model"), x = "", y = "Mean predicted count of birds on an average BBS route") +
# 
#     #geom_pointrange(data = indt, aes(x = rYear, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
#     geom_line(data = indt, aes(x = rYear, y = index,group = stratcode),colour = grey(0.4)) +#, colour = stratcode)) +
#     geom_ribbon(data = indt, aes(x = rYear, ymin = lci, ymax = uci,group = stratcode),fill = grey(0.3),alpha = 0.12)+ #, fill = stratcode , fill = grey(0.8)
#      geom_text_repel(data = labsbcr, aes(x = rYear, y = index,label = bcr),
#                      size = 2)+
#      facet_geo(facets = ~ prst,grid = stprovfacet,scales = "free_y", label = "code")+
#   scale_x_continuous(limits = c(1966, 2017), oob = squish, breaks = c(1970,2017),minor_breaks = c(1980,1990,2000,2010))
# pdf(file = paste0(sp," geofacet ",m," trajectory plot.pdf"),
#     height = 8.5,width = 11)
#   print(ptraj)
#   dev.off()
#     # for(p in stprovfacet$code)){
#   #   
#   # }#p
#   
#   
# }#m
# 
# #}#tempoarry end sp loop

#grid_design(data = stprov)

# write.csv(table(strat.cent[,c("ycol","xcol")]),paste("centroid.test",yb,xb,".csv"))
# ############# facet plot for strata
############# facet plot for strata
############# facet plot for strata
############# facet plot for strata
#########
cv = read.csv("bugs_loocv.csv")
for(m in c("Standard_loocv",
           "GAM_loocv",
           "GAMYE_loocv",
           "FD_loocv")){
  remt = which(abs(cv[,m]) > quantile(abs(cv[,m]),probs = 0.999))
  reminft = which(is.infinite(cv[,m]))
  if(m == "Standard_loocv"){
    rem = remt
    reminf = reminft
  }else{
    rem = unique(c(rem,remt))
    reminf = unique(c(reminf,reminft))
    
  }
}
### if comparing the 99.9% center of the distributions
#cv = cv[-rem,]
if(length(reminf) > 0){
  cv = cv[-reminf,]
}
  #cv = cv[-which(cv$rYear < 1980),]
plotstrat = F

pdf(paste0("All continental trajectories.pdf"),
    height = 11,
    width = 8.5)
# plot(x = (cv$count), y = cv[,cmpn],
#      main = cmpn)


tmp = index.cont#[which(index.cont$Model %in% c(m1,m2)),]
data_summary <- data.frame(Year = tmp$rYear,
                           Index = tmp$Index,
                           Index_neg = tmp$Index_neg,
                           Index_pos = tmp$Index_pos,
                           Model = tmp$Model)

data_summary$Model <- factor(data_summary$Model, levels = models,ordered = T)
by.yN = tapply(cv[,"count"],cv[,c("rYear")],mean,na.rm = T)
by.yNq1 = tapply(cv[,"count"],cv[,c("rYear")],quantile,probs = 0.25,na.rm = T)
by.yNq3 = tapply(cv[,"count"],cv[,c("rYear")],quantile,probs = 0.75,na.rm = T)


N = by.yN
Nl =  by.yNq1
Nu =by.yNq3

Ndf = data.frame(mean = as.numeric(N),
                 lci = as.numeric(Nl),
                 uci = as.numeric(Nu),
                 Year = 1966:2016)

# Ndf = loobyyear
# Ndf$Year = Ndf$rYear


###################################################
# Plot all 4 trajectories
###################################################

ptraj <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste("All Continental Population Trajectories"), x = "Year", y = "Index") + 
  geom_point(data = Ndf, aes(x = Year, y = mean),alpha = 0.3) +
  geom_line(data = data_summary, aes(x = Year, y = Index, colour = Model)) +
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  geom_ribbon(data = data_summary, aes(x = Year, ymin = Index_neg, ymax = Index_pos, fill = Model), alpha = 0.12)+
  coord_cartesian(ylim = c(0, max(Ndf$mean)),  xlim = c(1966, 2016))
print(ptraj)
dev.off()


#########################
#plotting the loocv scores
#########################


grpvar = "rYear"

for(grpvar in c("rYear","strat.name")){


for(m in models){
  cl = paste0(m,"_loocv")
luni = function(x){length(unique(x))}
byvar = tapply(-2*(cv[,cl]),cv[,grpvar],mean)
sdbyvar = tapply(-2*(cv[,cl]),cv[,grpvar],sd)
nrtsbyvar = tapply(cv[,"route"],cv[,grpvar],luni)
ncountsbyvar = tapply(cv[,"route"],cv[,grpvar],length)
nyrsbyvar = tapply(cv[,"year"],cv[,grpvar],luni)
skewbyvar = tapply(cv[,"count"],cv[,grpvar],skewness)

tmp = data.frame(var = names(byvar),
                        loocv = byvar,
                        loocv_se = sdbyvar/sqrt(ncountsbyvar),
                 loocv_lci = byvar-(sdbyvar/sqrt(ncountsbyvar)*1.96),
                 loocv_uci = byvar+(sdbyvar/sqrt(ncountsbyvar)*1.96),
                        nrts = nrtsbyvar,
                        ncounts = ncountsbyvar,
                        nyrs = nyrsbyvar,
                        skew = skewbyvar,
                      model = m)
tmp2 = cv[,c("count","strat","obser","year","firstyr","strat.name","route","rYear",cl)]
names(tmp2)[which(names(tmp2) == cl)] <- "loocv"
names(tmp2)[which(names(tmp2) == grpvar)] <- "var"
tmp2$loocv = -2*tmp2$loocv
tmp2$model = m
if(m == models[1]){
  comp.var = tmp
  stackcv = tmp2
}else{
  comp.var = rbind(comp.var,tmp)
  stackcv = rbind(stackcv,tmp2)
}

}#m

  comp.var$model = factor(comp.var$model,levels = models,ordered = T)
  
  stackcv$model = factor(stackcv$model,levels = models,ordered = T)
  
  if(grpvar == "strat.name"){
    comp.var = comp.var[rev(order(comp.var$ncounts)),]
    comp.var$var = factor(comp.var$var,ordered = T)
    stackcv$var = factor(stackcv$var,levels = levels(comp.var$var),ordered = T) 
    
  }
  
  
  pdf(paste(grpvar,"summary loocv scores.pdf"),
      width = 12,height = 5)
  
  varp = ggplot(data = comp.var,aes(x = var,y = loocv,colour = model))+
    geom_pointrange(aes(ymin = loocv_lci,ymax = loocv_uci),size = 0.15,position = position_dodge(width = 0.5))+
    scale_y_continuous(limits = c(0,max(comp.var$loocv_uci)))+
    #scale_fill_brewer(palette="Dark2")+
    scale_color_brewer(palette="Dark2")+
    labs(title = paste(sp,"loocv by",grpvar))
 
  
  print(varp) 
  dev.off()
  
  
  pdf(paste(grpvar,"violin loocv scores.pdf"),
            width = 12,height = 5)
  
  #stackcv$var = factor(stackcv$var,ordered = T)
  varp = ggplot(data = stackcv,aes(x = var,y = loocv,colour = model))+
    geom_violin(scale = "area",position = position_dodge(width = 0.5))+
    #scale_y_continuous(limits = c(0,max(comp.var$loocv_uci)))+
    scale_fill_brewer(palette="Dark2")+
    scale_color_brewer(palette="Dark2")+
    labs(title = paste(sp,"loocv by",grpvar))
  print(varp) 
 dev.off()
 
  
  pdf(paste(grpvar,"dotplot loocv scores.pdf"),
      width = 12,height = 5)
  
  varp = ggplot(data = stackcv,aes(x = var,y = loocv,colour = model))+
    geom_point(position = position_dodge(width = 0.5),size = 0.15)+
    #scale_y_continuous(limits = c(0,max(comp.var$loocv_uci)))+
    scale_fill_brewer(palette="Dark2")+
    scale_color_brewer(palette="Dark2")+
    labs(title = paste(sp,"loocv by",grpvar))
  print(varp) 
  dev.off()
  
  
  pdf(paste(grpvar,"dotplot loocv by count.pdf"),
      width = 12,height = 5)
  
  varp = ggplot(data = stackcv,aes(x = count,y = loocv,colour = model))+
    geom_point(position = position_dodge(width = 0.5),size = 0.15)+
    #scale_y_continuous(limits = c(0,max(comp.var$loocv_uci)))+
    facet_wrap(~model)+
    scale_fill_brewer(palette="Dark2")+
    scale_color_brewer(palette="Dark2")+
    labs(title = paste(sp,"loocv by",grpvar))
  print(varp) 
  dev.off()
  
  
  
}#grpvar





}#temp end to sp loop

for(m1 in c("GAM","GAMYE","FD")){
  for(m2 in c("Standard","FD","GAMYE")){
    if(m1 == m2){next}
    cmpn = paste0(m1,"-",m2)
    cmpnv = paste0(m1,"-",m2,"v")
cv[,cmpn] = -2*(cv[,paste0(m1,"_loocv")])-(-2*(cv[,paste0(m2,"_loocv")]))
cv[,cmpnv] = exp(cv[,paste0(m1,"_loocv")])/exp(cv[,paste0(m2,"_loocv")])/(1+(exp(cv[,paste0(m1,"_loocv")])/exp(cv[,paste0(m2,"_loocv")])))


luni = function(x){length(unique(x))}
bystrat = tapply(cv[,cmpn],cv[,"strat.name"],mean)
sdbystrat = tapply(cv[,cmpn],cv[,"strat.name"],sd)
nrtsbystrat = tapply(cv[,"route"],cv[,"strat.name"],luni)
ncountsbystrat = tapply(cv[,"route"],cv[,"strat.name"],length)
nyrsbystrat = tapply(cv[,"year"],cv[,"strat.name"],luni)
skewbystrat = tapply(cv[,"count"],cv[,"strat.name"],skewness)

comp.strat = data.frame(strat = names(bystrat),
                  cmpn = bystrat,
                  cmpn_se = sdbystrat/sqrt(ncountsbystrat),
                  nrts = nrtsbystrat,
                  ncounts = ncountsbystrat,
                  nyrs = nyrsbystrat,
                  skew = skewbystrat)
names(comp.strat)[c(2,3)] <- c(cmpn,paste0(cmpn,"_se"))

byyear = tapply(cv[,cmpn],cv[,"rYear"],mean)
skewy = tapply(cv[,"count"],cv[,"rYear"],skewness)
sdbyyear = tapply(cv[,cmpn],cv[,"rYear"],sd)
nrtsbyyear = tapply(cv[,"route"],cv[,"rYear"],luni)

loobyyear = data.frame(loo = byyear,
                       lci = byyear-(2*(sdbyyear/sqrt(nrtsbyyear))),
                       uci = byyear+(2*(sdbyyear/sqrt(nrtsbyyear))),
                       rYear = as.integer(names(byyear)),
                       skew = skewy,
                       scale = "all")


############################# strata by year loocv complilation


bystrat.y = tapply(cv[,cmpn],cv[,c("strat.name","rYear")],mean,na.rm = T)
sdbystrat.y = tapply(cv[,cmpn],cv[,c("strat.name","rYear")],sd,na.rm = T)
skewbystrat.y = tapply(cv[,"count"],cv[,c("strat.name","rYear")],skewness,na.rm = T)
#nrtsbystrat.y = tapply(cv[,"route"],cv[,c("strat.name","rYear")],luni)
ncountsbystrat.y = tapply(cv[,"route"],cv[,c("strat.name","rYear")],length)


bystrat.yN = tapply(cv[,"count"],cv[,c("strat.name","rYear")],median,na.rm = T)
bystrat.yNq1 = tapply(cv[,"count"],cv[,c("strat.name","rYear")],quantile,probs = 0.25,na.rm = T)
bystrat.yNq3 = tapply(cv[,"count"],cv[,c("strat.name","rYear")],quantile,probs = 0.75,na.rm = T)

by.yN = tapply(cv[,"count"],cv[,c("rYear")],mean,na.rm = T)
by.yNq1 = tapply(cv[,"count"],cv[,c("rYear")],quantile,probs = 0.25,na.rm = T)
by.yNq3 = tapply(cv[,"count"],cv[,c("rYear")],quantile,probs = 0.75,na.rm = T)

##############################
## strata trajectory plots
if(plotstrat){
pdf(paste0("strata level comparison plots ",cmpn,".pdf"),
    height = 11,
    width = 8.5)
for(st in unique(index$strat.name)){
  
  tmp = index[which(index$strat.name == st & index$model %in% c(m1,m2)),]
data_summary <- data.frame(Year = tmp$rYear,
                           Index = tmp$index,
                           Index_neg = tmp$lci,
                           Index_pos = tmp$uci,
                           Model = tmp$model)

data_summary$Model <- factor(data_summary$Model, levels = models)


N = bystrat.yN[st,]
Nl =  bystrat.yNq1[st,]
Nu =bystrat.yNq3[st,]

Ndf = data.frame(mean = as.numeric(N),
                   lci = as.numeric(Nl),
                   uci = as.numeric(Nu),
                   Year = 1966:2016)

###################################################
# Plot
###################################################

ptraj <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste("Population trajectories", st), x = "Year", y = "Index") + 
  geom_pointrange(data = Ndf, aes(x = Year, y = mean,ymin = lci,ymax = uci),alpha = 0.3) +
  geom_line(data = data_summary, aes(x = Year, y = Index, colour = Model)) +
  geom_ribbon(data = data_summary, aes(x = Year, ymin = Index_neg, ymax = Index_pos, fill = Model), alpha = 0.12)+
  scale_y_continuous(limits = c(0, max(data_summary$Index_pos)), expand = c(0,0), oob = squish)+
  scale_x_continuous(limits = c(1966, 2016), oob = squish)





tmpcv = cv[which(cv$strat.name == st),c("rYear",cmpn)]
names(tmpcv)[2] = "comp"

tmploo = bystrat.y[st,]
tmpsk = skewbystrat.y[st,]
tmploosd = sdbystrat.y[st,]
tmplool = tmploo-(0.1*tmploosd)
tmplool[which(is.na(tmplool))] <- tmploo[which(is.na(tmplool))]
tmploou = tmploo+(0.1*tmploosd)
tmploou[which(is.na(tmploou))] <- tmploo[which(is.na(tmploou))]
tmpdf = data.frame(loo = as.numeric(tmploo),
                   lci = as.numeric(tmplool),
                   uci = as.numeric(tmploou),
                   rYear = 1966:2016,
                   skew = tmpsk,
                   scale = "strat")
tmpdf = rbind(tmpdf,loobyyear)
#nrtsbystrat.y[st,]

ploobox <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste(cmpn,"looIC comparison", st), x = "Year", y = "Dif looIC") +
  annotate("text",hjust = "left",label = paste("favours",c(m1,m2)),
           y = c(quantile(tmpcv$comp,0.05),quantile(tmpcv$comp,0.95)),x = rep(1966,2))+
  geom_boxplot(data = tmpcv,varwidth = T,aes(x = rYear,y = comp, group = rYear),fill = "grey",colour = grey(0.2))+
  geom_hline(aes(yintercept = 0))+
  geom_col(position = "dodge")+
  geom_linerange(data = tmpdf, 
                 position = position_dodge(width = 0.4),alpha = 0.7,
                 aes(x = rYear, ymin = lci, ymax = uci, group = scale,colour = scale))+
  geom_point(data = tmpdf, 
             position = position_dodge(width = 0.4),alpha = 0.7,
             aes(x = rYear,y = loo, group = scale,colour = scale))+
  scale_y_continuous(limits = c(quantile(tmpcv$comp,0.02),quantile(tmpcv$comp,0.98)), expand = c(0,0), oob = squish)+
  scale_x_continuous(limits = c(1966, 2016), oob = squish)



# 
# ploo <- ggplot() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top") + 
#   labs(title = paste(cmpn,"looIC mean +- SE", st), x = "Year", y = "Dif looIC") +
#   annotate("text",hjust = "left",label = paste("favours",c(m1,m2)),
#            y = c(-1.1,1.1),x = rep(1966,2))+
#   geom_col(position = "dodge")+
#   geom_linerange(data = tmpdf, 
#                   position = position_dodge(width = 0.4),
#                   aes(x = rYear, ymin = lci, ymax = uci, group = scale,colour = scale))+
#   geom_point(data = tmpdf, 
#                   position = position_dodge(width = 0.4),
#                   aes(x = rYear,y = loo, group = scale,colour = scale))+
#   geom_hline(aes(yintercept = 0))+
# scale_y_continuous(limits = c(-1.2,1.2), expand = c(0,0), oob = squish)+
#   scale_x_continuous(limits = c(1966, 2016), oob = squish)



tmpncnt = ncountsbystrat.y[st,]
ncts = data.frame(rYear = 1966:2016,
                  ncnts = tmpncnt)

pn <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste("number of counts by year", st), x = "Year", y = "Number of counts") +
  geom_col(data = ncts,aes(x = rYear,y = ncnts))+
  scale_x_continuous(limits = c(1966, 2016), oob = squish)

psk <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste("skewness of loocv by year", st), x = "Year", y = "Skewness") +
  geom_point(data = tmpdf,aes(x = skew,y = loo))
  # scale_x_continuous(limits = c(1966, 2016), oob = squish)

plts = list(ptraj,ploobox,pn,psk)
lyt = matrix(ncol = 1,c(1,1,2,2,2,3,4))
multiplot(plotlist = plts,layout = lyt)

}## end strata loop

dev.off()
}#end if plot strata


pdf(paste0(cmpn,"overall comparisons.pdf"),
    height = 11,
    width = 8.5)
# plot(x = (cv$count), y = cv[,cmpn],
#      main = cmpn)


tmp = index.cont[which(index.cont$Model %in% c(m1,m2)),]
data_summary <- data.frame(Year = tmp$rYear,
                           Index = tmp$Index,
                           Index_neg = tmp$Index_neg,
                           Index_pos = tmp$Index_pos,
                           Model = tmp$Model)

data_summary$Model <- factor(data_summary$Model, levels = models)


N = by.yN
Nl =  by.yNq1
Nu =by.yNq3

Ndf = data.frame(mean = as.numeric(N),
                 lci = as.numeric(Nl),
                 uci = as.numeric(Nu),
                 Year = 1966:2016)

# Ndf = loobyyear
# Ndf$Year = Ndf$rYear


###################################################
# Plot
###################################################

ptraj <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste("Population trajectories continental"), x = "Year", y = "Index") + 
  geom_point(data = Ndf, aes(x = Year, y = mean),alpha = 0.3) +
  geom_line(data = data_summary, aes(x = Year, y = Index, colour = Model)) +
  geom_ribbon(data = data_summary, aes(x = Year, ymin = Index_neg, ymax = Index_pos, fill = Model), alpha = 0.12)+
coord_cartesian(ylim = c(0, max(Ndf$mean)),  xlim = c(1966, 2016))





tmpcv = cv[,c("rYear",cmpn)]
names(tmpcv)[2] = "comp"


tmpdf = loobyyear

# ploobox <- ggplot() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top") + 
#   labs(title = paste(cmpn,"looIC comparison"), x = "Year", y = "Dif looIC") +
#   annotate("text",hjust = "left",label = paste("favours",c(m1,m2)),
#            y = c(quantile(tmpcv$comp,0.05),quantile(tmpcv$comp,0.95)),x = rep(1966,2))+
#   geom_boxplot(data = tmpcv,varwidth = T,aes(x = rYear,y = comp, group = rYear),fill = "grey",colour = grey(0.2))+
#   geom_hline(aes(yintercept = 0))+
#   geom_col(position = "dodge")+
#   geom_linerange(data = tmpdf, 
#                  position = position_dodge(width = 0.4),alpha = 0.7,colour = "red",
#                  aes(x = rYear, ymin = lci, ymax = uci))+
#   geom_point(data = tmpdf, 
#              position = position_dodge(width = 0.4),alpha = 0.7,colour = "red",
#              aes(x = rYear,y = loo))+
#   coord_cartesian(ylim = c(quantile(tmpcv$comp,0.02),quantile(tmpcv$comp,0.98)),  xlim = c(1966, 2016))


##### means and ses plot

ploolrng <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = paste(cmpn,"looIC comparison"), x = "Year", y = "Dif looIC") +
  annotate("text",hjust = "left",label = paste("favours",c(m1,m2)),
           y = c(quantile(tmpdf$lci,0.01),quantile(tmpdf$uci,0.99)),
           x = rep(1990,2))+
  #geom_boxplot(data = tmpcv,varwidth = T,aes(x = rYear,y = comp, group = rYear),fill = "grey",colour = grey(0.2))+
  geom_hline(aes(yintercept = 0))+
  geom_col(position = "dodge")+
  geom_linerange(data = tmpdf, 
                 position = position_dodge(width = 0.4),alpha = 0.7,colour = "red",
                 aes(x = rYear, ymin = lci, ymax = uci))+
  geom_point(data = tmpdf, 
             position = position_dodge(width = 0.4),alpha = 0.7,colour = "red",
             aes(x = rYear,y = loo))+
  coord_cartesian(ylim = c(min(c(quantile(tmpdf$lci,0.01),min(tmpdf$loo))),max(c(quantile(tmpdf$uci,0.99),max(tmpdf$loo)))),  xlim = c(1966, 2016))









psky <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = "by year", y = paste(cmpn,"loocv"), x = "Skewness") +
  geom_hline(aes(yintercept = 0))+
  geom_point(data = tmpdf,aes(x = skew,y = loo))+
  geom_smooth(method = lm,data = tmpdf,aes(x = skew,y = loo),alpha = 0.2)
  





  tmpsk = comp.strat[,c("strat",cmpn,"skew")]
names(tmpsk)[2] = "loo"
  psks <- ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top") + 
  labs(title = "by strata", y = paste(cmpn,"loocv"), x = "Skewness") +
    geom_hline(aes(yintercept = 0))+
  geom_point(data = tmpsk,aes(x = skew,y = loo))+
    geom_smooth(method = lm,data = tmpsk,aes(x = skew,y = loo),alpha = 0.2)
  
# 
# ploo <- ggplot() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top") + 
#   labs(title = paste(cmpn,"looIC mean +- SE", st), x = "Year", y = "Dif looIC") +
#   annotate("text",hjust = "left",label = paste("favours",c(m1,m2)),
#            y = c(-1.1,1.1),x = rep(1966,2))+
#   geom_col(position = "dodge")+
#   geom_linerange(data = tmpdf, 
#                   position = position_dodge(width = 0.4),
#                   aes(x = rYear, ymin = lci, ymax = uci, group = scale,colour = scale))+
#   geom_point(data = tmpdf, 
#                   position = position_dodge(width = 0.4),
#                   aes(x = rYear,y = loo, group = scale,colour = scale))+
#   geom_hline(aes(yintercept = 0))+
# scale_y_continuous(limits = c(-1.2,1.2), expand = c(0,0), oob = squish)+
#   scale_x_continuous(limits = c(1966, 2016), oob = squish)


# 
# tmpncnt = ncountsbystrat.y[st,]
# ncts = data.frame(rYear = 1966:2016,
#                   ncnts = tmpncnt)
# 
# pn <- ggplot() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top") + 
#   labs(title = paste("number of counts by year", st), x = "Year", y = "Number of counts") +
#   geom_col(data = ncts,aes(x = rYear,y = ncnts))+
#   scale_x_continuous(limits = c(1966, 2016), oob = squish)

# plts = list(ptraj,ploobox,ploolrng,psks,psky)
# lyt = matrix(ncol = 2,c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,5,4,5),byrow = T)
# multiplot(plotlist = plts,layout = lyt)

plts = list(ptraj,ploolrng,psks,psky)
lyt = matrix(ncol = 2,c(1,1,1,1,1,1,2,2,2,2,2,2,3,4,3,4),byrow = T)
multiplot(plotlist = plts,layout = lyt)





















plot(comp.strat[,cmpn],x = sqrt(comp.strat[,"ncounts"]),
     main = cmpn,
     ylim = c(-0.5,0.5),
     xlim = c(0,50))
arrows(x0 = sqrt(comp.strat$ncounts),
       x1 = sqrt(comp.strat$ncounts),
       y0 = comp.strat[,cmpn]-comp.strat[,paste0(cmpn,"_se")],
       y1 = comp.strat[,cmpn]+comp.strat[,paste0(cmpn,"_se")],
       length = 0,
       col = grey(0.8))
abline(h = 0)
extrem = which(abs(comp.strat[,cmpn]) > quantile(abs(comp.strat[,cmpn]),0.9))
text(comp.strat[extrem,"strat"],
     x = sqrt(comp.strat[extrem,"ncounts"]),
     y = comp.strat[extrem,cmpn],
     pos = 4,
     cex = 0.6)

### above suggests that regions that favour the standard model have highly skewed distributions
### if the means are different from the medians, they're more likely to be positive
# x11()
# plot(cv[,cmpn],x = sqrt(cv[,"count"]),
#      main = cmpn)
# sm = loess(cv$cmpn~sqrt(cv$count),data = cv)
# smpl = predict(sm,newdata = data.frame(count = seq(0,10000,length.out = 1000)))
# lines(smpl,x = sqrt(seq(0,10000,length.out = 1000)),col = "red")

# "Nova Scotia Prince Edward Island-BCR14"
# "CA-NSPE-14"

mp = readOGR(dsn = "c:/bbs/2017/input", layer = "BBS_Analytical_Strata_2015_AlbersNA")
regs = read.csv("c:/bbs/2017/input/Strata_US_Can_new_strata_names.csv",stringsAsFactors = F)
regs[which(regs$prov == "NS"),"prov"] <- "NSPE"
sabrev = read.csv("c:/bbs/2017/input/state abrev.csv",stringsAsFactors = F)
# sabrev[which(sabrev$code == "NS"),]
regs = merge(regs,sabrev,by.x = c("country","prov"),
             by.y = c("country","code"))

regs$strat = paste0(regs$State,"-BCR",regs$bcr)
regs[which(regs$strat == "Nova Scotia-BCR14"),"strat"] <- "Nova Scotia Prince Edward Island-BCR14"

comp.strat = merge(comp.strat,regs[,c("strat","St_12")],all.x = T)

mp@data = merge(mp@data,comp.strat,by = "St_12",all.x = T)
mp@data$col =  transp.func("purple",0.5)
mp@data$col[which(mp@data[,cmpn] > 0)] = transp.func("darkorange",0.5)
mp@data$col[which(is.na(mp@data[,cmpn]))] = "transparent"

mp@data$col[which(mp@data[,cmpn]+mp@data[,paste0(cmpn,"_se")] < 0)] = "purple"

mp@data$col[which(mp@data[,cmpn]-mp@data[,paste0(cmpn,"_se")] > 0)] = "darkorange"
mp@data$col[which(is.na(mp@data[,cmpn]))] = "transparent"


cents = coordinates(mp)


plot(mp,col = mp@data$col,border = grey(0.7),
     main = cmpn)
text(x = cents[,1],
     y = cents[,2],
     mp@data$St_12,
     cex = 0.3)




# for(tt in c(1966,1970,2006)){
#   ye = 2016
#     ys = tt
#     ysr = c((tt-2):(tt+2))
#     yer = c(2014:2016)
#     
#     cvtmp = cv[which(cv$rYear %in% c(ye,ys)),]
# 
# }



# par(mar = c(4,12,1,1))
# boxplot(cv[,cmpn]~cv[,"strat.name"],
#         ylim = c(-2,2),
#         horizontal = T,
#         pars = list(las = 1,
#                     yaxs = "i",
#                     cex.axis = 0.5),
#         main = cmpn)
# abline(v = 0)
# points(x = comp.strat[,cmpn],
#        y = 1:nrow(comp.strat),
#        col = "red")

dev.off()

pdf(paste0("trajectory comparison ",cmpn,".pdf"))

print(ptraj)

dev.off()


  }
}
 



for(tt in c(1966,1970,2006)){
  ye = 2016
  ys = tt
  ysr = c((tt-2):(tt+2))
  if(tt == 1966){ysr = c(1966:1970)}
  yer = c(2012:2016)
  yall = c(tt:2016)
  

for(m1 in c("GAM","GAMYE","FD")){
  for(m2 in c("Standard","FD","GAMYE")){
    if(m1 == m2){next}
    cmpn = paste0(m1,"-",m2)
    cmpnv = paste0(m1,"-",m2,"v")
    tmpo = data.frame(comp = rep(cmpn,6),stringsAsFactors = F)
    
    cvyeys = cv[which(cv$rYear %in% c(ye,ys)),cmpn]
    cvvyeys = cv[which(cv$rYear %in% c(ye,ys)),cmpnv]
    cvys = cv[which(cv$rYear %in% c(ys)),cmpn]
    cvvys = cv[which(cv$rYear %in% c(ys)),cmpnv]
    cvyer = cv[which(cv$rYear %in% c(yer)),cmpn]
    cvvyer = cv[which(cv$rYear %in% c(yer)),cmpnv]
    cvysr = cv[which(cv$rYear %in% c(ysr)),cmpn]
    cvvysr = cv[which(cv$rYear %in% c(ysr)),cmpnv]
    cvyerysr = cv[which(cv$rYear %in% c(ysr,yer)),cmpn]
    cvvyerysr = cv[which(cv$rYear %in% c(ysr,yer)),cmpnv]
    cvyall = cv[which(cv$rYear %in% c(yall)),cmpn]
    cvvyall = cv[which(cv$rYear %in% c(yall)),cmpnv]
    tmplst = list(cvyeys = cvyeys,
                  cvvyeys = cvvyeys,
                  cvys = cvys,
                  cvvys = cvvys,
                  cvyer = cvyer,
                  cvvyer = cvvyer,
                  cvysr = cvysr,
                  cvvysr = cvvysr,
                  cvyerysr = cvyerysr,
                  cvvyerysr = cvvyerysr,
                  cvyall = cvyall,
                  cvvyall = cvvyall)
    j = 0
    for(ccc in c("yeys","ys","yer","ysr","yerysr","yall")){
      j = j+1
    tmpo[j,paste0("m_cv")] = mean(tmplst[[paste0("cv",ccc)]])
    tmpo[j,paste0("lci_cv")] = mean(tmplst[[paste0("cv",ccc)]])-(1.96*(sd(tmplst[[paste0("cv",ccc)]])/sqrt(length(tmplst[[paste0("cv",ccc)]]))))
    tmpo[j,paste0("uci_cv")] = mean(tmplst[[paste0("cv",ccc)]])+(1.96*(sd(tmplst[[paste0("cv",ccc)]])/sqrt(length(tmplst[[paste0("cv",ccc)]]))))
    
    tmpo[j,paste0("sum_cv")] = sum(tmplst[[paste0("cv",ccc)]])
    tmpo[j,paste0("sumlci_cv")] = sum(tmplst[[paste0("cv",ccc)]])-(1.96*(sd(tmplst[[paste0("cv",ccc)]])*sqrt(length(tmplst[[paste0("cv",ccc)]]))))
    tmpo[j,paste0("sumuci_cv")] = sum(tmplst[[paste0("cv",ccc)]])+(1.96*(sd(tmplst[[paste0("cv",ccc)]])*sqrt(length(tmplst[[paste0("cv",ccc)]]))))
    
    tmpo[j,paste0("m_cvv")] = mean(tmplst[[paste0("cvv",ccc)]])
    tmpo[j,paste0("lci_cvv")] = mean(tmplst[[paste0("cvv",ccc)]])-(1.96*(sd(tmplst[[paste0("cvv",ccc)]])/sqrt(length(tmplst[[paste0("cv",ccc)]]))))
    tmpo[j,paste0("uci_cvv")] = mean(tmplst[[paste0("cvv",ccc)]])+(1.96*(sd(tmplst[[paste0("cvv",ccc)]])/sqrt(length(tmplst[[paste0("cv",ccc)]]))))
    tmpo[j,paste0("pM1_cvv")] = length(which(tmplst[[paste0("cvv",ccc)]]>0.5))/length((tmplst[[paste0("cvv",ccc)]]))
    tmpo[j,paste0("N_cvv")] = length((tmplst[[paste0("cvv",ccc)]]))
    tmpo[j,"yrt"] = ccc
    }
    tmpo$startyear = tt
    
    
    if(m1 == "GAM" & m2 == "Standard" & tt == 1966){
      trendyears.summary = tmpo
      
    }else{
      trendyears.summary = rbind(trendyears.summary,tmpo)
    }
    
    
    

}}}

 write.csv(trendyears.summary,paste0(sp," trend year summary.csv"),row.names = F)
}#end species loop
# 
# 
# ptraj <- ggplot() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top") + 
#   labs(title = paste("Population trajectories continental"), x = "Year", y = "Index") + 
#   geom_point(data = Ndf, aes(x = Year, y = mean),alpha = 0.3) +
#   geom_line(data = data_summary, aes(x = Year, y = Index, colour = Model)) +
#   geom_ribbon(data = data_summary, aes(x = Year, ymin = Index_neg, ymax = Index_pos, fill = Model), alpha = 0.12)+
#   scale_y_continuous(limits = c(0, max(Ndf$mean)), expand = c(0,0), oob = squish)+
#   scale_x_continuous(limits = c(1966, 2016), oob = squish)
# 
# 
# 
# 
# 
# tmpcv = cv[,c("rYear",cmpn)]
# names(tmpcv)[2] = "comp"
# 
# 
# tmpdf = trendyears.summary
# 
# ploobox <- ggplot() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = "top") + 
#   labs(title = paste(cmpn,"looIC comparison"), x = "Year", y = "Dif looIC") +
#   annotate("text",hjust = "left",label = paste("favours",c("1st model","2nd model")),
#            y = c(max(tmpdf$lci_cv),max(tmpdf$uci_cv)),x = rep(1966,2))+
#   geom_hline(aes(yintercept = 0))+
#   geom_col(position = "dodge")+
#   geom_linerange(data = tmpdf, 
#                  position = position_dodge(width = 0.4),alpha = 0.7,
#                  aes(x = comp, ymin = lci_cv, ymax = uci_cv))+
#   geom_point(data = tmpdf, 
#              position = position_dodge(width = 0.4),alpha = 0.7,
#              aes(x = comp,y = m_cv))
#   # scale_y_continuous(limits = c(quantile(tmpcv$comp,0.02),quantile(tmpcv$comp,0.98)), expand = c(0,0), oob = squish)+
#   # scale_x_continuous(limits = c(1966, 2016), oob = squish)
# 
# 
# 
# 
# 
# 
# plts = list(ptraj,ploobox,psks,psky)
# lyt = matrix(ncol = 2,c(1,1,1,1,1,1,2,2,2,2,2,2,3,4,3,4),byrow = T)
# multiplot(plotlist = plts,layout = lyt)
# 









####################################################













