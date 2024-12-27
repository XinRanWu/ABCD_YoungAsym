load("E:/Acdamic/R_space/YoungAsym.RData")
require("pacman")
pacman::p_load(
  cluster,factoextra,psych,ggradar,pheatmap,ggplot2,ggiraph,
  ggiraphExtra,RColorBrewer,tableone,glmnet,rgl,ggseg,dplyr,
  reshape2,Rmisc,ggpubr,export,Rcpp,ggridges,stringr,ggfun,
  gghalves,utils,readr,parallel,lmerTest,foreach,tidyr,ggrepel,
  ggsci,mgcv,doParallel
)
################################ Figure. 1 ####################################

DK.info.index <- c(31,13,9,21,27,25,19,29,15,23,1,24,4,30,26,11,6,2,5,22,16,14,10,20,12,7,8,18,32,17,3,28,33,34)

DK.Lobe <- c(rep(c("Frontal","Frontal","Frontal","Frontal","Frontal","Frontal","Frontal","Frontal","Frontal",
                   "Parietal","Parietal","Cingulate","Cingulate","Cingulate","Cingulate","Parietal","Parietal",
                   "Parietal","Parietal","Parietal","Occipital","Occipital","Occipital","Occipital","Temporal",	
                   "Temporal","Temporal","Temporal","Temporal","Temporal","Temporal","Temporal","Temporal","Insula",
                   "Subcort","Subcort","Subcort","Subcort","Subcort","Subcort","Subcort"),2),"Subcort")


DK.info <- readxl::read_xlsx("E:/Acdamic/Acdamic_Toolset/hansen_receptors-main/data/lausanne/region_info_scale033.xlsx")
DK.info <- DK.info[,-1]
DK.info$lobe <- DK.Lobe
DK.info.cortex <- DK.info[DK.info$structure=="cortex",]
DK.info.cortex <- DK.info.cortex[c(DK.info.index,DK.info.index+34),]
yeo_7.id <- order(DK.info.cortex$von_economo[1:34])


################################ Figure. 2 ####################################


# data$brainregion <- gsub("Vol_","",data$brainregion)
data <- 
  aggregate(. ~ AgeGroup, Young_Asym[,c(86,76:85)], mean) %>% 
  gather(variable, value, -AgeGroup)
# data$variable <- factor(data$variable, 
#levels = names(Young_Asym[,c(76:85)]),
#ordered = TRUE,)
ggplot(data, aes(x = AgeGroup, y = value, group = variable, colour = variable)) +
  geom_line() + geom_point() + theme_minimal() + 
  scale_color_brewer(palette = "Paired") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid")#  + scale_color_brewer()


data <- 
  aggregate(. ~ AgeGroup, Young_Asym[,c(86,76:85)], mean) %>% 
  gather(variable, value, -AgeGroup)


# fit.abcdAsymBL <- data.frame()
# for (i in 8:85){
#      fit.ttest <- t.test(Young_Asym[Young_Asym$AgeGroup=="10",i])
#      fit.abcdAsymBL[i-7,"brain_var"] <- names(Young_Asym)[i]
#      fit.abcdAsymBL[i-7,c("t", "df", "p")] <- 
#        fit.ttest[c("statistic","parameter","p.value")]
# }
# 
# fit.abcdAsymBL[fit.abcdAsymBL$t > 0 & fit.abcdAsymBL$p<0.05,]
# which(fit.abcdAsymBL$t > 0 & fit.abcdAsymBL$p<0.05)



c(1,6,7,9,11,14,15,16,24,25,26,27,29,30,31,32,34)
c(36,37,45,49,50,51,52,53,54,55,56,57,59,60,61,64,65)
c(69,71,74,78)

c(2,3,4,5,10,12,13,17,18,20,21,22,28,33)
c(35,38,39,40,41,42,44,46,47,48,58,63,66,67,68)
c(70,72,73,75,76,77)




Young_Asym$Race[Young_Asym$Race!="White"&
                  Young_Asym$Race!="Black"&
                  !is.na(Young_Asym$Race)] <- "Other/Mixed"
Young_Asym$Race <- as.character(Young_Asym$Race)

Young_Freesurfer$Race[Young_Freesurfer$Race!="White"&
                  Young_Freesurfer$Race!="Black"&
                    !is.na(Young_Freesurfer$Race)] <- "Other/Mixed"
Young_Freesurfer$Race <- as.character(Young_Freesurfer$Race)



# ?????켣
fit.YoungAsym.Age = data.frame()
fit.YoungAsym.Gender = data.frame()
fit.YoungAsym.Handedness = data.frame()
fit.YoungAsym.ITV = data.frame()
# fit.YoungAsym.Mean = data.frame()
for (i in 10:87){
  f = as.formula(paste(names(Young_Asym)[i],
                       " ~ 1 + age + gender + Handedness + Race +",
                       "Vol_EstimatedTotalIntraCranialVol + ",
                       "(1|Dataset:ImagingCentreID) + (1|SubID)",sep = ""))
  
  res.lmer = lmerTest::lmer(f,Young_Asym)
  sum.lmer = summary(res.lmer)
  
  fit.YoungAsym.Age[i-9,c("demovar")] <- "age"
  fit.YoungAsym.Age[i-9,c("brainregion")] = names(Young_Asym)[i]
  fit.YoungAsym.Age[i-9,c("estimate","df", "t", "p")] =
    sum.lmer$coefficients["age", c("Estimate","df", "t value", "Pr(>|t|)")]
  
  fit.YoungAsym.Gender[i-9,c("demovar")] <- "gender"
  fit.YoungAsym.Gender[i-9,c("brainregion")] = names(Young_Asym)[i]
  fit.YoungAsym.Gender[i-9,c("estimate","df", "t", "p")] =
    sum.lmer$coefficients["genderM", c("Estimate","df", "t value", "Pr(>|t|)")]
  
  fit.YoungAsym.Handedness[i-9,c("demovar")] <- "Handedness"
  fit.YoungAsym.Handedness[i-9,c("brainregion")] = names(Young_Asym)[i]
  fit.YoungAsym.Handedness[i-9,c("estimate","df", "t", "p")] =
    sum.lmer$coefficients["Handedness", c("Estimate","df", "t value", "Pr(>|t|)")]
  
  fit.YoungAsym.ITV[i-9,c("demovar")] <- "ITV"
  fit.YoungAsym.ITV[i-9,c("brainregion")] = names(Young_Asym)[i]
  fit.YoungAsym.ITV[i-9,c("estimate","df", "t", "p")] =
    sum.lmer$coefficients["Vol_EstimatedTotalIntraCranialVol", c("Estimate", "df","t value", "Pr(>|t|)")]
  
  print(paste(names(Young_Asym)[i],' ~ Age + Sex + Handedness + ITV + Race', sep = ""))
}

fit.YoungAsym.DemoVar <- rbind(fit.YoungAsym.Age,fit.YoungAsym.Gender,
                               fit.YoungAsym.Handedness,fit.YoungAsym.ITV)
rm(fit.YoungAsym.Age,fit.YoungAsym.Gender,
   fit.YoungAsym.Handedness,fit.YoungAsym.ITV)

rm(list = c("res.lmer","sum.lmer"))

# p.fdr <- p.adjust(fit.YoungAsym.Age$p,method = "fdr")

## Interaction of Gender * Age
# ?????켣
fit.YoungAsym.AgeGender = data.frame()
for (i in 8:87){
  f = as.formula(paste(names(Young_Asym)[i],
                       " ~ age + gender + age*gender + Handedness + ",
                       "Vol_EstimatedTotalIntraCranialVol + ",
                       "(1|Dataset:ImagingCentreID) + (1|SubID)",sep = ""))
  res.lmer = lmerTest::lmer(f,Young_Asym)
  sum.lmer = summary(res.lmer)
  fit.YoungAsym.AgeGender[i-7,c("brainregion")] = names(Young_Asym)[i]
  fit.YoungAsym.AgeGender[i-7,c("estimate","df", "t", "p")] =
    sum.lmer$coefficients["age:genderM", c("Estimate","df", "t value", "Pr(>|t|)")]
  print(paste(names(IMAGEN_Asym)[i],' ~ Age : Sex', sep = ""))
}


fit.YoungAsym.DemoVar$Modal <- NA
fit.YoungAsym.DemoVar$Modal[which(str_detect(fit.YoungAsym.DemoVar$brainregion,"CT_"))] <- "CT"
fit.YoungAsym.DemoVar$Modal[which(str_detect(fit.YoungAsym.DemoVar$brainregion,"CSA_"))] <- "CSA"
fit.YoungAsym.DemoVar$Modal[which(str_detect(fit.YoungAsym.DemoVar$brainregion,"Vol_"))] <- "Vol"

behav_dic <- unique(fit.YoungAsym.DemoVar$demovar)
# behav_dic <- behav_dic[-1]
fit.YoungAsym.DemoVar.fdr <- fit.YoungAsym.DemoVar
for (i in 1:length(behav_dic)){
  fit.YoungAsym.DemoVar.fdr$pval[fit.YoungAsym.DemoVar.fdr$demovar == behav_dic[i]] <- 
    p.adjust(fit.YoungAsym.DemoVar.fdr[fit.YoungAsym.DemoVar.fdr$demovar == behav_dic[i],"pval"],method = "fdr")
}

fit.YoungAsym.DemoVar.Vol <- fit.YoungAsym.DemoVar.fdr[which(fit.YoungAsym.DemoVar.fdr$Modal=="Vol"),]
fit.YoungAsym.DemoVar.Vol[fit.YoungAsym.DemoVar.Vol$p<0.05,]
fit.YoungAsym.DemoVar.fdr <- fit.YoungAsym.DemoVar.fdr[-1*which(fit.YoungAsym.DemoVar.fdr$Modal=="Vol"),]
# fit.YoungAsym.DemoVar.fdr <- fit.YoungAsym.DemoVar.fdr[-1*which(fit.YoungAsym.DemoVar.fdr$demovar=="age"),]

# 600 300
# 800 300
# PosAgeEff9
# t[which(t==0)] <- NA
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic)),
                  t = fit.YoungAsym.DemoVar.fdr$t*(fit.YoungAsym.DemoVar.fdr$p<0.05), 
                  behav = fit.YoungAsym.DemoVar.fdr$demovar,
                  modal = fit.YoungAsym.DemoVar.fdr$Modal,
                  hemi = rep(c(rep("left", 34*length(behav_dic))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic)),
                     t = NA,
                     modal = rep(c("CSA","CT"),each = 34*length(behav_dic)*2),
                     behav = rep(behav_dic,each = 2),
                     hemi = rep(c(rep("left", 2*length(behav_dic))), 1))
data = rbind(data,data.NA);

for (i in c("age","gender","Handedness","ITV")){
  for (j in c('CSA','CT')){
    data2 <- data[data$behav==i&data$modal==j,]
    clim = ceiling( max(abs(data2$t),na.rm = T)*10)/10
    data2 %>%
      group_by(modal) %>% ggseg(atlas = dk,colour = "black",hemisphere = "left",
                                mapping = aes(fill = t),position = "stacked",size = .5) +
      theme(legend.position = "right") +
      scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                           na.value = "lightgrey", breaks = seq(-1*clim, 1, clim),limits = c(-1*clim, clim)) +
      theme_brain() + theme(text = element_text(family = "Arial", size = 15))
    
    ggsave(file = paste('DemoVar-',i,'-',j,'.png',sep = ''),
           device = "png",width = 160,height = 60,units = c("mm"))
  }
}

clim = 28
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_grid(modal~behav) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 3),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))
## 1100 550




YYoung_Asym$EstimatedTotalIntraCranialVol <- Young_Asym$EstimatedTotalIntraCranialVol/100000
residuals <- matrix(nrow = nrow(Young_Asym), ncol = 78)
for (i in 1:78){
  f <- as.formula(paste(names(Young_Asym)[i+9], "~ gender + Handedness + (1|Dataset:ImagingCentreID)", sep = ""))
  model <- lmer(f, data = Young_Asym)
  mm <- model.matrix(model)
  fixed_effects <- fixef(model) 
  predicted_vals <- mm %*% fixed_effects
  residuals[, i] <- Young_Asym[, i+9] - predicted_vals
}

for (i in 1:78){
  f <- as.formula(paste(names(Young_Asym)[i+9], "~ gender + Handedness + ImagingCentreID", sep = ""))
  model <- lm(f, data = Young_Asym)
  residuals[, i] <- residuals(model)
}
#  + EstimatedTotalIntraCranialVol

residuals <- data.frame(residuals)
names(residuals) <- names(Young_Asym)[10:87]
Young_AsymResiduals <- cbind(Young_Asym[,1:9],residuals)



g_sumTb <- Young_Asym[,c(Residuals8,10:87)] %>% 
  group_by(AgeGroup) %>% 
  summarise_all(list(mean = mean, sd = sd),na.rm = T)

df_m

data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),1),
                  r = colMeans(Young_Asym[Young_Asym$AgeGroup==19,c((1:34)+9)]),hemi = rep(c(rep("left", 34)), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), 1),
                     r = NA,hemi = rep(c(rep("left", 2)), 1))
data = rbind(data,data.NA);

max(data$r,na.rm = T)
min(data$r,na.rm = T)
# clim = .08
clim = .31
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = r),position = "stacked",size = .3) +
  theme(legend.position = "right") + 
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 3),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))






ean <- Young_sumTb[,c(1,2:79)] %>% gather(variable, value, -AgeGroup)
df_sd <- Young_sumTb[,c(1,80:157)] %>% gather(variable, value, -AgeGroup)
df_mean$variable = rep(names(Young_Asym[,c(10:87)]),each = 6)
names(df_mean)[3] <- "MeanVal"
df_sd$variable = rep(names(Young_Asym[,c(10:87)]),each = 6)
names(df_sd)[3] <- "SDVal"
df <- merge(df_mean,df_sd)
df$Upper <- df$MeanVal + df$SDVal
df$Lower <- df$MeanVal - df$SDVal
# df <- cbind(df[c(1:34,83:116,165:198),],rep(DK.info.cortex[1:34,],3))
df2 <- cbind(df[which(str_detect(df$variable,"CT")),],rep(DK.info.cortex[1:34,],5))


ggplot(data=df2[,c(1,2,3, 4,12,14)],
       aes(x = AgeGroup,y = MeanVal, group = variable ,color = lobe)) +
  # geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  scale_color_nejm() + 
  geom_line() +
  geom_errorbar(aes(ymin = MeanVal - SDVal, ymax = MeanVal + SDVal), width = 0.2)+
  theme_minimal()

# rm(df_mean,df_sd)
ggplot(data=df2[,c(1,2,3, 4,12,14)],
       aes(x = AgeGroup,y = MeanVal, group = variable ,color = lobe)) +
  # geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  scale_color_nejm() + 
  geom_smooth(method = "lm", formula = y ~ x, se = F) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  theme_minimal()


ggplot(data=df2[,c(1,2,3,12,14)],aes(x = AgeGroup,y = MeanVal, group = variable ,color = lobe)) +
  # geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  geom_point() + geom_line() + scale_color_nejm() + 
  geom_smooth(aes(group = 1),method = "lm", formula = y ~ x) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  theme_minimal()

# ????Э??????

abcd_asym.Norm.bl <- Young_Asym[Young_Asym$AgeGroup==10,]
abcd_asym.Norm.fu2 <- Young_Asym[Young_Asym$AgeGroup==12,]

subid_bl2fu2 <- intersect(abcd_asym.Norm.bl$SubID,abcd_asym.Norm.fu2$SubID)
bl2fu2 <- which(abcd_asym.Norm.bl$SubID %in% subid_bl2fu2)
fu22bl <- which(abcd_asym.Norm.fu2$SubID %in% subid_bl2fu2)
abcd_asym.Norm.delta <- cbind(abcd_asym.Norm.bl[bl2fu2,1:9],
                              (abcd_asym.Norm.fu2[fu22bl,10:87] - abcd_asym.Norm.bl[bl2fu2,10:87])/
                                (abcd_asym.Norm.fu2$age[fu22bl] - abcd_asym.Norm.bl$age[bl2fu2]))
abcd_asym.Norm.delta$age <- abcd_asym.Norm.fu2$age[fu22bl]
abcd_asym.Norm.delta$age_gap <- (abcd_asym.Norm.fu2$age[fu22bl] - abcd_asym.Norm.bl$age[bl2fu2])
# names(abcd_asym.Norm.delta)[1] <- "SubID"

abcd_asym.Norm.delta <- abcd_asym.Norm.delta[,c(1,2,88,3:87)]


IMAGEN_Asym.bl <- Young_Asym[Young_Asym$AgeGroup==14,]
IMAGEN_Asym.fu1 <- Young_Asym[Young_Asym$AgeGroup==19,]
IMAGEN_Asym.fu2 <- Young_Asym[Young_Asym$AgeGroup==22,]

subid_bl2fu2 <- intersect(IMAGEN_Asym.bl$SubID,IMAGEN_Asym.fu2$SubID)
bl2fu2 <- which(IMAGEN_Asym.bl$SubID %in% subid_bl2fu2)
fu22bl <- which(IMAGEN_Asym.fu2$SubID %in% subid_bl2fu2)
IMAGEN_Asym.delta <- cbind(IMAGEN_Asym.bl[bl2fu2,1:9],
                           (IMAGEN_Asym.fu2[fu22bl,10:87] - IMAGEN_Asym.bl[bl2fu2,10:87])/
                             (IMAGEN_Asym.fu2$age[fu22bl] - IMAGEN_Asym.bl$age[bl2fu2]))
IMAGEN_Asym.delta$age <- IMAGEN_Asym.fu2$age[fu22bl]
IMAGEN_Asym.delta$age_gap <- (IMAGEN_Asym.fu2$age[fu22bl] - IMAGEN_Asym.bl$age[bl2fu2])
IMAGEN_Asym.delta <- IMAGEN_Asym.delta[,c(1,2,88,3:87)]


subid_bl2fu1 <- intersect(IMAGEN_Asym.bl$SubID,IMAGEN_Asym.fu1$SubID)
bl2fu1 <- which(IMAGEN_Asym.bl$SubID %in% subid_bl2fu1)
fu12bl <- which(IMAGEN_Asym.fu1$SubID %in% subid_bl2fu1)
IMAGEN_Asym.delta <- cbind(IMAGEN_Asym.bl[bl2fu1,1:9],
                           (IMAGEN_Asym.fu1[fu12bl,10:87] - IMAGEN_Asym.bl[bl2fu1,10:87])/
                             (IMAGEN_Asym.fu1$age[fu12bl] - IMAGEN_Asym.bl$age[bl2fu1]))
IMAGEN_Asym.delta$age <- IMAGEN_Asym.fu1$age[fu12bl]
IMAGEN_Asym.delta$age_gap <- (IMAGEN_Asym.fu1$age[fu12bl] - IMAGEN_Asym.bl$age[bl2fu1])
IMAGEN_Asym.delta <- IMAGEN_Asym.delta[,c(1,2,88,3:87)]
IMAGEN_Asym.delta = na.omit(IMAGEN_Asym.delta)


subid_bl2fu2 <- intersect(IMAGEN_Asym.bl$SubID,IMAGEN_Asym.fu2$SubID)
bl2fu2 <- which(IMAGEN_Asym.bl$SubID %in% subid_bl2fu2)
fu22bl <- which(IMAGEN_Asym.fu2$SubID %in% subid_bl2fu2)
IMAGEN_Asym.delta <- cbind(IMAGEN_Asym.bl[bl2fu2,1:9],
                           (IMAGEN_Asym.fu2[fu22bl,10:87] - IMAGEN_Asym.bl[bl2fu2,10:87])/
                             (IMAGEN_Asym.fu2$age[fu22bl] - IMAGEN_Asym.bl$age[bl2fu2]))
IMAGEN_Asym.delta$age <- IMAGEN_Asym.fu2$age[fu22bl]
IMAGEN_Asym.delta$age_gap <- (IMAGEN_Asym.fu2$age[fu22bl] - IMAGEN_Asym.bl$age[bl2fu2])
IMAGEN_Asym.delta <- IMAGEN_Asym.delta[,c(1,2,88,3:87)]




rMat = cor(x = as.matrix(abcd_asym.Norm.delta[,c(11:88)]),use = "pairwise")
rownames(rMat)<-gsub("_"," ",rownames(rMat))
colnames(rMat)<-gsub("_"," ",colnames(rMat))
pheatmap(rMat,
         legend = TRUE,display_numbers = FALSE,
         cluster_rows = F, cluster_cols = F,fontsize = 7,
         breaks = unique(c(seq(-0.8,0.8, length=100))),
         color = colorRampPalette(c("blue", "white", "red"))(100)
)
abcd_CTDeltaDC = rowSums(rMat[35:68,35:68])
abcd_CSADeltaDC = rowSums(rMat[1:34,1:34])
abcd_Cor = diag(rMat[1:34,35:68])

rMat = cor(x = as.matrix(IMAGEN_Asym.delta[,c(11:88)]),use = "pairwise.complete.obs")
rownames(rMat)<-gsub("_"," ",rownames(rMat))
colnames(rMat)<-gsub("_"," ",colnames(rMat))
pheatmap(rMat,
         legend = TRUE,display_numbers = FALSE,
         cluster_rows = F, cluster_cols = F,fontsize = 7,
         breaks = unique(c(seq(-0.8,0.8, length=100))),
         color = colorRampPalette(c("blue", "white", "red"))(100)
)
# 800 750
IMAGEN_Cor = diag(rMat[1:34,35:68])

IMAGEN_CTDeltaDC = rowSums(rMat[35:68,35:68])
IMAGEN_CSADeltaDC = rowSums(rMat[1:34,1:34])

data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),1),
                  r = IMAGEN_Cor,hemi = rep(c(rep("left", 34)), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), 1),
                     r = NA,hemi = rep(c(rep("left", 2)), 1))
data = rbind(data,data.NA);

clim = .8
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = r),position = "stacked",size = .3) +
  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 3),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))
## 1100 550


DC <- abcd_CSADeltaDC# abcd_CTDC#  IMAGEN_CTDeltaDC
clim1 <- min(DC);clim2 <- max(DC);if (clim1<0){clim1<- -1*clim2}else{clim1 <-0}
# clim1 <- 0;clim2 <- 2.1
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),1),DC = DC,hemi = rep(c(rep("left", 34)), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'),1),DC = NA,hemi = rep(c(rep("left", 2)), 1));data = rbind(data,data.NA)
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",mapping = aes(fill = DC),position = "stacked",size = .3) +
  theme(legend.position = "right") + scale_fill_gradientn(colours = brewer.pal(11, "Reds"),
                                                          na.value = "lightgrey",breaks = seq(clim1, clim2, 1),limits = c(clim1, clim2)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))

data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",mapping = aes(fill = DC),position = "stacked",size = .3) +
  theme(legend.position = "right") + scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:-1:1],
                                                          na.value = "lightgrey",breaks = seq(clim1, clim2, 1),limits = c(clim1, clim2)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))
# 550  400

## 3. ???? x ????ģʽ????ϵͼ

label <- c("BST","CACg","CMF","Cu","En","Fu","IP","IT","IstCg","LO","LOrF","Lg",
           "MOrF","MT","PaH","PaC","Op","Or","Tr","PerCa","PoC","PoCg","PreC","PreCu",
           "RoACg","RoMF","SF","SP","ST","SM","FPol","TPol","TrT","Ins")
  
# gsub("CSA_","",names(Young_Asym)[8:41]) 
fit.YoungAsym.Age <- fit.YoungAsym.DemoVar[fit.YoungAsym.DemoVar$demovar=="age",]
data <- data.frame(
  CSA_BL = colMeans(Young_Asym[Young_Asym$AgeGroup==10,10:43],na.rm = T),
  CT_BL = colMeans(Young_Asym[Young_Asym$AgeGroup==10,44:77],na.rm = T),
  CSA_Delta = fit.YoungAsym.Age$t[1:34],
  CT_Delta = fit.YoungAsym.Age$t[35:68],
  label = label,
  yeo_7 = DK.info.cortex$yeo_7[1:34],
  lobe = DK.info.cortex$lobe[1:34]
)

library(gridExtra)
# 600 550
p1 <- ggplot(data, aes(x = CSA_BL, y = CSA_Delta, label = label, color = lobe)) +
  geom_point() +
  geom_smooth(aes(group = 1),method = "lm", se = TRUE, color = "red") + 
  geom_text_repel() +
  xlab("CSA Asym Baseline") +
  ylab("CSA Asym Age-Effect") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_color_npg() + 
  theme_minimal()

p2 <- ggplot(data, aes(x = CT_BL, y = CT_Delta, label = label, color = lobe)) +
  geom_point() +
  geom_smooth(aes(group = 1),method = "lm", se = TRUE, color = "red") + 
  geom_text_repel() +
  xlab("CT Asym Baseline") +
  ylab("CT Asym Age-Effect") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_color_npg() + 
  theme_minimal()

p3 <- ggplot(data, aes(x = CSA_Delta, y = CT_Delta, label = label, color = lobe)) +
  geom_point() +
  geom_smooth(aes(group = 1),method = "lm", se = TRUE, color = "red") + 
  geom_text_repel() +
  xlab("CSA Asym Age-Effect") +
  ylab("CT Asym Age-Effect") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_color_npg() + 
  theme_minimal()

p4 <- ggplot(data, aes(x = CSA_BL, y = CT_BL, label = label, color = lobe)) +
  geom_point() +
  geom_smooth(aes(group = 1),method = "lm", se = TRUE, color = "red") + 
  geom_text_repel() +
  xlab("CSA Asym Baseline") +
  ylab("CT Asym Baseline") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  scale_color_npg() + 
  theme_minimal()

grid.arrange(p1 + coord_cartesian(xlim = c(-.2, .2), ylim = c(-11, 11)), 
             p2 + coord_cartesian(xlim = c(-.04, .04), ylim = c(-28, 28)), 
             p3 + coord_cartesian(xlim = c(-11, 11), ylim = c(-28, 28)), 
             p4 + coord_cartesian(xlim = c(-.2, .2), ylim = c(-.04, .04)),
             nrow = 2, ncol = 2)


data2 <- data.frame(
  Vol_BL = colMeans(Young_Asym[which(Young_Asym$AgeGroup==10),78:87]),
  Vol_Delta = fit.YoungAsym.Age$t[69:78],
  label = gsub("Vol_","",names(Young_Asym)[78:87])
)


ggplot(data2, aes(x = Vol_BL, y = Vol_Delta, label = label, color = label)) +
  geom_smooth(aes(group = 1),method = "lm", se = TRUE, color = "red") + 
  geom_point() +
  geom_text_repel() +
  scale_color_brewer(palette = "Paired") + 
  xlab("Vol Asym Baseline") +
  ylab("Vol Asym Age-Effect") + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  theme_minimal()




# ??????Ҫ?Ŀ?
library(ggplot2)
library(ggsci)
# ????ʾ?????ݿ?
data <- data.frame(
  lobe = rep(DK.info.cortex$lobe, 4),
  tval = c(fit.YoungFreesurfer.Age$t[1:136]),
  modal = c(rep("CSA",68),rep("CT",68)),
  hemi = c(rep("left",34),rep("right",34),rep("left",34),rep("right",34)),
  roi = rep(label,4)
)

subCortROI = gsub("Vol_","",fit.YoungFreesurfer.Age$brainregion[137:156])
subCortROI = gsub("Left_","",subCortROI)
subCortROI = gsub("Right_","",subCortROI)

data2 <- data.frame(
  tval = c(fit.YoungFreesurfer.Age$t[137:156]),
  hemi = c(rep("left",10),rep("right",10)),
  roi = subCortROI
)

data2 <- data2[-1*which(data2$roi=="Cerebellum_Cortex"),]

dataCSA = data[data$modal=="CSA",]
dataCT = data[data$modal=="CT",]
# tdiff = dataCSA[1:34,"tval"] - data[35:68,"tval"]
# o = dataCSA$roi[order(DK.info.cortex$lobe[1:34][34:-1:1])]
dataCSA$roi <- factor(dataCSA$roi,ordered = TRUE,levels = label[34:1])
dataCT$roi <- factor(dataCT$roi,ordered = TRUE,levels = label[34:1])
# ?????Գ?????ͼ
ggplot(dataCSA, mapping = aes(x = roi, y = tval, fill = hemi)) + 
  #facet_grid( ~ modal) + 
  geom_bar(stat = "identity",alpha = .9,position = "dodge") + 
  scale_fill_npg() + 
  theme_minimal() + coord_flip()

ggplot(dataCT, mapping = aes(x = roi, y = tval, fill = hemi)) + 
  #facet_grid( ~ modal) + 
  geom_bar(stat = "identity",alpha = .9,position = "dodge") + 
  scale_fill_npg() + 
  theme_minimal() + coord_flip()
# 270 550

data2$roi <- factor(data2$roi,ordered = TRUE,levels = subCortROI[c(10:3,1)])
ggplot(data2, mapping = aes(x = roi, y = tval, fill = hemi)) + 
  #facet_grid( ~ modal) + 
  geom_bar(stat = "identity",alpha = .9,position = "dodge") + 
  scale_fill_npg() + 
  theme_minimal() + coord_flip()
#320 220
################################ Figure. 4 ####################################
## 4. Ԫ??????Meta?????ʡ???????

id = c(45,52,65,47,53,51,60,61,37,57,49,59,36,56,43,55,64,62,
       41,58,38,54,44,46,40,50,39,66,42,48,35,63,67,68,11,18,
       31,13,19,17,26,27,3,23,15,25,2,22,9,21,30,28,7,24,4,20,
       10,12,6,16,5,32,8,14,1,29,33,34)
library(wordcloud2)

CT_tAsym = fit.YoungAsym.Age$t[35:68]
CSA_tAsym = fit.YoungAsym.Age$t[1:34]

tmap = c(CT_tAsym,CT_tAsym)
# tmap = c(CSA_tAsym,CSA_tAsym)
x = cor(Neurosynth_maps[1:68,2:124],tmap[id]);rm(dataWeight)
dataWeight <- data.frame(word = names(Neurosynth_maps[,2:124]),weight = abs(x))
dataWeight$Symbol = NA;dataWeight$Symbol[which(x>0)]<- "red";dataWeight$Symbol[which(x<0)]<- "skyblue";
dataWeight<-dataWeight[-1*which(dataWeight$weight<0.15),]
wordcloud2(dataWeight, color = dataWeight$Symbol, fontWeight = "bold",rotateRatio = 0,size = .25)

dataWeight <- data.frame(word = names(Neurosynth_maps[,which(x>0)+1]),weight = abs(x)[which(x>0)])
dataWeight$Symbol<- "red";
wordcloud2(dataWeight, color = dataWeight$Symbol, fontWeight = "bold",rotateRatio = 0,size = .15)

dataWeight <- data.frame(word = names(Neurosynth_maps[,which(x<0)+1]),weight = abs(x)[which(x<0)])
dataWeight$Symbol<- "skyblue";
wordcloud2(dataWeight, color = dataWeight$Symbol, fontWeight = "bold",rotateRatio = 0,size = .15)


tmap = c(CT_tAsym,CT_tAsym)
# tmap = c(CSA_tAsym,CSA_tAsym)
x = cor(Receptors_maps[1:68,2:19],tmap[id]);rm(dataWeight)
dataWeight <- data.frame(word = names(Receptors_maps[,2:19]),weight = abs(x))
dataWeight$Symbol = NA;dataWeight$Symbol[which(x>0)]<- "red";dataWeight$Symbol[which(x<0)]<- "skyblue";
dataWeight<-dataWeight[-1*which(dataWeight$weight<0.15),]
wordcloud2(dataWeight, color = dataWeight$Symbol, fontWeight = "bold",rotateRatio = 0,size = .25)

Receptors_names
names(Receptors_maps)

tmap = c(CT_tAsym,CT_tAsym)
# tmap = c(CSA_tAsym,CSA_tAsym)
x = cor(Atrophy_maps[1:68,2:13],tmap[id]);rm(dataWeight)
dataWeight <- data.frame(word = names(Atrophy_maps[,2:13]),weight = abs(x))
dataWeight$Symbol = NA;dataWeight$Symbol[which(x>0)]<- "red";dataWeight$Symbol[which(x<0)]<- "skyblue";
dataWeight<-dataWeight[-1*which(dataWeight$weight<0.15),]
wordcloud2(dataWeight, color = dataWeight$Symbol, fontWeight = "bold",rotateRatio = 0,size = .25)


tmap = c(CT_tAsym,CT_tAsym)
# tmap = c(CSA_tAsym,CSA_tAsym)
x = cor(parcelExpressionLhDK[1:68,2:10028],tmap[id],use = "pairwise");rm(dataWeight)
dataWeight <- data.frame(word = names(parcelExpressionLhDK[,2:10028]),weight = abs(x))
dataWeight$Symbol = NA;dataWeight$Symbol[which(x>0)]<- "red";dataWeight$Symbol[which(x<0)]<- "skyblue";
dataWeight<-dataWeight[-1*which(dataWeight$weight<0.35),];dataWeight$weight <- dataWeight$weight - 0.35
wordcloud2(dataWeight, color = dataWeight$Symbol, fontWeight = "bold",rotateRatio = 0,size = .3)
dataWeight


################################ Figure. 5 ####################################
## 5. ??Ϊ????֪??֢״?????ط?????abcd+imagen??-PCA

abcd_asym.Norm <- abcd_asym.Norm[,-1*which(str_detect(names(abcd_asym.Norm),"_vol_cdk_"))]
abcd_asym.Norm.bl <- abcd_asym.Norm.bl[,-1*which(str_detect(names(abcd_asym.Norm.bl),"_vol_cdk_"))]
abcd_asym.Norm.fu2 <- abcd_asym.Norm.fu2[,-1*which(str_detect(names(abcd_asym.Norm.fu2),"_vol_cdk_"))]
abcd_asym.Norm.delta <- abcd_asym.Norm.delta[,-1*which(str_detect(names(abcd_asym.Norm.delta),"_vol_cdk_"))]




abcd_AsymVar <- merge(abcd_covaribles[,c(1:24)],abcd_multiVar,all = F)
abcd_AsymVar$race.4level[abcd_AsymVar$race.4level!="White"&
               abcd_AsymVar$race.4level!="Black"&
               !is.na(abcd_AsymVar$race.4level)] <- "Other/Mixed"

abcd_AsymVar$race.4level <- as.character(abcd_AsymVar$race.4level)
abcd_AsymVar.bl<-abcd_AsymVar[abcd_AsymVar$eventname=="baseline_year_1_arm_1",]
abcd_AsymVar.fu2<-abcd_AsymVar[abcd_AsymVar$eventname=="2_year_follow_up_y_arm_1",]
abcd_AsymVar.fu2<-abcd_AsymVar.fu2[,c(-9,-26,-27,-31,-33, -37)]
abcd_AsymVar.fu2<-merge(abcd_AsymVar.bl[,c(1,9)],abcd_AsymVar.fu2)
abcd_AsymVar.bl$ehi_y_ss_scoreb[abcd_AsymVar.bl$ehi_y_ss_scoreb>1]<-0




abcd_AsymVar.bl <- merge(abcd_AsymVar.bl, 
                         Young_Asym[Young_Asym$AgeGroup==10,c(1,10:87)],
                         by.x = "subjectkey",by.y = "SubID")

formulas <- list()
behav_dic <- names(abcd_multiVar)[c(3:12,18:37)]
brain_dic <- names(Young_Asym)[c(10:87)]
n<-0
for (j in 1:length(behav_dic)) {
  for (i in 1:length(brain_dic)) {
    n<-n+1
    formulas[[n]] <-  as.formula(
      paste(brain_dic[i],"~1+",behav_dic[j],
            "+sex+interview_age+ehi_y_ss_scoreb+",
            "race.6level+hisp+Pubertal+smri_vol_scs_intracranialv+",
            "(1|site_id_l:rel_family_id)",sep = ""))
  }
}

num_models <- length(formulas)
fit.abcdBehav <- vector("list", length = num_models)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
fit.abcdBehav <- 
  foreach(i = 1:num_models, .combine = "rbind") %dopar% {
    formula <- formulas[[i]]
    y <- all.vars(formula)[1]
    x <- all.vars(formula)[2]
    print(paste(y,' ~ ',x, sep = ""))
    mixed_model <- lmerTest::lmer(formula,data = abcd_AsymVar.bl)
    Estimate <- coef(summary(mixed_model))[x, "Estimate"]
    tval <- coef(summary(mixed_model))[x, "t value"]
    pval <- coef(summary(mixed_model))[x, "Pr(>|t|)"]
    df <- coef(summary(mixed_model))[x, "df"]
    data.frame(x = x, y = y, Estimate = Estimate,tval = tval, pval = pval, df = df)
  }

stopCluster(cl)
# fit.abcdBehav <- unlist(fit.abcdBehav)
fit.abcdBehav$Modal <- NA
fit.abcdBehav$Modal[which(str_detect(fit.abcdBehav$y,"CT_"))] <- "CT"
fit.abcdBehav$Modal[which(str_detect(fit.abcdBehav$y,"CSA_"))] <- "CSA"
# fit.abcdBehav$Modal[which(str_detect(fit.abcdBehav$y,""))] <- "CV"
fit.abcdBehav$Modal[which(str_detect(fit.abcdBehav$y,"Vol_"))] <- "Vol"
fit.abcdBehav.Figs <- fit.abcdBehav[which(fit.abcdBehav$Modal=="CT"),]
# 2500 x 800

data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic)),
                  t = fit.abcdBehav.Figs$t,
                  stringsAsFactors = FALSE, 
                  Modal = fit.abcdBehav.Figs$x,
                  hemi = rep(c(rep("left", 34*length(behav_dic))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic)),
                     t = NA,
                     stringsAsFactors = FALSE, 
                     Modal = rep(behav_dic,2),
                     hemi = rep(c(rep("left", 2*length(behav_dic))), 1))
data = rbind(data,data.NA)
clim = 6
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_wrap( ~ Modal) + 
  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))

fit.abcdBehav.fdr <- fit.abcdBehav
for (i in 1:length(behav_dic)){
  fit.abcdBehav.fdr$pval[fit.abcdBehav.fdr$x == behav_dic[i]] <- 
    p.adjust(fit.abcdBehav.fdr[fit.abcdBehav.fdr$x == behav_dic[i],"pval"],method = "fdr")
  if (sum(p.adjust(fit.abcdBehav[fit.abcdBehav$x == behav_dic[i],"pval"],method = "fdr")<0.05)==0){
    fit.abcdBehav.fdr <- fit.abcdBehav.fdr[-1*which(fit.abcdBehav.fdr$x == behav_dic[i]),]
  }
}

fit.abcdBehav.Vol <- fit.abcdBehav.fdr[which(fit.abcdBehav.fdr$Modal=="Vol"),]
fit.abcdBehav.Vol[fit.abcdBehav.Vol$pval<0.05,]

fit.abcdBehav.Figs <- fit.abcdBehav.fdr[which(fit.abcdBehav.fdr$Modal=="CT"),]
behav_dic.fdr <- unique(fit.abcdBehav.Figs$x)
# 2500 x 800
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic.fdr)),
                  t = fit.abcdBehav.Figs$tval*(fit.abcdBehav.Figs$pval<0.05), Modal = fit.abcdBehav.Figs$x,
                  hemi = rep(c(rep("left", 34*length(behav_dic.fdr))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic.fdr)),
                     t = NA,Modal = rep(behav_dic.fdr,each = 2),
                     hemi = rep(c(rep("left", 2*length(behav_dic.fdr))), 1))
data = rbind(data,data.NA);clim = 6
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_wrap( ~ Modal) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))
## 1100 500


# 1800 * 2400

abcd_AsymEnvir.bl <- merge(abcd_covaribles[,c(1:24)],abcd_enviromentalVar,all = F)
abcd_AsymEnvir.bl <- merge(abcd_AsymEnvir.bl, abcd_asym.Norm.bl)
abcd_AsymEnvir.bl$ehi_y_ss_scoreb[abcd_AsymEnvir.bl$ehi_y_ss_scoreb>1]<-0
behav_dic <- names(abcd_enviromentalVar)[c(5:13,19:21,27:51)]
brain_dic <- names(abcd_asym.Norm)[3:ncol(abcd_asym.Norm)]
fit.abcdEnvir <- data.frame() # +smri_vol_scs_intracranialv

formulas <- list();n<-0
for (j in 1:length(behav_dic)) {for (i in 1:length(brain_dic)) {n<-n+1
formulas[[n]] <-  as.formula(paste(brain_dic[i],"~1+",behav_dic[j],
                                   "+sex+interview_age+ehi_y_ss_scoreb+",
                                   "race.6level+hisp+Pubertal+smri_vol_scs_intracranialv+",
                                   "(1|site_id_l:rel_family_id)",sep = ""))}}

num_models <- length(formulas)
fit.abcdEnvir <- vector("list", length = num_models)
cl <- makeCluster(detectCores());registerDoParallel(cl)
fit.abcdEnvir <- 
  foreach(i = 1:num_models, .combine = "rbind") %dopar% {
    formula <- formulas[[i]]
    y <- all.vars(formula)[1]
    x <- all.vars(formula)[2]
    print(paste(y,' ~ ',x, sep = ""))
    mixed_model <- lmerTest::lmer(formula,data = abcd_AsymEnvir.bl)
    Estimate <- coef(summary(mixed_model))[x, "Estimate"]
    tval <- coef(summary(mixed_model))[x, "t value"]
    pval <- coef(summary(mixed_model))[x, "Pr(>|t|)"]
    df <- coef(summary(mixed_model))[x, "df"]
    data.frame(x = x, y = y, Estimate = Estimate,tval = tval, pval = pval, df = df)
  }

stopCluster(cl)
# fit.abcdBehav <- unlist(fit.abcdBehav)
fit.abcdEnvir$Modal <- NA
fit.abcdEnvir$Modal[which(str_detect(fit.abcdEnvir$y,"_thick_"))] <- "CT"
fit.abcdEnvir$Modal[which(str_detect(fit.abcdEnvir$y,"_area_"))] <- "CSA"
fit.abcdEnvir$Modal[which(str_detect(fit.abcdEnvir$y,"_vol_cdk_"))] <- "CV"
fit.abcdEnvir$Modal[which(str_detect(fit.abcdEnvir$y,"_vol_scs_"))] <- "Vol"
fit.abcdEnvir.Figs <- fit.abcdEnvir[which(fit.abcdEnvir$Modal=="CT"),]
# 2500 x 800
fit.abcdEnvir[fit.abcdEnvir$pval<(0.05/76),]



fit.abcdEnvir.fdr <- fit.abcdEnvir
for (i in 1:length(behav_dic)){
  fit.abcdEnvir.fdr$pval[fit.abcdEnvir.fdr$x == behav_dic[i]] <- 
    p.adjust(fit.abcdEnvir.fdr[fit.abcdEnvir.fdr$x == behav_dic[i],"pval"],method = "fdr")
  if (sum(p.adjust(fit.abcdEnvir[fit.abcdEnvir$x == behav_dic[i],"pval"],method = "fdr")<0.05)==0){
    fit.abcdEnvir.fdr <- fit.abcdEnvir.fdr[-1*which(fit.abcdEnvir.fdr$x == behav_dic[i]),]
  }
}
fit.abcdEnvir.fdr <- fit.abcdEnvir.fdr[-1*which(fit.abcdEnvir.fdr$Modal=="CV"),]

fit.abcdEnvir.Vol <- fit.abcdEnvir.fdr[which(fit.abcdEnvir.fdr$Modal=="Vol"),]
fit.abcdEnvir.Vol[fit.abcdEnvir.Vol$pval<0.05,]

fit.abcdEnvir.Figs <- fit.abcdEnvir.fdr[which(fit.abcdEnvir.fdr$Modal=="CSA"),]
behav_dic.fdr <- unique(fit.abcdEnvir.Figs$x)
# 2500 x 800
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic.fdr)),
                  t = fit.abcdEnvir.Figs$tval*(fit.abcdEnvir.Figs$pval<0.05), Modal = fit.abcdEnvir.Figs$x,
                  hemi = rep(c(rep("left", 34*length(behav_dic.fdr))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic.fdr)),
                     t = NA,Modal = rep(behav_dic.fdr,2),
                     hemi = rep(c(rep("left", 2*length(behav_dic.fdr))), 1))
data = rbind(data,data.NA);clim = 6
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_wrap( ~ Modal) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))

# no environmental vars significant



abcd_AsymVar.delta <- merge(abcd_covaribles[,c(1:24)],abcd_multiVar,all = F)
abcd_AsymVar.delta <- merge(abcd_AsymVar.delta, abcd_asym.Norm.delta)
abcd_AsymVar.delta <- abcd_AsymVar.delta[,c(-1*c(6:10,24))]
abcd_AsymVar.delta <- merge(abcd_AsymVar.delta,abcd_covaribles_BL[,c(1,6:10,24)])
abcd_AsymVar.delta$ehi_y_ss_scoreb[abcd_AsymVar.delta$ehi_y_ss_scoreb>1]<-0

formulas <- list()
behav_dic <- names(abcd_multiVar)[c(3:4,7:9,11:47)]
brain_dic <- names(abcd_asym.Norm)[3:ncol(abcd_asym.Norm)]
n<-0
for (j in 1:length(behav_dic)) {for (i in 1:length(brain_dic)) {n<-n+1
formulas[[n]] <-  as.formula(paste(brain_dic[i],"~1+",behav_dic[j],"+sex+interview_age+ehi_y_ss_scoreb+",
                                   "race.6level+hisp+Pubertal+smri_vol_scs_intracranialv+",
                                   "(1|site_id_l:rel_family_id)",sep = "")) }}

cl <- makeCluster(detectCores());registerDoParallel(cl)
num_models <- length(formulas)
fit.abcdDeltaBehav <- vector("list", length = num_models)
fit.abcdDeltaBehav <- 
  foreach(i = 1:num_models, .combine = "rbind",.errorhandling = "remove") %dopar% {
    formula <- formulas[[i]]
    y <- all.vars(formula)[1]
    x <- all.vars(formula)[2]
    print(paste(y,' ~ ',x, sep = ""))
    mixed_model <- lmerTest::lmer(formula,data = abcd_AsymVar.delta)
    Estimate <- coef(summary(mixed_model))[x, "Estimate"]
    tval <- coef(summary(mixed_model))[x, "t value"]
    pval <- coef(summary(mixed_model))[x, "Pr(>|t|)"]
    df <- coef(summary(mixed_model))[x, "df"]
    data.frame(x = x, y = y, Estimate = Estimate,tval = tval, pval = pval, df = df)
  }

# stopCluster(cl)
# fit.abcdBehav <- unlist(fit.abcdBehav)
fit.abcdDeltaBehav$Modal <- NA
fit.abcdDeltaBehav$Modal[which(str_detect(fit.abcdDeltaBehav$y,"_thick_"))] <- "CT"
fit.abcdDeltaBehav$Modal[which(str_detect(fit.abcdDeltaBehav$y,"_area_"))] <- "CSA"
fit.abcdDeltaBehav$Modal[which(str_detect(fit.abcdDeltaBehav$y,"_vol_cdk_"))] <- "CV"
fit.abcdDeltaBehav$Modal[which(str_detect(fit.abcdDeltaBehav$y,"_vol_scs_"))] <- "Vol"

# 2500 x 800
fit.abcdDeltaBehav.fdr <- fit.abcdDeltaBehav
for (i in 1:length(behav_dic)){
  fit.abcdDeltaBehav.fdr$pval[fit.abcdDeltaBehav.fdr$x == behav_dic[i]] <- 
    p.adjust(fit.abcdDeltaBehav.fdr[fit.abcdDeltaBehav.fdr$x == behav_dic[i],"pval"],method = "fdr")
  if (sum(p.adjust(fit.abcdDeltaBehav[fit.abcdDeltaBehav$x == behav_dic[i],"pval"],method = "fdr")<0.05)==0){
    fit.abcdDeltaBehav.fdr <- fit.abcdDeltaBehav.fdr[-1*which(fit.abcdDeltaBehav.fdr$x == behav_dic[i]),]
  }
}
fit.abcdDeltaBehav.fdr <- fit.abcdDeltaBehav.fdr[-1*which(fit.abcdDeltaBehav.fdr$Modal=="CV"),]

fit.abcdDeltaBehav.Vol <- fit.abcdDeltaBehav.fdr[which(fit.abcdDeltaBehav.fdr$Modal=="Vol"),]
fit.abcdDeltaBehav.Vol[fit.abcdDeltaBehav.Vol$pval<0.05,]

fit.abcdDeltaBehav.Figs <- fit.abcdDeltaBehav.fdr[which(fit.abcdDeltaBehav.fdr$Modal=="CSA"),]
behav_dic.fdr <- unique(fit.abcdDeltaBehav.Figs$x)
# 2500 x 800
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic.fdr)),
                  t = fit.abcdDeltaBehav.Figs$tval*(fit.abcdDeltaBehav.Figs$pval<0.05), Modal = fit.abcdDeltaBehav.Figs$x,
                  hemi = rep(c(rep("left", 34*length(behav_dic.fdr))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic.fdr)),
                     t = NA,Modal = rep(behav_dic.fdr,2),
                     hemi = rep(c(rep("left", 2*length(behav_dic.fdr))), 1))
data = rbind(data,data.NA);clim = 6
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_wrap( ~ Modal) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))














abcd_AsymEnvir.delta <- merge(abcd_covaribles[,c(1:24)],abcd_enviromentalVar,all = F)
abcd_AsymEnvir.delta <- merge(abcd_AsymEnvir.delta, abcd_asym.Norm.delta)
abcd_AsymEnvir.delta <- abcd_AsymEnvir.delta[,c(-1*c(6:10,24))]
abcd_AsymEnvir.delta <- merge(abcd_AsymEnvir.delta,abcd_covaribles_BL[,c(1,6:10,24)])
abcd_AsymEnvir.delta$ehi_y_ss_scoreb[abcd_AsymEnvir.delta$ehi_y_ss_scoreb>1]<-0

behav_dic <- names(abcd_enviromentalVar)[c(5:13,19:21,27:51)]
brain_dic <- names(abcd_asym.Norm.delta)[3:ncol(abcd_asym.Norm)]
fit.abcdDeltaEnvir <- data.frame() # +smri_vol_scs_intracranialv
formulas <- list();n<-0
for (j in 1:length(behav_dic)) {for (i in 1:length(brain_dic)) {n<-n+1
formulas[[n]] <-  as.formula(paste(brain_dic[i],"~1+",behav_dic[j],"+sex+interview_age+ehi_y_ss_scoreb+",
                                   "race.6level+hisp+Pubertal+smri_vol_scs_intracranialv+",
                                   "(1|site_id_l:rel_family_id)",sep = ""))}}

num_models <- length(formulas)
fit.abcdDeltaEnvir <- vector("list", length = num_models)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
fit.abcdDeltaEnvir <- 
  foreach(i = 1:num_models, .combine = "rbind",.errorhandling = "pass") %dopar% {
    formula <- formulas[[i]]
    y <- all.vars(formula)[1]
    x <- all.vars(formula)[2]
    print(paste(y,' ~ ',x, sep = ""))
    mixed_model <- lmerTest::lmer(formula,data = abcd_AsymEnvir.delta)
    Estimate <- coef(summary(mixed_model))[x, "Estimate"]
    tval <- coef(summary(mixed_model))[x, "t value"]
    pval <- coef(summary(mixed_model))[x, "Pr(>|t|)"]
    df <- coef(summary(mixed_model))[x, "df"]
    data.frame(x = x, y = y, Estimate = Estimate,tval = tval, pval = pval, df = df)
  }

stopCluster(cl)
# fit.abcdBehav <- unlist(fit.abcdBehav)
fit.abcdDeltaEnvir$Modal <- NA
fit.abcdDeltaEnvir$Modal[which(str_detect(fit.abcdDeltaEnvir$y,"_thick_"))] <- "CT"
fit.abcdDeltaEnvir$Modal[which(str_detect(fit.abcdDeltaEnvir$y,"_area_"))] <- "CSA"
fit.abcdDeltaEnvir$Modal[which(str_detect(fit.abcdDeltaEnvir$y,"_vol_cdk_"))] <- "CV"
fit.abcdDeltaEnvir$Modal[which(str_detect(fit.abcdDeltaEnvir$y,"_vol_scs_"))] <- "Vol"
fit.abcdDeltaEnvir.Figs <- fit.abcdDeltaEnvir[which(fit.abcdDeltaEnvir$Modal=="CT"),]
# 2500 x 800


IMAGEN_cantab <- merge(IMAGEN_cantab_BL,IMAGEN_cantab_FU2,all = TRUE)
IMAGEN_cantab <- merge(IMAGEN_cantab,IMAGEN_cantab_FU3,all = TRUE)


IMAGEN_dawba_BL$AgeGroup = 14
IMAGEN_dawba_FU1$AgeGroup = 17
IMAGEN_dawba_FU2$AgeGroup = 19
IMAGEN_dawba_FU3$AgeGroup = 22  

IMAGEN_dawba <- merge(IMAGEN_dawba_BL,IMAGEN_dawba_FU1,all = TRUE)
IMAGEN_dawba <- merge(IMAGEN_dawba,IMAGEN_dawba_FU2,all = TRUE)
IMAGEN_dawba <- merge(IMAGEN_dawba,IMAGEN_dawba_FU3,all = TRUE)

IMAGEN_LEQ <- merge(IMAGEN_LEQ_BL,IMAGEN_LEQ_FU1,all = TRUE)
IMAGEN_LEQ <- merge(IMAGEN_LEQ,IMAGEN_LEQ_FU2,all = TRUE)
IMAGEN_LEQ <- merge(IMAGEN_LEQ,IMAGEN_LEQ_FU3,all = TRUE)


IMAGEN_NI_DATA_BL$AgeGroup = 14
IMAGEN_NI_DATA_FU2$AgeGroup = 19
IMAGEN_NI_DATA_FU3$AgeGroup = 22  

IMAGEN_NI_DATA <- merge(IMAGEN_NI_DATA_BL,IMAGEN_NI_DATA_FU2,all = TRUE)
IMAGEN_NI_DATA <- merge(IMAGEN_NI_DATA,IMAGEN_NI_DATA_FU3,all = TRUE)

IMAGEN_PBQ <- merge(IMAGEN_PBQ_BL,IMAGEN_PBQ_FU1,all = TRUE)

IMAGEN_DAWBA <- IMAGEN_dawba[,c("PSC2","AgeGroup","age","gender",
                                SDQNames[4:length(SDQNames)],
                                DAWBANames[4:length(DAWBANames)])]
# IMAGEN_CANTAB <- IMAGEN_cantab

NINames <- c("WISCIV_DigitSpan_backward","WISCIV_DigitSpan_forward","WISCIV_DigitSpan_longest_backward",
             "WISCIV_DigitSpan_longest_forward","WISCIV_MatrixReasoning","WISCIV_Similarities","WISCIV_Vocabulary")

IMAGEN_WISC <- IMAGEN_NI_DATA[,c("User code","AgeGroup",NINames)]
IMAGEN_WISC[IMAGEN_WISC==-1]<-NA
names(IMAGEN_WISC)[1] <- "PSC2"
IMAGEN_WISC$PSC2 <- gsub("-I","",IMAGEN_WISC$PSC2)

IMAGEN_multiVar <- merge(IMAGEN_DAWBA,IMAGEN_WISC)
IMAGEN_multiVar <- IMAGEN_multiVar[, colSums(is.na(IMAGEN_multiVar)) != nrow(IMAGEN_multiVar)]
IMAGEN_multiVar <- IMAGEN_multiVar[, apply(IMAGEN_multiVar, 2, function(x) sd(x, na.rm = TRUE) != 0)]
# "ni_height","ni_mass"

# CANTAB:PRT	Working Memory
# CANTAB:SRMT	Spatial Memory
# CANTAB:AGNG	Emotional bias
# CANTAB:CGT	Risk - taking, decision - making
# CANTAB:RVP	Attention
# Morphed Faces Task (Pollak & Kistler, 2002)	Face recognition
# Emotional Dot - Probe (MacLeod et al., 1986)	Attentional bias to emotions
# Passive Avoidance Learning Paradigm (PALP, Arnett & Newman, 2000)	Behavioural inhibition
# WISC:Vocabulary Similarities	Verbal IQ
# WISC:Block Design, Matrix Reasoning, Digit Span	Non - verbal IQ

fit.IMAGENBehav <- data.frame()
names(IMAGEN_multiVar)[1] <- "SubID"
IMAGEN_multiVar_BL <- IMAGEN_multiVar[which(IMAGEN_multiVar$AgeGroup=='14'),]
IMAGEN_multiVar_BL <- IMAGEN_multiVar_BL[, colSums(is.na(IMAGEN_multiVar_BL))!= nrow(IMAGEN_multiVar_BL)]
IMAGEN_multiVar_BL <- IMAGEN_multiVar_BL[, apply(IMAGEN_multiVar_BL, 2, function(x) sd(x, na.rm = TRUE) != 0)]

IMAGEN_AsymBehav_BL = merge(
  IMAGEN_multiVar_BL[,c(1,2,4:ncol(IMAGEN_multiVar_BL))],
  Young_Asym[Young_Asym$AgeGroup==14,])

IMAGEN_AsymBehav_BL$Vol_EstimatedTotalIntraCranialVol = 
  IMAGEN_AsymBehav_BL$Vol_EstimatedTotalIntraCranialVol/10000
brain_dic <- names(Young_Asym)[10:87]
behav_dic <- names(IMAGEN_multiVar_BL)[4:ncol(IMAGEN_multiVar_BL)]



n = 0
for (j in 1:length(behav_dic)) {
  for (i in 1:length(brain_dic)) {
    n = n + 1
    f = as.formula(
      paste(brain_dic[i],
            " ~ age + gender + Handedness + Race +
            Vol_EstimatedTotalIntraCranialVol + 
            (1|ImagingCentreID) + ",behav_dic[j],sep = "")
    )
    res.lmer = lmerTest::lmer(f, IMAGEN_AsymBehav_BL)
    sum.lmer = summary(res.lmer)
    fit.IMAGENBehav[n, c("x")] = behav_dic[j]
    fit.IMAGENBehav[n, c("y")] = brain_dic[i]
    fit.IMAGENBehav[n, c("estimate","df", "tval", "pval")] =
      sum.lmer$coefficients[behav_dic[j], 
                            c("Estimate","df", "t value", "Pr(>|t|)")]
    print(paste(brain_dic[i],' ~ ',behav_dic[j], sep = ""))
  }
}

rm(list = c("res.lmer","sum.lmer"))
fit.IMAGENBehav$Modal <- NA
fit.IMAGENBehav$Modal[which(str_detect(fit.IMAGENBehav$y,"CT_"))] <- "CT"
fit.IMAGENBehav$Modal[which(str_detect(fit.IMAGENBehav$y,"CSA_"))] <- "CSA"
fit.IMAGENBehav$Modal[which(str_detect(fit.IMAGENBehav$y,"Vol_"))] <- "Vol"
fit.IMAGENBehav.Figs <- fit.IMAGENBehav[which(fit.IMAGENBehav$Modal=="CT"),]
# 2500 x 800 *(fit.IMAGENBehav.Figs$p<0.05/34)
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic)),
                  t = fit.IMAGENBehav.Figs$tval,behav = fit.IMAGENBehav.Figs$x,
                  hemi = rep(c(rep("left", 34*length(behav_dic))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic)),
                     t = NA,behav = rep(behav_dic,each = 2),
                     hemi = rep(c(rep("left", 2*length(behav_dic))), 1))
data = rbind(data,data.NA);clim = 4
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .5) +
  facet_wrap( ~ behav) + 
  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))


fit.IMAGENBehav.fdr <- fit.IMAGENBehav
for (i in 1:length(behav_dic)){
  fit.IMAGENBehav.fdr$pval[fit.IMAGENBehav.fdr$x == behav_dic[i]] <- 
    p.adjust(fit.IMAGENBehav.fdr[fit.IMAGENBehav.fdr$x == behav_dic[i],"pval"],method = "fdr")
  if (sum(p.adjust(fit.IMAGENBehav[fit.IMAGENBehav$x == behav_dic[i],"pval"],method = "fdr")<0.05)==0){
    fit.IMAGENBehav.fdr <- fit.IMAGENBehav.fdr[-1*which(fit.IMAGENBehav.fdr$x == behav_dic[i]),]
  }
}

fit.IMAGENBehav.Vol <- fit.IMAGENBehav.fdr[which(fit.IMAGENBehav.fdr$Modal=="Vol"),]
fit.IMAGENBehav.Vol[fit.IMAGENBehav.Vol$pval<0.05,]

fit.IMAGENBehav.Figs <- fit.IMAGENBehav.fdr[which(fit.IMAGENBehav.fdr$Modal=="CT"),]
behav_dic.fdr <- unique(fit.IMAGENBehav.Figs$x)
# 2500 x 800
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic.fdr)),
                  t = fit.IMAGENBehav.Figs$tval*(fit.IMAGENBehav.Figs$pval<0.05), Modal = fit.IMAGENBehav.Figs$x,
                  hemi = rep(c(rep("left", 34*length(behav_dic.fdr))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic.fdr)),
                     t = NA,Modal = rep(behav_dic.fdr,2),
                     hemi = rep(c(rep("left", 2*length(behav_dic.fdr))), 1))
data = rbind(data,data.NA);clim = 6
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_wrap( ~ Modal) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))
# 850 350



IMAGEN_Asym_BL <- IMAGEN_Asym.Match[which(IMAGEN_Asym.Match$AgeGroup=="14"),]
IMAGEN_Asym_FU2 <- IMAGEN_Asym.Match[which(IMAGEN_Asym.Match$AgeGroup=="19"),]
IMAGEN_Asym_FU3 <- IMAGEN_Asym.Match[which(IMAGEN_Asym.Match$AgeGroup=="22"),]
SubID.bl2fu2 <- intersect(IMAGEN_Asym_BL$SubID,IMAGEN_Asym_FU2$SubID)
bl2fu2 <- which(IMAGEN_Asym_BL$SubID %in% SubID.bl2fu2)
fu22bl <- which(IMAGEN_Asym_FU2$SubID %in% SubID.bl2fu2)
IMAGEN_Asym.bl2fu2delta = 
  (IMAGEN_Asym_FU2[fu22bl,7:84] - IMAGEN_Asym_BL[bl2fu2,7:84])/
  (IMAGEN_Asym_FU2[fu22bl,2] - IMAGEN_Asym_BL[bl2fu2,2])
IMAGEN_Asym.bl2fu2delta <- cbind(IMAGEN_Asym_BL[bl2fu2,c(1:6)],IMAGEN_Asym.bl2fu2delta)
IMAGEN_Asym.bl2fu2delta$ageGap <- IMAGEN_Asym_FU2[fu22bl,2] - IMAGEN_Asym_BL[bl2fu2,2]
IMAGEN_Asym.bl2fu2delta <- IMAGEN_Asym.bl2fu2delta[-1*which(IMAGEN_Asym.bl2fu2delta$ageGap==0),]
IMAGEN_Asym.bl2fu2delta <- IMAGEN_Asym.bl2fu2delta[-1*which(is.na(IMAGEN_Asym.bl2fu2delta$ageGap)),]


fit.IMAGENDeltaBehav <- data.frame()
IMAGEN_multiVar_FU2 <- IMAGEN_multiVar[which(IMAGEN_multiVar$AgeGroup=='19'),]
IMAGEN_multiVar_FU2 <- IMAGEN_multiVar_FU2[, colSums(is.na(IMAGEN_multiVar_FU2)) != nrow(IMAGEN_multiVar_FU2)]
IMAGEN_multiVar_FU2 <- IMAGEN_multiVar_FU2[, apply(IMAGEN_multiVar_FU2, 2, function(x) sd(x, na.rm = TRUE) != 0)]
IMAGEN_multiVar_FU3 <- IMAGEN_multiVar[which(IMAGEN_multiVar$AgeGroup=='22'),]
IMAGEN_multiVar_FU3 <- IMAGEN_multiVar_FU3[, colSums(is.na(IMAGEN_multiVar_FU3)) != nrow(IMAGEN_multiVar_FU3)]
IMAGEN_multiVar_FU3 <- IMAGEN_multiVar_FU3[, apply(IMAGEN_multiVar_FU3, 2, function(x) sd(x, na.rm = TRUE) != 0)]

IMAGEN_AsymDeltaBehav_FU2 = merge(IMAGEN_multiVar_FU2[,c(1,2,4:ncol(IMAGEN_multiVar_FU2))],
                                  IMAGEN_Asym.bl2fu2delta[,c(1,3:ncol(IMAGEN_Asym.bl2fu2delta))])
IMAGEN_AsymDeltaBehav_FU2$Vol_EstimatedTotalIntraCranialVol = 
  IMAGEN_AsymDeltaBehav_FU2$Vol_EstimatedTotalIntraCranialVol/10000
brain_dic <- names(IMAGEN_Asym.Match)[7:84]
behav_dic <- names(IMAGEN_multiVar_FU2)[4:ncol(IMAGEN_multiVar_FU2)]
n = 0
for (j in 1:length(behav_dic)) {
  for (i in 1:length(brain_dic)) {
    n = n + 1
    f = as.formula(
      paste(brain_dic[i],
            " ~ age + gender + Handedness + Vol_EstimatedTotalIntraCranialVol + ",
            "(1|ImagingCentreID) + ",behav_dic[j],sep = "")
    )
    res.lmer = lmerTest::lmer(f, IMAGEN_AsymDeltaBehav_FU2)
    sum.lmer = summary(res.lmer)
    fit.IMAGENDeltaBehav[n, c("x")] = behav_dic[j]
    fit.IMAGENDeltaBehav[n, c("y")] = brain_dic[i]
    fit.IMAGENDeltaBehav[n, c("estimate","df", "tval", "pval")] =
      sum.lmer$coefficients[behav_dic[j], 
                            c("Estimate","df", "t value", "Pr(>|t|)")]
    print(paste(brain_dic[i],' ~ ',behav_dic[j], sep = ""))
    
  }
}

rm(list = c("res.lmer","sum.lmer"))
fit.IMAGENDeltaBehav$Modal <- NA
fit.IMAGENDeltaBehav$Modal[which(str_detect(fit.IMAGENDeltaBehav$brainregion,"CT_"))] <- "CT"
fit.IMAGENDeltaBehav$Modal[which(str_detect(fit.IMAGENDeltaBehav$brainregion,"CSA_"))] <- "CSA"
fit.IMAGENDeltaBehav$Modal[which(str_detect(fit.IMAGENDeltaBehav$brainregion,"Vol_"))] <- "Vol"
fit.IMAGENDeltaBehav.Figs <- fit.IMAGENDeltaBehav[which(fit.IMAGENDeltaBehav$Modal=="CSA"),]
# 2500 x 800 *(fit.IMAGENDeltaBehav.Figs$p<0.05/34)
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic)),
                  t = fit.IMAGENDeltaBehav.Figs$t,behav = fit.IMAGENDeltaBehav.Figs$behaviorvar,
                  hemi = rep(c(rep("left", 34*length(behav_dic))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic)),
                     t = NA,behav = rep(behav_dic,each = 2),
                     hemi = rep(c(rep("left", 2*length(behav_dic))), 1))
data = rbind(data,data.NA);data$t <- as.numeric(data$t);clim = 4
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .5) +
  facet_wrap( ~ behav) + 
  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))


fit.IMAGENDeltaBehav.fdr <- fit.IMAGENDeltaBehav
for (i in 1:length(behav_dic)){
  fit.IMAGENDeltaBehav.fdr$p[fit.IMAGENDeltaBehav.fdr$behaviorvar == behav_dic[i]] <- 
    p.adjust(fit.IMAGENDeltaBehav.fdr[fit.IMAGENDeltaBehav.fdr$behaviorvar == behav_dic[i],"p"],method = "fdr")
  if (sum(p.adjust(fit.IMAGENDeltaBehav[fit.IMAGENDeltaBehav$behaviorvar == behav_dic[i],"p"],method = "fdr")<0.05)==0){
    fit.IMAGENDeltaBehav.fdr <- fit.IMAGENDeltaBehav.fdr[-1*which(fit.IMAGENDeltaBehav.fdr$behaviorvar == behav_dic[i]),]
  }
}

fit.IMAGENDeltaBehav.Vol <- fit.IMAGENDeltaBehav.fdr[which(fit.IMAGENDeltaBehav.fdr$Modal=="Vol"),]
fit.IMAGENDeltaBehav.Vol[fit.IMAGENDeltaBehav.Vol$p<0.05,]

fit.IMAGENDeltaBehav.Figs <- fit.IMAGENDeltaBehav.fdr[which(fit.IMAGENDeltaBehav.fdr$Modal=="CSA"),]
behav_dic.fdr <- unique(fit.IMAGENDeltaBehav.Figs$behaviorvar)
# 2500 x 800
data = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),length(behav_dic.fdr)),
                  t = fit.IMAGENDeltaBehav.Figs$t*(fit.IMAGENDeltaBehav.Figs$p<0.05), Modal = fit.IMAGENDeltaBehav.Figs$behaviorvar,
                  hemi = rep(c(rep("left", 34*length(behav_dic.fdr))), 1))
data.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(behav_dic.fdr)),
                     t = NA,Modal = rep(behav_dic.fdr,2),
                     hemi = rep(c(rep("left", 2*length(behav_dic.fdr))), 1))
data = rbind(data,data.NA);clim = 6
data %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = t),position = "stacked",size = .3) +
  facet_wrap( ~ Modal) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 1),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))






## PLS
library(pls)
library(permimp)
# ׼??????
X <- matrix(...)  # ??????�� X ?????ݣ???????һ???????????ݿ???ÿ?д???һ????��
Y <- matrix(...)  # ??Ӧ??�� Y ?????ݣ???????һ???????????ݿ???ÿ?д???һ????��
C <- matrix(...)  # Э??�� C ?????ݣ???????һ???????????ݿ???ÿ?д???һ????��
# ִ?? PLSR ????
pls_model <- plsr(X, Y, scale = TRUE, x_scores = X, y_scores = Y, covariates = C)
# ִ???û?????
perm_result <- permimp(pls_model, nperm = 1000)  # ?????û???????????ʾ??Ϊ1000??
# ?鿴?û?????????
summary(perm_result)
scores(perm_result)
# ??ȡ??��??Ҫ??
loadings <- loadings(pls_model)
# ???չ??׶Ƚ???????
sorted_loadings <- sort(abs(loadings), decreasing = TRUE)
# ????Ȩ?????ߵı?��
top_variables <- names(sorted_loadings)[1:5]  # ѡ??ǰ5????��?????Ը?????Ҫ???е???
print(top_variables)




################################ Figure. 6 ####################################
## 6. ???????ط?????abcd?????б���-PCA

################################ Figure. 7 ####################################
## 7. GWAS????????????ƫ/??ƫ??????????/???٣? ?????б���?໥??��??????
GWAS_Path = "E:/Acdamic/Article_ABCD_Asymmetry/ABCD_GWAS/FUMA_gene2func114166/"
# 161 165 166 167
GS <- read.cVol(paste(GWAS_Path,"GS.txt",sep = ""),sep = "\t")
Pthres = 1e-10
# GSFigs <- rbind(GS[GS$adjP<Pthres&GS$Category=="GO_bp",1:6],
#                 GS[GS$adjP<Pthres&GS$Category=="GO_cc",1:6],
#                 GS[GS$adjP<Pthres&GS$Category=="GO_mf",1:6])


GSFigs <- rbind(GS[GS$Category=="GO_bp",1:6],
                GS[GS$Category=="GO_cc",1:6],
                GS[GS$Category=="GO_mf",1:6])
GSFigs <- GSFigs[order(GSFigs$adjP)[1:30],]
GSFigs <- GSFigs[-1*which(GSFigs$adjP>Pthres),]

# GS[GS$adjP<Pthres&GS$Category=="KEGG",1:6]
# GS[GS$adjP<Pthres&GS$Category=="GWAScatalog",1:6]
library(wordcloud2)
rm(dataWeight)
dataWeight <- data.frame(word = GSFigs$GeneSet,weight = -1*log(GSFigs$adjP)/200,category = GSFigs$Category)
dataWeight$col = NA;
dataWeight$col[which(dataWeight$category=="GO_bp")]<- pal_npg("nrc")(3)[1];
dataWeight$col[which(dataWeight$category=="GO_cc")]<- pal_npg("nrc")(3)[2];
dataWeight$col[which(dataWeight$category=="GO_mf")]<- pal_npg("nrc")(3)[3];
wordcloud2(dataWeight, color = dataWeight$col, rotateRatio = 0,size = .1)

DEG <- read.cVol(paste(GWAS_Path,"gtex_v8_ts_general_DEG.txt",sep = ""),sep = "\t")
d <- -1*log(DEG$adjP)
d[which(DEG$adjP>0.05)] <- NA
DEG.Mat <- matrix(data = d,nrow = 30,ncol = 3)
row.names(DEG.Mat) <- DEG$GeneSet[1:30]
colnames(DEG.Mat) <- unique(DEG$Category)

pheatmap(t(DEG.Mat[order(DEG.Mat[,3],decreasing = TRUE),]),cluster_cols = FALSE,cluster_rows = FALSE)



DEG <- read.cVol(paste(GWAS_Path,"gtex_v8_ts_DEG.txt",sep = ""),sep = "\t")
d <- -1*log(DEG$adjP)
d[which(DEG$adjP>0.05)] <- NA
DEG.Mat <- matrix(data = d,nrow = 54,ncol = 3)
row.names(DEG.Mat) <- DEG$GeneSet[1:54]
colnames(DEG.Mat) <- unique(DEG$Category)

pheatmap(t(DEG.Mat[order(DEG.Mat[,3],decreasing = TRUE),]),cluster_cols = FALSE,cluster_rows = FALSE)

################################ Figure. 8 ####################################
## 8. ???۽??ͣ????ܴ??ڵ???ǰ???????ͺ?????



#################################### SM #######################################
## ??¼
############################## SM Figure. 1 ###################################
## 1. Ƥ????????????ƫ?໯????????
############################## SM Figure. 2 ###################################
## 2. ?˿?ѧ???ع?��???????Ա????ഺ?ڡ??Ա????佻??ЧӦ?????֣?
############################## SM Figure. 3 ###################################
## 3. ?????ݿ?????????????????
############################## SM Figure. 4 ###################################
## 4. imagen????????
############################## SM Figure. 5 ###################################
## 5. ??????ЧӦ????
library(gamm4)
# library(ade4)
# ѭ??????ÿһ?е?ģ??
models <- list()
for (i in 8:85) {
  # ??????ʽ
  formula <- paste(names(Young_Asym)[i], 
                   "~ s(age) + gender + Handedness + Vol_EstimatedTotalIntraCranialVol")
  model <- gamm4(as.formula(formula), data = Young_Asym,random= ~ (1|ImagingCentreID))
  models[[i-7]] <- model
}
# ????Ԥ??????
plot(models[[10]]$gam)
# ѭ??Ԥ???ͻ?ͼ
for (i in 1:78) {
  p <- ggplot(aes(x = models[[i]]$gam$model$age, y = models[[i]]$gam$y)) +
    geom_line() +
    xlab("Age") +
    ylab(paste("Outcome", i)) +
    ggtitle(paste("Outcome", i, "by Age"))
  print(p)
}


############################## SM Figure. 6 ###################################
## 6. ???????ط?????????Ԥ??δ��???????켣Ԥ??δ��??


fit.YoungFreesurfer.Age <- data.frame()
for (i in 10:167){
  f = as.formula(paste(names(Young_Freesurfer)[i],
                       " ~ age + gender + Handedness + ",
                       "Vol_EstimatedTotalIntraCranialVol + ",
                       "(1|Dataset:ImagingCentreID) + (1|SubID)",sep = ""))
  res.lmer = lmerTest::lmer(f,Young_Freesurfer)
  sum.lmer = summary(res.lmer)
  fit.YoungFreesurfer.Age[i-9,c("y")] = names(Young_Freesurfer)[i]
  fit.YoungFreesurfer.Age[i-9,c("estimate","df", "t", "p")] = 
    sum.lmer$coefficients["age", c("Estimate","df", "t value", "Pr(>|t|)")]
  print(paste(names(Young_Freesurfer)[i],' ~ Age', sep = ""))
}
rm(list = c("res.lmer","sum.lmer"))

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
dk_id.NA <- formatC(c(1,5), width = 4, flag = '0')

data = data.frame(roi = rep(dk_id, 1),t = fit.YoungFreesurfer.Age$t[1:(68*2)],
                  stringsAsFactors = FALSE, Modal = c(rep("CSA", 68),rep("CT", 68)),
                  hemi = rep(c(rep("left", 34),rep("right", 34)), 2));
clim = 100
data[data$Modal=="CT",] %>%
  group_by(Modal) %>% 
  ggseg(atlas = dk,colour = "black",# hemisphere = "both",
        mapping = aes(fill = t),position = "stacked",size = .5) +
  facet_wrap( ~ Modal, ncol = 1) + 
  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey", breaks = seq(-1*clim, 1, clim),limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))

# 660 450
############################## SM Figure. 7 ###################################
## 7. GWAS Manhattan ͼ

## GCTA
# ?????ļ?·??
genotype_file <- "genotype.bed"  # ?滻Ϊʵ?ʵĻ??????????ļ?·????????
phenotype_file <- "phenotype.txt"  # ?滻Ϊʵ?ʵı????????ļ?·????????
# ????GCTA????
gcta_cmd <- paste("gcta64 --bfile", genotype_file, "--pheno", phenotype_file,
                  "--make-grm --out output")
# ????GCTA????
system(gcta_cmd)




library(tableone)
CreateTableOne(
               data = Young_Asym[,c(2:3,5:7,86,87)],
               factorVars = c("AgeGroup"))

print(Young_grouped)

res <- Young_Asym %>%
  group_by(AgeGroup,gender) %>%           # ??Age??Gender?н??з???
  summarise(
    Count = n()
  )


