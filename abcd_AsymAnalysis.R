# ABCD Plots
# bl 2y 4y & Age Effect & Age Effect (Raw Measures)
# Correlation between Patterns
# Set a different CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# List of packages to install and load
packages <- c(
  "cluster", "factoextra", "psych", "ggradar", "pheatmap", "ggplot2",
  "ggiraph", "ggiraphExtra", "RColorBrewer", "tableone",
  "glmnet", "rgl", "ggseg", "dplyr", "gamm4", "viridis", "reshape2",
  "Rmisc", "ggpubr", "Rcpp", "ggridges", "stringr", "ggfun","gridExtra",
  "tidyr", "doParallel", "gghalves", "utils", "readr", "parallel",
  "lmerTest", "foreach", "ggsci", "longCombat", "circlize","qvalue",
  "ggrepel","VennDiagram","venn","snpsettest","wordcloud2","matrixStats",
  "UpSetR","ggalt"
)

# Install and load packages
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


# 读取数据
data <- read_csv("E:/Acdamic/Article_ABCD_Asymmetry/abcd_AsymPheno.csv",na = "NaN")
data$eventname <- recode(data$eventname, 
                         'baseline_year_1_arm_1' = 'bl', 
                         '2_year_follow_up_y_arm_1' = '2y', 
                         '4_year_follow_up_y_arm_1' = '4y')
data$rel_family_id <- as.character(data$rel_family_id)
data$rel_birth_id <- as.character(data$rel_birth_id)
data$school_id <- as.character(data$school_id)
data$district_id <- as.character(data$district_id)
data$smri_vol_scs_intracranialv <- data$smri_vol_scs_intracranialv/10000
ehi_y_ss_scoreb <- data$ehi_y_ss_scoreb
data$ehi_y_ss_scoreb[ehi_y_ss_scoreb==1] <- 3
data$ehi_y_ss_scoreb[ehi_y_ss_scoreb==2] <- 1
data$ehi_y_ss_scoreb[ehi_y_ss_scoreb==3] <- 2


data$eventname <- factor(data$eventname, levels = c("bl", "2y", "4y"))

ggplot(data, aes(x = eventname, y = smri_area_cdk_lobfr, group = src_subject_id)) +
  geom_line(alpha = .1) +
  geom_point() +
  labs(title = "Individual Trajectories Over Time",
       x = "Time",
       y = "Value",
       color = "Subject ID") +
  theme_minimal() +
  theme(legend.position = "none")  # 隐藏图例（如果被试太多，可以保留图例）


data_long <- melt(data[,c(1,2,16:95)], id.vars = c("src_subject_id", "eventname"))

# 计算每个时间点和每个变量的均值和标准误
summary_df <- aggregate(value ~ eventname + variable, data_long, function(x) {
  mean_value <- mean(x, na.rm = TRUE)
  se_value <- sd(x, na.rm = TRUE) / sqrt(length(x))
  c(mean = mean_value, se = se_value)
})

# 分离均值和标准误
summary_df <- data.frame(
  time = summary_df$eventname,
  variable = summary_df$variable,
  mean_value = summary_df$value[, "mean"],
  se_value = summary_df$value[, "se"]
)

# 检查 summary_df 结构
head(summary_df)

# 绘制带误差条的均值图
ggplot(summary_df, aes(x = time, y = mean_value, color = variable, group = variable)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  labs(title = "Mean Value with Error Bars at Different Time Points",
       x = "Time",
       y = "Mean Value") +
  theme_minimal() +
  theme(legend.position = "none") +

  facet_wrap(~variable, scales = "free_y")






data_long <- melt(data_raw[, c(1, 2, 16:95, grep("lh$", names(data_raw)), grep("rh$", names(data_raw)))],
                  id.vars = c("src_subject_id", "eventname"))

# 使用 aggregate 函数计算每个事件和变量的均值和标准误
summary_df <- aggregate(value ~ eventname + variable, data_long, function(x) {
  mean_value <- mean(x, na.rm = TRUE)
  se_value <- sd(x, na.rm = TRUE) / sqrt(length(x))
  data.frame(mean_value = mean_value, se_value = se_value)
})

summary_mean <- aggregate(value ~ eventname + variable, data_long, mean)
summary_se <- aggregate(value ~ eventname + variable, data_long, sd)

# 移除不再需要的 value 列
summary_df$mean <- summary_mean$value
summary_df$se <- summary_se$value

summary_df$hemi <- ifelse(grepl("lh$", summary_df$variable), "lh", 
                          ifelse(grepl("rh$", summary_df$variable), "rh", NA))

# 去除 variable 后缀的 lh 和 rh
summary_df$variable <- gsub("(lh|rh)$", "", summary_df$variable)
summary_df$value <- NA
summary_df <- summary_df[!is.na(summary_df$hemi),]

# 绘制带误差条的均值图
ggplot(summary_df, aes(x = eventname, y = mean, color = variable, group = hemi, shape = hemi)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  labs(title = "Mean Value with Error Bars at Different Time Points",
       x = "Time",
       y = "Mean Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~gsub("(lh|rh)$", "", variable), scales = "free_y")

ggplot(summary_df, aes(x = eventname, y = mean, color = variable, group = hemi, shape = hemi)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  labs(title = "Mean Value with Error Bars at Different Time Points",
       x = "Time",
       y = "Mean Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~gsub("(lh|rh)$", "", variable), scales = "free_y")





data_residual <- data
for (i in names(data)[16:95]){
  model <- lmer(as.formula(paste0(i, '~sex+Pubertal+Race+ehi_y_ss_scoreb+',
                                  'smri_vol_scs_intracranialv+(1|site_id_l:rel_family_id)+',
                                  '(1|src_subject_id)')), data = data_residual)
  fixed_effects <- predict(model, re.form = NA)
  data_residual[, i] <- data_residual[, i] - fixed_effects
  print(paste('Remove Covariables Effect: ',i ,sep = ""))
}


# x = colMeans(data_residual[data_residual$eventname=='bl',16:49])
x = as.double(colMeans(data[data$eventname=='bl',c(16:49)]));xmodel = 'CT'; hemi = 'left'
clim = ceiling(max(abs(x))*100)/100
plotDat.noNA = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'), 1),
                          val = x, model = rep(c(rep(xmodel, length(x), each = 1)), 1),
                          hemi = rep(c(rep(hemi, length(x))), 1))
plotDat.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(x)),
                        val = NA,model = rep(xmodel,each = 2),
                        hemi = rep(c(rep(hemi, 2*length(x))), 1))
plotDat = rbind(plotDat.noNA,plotDat.NA);
plotDat %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = val),position = "stacked",size = .65) +
  facet_wrap(  ~ model) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 3),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))

# 1032.730 744.545

# data = data[-1*which(data$eventname=="4y"),]
xList = c("interview_age","sex","ehi_y_ss_scoreb","smri_vol_scs_intracranialv")
xList2 = c("interview_age","sexM","ehi_y_ss_scoreb","smri_vol_scs_intracranialv")
fit.demovar <- data.frame()
n = 0
for (i in names(data)[16:95]){
  f = as.formula(paste(i,'~interview_age+sex+Pubertal+Race+ehi_y_ss_scoreb+',
                       'smri_vol_scs_intracranialv+(1|site_id_l:rel_family_id)+',
                       '(1|src_subject_id)',sep = "")
                 )
  fit = lmerTest::lmer(f,data);fit.summary = summary(fit)
  for (j in 1:length(xList)){
    n = n + 1
    fit.demovar[n,c("x")] <- xList[j]
    fit.demovar[n,c("y")] = i
    fit.demovar[n,c("estimate","df", "tval", "pval")] =
      fit.summary$coefficients[xList2[j], c("Estimate","df", "t value", "Pr(>|t|)")]
  }
  print(paste('Linear Mixed Model: ',i,'~sex+Race+ehi_y_ss_scoreb+Pubertal+',
              'smri_vol_scs_intracranialv', sep = ""))
}


xList = c("interview_age")
fit.demovar.raw <- data.frame()
n = 0
for (i in names(data)[16:95]){
  f = as.formula(paste(i,'~interview_age+sex+Pubertal+Race+ehi_y_ss_scoreb+',
                       'smri_vol_scs_intracranialv+(1|site_id_l:rel_family_id)+',
                       '(1|src_subject_id)',sep = "")
  )
  fit = lmerTest::lmer(f,data);fit.summary = summary(fit)
  for (j in 1:length(xList)){
    n = n + 1
    fit.demovar.raw[n,c("x")] <- xList[j]
    fit.demovar.raw[n,c("y")] = i
    fit.demovar.raw[n,c("estimate","df", "tval", "pval")] =
      fit.summary$coefficients[xList2[j], c("Estimate","df", "t value", "Pr(>|t|)")]
  }
  print(paste('Linear Mixed Model: ',i,'~sex+Race+ehi_y_ss_scoreb+Pubertal+',
              'smri_vol_scs_intracranialv', sep = ""))
}







# x = fit.demovar[fit.demovar$x=="interview_age","tval"][34+c(1:34)];
xid = c(1:34)+34
# x0 = fit.demovar[fit.demovar$x=="interview_age","tval"]*
#   (p.adjust(fit.demovar[fit.demovar$x=="interview_age","pval"],method = "fdr")<0.05)
v = "interview_age"
x0 = fit.demovar[fit.demovar$x==v,"tval"]*
  (qvalue(p = fit.demovar[fit.demovar$x==v,"pval"])$qvalues<0.05)

x = x0[xid]
xmodel = 'mu'; hemi = 'left'
clim = ceiling(max(abs(x))*100)/100
plotDat.noNA = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'), 1),
                          val = x, model = rep(c(rep(xmodel, length(x), each = 1)), 1),
                          hemi = rep(c(rep(hemi, length(x))), 1))
plotDat.NA = data.frame(roi = rep(formatC(c(1,5),width = 4,flag = '0'), length(x)),
                        val = NA,model = rep(xmodel,each = 2),
                        hemi = rep(c(rep(hemi, 2*length(x))), 1))
plotDat = rbind(plotDat.noNA,plotDat.NA);
plotDat %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = val),position = "stacked",size = .65) +
  facet_wrap(  ~ model) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
                       na.value = "lightgrey",breaks = seq(-1*clim, clim, 3),
                       limits = c(-1*clim, clim)) +
  theme_brain() + theme(text = element_text(size = 15))

## GWAS
signif_snps_AsymDelta <-
  read.csv("E:/Acdamic/Article_ABCD_Asymmetry/ABCD_GWAS/signif_snps_AsymDelta.txt",sep = " ")
signif_snps_Asym <-
  read.csv("E:/Acdamic/Article_ABCD_Asymmetry/ABCD_GWAS/signif_snps_Asym.txt",sep = " ")

signif_snps_Asym <- signif_snps_Asym[signif_snps_Asym$P<5e-8,]
signif_snps_AsymDelta <- signif_snps_AsymDelta[signif_snps_AsymDelta$P<5e-8,]



n = 0
GeneMaps_List <- list()
for (i in unique(signif_snps_Asym$VAR)) {
  n = n + 1
  hsumstats <- subset(signif_snps_Asym, VAR == i)
  names(hsumstats)[c(1:3, 12)] <- c("chr", "pos", "id", "p")
  if (sum(duplicated(hsumstats))>0){
    hsumstats <- hsumstats[-1*which(duplicated(hsumstats)),]
  }
  snp_sets <- map_snp_to_gene(hsumstats, gene.curated.GRCh37)
  GeneMaps <- snp_sets[["map"]]
  for (j in 1:nrow(GeneMaps)) {
    GeneMaps$var.id[j] <- i
    GeneMaps$gene.name[j] <- gene.curated.GRCh37[
      gene.curated.GRCh37$gene.id == GeneMaps$gene.id[j],]$gene.name
  }
  GeneMaps_List[[n]] <- GeneMaps
}

GeneNames_List <- list()
for (i in 1:length(GeneMaps_List)){
  GeneNames_List[[i]] <- unique(GeneMaps_List[[i]][,c("var.id","gene.name")])
}

GeneNames_Df <- do.call(rbind, GeneNames_List)
GeneNames_Df <- GeneNames_Df[-1*which(is.na(GeneNames_Df$gene.name)),]
GeneMaps_Df <- do.call(rbind, GeneMaps_List)
GeneMaps_Df<-na.omit(GeneMaps_Df)


write.csv(GeneNames_Df,file = "GeneNames_Asym_GRCh37.csv")
write.csv(GeneMaps_Df,file = "GeneMaps_Asym_GRCh37.csv")


n = 0
GeneMaps_List <- list()
for (i in unique(signif_snps_AsymDelta$VAR)) {
  n = n + 1
  hsumstats <- subset(signif_snps_AsymDelta, VAR == i)
  names(hsumstats)[c(1:3, 12)] <- c("chr", "pos", "id", "p")
  if (sum(duplicated(hsumstats))>0){
    hsumstats <- hsumstats[-1*which(duplicated(hsumstats)),]
  }
  snp_sets <- map_snp_to_gene(hsumstats, gene.curated.GRCh37)
  GeneMaps <- snp_sets[["map"]]
  for (j in 1:nrow(GeneMaps)) {
    GeneMaps$var.id[j] <- i
    GeneMaps$gene.name[j] <- gene.curated.GRCh37[
      gene.curated.GRCh37$gene.id == GeneMaps$gene.id[j],]$gene.name
  }
  GeneMaps_List[[n]] <- GeneMaps
}

GeneNames_List <- list()
for (i in 1:length(GeneMaps_List)){
  GeneNames_List[[i]] <- unique(GeneMaps_List[[i]][,c("var.id","gene.name")])
}

GeneNames_Df <- do.call(rbind, GeneNames_List)
GeneNames_Df <- GeneNames_Df[-1*which(is.na(GeneNames_Df$gene.name)),]
GeneMaps_Df <- do.call(rbind, GeneMaps_List)
GeneMaps_Df<-na.omit(GeneMaps_Df)

write.csv(GeneNames_Df,file = "GeneNames_AsymDelta_GRCh37.csv")
write.csv(GeneMaps_Df,file = "GeneMaps_AsymDelta_GRCh37.csv")



# GWAS (CSA / CT * bl / delta)
# Enrichment (Bubbles / BrainSpan)

GeneMaps_Asym_GRCh37 <- read_csv("GeneMaps_Asym_GRCh37.csv", 
                                 col_types = cols(...1 = col_skip()))

GeneMaps_AsymDelta_GRCh37 <- read_csv("GeneMaps_AsymDelta_GRCh37.csv", 
                                      col_types = cols(...1 = col_skip()))

GeneSubset = GeneMaps_AsymDelta_GRCh37[,c("var.id","gene.name")]
GeneSubset = unique(GeneSubset)
GeneSubset = na.omit(GeneSubset)

GeneCode= unique(GeneSubset$var.id)
GeneOverlapMat <- matrix(data = NA,nrow = length(GeneCode),ncol = length(GeneCode))
GeneOverlapList = list()
GeneOverlapSize = c()

GeneSubsetAsym = GeneMaps_Asym_GRCh37[,c("var.id","gene.name")]
GeneSubsetAsym = unique(GeneSubsetAsym);GeneSubsetAsym = na.omit(GeneSubsetAsym)

GeneSubsetDeltaAsym = GeneMaps_AsymDelta_GRCh37[,c("var.id","gene.name")]
GeneSubsetDeltaAsym = unique(GeneSubsetDeltaAsym);GeneSubsetDeltaAsym = na.omit(GeneSubsetDeltaAsym)

# install.packages("venn")

sets <- list(
  CSA_Asym = unique(GeneSubsetAsym$gene.name[(GeneSubsetAsym$var.id<35)]),
  CT_Asym = unique(GeneSubsetAsym$gene.name[(GeneSubsetAsym$var.id>34)&(GeneSubsetAsym$var.id<69)]),
  Vol_Asym = unique(GeneSubsetAsym$gene.name[(GeneSubsetAsym$var.id>68)]),
  CSA_DeltaAsym = unique(GeneSubsetDeltaAsym$gene.name[(GeneSubsetDeltaAsym$var.id<35)]),
  CT_DeltaAsym = unique(GeneSubsetDeltaAsym$gene.name[(GeneSubsetDeltaAsym$var.id>34)&(GeneSubsetDeltaAsym$var.id<69)]),
  Vol_DeltaAsym = unique(GeneSubsetDeltaAsym$gene.name[(GeneSubsetDeltaAsym$var.id>68)])
)

sets <- sets[6:1]

venn(sets,
     zcolor='style', 
     opacity = 0.3,  
     box = F,        
     ilcs = 0.5,     
     sncs = 1        
)

# devtools::install_github("hms-dbmi/UpSetR")

binary_matrix <- function(sets) {
  unique_elements <- unique(unlist(sets))
  binary_data <- sapply(sets, function(x) as.integer(unique_elements %in% x))
  rownames(binary_data) <- unique_elements
  return(binary_data)
}

data_matrix <- binary_matrix(sets)

# UpSet Plots
upset(as.data.frame(data_matrix), 
      sets = names(sets), 
      main.bar.color = "steelblue",
      sets.bar.color = "skyblue",
      mb.ratio = c(0.55, 0.45),
      order.by = "degree",
      decreasing = FALSE,
      text.scale = 2,
      keep.order = TRUE
)


c("Lateral_Ventricle","Inf_Lat_Vent","Thalamus_Proper","Caudate","Putamen",
  "Pallidum","Hippocampus","Amygdala","Accumbens_area","VentralDC")

rm(SNP.Num)
SNP.Num = rep(0,78)
for (i in 1:78){
  SNP.Num[i] = nrow(signif_snps_Asym[signif_snps_Asym$VAR==i,])
}

SNP.Num

SNP.Num = rep(0,78)
for (i in 1:78){
  SNP.Num[i] = nrow(signif_snps_AsymDelta[signif_snps_AsymDelta$VAR==i,])
}
SNP.Num


rm(Gene.Num)
Gene.Num = rep(0,78)
for (i in 1:78){
  Gene.Num[i] = length(unique(GeneMaps_Asym_GRCh37[GeneMaps_Asym_GRCh37$var.id==i,]$gene.name))
}

Gene.Num

Gene.Num = rep(0,78)
for (i in 1:78){
  Gene.Num[i] = length(unique(GeneMaps_AsymDelta_GRCh37[GeneMaps_AsymDelta_GRCh37$var.id==i,]$gene.name))
}
Gene.Num


plotDat = data.frame(roi = rep(formatC(c(2:4, 6:36),width = 4,flag = '0'),2),
                  genenum =Gene.Num[1:68], model = rep(c("CSA","CT"),each = 34),hemi = c(rep("left", 68)))
plotDat.NA = data.frame(roi = formatC(c(1,5),width = 4,flag = '0'),
                     genenum = NA,model = rep(c("CSA","CT"),each = 2),hemi = c(rep("left", 4)))
plotDat = rbind(plotDat,plotDat.NA);clim = 420 # max(plotDat$genenum,na.rm = T)
plotDat %>% 
  ggseg(atlas = dk,colour = "black",hemisphere = "left",
        mapping = aes(fill = genenum),position = "stacked",size = .5) +
  theme(legend.position = "right") +
  facet_grid(model ~.) + 
  scale_fill_gradientn(colours = brewer.pal(9, "Reds"),
                       na.value = "lightgrey",breaks = seq(0, clim, 100),
                       limits = c(0, clim)) +
  theme_brain() + theme(text = element_text(size = 15))



geneList = unique(GeneMaps_Asym_GRCh37$gene.name)
write.csv(geneList,file = "GeneNames_Asym_GRCh37.csv")

geneList = unique(GeneMaps_AsymDelta_GRCh37$gene.name)
write.csv(geneList,file = "GeneNames_AsymDelta_GRCh37.csv")


# install.packages("SAIGE")

GS = read.csv(file = paste("E:/Acdamic/Article_ABCD_Asymmetry/ABCD_GWAS/",
                           "FUMA_gene2func494760/GS.txt",sep = ""),sep = "\t",quote = "  ")


GS$N_genes[GS$Category=="GO_bp"]#size
GS$p[GS$Category=="GO_bp"]#x pos
# color
GS$GeneSet = stringr::str_to_title((gsub("GOCC ","",gsub("GOBP ","",gsub("_"," ",GS$GeneSet)))))
GS$N_overlap = GS$N_overlap/57241
GS$logP = -1*log10(GS$p)

GS_bp <- GS[GS$Category=="GO_bp",];
GS_bp$GeneSet <- factor(GS_bp$GeneSet,ordered = T,levels = GS_bp$GeneSet[order(GS_bp$logP)])
GS_cc <- GS[GS$Category=="GO_cc",];
GS_cc$GeneSet <- factor(GS_cc$GeneSet,ordered = T,levels = GS_cc$GeneSet[order(GS_cc$logP)])

GS_go <- GS[(GS$Category=="GO_bp")|(GS$Category=="GO_cc")|(GS$Category=="GO_mf"),];
GS_go_top_30 <- GS_go %>%
  arrange(p) %>%
  slice(1:30)

ggplot(GS_bp,aes(N_overlap,GeneSet),size=N_genes,color=logP)+
  geom_point(aes(size=N_genes,color=logP))+
  scale_color_gradient(low="gray",high = "red")+ 
  labs(color=expression(p),size="N_genes",
       x="RichFactor",y="Pathway name",
       title="Pathway enrichment") + 
  theme_minimal()

ggplot(GS_cc,aes(N_overlap,GeneSet),size=N_genes,color=logP)+
  geom_point(aes(size=N_genes,color=logP))+
  scale_color_gradient(low="gray",high = "red")+ 
  labs(color=expression(p),size="N_genes",
       x="RichFactor",y="Pathway name",
       title="Pathway enrichment")+
  theme_minimal()


ggplot(GS_go_top_30, aes(x = N_overlap, y = GeneSet, size = N_genes, color = Category)) +
  geom_point() +
  scale_color_manual(values = c("GO_bp" = "red", "GO_cc" = "blue", "GO_mf" = "green")) +
  labs(color = "Category", size = "N_genes",
       x = "RichFactor", y = "Pathway name",
       title = "Pathway enrichment") +
  theme_minimal()



# 600 620
ggplot(GS_go_top_30, aes(x = N_overlap, y = GeneSet, size = N_genes, color = logP, shape = Category)) +
  geom_point(aes(size=N_genes,color=logP)) +
  scale_color_gradient(low = "gray", high = "red") + # 颜色从灰色到红色
  scale_shape_manual(values = c("GO_bp" = 16, "GO_cc" = 18, "GO_mf" = 17)) + # 手动指定形状，23=菱形，21=圆形，24=三角形
  labs(color = "logP", size = "N_genes", shape = "Category",
       x = "RichFactor", y = "Pathway name",
       title = "Pathway enrichment") +
  theme_minimal()


gtex = read.csv(file = paste("E:/Acdamic/Article_ABCD_Asymmetry/ABCD_GWAS/",
                             "FUMA_gene2func494760/gtex_v8_ts_avg_normTPM_exp.txt",sep = ""),sep = "\t",quote = "  ")
gtex = read.csv(file = paste("E:/Acdamic/Article_ABCD_Asymmetry/ABCD_GWAS/",
                             "FUMA_gene2func494760/gtex_v8_ts_general_DEG.txt",sep = ""),sep = "\t",quote = "  ")

gt<-gtex[gtex$Category=="DEG.twoside",c(2:6)]
gt$GeneSet<- factor(gt$GeneSet,ordered = T,levels = gt$GeneSet[order(-1*log10(gt$p))])
ggplot(gt,aes(x =GeneSet,y=N_overlap,fill=-1*log10(p))) + 
  geom_bar(width = 1, stat = "identity") + 
  theme_minimal() + # coord_cartesian(ylim = c(-20,120)) + 
  coord_polar()
# theme(axis.text.x = element_text(angle = 95, hjust = 1)) + 

gtex$Category <- factor(gtex$Category,ordered = T,
                        levels = c("DEG.up","DEG.down","DEG.twoside"))
gtex$GeneSet<- factor(gtex$GeneSet,ordered = T,
                      levels = gt$GeneSet[order(-1*log10(gt$p),decreasing = T)])
#  + scale_fill_brewer(palette = "Set1")
ggplot(gtex,aes(x =GeneSet,y=N_overlap,fill=p< 0.001)) + 
  facet_grid(rows = "Category") + 
  geom_bar(width = .9, stat = "identity") + 
  theme_minimal() + 
  scale_fill_manual(values = c("grey", "red"), 
                    labels = c("p >= 0.001", "p < 0.001")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 40))



## Import Brain span data 
## / START

path = 'E:/Acdamic/Acdamic_Toolset/temporal-brain-expression-master/genes_matrix_csv/';
rows_metadata <- read_csv(paste(path,"rows_metadata.csv",sep = ""))
columns_metadata <- read_csv(paste(path,"columns_metadata.csv",sep = ""))
expression_matrix <- read_csv(paste(path,"expression_matrix.csv",sep = ""), col_names = FALSE)
expression_matrix <- expression_matrix [,-1]

# change var name gene_symbol to GENE
rows_metadata <- rows_metadata %>% select(., 1:3, GENE = gene_symbol, 5)

# merge rows_metadata with expression
dat <- cbind(rows_metadata, expression_matrix)

## Check overall sample size 
samplecount <- columns_metadata %>% group_by(donor_id) %>% count()# %>% arrange(desc(n))

## Check sample size for each of the available tissue * time point 
structure <- columns_metadata %>% group_by(structure_name) %>% count()
timepoint <- columns_metadata %>% group_by(age) %>% count()

## Recode Age into Weeks
columns_metadata <- 
  columns_metadata %>% 
  mutate(Weeks = ifelse(age =="8 pcw",8,
                        ifelse(age =="9 pcw",9, 
                               ifelse(age =="12 pcw",12,
                                      ifelse(age =="13 pcw",13,
                                             ifelse(age =="16 pcw",16,
                                                    ifelse(age =="17 pcw",17,
                                                           ifelse(age =="19 pcw",19,
                                                                  ifelse(age =="21 pcw",21,
                                                                         ifelse(age =="24 pcw",24,
                                                                                ifelse(age =="25 pcw",25,
                                                                                       ifelse(age =="26 pcw",26,
                                                                                              ifelse(age =="35 pcw",35,
                                                                                                     ifelse(age =="37 pcw",37,
                                                                                                            ifelse(age =="4 mos",53,
                                                                                                                   ifelse(age =="10 mos",77,
                                                                                                                          ifelse(age =="1 yrs",89,
                                                                                                                                 ifelse(age =="2 yrs",141,
                                                                                                                                        ifelse(age =="3 yrs",193,
                                                                                                                                               ifelse(age =="4 yrs",245,
                                                                                                                                                      ifelse(age =="8 yrs",453,
                                                                                                                                                             ifelse(age =="11 yrs",609,
                                                                                                                                                                    ifelse(age =="13 yrs",713,
                                                                                                                                                                           ifelse(age =="15 yrs",817,
                                                                                                                                                                                  ifelse(age =="18 yrs",973,
                                                                                                                                                                                         ifelse(age =="19 yrs",1025,
                                                                                                                                                                                                ifelse(age =="21 yrs",1129,
                                                                                                                                                                                                       ifelse(age =="23 yrs",1233,
                                                                                                                                                                                                              ifelse(age =="30 yrs",1597,
                                                                                                                                                                                                                     ifelse(age =="36 yrs",1909,
                                                                                                                                                                                                                            ifelse(age =="37 yrs",1961,
                                                                                                                                                                                                                                   ifelse(age =="40 yrs",2117, "no"))))))))))))))))))))))))))))))))

columns_metadata <- columns_metadata %>% 
  mutate(Sex = ifelse(gender == "M", 1, ifelse(gender == "F", 0, "NA")))


struct.info = unique(timepoint[,c(6:8)])
struct.id.rm <- struct.info$structure_acronym[
  struct.info$structure_id %in% 
    c(10361,10550,10665,10552,10391,10551,10294,10333,10656,10657,10398)]

df.expression_matrix.clean.melt <- 
  df.expression_matrix.clean.melt[!(
    df.expression_matrix.clean.melt$structure_acronym %in% 
      struct.id.rm),]
names(geneList)[1]<-"GENE"
geneList <- merge(geneList,dat[,c(3,4)])


df <- df.expression_matrix.clean.melt[
  which(df.expression_matrix.clean.melt$ensembl_gene_id %in% geneList$ensembl_gene_id),]

df <- df[!(df$structure_acronym %in% struct.id.rm),]
df <- df[(df$ensembl_gene_id %in% geneList$ensembl_gene_id),]

df$age <- factor(df$age,ordered = T, levels = unique(df$age))
# df$Weeks = as.double(df$Weeks)
df <- na.omit(df)
df$structure_acronym <- factor(df$structure_acronym,levels = unique(df$structure_acronym))

df <- df[!df$structure_acronym=="M1C-S1C",]
df <- df[!df$structure_acronym=="PCx",]
df <- df[!df$structure_acronym=="Ocx",]

# STC MFC DFC OFC ITC VFC TCx A1C V1C M1C IPC S1C
df$lobe = NA
df$lobe[df$structure_acronym=="DFC"|df$structure_acronym=="MFC"|
          df$structure_acronym=="OFC"|df$structure_acronym=="VFC"] <- "Fr"
df$lobe[df$structure_acronym=="STC"|df$structure_acronym=="ITC"|
          df$structure_acronym=="A1C"] <- "Tr"
df$lobe[df$structure_acronym=="M1C"|df$structure_acronym=="S1C"] <- "SeMo"
df$lobe[df$structure_acronym=="IPC"] <- "Pr"
df$lobe[df$structure_acronym=="V1C"] <- "Or"


g = ggplot(data = df, aes(x = age, y = value, 
                          group = structure_acronym,
                          color = lobe)) +
  geom_smooth(method = "loess",se = T, alpha = 0.03) + # geom_point() + # scale_color_brewer() +
  labs(title = "Brain Span", x = "Age", y = "Exp", color = "lobe") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 95, hjust = 1)) + coord_cartesian(ylim = c(1.7, 2.3))
g + scale_color_jama()
brewer.pal(n = 5, name = "RdYlBu")

g + scale_color_viridis_d()
g + scale_colour_brewer(palette = "RdYlBu")

g = ggplot(data = df, aes(x = age, y = value, 
                          group = structure_acronym,
                          color = structure_acronym)) +
  geom_smooth(method = "loess",se = T, alpha = 0.05) + # geom_point() + # scale_color_brewer() +
  labs(title = "Brain Span", x = "Age", y = "Exp", color = "lobe") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 95, hjust = 1)) + coord_cartesian(ylim = c(1.7, 2.3))
g


############################################################
# 11 regions: violin plot (high IQ vs low IQ)


data_raw <- read.csv("E:/Acdamic/Article_ABCD_Asymmetry/abcd_fs.csv",na = "",sep = ";")
data_raw$eventname <- recode(data_raw$eventname, 
                         'baseline_year_1_arm_1' = 'bl', 
                         '2_year_follow_up_y_arm_1' = '2y', 
                         '4_year_follow_up_y_arm_1' = '4y')
data_raw$rel_family_id <- as.character(data_raw$rel_family_id)
data_raw$rel_birth_id <- as.character(data_raw$rel_birth_id)
data_raw$school_id <- as.character(data_raw$school_id)
data_raw$district_id <- as.character(data_raw$district_id)
data_raw$smri_vol_scs_intracranialv <- data_raw$smri_vol_scs_intracranialv/10000
ehi_y_ss_scoreb <- data_raw$ehi_y_ss_scoreb
data_raw$ehi_y_ss_scoreb[ehi_y_ss_scoreb==1] <- 3
data_raw$ehi_y_ss_scoreb[ehi_y_ss_scoreb==2] <- 1
data_raw$ehi_y_ss_scoreb[ehi_y_ss_scoreb==3] <- 2



library(tidyverse)
library(lme4)


brain_vars <- c("smri_area_cdk_lobfr","smri_thick_cdk_fusiform","smri_thick_cdk_ihcate",
                "smri_thick_cdk_parahpal","smri_thick_cdk_sufr","smri_vol_scs_tp","smri_vol_scs_pallidum")
brain_vars_hemi <- paste(rep(brain_vars,each = 2),rep(c("lh","rh"),2),sep = "")
brain_vars_hemi


data_raw_bl <- data_raw[data_raw$eventname=='bl',]
data_raw_bl <- cbind(data_raw_bl[,c(1:4,8,12:15,185,176)],data_raw_bl[,brain_vars_hemi])
data_raw_bl <- na.omit(data_raw_bl)
# 
# 
# formula <- "~ interview_age + Pubertal + sex + Race + smri_vol_scs_intracranialv + ehi_y_ss_scoreb + (1 | site_id_l:rel_family_id)"
# get_residuals <- function(df,variable,formula) {
#   model <- lmerTest::lmer(as.formula(paste0(variable, formula)), data = df)
#   residuals <- resid(model) + fixef(model)[1]
#   return(residuals)
# }
# 
# for (var in brain_vars_hemi) {
#   data_raw_bl[[paste0(var, "_resid")]] <- get_residuals(data_raw_bl,var,formula)
# }
# 
# data_raw_bl <- data_raw_bl[,-1*c(2:11,12:33)]
# names(data_raw_bl) <- gsub("_resid","",names(data_raw_bl))
# 

# 分为前30%和后30%两个组
data_raw_bl <- data_raw_bl %>%
  mutate(
    group = case_when(
      nihtbx_cryst_uncorrected <= quantile(nihtbx_cryst_uncorrected, 0.3) ~ "Low",
      nihtbx_cryst_uncorrected >= quantile(nihtbx_cryst_uncorrected, 0.7) ~ "High",
      TRUE ~ "Middle"
    )
  )

data_raw_bl <- data_raw_bl[-1*which(data_raw_bl$group=="Middle"),]

data_raw_left <- data_raw_bl[,paste(brain_vars,"lh",sep = "")];names(data_raw_left) <- brain_vars;data_raw_left$hemi <- "left"
data_raw_right <- data_raw_bl[,paste(brain_vars,"rh",sep = "")];names(data_raw_right) <- brain_vars;data_raw_right$hemi <- "right"
data_raw_stack <- rbind(cbind(data_raw_bl[,c(1:11,26)],data_raw_left),
                        cbind(data_raw_bl[,c(1:11,26)],data_raw_right))
data_raw_stack <- na.omit(data_raw_stack)



data_raw_stack$group <- factor(data_raw_stack$group, 
                               levels = c("Low", "High"), 
                               ordered = TRUE)
data_raw_stack$hemi <- factor(data_raw_stack$hemi, 
                               levels = c("right", "left"), 
                               ordered = TRUE)

data_raw_stack2 <- data_raw_stack
data_raw_stack2$group <- NA
data_raw_stack2$hemi <- NA

data_raw_stack2$group[which(data_raw_stack$group=="High")] <- 1
data_raw_stack2$group[which(data_raw_stack$group=="Low")] <- 0
data_raw_stack2$hemi[which(data_raw_stack$hemi=="left")] <- 1
data_raw_stack2$hemi[which(data_raw_stack$hemi=="right")] <- 0


fit.interact <- data.frame(x = character(),
                           y = character(),
                           estimate = numeric(),
                           df = numeric(),
                           tval = numeric(),
                           pval = numeric(),
                           stringsAsFactors = FALSE)

n = 0

# 循环拟合所有变量的模型
for (i in 1:length(brain_vars)) {
  formula <- as.formula(paste(brain_vars[i], "~ group * hemi + (1|src_subject_id)"))
  model <- lmer(formula, data = data_raw_stack)
  fit.summary <- summary(model)
  
  # 提取主效应和交互效应的统计量
  terms <- c("group.L", "hemi.L", "group.L:hemi.L")
  
  for (term in terms) {
    # 提取系数和p值
    term_info <- fit.summary$coefficients[term, ]
    n = n + 1
    
    fit.interact[n, "y"] <- brain_vars[i]
    fit.interact[n, "x"] <- term
    fit.interact[n, "estimate"] <- term_info["Estimate"]
    fit.interact[n, "df"] <- term_info["df"]
    fit.interact[n, "tval"] <- term_info["t value"] # 近似F值
    fit.interact[n, "pval"] <- term_info["Pr(>|t|)"]
  }
}

fit.interact
write.csv(fit.interact,"fit.interact5.csv")


Vio.plotDat <- data_raw_bl %>%
  gather(key = "variable", value = "residual", starts_with("smri_")) %>%
  mutate(
    measure = sub("(.*)(lh|rh)$", "\\1", variable),
    hemi = sub(".*(lh|rh)$", "\\1", variable),
    measure_type = case_when(
      grepl("area", measure) ~ "Area",
      grepl("thick", measure) ~ "Thickness",
      grepl("vol", measure) ~ "Volume",
      TRUE ~ NA_character_
    )
  )
Vio.plotDat <- Vio.plotDat[Vio.plotDat$variable!="smri_vol_scs_intracranialv",]


Vio.plotDat$hemi_vio <- gsub("h","",Vio.plotDat$hemi)

Vio.plotDat <- Vio.plotDat[Vio.plotDat$group!="Middle",]
ggplot(Vio.plotDat, 
       aes(x = group, y = residual, fill = group, alpha = hemi)) +
  geom_half_violin(data = Vio.plotDat[Vio.plotDat$hemi=="lh",],side = "l") +
  geom_half_violin(data = Vio.plotDat[Vio.plotDat$hemi=="rh",],side = "r") +
  scale_alpha_manual(values = c("lh" = 1, "rh" = 0.5)) +
  geom_half_boxplot(data = Vio.plotDat[Vio.plotDat$hemi=="lh",],side = "l",fill = "white", alpha = 1, outlier.shape = NA,show.legend = T) +
  geom_half_boxplot(data = Vio.plotDat[Vio.plotDat$hemi=="rh",],side = "r",fill = "white", alpha = 1, outlier.shape = NA,show.legend = T) +
  scale_fill_manual(values = c("Low" = "red", "High" = "blue")) +
  facet_wrap(~measure, scales = "free_y") +
  labs(
    title = "Residuals of Brain Regions by Hemisphere and Group",
    x = "Brain Region",
    y = "Residuals",
    fill = "Group",
    alpha = "Hemi"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  stat_summary(data = Vio.plotDat[Vio.plotDat$hemi=="lh",],fun = mean, geom = "line", aes(group = measure), alpha = 1, size = .5) +
  stat_summary(data = Vio.plotDat[Vio.plotDat$hemi=="rh",],fun = mean, geom = "line", aes(group = measure), alpha = .5, size = .5) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2)
