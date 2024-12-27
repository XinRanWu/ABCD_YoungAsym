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

############################################################
# 8 regions: violin plot (high IQ vs low IQ)


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
