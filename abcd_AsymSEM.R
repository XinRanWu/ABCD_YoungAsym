# ABCD CLPM (Asymmetry & Cognition)
# 加载需要的包
library(dplyr)
library(tidyr)
library(lavaan)
library(lme4)
library(readr)

mri_y_tmf <- read_csv("E:/Acdamic/Data_ABCD/abcd-data-release-5.0/mri_y_tmf.csv")
mri_y_tmf_cc <- mri_y_tmf[,c(1,2,seq(21,210,42))]
mri_y_tmf_cc$eventname <- recode(mri_y_tmf_cc$eventname, 
                         'baseline_year_1_arm_1' = 'bl', 
                         '2_year_follow_up_y_arm_1' = '2y', 
                         '4_year_follow_up_y_arm_1' = '4y')

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

data <- data[rowSums(is.na(data[,c("interview_age","sex","Race","Pubertal","ehi_y_ss_scoreb",
                                   "smri_vol_scs_intracranialv","site_id_l","rel_family_id")]))==0,]

data <- merge(data,abcd_brainPheno[,c(1,2,97:125)])

############################################################

# 1. 模型拟合效果
# 模型拟合效果可以通过一系列拟合度指标来评估。主要关注以下几个拟合度指标：
# 
# Chi-Square Test of Model Fit (χ²)：评估模型拟合数据的好坏。一般来说，较小的χ²值表明模型拟合良好，但χ²值对样本量非常敏感。
# Root Mean Square Error of Approximation (RMSEA)：RMSEA ≤ 0.06 表示模型拟合良好；0.06 < RMSEA ≤ 0.08 表示拟合可接受；RMSEA > 0.08 表示模型拟合不佳。
# Comparative Fit Index (CFI) 和 Tucker-Lewis Index (TLI)：这两个指标一般在0到1之间，越接近1表示模型拟合越好。通常，CFI和TLI ≥ 0.95表示模型拟合良好；0.90 ≤ CFI和TLI < 0.95表示拟合可以接受。
# Standardized Root Mean Square Residual (SRMR)：SRMR ≤ 0.08 表示模型拟合良好。
# 2. 路径系数和中介效应
# 除了拟合度指标，还需要关注路径系数和中介效应的显著性。具体来看：
# 
# 路径系数（Path Coefficients）：查看每条路径（例如X -> M, M -> Y, X -> Y）的标准化系数和P值。P值小于0.05通常表示路径系数显著。
# 中介效应：主要关注间接效应的估计值及其显著性。间接效应通常通过路径系数的乘积表示（例如a*b）。可以通过标准误和P值来判断间接效应是否显著。

## Mediation of Asymmetry

# 定义和拟合中介模型
model_template <- '
  # X -> M 
  M_bl ~ a*X_bl
  
  # M -> Y 
  Y_2y ~ b*M_bl
  
  # X -> Y direct path
  Y_2y ~ c*X_bl
  
  # indirect path (a*b)
  indirect := a * b
  
  # total path
  total := c + (a * b)
'

results <- data.frame(
  X = character(),M = character(),Y = character(),
  ChiSquare = numeric(),df = numeric(),pvalue = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_pvalue = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric(),
  X_to_M_Estimate = numeric(),X_to_M_StdErr = numeric(),X_to_M_zvalue = numeric(),X_to_M_pvalue = numeric(),
  M_to_Y_Estimate = numeric(),M_to_Y_StdErr = numeric(),M_to_Y_zvalue = numeric(),M_to_Y_pvalue = numeric(),
  X_to_Y_Estimate = numeric(),X_to_Y_StdErr = numeric(),X_to_Y_zvalue = numeric(),X_to_Y_pvalue = numeric(),
  indirect_Estimate = numeric(),indirect_StdErr = numeric(),indirect_zvalue = numeric(),indirect_pvalue = numeric(),
  total_Estimate = numeric(),total_StdErr = numeric(),total_zvalue = numeric(),total_pvalue = numeric()
)

id = 0
for (Yvar in c('cbcl_scr_syn_totprob_r')){
  for (Xvar in c("latent_factor_ss_general_ses")){
    for (Mvar in names(abcd_brainPheno)[17:96]){
      
      id =id + 1
      
      data_residual <- abcd_brainPheno[,c("src_subject_id","eventname","interview_age","Pubertal",
                                          "sex","Race","smri_vol_scs_intracranialv","ehi_y_ss_scoreb",
                                          "site_id_l","rel_family_id",Xvar, Mvar, Yvar)]
      
      data_residual$eventname <- as.character(data_residual$eventname)
      
      colnames(data_residual) <- sub(Yvar, "Y", colnames(data_residual))
      colnames(data_residual) <- sub(Xvar, "X", colnames(data_residual))
      colnames(data_residual) <- sub(Mvar, "M", colnames(data_residual))
      
      
      data_residual_bl <- na.omit(data_residual[data_residual$eventname=='bl',])
      data_residual_2y <- na.omit(data_residual[data_residual$eventname=='2y',-11])
      
      data_residual_bl[, "M"] <- residuals(lmer(M~interview_age+Pubertal+sex+Race+
                                                  smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                  (1|site_id_l:rel_family_id),data = data_residual_bl))
      data_residual_2y[, "M"] <- residuals(lmer(M~interview_age+Pubertal+sex+Race+
                                                  smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                  (1|site_id_l:rel_family_id),data = data_residual_2y))
      data_residual_2y$X <- NA
      
      data_residual <- rbind(data_residual_bl,data_residual_2y)
      
      duplicates <- data_residual %>%
        group_by(src_subject_id, eventname) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n > 1L)
      
      data_residual <- data_residual %>%
        group_by(src_subject_id, eventname) %>%
        summarise(across(c(X, M, Y), mean, na.rm = TRUE), .groups = "drop")
      
      data_residual <- data_residual %>%
        mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
      
      # 为lavaan包准备数据                                                 
      data_wide <- data_residual %>%
        dplyr::select(src_subject_id, eventname, X, M, Y) %>%
        pivot_wider(names_from = eventname, values_from = c(X, M, Y))
      
      # 确保所有变量都存在并且没有缺失值
      required_vars <- c("M_bl", "Y_2y", "X_bl")
      
      data_wide_filtered <- data_wide %>%
        filter(complete.cases(select(., all_of(required_vars))))
      
      fit <- sem(model_template, data = data_wide_filtered)
      
      # summary(fit, fit.measures = TRUE, standardize = TRUE)
      
      # 提取拟合度指标
      fit_measures <- fitMeasures(fit)
      chi_square <- fit_measures["chisq"]
      df <- fit_measures["df"]
      pvalue <- fit_measures["pvalue"]
      rmsea <- fit_measures["rmsea"]
      rmsea_ci <- fit_measures["rmsea.ci.lower"] + fit_measures["rmsea.ci.upper"]
      rmsea_pvalue <- fit_measures["rmsea.pvalue"]
      cfi <- fit_measures["cfi"]
      tli <- fit_measures["tli"]
      srmr <- fit_measures["srmr"]
      
      # 提取路径系数和中介效应
      path_estimates <- parameterEstimates(fit, standardized = TRUE)
      x_to_m <- path_estimates[path_estimates$lhs == "M_bl" & path_estimates$rhs == "X_bl", ]
      m_to_y <- path_estimates[path_estimates$lhs == "Y_2y" & path_estimates$rhs == "M_bl", ]
      x_to_y <- path_estimates[path_estimates$lhs == "Y_2y" & path_estimates$rhs == "X_bl", ]
      indirect <- path_estimates[path_estimates$lhs == "indirect", ]
      total <- path_estimates[path_estimates$lhs == "total", ]
      
      # 将结果添加到数据框
      results <- rbind(results, data.frame(
        X = Xvar,M = Mvar,Y = Yvar,
        ChiSquare = chi_square,df = df,pvalue = pvalue,
        RMSEA = rmsea,RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],RMSEA_pvalue = rmsea_pvalue,
        CFI = cfi,TLI = tli,SRMR = srmr,
        X_to_M_Estimate = x_to_m$est,X_to_M_StdErr = x_to_m$se,X_to_M_zvalue = x_to_m$z,X_to_M_pvalue = x_to_m$pvalue,
        M_to_Y_Estimate = m_to_y$est,M_to_Y_StdErr = m_to_y$se,M_to_Y_zvalue = m_to_y$z,M_to_Y_pvalue = m_to_y$pvalue,
        X_to_Y_Estimate = x_to_y$est,X_to_Y_StdErr = x_to_y$se,X_to_Y_zvalue = x_to_y$z,X_to_Y_pvalue = x_to_y$pvalue,
        indirect_Estimate = indirect$est,indirect_StdErr = indirect$se,indirect_zvalue = indirect$z,indirect_pvalue = indirect$pvalue,
        total_Estimate = total$est,total_StdErr = total$se,total_zvalue = total$z,total_pvalue = total$pvalue
      ))
      
      print(paste('Mediation: ',Xvar,' - ',Mvar,' - ',Yvar, sep = ""))
      
    }
  }
}


results[which((results$total_pvalue<0.05)&
                (results$indirect_pvalue<0.05)&
                (results$X_to_M_pvalue<0.05)&
                (results$M_to_Y_pvalue<0.05)),
        c(1:3,13:33)]

## Mediation of Asymmetry changing rate

############################################################

# 常用的拟合指标及其参考值:
# CFI (Comparative Fit Index, 比较拟合指数):CFI ≥ 0.95 很好的拟合。0.90 ≤ CFI < 0.95 可以接受的拟合。CFI < 0.90 差拟合。
# TLI (Tucker-Lewis Index, 塔克-刘易斯指数):TLI ≥ 0.95 好拟合。、0.90 ≤ TLI < 0.95 可以接受的拟合。TLI < 0.90 差拟合。
# RMSEA (Root Mean Square Error of Approximation, 近似均方根残差):
# RMSEA ≤ 0.05 好拟合。0.05 < RMSEA ≤ 0.08 可以接受的拟合。0.08 < RMSEA ≤ 0.10 较差的拟合。RMSEA > 0.10 非常差的拟合。 
# 通常情况下,如果一个模型的CFI和TLI大于0.95,且RMSEA小于0.06,则可以认为该模型拟合数据很好。但是,这些只是一般的经验参考值,在实际应用中还需要结合理论和数据的具体情况进行综合判断。
# 注意事项:CFI和TLI越接近1越好,RMSEA越接近0越好。SRMR也可以一并考虑,模型参数的显著性、理论依据以及交叉验证也很重要。
# 对于复杂模型,适度的拟合可能已经可以接受。
# 
# 参数估计值:包括回归系数、协方差/相关性等的标准化和非标准化估计值。
# 对于回归系数,正值表示正向影响,负值表示负向影响。
# 标准误差:
# 每个参数估计值对应的标准误差。
# z值和p值:z值是参数估计值除以标准误差的比值,反映了该参数是否显著不同于0。p值是基于z值计算的概率值,小于0.05表示在5%的水平上显著。
# 参数约束:一些参数可能被约束为相等,反映特定的理论假设。

# 路径A_2y ~ c1*B_bl表示基线时B影响2年后的A(B→A)
# 如果该路径系数c1显著不等于0,则支持B影响A
# 路径B_2y ~ d1*A_bl表示基线时A影响2年后的B(A→B)
# 如果该路径系数d1显著不等于0,则支持A影响B
# 通常,如果A→B和B→A的交叉滞后影响路径均显著,则表明A和B之间可能存在相互影响。如果只有单向路径显著,则支持单向影响的假设。

## 3 time points CLPM

# Good fit for a CLPM model
# CFI> 0.95, TLI> 0.95, RMSEA<0.06 and SRMR< 0.08.

# 初始化结果数据框
results_fit <- list()

results_fit_indices <- data.frame(
  A = character(),B = character(),
  ChiSquare = numeric(),DF = numeric(),P_value = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_P_value = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric()
)

results_path_estimates <- data.frame(
  A = character(),B = character(),
  lhs = character(),op = character(),rhs = character(),
  est = numeric(),se = numeric(),z = numeric(),pvalue = numeric(),std.lv = numeric(),std.all = numeric()
)

model_syntax <- "
  # 自回归路径
  A_2y ~ a1*A_bl
  A_4y ~ a2*A_2y
  
  B_2y ~ b1*B_bl  
  B_4y ~ b2*B_2y
  
  # 交叉滞后路径
  A_2y ~ c1*B_bl
  A_4y ~ c2*B_2y
  
  B_2y ~ d1*A_bl
  B_4y ~ d2*A_2y
  
  # 同时相关
  A_bl ~~ B_bl
  A_2y ~~ B_2y
  A_4y ~~ B_4y
"


# 循环构建三个双向关系模型
id = 0

for (Bvar in c('sds_p_ss_total','nihtbx_cryst_uncorrected','cbcl_scr_syn_attention_r')) {
  
  for (Avar in names(data)[16:95]) {
    
    id = id + 1;
    
    data_residual <- data
    
    colnames(data_residual) <- sub(Bvar, "B", colnames(data_residual))
    colnames(data_residual) <- sub(Avar, "A", colnames(data_residual))
    
    data_residual_bl <- data_residual[data_residual$eventname=='bl',]
    data_residual_2y <- data_residual[data_residual$eventname=='2y',]
    data_residual_4y <- data_residual[data_residual$eventname=='4y',]
    
    data_residual_bl[, "A"] <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                                smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                (1|site_id_l:rel_family_id),data = data_residual_bl))
    data_residual_2y[, "A"] <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                                smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                (1|site_id_l:rel_family_id),data = data_residual_2y))
    data_residual_4y[, "A"] <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                                smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                (1|site_id_l:rel_family_id),data = data_residual_4y))
    
    data_residual <- rbind(data_residual_bl,data_residual_2y,data_residual_4y)
    
    duplicates <- data_residual %>%
      group_by(src_subject_id, eventname) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1L)
    
    data_residual <- data_residual %>%
      group_by(src_subject_id, eventname) %>%
      summarise(across(c(A, B), mean, na.rm = TRUE), .groups = "drop")
    
    data_residual <- data_residual %>%
      mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
    
    # 为lavaan包准备数据                                                 
    data_wide <- data_residual %>%
      dplyr::select(src_subject_id, eventname, A, B) %>%
      pivot_wider(names_from = eventname, 
                  values_from = c(A, B))
    
    
    # 确保数据集不含缺失值的行
    required_vars <- c("A_bl", "A_2y", "A_4y", "B_bl", "B_2y", "B_4y")
    data_wide_filtered <- data_wide %>%
      filter(complete.cases(select(., all_of(required_vars))))
    
    data_wide_standardized <- data_wide_filtered %>%
      mutate(across(all_of(required_vars), scale))
    
    fit <- lavaan::sem(model_syntax, data = data_wide_standardized)
    
    fit_measures <- fitMeasures(fit)
    parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
    
    fit_indices <- data.frame(
      A = Avar,B = Bvar,
      ChiSquare = fit_measures["chisq"],
      DF = fit_measures["df"],
      P_value = fit_measures["pvalue"],
      RMSEA = fit_measures["rmsea"],
      RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
      RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
      RMSEA_P_value = fit_measures["rmsea.pvalue"],
      CFI = fit_measures["cfi"],
      TLI = fit_measures["tli"],
      SRMR = fit_measures["srmr"]
    )
    
    results_fit[[id]] <- fit 
    
    results_fit_indices <- rbind(results_fit_indices, fit_indices)
    
    # 提取路径估计
    parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
    
    # 将路径估计添加到 results_path_estimates 数据框
    path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                              "est", "se", "z", "pvalue", "std.lv", "std.all")]
    path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                       B = rep(Bvar,nrow(path_estimates))),
                            path_estimates)
    
    results_path_estimates <- rbind(results_path_estimates, path_estimates)
    
    print(paste('3-Wave Cross-Lagging Panel Model: ',Avar,' - ',Bvar, sep = ""))
    
  }
}


results_path_estimates[(results_path_estimates$lhs=='A_2y')&
                         (results_path_estimates$rhs=='B_bl')&
                         (results_path_estimates$pvalue<0.05),]

results_path_estimates[(results_path_estimates$lhs=='B_2y')&
                         (results_path_estimates$rhs=='A_bl')&
                         (results_path_estimates$pvalue<0.05),]

results_path_estimates[(results_path_estimates$lhs=='A_4y')&
                         (results_path_estimates$rhs=='B_2y')&
                         (results_path_estimates$pvalue<0.05),]

results_path_estimates[(results_path_estimates$lhs=='B_4y')&
                         (results_path_estimates$rhs=='A_2y')&
                         (results_path_estimates$pvalue<0.05),]



results_path_estimates[(results_path_estimates$A=='smri_area_cdk_cuneus')&
                         (results_path_estimates$B=='nihtbx_cryst_uncorrected'),]


## 2 wave CLPM

# 初始化结果数据框
results_fit <- list()

results_fit_indices <- data.frame(
  A = character(),B = character(),
  ChiSquare = numeric(),DF = numeric(),P_value = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_P_value = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric()
)

results_path_estimates <- data.frame(
  A = character(),B = character(),
  lhs = character(),op = character(),rhs = character(),
  est = numeric(),se = numeric(),z = numeric(),pvalue = numeric(),std.lv = numeric(),std.all = numeric()
)

model_syntax <- "
  # 自回归路径
  A_2y ~ a1*A_bl
  B_2y ~ b1*B_bl  
  
  # 交叉滞后路径
  A_2y ~ c1*B_bl
  B_2y ~ d1*A_bl
  
  # 同时相关
  A_bl ~~ B_bl
  A_2y ~~ B_2y
"

id = 0

for (Bvar in c('nihtbx_cryst_uncorrected','cbcl_scr_syn_totprob_r')) {
  for (Avar in names(abcd_brainPheno)[17:96]) {
    id = id + 1;
    
    data_residual <- na.omit(abcd_brainPheno[,c('src_subject_id','eventname','interview_age','Pubertal','sex',
                                                'Race','smri_vol_scs_intracranialv','ehi_y_ss_scoreb','site_id_l',
                                                'rel_family_id',Avar,Bvar)])
    
    colnames(data_residual) <- sub(Bvar, "B", colnames(data_residual))
    colnames(data_residual) <- sub(Avar, "A", colnames(data_residual))
    
    data_residual_bl <- data_residual[data_residual$eventname=='bl',]
    data_residual_2y <- data_residual[data_residual$eventname=='2y',]
    
    data_residual_bl$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                           smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                           (1|site_id_l:rel_family_id),data = data_residual_bl))
    data_residual_2y$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                           smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                           (1|site_id_l:rel_family_id),data = data_residual_2y))
    
    data_residual <- rbind(data_residual_bl,data_residual_2y)
    
    duplicates <- data_residual %>%
      group_by(src_subject_id, eventname) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1L)
    
    data_residual <- data_residual %>%
      group_by(src_subject_id, eventname) %>%
      summarise(across(c(A, B), mean, na.rm = TRUE), .groups = "drop")
    
    data_residual <- data_residual %>%
      mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
    
    # 为lavaan包准备数据                                                 
    data_wide <- data_residual %>%
      dplyr::select(src_subject_id, eventname, A, B) %>%
      pivot_wider(names_from = eventname, 
                  values_from = c(A, B))
    
    
    # 确保数据集不含缺失值的行
    names(data_wide)[2:5] <- c("A_bl", "A_2y", "B_bl", "B_2y")
    required_vars <- c("A_bl", "A_2y", "B_bl", "B_2y")
    data_wide_filtered <- data_wide %>%
      filter(complete.cases(select(., all_of(required_vars))))
    
    data_wide_standardized <- data_wide_filtered %>%
      mutate(across(all_of(required_vars), scale))
    
    fit <- lavaan::sem(model_syntax, data = data_wide_standardized)
    
    fit_measures <- fitMeasures(fit)
    parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
    
    fit_indices <- data.frame(
      A = Avar,B = Bvar,
      ChiSquare = fit_measures["chisq"],
      DF = fit_measures["df"],
      P_value = fit_measures["pvalue"],
      RMSEA = fit_measures["rmsea"],
      RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
      RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
      RMSEA_P_value = fit_measures["rmsea.pvalue"],
      CFI = fit_measures["cfi"],
      TLI = fit_measures["tli"],
      SRMR = fit_measures["srmr"]
    )
    
    results_fit[[id]] <- fit 
    
    results_fit_indices <- rbind(results_fit_indices, fit_indices)
    
    # 提取路径估计
    parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
    
    # 将路径估计添加到 results_path_estimates 数据框
    path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                              "est", "se", "z", "pvalue", "std.lv", "std.all")]
    path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                       B = rep(Bvar,nrow(path_estimates))),
                            path_estimates)
    
    results_path_estimates <- rbind(results_path_estimates, path_estimates)
    
    print(paste('2-Wave Cross-Lagging Panel Model: ',Avar,' - ',Bvar, sep = ""))
    
  }
}

results_path_estimates[(results_path_estimates$lhs=='A_2y')&
                         (results_path_estimates$rhs=='B_bl')&
                         (results_path_estimates$pvalue<0.05),]

results_path_estimates[(results_path_estimates$lhs=='B_2y')&
                         (results_path_estimates$rhs=='A_bl')&
                         (results_path_estimates$pvalue<0.05),]


ggseg.hemi(results_path_estimates$z[(results_path_estimates$B=="cbcl_scr_syn_totprob_r")&
                                      (results_path_estimates$lhs=='A_2y')&
                                      (results_path_estimates$rhs=='B_bl')][1:68],c(rep("CSA",34),rep("CT",34)),7)
ggseg.hemi(results_path_estimates$z[(results_path_estimates$B=="cbcl_scr_syn_totprob_r")&
                                      (results_path_estimates$lhs=='B_2y')&
                                      (results_path_estimates$rhs=='A_bl')][1:68],c(rep("CSA",34),rep("CT",34)),7)


# Bvar in c('nihtbx_cryst_uncorrected','cbcl_scr_syn_totprob_r')) {
# for (Avar in names(abcd_brainPheno)[17:96]
B_to_A_id = which((results_path_estimates$B=="cbcl_scr_syn_totprob_r")&
                    (results_path_estimates$lhs=='A_2y')&
                    (results_path_estimates$rhs=='B_bl'))

A_to_B_id = which((results_path_estimates$B=="cbcl_scr_syn_totprob_r")&
                    (results_path_estimates$lhs=='B_2y')&
                    (results_path_estimates$rhs=='A_bl'))


B_to_A_id = which((results_path_estimates$B=="nihtbx_cryst_uncorrected")&
                    (results_path_estimates$lhs=='A_2y')&
                    (results_path_estimates$rhs=='B_bl'))

A_to_B_id = which((results_path_estimates$B=="nihtbx_cryst_uncorrected")&
                    (results_path_estimates$lhs=='B_2y')&
                    (results_path_estimates$rhs=='A_bl'))

x = results_path_estimates$z[A_to_B_id]*(qvalue::qvalue(results_path_estimates$pvalue[A_to_B_id])$qvalue<0.05)
ggseg.hemi(x[1:68],c(rep("CSA",34),rep("CT",34)),7)

x = results_path_estimates$z[B_to_A_id]*(qvalue::qvalue(results_path_estimates$pvalue[B_to_A_id])$qvalue<0.05)
ggseg.hemi(x[1:68],c(rep("CSA",34),rep("CT",34)),7)

x = results_path_estimates$z[A_to_B_id]*(results_path_estimates$pvalue[A_to_B_id]<0.05)
ggseg.hemi(x[1:68],c(rep("CSA",34),rep("CT",34)),7)

x = results_path_estimates$z[B_to_A_id]*(results_path_estimates$pvalue[B_to_A_id]<0.05)
ggseg.hemi(x[1:68],c(rep("CSA",34),rep("CT",34)),7)

############################################################

save.image('abcd_AsymSEM.RData')

############################################################

## Plot
# 加载需要的包
library(semPlot)

# 假设您已经估计了模型并将结果保存在fit对象中

# 绘制路径图
semPath(fit,
        layout = 'tree2',  # 布局类型
        edge.label.cex = 0.6, # 路径系数字体大小
        sizeMan = 8, # 潜在变量节点大小  
        sizeXlat = 8, # 标识符节点大小
        residuals = FALSE, # 不显示残差
        residuals.label = FALSE, 
        rotation = 2, # 旋转角度
        layout.eta = 1.5, # 调整布局参数
        layout.tau = 1.5,
        edge.color = 'black', # 路径线条颜色
        hashingDeptor = 'c34', #单箭头箭尾类型
        hashingDeptee = 'boxblank',#单箭头箭头类型
        whatGradients = c(1, 2), # 渐变色
        edge.label.color = 'black', # 路径系数颜色
        label.font = 2.5, # 标签字体大小
        nCharNodes = 30  # 单词截断长度
)

# 保存路径图为PNG文件
png(file = "path_diagram.png", width = 8, height = 6, units = "in", res = 300)
semPath(fit, ...)  # 替换为您之前使用的semPath函数及其参数
dev.off()











################################################################################
id = 0

Aid = 0
Bid = 0

for (Bid in 1:length(names(data)[16:95])) {
  
  for (Aid in 1:length(names(data)[16:95])) {
    
    if (Aid > Bid){
      id = id + 1;
      Avar = names(data)[16:95][Aid]
      Bvar = names(data)[16:95][Bid]
      
      data_residual <- data
      
      colnames(data_residual) <- sub(Bvar, "B", colnames(data_residual))
      colnames(data_residual) <- sub(Avar, "A", colnames(data_residual))
      
      data_residual_bl <- data_residual[data_residual$eventname=='bl',]
      data_residual_2y <- data_residual[data_residual$eventname=='2y',]
      
      data_residual_bl$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                             smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                             (1|site_id_l:rel_family_id),data = data_residual_bl))
      data_residual_2y$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                             smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                             (1|site_id_l:rel_family_id),data = data_residual_2y))
      
      data_residual <- rbind(data_residual_bl,data_residual_2y)
      
      duplicates <- data_residual %>%
        group_by(src_subject_id, eventname) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n > 1L)
      
      data_residual <- data_residual %>%
        group_by(src_subject_id, eventname) %>%
        summarise(across(c(A, B), mean, na.rm = TRUE), .groups = "drop")
      
      data_residual <- data_residual %>%
        mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
      
      # 为lavaan包准备数据                                                 
      data_wide <- data_residual %>%
        dplyr::select(src_subject_id, eventname, A, B) %>%
        pivot_wider(names_from = eventname, 
                    values_from = c(A, B))
      
      
      # 确保数据集不含缺失值的行
      required_vars <- c("A_bl", "A_2y", "B_bl", "B_2y")
      data_wide_filtered <- data_wide %>%
        filter(complete.cases(select(., all_of(required_vars))))
      
      data_wide_standardized <- data_wide_filtered %>%
        mutate(across(all_of(required_vars), scale))
      
      fit <- lavaan::sem(model_syntax, data = data_wide_standardized)
      
      fit_measures <- fitMeasures(fit)
      parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
      
      fit_indices <- data.frame(
        A = Avar,B = Bvar,
        ChiSquare = fit_measures["chisq"],
        DF = fit_measures["df"],
        P_value = fit_measures["pvalue"],
        RMSEA = fit_measures["rmsea"],
        RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
        RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
        RMSEA_P_value = fit_measures["rmsea.pvalue"],
        CFI = fit_measures["cfi"],
        TLI = fit_measures["tli"],
        SRMR = fit_measures["srmr"]
      )
      
      results_fit[[id]] <- fit 
      
      results_fit_indices <- rbind(results_fit_indices, fit_indices)
      
      # 提取路径估计
      parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
      
      # 将路径估计添加到 results_path_estimates 数据框
      path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                                "est", "se", "z", "pvalue", "std.lv", "std.all")]
      path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                         B = rep(Bvar,nrow(path_estimates))),
                              path_estimates)
      
      results_path_estimates <- rbind(results_path_estimates, path_estimates)
      
      print(paste('2-Wave Cross-Lagging Panel Model: ',Avar,' - ',Bvar, sep = ""))
      
    }
    
  }
}

results_path_estimates[(results_path_estimates$lhs=='A_2y')&
                         (results_path_estimates$rhs=='B_bl')&
                         (results_path_estimates$pvalue<5^-80),]

results_path_estimates[(results_path_estimates$lhs=='B_2y')&
                         (results_path_estimates$rhs=='A_bl')&
                         (results_path_estimates$pvalue<5^-80),]


################################################################################
## 2 wave CLPM

# 初始化结果数据框
results_fit <- list()

results_fit_indices <- data.frame(
  A = character(),B = character(),
  ChiSquare = numeric(),DF = numeric(),P_value = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_P_value = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric()
)

results_path_estimates <- data.frame(
  A = character(),B = character(),
  lhs = character(),op = character(),rhs = character(),
  est = numeric(),se = numeric(),z = numeric(),pvalue = numeric(),std.lv = numeric(),std.all = numeric()
)

model_syntax <- "
  # 自回归路径
  A_2y ~ a1*A_bl
  B_2y ~ b1*B_bl  
  
  # 交叉滞后路径
  A_2y ~ c1*B_bl
  B_2y ~ d1*A_bl
  
  # 同时相关
  A_bl ~~ B_bl
  A_2y ~~ B_2y
"

id = 0

for (Bvar in names(mri_y_tmf_cc)[3:7]) {
  
  for (Avar in names(data)[16:95]) {
    
    id = id + 1;
    
    data_residual <- data
    
    colnames(data_residual) <- sub(Bvar, "B", colnames(data_residual))
    colnames(data_residual) <- sub(Avar, "A", colnames(data_residual))
    
    data_residual_bl <- data_residual[data_residual$eventname=='bl',]
    data_residual_2y <- data_residual[data_residual$eventname=='2y',]
    
    data_residual_bl$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                           smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                           (1|site_id_l:rel_family_id),data = data_residual_bl))
    data_residual_2y$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                           smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                           (1|site_id_l:rel_family_id),data = data_residual_2y))
    
    data_residual <- rbind(data_residual_bl,data_residual_2y)
    
    duplicates <- data_residual %>%
      group_by(src_subject_id, eventname) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1L)
    
    data_residual <- data_residual %>%
      group_by(src_subject_id, eventname) %>%
      summarise(across(c(A, B), mean, na.rm = TRUE), .groups = "drop")
    
    data_residual <- data_residual %>%
      mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
    
    # 为lavaan包准备数据                                                 
    data_wide <- data_residual %>%
      dplyr::select(src_subject_id, eventname, A, B) %>%
      pivot_wider(names_from = eventname, 
                  values_from = c(A, B))
    
    
    # 确保数据集不含缺失值的行
    required_vars <- c("A_bl", "A_2y", "B_bl", "B_2y")
    data_wide_filtered <- data_wide %>%
      filter(complete.cases(select(., all_of(required_vars))))
    
    data_wide_standardized <- data_wide_filtered %>%
      mutate(across(all_of(required_vars), scale))
    
    fit <- lavaan::sem(model_syntax, data = data_wide_standardized)
    
    fit_measures <- fitMeasures(fit)
    parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
    
    fit_indices <- data.frame(
      A = Avar,B = Bvar,
      ChiSquare = fit_measures["chisq"],
      DF = fit_measures["df"],
      P_value = fit_measures["pvalue"],
      RMSEA = fit_measures["rmsea"],
      RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
      RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
      RMSEA_P_value = fit_measures["rmsea.pvalue"],
      CFI = fit_measures["cfi"],
      TLI = fit_measures["tli"],
      SRMR = fit_measures["srmr"]
    )
    
    results_fit[[id]] <- fit 
    
    results_fit_indices <- rbind(results_fit_indices, fit_indices)
    
    # 提取路径估计
    parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
    
    # 将路径估计添加到 results_path_estimates 数据框
    path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                              "est", "se", "z", "pvalue", "std.lv", "std.all")]
    path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                       B = rep(Bvar,nrow(path_estimates))),
                            path_estimates)
    
    results_path_estimates <- rbind(results_path_estimates, path_estimates)
    
    print(paste('2-Wave Cross-Lagging Panel Model: ',Avar,' - ',Bvar, sep = ""))
    
  }
}

results_path_estimates[(results_path_estimates$lhs=='A_2y')&
                         (results_path_estimates$rhs=='B_bl')&
                         (results_path_estimates$pvalue<5^-6),]

results_path_estimates[(results_path_estimates$lhs=='B_2y')&
                         (results_path_estimates$rhs=='A_bl')&
                         (results_path_estimates$pvalue<5^-6),]


################################################################################
## 3 var 2 wave CLPM

# 初始化结果数据框
results_fit <- list()

results_fit_indices <- data.frame(
  A = character(),B = character(),C = character(),
  ChiSquare = numeric(),DF = numeric(),P_value = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_P_value = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric()
)

results_path_estimates <- data.frame(
  A = character(),B = character(),C = character(),
  lhs = character(),op = character(),rhs = character(),
  est = numeric(),se = numeric(),z = numeric(),pvalue = numeric(),std.lv = numeric(),std.all = numeric()
)

# 定义CLPM模型
model_syntax <- '
  # 自回归路径
  A_2y ~ A_bl
  B_2y ~ B_bl
  C_2y ~ C_bl
  
  # 跨滞路径
  A_2y ~ B_bl + C_bl
  B_2y ~ A_bl + C_bl
  C_2y ~ A_bl + B_bl
  
  # 协方差
  A_bl ~~ B_bl
  A_bl ~~ C_bl
  B_bl ~~ C_bl
  
  A_2y ~~ B_2y
  A_2y ~~ C_2y
  B_2y ~~ C_2y
'

id = 0

for (Bvar in c("dmdtifp1_pc1")) {
  
  for (Avar in names(abcd_brainPheno)[17:96]) {
    
    for (Cvar in  c('nihtbx_cryst_uncorrected',
                    'nihtbx_fluidcomp_uncorrected',
                    'nihtbx_totalcomp_uncorrected',
                    'cbcl_scr_syn_totprob_r')) {
      
      id = id + 1;
      
      data_residual <- abcd_brainPheno
      
      colnames(data_residual) <- sub(Bvar, "B", colnames(data_residual))
      colnames(data_residual) <- sub(Avar, "A", colnames(data_residual))
      colnames(data_residual) <- sub(Cvar, "C", colnames(data_residual))
      
      data_residual_bl <- data_residual[data_residual$eventname=='bl',]
      data_residual_2y <- data_residual[data_residual$eventname=='2y',]
      
      data_residual_bl$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                             smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                             (1|site_id_l:rel_family_id),data = data_residual_bl))
      data_residual_2y$A <- residuals(lmer(A~interview_age+Pubertal+sex+Race+
                                             smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                             (1|site_id_l:rel_family_id),data = data_residual_2y))
      
      data_residual_bl$B <- residuals(lmer(B~interview_age+Pubertal+sex+Race+
                                             smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                             (1|site_id_l:rel_family_id),data = data_residual_bl))
      data_residual_2y$B <- residuals(lmer(B~interview_age+Pubertal+sex+Race+
                                             smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                             (1|site_id_l:rel_family_id),data = data_residual_2y))
      
      data_residual <- rbind(data_residual_bl,data_residual_2y)
      
      duplicates <- data_residual %>%
        group_by(src_subject_id, eventname) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n > 1L)
      
      data_residual <- data_residual %>%
        group_by(src_subject_id, eventname) %>%
        summarise(across(c(A, B), mean, na.rm = TRUE), .groups = "drop")
      
      data_residual <- data_residual %>%
        mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
      
      # 为lavaan包准备数据                                                 
      data_wide <- data_residual %>%
        dplyr::select(src_subject_id, eventname, A, B, C) %>%
        pivot_wider(names_from = eventname, 
                    values_from = c(A, B, C))
      
      
      # 确保数据集不含缺失值的行
      required_vars <- c("A_bl", "A_2y", "B_bl", "B_2y","C_bl", "C_2y")
      data_wide_filtered <- data_wide %>%
        filter(complete.cases(select(., all_of(required_vars))))
      
      data_wide_standardized <- data_wide_filtered %>%
        mutate(across(all_of(required_vars), scale))
      
      fit <- lavaan::sem(model_syntax, data = data_wide_standardized)
      
      fit_measures <- fitMeasures(fit)
      parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
      
      fit_indices <- data.frame(
        A = Avar,B = Bvar,C = Cvar,
        ChiSquare = fit_measures["chisq"],
        DF = fit_measures["df"],
        P_value = fit_measures["pvalue"],
        RMSEA = fit_measures["rmsea"],
        RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
        RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
        RMSEA_P_value = fit_measures["rmsea.pvalue"],
        CFI = fit_measures["cfi"],
        TLI = fit_measures["tli"],
        SRMR = fit_measures["srmr"]
      )
      
      results_fit_indices <- rbind(results_fit_indices, fit_indices)
      
      # 提取路径估计
      parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
      
      # 将路径估计添加到 results_path_estimates 数据框
      path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                                "est", "se", "z", "pvalue", "std.lv", "std.all")]
      path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                         B = rep(Bvar,nrow(path_estimates)),
                                         C = rep(Cvar,nrow(path_estimates))),
                              path_estimates)
      
      results_path_estimates <- rbind(results_path_estimates, path_estimates)
      
      print(paste('2-Wave Cross-Lagging Panel Model: ',Avar,' - ',Bvar,' - ',Cvar, sep = ""))
    }
  }
}

results_path_estimates[(results_path_estimates$lhs=='A_2y')&
                         (results_path_estimates$rhs=='B_bl')&
                         (results_path_estimates$pvalue<0.05),]

results_path_estimates[(results_path_estimates$lhs=='B_2y')&
                         (results_path_estimates$rhs=='A_bl')&
                         (results_path_estimates$pvalue<5^-6),]

############################################################
## SEM

# 初始化结果数据框

abcd_brainPheno.bl <- abcd_brainPheno[abcd_brainPheno$eventname=="bl",]
abcd_brainPheno.2y <- abcd_brainPheno[abcd_brainPheno$eventname=="2y",]
abcd_brainPheno.bl$dmdtifp1_pc1 <- prcomp(apply(abcd_brainPheno.bl[, 97:100], 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)), center = TRUE, scale. = TRUE)$x[, 1]
abcd_brainPheno.bl$dmdtifp1_pc1[is.na(abcd_brainPheno.bl$dmdtifp1_19)] <- NA

abcd_brainPheno.2y$dmdtifp1_pc1 <- prcomp(apply(abcd_brainPheno.2y[, 97:100], 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)), center = TRUE, scale. = TRUE)$x[, 1]
abcd_brainPheno.2y$dmdtifp1_pc1[is.na(abcd_brainPheno.2y$dmdtifp1_19)] <- NA

abcd_brainPheno = rbind(abcd_brainPheno.bl,abcd_brainPheno.2y)

results_fit_indices <- data.frame(
  A = character(),B = character(),C = character(),D = character(),
  ChiSquare = numeric(),DF = numeric(),P_value = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_P_value = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric()
)

results_path_estimates <- data.frame(
  A = character(),B = character(),C = character(),D = character(),
  lhs = character(),op = character(),rhs = character(),
  est = numeric(),se = numeric(),z = numeric(),pvalue = numeric(),std.lv = numeric(),std.all = numeric()
)

sem_model <- '
  # 定义路径
  B_bl ~ A_bl
  C_2y ~ A_bl + B_bl
  D_2y ~ A_bl + B_bl + C_2y
'

id = 0
# "nihtbx_cryst_uncorrected"
for (Avar in c("latent_factor_ss_perinatal")){
  for (Cvar in c(names(abcd_brainPheno)[17:96])){
    for (Bvar in c("dmdtifp1_pc1")){
      for (Dvar in c("cbcl_scr_syn_totprob_r")){
        
        id = id + 1;
        
        sem_data <- abcd_brainPheno[,c("src_subject_id","eventname","interview_age","Pubertal",
                                       "sex","Race","smri_vol_scs_intracranialv","ehi_y_ss_scoreb",
                                       "site_id_l","rel_family_id",Avar, Bvar, Cvar, Dvar)]
        # 重命名列
        names(sem_data) <- c("src_subject_id","eventname","interview_age","Pubertal",
                             "sex","Race","smri_vol_scs_intracranialv","ehi_y_ss_scoreb",
                             "site_id_l","rel_family_id","A", "B", "C", "D")
        # 在使用 pivot_wider 之前先去除重复的组合
        sem_data <- sem_data %>%
          distinct(src_subject_id, eventname, .keep_all = TRUE)
        
        sem_data_residual <- sem_data
        sem_data_residual_bl <- na.omit(sem_data_residual[sem_data_residual$eventname=='bl',])
        sem_data_residual_2y <- na.omit(sem_data_residual[sem_data_residual$eventname=='2y',-11])
        
        sem_data_residual_bl$B <- residuals(lmer(B~interview_age+Pubertal+sex+Race+
                                                   smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_bl))
        sem_data_residual_2y$B <- residuals(lmer(B~interview_age+Pubertal+sex+Race+
                                                   smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_2y))
        
        sem_data_residual_bl$C <- residuals(lmer(C~interview_age+Pubertal+sex+Race+
                                                   ehi_y_ss_scoreb+
                                                   smri_vol_scs_intracranialv+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_bl))
        sem_data_residual_2y$C <- residuals(lmer(C~interview_age+Pubertal+sex+Race+
                                                   ehi_y_ss_scoreb+
                                                   smri_vol_scs_intracranialv+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_2y))
        
        sem_data_residual_2y$A <- NA
        
        sem_data_residual <- rbind(sem_data_residual_bl,sem_data_residual_2y)
        
        # 使用 pivot_wider 将数据转换为宽格式
        sem_data_wide <- sem_data_residual %>%
          select(src_subject_id, eventname, A, B, C, D) %>%
          pivot_wider(names_from = eventname, values_from = c(A, B, C, D))
        
        # 检查 eventname 列中的唯一值
        unique_eventnames <- unique(sem_data$eventname)
        
        # 确保数据集不含缺失值的行
        sem_data_wide <- sem_data_wide[,c("src_subject_id","A_bl", "B_bl", "C_2y", "D_2y")]
        sem_data_wide <- na.omit(as.data.frame(sem_data_wide))
        
        # 再次拟合模型
        fit <- sem(sem_model, data = sem_data_wide)
        
        # 查看模型拟合结果
        # summary(fit, fit.measures = TRUE, standardized = TRUE)
        fit_measures <- fitMeasures(fit)
        parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
        
        fit_indices <- data.frame(
          A = Avar,B = Bvar,C = Cvar,D = Cvar,
          ChiSquare = fit_measures["chisq"],
          DF = fit_measures["df"],
          P_value = fit_measures["pvalue"],
          RMSEA = fit_measures["rmsea"],
          RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
          RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
          RMSEA_P_value = fit_measures["rmsea.pvalue"],
          CFI = fit_measures["cfi"],
          TLI = fit_measures["tli"],
          SRMR = fit_measures["srmr"]
        )
        
        results_fit_indices <- rbind(results_fit_indices, fit_indices)
        
        # 提取路径估计
        parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
        
        # 将路径估计添加到 results_path_estimates 数据框
        path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                                  "est", "se", "z", "pvalue", "std.lv", "std.all")]
        path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                           B = rep(Bvar,nrow(path_estimates)),
                                           C = rep(Cvar,nrow(path_estimates)),
                                           D = rep(Dvar,nrow(path_estimates))),
                                path_estimates)
        
        results_path_estimates <- rbind(results_path_estimates, path_estimates)
        
        print(paste('2-Wave SEM: ',Avar,' - ',Bvar,' - ',Cvar,' - ',Dvar, sep = ""))
        
      }
    }
  }
}
rm(sem_data)
rm(sem_data_clean)
rm(sem_data_wide)
rm(sem_data_scaled)
rm(sem_data_residual_2y)
rm(sem_data_residual_bl)

results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                         (results_path_estimates$rhs=='C_bl')&
                         (results_path_estimates$pvalue<0.05),]

min(results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                             (results_path_estimates$rhs=='C_2y'),]$pvalue)

min(results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                             (results_path_estimates$rhs=='B_bl'),]$pvalue)

min(results_path_estimates[(results_path_estimates$lhs=='B_bl')&
                             (results_path_estimates$rhs=='A_bl'),]$pvalue)

ggseg.hemi(results_path_estimates[(results_path_estimates$lhs=='B_bl')&
                                    (results_path_estimates$rhs=='A_bl'),]$z[1:68],c(rep("CSA",34),rep("CT",34)),15)

ggseg.hemi(
  results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                           (results_path_estimates$rhs == 'A_bl'), ]$z[1:68] *
    (results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                              (results_path_estimates$rhs ==
                                 'A_bl'), ]$pvalue[1:68] < 0.005),
  c(rep("CSA", 34), rep("CT", 34)),7)

ggseg.hemi(
  results_path_estimates[(results_path_estimates$lhs == 'D_2y') &
                           (results_path_estimates$rhs == 'C_2y'), ]$z[1:68] *
    (results_path_estimates[(results_path_estimates$lhs == 'D_2y') &
                              (results_path_estimates$rhs ==
                                 'C_2y'), ]$pvalue[1:68] < 0.005),
  c(rep("CSA", 34), rep("CT", 34)),7)


ggseg.hemi(
  results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                           (results_path_estimates$rhs == 'B_bl'), ]$z[1:68] *
    (results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                              (results_path_estimates$rhs ==
                                 'B_bl'), ]$pvalue[1:68] < 0.005),
  c(rep("CSA", 34), rep("CT", 34)),15)

ggseg.hemi(results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                                    (results_path_estimates$rhs=='C_2y'),]$z[1:68],c(rep("CSA",34),rep("CT",34)),5)

ggseg.hemi(results_path_estimates[(results_path_estimates$lhs=='C_2y')&
                                    (results_path_estimates$rhs=='B_bl'),]$z[1:68],c(rep("CSA",34),rep("CT",34)),15)


############################################################
## SEM

# 初始化结果数据框

results_fit_indices <- data.frame(
  A = character(),B = character(),C = character(),D = character(),
  ChiSquare = numeric(),DF = numeric(),P_value = numeric(),
  RMSEA = numeric(),RMSEA_CI_Lower = numeric(),RMSEA_CI_Upper = numeric(),RMSEA_P_value = numeric(),
  CFI = numeric(),TLI = numeric(),SRMR = numeric()
)

results_path_estimates <- data.frame(
  A = character(),B = character(),C = character(),D = character(),
  lhs = character(),op = character(),rhs = character(),
  est = numeric(),se = numeric(),z = numeric(),pvalue = numeric(),std.lv = numeric(),std.all = numeric()
)

sem_model <- '
  # 定义路径
  B_bl ~ A_bl
  C_2y ~ A_bl + B_bl
  D_4y ~ A_bl + B_bl + C_2y
'

id = 0

for (Avar in c("latent_factor_ss_perinatal")){
  for (Cvar in c(names(abcd_brainPheno)[17:96])){
    for (Bvar in c("dmdtifp1_pc1")){
      for (Dvar in c("nihtbx_cryst_uncorrected")){
        
        id = id + 1;
        
        sem_data <- abcd_brainPheno[,c("src_subject_id","eventname","interview_age","Pubertal",
                                       "sex","Race","smri_vol_scs_intracranialv","ehi_y_ss_scoreb",
                                       "site_id_l","rel_family_id",Avar, Bvar, Cvar, Dvar)]
        # 重命名列
        names(sem_data) <- c("src_subject_id","eventname","interview_age","Pubertal",
                             "sex","Race","smri_vol_scs_intracranialv","ehi_y_ss_scoreb",
                             "site_id_l","rel_family_id","A", "B", "C", "D")
        # 在使用 pivot_wider 之前先去除重复的组合
        sem_data <- sem_data %>%
          distinct(src_subject_id, eventname, .keep_all = TRUE)
        
        sem_data_residual <- sem_data
        sem_data_residual_bl <- na.omit(sem_data_residual[sem_data_residual$eventname=='bl',])
        sem_data_residual_2y <- na.omit(sem_data_residual[sem_data_residual$eventname=='2y',-11])
        sem_data_residual_4y <- sem_data_residual[(sem_data$eventname=='4y')&(!is.na(sem_data$D)),]
        
        sem_data_residual_bl$B <- residuals(lmer(B~interview_age+Pubertal+sex+Race+
                                                   smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_bl))
        sem_data_residual_2y$B <- residuals(lmer(B~interview_age+Pubertal+sex+Race+
                                                   smri_vol_scs_intracranialv+ehi_y_ss_scoreb+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_2y))
        
        sem_data_residual_bl$C <- residuals(lmer(C~interview_age+Pubertal+sex+Race+
                                                   ehi_y_ss_scoreb+
                                                   smri_vol_scs_intracranialv+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_bl))
        sem_data_residual_2y$C <- residuals(lmer(C~interview_age+Pubertal+sex+Race+
                                                   ehi_y_ss_scoreb+
                                                   smri_vol_scs_intracranialv+
                                                   (1|site_id_l:rel_family_id),data = sem_data_residual_2y))
        
        sem_data_residual_2y$A <- NA
        
        sem_data_residual <- rbind(sem_data_residual_bl,sem_data_residual_2y,sem_data_residual_4y)
        
        # 使用 pivot_wider 将数据转换为宽格式
        sem_data_wide <- sem_data_residual %>%
          select(src_subject_id, eventname, A, B, C, D) %>%
          pivot_wider(names_from = eventname, values_from = c(A, B, C, D))
        
        # 检查 eventname 列中的唯一值
        unique_eventnames <- unique(sem_data$eventname)
        
        # 确保数据集不含缺失值的行
        sem_data_wide <- sem_data_wide[,c("src_subject_id","A_bl", "B_bl", "C_2y", "D_4y")]
        sem_data_wide <- na.omit(as.data.frame(sem_data_wide))
        
        # 再次拟合模型
        fit <- sem(sem_model, data = sem_data_wide)
        
        # 查看模型拟合结果
        # summary(fit, fit.measures = TRUE, standardized = TRUE)
        fit_measures <- fitMeasures(fit)
        parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
        
        fit_indices <- data.frame(
          A = Avar,B = Bvar,C = Cvar,D = Cvar,
          ChiSquare = fit_measures["chisq"],
          DF = fit_measures["df"],
          P_value = fit_measures["pvalue"],
          RMSEA = fit_measures["rmsea"],
          RMSEA_CI_Lower = fit_measures["rmsea.ci.lower"],
          RMSEA_CI_Upper = fit_measures["rmsea.ci.upper"],
          RMSEA_P_value = fit_measures["rmsea.pvalue"],
          CFI = fit_measures["cfi"],
          TLI = fit_measures["tli"],
          SRMR = fit_measures["srmr"]
        )
        
        results_fit_indices <- rbind(results_fit_indices, fit_indices)
        
        # 提取路径估计
        parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
        
        # 将路径估计添加到 results_path_estimates 数据框
        path_estimates <- parameter_estimates[, c("lhs", "op", "rhs", 
                                                  "est", "se", "z", "pvalue", "std.lv", "std.all")]
        path_estimates <- cbind(data.frame(A = rep(Avar,nrow(path_estimates)),
                                           B = rep(Bvar,nrow(path_estimates)),
                                           C = rep(Cvar,nrow(path_estimates)),
                                           D = rep(Dvar,nrow(path_estimates))),
                                path_estimates)
        
        results_path_estimates <- rbind(results_path_estimates, path_estimates)
        
        print(paste('2-Wave SEM: ',Avar,' - ',Bvar,' - ',Cvar,' - ',Dvar, sep = ""))
        
      }
    }
  }
}
rm(sem_data)
rm(sem_data_clean)
rm(sem_data_wide)
rm(sem_data_scaled)
rm(sem_data_residual_2y)
rm(sem_data_residual_bl)

results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                         (results_path_estimates$rhs=='C_bl')&
                         (results_path_estimates$pvalue<0.05),]

min(results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                             (results_path_estimates$rhs=='C_2y'),]$pvalue)

min(results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                             (results_path_estimates$rhs=='B_bl'),]$pvalue)

min(results_path_estimates[(results_path_estimates$lhs=='B_bl')&
                             (results_path_estimates$rhs=='A_bl'),]$pvalue)

ggseg.hemi(results_path_estimates[(results_path_estimates$lhs=='B_bl')&
                                    (results_path_estimates$rhs=='A_bl'),]$z[1:68],c(rep("CSA",34),rep("CT",34)),15)

ggseg.hemi(
  results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                           (results_path_estimates$rhs == 'A_bl'), ]$z[1:68] *
    (results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                              (results_path_estimates$rhs ==
                                 'A_bl'), ]$pvalue[1:68] < 0.005),
  c(rep("CSA", 34), rep("CT", 34)),7)

ggseg.hemi(
  results_path_estimates[(results_path_estimates$lhs == 'D_2y') &
                           (results_path_estimates$rhs == 'C_2y'), ]$z[1:68] *
    (results_path_estimates[(results_path_estimates$lhs == 'D_2y') &
                              (results_path_estimates$rhs ==
                                 'C_2y'), ]$pvalue[1:68] < 0.005),
  c(rep("CSA", 34), rep("CT", 34)),7)


ggseg.hemi(
  results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                           (results_path_estimates$rhs == 'B_bl'), ]$z[1:68] *
    (results_path_estimates[(results_path_estimates$lhs == 'C_2y') &
                              (results_path_estimates$rhs ==
                                 'B_bl'), ]$pvalue[1:68] < 0.005),
  c(rep("CSA", 34), rep("CT", 34)),15)

ggseg.hemi(results_path_estimates[(results_path_estimates$lhs=='D_2y')&
                                    (results_path_estimates$rhs=='C_2y'),]$z[1:68],c(rep("CSA",34),rep("CT",34)),5)

ggseg.hemi(results_path_estimates[(results_path_estimates$lhs=='C_2y')&
                                    (results_path_estimates$rhs=='B_bl'),]$z[1:68],c(rep("CSA",34),rep("CT",34)),15)
