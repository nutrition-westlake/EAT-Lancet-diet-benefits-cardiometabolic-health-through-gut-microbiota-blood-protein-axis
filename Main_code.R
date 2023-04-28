####Main codes for Adherence to the EAT-Lancet diet benefits cardiometabolic health through gut microbiota-blood protein axis (for questions, please contact Kui Deng, dengkui@westlake.edu.cn or dengkui_stat@163.com)
#1.Adherence to the EAT-Lancet diet was prospectively associated with specific circulating proteins
#2.Serum proteomic biomarkers of the EAT-Lancet diet were associated with cardiometabolic health
#3.Serum proteins mediated the association between the EAT-Lancet diet and cardiometabolic health
#4.Gut microbiota linked the associations between the EAT-Lancet diet and its related serum proteins


##All codes were executed using R


###1.Adherence to the EAT-Lancet diet was prospectively associated with specific circulating proteins
##LASSO model was used in the discovary cohort to select EAT-Lancet diet-related proteins
library(glmnet)
lambda_min_zong <- c()

for (i in 1:10){
  cat("######",i,"#######\n")
  set.seed(i)
  lasso_cv <- cv.glmnet(x=scale(as.matrix(data_protein123_protein_zong_follow_new[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3),])),y=factor(data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3)]),alpha=1,family="binomial")
  lambda_min_zong <- c(lambda_min_zong,lasso_cv$lambda.min)
}

lasso_model <- glmnet(x=scale(as.matrix(data_protein123_protein_zong_follow_new[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3),])),y=factor(data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3)]),alpha=1,family="binomial",lambda=min(lambda_min_zong))
lasso_model_beta <- as.matrix(lasso_model$beta)
lasso_model_beta_select <- lasso_model_beta[lasso_model_beta[,1]!=0,]
Protein_select <- names(lasso_model_beta_select)

##Multivariable linear regression  was used to further examine the prospective associations between the EAT-Lancet score (highest vs. lowest tertile) and selected serum proteins by LASSO, adjusted for potential confounders
#confounders included age, sex, BMI, smoking status, alcohol status, education, income, physical activity, total energy intake, T2D status, hypertension status, dyslipidemia medication usage, time interval, inner sequencing batch, and corresponding baseline protein abundance
Protein_longitudinal123_coef <- c()
for (i in Protein_select){
  cat("##########",i,"#############\n")
  data_fit_temp <- data.frame(protein_base=data_protein123_protein_zong_base_new[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3),i],protein_follow=data_protein123_protein_zong_follow_new[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3),i],data_protein123_phenotype_zong_base_new[data_protein123_phenotype_zong_base_new$Eat_Lancet_Score_Group %in% c(1,3),],stringsAsFactors = F)
  if (i %in% c("q9brj2","p01023")){
    data_fit_temp$protein_base <- log(data_fit_temp$protein_base)
    
    data_fit_temp$protein_follow <- log(data_fit_temp$protein_follow)
  }
  fit_temp <- lm(scale(protein_follow)~factor(Eat_Lancet_Score_Group)+scale(protein_base)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+Time_diff,data=data_fit_temp)
  fit_temp_summary <- summary(fit_temp)
  CI_temp <- confint(fit_temp)[2,]
  coef <- fit_temp_summary$coefficients[2,]
  coef <- c(coef,CI_temp)
  Protein_longitudinal123_coef <- rbind(Protein_longitudinal123_coef,coef)
}

Protein_longitudinal123_coef <- as.data.frame(Protein_longitudinal123_coef)
rownames(Protein_longitudinal123_coef) <- Protein_select

Protein_longitudinal123_coef$FDR <- p.adjust(Protein_longitudinal123_coef$`Pr(>|t|)`,method="fdr")

Protein_longitudinal123_coef_select <- Protein_longitudinal123_coef[Protein_longitudinal123_coef$FDR < 0.05,]

##Constructing Lancet-protein index 
Lancet_protein_index_phase123_follow <- scale(data_protein123_protein_zong_follow_new[,rownames(Protein_longitudinal123_coef_select)]) %*% sign(Protein_longitudinal123_coef_select$Estimate)
Lancet_protein_index_phase123_base <- scale(data_protein123_protein_zong_base_new[,rownames(Protein_longitudinal123_coef_select)]) %*% sign(Protein_longitudinal123_coef_select$Estimate)

#Examining the prospective association between the EAT-Lancet score (highest vs. lowest tertile) and Lancet-protein index using multivariable linear regression adjusting for the same confounders as above
data_fit_temp <- data.frame(Index_follow=Lancet_protein_index_phase123_follow,index_base=Lancet_protein_index_phase123_base,data_protein123_phenotype_zong_base_new,stringsAsFactors = F)
data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_Score_Group %in% c(1,3),]

fit_phase123_temp <- lm(scale(Index_follow)~factor(Eat_Lancet_Score_Group)+scale(index_base)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+Time_diff,data=data_fit_temp)

##Sensitivity analysis by linear mixed-effect model including all repeated-measured proteomics data and accounting for within-person correlation
#For selected proteins
Protein_longitudinal123_coef_repeat <- c()
for (i in rownames(Protein_longitudinal123_coef_select)){
  cat("##########",i,"#############\n")
  
  data_fit_temp_follow <- data.frame(protein_follow=data_protein123_phenotype_zong_repeat_follow[,i],data_protein123_phenotype_zong_repeat_follow,stringsAsFactors = F)
  data_fit_temp_base <- data.frame(protein_base=data_protein123_phenotype_zong_repeat_base[,i],id=data_protein123_phenotype_zong_repeat_base$id,stringsAsFactors = F)
  data_fit_temp <- merge(data_fit_temp_base,data_fit_temp_follow,by="id")
  data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_Score_Group %in% c(1,3),]
  
  
  if (i %in% c("q9brj2","p01023")){
    data_fit_temp$protein_base <- log(data_fit_temp$protein_base)
    
    data_fit_temp$protein_follow <- log(data_fit_temp$protein_follow)
  }
  
  library(lme4)
  fit_temp <- lmer(scale(protein_follow)~factor(Eat_Lancet_Score_Group)+scale(protein_base)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+(1|id),data=data_fit_temp)
  
  fit_temp_summary <- summary(fit_temp)
  CI_temp <- confint(fit_temp)[4,]
  
  coef <- fit_temp_summary$coefficients[2,]
  coef <- c(coef,CI_temp)
  Protein_longitudinal123_coef_repeat <- rbind(Protein_longitudinal123_coef_repeat,coef)
}

Protein_longitudinal123_coef_repeat <- as.data.frame(Protein_longitudinal123_coef_repeat)
rownames(Protein_longitudinal123_coef_repeat) <- rownames(Protein_longitudinal123_coef_select)

#For Lancet-protein index 
data_protein123_phenotype_zong_repeat$Lancet_protein_index_phase123 <- scale(data_protein123_phenotype_zong_repeat[,rownames(Protein_longitudinal123_coef_select)]) %*% sign(Protein_longitudinal123_coef_select$Estimate)

library(lme4)
fit_phase123_temp_repeat <- lmer(scale(Index_follow)~factor(Eat_Lancet_Score_Group)+scale(Index_base)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+(1|id),data=data_fit_temp)

##All above analyses were performed in the validation cohort (inner sequencing batch was not adjusted)

##Meta the results from the discovery and validation cohorts
#For example
data_discovery <- Protein_longitudinal123_coef_select
data_validation <- Protein_coef_validation

library(metafor)
Meta_discovery_validation <- c()
for (i in 1:nrow(data_discovery)){
  cat("#########",i,"###########\n")
  data_discovery_temp <- data_discovery[i,c("Estimate","Std. Error")]
  data_validation_temp  <-data_validation[i,c("Estimate","Std. Error")]
  
  data_test_temp <- rbind(data_discovery_temp,data_validation_temp)
  meta_result_temp <- rma(yi=Estimate,sei=`Std. Error`,data=data_test_temp,method="FE")
  result_temp <- data.frame(protein=rownames(data_discovery)[i],B_discovery=data_discovery$Estimate[i],SE_discovery=data_discovery$`Std. Error`[i],p_discovery=data_discovery$`Pr(>|t|)`[i],
                            B_validation=data_validation$Estimate[i],SE_validation=data_validation$`Std. Error`[i],p_validation=data_validation$`Pr(>|t|)`[i],
                            I2=meta_result_temp$I2,heterogeneity=meta_result_temp$QEp,B_meta=meta_result_temp$beta,SE_meta=meta_result_temp$se,CI_lower_meta=meta_result_temp$ci.lb,CI_upper_meta=meta_result_temp$ci.ub,p_meta=meta_result_temp$pval,stringsAsFactors=F
  )
  Meta_discovery_validation <- rbind(Meta_discovery_validation,result_temp)
}

Meta_discovery_validation$FDR_meta <- p.adjust(Meta_discovery_validation$p_meta,method="fdr")


###2.Serum proteomic biomarkers of the EAT-Lancet diet were associated with cardiometabolic health
##Linear mixed-effect model was used to examine the associations between EAT-Lancet diet-related proteins and cardiometabolic risk factors 
#Discovery cohort
traits <- c("Glu_new","Ins_new","HOMA_IR_new","HbA1c_new","TC_new","TG_new","HDL_new","LDL_new","SBP_new","DBP_new","BMI_new")

Protein_selected <- rownames(Protein_longitudinal123_coef_select)

library(lme4)
library(lmerTest)
Protein_trait_result_discovery <- c()
for (i in Protein_selected){
  data_result_temp <- c()
  for (j in traits){
    data_temp <- data.frame(traits=data_protein_phenotype[,j],
                            protein=data_protein_phenotype[,i],
                            data_protein_phenotype,stringsAsFactors = F)
    if (j %in% c("BMI_new")){
      fit_temp <- lmer(scale(traits)~scale(protein)+age+sex+smoke+alc+factor(edu3)+factor(income4)+MET+energy+factor(phase)+(1|id),data=data_temp[data_temp$batch==1,])
    } else {
      fit_temp <- lmer(scale(traits)~scale(protein)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+factor(phase)+(1|id),data=data_temp[data_temp$batch==1,])
    }
    
    fit_temp_summary <- summary(fit_temp)
    
    
    coef_temp <- fit_temp_summary$coefficients[2,]
    CI_temp <- confint(fit_temp)[4,]
    
    coef_temp <- c(protein=i,traits=j,coef_temp,CI_temp)
    data_result_temp <- rbind(data_result_temp,coef_temp)
    
  }
  data_result_temp <- as.data.frame(data_result_temp)
  data_result_temp[,-c(1:2)] <- apply(data_result_temp[,-c(1:2)],2,as.numeric)
  data_result_temp$FDR <- p.adjust(data_result_temp$`Pr(>|t|)`,method="fdr")
  
  Protein_trait_result_discovery <- rbind(Protein_trait_result_discovery,data_result_temp)
  
  
}

#Validation cohort
Protein_trait_result_validation <- c()
for (i in Protein_selected){
  data_result_temp <- c()
  for (j in traits){
    data_temp <- data.frame(traits=data_protein_phenotype[,j],
                            protein=data_protein_phenotype[,i],
                            data_protein_phenotype,stringsAsFactors = F)
    
    if (j %in% c("BMI_new")){
      fit_temp <- lmer(scale(traits)~scale(protein)+age+sex+smoke+alc+factor(edu3)+factor(income4)+MET+energy+(1|id),data=data_temp[data_temp$batch==2,])
    } else {
      fit_temp <- lmer(scale(traits)~scale(protein)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+(1|id),data=data_temp[data_temp$batch==2,])
    }
    fit_temp_summary <- summary(fit_temp)
    coef_temp <- fit_temp_summary$coefficients[2,]
    CI_temp <- confint(fit_temp)[4,]
    
    coef_temp <- c(protein=i,traits=j,coef_temp,CI_temp)
    data_result_temp <- rbind(data_result_temp,coef_temp)
    
  }
  data_result_temp <- as.data.frame(data_result_temp)
  data_result_temp[,-c(1:2)] <- apply(data_result_temp[,-c(1:2)],2,as.numeric)
  data_result_temp$FDR <- p.adjust(data_result_temp$`Pr(>|t|)`,method="fdr")
  
  Protein_trait_result_validation <- rbind(Protein_trait_result_validation,data_result_temp)
}

#Meta the results from the discovery and validation cohorts
data_discovery <- Protein_trait_result_discovery
data_validation <- Protein_trait_result_validation

library(metafor)
Meta_protein_traits_discovery_validation <- c()
for (i in 1:nrow(data_discovery)){
  cat("#########",i,"###########\n")
  data_discovery_temp <- data_discovery[i,c("Estimate","Std. Error")]
  data_validation_temp  <-data_validation[i,c("Estimate","Std. Error")]
  
  data_test_temp <- rbind(data_discovery_temp,data_validation_temp)
  meta_result_temp <- rma(yi=Estimate,sei=`Std. Error`,data=data_test_temp,method="FE")
  result_temp <- data.frame(protein=data_discovery$protein[i],traits=data_discovery$traits[i],
                            B_discovery=data_discovery$Estimate[i],CI_lower_discovery=data_discovery$`2.5 %`[i],CI_upper_discovery=data_discovery$`97.5 %`[i],p_discovery=data_discovery$FDR_zong[i],
                            B_validation=data_validation$Estimate[i],CI_lower_validation=data_validation$`2.5 %`[i],CI_upper_validation=data_validation$`97.5 %`[i],p_validation=data_validation$`Pr(>|t|)`[i],
                            B_meta=meta_result_temp$beta,CI_lower_meta=meta_result_temp$ci.lb,CI_upper_meta=meta_result_temp$ci.ub,p_meta=meta_result_temp$pval,heterogeneity=meta_result_temp$QEp,stringsAsFactors=F
  )
  Meta_protein_traits_discovery_validation <- rbind(Meta_protein_traits_discovery_validation,result_temp)
}

Meta_protein_traits_discovery_validation$FDR_meta <- p.adjust(Meta_protein_traits_discovery_validation$p_meta,method="fdr")

##Similar codes were used for Lancet-protein index-cardiometabolic risk factor associations

##Logistic model was used to examine prospective associations between baseline EAT-Lancet diet-related proteins and incident cardiometabolic diseases 
#Taking protein-type 2 diabetas association as an example
data_phenotype_disease_NL_protein <- data_phenotype_disease[data_protein_phenotype_NL$id,]

#For discovery
Result_protein_T2D_discovery <- c()
for (protein_temp in rownames(Protein_longitudinal123_coef_select)){
  data_fit_temp <- data.frame(protein=data_protein_phenotype_NL[,protein_temp],
                              batch=data_protein_phenotype_NL$batch,
                              phase=data_protein_phenotype_NL$phase,
                              data_phenotype_disease_NL_protein,stringsAsFactors = F
  )
  data_fit_temp <- data_fit_temp[-which(data_fit_temp$DM_diagnosis_v1==1),]
  T2D_model_protein <- glm(DM_diagnosis_follow~scale(protein)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+factor(phase),data=data_fit_temp[data_fit_temp$batch==1,],family=binomial(link="logit"))
  T2D_model_protein_summary <- summary(T2D_model_protein)
  OR_summary <- exp(c(T2D_model_protein_summary$coefficients[2,1],confint(T2D_model_protein)[2,]))
  
  names(OR_summary) <- c("OR","OR_lower","OR_upper")
  
  coef_temp <- c(T2D_model_protein_summary$coefficients[2,],OR_summary)
  Result_protein_T2D_discovery <- rbind(Result_protein_T2D_discovery,coef_temp)
  
}

rownames(Result_protein_T2D_discovery) <-rownames(Protein_longitudinal123_coef_select)
Result_protein_T2D_discovery <- as.data.frame(Result_protein_T2D_discovery)

#For validation
Result_protein_T2D_validation <- c()
for (protein_temp in rownames(Protein_longitudinal123_coef_select)){
  data_fit_temp <- data.frame(protein=data_protein_phenotype_NL[,protein_temp],
                              batch=data_protein_phenotype_NL$batch,
                              data_phenotype_disease_NL_protein,stringsAsFactors = F
  )
  data_fit_temp <- data_fit_temp[-which(data_fit_temp$DM_diagnosis_v1==1),]
  T2D_model_protein <- glm(DM_diagnosis_follow~scale(protein)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy,data=data_fit_temp[data_fit_temp$batch==2,],family=binomial(link="logit"))
  T2D_model_protein_summary <- summary(T2D_model_protein)
  OR_summary <- exp(c(T2D_model_protein_summary$coefficients[2,1],confint(T2D_model_protein)[2,]))
  
  names(OR_summary) <- c("OR","OR_lower","OR_upper")
  
  coef_temp <- c(T2D_model_protein_summary$coefficients[2,],OR_summary)
  Result_protein_T2D_validation <- rbind(Result_protein_T2D_validation,coef_temp)
  
}

rownames(Result_protein_T2D_validation) <-rownames(Protein_longitudinal123_coef_select)
Result_protein_T2D_validation <- as.data.frame(Result_protein_T2D_validation)

#Meta the results from the discovery and validation cohorts
data_discovery <- Result_protein_T2D_discovery
data_validation <- Result_protein_T2D_validation

library(metafor)
Meta_protein_T2D <- c()
for (i in 1:nrow(data_discovery)){
  cat("#########",i,"###########\n")
  data_discovery_temp <- data_discovery[i,c("Estimate","Std. Error")]
  data_validation_temp  <-data_validation[i,c("Estimate","Std. Error")]
  
  data_test_temp <- rbind(data_discovery_temp,data_validation_temp)
  meta_result_temp <- rma(yi=Estimate,sei=`Std. Error`,data=data_test_temp,method="FE")
  result_temp <- data.frame(Protein=rownames(data_discovery)[[i]],OR_discovery=data_discovery$OR[i],OR_lower_discovery=data_discovery$OR_lower[i],OR_upper_discovery=data_discovery$OR_upper[i],SE_discovery=data_discovery$`Std. Error`[i],p_discovery=data_discovery$`Pr(>|z|)`[i],
                            OR_validation=data_validation$OR[i],OR_lower_validation=data_validation$OR_lower[i],OR_upper_validation=data_validation$OR_upper[i],SE_validation=data_validation$`Std. Error`[i],p_validation=data_validation$`Pr(>|z|)`[i],
                            heterogeneity=meta_result_temp$QEp,B_meta=meta_result_temp$beta,SE_meta=meta_result_temp$se,CI_lower_meta=meta_result_temp$ci.lb,CI_upper_meta=meta_result_temp$ci.ub,p_meta=meta_result_temp$pval,stringsAsFactors=F
  )
  Meta_protein_T2D <- rbind(Meta_protein_T2D,result_temp)
}

Meta_protein_T2D$FDR_meta <- p.adjust(Meta_protein_T2D$p_meta,method="fdr")

##Similar codes were used for other cardiometabolic diseases (hypertension, dyslipidemia, and metabolic syndrome) and Lancet-protein index

###3.Serum proteins mediated the association between the EAT-Lancet diet and cardiometabolic health
##For the mediation: Eat-Lancet Diet to proteins to cardiometabolic risk factors
#For discovery
Protein_traits_mediation_result_discovery <- c()
for (i in c(1:nrow(Meta_protein_traits_discovery_validation_select))){
  cat("#############",i,"#############\n")
  trait_temp <- Meta_protein_traits_discovery_validation_select$traits[i]
  protein_temp <- Meta_protein_traits_discovery_validation_select$protein[i]
  
  data_trait_v3 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v3",sep="")]
  data_trait_v4 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v4",sep="")]
  
  
  data_protein_base <- data_protein_phenotype_NL_F2_NL[,protein_temp]
  data_protein_F2 <- data_protein_phenotype_NL_F2_F2[,protein_temp]
  
  
  if (protein_temp %in% c("q9brj2","p01023")){
    data_protein_base <- log(data_protein_base)
    data_protein_F2 <- log(data_protein_F2)
  }
  
  
  data_fit_temp <- data.frame(Protein_base=scale(data_protein_base),
                              Protein_F2 = scale(data_protein_F2),
                              Eat_Lancet_group=data_protein_phenotype_NL_F2_NL$Eat_Lancet_Score_Group,
                              trait_v3=scale(data_trait_v3),trait_v4=scale(data_trait_v4),data_protein_phenotype_NL_F2_NL,stringsAsFactors = F)
  
  data_fit_temp <- data_fit_temp[data_fit_temp$batch==1,]
  data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_group %in% c(1,3),]
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==1] <- 0
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==3] <- 1
  data_fit_temp <- data_fit_temp[!(is.na(data_fit_temp$trait_v3) | is.na(data_fit_temp$trait_v4)),]
  
    if (trait_temp %in% c("BMI")){
    Med_fit_discovery <- lm(Protein_F2~Eat_Lancet_group+Protein_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+Time_diff,data=data_fit_temp)
    Out_fit_discovery <- lm(trait_v4~Protein_F2+Eat_Lancet_group+trait_v3+age+sex+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3+factor(phase),data=data_fit_temp)
    
    
  } else {
    Med_fit_discovery <- lm(Protein_F2~Eat_Lancet_group+Protein_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+Time_diff,data=data_fit_temp)
    Out_fit_discovery <- lm(trait_v4~Protein_F2+Eat_Lancet_group+trait_v3+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3+factor(phase),data=data_fit_temp)
  }
  
  library(mediation)
  set.seed(1234)
  Med_discovery <- mediate(Med_fit_discovery, Out_fit_discovery, treat = "Eat_Lancet_group",mediator = "Protein_F2")
  
  Result_mediation_temp <- c(P_mediation=Med_discovery$d.avg.p,Prop_mediation=Med_discovery$n.avg)
  Name_temp <- c(X="Eat_Lancet_Score_Group",Mediator=protein_temp,Y=trait_temp)
  
  Result_mediation_temp <- c(Name_temp,Result_mediation_temp)
  Protein_traits_mediation_result_discovery <- rbind(Protein_traits_mediation_result_discovery,Result_mediation_temp)
}

Protein_traits_mediation_result_discovery <- as.data.frame(Protein_traits_mediation_result_discovery,stringsAsFactors = F)
Protein_traits_mediation_result_discovery[,-c(1:3)] <- apply(Protein_traits_mediation_result_discovery[,-c(1:3)],2,as.numeric)
Protein_traits_mediation_result_discovery_select <- Protein_traits_mediation_result_discovery[Protein_traits_mediation_result_discovery$P_mediation<0.05,]

#For validation
Protein_traits_mediation_result_validation <- c()
for (i in c(1:nrow(Protein_traits_mediation_result_discovery_select))){
  cat("#############",i,"#############\n")
  trait_temp <- Protein_traits_mediation_result_discovery_select$Y[i]
  protein_temp <- Protein_traits_mediation_result_discovery_select$Mediator[i]
  
  data_trait_v3 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v3",sep="")]
  data_trait_v4 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v4",sep="")]
  
  
  data_protein_base <- data_protein_phenotype_NL_F2_NL[,protein_temp]
  data_protein_F2 <- data_protein_phenotype_NL_F2_F2[,protein_temp]
  
  
  if (protein_temp %in% c("q9brj2","p01023")){
    data_protein_base <- log(data_protein_base)
    data_protein_F2 <- log(data_protein_F2)
  }
  
  
  data_fit_temp <- data.frame(Protein_base=scale(data_protein_base),
                              Protein_F2 = scale(data_protein_F2),
                              Eat_Lancet_group=data_protein_phenotype_NL_F2_NL$Eat_Lancet_Score_Group,
                              trait_v3=scale(data_trait_v3),trait_v4=scale(data_trait_v4),data_protein_phenotype_NL_F2_NL,stringsAsFactors = F)
  
  data_fit_temp <- data_fit_temp[data_fit_temp$batch==2,]
  data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_group %in% c(1,3),]
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==1] <- 0
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==3] <- 1
  data_fit_temp <- data_fit_temp[!(is.na(data_fit_temp$trait_v3) | is.na(data_fit_temp$trait_v4)),]
  
  if (trait_temp %in% c("BMI")){
    Med_fit_validation <- lm(Protein_F2~Eat_Lancet_group+Protein_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+Time_diff,data=data_fit_temp)
    Out_fit_validation <- lm(trait_v4~Protein_F2+Eat_Lancet_group+trait_v3+age+sex+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3,data=data_fit_temp)
    
    
  } else {
    Med_fit_validation <- lm(Protein_F2~Eat_Lancet_group+Protein_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+Time_diff,data=data_fit_temp)
    Out_fit_validation <- lm(trait_v4~Protein_F2+Eat_Lancet_group+trait_v3+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3,data=data_fit_temp)
  }
  
  library(mediation)
  set.seed(1234)
  Med_validation <- mediate(Med_fit_validation, Out_fit_validation, treat = "Eat_Lancet_group",mediator = "Protein_F2")
  
  Result_mediation_temp <- c(P_mediation=Med_validation$d.avg.p,Prop_mediation=Med_validation$n.avg)
  Name_temp <- c(X="Eat_Lancet_Score_Group",Mediator=protein_temp,Y=trait_temp)
  
  Result_mediation_temp <- c(Name_temp,Result_mediation_temp)
  Protein_traits_mediation_result_validation <- rbind(Protein_traits_mediation_result_validation,Result_mediation_temp)
  
}

Protein_traits_mediation_result_validation <- as.data.frame(Protein_traits_mediation_result_validation)


##For the mediation: Eat-Lancet Diet to Lancet-protein index to cardiometabolic risk factors
#For discovery
LancetProtein_traits_mediation_result_discovery <- c()
for (i in c(1:nrow(Meta_lancetScore_traits_discovery_validation_select))){
  cat("#############",i,"#############\n")
  trait_temp <- Meta_lancetScore_traits_discovery_validation_select$traits[i]
  
  data_trait_v3 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v3",sep="")]
  data_trait_v4 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v4",sep="")]
  
  data_fit_temp <- data.frame(Index_base=scale(data_protein_phenotype_NL_F2_NL$Eat_Lancet_protein_index_base),
                              Index_F2 = scale(data_protein_phenotype_NL_F2_NL$Eat_Lancet_protein_index_F2),
                              Eat_Lancet_group=data_protein_phenotype_NL_F2_NL$Eat_Lancet_Score_Group,
                              trait_v3=scale(data_trait_v3),trait_v4=scale(data_trait_v4),data_protein_phenotype_NL_F2_NL,stringsAsFactors = F)
  
  data_fit_temp <- data_fit_temp[data_fit_temp$batch==1,]
  data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_group %in% c(1,3),]
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==1] <- 0
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==3] <- 1
  data_fit_temp <- data_fit_temp[!(is.na(data_fit_temp$trait_v3) | is.na(data_fit_temp$trait_v4)),]
  
  if (trait_temp %in% c("BMI")){
    Med_fit_discovery <- lm(Index_F2~Eat_Lancet_group+Index_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+Time_diff,data=data_fit_temp)
    Out_fit_discovery <- lm(trait_v4~Index_F2+Eat_Lancet_group+trait_v3+age+sex+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3+factor(phase),data=data_fit_temp)
    
    
  } else {
    Med_fit_discovery <- lm(Index_F2~Eat_Lancet_group+Index_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase)+Time_diff,data=data_fit_temp)
    Out_fit_discovery <- lm(trait_v4~Index_F2+Eat_Lancet_group+trait_v3+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3+factor(phase),data=data_fit_temp)
  }
  library(mediation)
  set.seed(1234)
  Med_discovery <- mediate(Med_fit_discovery, Out_fit_discovery, treat = "Eat_Lancet_group",mediator = "Index_F2")
  
  Result_mediation_temp <- c(P_mediation=Med_discovery$d.avg.p,Prop_mediation=Med_discovery$n.avg)
  Name_temp <- c(X="Eat_Lancet_Score_Group",Mediator="Lancet_Protein_Index",Y=trait_temp)
  
  Result_mediation_temp <- c(Name_temp,Result_mediation_temp)
  LancetProtein_traits_mediation_result_discovery <- rbind(LancetProtein_traits_mediation_result_discovery,Result_mediation_temp)
  
}
LancetProtein_traits_mediation_result_discovery <- as.data.frame(LancetProtein_traits_mediation_result_discovery)
LancetProtein_traits_mediation_result_discovery[,-c(1:3)] <- apply(LancetProtein_traits_mediation_result_discovery[,-c(1:3)],2,as.numeric)

LancetProtein_traits_mediation_result_discovery_select <- LancetProtein_traits_mediation_result_discovery[LancetProtein_traits_mediation_result_discovery$P_mediation < 0.05,]


#For Validation
LancetProtein_traits_mediation_result_validation <- c()

for (i in c(1:nrow(LancetProtein_traits_mediation_result_discovery_select))){
  cat("#############",i,"#############\n")
  trait_temp <- LancetProtein_traits_mediation_result_discovery_select$Y[i]
  
  data_trait_v3 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v3",sep="")]
  data_trait_v4 <- data_protein_phenotype_NL_F2_NL[,paste(trait_temp,"_v4",sep="")]
  
  data_fit_temp <- data.frame(Index_base=scale(data_protein_phenotype_NL_F2_NL$Eat_Lancet_protein_index_base),
                              Index_F2 = scale(data_protein_phenotype_NL_F2_NL$Eat_Lancet_protein_index_F2),
                              Eat_Lancet_group=data_protein_phenotype_NL_F2_NL$Eat_Lancet_Score_Group,
                              trait_v3=scale(data_trait_v3),trait_v4=scale(data_trait_v4),data_protein_phenotype_NL_F2_NL,stringsAsFactors = F)
  
  data_fit_temp <- data_fit_temp[data_fit_temp$batch==2,]
  data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_group %in% c(1,3),]
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==1] <- 0
  data_fit_temp$Eat_Lancet_group[data_fit_temp$Eat_Lancet_group==3] <- 1
  data_fit_temp <- data_fit_temp[!(is.na(data_fit_temp$trait_v3) | is.na(data_fit_temp$trait_v4)),]
  
  if (trait_temp %in% c("BMI")){
    Med_fit_validation <- lm(Index_F2~Eat_Lancet_group+Index_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+Time_diff,data=data_fit_temp)
    Out_fit_validation <- lm(trait_v4~Index_F2+Eat_Lancet_group+trait_v3+age+sex+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3,data=data_fit_temp)
    
    
  } else {
    Med_fit_validation <- lm(Index_F2~Eat_Lancet_group+Index_base+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+Time_diff,data=data_fit_temp)
    Out_fit_validation <- lm(trait_v4~Index_F2+Eat_Lancet_group+trait_v3+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+Time_diff_F2_F3,data=data_fit_temp)
  }
  library(mediation)
  set.seed(1234)
  Med_validation <- mediate(Med_fit_validation, Out_fit_validation, treat = "Eat_Lancet_group",mediator = "Index_F2")
  
  
  Result_mediation_temp <- c(P_mediation=Med_validation$d.avg.p,Prop_mediation=Med_validation$n.avg)
  Name_temp <- c(X="Eat_Lancet_Score_Group",Mediator="Lancet_Protein_Index",Y=trait_temp)
  
  Result_mediation_temp <- c(Name_temp,Result_mediation_temp)
  LancetProtein_traits_mediation_result_validation <- rbind(LancetProtein_traits_mediation_result_validation,Result_mediation_temp)
  
}

LancetProtein_traits_mediation_result_validation[,-c(1:3)] <- apply(LancetProtein_traits_mediation_result_validation[,-c(1:3)],2,as.numeric)


###4.Gut microbiota linked the associations between the EAT-Lancet diet and its related serum proteins
##LASSO model was used to select EAT-Lancet diet-related microbes
library(glmnet)
lambda_min_zong_microbiome <- c()
for (i in 1:10){
  cat("######",i,"#######\n")
  set.seed(i)
  lasso_cv_microbiome <- cv.glmnet(x=scale(as.matrix(data_microbiome_final_lasso)),y=factor(data_microbiome_phenotype_final_lasso$Eat_Lancet_Score_Group),alpha=1,family="binomial")
  lambda_min_zong_microbiome <- c(lambda_min_zong_microbiome,lasso_cv_microbiome$lambda.min)
}

lasso_microbiome_model <- glmnet(x=scale(as.matrix(data_microbiome_final_lasso)),y=factor(data_microbiome_phenotype_final_lasso$Eat_Lancet_Score_Group),alpha=1,family="binomial",lambda=min(lambda_min_zong_microbiome))
lasso_microbiome_model_beta <- as.matrix(lasso_microbiome_model$beta)
lasso_microbiome_model_beta_select <- lasso_microbiome_model_beta[lasso_microbiome_model_beta[,1]!=0,]
microbiome_select <- names(lasso_microbiome_model_beta_select)

##Multivariable linear regression was used to further examine the prospective associations between the EAT-Lancet score (highest vs. lowest tertile) and selected gut microbes by LASSO, adjusted for potential confounders
#confounders included age, sex, BMI, smoking status, alcohol status, education, income, physical activity, total energy intake, T2D status, hypertension status, dyslipidemia medication usage, time interval, and Bristol stool score
microbiome_linear_coef <- c()
for (i in microbiome_select){
  cat("#######",i,"#########\n")
  data_fit_temp <- data.frame(microbiome=data_microbiome_final_CLR[,i],data_microbiome_phenotype_final,stringsAsFactors = F)
  fit_temp <- lm(scale(microbiome)~factor(Eat_Lancet_Score_Group)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+bristol_scale+Time_diff,data=data_fit_temp)
  fit_temp_summary <- summary(fit_temp)
  CI_temp <- confint(fit_temp)[3,]
  coef_temp <- c(fit_temp_summary$coefficients[3,],CI_temp)
  microbiome_linear_coef <- rbind(microbiome_linear_coef,coef_temp)
}

microbiome_linear_coef <- as.data.frame(microbiome_linear_coef)
rownames(microbiome_linear_coef) <- microbiome_select
microbiome_linear_coef$FDR <- p.adjust(microbiome_linear_coef$`Pr(>|t|)`,method="fdr")

##Sensitivity analysis by linear mixed-effect model
microbiome_linear_repeat_coef <- c()
for (i in microbiome_select){
  cat("#######",i,"#########\n")
  library(lme4)
  data_fit_temp <- data.frame(microbiome=data_microbiome_final_repeat_CLR[,i],data_microbiome_phenotype_final_repeat,stringsAsFactors = F)
  data_fit_temp <- data_fit_temp[data_fit_temp$Eat_Lancet_Score_Group %in% c(1,3),]
  fit_temp <- lmer(scale(microbiome)~factor(Eat_Lancet_Score_Group)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+bristol_scale+(1|id),data=data_fit_temp)
  fit_temp_summary <- summary(fit_temp)
  CI_temp <- confint(fit_temp)[4,]
  coef_temp <- c(fit_temp_summary$coefficients[2,],CI_temp)
  microbiome_linear_repeat_coef <- rbind(microbiome_linear_repeat_coef,coef_temp)
}

microbiome_linear_repeat_coef <- as.data.frame(microbiome_linear_repeat_coef)
rownames(microbiome_linear_repeat_coef) <- microbiome_select
microbiome_linear_repeat_coef$FDR <- p.adjust(microbiome_linear_repeat_coef$`Pr(>|t|)`,method="fdr")

##The associations between EAT-Lancet diet-related microbes and proteins
Microbe_identified <- rownames(microbiome_linear_coef)
Protein_identified <- rownames(Protein_longitudinal123_coef_select)

Microbe_protein_association_result_new <- c()
for (i in Microbe_identified){
  data_result_temp <- c()
  for (j in Protein_identified){
    cat("##########",i,"#########",j,"###\n")
    data_fit_temp <- data.frame(Protein=data_protein_zong_new[,j],
                                Microbe=data_microbiome_zong_new_CLR[,i],
                                data_phenotype_protein_microbiome_new,stringsAsFactors = F)
    
    if (i %in% c("q9brj2","p01023")){
      data_fit_temp$Protein <- log(data_fit_temp$Protein)
      
    }
    
    fit_model <- lm(scale(Protein)~scale(Microbe)+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase),data=data_fit_temp)
    fit_model_summary <- summary(fit_model)
    coef <- fit_model_summary$coefficients[2,]
    CI_temp <- confint(fit_model)[2,]
    coef <- c(Microbe=i,Protein=j,coef,CI_temp)
    data_result_temp <- rbind(data_result_temp,coef)
    
  }
  data_result_temp <- as.data.frame(data_result_temp)
  data_result_temp[,-c(1:2)] <- apply(data_result_temp[,-c(1:2)],2,as.numeric)
  
  Microbe_protein_association_result_new <- rbind(Microbe_protein_association_result_new,data_result_temp)
  
}
Microbe_protein_association_result_new$FDR <- p.adjust(Microbe_protein_association_result_new$`Pr(>|t|)`,method="fdr")
Microbe_protein_association_result_new_select <- Microbe_protein_association_result_new[Microbe_protein_association_result_new$FDR < 0.05,]

##Mediation analysis for EAT-Lancet diet-Rothia mucilaginosa-A2M association
# p01023 --- A2M  s49 --- Rothia mucilaginosa 
data_fit_temp_new <- data.frame(Protein=data_protein_zong_new[,"p01023"],
                                 Microbe=data_microbiome_zong_new_CLR[,"s49"],
                                 data_phenotype_protein_microbiome_new,stringsAsFactors = F)

data_fit_temp_new$Protein <- log(data_fit_temp_new$Protein)
  

data_fit_temp_new$Protein <- scale(data_fit_temp_new$Protein)
data_fit_temp_new$Microbe <- scale(data_fit_temp_new$Microbe)


library(mediation)
data_fit_temp_new <- data_fit_temp_new[data_fit_temp_new$Eat_Lancet_Score_Group %in% c(1,3),]

data_fit_temp_new$Eat_Lancet_Score_Group[data_fit_temp_new$Eat_Lancet_Score_Group==1] <- 0
data_fit_temp_new$Eat_Lancet_Score_Group[data_fit_temp_new$Eat_Lancet_Score_Group==3] <- 1


Med_fit_temp <- lm(Microbe~Eat_Lancet_Score_Group+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+bristol_scale,data=data_fit_temp_new)
Out_fit_temp  <- lm(Protein~Microbe+Eat_Lancet_Score_Group+age+sex+BMI+smoke+alc+factor(edu3)+factor(income4)+MET+energy+DM_diagnosis_v1+Hyper_diagnosis_v1+dyslipi_med_new+factor(phase),data=data_fit_temp_new)

library(mediation)
set.seed(1234)
Med_s49_p01023_new <- mediate(Med_fit_temp, Out_fit_temp, treat = "Eat_Lancet_Score_Group",mediator = "Microbe")
summary(Med_s49_p01023_new)


