library(xgboost)
library(pROC)
library(dplyr)
library(Matrix)
library(ggplot2)

load("saved_variables.RData")

# trainset
train_y<-ifelse(train_sample_use$sample_type=="AML",1,0)
dtrain <- xgb.DMatrix(data = as.matrix(train_sample_use[,1:(ncol(train_sample_use)-2)]), label = train_sample_use$sample_type_bin)

# valiset
vali_y<-ifelse(vali_sample_use$sample_type=="AML",1,0)
dvali  <- xgb.DMatrix(data = as.matrix(vali_sample_use[,1:(ncol(vali_sample_use)-2)]), label = vali_sample_use$sample_type_bin)


# 使用已有模型进行预测
train_preds <- predict(bst, dtrain)
vali_preds <- predict(bst, dvali)

# 基于预测的概率获取预测标签
train_labels <- case_when(
  train_preds > 0.5  ~ 1,
  train_preds < 0.5  ~ 0,
  TRUE               ~ 0.5
)

vali_labels <- case_when(
  vali_preds > 0.5  ~ 1,
  vali_preds < 0.5  ~ 0,
  TRUE              ~ 0.5
)

# 创建混淆矩阵
xgb.cf_train <- caret::confusionMatrix(as.factor(train_labels), as.factor(train_sample_use$sample_type_bin))
xgb.cf_vali <- caret::confusionMatrix(as.factor(vali_labels), as.factor(vali_sample_use$sample_type_bin))

# 打印混淆矩阵
print(xgb.cf_train)
print(xgb.cf_vali)

# 使用原始的预测概率计算ROC对象
train_roc_obj <- roc(train_sample_use$sample_type_bin, train_preds)
vali_roc_obj <- roc(vali_sample_use$sample_type_bin, vali_preds)

# 获取训练集和验证集的AUC值，并为图例生成相应的标签
train_auc_label <- sprintf("Train AUC: %.2f", auc(train_roc_obj))
vali_auc_label <- sprintf("Test AUC: %.2f", auc(vali_roc_obj))

# 创建ROC曲线
roc_plot <- ggplot() + 
  geom_line(data = data.frame(sensitivities = train_roc_obj$sensitivities, specificities = train_roc_obj$specificities), 
            aes(x = sensitivities, y = specificities, colour = train_auc_label), linewidth = 2) +
  geom_line(data = data.frame(sensitivities = vali_roc_obj$sensitivities, specificities = vali_roc_obj$specificities), 
            aes(x = sensitivities, y = specificities, colour = vali_auc_label), linewidth = 2) +
  ggtitle("ROC Curves") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 2),
        axis.ticks = element_line(linewidth = 1, linetype = 1),
        plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 2, hjust = 0.5),
        axis.text = element_text(size = 15, face = "bold", colour = "black"),
        legend.title = element_blank(),
        legend.key.width = unit(0.7, "cm"),
        legend.position = "bottom") +
  scale_color_manual(values = c("#00448C", "#EE0000")) +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
roc_plot

save(list = ls(),file = "malignant_data.Rdata")

