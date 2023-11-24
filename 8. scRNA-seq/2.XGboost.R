library(xgboost)
library(pROC)
library(randomForest)
library(dplyr)
library(Matrix)
library(ggplot2)
load("xgboost_model_for_AML_scRNA-seq_191genes.Rdata")

# 1. 数据准备
# 将因变量转为二进制

# load("train_vali_LS_matrix.Rdata")
genes<-read.table("AML_signature_diff.txt",sep = "\t",header = T)[,1]
genes <- intersect(genes[genes%in%rownames(train_matrix)],genes[genes%in%rownames(vali_matrix)])
train_sample_use <- as.data.frame(t(train_matrix[genes,]))
train_sample_use$sample_type <- as.character(train_scRNA_obj@meta.data$sample_type)

vali_sample_use <- as.data.frame(t(vali_matrix[genes,]))
vali_sample_use$sample_type <- as.character(vali_scRNA_obj@meta.data$sample_type)

train_sample_use$sample_type_bin <- as.numeric(train_sample_use$sample_type == "AML")
vali_sample_use$sample_type_bin <- as.numeric(vali_sample_use$sample_type == "AML")

# vali_sample_use
# 为XGBoost准备数据
dtrain <- xgb.DMatrix(data = as.matrix(train_sample_use[,1:(ncol(train_sample_use)-2)]), label = train_sample_use$sample_type_bin)
dvali  <- xgb.DMatrix(data = as.matrix(vali_sample_use[,1:(ncol(vali_sample_use)-2)]), label = vali_sample_use$sample_type_bin)

# 2. 训练XGBoost模型
params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  max_depth = 3,
  eta = 0.01,
  colsample_bytree = 0.5
)
bst <- xgb.train(params, dtrain, nrounds = 1000, watchlist = list(train=dtrain, vali=dvali))
importance_matrix <- xgb.importance(feature_names = bst$feature_names, model = bst)

# 数据预处理以适应堆积柱状图
melted_data <- reshape2::melt(importance_matrix, id.vars = "Feature")

# 使用ggplot2来创建堆积柱状图
ggplot(melted_data, aes(x=Feature, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#03478A", "#E80809", "#44B543")) +
  ggtitle("Feature Importance on Multiple Metrics") +
  xlab("Features") +
  ylab("Metric Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = importance_matrix$Feature)  # 对x轴进行排序，如果需要

# 3. 使用模型进行预测
train_preds <- predict(bst, dtrain)
train_labels <- ifelse(train_preds >= 0.5, 1, 0)

vali_preds <- predict(bst, dvali)
vali_labels <- ifelse(vali_preds >= 0.5, 1, 0)

train_roc_obj <- roc(train_sample_use$sample_type_bin, train_preds)
vali_roc_obj <- roc(vali_sample_use$sample_type_bin, vali_preds)

# 获取训练集和验证集的AUC值，并为图例生成相应的标签
train_auc_label <- sprintf("Train AUC: %.2f", auc(train_roc_obj))
vali_auc_label <- sprintf("Test AUC: %.2f", auc(vali_roc_obj))
roc_plot <- ggplot() + 
  geom_line(data = data.frame(sensitivities = train_roc_obj$sensitivities, specificities = train_roc_obj$specificities), 
            aes(x = sensitivities, y = specificities, colour = train_auc_label), size = 2) +
  geom_line(data = data.frame(sensitivities = vali_roc_obj$sensitivities, specificities = vali_roc_obj$specificities), 
            aes(x = sensitivities, y = specificities, colour = vali_auc_label), size = 2) +
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
ggsave(filename = "roc_plot.pdf",plot = roc_plot,width = 6,height = 6.5)
save(list = ls(),file = "xgboost_model_for_AML_scRNA-seq_191genes.Rdata")


load("xgboost_model_for_AML_scRNA-seq_191genes.Rdata")
load("umap.Rdata")

train_umap_data<-cbind(train_umap,train_preds,pred = train_labels,orig = train_sample_use$sample_type_bin,pred_score = train_preds)%>%as.data.frame()
vali_umap_data<-cbind(vali_umap,vali_preds,pred = vali_labels,orig =vali_sample_use$sample_type_bin,pred_score = vali_preds)%>%as.data.frame()

library(ggplot2)
library(viridis)
library(dplyr)

# 为train_umap_data预处理
train_umap_data <- train_umap_data %>%
  mutate(
    shape = ifelse(orig == 0, "triangle", "circle"),
    color = ifelse(pred == orig, 
                   ifelse(pred == 0, "#034788", "#E60808"), 
                   "#E0E0E0")
  )
train_umap_data$color <- factor(train_umap_data$color,levels = c("#034788","#E60808","#E0E0E0"))
train_umap_data <- train_umap_data[order(train_umap_data$pred_score,decreasing = T),]

# 为vali_umap_data预处理
vali_umap_data <- vali_umap_data %>%
  mutate(
    shape = ifelse(orig == 0, "triangle", "circle"),
    color = ifelse(pred == orig, 
                   ifelse(pred == 0, "#034788", "#E60808"), 
                   "#E0E0E0")
  )
vali_umap_data$color <- factor(vali_umap_data$color,levels = c("#034788","#E60808","#E0E0E0"))
vali_umap_data <- vali_umap_data[order(vali_umap_data$pred_score,decreasing = T),]


library(ggplot2)
library(gridExtra)
plot1 <-ggplot(train_umap_data, aes(x=UMAP_1, y=UMAP_2, shape=shape, color=color)) +
  geom_point(size=3, alpha=0.8) +
  scale_shape_manual(name="Original Type", values=c("triangle"=17, "circle"=16), labels=c("AML", "Normal")) +
  scale_color_manual(name="Predicted Type",
                     values=c("#034788","#E60808","#E0E0E0"),
                     labels=c("Predicted Normal", "Predicted AML", "Mismatched")) +
  theme_minimal() +
  labs(title="Train Data: Original vs Predicted")
# plot1
plot2 <- ggplot(vali_umap_data, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(shape=shape, color=color), size=3, alpha=0.8) +  
  scale_shape_manual(name="Original Type", values=c("triangle"=17, "circle"=16), labels=c("AML", "Normal")) +
  scale_color_manual(name="Predicted Type",
                     values=c("#034788","#E60808","#E0E0E0"),
                     labels=c("Predicted Normal", "Predicted AML", "Mismatched")) +
  theme_minimal() +
  labs(title="Test Data: Original vs Predicted")

plot3 <- ggplot(train_umap_data, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(shape=shape, color=pred_score), size=3, alpha=0.5) +  
  scale_shape_manual(name="Original Type", values=c("triangle"=17, "circle"=16), labels=c("AML", "Normal")) +
  scale_color_viridis(name="Pred Score") +
  theme_minimal() +
  labs(title="Train Data: Original Type with Pred Score")

plot4 <- ggplot(vali_umap_data, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(shape=shape, color=pred_score), size=3, alpha=0.5) +  
  scale_shape_manual(name="Original Type", values=c("triangle"=17, "circle"=16), labels=c("AML", "Normal")) +
  scale_color_viridis(name="Pred Score") +
  theme_minimal() +
  labs(title="Test Data: Original Type with Pred Score")

final_plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
ggsave(filename = "merge.pdf",plot = final_plot,height = 6,width = 10)
save(list = ls(),file = "tem_save_for_plot.Rdata")
