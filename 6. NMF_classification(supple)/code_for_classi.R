library(NMF)
library(pheatmap)
library(survival)
library(survminer)
library(limma)
setwd("../7. AML_machine_learning/")
# load("classification_data.Rdata")
# load("beatAML_all_data.Rdata")
load("classification_data_survival_2class.Rdata")
# 数据预处理
# 计算每个基因的方差
gene_variances <- apply(beatAML_FPKM, 1, var)

# 保留方差不为0的基因
beatAML_FPKM_filtered <- beatAML_FPKM[gene_variances > 0.5, ]
beatAML_FPKM_filtered[1:5, 1:5]

expr <- as.matrix(beatAML_FPKM_filtered)
expr <- expr[rownames(expr) %in% AML_signature,]

# 使用NMF进行一致性聚类
set.seed(1234)
res <- nmf(expr, rank = 2:10, nrun = 50, method = 'brunet')

library(ConsensusClusterPlus)

library(NbClust)
pdf(file = "CDF.pdf",width = 6,height = 6)
results <- ConsensusClusterPlus(
  d = expr, 
  maxK = 10, 
  reps = 50, 
  pItem = 0.8, 
  pFeature = 1, 
  clusterAlg = "km",
  title = "AML_cluster",
  distance = "pearson"
)
dev.off()
results <- ConsensusClusterPlus(
  d = expr, 
  maxK = 10, 
  reps = 50, 
  pItem = 0.8, 
  pFeature = 1, 
  clusterAlg = "km",
  title = "AML_cluster",
  distance = "pearson"
)

final_rank <- 2
res_final <- nmf(expr, rank = final_rank, nrun = 100, method = 'brunet')

# 获取样本的聚类结果
H <- coef(res_final)
clusters <- apply(H, 2, which.max)
table(clusters)
beatAML_cli$cluster <- clusters[as.character(beatAML_cli$patient)]

# 获取聚类结果(图)
consensusmap(res_final)
pdf("beatAML_consensusmap_output.pdf")
consensusmap(res_final)
dev.off()

# # 为每个组构建生存函数

subset_data <- beatAML_cli
km_fit_subset <- survfit(Surv(futime, fustat) ~ cluster, data = subset_data)
pdf(file = "beatAML_survival_plot_NMF.pdf",width = 7,height = 7)
ggsurvplot(
  km_fit_subset,
  data = subset_data,
  risk.table = T,
  pval = TRUE,
  conf.int = TRUE,
  title = "Kaplan-Meier Survival Curve",
  xlab = "Time(days)",
  ylab = "Survival Probability",
  palette = c("#E60808", "#034788")  # 根据cluster值分配颜色
)
dev.off()

design <- model.matrix(~ 0 + factor(clusters))
colnames(design) <- c("Cluster1", "Cluster2")
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(Cluster2-Cluster1, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
sig_genes <- topTable(fit2, adjust = "fdr", n = Inf)
threshold <- 1
filtered_genes <- sig_genes[abs(sig_genes$logFC) > threshold,]
sig_genes <- filtered_genes[filtered_genes$adj.P.Val < 0.05, ]
selected_expr <- expr[rownames(sig_genes),]
z_scores <- t(apply(selected_expr, 1, function(x) (x - mean(x)) / sd(x)))
sorted_samples <- names(clusters)[order(clusters)]
z_scores[z_scores >= 2] <- 2
z_scores[z_scores <= -2] <- -2
pdf(file = "beatAML_heatmap.pdf",width = 7,height = 5)
pheatmap(
  z_scores[, sorted_samples],
  cluster_rows = T,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_col = data.frame(Cluster = factor(clusters[sorted_samples], levels = 1:2))
)
dev.off()

gene_cluster = hclust(dist(z_scores))
# 假设你希望将基因划分为k个类
k = 2 # 你可以根据需要更改这个值
gene_clusters = cutree(gene_cluster, k)
clustered_genes_list = split(rownames(z_scores), gene_clusters)

library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

# 创建一个空的数据框，以便于后续的数据填充
result_df <- data.frame(gene = character(), group = numeric())

# 遍历列表
for (group_name in names(clustered_genes_list)) {
  temp_df <- data.frame(
    gene = clustered_genes_list[[group_name]],
    group = group_name
  )
  result_df <- rbind(result_df, temp_df)
}
write.table(result_df,file = "beatAML_group_result.xls",sep = "\t",quote = F,row.names = F,col.names = T)
