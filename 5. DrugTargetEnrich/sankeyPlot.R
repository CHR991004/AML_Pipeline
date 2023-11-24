# drugbank_data_preprocess

# Drug_target_DGIdb <- read.table("Drug_target_DGIdb.xls",sep = "\t",header = T,quote = "")
# Drug_target_DGIdb$drug_name <- str_to_lower(Drug_target_DGIdb$drug_name)
# Drug_target_DGIdb <- Drug_target_DGIdb[Drug_target_DGIdb$genes!="",]%>%
#   unique()
# Drug_target_Drugbank <- read.table("Drug_target_Drugbank.xls",sep = "\t",header = T,quote = "")
# Drug_target_Drugbank$drug_name <- str_to_lower(Drug_target_Drugbank$drug_name)
# Drug_target_Drugbank <- unique(Drug_target_Drugbank)
# 
# write.table(Drug_target_DGIdb,"Drug_target_DGIdb_usage.xls",sep = "\t",row.names = F,col.names = T,quote = F)
# write.table(Drug_target_Drugbank,"Drug_target_Drugbank_usage.xls",sep = "\t",row.names = F,col.names = T,quote = F)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

load("preprocessed_data2.Rdata")
Drug_target_DGIdb <- read.table("Drug_target_DGIdb_usage.xls",sep = "\t",header = T,quote = "")
Drug_target_Drugbank <- read.table("Drug_target_Drugbank_usage.xls",sep = "\t",header = T,quote = "")

input_genes <- read.table("selected_genes.txt",sep = "\t",header = F)

# 初始化一个空的字符向量来存储分割后的结果
split_genes <- character()

# 遍历 input_genes 的每一行
for(i in seq_len(nrow(input_genes))) {
  # 使用 strsplit 函数分割字符串，然后转换为字符向量
  temp <- unlist(strsplit(as.character(input_genes$V1[i]), ", "))
  # 将结果合并到 split_genes 字符向量中
  split_genes <- c(split_genes, temp)
}

split_genes <- unique(split_genes)

input_drugs <- read.table("selected_drugs.txt",sep = "\t",header = F)

# 初始化一个空的字符向量来存储分割后的结果
split_drugs <- character()

# 遍历 input_genes 的每一行
for(i in seq_len(nrow(input_drugs))) {
  # 使用 strsplit 函数分割字符串，然后转换为字符向量
  temp <- unlist(strsplit(as.character(input_drugs$V1[i]), ", "))
  # 将结果合并到 split_genes 字符向量中
  split_drugs <- c(split_drugs, temp)
}

split_drugs <- unique(split_drugs)

# input_genes <- c("FLT3","RUNX1","EGFR")


# 从Drug_target_Drugbank数据库中筛选出输入基因对应的数据
filtered_Drugbank <- Drug_target_Drugbank %>%
  filter(genes %in% split_genes) %>%
  mutate(source = "Drugbank")

# 从Drug_target_DGIdb数据库中筛选出输入基因对应的数据
filtered_DGIdb <- Drug_target_DGIdb %>%
  filter(genes %in% split_genes) %>%
  mutate(source = "DGIdb")

# 将两个筛选后的数据集合并成一个汇总表
combined_data <- bind_rows(filtered_Drugbank, filtered_DGIdb)

# 定义获取ATC编码的函数
get_atc_codes <- function(drug_name, atc_list) {
  codes <- character()
  
  # 遍历 atc_list
  for (atc_code in names(atc_list)) {
    # 检查药物名称是否在当前ATC编码的列表中
    if (drug_name %in% atc_list[[atc_code]]) {
      codes <- c(codes, atc_code)
    }
  }
  
  return(if (length(codes) > 0) paste(codes, collapse = ",") else NA)
}

# 应用函数到 combined_data 的 drug_name 列以获取ATC编码
combined_data$ATC_codes <- sapply(combined_data$drug_name, get_atc_codes, atc_list = atc_list_unique)
combined_data <- combined_data %>%
  # 使用mutate将ATC_codes的NA替换为"unknown"
  mutate(ATC_codes = ifelse(is.na(ATC_codes), "unknown", ATC_codes)) %>%
  # 使用separate_rows将逗号分隔的ATC_codes分开
  separate_rows(ATC_codes, sep = ",")

# 使用ATC筛选
# selected_atc_code <- read.table("selected_atc.txt",sep = "\t",header = F)[,1]

# combined_data <- combined_data[combined_data$ATC_codes%in%selected_atc_code,]

# 使用drugname筛选
combined_data <- combined_data[combined_data$drug_name%in%split_drugs,]

# 加载必要的库
library(ggalluvial)
library(ggplot2)
combined_data$genes <- as.factor(combined_data$genes)
combined_data$drug_name <- as.factor(combined_data$drug_name)
combined_data$source <- as.factor(combined_data$source)
combined_data$ATC_codes <- as.factor(combined_data$ATC_codes)
# 统计每个组合的频率
data_for_plot <- combined_data %>%
  count(genes, drug_name, ATC_codes, name = "Freq")

# 绘制桑基图
sankey_plot <- ggplot(data_for_plot,
       aes(y = Freq,
           axis1 = genes, axis2 = drug_name, axis3 = ATC_codes)) +
  geom_alluvium(aes(fill = genes),
                decreasing = FALSE) +
  geom_stratum(aes(fill = genes), decreasing = FALSE) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(length(unique(data_for_plot$genes)))) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3, angle = 0,
            decreasing = FALSE) +
  # 自定义X和Y轴上的标签
  scale_x_continuous(breaks = 1:3,
                     labels = c("genes", "drugs", "ATC_code"),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  # 设置图形主题样式
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Group", y = "Frequency")
sankey_plot
# write.table(combined_data,file = "data_for_sankey.xls",quote = F,row.names = F,col.names = T,sep = "\t")
# ggsave(filename = "sankey_output.pdf",plot = sankey_plot,width = 8,height = 10)
