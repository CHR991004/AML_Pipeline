AML_signature <- read.table("./Enrichment_analysis/input/AML_signature.txt",sep = "\t",header = T)

# 检索出的基因富集分析
getOption("clusterProfiler.download.method")
# install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto')
library("clusterProfiler")
library("org.Hs.eg.db")
library("topGO")
library("DOSE")
library("ggplot2")
library("ggsci")

Sys.setlocale(category = "LC_ALL",locale = "English_United States.1252")
genes=as.vector(AML_signature[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)

ALL <- enrichGO(entrezIDs, OrgDb = org.Hs.eg.db ,ont="ALL",
                pvalueCutoff = 0.05,pAdjustMethod = 'none',
                qvalueCutoff = 2)
GO_results<-ALL@result
BP<-subset(GO_results,ONTOLOGY=="BP")[1:5,]
CC<-subset(GO_results,ONTOLOGY=="CC")[1:5,]
MF<-subset(GO_results,ONTOLOGY=="MF")[1:5,]
GO_data<-rbind(BP,CC,MF)
GO_data$Description<-factor(GO_data$Description,levels=GO_data$Description)

GO_plot=ggplot(data = GO_data, mapping = aes(x=Description,y=Count,fill=ONTOLOGY))+
  geom_bar(stat="identity")+
  coord_flip()+
  ggtitle("GO plot")+
  scale_x_discrete(limits=rev(levels(GO_data$Description)))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=12,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", linewidth = 3, linetype="solid"),
        axis.ticks = element_line(linewidth =  2),
        plot.title = element_text(size = 30, face = "bold"))+
  scale_fill_lancet()
GO_plot
write.table(GO_results,"GO_enrichment_results.txt",sep = "\t",row.names = F,col.names = T,quote = F)
ggsave(GO_plot, filename = "GO.pdf", width = 8, height = 7)

KEGG_results <- enrichKEGG(gene = entrezIDs, organism = "hsa", pvalueCutoff = 2, qvalueCutoff = 2)
KEGG_data <- KEGG_results@result[1:10,]
KEGG_data$Description<-factor(KEGG_data$Description,levels=KEGG_data[order(KEGG_data$pvalue,decreasing = T),]$Description)
KEGG_data$pvalue<--log10(KEGG_data$pvalue)
KEGG_plot <- ggplot(KEGG_data, aes(x=Description, y=Count,fill=pvalue)) +
  geom_bar(stat="identity", width=0.8) + # width可设置条形图宽度
  coord_flip() +
  xlab("") +
  ggtitle("KEGG plot",)+
  theme_bw()+
  scale_fill_gradient(name="-log10(P)",low = "#e2e071",  #设置p值low 和 high的颜色,最终是两种颜色渐变
                      high = "#7f3185")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.title = element_text(size = 15,face = "bold",colour = "black"), # 可以设置标题字体大小、是否加粗、颜色
        axis.text = element_text(size=15,face = "bold",colour = "black"), #同理
        panel.border = element_rect(fill=NA,color="black", linewidth=3, linetype="solid"), # 外围加框，fill为填充色
        axis.ticks = element_line(linewidth = 2),
        plot.title = element_text(size = 30, face = "bold")) # 坐标轴线条加粗
KEGG_plot
ggsave("KEGG_plot.pdf",KEGG_plot,width = 8, height = 7)
write.table(KEGG_results,"KEGG_enrichment_results.txt",sep = "\t",row.names = F,col.names = T,quote = F)


library(simplifyEnrichment)
go_id <- GO_results[,"ID"]
mat <-  GO_similarity(go_id, 
                      ont =  c("BP", "CC", "MF"),
                      db = 'org.Hs.eg.db',
                      measure = "Rel",
                      remove_orphan_terms = F)
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
cl = binary_cut(mat)
export_to_shiny_app(mat, cl)

# df = simplifyGO(mat)
save(list = ls(),file = "enrichment_reuslts.Rdata")
load("enrichment_reuslts.Rdata")
library(org.Hs.eg.db)
library(clusterProfiler)
setReadable(KEGG_results,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
KEGG_selected <- read.table("Enrichment_analysis/output/KEGG_select.txt",sep = "\t",header = F)[,1]
KEGG_data_selected <-KEGG_results[KEGG_results$Description%in%KEGG_selected,]
KEGG_data_selected $pvalue<--log10(KEGG_data_selected $pvalue)
KEGG_data_selected <- KEGG_data_selected %>%
  arrange(pvalue)
# 将Description转化为一个因子，其水平按照pvalue的顺序
KEGG_data_selected$Description <- factor(KEGG_data_selected$Description, levels = KEGG_data_selected$Description)

KEGG_plot <- ggplot(KEGG_data_selected, aes(x=Description, y=Count,fill=pvalue)) +
  geom_bar(stat="identity", width=0.8) + # width可设置条形图宽度
  coord_flip() +
  xlab("") +
  # ggtitle("KEGG plot",)+
  theme_bw()+
  scale_fill_gradient(name="-log10(P)",low = "#e2e071",  #设置p值low 和 high的颜色,最终是两种颜色渐变
                      high = "#7f3185")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.title = element_text(size = 12,colour = "black"), # 可以设置标题字体大小、是否加粗、颜色
        axis.text = element_text(size=12,colour = "black"), #同理
        panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"), # 外围加框，fill为填充色
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(size=12,  colour="black"),
        legend.title = element_text(size=12,  colour="black")) # 坐标轴线条加粗
KEGG_plot
ggsave("./Enrichment_analysis/output/KEGG_plot.pdf",KEGG_plot,width = 7, height = 6)
