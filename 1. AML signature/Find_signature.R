
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape)
library(corrplot)
library(ggpubr)
library(do)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

rm(list = ls())
setwd("./AML_signature/")

blood_output <- read.table("multidata_genematrix.txt",sep = "\t",header = T)

data <- blood_output[,2:(ncol(blood_output)-1)]
plot_matrix <- matrix(nrow = ncol(data),ncol = 3,dimnames = list(NULL,c("name","Upregulated","Downregulated")))
for (x in 1:ncol(data)) {
  data_name <- colnames(data)[x]
  non_na_data <- na.omit(data[,x])
  Upregulated <- length(non_na_data[non_na_data>0])
  Downregulated <- length(non_na_data[non_na_data<0])
  plot_matrix[x,] <- c(data_name,Upregulated,Downregulated)
}
plot_matrix <- as.data.frame(plot_matrix)
plot_matrix <- melt(plot_matrix,id.vars = c("name"))

plot_matrix$value <- as.numeric(plot_matrix$value) 

pl <- ggplot(data=plot_matrix, aes(x=reorder(name,-value), y=value)) +
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values=c("#005187","#e5082c"))+
  scale_y_continuous(breaks=c(1000,2000,3000),
                     labels=c("1000", "2000", "3000"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  labs(x="",y="# of reported DE genes",title = "Reanalysis")+
  theme(text = element_text(family = "Arial",face = "bold"))
pl

ggsave(pl, filename = "AML_figure_1c.pdf", device = cairo_pdf, 
       width = 8, height = 7, units = "in")


multidata <- read.table("multidata_genematrix.txt",sep = "\t",header = T,check.names = F)[,c(1:27)]
length_compute <- function(x){
  length(x[!is.na(x)])-1
}

genematrix_count <- apply(multidata,1,length_compute)
multidata <- cbind(multidata,count=genematrix_count)


point_plot_data <- as.data.frame(table(multidata$count))
point_plot_data[,1] <- as.numeric(point_plot_data[,1])
supple_Fig1D <- ggplot(data = point_plot_data, aes(x=Var1,y=Freq)) +
  geom_smooth(mapping = aes(x=Var1,y=Freq),color="#4baf49",
              method = "loess",se = F,span=1,formula = y~I(1/x),size=1.5) +
  geom_point(shape=1,size=3,color="#4baf49",stroke=2)+
  labs(x="number of overlaps across studies",y="DE genes")+
  ylim(0,6500)+
  scale_x_continuous(breaks = as.numeric(point_plot_data$Var1))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=15,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks = element_line(size = 2))
supple_Fig1D
ggsave(supple_Fig1D, filename = "AML_supple Fig1D.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")

rm(list = ls())
Fig1F_data <- read.table("data_nocount.txt",sep = "\t",header = T,check.names = F,row.names = 1)

heatmap_data<-matrix(data = NA,nrow = ncol(Fig1F_data),
                     ncol = ncol(Fig1F_data),
                     dimnames = list(colnames(Fig1F_data),colnames(Fig1F_data)))

two_inter <- list()
for (i in colnames(Fig1F_data)) {
  for (x in colnames(Fig1F_data)) {
    two_inter[[i]][[x]] <- rownames(na.omit(Fig1F_data[,c(i,x)]))
    heatmap_data[i,x] <- nrow(na.omit(Fig1F_data[,c(i,x)]))
  }
}
heatmap_data[upper.tri(heatmap_data)]<-NA

fraction <- c()
for (n in names(two_inter)) {
  inter_names <- names(two_inter)[names(two_inter)!=n]
  rname <- c()
  for (n2 in inter_names) {
    rname <- unique(c(rname,two_inter[[n]][[n2]]))
  }
  fraction <- c(fraction,length(rname)/length(Fig1F_data[,n][!is.na(Fig1F_data[,n])]))
}
heatmap_data_gg<-reshape2::melt(heatmap_data)
heatmap_plot<-	ggplot(heatmap_data_gg,aes(x=Var1,y=Var2,fill=log2(value)))+
  geom_tile(color = "white",na.rm = T,size=0.5)+
  scale_fill_gradient2(low = "#800080", high = "#faed5d", mid = "#65aeb0", 
                       midpoint = 6,  space = "Lab",na.value = "white",guide = "colourbar") +
  labs(x="",y="",title = "Re-analysis")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        panel.border = element_rect(fill = NA,color = "black",size = 2),
        axis.ticks = element_line(size = 1,linetype = 1),
        legend.text = element_blank(),
        plot.title = element_text(size=18,face = "bold",colour = "black",vjust = 2,hjust = 0.5),
        axis.text = element_text(size=15,face = "bold",colour = "black"),
        legend.title = element_blank(),legend.key.width = unit(0.7,"cm") 
  )+
  coord_fixed(1)+
  geom_text(aes(label=value),colour= "white",size=2.5)
heatmap_plot

fraction_data <- as.data.frame(cbind(name=rep("col",length(fraction)),fraction,name_true = colnames(Fig1F_data)))
fraction_data$fraction <- as.numeric(fraction_data$fraction)
fraction_data$name_true<-factor(fraction_data$name_true,levels = fraction_data$name_true)

inter_ratio <- ggplot(data = fraction_data,aes(x=name,y=name_true,fill=fraction))+
  geom_tile(color = "white",size = 0.5)+
  theme_bw()+
  scale_fill_gradient2(low = "#800080", high = "#faed5d", mid = "#65aeb0", 
                       midpoint = 0.5,  space = "Lab",na.value = "white")+
  coord_fixed(ratio=1)+
  geom_text(aes(label=signif(as.numeric(as.character(fraction)), 2)),colour="white")+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none",
        panel.border = element_rect(fill = NA,color = "black",size = 2),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_line(size = 2),
        text = element_text(color = "white")
  )

inter_ratio
all_plot <- ggarrange(heatmap_plot, inter_ratio,
                      ncol = 2, nrow = 1, 
                      widths = c(12, 1), heights = 9,
                      common.legend = TRUE,hjust = 0)	
ggsave(all_plot, filename = "supple Fig1F.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")

rm(list = ls())

sample1 <- sample(1:26, size = 12)
sample2 <- sample(1:26, size = 12)
sample3 <- sample(1:26, size = 12)
sample4 <- sample(1:26, size = 12)
sample_list <- rbind(sample1,sample2,sample3,sample4)
write.table(sample_list,"sample_list.txt",sep = "\t",quote = F,row.names = T,col.names = F)
Fig1E_data <- read.table("data_nocount.txt",sep = "\t",header = T,check.names = F,row.names = 1)[,sample4]

length_compute <- function(x){
  length(x[!is.na(x)])
}
GMC_count_row <- function(one_col,rt){
  count_all <- rt[,"count"]
  one_col <- count_all[!is.na(one_col)]
  count_mean <- mean(one_col)
}
GMC_count_all <- function(data_set,Fig1E_data,	GMC_count_row){
  rt <- Fig1E_data[,data_set]
  rt	<- cbind(rt, count=Fig1E_data$count)
  if ((ncol(rt)>=3)) {
    row_gene_mean<-mean(apply(rt[,1:(ncol(rt)-1)],2,FUN = GMC_count_row,rt=rt))
  }else{
    row_gene_mean<-apply(as.matrix(rt[,1]),2,FUN = GMC_count_row,rt=rt)
  }
  
}

CV_all <- function(data_set,Fig1E_data){
  rt <- as.matrix(Fig1E_data[,data_set])
  delete_row <- which(apply(rt,1,function(x) all(is.na(x))))
  if (length(delete_row)>0) {
    rt<-rt[-delete_row,]
  }
  gene_EXP <- as.vector(t(rt))
  gene_EXP <- gene_EXP[!is.na(gene_EXP)]
  rt_CV <- sd(gene_EXP)/mean(gene_EXP)
}

genematrix_count <- apply(Fig1E_data,1,length_compute)
Fig1E_data <- cbind(Fig1E_data,count=genematrix_count)
Fig1E_data <- Fig1E_data[Fig1E_data$count>1,]
Fig1E_CV_data <- Fig1E_data
Fig1E_CV_data[is.na(Fig1E_CV_data)]<-0
result_all <- list()

for (i in 1:11) {
  data_set <- t(combn(colnames(Fig1E_data)[-13],m = i))
  CV_result <- apply(data_set,1,FUN = CV_all,Fig1E_data=Fig1E_data)
  GMC_count_result <- apply(data_set,MARGIN = 1,FUN = GMC_count_all,Fig1E_data=Fig1E_data,GMC_count_row=GMC_count_row)
  plot_data_single <- cbind(CV=CV_result,GMC=GMC_count_result,type=i)
  result_all[[i]] <- plot_data_single
}
plot_data <- as.data.frame(do.call(what = rbind,args = result_all))
plot_data[plot_data<0]<-abs(plot_data[plot_data<0])
plot_data$GMC2<-(plot_data$GMC-min(plot_data$GMC))/(max(plot_data$GMC)-min(plot_data$GMC))*12

plot_data_article <- matrix(nrow = 11,ncol = 3,dimnames = list(NULL,c("GMC","CV","type")))
plot_data_article <- as.data.frame(plot_data_article)

supple_Fig1D<-ggplot(data = plot_data, aes(x=GMC,y=CV)) +
  ylim(0,10)+
  geom_point(aes(x=GMC2,y=CV,colour=factor(type)),shape=1,size=3,stroke=2)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  # geom_smooth(data = plot_data,mapping = aes(x=GMC,y=CV,colour=factor(type)),
  # 												method = "loess",se = F,span=1,formula = y~I(x),size=1.5)+
  scale_color_manual(values = c("1"="#2c3c97",
                                "2"="#0f75b6",
                                "3"="#4faed3",
                                "4"="#95d8e7",
                                "5"="#d6f1f7",
                                "6"="#faf5be",
                                "7"="#ffdf89",
                                "8"="#ffad5b",
                                "9"="#ff6336",
                                "10"="#fa2b1d",
                                "11"="#bf1625"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        panel.border = element_rect(fill = NA,color = "black",size = 2),
        axis.ticks = element_line(size = 1,linetype = 1),
        plot.title = element_text(size=18,face = "bold",colour = "black",vjust = 2,hjust = 0.5),
        axis.text = element_text(size=15,face = "bold",colour = "black"),
        legend.title = element_blank(),legend.key.width = unit(0.7,"cm") 
  )

ggsave(supple_Fig1D, filename = "supple Fig1D_list4.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")

x = 1:30
y = 1-1/x
Fig1F_data <- as.data.frame(cbind(x,y))
Fig1F <- ggplot(data = Fig1F_data, aes(x=x,y=y)) +
  geom_smooth(mapping = aes(x=x,y=y),color="#e5082c",
              method = "loess",se = F,span=1,formula = y~I(1/x),size=2) +
  geom_point(shape=1,size=4,color="#000000",stroke = 1.5)+
  labs(x="number of studys",y="resolution")+
  ylim(0,1)+
  scale_x_continuous(breaks=c(5,10,15,20,25,30),
                     labels=c("5", "10", "15","20","25","30"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=15,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=3,
                                    linetype="solid"),
        axis.ticks = element_line(size = 2))+
  geom_segment(data=Fig1F_data,aes(x=26,y=-Inf,xend=x[26],yend=y[26]),
               col="#ff8d00",size=1,lty=8)+
  geom_segment(data=Fig1F_data,aes(x=-Inf,y=y[26],xend=26,yend=y[26]),
               col="#ff8d00",size=1,lty=8)

Fig1F
ggsave(Fig1F, filename = "Fig1F.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")

Fig1E_data <- read.table("multidata_genematrix.txt",
                         sep = "\t",header = T,check.names = F,row.names = 1)
AML_signature <- rownames(Fig1E_data)[Fig1E_data$count>=6]
Fig1E_plot_data <- matrix(nrow = 26,ncol = 2,
                          dimnames = list(colnames(Fig1E_data)[1:26],
                                          c("similarity","gene in SR")))
for (research_id in 1:26) {
  single_research <- rownames(Fig1E_data)[!is.na(Fig1E_data[,research_id])]
  Fig1E_plot_data[research_id,"gene in SR"] <- length(single_research)
  inter_gene <- intersect(AML_signature,single_research)
  Fig1E_plot_data[research_id,"similarity"] <- 1-(length(inter_gene)/length(single_research))
}

Fig1E_plot_data[,1]<-scales::rescale(Fig1E_plot_data[,1],to = c(0.1,1))
Fig1E_plot_data_scaled <- cbind(name=rownames(Fig1E_plot_data),Fig1E_plot_data)%>%as.data.frame()
Fig1E_plot_data_scaled$name<- factor(Fig1E_plot_data_scaled$name,levels = Fig1E_plot_data_scaled$name)
Fig1E_plot_data_scaled$`gene in SR`<-as.numeric(Fig1E_plot_data_scaled$`gene in SR`)

Fig1E_plot<-ggplot(Fig1E_plot_data_scaled, aes(x=name, y=similarity))+
  geom_point(aes(x=name,y=similarity, size=`gene in SR`),shape=1, stroke = 1)+
  coord_polar()+
  geom_segment(aes(y=0, xend=name, yend=similarity))+
  theme_bw()+
  scale_size(range = c(3,10))

ggsave(Fig1E_plot, filename = "Fig1E_plot.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")

rm(list = ls())
Fig2A_data <- read.table("multidata_genematrix.txt",
                         sep = "\t",header = T,check.names = F,row.names = 1)
naomit_mean<-function(x){
  x<-x[!is.na(x)]%>%mean()
}

Fig2A_data$logFC<-apply(Fig2A_data[,1:26],1,naomit_mean)
AML_signature_matrix <- Fig2A_data[Fig2A_data$count>=6,]
diffAll <- cbind(genes=rownames(AML_signature_matrix),logFC=AML_signature_matrix$logFC)%>%as.data.frame()

entrez_list <- as.data.frame(t(as.data.frame(mget(diffAll$genes, org.Hs.egSYMBOL2EG, ifnotfound=NA),row.names = c("entrez1","entrez2"))))
entrez_list <- cbind(genes=rownames(entrez_list),entrez_list)
entrez_list$genes<-str_replace(entrez_list$genes,"\\.","-")
entrez_logFC <- merge(entrez_list,diffAll,by="genes")%>%
  unite(., "entrez", entrez1,entrez2,sep = ",")%>%
  separate_rows(.,entrez,sep = ",")%>%.[!duplicated(.$entrez),]

AML_signature <- entrez_logFC$genes
AML_signature_up <- entrez_logFC$genes[entrez_logFC$logFC>0]
AML_signature_down <- entrez_logFC$genes[entrez_logFC$logFC<0]

AML_entrez_up <- entrez_logFC[entrez_logFC$logFC>0,"entrez"]$entrez
AML_entrez_down <- entrez_logFC[entrez_logFC$logFC<0,"entrez"]$entrez

GO_up <- enrichGO(AML_entrez_up, OrgDb = org.Hs.eg.db ,ont="ALL",pvalueCutoff = 0.05,pAdjustMethod = 'fdr',qvalueCutoff = 0.1)
GO_down <- enrichGO(AML_entrez_down, OrgDb = org.Hs.eg.db ,ont="ALL",pvalueCutoff = 0.05,pAdjustMethod = 'fdr',qvalueCutoff = 0.1)

GO_up_All <- GO_up@result
GO_down_All <- GO_down@result

CC_inner <- intersect(GO_up_All[GO_up_All$ONTOLOGY=="CC","Description"],GO_down_All[GO_down_All$ONTOLOGY=="CC","Description"])
GO_up <- GO_up_All[GO_up_All$Description%in%CC_inner,]
GO_down <- GO_down_All[GO_down_All$Description%in%CC_inner,]
GO_up<-cbind(GO_up,type=rep("up",nrow(GO_up)))
GO_down<-cbind(GO_down,type=rep("down",nrow(GO_up)))
GO_down$Count<--(GO_down$Count)
GO_all<-rbind(GO_up,GO_down)
GO_plot <- ggplot(GO_all,aes(x = reorder(Description,-Count),y = Count))+
  geom_bar(aes(fill = type),stat = "identity", width=0.5) + 
  coord_flip()+
  scale_fill_manual(values=c("#005187","#e5082c"))+
  geom_text(aes(x=Description,y=Count,label=as.character(abs(Count))),size=5)+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=12,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks = element_line(size = 2))+
  labs(x="",y="Number of genes",title = "")
GO_plot
ggsave(GO_plot, filename = "Fig2A_plot.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")

KEGG_up <- enrichKEGG(AML_entrez_up, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.1)
KEGG_down <- enrichKEGG(AML_entrez_down, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.1)
KEGG_up_All <- KEGG_up@result
KEGG_down_All <- KEGG_down@result
KEGG_inner <- intersect(KEGG_up_All$Description,KEGG_down_All$Description)
KEGG_up <- KEGG_up_All[KEGG_up_All$Description%in%KEGG_inner,]
KEGG_down <- KEGG_down_All[KEGG_down_All$Description%in%KEGG_inner,]
KEGG_all <- rbind()


ggplot(GO_all,aes(x = reorder(Description,-Count),y = Count))+
  geom_bar(aes(fill = type),stat = "identity", width=0.5) + 
  coord_flip()+
  scale_fill_manual(values=c("#005187","#e5082c"))+
  geom_text(aes(x=Description,y=Count,label=as.character(abs(Count))),size=5)+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=12,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks = element_line(size = 2))+
  labs(x="",y="Number of genes",title = "")

write.table(entrez_logFC,"AML_signature.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(KEGG_up,"KEGG_up.txt",sep = "\t",quote = F,row.names = F,col.names = T)

