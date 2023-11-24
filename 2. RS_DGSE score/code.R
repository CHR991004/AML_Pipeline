Fig4A_data <- read.table("multidata_genematrix.txt",
                         sep = "\t",header = T,check.names = F)
AML_signature <- read.table("AML_signature_diff.txt",
                            sep = "\t",header = T,check.names = F)

rank_score<-function(x){
  ifelse(x==1,scn<--1,scn<-x)
  return(scn)
}

rank_score_list <- list()
for (test in 1:26) {
  test_data <- Fig4A_data[,c(1,test+1,28)]
  score_data <- inner_join(test_data,AML_signature,by="genes")%>%na.omit()
  rank_score_list[[colnames(Fig4A_data)[test+1]]] <- sum(apply(as.data.frame(score_data[,3]),1,rank_score))
}
rank_score_list[["Aging_signiture"]]<-sum(Fig4A_data[Fig4A_data$genes%in%AML_signature$genes,"count"])

Fig4A_data <- read.table("data_nocount.txt",
                         sep = "\t",header = T,check.names = F)
AML_signature <- read.table("AML_signature_diff.txt",
                            sep = "\t",header = T,check.names = F)
AML_signature$rank<-rank(AML_signature$logFC,ties.method = "first",)
EnSc_fixed_list <- list()
for (test in 1:26) {
  test_data <- Fig4A_data[,c(1,test+1)]%>%na.omit()
  score_data <- intersect(test_data$genes,Fig4A_data$genes)
  length_test <- length(score_data)
  AST_data <- full_join(test_data,AML_signature,by="genes")
  length_AST <- nrow(AST_data)
  length_AS <- nrow(AML_signature)
  length_AL <- nrow(Fig4A_data)
  EnSc_fixed_list[[test]] <- (length_AST-((length_test*length_AS)/length_AL))/length_test
}

test_data<-AML_signature
score_data <- intersect(test_data$genes,Fig4A_data$genes)
length_test <- length(score_data)
AST_data <- full_join(test_data,AML_signature,by="genes")
length_AST <- nrow(AST_data)
length_AS <- nrow(AML_signature)
length_AL <- nrow(Fig4A_data)
EnSc_fixed_list[["Aging_signiture"]]<-(length_AST-((length_test*length_AS)/length_AL))/length_test

Fig4A_data <- read.table("multidata_genematrix.txt",
                         sep = "\t",header = T,check.names = F)
AML_signature <- read.table("AML_signature_diff.txt",
                            sep = "\t",header = T,check.names = F)

for (i in 1:26) {
  AML_list <- Fig4A_data[,c(1,28)]
  EnSc_matrix <- matrix(nrow = nrow(AML_list),ncol = 3,dimnames = list(NULL,c("num","score","gene")))
  test_data <- Fig4A_data[,c(1,i+1)]%>%na.omit()
  ref_list = test_data[,2]
  names(ref_list) = test_data$genes
  
  for (group_scale in 1:nrow(AML_list)) {
    EnSc_matrix[group_scale,3] <- AML_list$genes[group_scale]
    EnSc_matrix[group_scale,1]<-group_scale
    AML_group <- AML_list[1:group_scale,]
    length_test <- length(intersect(AML_list$genes,test_data$genes))
    length_AL <- nrow(AML_list)
    length_AS <- group_scale
    score_data <- inner_join(test_data,AML_group,by="genes")
    length_AST <- nrow(score_data)
    EnSc_matrix[group_scale,2] <- (length_AST-((length_test*length_AS)/length_AL))/length_test
  }
  EnSc_matrix<-as.data.frame(EnSc_matrix)
  EnSc_matrix$inset <- rep(0, nrow(Fig4A_data))
  EnSc_matrix$inset[which(EnSc_matrix$gene %in% test_data$gene)] <- EnSc_matrix$num[which(EnSc_matrix$gene %in% test_data$gene)]
  EnSc_matrix$num<-as.numeric(EnSc_matrix$num)
  EnSc_matrix$score<-as.numeric(EnSc_matrix$score)
  EnSc_matrix$inset<-as.numeric(EnSc_matrix$inset)
  
  GSEA_plot <- ggplot(data = EnSc_matrix,aes(x = num,y = score))+
    geom_line(size = 1.2,color = "#005187")+
    geom_point(aes(x=inset,y=min(score)),shape = 108, size = 10, color = "#e51a2f")+
    labs(subtitle = paste0("Enrch.Score = ",signif(max(EnSc_matrix$score),digits = 3)),
         title = colnames(test_data)[2],x="Number of genes",y="Enrichment score")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    theme(panel.border = element_rect(fill = NA,color = "black",size = 2),
          axis.ticks = element_line(size = 1,linetype = 1),
          plot.title = element_text(size=18,face = "bold",colour = "black",vjust = 2,hjust = 0.5),
          axis.text = element_text(size=15,face = "bold",colour = "black"),
          legend.title = element_blank(),legend.key.width = unit(0.7,"cm"),
          axis.title=element_text(size=18,colour = "black"),
          plot.subtitle=element_text(size=15, color="black")
    )
  GSEA_plot
  ggsave(GSEA_plot, filename = paste0(colnames(test_data)[2],"_Enrichment_plot.PDF"), device = cairo_pdf, 
         width = 8, height = 7, units = "in")
}

library(ggplot2)
library(ggrepel)
lineplot<-openxlsx::read.xlsx("DGSE_score.xlsx",sheet = 1)
lineplot<-lineplot[order(lineplot$score,decreasing = T),]
p<-ggplot(lineplot, aes(x = score,y = value,label=data)) +
  geom_bar(mapping = aes(x = score),width = 0.002,stat = "identity",fill="#005187") +
  geom_point(aes(x=score,y=value),colour="#e3192f",size=1.5)+
  theme(legend.position = "top")+
  geom_text_repel(data=lineplot,box.padding   = 0.35, 
                  point.padding = 1,nudge_y = 0.5,
                  segment.color = 'black',max.overlaps = 10,vjust=0.2)+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=12,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks = element_line(size = 2))
ggsave(p, filename = "DGSE_score_count.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")
rank_score_matrix<-do.call(rbind,rank_score_list)
rs_line_plot<-cbind(data=rownames(rank_score_matrix),
                    score=rank_score_matrix[,1],value=rep(1,nrow(rank_score_matrix)))

rs_line_plot<-as.data.frame(rs_line_plot)
rs_line_plot[,3]<-as.numeric(rs_line_plot[,3])
rs_line_plot[,2]<-as.numeric(rs_line_plot[,2])
rs_line_plot <- rs_line_plot[-27,]
p<-ggplot(rs_line_plot, aes(x = score,y = value,label=data)) +
  ylim(0,2)+
  geom_bar(mapping = aes(x = score),width = 3,stat = "identity",fill="#005187") +
  geom_point(aes(x=score,y=value),colour="#e3192f",size=1.5)+
  theme(legend.position = "top")+
  geom_text_repel(data=rs_line_plot,box.padding   = 0.35, 
                  point.padding = 1,nudge_y = 0.5,
                  segment.color = 'black',max.overlaps = 20,vjust=0.2)+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.title = element_text(size = 15,face = "bold",colour = "black"),
        axis.text = element_text(size=12,face = "bold",colour = "black"),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks = element_line(size = 2))
ggsave(p, filename = "rank_score_count.PDF", device = cairo_pdf, 
       width = 8, height = 7, units = "in")
write.table(rs_line_plot,"rs_score.txt",sep = "\t",row.names = F,col.names = T,quote = F)

