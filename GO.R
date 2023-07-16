library("topGO")
library("Rgraphviz")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("stringr")
df<-readxl::read_excel("targets2.xlsx") #target基因列表
while(count <= 10){
  for (i in 1:length(rownames(df))){
    gene <- as.character(df[i, "gene"])
    gene_vector <- unlist(strsplit(gene, ","))
    if (length(gene_vector) <= 3){
      next
    }else{
      gene <- bitr(gene_vector[1:36],
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
      
      ego <-enrichGO(
        gene = gene$ENTREZID,   ##ENTREZID
        keyType = "ENTREZID", ##输入基因类型
        OrgDb   = org.Hs.eg.db,##导入背景基因
        ont     ="all",        ##GO的种类，BP,CC,MF
        pAdjustMethod = "BH",   ##矫正的类型
        pvalueCutoff   = 0.05,  ##P值过滤值
        readable       = TRUE)
      barplot(ego,
              drop = TRUE,
              showCategory = 6, ##显示前10个GO term
              split="ONTOLOGY")+
        scale_y_discrete(labels=function(x) str_wrap(x,width = 80))+
        facet_grid(ONTOLOGY~.,scale = 'free')
      
      ggsave(paste("motif",".png",sep=as.character(i)),height = 5,width = 10)
      count = count +1
    }
  }
}
