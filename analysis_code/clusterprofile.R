library(DOSE)
library(dplyr)
library(patchwork)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(enrichplot)
library(ggplot2)
library(cowplot)
options(repr.plot.width=15,repr.plot.height=10)

Args<-commandArgs()
deg_path = Args[6]
output_path = Args[7]
pvalue_cut = as.numeric(Args[8])
fc_cut = as.numeric(Args[9])
geno = Args[10]
mode = Args[11]


GO_up = paste0(output_path,'/', mode, '_GO_up.csv')
#GO_up_pdf = paste0(output_path,'/GO_up.pdf')

GO_down = paste0(output_path,'/', mode,'_GO_down.csv')
#GO_down_pdf = paste0(output_path,'/GO_down.pdf')

KEGG_up = paste0(output_path,'/', mode,'_KEGG_up.csv')
#KEGG_up_pdf = paste0(output_path,'/KEGG_up.pdf')

KEGG_down = paste0(output_path,'/', mode,'_KEGG_down.csv')
#KEGG_down_pdf = paste0(output_path,'/KEGG_down.pdf')

DEG_result <- read.table(deg_path,
                         header = TRUE, sep = ",",row.names=1)

##go
go <- function(DEG_result, ref){
    ## 转 ensg
    if (ref == 'mm10'){
        ref_org = org.Mm.eg.db
    }
    else{
        ref_org = org.Hs.eg.db
    }
    
    gene.df <- bitr(rownames(DEG_result), fromType = "SYMBOL",
        toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
        OrgDb = ref_org)
    ego <- enrichGO(gene          = gene.df$ENSEMBL,
                keyType = "ENSEMBL",
                OrgDb         = ref_org,
                ont           = "BP", ##全部查看
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
        readable      = TRUE)
    
    return (ego)
}

kegg <- function(DEG_result, ref){
    ## 转 ensg
    if (ref == 'mm10'){
        ref_org = org.Mm.eg.db
        ref_organ = 'mmu'
    }
    else{
        ref_org = org.Hs.eg.db
        ref_organ = 'hsa'
    }
    
    gene.df <- bitr(rownames(DEG_result), fromType = "SYMBOL",
        toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
        OrgDb = ref_org)
    
    kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = ref_organ, #mmu 为鼠， hsa为人
                 pvalueCutoff = 1)
    
    return (kk)
}
up = DEG_result[which(DEG_result['log2FoldChange'] > fc_cut),]
up = up[which(up['pvalue'] < pvalue_cut),]
down = DEG_result[which(DEG_result['log2FoldChange'] < -fc_cut),]
down = down[which(down['pvalue'] < pvalue_cut),]

ego_up <- go(up, geno)
ego_down <- go(down, geno)

if (length(ego_up) != 0){
    write.table(head(ego_up, 20), file = GO_up,  sep = ",", row.names = TRUE,
            col.names = TRUE)
    #pdf(file =GO_up_pdf)
    #dotplot(ego_up, x = "GeneRatio", color = 'p.adjust', showCategory = 5, title="GO_UP") #点图，按富集的数从大到小的
    #dev.off()
}

if (length(ego_down) != 0){
    write.table(head(ego_down, 20), file = GO_down,  sep = ",", row.names = TRUE,
            col.names = TRUE)
    #pdf(file =GO_down_pdf)
    #dotplot(ego_down, x = "GeneRatio", color = 'p.adjust', showCategory = 5, title="GO_DOWN") #点图，按富集的数从大到小的
    #dev.off()
}

# kk_up <- kegg(up, geno)
# kk_down <- kegg(down, geno)
# if (length(kk_up) != 0){
#     write.table(head(kk_up, 20), file = KEGG_up,  sep = ",", row.names = TRUE,
#             col.names = TRUE)
#     #pdf(file =KEGG_up_pdf)
#     #dotplot(kk_up, x = "GeneRatio", color = 'p.adjust',  title="KEGG_UP") #点图，按富集的数从大到小的
#     #dev.off()
# }

# if (length(kk_down) != 0){
#     write.table(head(kk_down, 20), file = KEGG_down,  sep = ",", row.names = TRUE,
#             col.names = TRUE)
#     #pdf(file =KEGG_down_pdf)
#     #dotplot(kk_down, x = "GeneRatio", color = 'p.adjust',  title="KEGG_DOWN") #点图，按富集的数从大到小的
#     #dev.off()
# }

## GSEA 
marker = DEG_result
marker['SYMBOL'] <- row.names(marker)
marker = marker[,c('pvalue', 'log2FoldChange', 'SYMBOL')]
colnames(marker) <- c('p_val', 'avg_log2FC', 'SYMBOL')

data = marker
data_up <- subset(data, avg_log2FC>0.15)
data_down <- subset(data, avg_log2FC< -0.15)
data <- rbind(data_up, data_down)
data <- data %>% arrange(desc(avg_log2FC))

gene <- data$SYMBOL

if (geno == 'mm10'){
    gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 
}else{
    gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
}
 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
data_all <- data %>% 
 inner_join(gene,by="SYMBOL")
#dim(data_all)
data_all_sort <- data_all %>% arrange(desc(avg_log2FC))
geneList = data_all_sort$avg_log2FC 
names(geneList) <- data_all_sort$ENTREZID 

if (geno == 'mm10'){
    gse.GO <- gseGO(
         geneList, #geneList
         minGSSize=1,
         ont = "BP",  # 可选"BP"、"MF"和"CC"或"ALL"
         OrgDb = org.Mm.eg.db, #人 注释基因
         keyType = "ENTREZID",
         pvalueCutoff = 0.99,
         pAdjustMethod = "BH",#p值校正方法，
            exponent = 0.5,
    )
}else{
    gse.GO <- gseGO(
         geneList, #geneList
         minGSSize=1,
         ont = "BP",  # 可选"BP"、"MF"和"CC"或"ALL"
         OrgDb = org.Hs.eg.db, #人 注释基因
         keyType = "ENTREZID",
         pvalueCutoff = 0.99,
         pAdjustMethod = "BH",#p值校正方法，
            exponent = 0.5,
    )
}

gsea_df <- data.frame(gse.GO)
gsea_result_1 = paste0(output_path,'/', mode, '_gsea_result_1.csv')
gsea_result_2 = paste0(output_path,'/', mode, '_gsea_result_2.csv')

write.csv(gsea_df, file=gsea_result_1,
                  sep = ",", row.names = TRUE, col.names = TRUE)
write.csv(data_all_sort, file=gsea_result_2, 
         sep = ",", row.names = TRUE, col.names = TRUE)
## python annotation

## default save top 5 GSEA and all result
for(i in c(1,2,3,4,5)){
    top_go = gsea_df[i,'ID']
    pdf_path = paste0(output_path,'/', mode, '_gsea_top_', i,'.pdf')
    #print(top_go)
    #print(pdf_path)
    p1 <- gseaplot2(gse.GO,
         title = "",  #设置title
         top_go, 
         color="red", #线条颜色
         base_size = 18, #基础字体的大小
         subplots = 1:2, #展示上2部分
         pvalue_table = T) # 显示p值
    
    pdf(pdf_path, width=10, height=6)
    print(p1)
    dev.off()
}

