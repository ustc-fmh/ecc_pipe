library(DOSE)
library(dplyr)
library(patchwork)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(enrichplot)
options(repr.plot.width=15,repr.plot.height=10)

Args<-commandArgs()
deg_path = Args[6]
output_path = Args[7]
pvalue_cut = as.numeric(Args[8])
fc_cut = as.numeric(Args[9])
geno = Args[10]


GO_up = paste0(output_path,'/GO_up.csv')
# GO_up_pdf = paste0(output_path,'/GO_up.pdf')

GO_down = paste0(output_path,'/GO_down.csv')
# GO_down_pdf = paste0(output_path,'/GO_down.pdf')

KEGG_up = paste0(output_path,'/KEGG_up.csv')
# KEGG_up_pdf = paste0(output_path,'/KEGG_up.pdf')

KEGG_down = paste0(output_path,'/KEGG_down.csv')
# KEGG_down_pdf = paste0(output_path,'/KEGG_down.pdf')

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
                ont           = "ALL", ##全部查看
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
kk_up <- kegg(up, geno)
kk_down <- kegg(down, geno)

if (length(ego_up) != 0){
    write.table(head(ego_up, 20), file = GO_up,  sep = ",", row.names = TRUE,
            col.names = TRUE)
    # pdf(file =GO_up_pdf)
    # dotplot(ego_up, x = "GeneRatio", color = 'p.adjust',  title="GO_UP") #点图，按富集的数从大到小的
    # dev.off()
}

if (length(ego_down) != 0){
    write.table(head(ego_down, 20), file = GO_down,  sep = ",", row.names = TRUE,
            col.names = TRUE)
    # pdf(file =GO_down_pdf)
    # dotplot(ego_down, x = "GeneRatio", color = 'p.adjust',  title="GO_DOWN") #点图，按富集的数从大到小的
    # dev.off()
}

if (length(kk_up) != 0){
    write.table(head(kk_up, 20), file = KEGG_up,  sep = ",", row.names = TRUE,
            col.names = TRUE)
    # pdf(file =KEGG_up_pdf)
    # dotplot(kk_up, x = "GeneRatio", color = 'p.adjust',  title="KEGG_UP") #点图，按富集的数从大到小的
    # dev.off()
}

if (length(kk_down) != 0){
    write.table(head(kk_down, 20), file = KEGG_down,  sep = ",", row.names = TRUE,
            col.names = TRUE)
    # pdf(file =KEGG_down_pdf)
    # dotplot(kk_down, x = "GeneRatio", color = 'p.adjust',  title="KEGG_DOWN") #点图，按富集的数从大到小的
    # dev.off()
}