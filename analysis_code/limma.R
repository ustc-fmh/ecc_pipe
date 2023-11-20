library('limma')
library('edgeR')
Args<-commandArgs()
countdir = Args[6]
metadir = Args[7]

output_norm = paste0(Args[8],'/limma_norm_matrix.csv')
output_result = paste0(Args[8],'/limma_result.csv')

readcount=read.table(countdir,sep="\t", row.names=1, header=T)
coldata<-read.table(metadir, sep="\t", row.names=2, header=F) ## F
colnames(coldata) <- c('sample', 'tool')

##add
coldata<-coldata[,'sample']
exp=unique(coldata)[1]
ctr=unique(coldata)[2]
#分组矩阵design构建
design <- model.matrix(~0+factor(coldata)) 
colnames(design) <- levels(factor(coldata))
rownames(design) <- colnames(readcount)
## 表达矩阵DGEList构建与过滤低表达基因
dge <- DGEList(counts=readcount) 
# keep.exprs <- filterByExpr(dge,design=design) #过滤低表达基因
# dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 

dge <- calcNormFactors(dge)  #归一化基因表达分布,得到的归一化系数被用作文库大小的缩放系数
cont.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr), #比对顺序实验/对照
                             levels = design)
de <- voom(dge,design,plot=TRUE, normalize="quantile")
fit1 <- lmFit(de, design)               #线性拟合
fit2 <- contrasts.fit(fit1,cont.matrix) #统计检验
efit <- eBayes(fit2, trend=F)  #Apply empirical Bayes smoothing to the standard errors

tempDEG <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)  #padj值从小到大排列
DEG_limma_voom  <- na.omit(tempDEG)
colnames(DEG_limma_voom) = c('log2FoldChange', 'AveExpr', 't', 'pvalue', 'adj.P.Val', 'B')

write.table(de$E,output_norm,sep=",",row.names = T,quote = F)

resdata <-  merge(as.data.frame(DEG_limma_voom),as.data.frame(de$E),by="row.names",sort=FALSE)
write.table(resdata, output_result,sep=",",row.names = F,quote = F)
