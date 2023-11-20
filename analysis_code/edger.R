library('edgeR')
Args<-commandArgs()
countdir = Args[6]
metadir = Args[7]

output_norm = paste0(Args[8],'/edger_norm_matrix.csv')
output_result = paste0(Args[8],'/edger_result.csv')

readcount=read.table(countdir,sep="\t", row.names=1, header=T)
coldata<-read.table(metadir, sep="\t", row.names=2, header=F) ## F
colnames(coldata) <- c('sample', 'tool')
coldata<-coldata['sample']

## add 
group = coldata[,'sample']
y<-DGEList(counts=readcount, group=group)# 构建基因表达列表
#dim(y)
y<-calcNormFactors(y, method="TMM")#  method="RLE" method="upperquartile" 'TMM'
y<-estimateCommonDisp(y)
y<-estimateTagwiseDisp(y)
et<-exactTest(y)
colnames(et) = c('log2FoldChange', 'log2CPM', 'pvalue')

write.table(y$pseudo.counts,output_norm,sep=",",row.names = T,quote = F)

resdata <-  merge(as.data.frame(et),as.data.frame(y$pseudo.counts),by="row.names",sort=FALSE)
write.table(resdata,output_result,sep=",",row.names = F,quote = F)
