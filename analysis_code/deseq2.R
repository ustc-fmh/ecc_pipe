library(DESeq2)
Args<-commandArgs()
countdir = Args[6]
metadir = Args[7]

output_norm = paste0(Args[8],'/deseq2_norm_matrix.csv')
output_result = paste0(Args[8],'/deseq2_result.csv')

readcount=read.table(countdir,sep="\t", row.names=1, header=T)
coldata<-read.table(metadir, sep="\t", row.names=2, header=F) ## F
colnames(coldata) <- c('sample', 'tool')
coldata<-coldata['sample']

#colnames(readcount) <- rownames(coldata)
all(rownames(coldata) %in% colnames(readcount))
all(rownames(coldata) == colnames(readcount))
dds <- DESeqDataSetFromMatrix(readcount, coldata, design= ~ sample)

dds <- dds[rowSums(counts(dds)) > 4,]

dds2 <- DESeq(dds)
#res <- results(dds2)
res <- results(dds2, contrast=c("sample", unique(coldata[,'sample'])))
resOrdered <- res[order(res$pvalue),]
res05 <- results(dds2, alpha=0.05)
summary(res05)

normcount<-counts(dds2,normalize=TRUE)
write.table(normcount,output_norm,sep=",",row.names = T,quote = F)

log_norcount<-log2(normcount+1)
resdata <-  merge(as.data.frame(resOrdered),as.data.frame(log_norcount),by="row.names",sort=FALSE)
write.table(resdata,output_result,sep=",",row.names = F,quote = F)