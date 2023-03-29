## circlize.R

library(circlize)
Args<-commandArgs()
ecc_result_path = Args[6]
gene_anno_path = Args[7]
output_pdf_path = Args[8]

##ecc_result
ecc_result = read.table(ecc_result_path,
                        sep='\t', col.names=c('chr','start','end','count','id'))
ecc_result$name = paste0(ecc_result$chr,'_',ecc_result$start,'_',ecc_result$end)
ecc_result = ecc_result[, c('name','start', 'end','chr')]

gene_anno = read.table(gene_anno_path,
                        sep='\t', col.names=c('chr','start','end','count','id', 'gchr','gstart','gend','gid'))
gene_anno$name = paste0(gene_anno$chr,'_',gene_anno$start,'_',gene_anno$end)
gene_anno = gene_anno[, c('name','gstart','gend','gid')]


# Plot All gene in a circos
pdf(output_pdf_path)
circos.genomicInitialize(ecc_result, labels.cex = 1.25, sector.names = ecc_result$id) + # Step 1 -> circos.initialize
circos.par(points.overflow.warning=F) 
circos.track(ylim = c(0, 1), 
    bg.col = c("#FF000040", "#00FF0040", "#0000FF40"),
    bg.border = NA, track.height = 0.15) # Step2 -> track&modification repeat
# print(get.all.sector.index())

circos.track(ylim = c(0,1), bg.border = NA, track.height = 0.05)

for(name in ecc_result$name){
    subanno = gene_anno[gene_anno$name == name,]
    if(dim(subanno)[1] != 0){
    for(i in 1:dim(subanno)[1]){
        gene = subanno[i,]
        circos.arrow(gene$gstart, gene$gend, # is the position of chr(begin at start of eccDNA) but not eccDNA length(which begin at 0)
                    sector.index = name, # default is the last sector plot
                    arrow.head.length = mm_x(1), arrow.head.width = 0.90,
                    col='#000000') # default: arrow.head.length = mm_x(5), arrow.head.width = width*2
        circos.text((gene$gstart+gene$gend)/2, -0.7, gene$gid, # x,y,labels
                   facing='reverse.clockwise', niceFacing =T, adj=c(0,0), cex=0.7)
#         print(gene)
    }
    }
}
# print(get.all.track.index())
dev.off()

