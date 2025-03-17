library(karyoploteR)

#chrNames = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22")
chrNames = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")

x=read.table("in.table", header=F)
colnames(x)=c("chr", "start", "end")
x=x[x$chr %in% chrNames,]
x$chr=factor(x$chr, levels=chrNames)
sv=makeGRangesFromDataFrame(x)

cyto = read.table("chm13v2.0_cytobands_allchrs.bed.gz", header=F)
colnames(cyto) = c("chr", "start", "end", "name", "gieStain")
cyto=cyto[cyto$chr %in% chrNames,]
cyto$chr=factor(cyto$chr, levels=chrNames)
cyto=makeGRangesFromDataFrame(cyto, keep.extra.columns=TRUE)

png("fn.density.png", width=600, height=800)
kp = plotKaryotype(genome="T2T.CHM13v2.0", plot.type=1, main="FN Density", cytobands=cyto, chromosomes=chrNames)
kp = kpPlotDensity(kp, sv)
dev.off()
