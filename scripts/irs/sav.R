library(ggplot2)
library(scales)
library(grid)
library(reshape2)


cs = read.table("DEL.affy.pval.gz", header=T)

txtFontSize=10
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), plot.title=element_text(size=axisTtlFontSize))

# FDR by vaf and size
svType = substr(cs, 0, 3)
colNames=c("LOWERPVALUE")
cs$size = cs$END - cs$START

png("irs.png", width=600, height=600)
p = ggplot(data=cs, aes(x=LOWERPVALUE))
p = p + geom_histogram(bins=50)
p = p + scienceTheme
p
dev.off()

png("irs2.png", width=600, height=600)
p = ggplot(data=cs, aes(x=LOWERPVALUE, y=size))
p = p + geom_point(alpha=0.1)
p = p + scienceTheme
p = p + ylab("SV size")
p
dev.off()

# Overall FDR
nSiteDel=nrow(cs)
cs = cs[!is.na(cs$LOWERPVALUE),]
nEvalDel=nrow(cs)
nHighDel=sum(cs$LOWERPVALUE>0.5)
fdrDel=0
if (nEvalDel>0) { fdrDel=round((2*nHighDel)/nEvalDel, digits=4); }

print(paste0("Deletions (#n=", nSiteDel, ")"))
print(paste0("FDR(DEL)=", fdrDel, " (#n=", nEvalDel, ")"))
print(warnings())

