library(ggplot2)
library(reshape2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x = element_text(size = axisFontSize), strip.text.y = element_text(size = axisFontSize))


x=read.table("ext.tsv", header=F)
colnames(x) = c("ext", "dataset", "all", "shared", "unique", "novelty")
print(summary(x))
df = melt(x[, c("ext", "dataset", "shared", "unique")], id.vars=c("ext", "dataset"))
df$variable = factor(df$variable, levels=c("shared", "unique"))
df = df[order(c(df$variable), decreasing = T),]
print(df)

x = x[x$dataset=="1kG_ONT",]
png("novelty.png")
p = ggplot(data=x, aes(x=ext, y=novelty))
p = p + geom_line(color="black")
p = p + xlab("Basepair extension of insertion integration point")
p = p + ylab("Fraction of novel insertions compared to gnomAD")
p
dev.off()

