library(ggplot2)
library(reshape2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x = element_text(size = axisFontSize), strip.text.y = element_text(size = axisFontSize))

x = read.table("fdr.by_sample.tsv", header=T)
x = x[,c("sample","svtype","class","fdr","trf_fraction")]
df = melt(x, id.vars=c("sample", "svtype", "class"))
df[df$class=="SMALL",]$class = "<250bp"
df[df$class=="LARGE",]$class = ">=250bp"
df$class=factor(df$class, levels=c("<250bp", ">=250bp"))

png("bySample.png", width=600, height=400)
p = ggplot(data=df, aes(x=class, y=value))
p = p + geom_boxplot(aes(color=variable))
p = p + xlab("SV size")
p = p + ylab("Value") + ylim(0,1)
p = p + scienceTheme
p = p + scale_colour_manual(labels = c("FDR", "Tandem-repeat SV Fraction"), values=c("red", "blue"))
p = p + theme(legend.position="bottom")
p = p + facet_wrap(~svtype)
p = p + labs(color="")
p
dev.off()

