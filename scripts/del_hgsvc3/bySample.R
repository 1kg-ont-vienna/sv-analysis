library(ggplot2)
library(reshape2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x = element_text(size = axisFontSize), strip.text.y = element_text(size = axisFontSize))


x=read.table("stats.tsv", header=T, sep="\t")
x=x[x$dataset=="1kG_ONT",]
print(summary(x))
df = melt(x[,c("size","sample","FDR","TRfraction")], id.vars=c("size","sample"))
df[df$size=="small",]$size = "<250bp"
df[df$size=="large",]$size = ">=250bp"
df$size=factor(df$size, levels=c("<250bp", ">=250bp"))

png("del.png", width=1200, height=800)
p = ggplot(data=df, aes(x=size, y=value))
p = p + geom_boxplot(aes(color=variable))
p = p + xlab("Deletion size")
p = p + ylab("Value") + ylim(0,1)
p = p + scale_fill_brewer(palette = "Set2")
p = p + scienceTheme
p = p + theme(legend.position="bottom")
p = p + labs(color="")
p
dev.off()

