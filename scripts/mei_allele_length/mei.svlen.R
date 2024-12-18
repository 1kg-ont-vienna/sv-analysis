library(ggplot2)
library(reshape2)
library(scales)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x = element_text(size = axisFontSize), strip.text.y = element_text(size = axisFontSize))

x=read.table("ont.mei")
y=read.table("nygc.mei")
colnames(x)=c("chr", "pos", "svlen", "family")
colnames(y)=c("chr", "pos", "svlen", "family")
x$dataset="1kG_ONT"
y$dataset="Illumina"
df = rbind(x[,c("svlen","family","dataset")], y[,c("svlen","family","dataset")])

print(summary(df))
png("mei_svlen.png", width=1200, height=800)
p = ggplot(data=df, aes(x=svlen))
p = p + geom_freqpoly(aes(color=family, linetype=dataset), linewidth=1, bins=50)
p = p + xlab("MEI length (in bp)")
p = p + ylab("Number of ascertained MEIs (in log-scale)")
p = p + scale_x_continuous(labels=comma, limits=c(0,10000))
p = p + scale_y_log10()
p = p + scale_color_brewer(palette = "Dark2")
p = p + scienceTheme
p = p + theme(legend.position="right")
p = p + labs(color="Family", linetype="Dataset")
p
dev.off()
