library(ggplot2)
library(reshape2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x = element_text(size = axisFontSize), strip.text.y = element_text(size = axisFontSize))


x=read.table("stats.tsv", header=T)
x$sensitivity = x$shared / x$all
df = melt(x[,c("mei","sample","dataset","shared","unique")], id.vars=c("mei","sample","dataset"))
df$variable = factor(df$variable, levels=c("shared", "unique"))
df = df[order(c(df$variable), decreasing = T),]
print(df)

png("mei.png", width=1200, height=800)
p = ggplot(data=df, aes(x=sample, y=value))
p = p + geom_bar(aes(fill=variable), stat="identity", color="black")
p = p + xlab("Shared Samples")
p = p + ylab("Number of ascertained MEIs")
p = p + scale_fill_brewer(palette = "Set2")
p = p + facet_wrap(dataset~mei, scales='free')
p = p + scienceTheme
p = p + theme(legend.position="bottom")
p = p + labs(fill="MEI comparison")
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
dev.off()

s = read.table("samples.info", header=T, sep="\t", comment.char="$")
s$COVBIN="<=20x"
s[s$COVERAGE>20,]$COVBIN=">20x"
s$READN50BIN="<=20Kbp"
s[s$N50>20000,]$READN50BIN=">20Kbp"
df = merge(x, s)
print(df)


df = df[df$dataset=="HGSVC3",]
png("params.png", width=1200, height=800)
p = ggplot(data=df, aes(x=READN50BIN, y=sensitivity))
p = p + geom_boxplot(aes(fill=COVBIN), color="black")
p = p + xlab("N50 read length")
p = p + ylab("Shared fraction of MEIs\nwith whole-genome assemblies")
p = p + scale_fill_brewer(palette = "Dark2")
p = p + facet_wrap(~mei)
p = p + scienceTheme
p = p + theme(legend.position="bottom")
p = p + labs(fill="Coverage")
p
dev.off()
