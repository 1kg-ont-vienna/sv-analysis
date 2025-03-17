library(ggplot2)
library(scales)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisFontSize))

args=commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=F)
colnames(x)=c("id", "type", "size", "maf", "ac", "status")
x$type = factor(x$type)
x$status = factor(x$status)
print(summary(x))

png(paste0(args[1], ".size.png"), width=1000, height=400)
p = ggplot(data=x, aes(x=size))
p = p + geom_freqpoly(aes(group=status, color=status))
p = p + scale_x_log10(labels=comma, limits=c(min(x$size), max(x$size))) + facet_wrap(~type)
#p = p + scale_y_log10()
p = p + xlab("SV size (in bp)")
p = p + ylab("Count")
#p = p + ggtitle(args[1])
p = p + scienceTheme
for(i in c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000)) { p = p + geom_hline(yintercept=i, linetype="dashed", color="lightgrey"); }
p
dev.off()

png(paste0(args[1], ".maf.png"), width=1000, height=400)
p = ggplot(data=x, aes(x=maf))
p = p + geom_freqpoly(aes(group=status, color=status))
p = p + scale_x_continuous(labels=comma, limits=c(min(x$maf), 0.5)) + facet_wrap(~type)
#p = p + scale_y_log10()
p = p + xlab("Minor Allele Frequency (MAF)")
p = p + ylab("Count")
#p = p + ggtitle(args[1])
p = p + scienceTheme
for(i in c(1000, 2000, 3000, 4000, 5000)) { p = p + geom_hline(yintercept=i, linetype="dashed", color="lightgrey"); }
p
dev.off()
