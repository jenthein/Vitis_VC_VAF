library(ggplot2)
library(dplyr)

########################################################
## jens.theine@uni-bielefeld.de / V:07.05.2024
########################################################

################### config #############################

file = "calls_q30filteredvariants_biSNPs_AN__Chr_Pos_DP_VAF.tsv"
out = "/out_all/"

contigs = list("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrUn")
contig_of_interest = "chr14"

anno = "__VAF_Cov_plot"

########################################################


all <- read.csv(file, header=FALSE, sep = "\t", quote="")
colnames(all) <- c("tig", "pos", "valCov", "valVAF")

for(contig_of_interest in contigs){

df <- all[grep(contig_of_interest,all$tig),]
##df <- filter(df, between(pos, START, END))

df$valVAF = df$valVAF*100

fulloutpath = paste(out, contig_of_interest, anno, sep = "")


density<-ggplot()+ 
        geom_point(data=df, aes(x=pos, y=valVAF ), color="darkblue", size=0.001, alpha=0.6) +
        geom_point(data=df, aes(x=pos, y=valCov), color="darkred", size=0.1, alpha=0.6)+
        scale_y_continuous(limits=c(0, 100), sec.axis = sec_axis(~ ./100, name = "variant allele frequency (VAF, blue)"))+
        labs(x=paste0(contig_of_interest," [bp]"), y="read depth (DP, red)")

print(density)
ggsave(paste0(fulloutpath, "_alldata.pdf"), width=16, height=9, dpi = 150)
ggsave(paste0(fulloutpath, "_alldata.png"), width=16, height=9, dpi = 150)

densityAVG<-ggplot()+ 
	geom_smooth(data=df, aes(x=pos, y=valVAF ), colour="darkblue", size=1, method="loess", span=0.01, se=FALSE) +
	geom_smooth(data=df, aes(x=pos, y=valCov ), colour="darkred", size=1, method="loess", span=0.01, se=FALSE) +
        scale_y_continuous(limits=c(0, 100), sec.axis = sec_axis(~ ./100, name = "avg. VAF (blue)"))+
        labs(x=paste0(contig_of_interest," [bp]"), y="avg. DP (red)")


print(densityAVG)
ggsave(paste0(fulloutpath, "_avgdata.pdf"), width=16, height=9, dpi = 150)
ggsave(paste0(fulloutpath, "_avgdata.png"), width=16, height=9, dpi = 150)

densityallinone<-ggplot()+ 
        geom_point(data=df, aes(x=pos, y=valVAF ), color="darkblue", size=0.001, alpha=0.6) +
        geom_point(data=df, aes(x=pos, y=valCov), color="darkred", size=0.1, alpha=0.6)+
	geom_smooth(data=df, aes(x=pos, y=valVAF ), colour="blue", size=1, method="loess", span=0.01, se=FALSE, alpha=0.8) +
	geom_smooth(data=df, aes(x=pos, y=valCov ), colour="red", size=1, method="loess", span=0.01, se=FALSE, alpha=0.8) +
        scale_y_continuous(limits=c(0, 100), sec.axis = sec_axis(~ ./100, name = "variant allele frequency (VAF, blue)"))+
        labs(x=paste0(contig_of_interest," [bp]"), y="read depth (DP, red)")

print(densityallinone)
ggsave(paste0(fulloutpath, "_allinone.pdf"), width=16, height=9, dpi = 150)
ggsave(paste0(fulloutpath, "_allinone.png"), width=16, height=9, dpi = 150)

dev.off()

}
