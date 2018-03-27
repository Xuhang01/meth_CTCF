require(ggplot2)
args <- commandArgs(T)
probe_id <- args[1]
gene_id <- args[2]
beta <- as.numeric(args[3])
qvalue <- as.numeric(args[4])
df_pairs <- read.table(paste('./data/eQTL_input_icgc_donor_id/significant_eQTL_pairs/',gene_id,'_',probe_id,'.txt',sep=""),
                       header=T, sep="\t")
# ggplot(df_pairs, aes_string(x=probe_id,y=gene_id))+geom_point() + geom_smooth(method='lm',formula=y~x)
# ggsave(paste("./misc/pictures/gene_probe_pairs/",gene_id,"_",probe_id,".png",sep=""))


library(gridExtra)
hist_top <- ggplot(df_pairs, aes_string(probe_id))+geom_histogram(bins=30)
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter <- ggplot(df_pairs, aes_string(x=probe_id,y=gene_id))+geom_point() + geom_smooth(method='lm',formula=y~x) +
	scale_x_continuous(limits = c(0, 1)) + 
	annotate("text", x = 0.85, y = max(df_pairs[gene_id])*0.95, label = paste("beta=", format(beta,digits=2), "\n", "qvalue=", format(qvalue,digits=2)))
hist_right <- ggplot(df_pairs, aes_string(gene_id))+geom_histogram(bins=30) + coord_flip()
# grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))

g <- arrangeGrob(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
ggsave(paste("./misc/pictures/gene_probe_pairs/",gene_id,"_",probe_id,".png",sep=""), width=10, height=6, g)

