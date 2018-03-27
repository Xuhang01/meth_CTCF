library(ggplot2)

######## 1. donor_count plot
df_count <- read.table("./misc/stat/donor_count.txt")
df_count$proj_code <- rownames(df_count) 


count <- c(df_count$icgc_donor_id,df_count$icgc_specimen_id,df_count$icgc_sample_id)
type <- rep(c('icgc_donor_id','icgc_specimen_id','icgc_sample_id'),times=c(19,19,19))
proj_code <- rep(df_count$proj_code, times=3)

df_count <- data.frame(count, type, proj_code)
df_count$type <- factor(df_count$type, levels=c("icgc_donor_id","icgc_specimen_id","icgc_sample_id"))
ggplot(df_count, aes(x=proj_code, y=count, fill=type))+
  geom_col(position = "dodge")+
  geom_text(aes(label=count),position=position_dodge(width=1),vjust=-1, size=2, angle=30,hjust=0)+
  theme(axis.text.x = element_text(angle=50,hjust=1))
ggsave("./misc/pictures/sample_count.jpeg", width=10, height=6)


######## 2. density plot of methylation level
library(ggplot2)
density <- function(probe_id) {
  df_probe <- read.table(paste("./misc/stat/probe_id/", probe_id ,".txt",sep=''),head=T,sep="\t")
  df_ICGC <- df_probe
  df_ICGC$project_code <- "ICGC"
  df_total <- rbind(df_probe, df_ICGC)
  ggplot(df_total, aes_string(x=probe_id))+geom_density()+xlim(0,1)+facet_wrap( ~ project_code,ncol=5)
  ggsave(paste("./misc/pictures/probe_id/", probe_id ,".jpeg", sep=''),width=10,height=6)
}
probe_id_list <- array(c('cg00000029','cg00000165','cg00000236','cg00000289','cg00000292','cg00000321','cg00000363','cg00000622'))
apply(probe_id_list, 1, density)



######## 3. relation between mean and std
mean_and_std <- read.table('./misc/stat/mean_and_std.txt', header=T, sep='\t')
#### only ggplot2
library(ggplot2)
ggplot(mean_and_std,aes(x=mean_meth, y=std_meth))+geom_point(  size=2, colour='red',shape=16, alpha=0.6)
ggsave('./misc/pictures/mean_and_std.jpeg',width=10, height=6)
#### gridExtra
library(gridExtra)
hist_top <- ggplot(mean_and_std, aes(mean_meth))+geom_histogram()
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter <- ggplot(mean_and_std,aes(x=mean_meth, y=std_meth)) + 
  geom_point( size=1, colour='red',shape=16, alpha=0.6)
hist_right <- ggplot(mean_and_std, aes(std_meth))+geom_histogram()+coord_flip()
g <- arrangeGrob(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
ggsave('./misc/pictures/mean_and_std.gridExtra.jpeg',width=10, height=6, g)
#### ggExtra
library(ggExtra)
scatter <- ggplot(mean_and_std,aes(x=mean_meth, y=std_meth)) + 
  geom_point( size=1, colour='red',shape=16, alpha=0.6) + 
  theme_classic()
# ggExtra::ggMarginal(scatter, type="histogram")
ggsave('./misc/pictures/mean_and_std.ggExtra.jpeg',ggMarginal(scatter, type="histogram"),width=10, height=6)



######## 4. probe related to Island or Not (Kolmogorov-Smirnov test show there are difference between them)
library(ggplot2)
Island <- read.table("./misc/stat/Island_and_non_Island_comparison.txt",header=T,sep="\t")
# case_count comparsion
ggplot(Island, aes(x=is_Island, y=case_count))+geom_boxplot()+geom_jitter(alpha=0.5)
ggsave("./misc/pictures/case_count.compare_Island_non-Island.jpeg",width=8,height=8)
ks.test(subset(Island, is_Island=='Island')$case_count, subset(Island, is_Island != 'Island')$case_count)
#D = 0.083632, p-value = 0.005762
# mean_meth comparison
ggplot(Island, aes(x=is_Island, y=mean_meth))+geom_boxplot()+geom_jitter(alpha=0.4)
ks.test(subset(Island, is_Island=='Island')$mean_meth, subset(Island, is_Island != 'Island')$mean_meth)
ggsave("./misc/pictures/mean_meth.compare_Island_non-Island.jpeg",width=8,height=8)
# std_meth comparison
ggplot(Island, aes(x=is_Island, y=std_meth))+geom_boxplot()+geom_jitter(alpha=0.5)
ks.test(subset(Island, is_Island=='Island')$std_meth, subset(Island, is_Island != 'Island')$std_meth)
ggsave("./misc/pictures/std_meth.compare_Island_non-Island.jpeg",width=8,height=8)
# mean std relation according to the type of probe
hist_top <- ggplot(Island, aes(mean_meth))+geom_density(aes(colour=is_Island))+ theme(legend.position="none")
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter <- ggplot(Island, aes(x=mean_meth, y=std_meth))+geom_point(aes(colour=is_Island))
hist_left <- ggplot(Island, aes(std_meth))+geom_density(aes(colour=is_Island))+coord_flip()+ theme(legend.position="none")
g <- arrangeGrob(empty,hist_top, hist_right, scatter, ncol=2, nrow=2, widths=c(1,4), heights=c(1,4))
ggsave('./misc/pictures/mean_and_std.Island.jpeg',width=10, height=6, g)

######## position between probe and TSS
