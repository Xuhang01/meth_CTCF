require(ggplot2)

df_cis_eQTL_0 <- read.table("./data/eQTL_input_icgc_donor_idresults/results/0_cis_eQTL.txt",sep="\t",head=T)


qqnorm(df_cis_eQTL_0$p.value)

df$threshold <- with(df, ifelse(
  -log10(FDR)>100 & abs(beta) >0.01, "Green", ifelse(
    -log10(FDR)>100 & abs(beta) <0.01, "Red", ifelse(
      -log10(FDR)<100 & abs(beta) >0.01, "Orange", "Black"))))

ggplot(data=df, aes(x=beta, y=-log10(p.value))) + 
  geom_point(aes(colour = threshold))+ 
  scale_colour_manual(values=c("Green"="Green", "Red"="Red","Orange"="Orange","Black"="Black"))
