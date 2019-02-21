rm(list=ls(all=TRUE)) 
setwd("/Users/shding/Dropbox/")
load("real_minus_random_means_cons_div_181211.rda")

library(ggplot2)
pdf("phyloP_differences_cons_div_Mouse.pdf", useDingbats=F)
ggplot(real_minus_random_means_cons_div, aes(x=Cons_Div, y=Diff_contig)) + 
  geom_jitter(aes( col=subclass_color), width = 0.2, size=2) + scale_color_identity()  + 
  scale_x_discrete(limits = unique(real_minus_random_means_cons_div$Cons_Div)) + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5) + 
  ylim(0,0.14)
dev.off()	

