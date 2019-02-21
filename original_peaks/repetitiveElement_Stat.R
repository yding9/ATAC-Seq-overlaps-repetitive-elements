rm(list=ls(all=TRUE)) 
setwd("C:/Users/yi.ding/Desktop/Peaks_JM")
x = read.table("overlap_summary.txt", sep = " ", header = F)
percentages_x = x[,1] / x[,3] * 100
x = cbind(x, percentages_x)
overlap_matrix = matrix(percentages_x, nrow = 68, ncol =7)
overlap_matrix_colnames = c("Repetitive.elements", "LINE", "LTR", "Satellite", "Simple.repeat", "SINE", "DNA.transposon")
colnames(overlap_matrix) = overlap_matrix_colnames
rownames = readLines("rowName_R_68rows.txt")
rownames_overlapMatrix = as.vector(rownames)
rownames(overlap_matrix) = rownames_overlapMatrix
#heatmap(overlap_matrix)


x.top9832 = read.table("overlap_top9832_summary.txt", sep = " ", header = F)
percent_x.top9832 = x.top9832[,1] / x.top9832[,3] * 100
x.top9832 = cbind(x.top9832, percent_x.top9832)
top9832_overlap_matrix = matrix(percent_x.top9832, nrow = 68, ncol =7)
top9832_overlap_matrix_colnames = c("Top_Repetitive.elements", "Top_LINE", "Top_LTR", "Top_Satellite", "Top_Simple.repeat", "Top_SINE", "Top_DNA.transposon")
colnames(top9832_overlap_matrix) = top9832_overlap_matrix_colnames
rownames(top9832_overlap_matrix) = rownames_overlapMatrix
overlap_includeTop9832_matrix = cbind(overlap_matrix, top9832_overlap_matrix)
#write.table(overlap_includeTop9832_matrix, file = "overlap_matrix_includeTop9832peaks.xls", sep = "\t")

i = 1
coef_Allvs.Subset = numeric(0)
while (i < length(overlap_matrix_colnames) + 1){
  res = cor(overlap_includeTop9832_matrix[,overlap_matrix_colnames[i]],overlap_includeTop9832_matrix[, top9832_overlap_matrix_colnames[i]], method = "pearson" )
  coef_Allvs.Subset[[i]] = res
  i = i + 1
}
names(coef_Allvs.Subset) = overlap_matrix_colnames
coef_Allvs.Subset




#Install the latest version from GitHub as follow (recommended):
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library("ggpubr")
#hist(x$percentages)
#hist(x$percentages, breaks = 50, xlim = c(0,100))
#hist(x$percentages[grep("LINE", x$V2)], breaks = 50)
#hist(x$percentages[grep("LINE", x$V2)], breaks = 10, xlim = c(0, 100))
ggscatter(as.data.frame(overlap_includeTop9832_matrix), x = "Repetitive.elements", y = "Top_Repetitive.elements", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Overlap with Repetitive Elements", ylab = "Subset peaks overlap with Repetitive Elements")


# integrate mean and standard devitation from bootstrap results. The following matrix is calculated using std_mean_all.R script (written by myself).
bootstrap_matrix = read.table("STD_MEAN_includeShuffle_matrix.xls", sep = "\t")
overlap_matrix_AllInOne = cbind(overlap_includeTop9832_matrix, bootstrap_matrix)
View(overlap_matrix_AllInOne)
overlap_matrix_AllInOne[,c("Repetitive.elements", "repeatMaskerintsxns100xCounts_average")]
write.table(overlap_matrix_AllInOne, file = "overlap_matrix_AllInOne.xls", sep = "\t")
