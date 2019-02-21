#setwd("C:/Users/yi.ding/Desktop/Peaks_JM")
rm(list=ls())
setwd("/Users/shding/Dropbox/Peaks_JM")
library(reshape)
library(tibble)
library(dplyr)
library(stringr)
#read in data (had to pre-format it so it would work like this, this pre-format was done in excel by JohnM)
a = read.table("overlap_matrix_AllInOne_includeShuffle.tsv", header=T)
av_cols = c(5,7,9,11,13,15,17,19,21,23,25,27,29,31)
sd_cols = c(6,8,10,12,14,16,18,20,22,24,26,28,30,32)
#split into two for avs and sds
avs = a[,c(1:4,av_cols)]
sds = a[,c(1:4,sd_cols)]

#format the data
avs_melt = melt(data.frame(avs))    # melt function comes from the Reshape Package
sds_melt = melt(data.frame(sds))
avs_melt$variable = gsub("intsxns100xCounts","",avs_melt$variable)   # gsub() function repalces all matches of a string and returns a string vector of the same length. Unmatched substrings will returned unchanged.
sds_melt$variable = gsub("intsxns100xCounts","",sds_melt$variable)  # gsub(pattern, replacement, x, ignore.case = FALSE). x is a string or string vector.

#split into lists and reformat those lists
peakset_list_av = split(avs_melt,list(avs_melt$Peakset_type, avs_melt$variable, avs_melt$Species))
peakset_list_sd = split(sds_melt,list(sds_melt$Peakset_type, sds_melt$variable, sds_melt$Species))

####################################################################################################################################################################################
# after changing input file from overlap_matrix_AllInOne.tsv to overlap_matrix_AllInOne_includeShuffle.tsv, there are bugs in the following code. 
# The regular expression _.+$ doesn't work on shuffled_* groups. Need to fix the codes below.
####################################################################################################################################################################################
for (x in seq_along(peakset_list_av)) { # seq_along(list/vector): will get the length of the list and generate a vector from 1: length(list/vector)
  names(peakset_list_av[[x]])[6] = gsub("_","",str_match(peakset_list_av[[x]]$variable[1], regex("_.+$"))); # (peakset_list_av[[1]])[6]: the last column "value".
  # str_match(string, pattern): string is input vector (either a character vector or somethig coercible to one). Pattern uses regular expression. 
  # gsub: replace the "_" with "" (empty) in the output from str_match function. str_match DNATransposon_average, regex: regular expression. 
  #_.+$: _ dash symbol in the vector. . any character except \n. + matches at least 1 time. $ End of the string.
  peakset_list_av[[x]]$variable = gsub(str_match(peakset_list_av[[x]]$variable[1], regex("_.+$")),"",peakset_list_av[[x]]$variable);
}	
for (x in seq_along(peakset_list_sd)) { 
  names(peakset_list_sd[[x]])[6] = gsub("_","",str_match(peakset_list_sd[[x]]$variable[1], regex("_.+$")));
  peakset_list_sd[[x]]$variable = gsub(str_match(peakset_list_sd[[x]]$variable[1], regex("_.+$")),"",peakset_list_sd[[x]]$variable);
}

#now combine those lists back together and name them appropriately
# lapply(x, FUN, ...): return a list of the same length as x, each element of which is the result of applying FUN to the corresponding element of x. 
# So in the following case, peak_list_av is a list. lapply loops through each element in the list and apply as.tibble function to convert every element as tibble data type. 
peakset_list_av = lapply(peakset_list_av, as.tibble)
peakset_list_sd = lapply(peakset_list_sd, as.tibble)
peakset_list_comb = list()
# join two tibbles(tbls) together: by default, *_join() will do a natural join, using all variables with common names across the two tables (works fine without the right order).
# left_join(): return all rows from x, and all columns from x and y. Rows in x with no match in y will have NA values in the new columns. If there are multiple matches between x and y, all combinations of the matches are returned.
for (i in 1:length(peakset_list_av)){ 
  peakset_list_comb[[as.numeric(i)]] = left_join(peakset_list_av[[as.numeric(i)]], peakset_list_sd[[as.numeric(i)]])
}	
for (i in 1:length(peakset_list_comb)){ 
  names(peakset_list_comb)[as.numeric(i)] = paste0(peakset_list_comb[[as.numeric(i)]]$Species[1], "_", 
                                                   peakset_list_comb[[as.numeric(i)]]$variable[1], "_",
                                                   peakset_list_comb[[as.numeric(i)]]$Peakset_type[1])
}

#finally plot them all out to a pdf (can use cowplot grid_arrange and other things for this, but this is simple
pdf("allPlots_JKM_181119.pdf")
library(ggplot2)
for (i in names(peakset_list_comb)) {   
  p = ggplot(peakset_list_comb[[i]], aes(x=Cell_subclass, y = average)) + geom_col() +
    geom_errorbar(aes(ymin=average-SD, ymax=average+SD), width = 0.3) + labs(title = i, y = "Percent overlap") + ylim(0,80) +
    theme(axis.text.x = element_text(angle = 90, size=rel(0.6))); # specify the font size relative to the base_size included in themes such as theme_bw() using the rel() function.
  print(p)	
}	
dev.off(4)

# An aethetic is a visual property of the objects in your plot. Aesthetics include things like the size, the shape, or the color of your points. 
# For each aesthetic, you use aes() to associate (map) the name of the aesthetic with a variable to display.
# One common problem when creating ggplot2 graphics is to put the + in the wrong place: it has to come at the end of the line, not the start.
