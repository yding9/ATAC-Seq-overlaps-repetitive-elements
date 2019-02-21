###########################################################################################################################################
## Import and Plot the bootstrap results.
###########################################################################################################################################
rm(list=ls())
setwd("C:/Users/yi.ding/Dropbox/Overlap_RepeatElements_scripts/subclass_specific_peaks")
library(reshape)
library(tibble)
library(dplyr)
library(stringr)
library(ggplot2)
## import bootstrap data in excel obtained by using calculate_mean_std_from_bootstrap_specificPeaks.R
## The following is trying to reformat the data to include columns like Species, Name, Cell_subclass, Peakset_type.
new_allInOne <- read.table("specific_STD_MEAN_includeShuffle_matrix.xls", sep = "\t", header = T)
new_allInOne <- cbind(Name = rownames(new_allInOne),new_allInOne)
rownames(new_allInOne) <- NULL
list_of_strings <- stringr::str_split(new_allInOne$Name, "_")
Species <- character(length = length(list_of_strings))
Cell_subclass <- character(length = length(list_of_strings))
Peakset_type <- character(length = length(list_of_strings))
for (i in seq_along(list_of_strings)){
  if ( grepl("^h",list_of_strings[[i]][1]) ){
    Species_temp <- c("Human")
  } else{
    Species_temp <- c("Mouse")
  }
  Cell_subclass_temp <- c(substring(list_of_strings[[i]][1],2))
  if (Cell_subclass_temp == "Astrocyte"){
    Cell_subclass_temp <- "Astro"
  }
  if (Cell_subclass_temp == "Microglia"){
    Cell_subclass_temp <- "Micro"
  }
  if (Cell_subclass_temp == "OligodendrocyteOPC"){
    Cell_subclass_temp <- "OligoOPC"
  }
  Peakset_type_temp <- c(list_of_strings[[i]][2])
  if (Peakset_type_temp == "Con"){
    Peakset_type_temp <- "Conserved"
  }
  if (Peakset_type_temp == "Div"){
    Peakset_type_temp <- "Divergent"
  }
  Species[i] <- Species_temp
  Cell_subclass[i] <- Cell_subclass_temp
  Peakset_type[i] <- Peakset_type_temp
}

new_allInOne <- cbind(Name = new_allInOne$Name, Species, Cell_subclass, Peakset_type, new_allInOne[,c(2:ncol(new_allInOne))])

av_cols = c(seq(5,31,2))

avs = new_allInOne[,c(1:4,av_cols)]

# calculate enrichment (%overlap for real data/ %overlap for shuffed data) and create a new matrix only contains enrichment.
# generate an empty matrix to hold data in the for loop.
enrichment_matrix = matrix(nrow= nrow(avs), ncol = 7)
# generate an empty vector with characters to save column names in the for loop.
colnames_enrichment = vector(mode="character", length=7) 
for (i in 5:11){
  enrichment_matrix[,i-4] = avs[,i]/avs[,i+7]
  colnames_enrichment[i-4] = paste0(gsub("intsxns100xCounts_.*$","", colnames(avs)[i]))
}
colnames(enrichment_matrix) = colnames_enrichment
enrichment_df <- cbind(avs[,c(1:4)], enrichment_matrix)

## create new dataset without missing data and only contains conserved and divergent peak types.
enrichment_df_ConDiv <- as.tibble(filter(enrichment_df, Name != "mL56IT_Con" & Name != "hL56IT_Con" & Peakset_type != "All"))
enrichment_df_ConDiv$Cell_subclass <- as.character(enrichment_df_ConDiv$Cell_subclass)

## integrate color codes into enrichment_df
## this is to generate a hash map manually. It's better to use the left_join function.
# library(hash)
# color_hash <- hash(c('Astro','DL','L23','L4','L56IT','LAMP5','Micro','OligoOPC','PVALB','SST','VIP'), 
#                    c('#604f22','#256100', '#9cea95','#00b4bb', '#000C98','#e5938d', '#AAAAAA','#577d76','#e82c49','#E16115','#b200b8'))
# subclass_color <- character(length = nrow(enrichment_df_ConDiv))
# for (i in seq(1:nrow(enrichment_df_ConDiv))){
#   subclass_color[i] <- color_hash[[as.character(enrichment_df_ConDiv[i,3])]]
# }
# enrichment_df_ConDiv_colorCode <- cbind(subclass_color = subclass_color, enrichment_df_ConDiv)                            
#### left_join: must remove duplicate in the binding data frame!
load("C:/Users/yi.ding/Dropbox/Overlap_RepeatElements_scripts/subclass_specific_peaks/real_minus_random_means_cons_div_181211.rda")
color_code <- select(as.tibble(real_minus_random_means_cons_div), subclass_label, subclass_color)
color_code <-  distinct(color_code,subclass_label, .keep_all= TRUE)
enrichment_df_ConDiv_colorCode <- left_join(enrichment_df_ConDiv, color_code, by = c("Cell_subclass" = "subclass_label"))


#format the data using melt function.
enrichment_melt <- melt(data.frame(enrichment_df_ConDiv_colorCode))
colnames(enrichment_melt)[c(6,7)] <- c("Repetitive_element_type", "Enrichment_value")
plot_class <- paste0(as.character(enrichment_melt$Repetitive_element_type), "_", as.character(enrichment_melt$Peakset_type))
enrichment_melt <- cbind(enrichment_melt, plot_class = plot_class)
save(enrichment_melt, file = "enrichment_melt.rda")

##jitter plots using ggplot2 and pre-determined color codes.
pdf("Human_OverlapRE_Enrichment_Cons_Div.pdf", useDingbats=F)
ggplot(filter(enrichment_melt, Species == "Human"), 
            aes(x = fct_relevel(plot_class, "repeatMasker_Conserved", "repeatMasker_Divergent"), y = Enrichment_value)) +
  geom_jitter(aes( col = subclass_color), position=position_jitter(0.2)) + scale_color_identity() +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                geom = "crossbar", width = 0.5) + 
  labs(title="Human Peaks",x="Repetitive element types", y = "Enrichment") +
  theme(axis.text.x = element_text(angle = 90, size=rel(1), hjust =1, vjust = 0))
	

pdf("Mouse_OverlapRE_Enrichment_Cons_Div.pdf", useDingbats=F)
ggplot(filter(enrichment_melt, Species == "Mouse"), 
       aes(x = fct_relevel(plot_class, "repeatMasker_Conserved", "repeatMasker_Divergent"), y = Enrichment_value)) +
  geom_jitter(aes( col = subclass_color), position=position_jitter(0.2)) + scale_color_identity() +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5) + 
  labs(title="Mouse Peaks",x="Repetitive element types", y = "Enrichment") +
  theme(axis.text.x = element_text(angle = 90, size=rel(1), hjust =1, vjust = 0))
dev.off()	















## one question here: why the order in figure p1 doesn't match the data order in enrichment_melt data.frame.
(pp1 <- ggplot(filter(enrichment_melt, Species == "Human"), 
              aes(x = Repetitive_element_type, y = Enrichment_value)) +
    geom_jitter(aes( col = subclass_color), position=position_jitter(0.2)) + scale_color_identity() +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.5) + 
    labs(title="Human Peaks",x="Repetitive element types", y = "Enrichment") +
    theme(axis.text.x = element_text(angle = 90, size=rel(1), hjust =1, vjust = 0)))












