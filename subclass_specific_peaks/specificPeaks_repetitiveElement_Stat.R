rm(list=ls(all=TRUE)) 
setwd("/Users/shding/Dropbox/Overlap_RepeatElements_scripts/subclass_specific_peaks")

######################################################################################################################################
## import data without bootstrap.
######################################################################################################################################
x <-  read.table("specific_overlap_summary.txt", sep = " ", header = F)
# calculate the overlap percentage using column 1 and column 3.
percentages_x <-  x[,1] / x[,3] * 100
x <-  cbind(x, percentages_x)
#gsub("_peaks.*$","",x[,2])  # an example of using gsub and regular expression.
## The following part is to extract row names from column 1 of specific_overlap_summary.txt
list_strings <-  stringr::str_split(x[,2],"_")
row_names <- character(length(list_strings))
for (index in seq_along(list_strings)){
  if (grepl("olp",list_strings[[index]][5]) ){
    row_names[index] <- paste0(list_strings[[index]][1], "_All")
  } else {
    row_names[index] <- paste0(list_strings[[index]][1], "_", substring(list_strings[[index]][5],3))
  }
}
row_names_matrix <- row_names[1:68]
#unique(row_names)
#length(unique(row_names))
# convert the overlap percentage into a 68X7 matrix with row names and column names.
overlap_matrix = matrix(percentages_x, nrow = 68, ncol =7)
unique(sub(".*olp","",x[,2]))  # to check the order of column names.
overlap_matrix_colnames = c("Repetitive.elements", "LINE", "LTR", "Satellite", "Simple.repeat", "SINE", "DNA.transposon")
colnames(overlap_matrix) = overlap_matrix_colnames
rownames(overlap_matrix) = row_names_matrix

#heatmap(overlap_matrix)  # This function comes from base R and it is not very useful. Try ggplot2 instead.
## In order to use ggplot2 for heatmap, we have to melt the dataframe (matrix) to a univariate format.
# convert the rownames of the overlap_matrix to a proper column of the data frame.
overlap_df <- data.frame(cbind(Peak.type = rownames(overlap_matrix), overlap_matrix))
#If you want to then remove the original rownames:
rownames(overlap_df) <- NULL
df_heatmap <- melt(overlap_df,id.vars = "Peak.type")
names(df_heatmap)[2:3] <- c("RepeatElement_groups","overlap_percentage")
###This step is very important! Otherwise, Error: Discrete value supplied to continuous scale.
df_heatmap$overlap_percentage <- as.numeric(df_heatmap$overlap_percentage)

ggplot(df_heatmap,aes(x = RepeatElement_groups, y = Peak.type)) +
     geom_tile(aes(fill = overlap_percentage), color = "white") +
     scale_fill_gradient(low = "white", high = "steelblue") +
     ylab("Subclass-specific Peak Types ") +
     xlab("Repeat Elements Groups") +
     theme(legend.title = element_text(size = rel(0.3)),
                     legend.text = element_text(size = rel(0.3)),
                     plot.title = element_text(size=rel(1)),
                     axis.title=element_text(size= rel(0.5),face="bold"),
                     axis.text.y = element_text(size = rel(0.28)),
                     axis.text.x = element_text(size = rel(0.4),angle = 90, hjust = 1)) +
     labs(fill = "Overlap Percentage")


###########################################################################################################################################
## Import and Plot the bootstrap results.
###########################################################################################################################################
rm(list=ls())
setwd("/Users/shding/Dropbox/Overlap_RepeatElements_scripts/subclass_specific_peaks")
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
  Peakset_type_temp <- c(list_of_strings[[i]][2])
  Species[i] <- Species_temp
  Cell_subclass[i] <- Cell_subclass_temp
  Peakset_type[i] <- Peakset_type_temp
}
new_allInOne <- cbind(Name = new_allInOne$Name,Species, Cell_subclass, Peakset_type, new_allInOne[,c(2:ncol(new_allInOne))])
### Until this step, the right format of specific_STD_MEAN_includeShuffle_matrix is good for use in R.

library(reshape)
library(tibble)
library(dplyr)
library(stringr)
library(ggplot2)

av_cols = c(seq(5,31,2))
sd_cols = c(seq(6,32,2))
#split into two for avs and sds
avs = new_allInOne[,c(1:4,av_cols)]
sds = new_allInOne[,c(1:4,sd_cols)]
# calculate enrichment (%overlap for real data/ %overlap for shuffed data) and add to avs matrix.
# generate an empty matrix to hold data in the for loop.
enrichment_matrix = matrix(nrow= nrow(avs), ncol = 7)
# generate an empty vector with characters to save column names in the for loop.
colnames_enrichment = vector(mode="character", length=7) 
for (i in 5:11){
  enrichment_matrix[,i-4] = avs[,i]/avs[,i+7]
  colnames_enrichment[i-4] = paste0("enrichment_",gsub("intsxns100xCounts_.*$","", colnames(avs)[i]))
}
colnames(enrichment_matrix) = colnames_enrichment
avs = cbind(avs,enrichment_matrix)

### Keep the 68X7 matrix dimension before removing NA values (For the purpose to use the same code with non specific peaks).
## create new dataset without missing data
avs <- filter(avs, Name != "mL56IT_Con" & Name != "hL56IT_Con")
sds <- filter(sds, Name != "mL56IT_Con" & Name != "hL56IT_Con")
#avs <- na.omit(avs)
#sds <- na.omit(sds)

#format the data using melt function.
# melt function comes from the Reshape Package. Change multivariate format into univariate format (X1,X2,X3).
avs_melt = melt(data.frame(avs))
sds_melt = melt(data.frame(sds))
avs_melt$variable = gsub("intsxns100xCounts","",avs_melt$variable)   
# gsub() function repalces all matches of a string and returns a string vector of the same length. Unmatched substrings will returned unchanged.
sds_melt$variable = gsub("intsxns100xCounts","",sds_melt$variable)  
# gsub(pattern, replacement, x, ignore.case = FALSE). x is a string or string vector.

##################################################################################################################
# Plot: x axis(human_con, human_div, mouse_con, mouse_div). y axis (enrichment of overlapping with repeatMasker).
# Compute the average for the variable avs_melt$value grouped according to the Peakset_type, variable, and species.
### na.rm = TRUE is very important here!
avs_melt_aggregated = aggregate(avs_melt$value, 
                                by = list(avs_melt$Peakset_type, avs_melt$variable, avs_melt$Species),
                                FUN = mean, na.rm = TRUE)

names(avs_melt_aggregated) = c("Peakset_type","variable","Species","value")
# split the whole data.frame into a list of subset data.frames grouped according to variable types.
peakset_list_avs_melt_aggregated = split(avs_melt_aggregated, list(avs_melt_aggregated$variable))
# rename the Peakset_type to include information from Species (human or mouse) for better plotting outcomes.
peakset_list_avs_melt_aggregated[["enrichment_repeatMasker"]]$Peakset_type = 
  paste0(peakset_list_avs_melt_aggregated[["enrichment_repeatMasker"]]$Species,
         "_",peakset_list_avs_melt_aggregated[["enrichment_repeatMasker"]]$Peakset_type)
# use ggplot2 to generate a bar plot. geom_col() uses the heights of the bars to represent values in the data.
# geom_bar() makes the height of the bar proportional to the number of cases in each group. It counts the number of cases at each x position.
p_enrichRM_aggregate = ggplot(peakset_list_avs_melt_aggregated[["enrichment_repeatMasker"]], aes(x = Peakset_type, y = value)) + geom_col()+
  labs(title = "enrichment of peaks at repetitive elements", y = "enrichment at repetitive elements") +
  theme(axis.text.x = element_text(angle = 90, size = rel(1.3)))
print(p_enrichRM_aggregate)
##################################################################################################################

#### Research question: do any of the repetitive element types (LINEs, SINEs, LTRs, etc.) show positive enrichment (>1)?
repetitive_types = c("enrichment_repeatMasker","enrichment_LINE", 
                     "enrichment_LTR", "enrichment_Satellite", "enrichment_SimpleRepeat", 
                     "enrichment_SINE", "enrichment_DNATransposon")
# this is to initialize a list if you know how long your list is going to be. 
enrich_largerOne <- vector("list", 7)
# this is to initialize a vector to save the row numbers of each element in the list.
enrich_numbers = vector(mode="numeric", length = 7)
for (i in 1:7){
  enrich_largerOne[[i]] = avs[avs[repetitive_types[i]] > 1,
                              c("Name","Species","Cell_subclass", "Peakset_type", repetitive_types[i])];
  names(enrich_largerOne)[i] = repetitive_types[i];
  enrich_numbers[i] = nrow(enrich_largerOne[[i]]);
}
names(enrich_numbers) <- repetitive_types
# report the number of repetitive elements types with positive enrichment.
enrich_numbers
# report the detail individuals with positive enrichment.
enrich_largerOne
##################################################################################################################

### Research question: do any of the repetitive element types (LINEs, SINEs, etc.) show large differences between
#   conserved and divergent peak subsets?

enrich_matrix_new = melt(data.frame(cbind(avs[,1:4],avs[,repetitive_types])))
# add two new columns (Subclass_Species, Peaktype_variable) for better plotting performance.
enrich_matrix_new$Subclass_Species = paste0(enrich_matrix_new$Cell_subclass,"_",enrich_matrix_new$Species)
enrich_matrix_new$Peaktype_variable = paste0(enrich_matrix_new$variable, "_", enrich_matrix_new$Peakset_type)
# split the matrix into a list of matrix/dataframe based on column factors.
list_groupByPeakset_enrichment = split(enrich_matrix_new, list(enrich_matrix_new$Peakset_type, enrich_matrix_new$variable))
summary(list_groupByPeakset_enrichment)

pdf("SpecificPeaks_SubclassXRepeatElements.pdf")
for (i in names(list_groupByPeakset_enrichment)){
  p = ggplot(list_groupByPeakset_enrichment[[i]], aes(x = Subclass_Species, y = value, color = Species)) + geom_col(fill="darkgray") +
    labs(title = i, y = "Enrichment of overlap percent") + ylim(0,3) +
    theme(axis.text.x = element_text(angle = 90, size=rel(1)))
  print(p)
}   # here I include ylim(0,3) to give a constant scale across figures.
dev.off(4)


## plot conserved vs divergent peaks overlapped with 7 repeat element groups based on species (human and mouse).
list_groupByCellSubclass_enrichment = split(enrich_matrix_new, list(enrich_matrix_new$Cell_subclass, enrich_matrix_new$Species))
summary(list_groupByCellSubclass_enrichment)
pdf("SpecificPeaks_ConserveXDivergentXRepeatElementsXSpecies.pdf")
for (i in names(list_groupByCellSubclass_enrichment)){
  p2 = ggplot(list_groupByCellSubclass_enrichment[[i]], aes(x = Peaktype_variable, y = value, color = Peakset_type)) + geom_col(fill="darkgray") +
    labs(title = i, y = "Enrichment of overlap percent")  +
    theme(axis.text.x = element_text(angle = 90, size=rel(1)))
  print(p2)
}   # here I exclude ylim(0,3) to reflect the true values for each comparison in the same figure.
dev.off(4)

## plot peak enrichment merged across 11 subclasses based on species (human and mouse).
#names(list_groupByCellSubclass_enrichment)[1] ## To access the name of the first element in list_groupByCellSubclass_enrichment. 
enrich_matrix_agg1 = as.tibble(aggregate(enrich_matrix_new$value,
                                         by = list(enrich_matrix_new$Peaktype_variable, enrich_matrix_new$Species),
                                         FUN = mean, na.rm = TRUE))
names(enrich_matrix_agg1) = c("variable_peakType", "Species", "mean_acrossSubclass")
enrich_matrix_agg2 = as.tibble(aggregate(enrich_matrix_new$value,
                                         by = list(enrich_matrix_new$Peaktype_variable, enrich_matrix_new$Species),
                                         FUN = sd, na.rm = TRUE))
names(enrich_matrix_agg2) = c("variable_peakType", "Species", "std_acrossSubclass")
enrich_matrix_aggregate = left_join(enrich_matrix_agg1,enrich_matrix_agg2)
enrich_matrix_aggregate$variable_peakType_Species = paste0(enrich_matrix_aggregate$variable_peakType, "_", enrich_matrix_aggregate$Species)
pdf("Specific peak enrichment merged across cell subclasses.pdf")
p3 = ggplot(enrich_matrix_aggregate, aes(x = variable_peakType_Species, y = mean_acrossSubclass, color = Species)) + geom_col(fill="darkgray") +
  geom_errorbar(aes(ymin=mean_acrossSubclass-std_acrossSubclass, ymax= mean_acrossSubclass+std_acrossSubclass)) +
  labs(title = "peak enrichment merged across cell subclasses", y = "Enrichment of overlap percent")  +
  theme(axis.text.x = element_text(angle = 90, size=rel(1), hjust =1, vjust = 0))
print(p3)
dev.off(4)

