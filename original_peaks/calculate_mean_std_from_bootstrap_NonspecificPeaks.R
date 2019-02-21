# the original std_mean_all.R script has bugs. The output is not in the right order. 
#So do not use external rowname files.Extract the name information from input files directly.
### Change working directory as needed.If run on linux machine, and files and scripts are under the same directory, no need to setwd.
#setwd("C:/Users/yi.ding/Desktop/specific_peaks")
rm(list = ls(all=T))

# Calculate the standard devitation_mean matrix using the real data with bootstrap using 80% peaks.
file_loop = c("repeatMaskerintsxns100xCounts", "LINEintsxns100xCounts", "LTRintsxns100xCounts", "Satelliteintsxns100xCounts",
              "SimpleRepeatintsxns100xCounts", "SINEintsxns100xCounts", "DNATransposonintsxns100xCounts", 
              "shuffled_repeatMaskerintsxns100xCounts", "shuffled_LINEintsxns100xCounts", "shuffled_LTRintsxns100xCounts", 
              "shuffled_Satelliteintsxns100xCounts","shuffled_SimpleRepeatintsxns100xCounts", 
              "shuffled_SINEintsxns100xCounts", "shuffled_DNATransposonintsxns100xCounts")

num_peak_type = 68;
row_names <- character(num_peak_type)
col_names <- character(2 * length(file_loop))
result <- matrix(,nrow = num_peak_type, ncol= 2 * length(file_loop))
for (index in seq_along(file_loop)){
  if (index <= 7) {
    inFileNames <- sort(grep(list.files(path = ".", pattern = file_loop[index], full.names = T), pattern = 'shuffled_', inv = T, value = T));
  } else {
    inFileNames <- sort(list.files(path = ".", pattern = file_loop[index], full.names = T));
  }
  list_strings <- stringr::str_split(inFileNames,"_")
  if (length(inFileNames) != num_peak_type) {
    print("Incorrect number of files!!!");
  }
  for (i in seq_along(inFileNames)){
    # the difference btw nonspecific peaks and specific peaks is the index of string elements.
    if (index <= 7) {
      if (grepl(file_loop[index], list_strings[[i]][4])){
        row_name_temp <- paste0(substring(list_strings[[i]][1],3),"_All")
      } else {
        row_name_temp <- paste0(substring(list_strings[[i]][1],3),"_", substring(list_strings[[i]][4],3))
      }
    } else {
      if (grepl("shuffled", list_strings[[i]][4])){
        row_name_temp <- paste0(substring(list_strings[[i]][1],3),"_All")
      } else {
        row_name_temp <- paste0(substring(list_strings[[i]][1],3),"_", substring(list_strings[[i]][4],3))
      }
    }
    
    row_ind = 0
    if (index == 1) {
      row_names[i] = row_name_temp;
      row_ind = i;
    } else {
      row_ind = match(row_name_temp, row_names);
    }
    bootstrap_percent <-  as.numeric(readLines(inFileNames[i]))*100;
    std <-  sd(bootstrap_percent);
    average <-  mean(bootstrap_percent);
    std_mean_vector <-  c(average, std);
    result[row_ind, c(2 * index - 1,2 * index)] = std_mean_vector
  }
  col_names[c(2 * index - 1,2 * index)] = c(paste0(file_loop[index],"_average"), paste0(file_loop[index], "_STD"))
}
colnames(result) <- col_names
rownames(result) <- row_names
write.table(result, file = "Nonspecific_STD_MEAN_includeShuffle_matrix.xls", sep = "\t")

###################################################################################################################################################
### this block of code works, but the order of peak types are not well controlled. Have to use tibble later to leftjoin by Peak.type. 
# file_loop = c("repeatMaskerintsxns100xCounts", "LINEintsxns100xCounts", "LTRintsxns100xCounts", "Satelliteintsxns100xCounts",
#               "SimpleRepeatintsxns100xCounts", "SINEintsxns100xCounts", "DNATransposonintsxns100xCounts")
# 
# list_of_df <- vector("list", length(file_loop))
# for (index in seq_along(file_loop)){
#   inFileNames <- sort(grep(list.files(path = ".", pattern = file_loop[index], full.names = T), pattern = 'shuffled_', inv = T, value = T));
#   list_strings <- stringr::str_split(inFileNames,"_")
#   row_names <- character(length(inFileNames))
#   matrix_1 <- matrix(,nrow = length(inFileNames), ncol = 2)
#   for (i in seq_along(inFileNames)){
#     if (grepl(file_loop[index], list_strings[[i]][5])){
#       row_names[i] <- paste0(substring(list_strings[[i]][1],3),"_All")
#     } else {
#       row_names[i] <- paste0(substring(list_strings[[i]][1],3),"_", substring(list_strings[[i]][5],3))
#     }
#     bootstrap_percent <-  as.numeric(readLines(inFileNames[i]))*100;
#     std <-  sd(bootstrap_percent);
#     average <-  mean(bootstrap_percent);
#     std_mean_vector <-  c(average, std);
#     matrix_1[i,] <- std_mean_vector;
#     
#   }
#   colnames(matrix_1) <- c(paste0(file_loop[index],"_average"), paste0(file_loop[index], "_SD"));
#   df_mean_std <- data.frame(cbind(Peak.type = row_names, matrix_1))
#   list_of_df[[index]] <- df_mean_std
# }
###################################################################################################################################################

