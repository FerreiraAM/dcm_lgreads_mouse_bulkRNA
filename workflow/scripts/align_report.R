#' Extract STAR mapping categories
#'
#' Extract percentage of reads falling into different mapping categories from the STAR log file
#' \code{Log.final.out}.
#'
#' @param star_log_file Path to the \code{Log.final.out} STAR log file.
get_mapping_cats <- function(star_log_file) {
  
  # read log file
  star_log <- readLines(star_log_file)
  
  # extract total number of reads
  total_reads <- star_log[6]
  
  # extract lines containing desired data
  values <- c("nb_input_reads" = star_log[6],
                  "average_input_read_length" = star_log[7],
                  "nb_uniquely_mapped" = star_log[9]
                  # "perc_uniquely_mapped" = star_log[10],
                  # "perc_multi_mapped" = star_log[25],
                  # "perc_too_many_loci" = star_log[27],
                  # "umapped_mismatch" = star_log[29],
                  # "unmapped_short" = star_log[30],
                  # "nb_unmapped_other" = star_log[31]
  )
  
  # only retain numbers from strings
  total_reads <- gsub(total_reads, pattern = "[^0-9]", replacement = "")
  values <- gsub(values, pattern = "[^0-9.]", replacement = "")
  
  # get mapping categories and convert to factor with specified levels
  categories <- factor(names(values), levels = names(values))
  
  # create data.frame with values in mapping cats
  data.frame(category = categories, values = as.numeric(values))
  
}



# Log files 
log_files <- unlist(snakemake@input)
# Extract information
list_info <- lapply(X = log_files, FUN = function(one_log){
  # Extract samples info from the path
  break_path <- unlist(strsplit(one_log, "/"))
  sample <- break_path[4]
  experiment <- break_path[3]
  # Extract info from log file
  df_info <- get_mapping_cats(one_log)
  data.frame(df_info, "sample" = sample, "experiment" = experiment)
})
# As data.frame
table_nb_reads <- do.call(rbind, list_info)
# Save
write.table(table_nb_reads, file = unlist(snakemake@output), col.names = TRUE, row.names = FALSE)
