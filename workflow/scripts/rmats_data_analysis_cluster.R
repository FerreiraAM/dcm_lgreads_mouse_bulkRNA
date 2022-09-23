# rMATS data analysis
# save.image("rmats.rda")
# stop()

# Libraries
library("tidyverse")

# Functions
# Read files and add splice junction type
fun_read_SJ <- function(ls_file_names){
  lapply(ls_file_names, FUN = function(onefile){
    # Read file
    df <- read.table(onefile, header = TRUE)
    # Add splice junction type
    df$SJ_name <- unlist(strsplit(x = basename(onefile), "[.]"))[1]
    # Return
    df
  })
}
# Combine list of file in a data.frame
fun_comb_in_df <- function(ls_files, ls_file_names){
  # Rename
  names(ls_files) <- basename(ls_file_names)
  # Dataframe
  dplyr::bind_rows(ls_files)
}
# Basics stats from the files
fun_basic_stats <- function(df_files, fdr = 0.01){
  # Stats
  nb_event_per_SJtype <- df_files %>% 
    group_by(SJ_name) %>% 
    count(SJ_name)
  # Number of affected genes 
  nb_affected_genes <- length(unique(df_files$geneSymbol))
  # Number of gene affected by a significant SJ event
  df_files_sig <- dplyr::filter(df_files, FDR < fdr)
  nb_genes_sig <- length(unique(df_files_sig$geneSymbol))
  # Return
  list("nb_event_per_SJtype" = nb_event_per_SJtype,
       "nb_affected_genes" = nb_affected_genes,
       "nb_genes_sig" = nb_genes_sig)
}

# R636Q HOM vs WT
# List files
ls_R636Q_HOMvsWT_JC_filesname <- list.files(path = "rmats_directory/R636Q_HOMvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_R636Q_HOMvsWT_JC_files <- fun_read_SJ(ls_R636Q_HOMvsWT_JC_filesname)
# Convert to a data.frame
df_R636Q_HOMvsWT_JC_files <- fun_comb_in_df(ls_R636Q_HOMvsWT_JC_files, 
                                            ls_R636Q_HOMvsWT_JC_filesname)
# R636Q HET vs WT
# List files
ls_R636Q_HETvsWT_JC_filesname <- list.files(path = "rmats_directory/R636Q_HETvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_R636Q_HETvsWT_JC_files <- fun_read_SJ(ls_R636Q_HETvsWT_JC_filesname)
# Convert to a data.frame
df_R636Q_HETvsWT_JC_files <- fun_comb_in_df(ls_R636Q_HETvsWT_JC_files, 
                                            ls_R636Q_HETvsWT_JC_filesname)

# P635L
# HOM vs WT
# List files
ls_P635L_HOMvsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_HOMvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_P635L_HOMvsWT_JC_files <- fun_read_SJ(ls_P635L_HOMvsWT_JC_filesname)
# Convert to a data.frame
df_P635L_HOMvsWT_JC_files <- fun_comb_in_df(ls_P635L_HOMvsWT_JC_files, 
                                            ls_P635L_HOMvsWT_JC_filesname)
# HET vs WT
# List files
ls_P635L_HETvsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_HETvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_P635L_HETvsWT_JC_files <- fun_read_SJ(ls_P635L_HETvsWT_JC_filesname)
# Convert to a data.frame
df_P635L_HETvsWT_JC_files <- fun_comb_in_df(ls_P635L_HETvsWT_JC_files, 
                                            ls_P635L_HETvsWT_JC_filesname)

###### Overlap ########
# HOM vs WT
# Add a column in each dataframe
df_P635L_HOMvsWT_JC_files$mutant <- "P635L"
df_P635L_HOMvsWT_JC_files$comparison <- "HOMvsWT"
df_R636Q_HOMvsWT_JC_files$mutant <- "R636Q"
df_R636Q_HOMvsWT_JC_files$comparison <- "HOMvsWT"
df_P635L_HETvsWT_JC_files$mutant <- "P635L"
df_P635L_HETvsWT_JC_files$comparison <- "HETvsWT"
df_R636Q_HETvsWT_JC_files$mutant <- "R636Q"
df_R636Q_HETvsWT_JC_files$comparison <- "HETvsWT"
# Intersect all
intersect_all_genesymbol <- Reduce(intersect, list(unique(df_P635L_HOMvsWT_JC_files$geneSymbol),
                                                   unique(df_R636Q_HOMvsWT_JC_files$geneSymbol),
                                                   unique(df_P635L_HETvsWT_JC_files$geneSymbol),
                                                   unique(df_R636Q_HETvsWT_JC_files$geneSymbol)))
# Combine
df_all_JC_files <- rbind(df_P635L_HOMvsWT_JC_files,
                         df_R636Q_HOMvsWT_JC_files,
                         df_P635L_HETvsWT_JC_files,
                         df_R636Q_HETvsWT_JC_files)
# Filtering
# Common gene symbols
df_all_JC_files_common_genesymbols <- dplyr::filter(df_all_JC_files, 
                                                    geneSymbol %in% intersect_all_genesymbol)
# For each gene
ls_sameevent_pergene <- lapply(unique(df_all_JC_files_common_genesymbols$geneSymbol), 
                               FUN = function(onegene, df = df_all_JC_files_common_genesymbols){
                                 print(onegene)
                                 # Filter the dataframe
                                 df_of_interest <- dplyr::filter(df, geneSymbol == onegene)
                                 # SJ names
                                 SJ_name_of_interest <- unique(df_of_interest$SJ_name)
                                 # Check for each time of SJ
                                 # SE
                                 # SE: exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
                                 # The inclusion form includes the target exon (exonStart_0base, exonEnd)
                                 if("SE" %in% SJ_name_of_interest){
                                   # Filter
                                   SE_test <- dplyr::filter(df_of_interest, SJ_name == "SE")
                                   # Count the number of unique combination of coordinates
                                   # nb_find = 4 means that it was found 4 times (in the 4 conditions)
                                   SE_nb_SJ <- SE_test %>% 
                                     group_by(exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE) %>% 
                                     summarize(nb_find = n())
                                   # Filter nb_find == 4
                                   SE_nb_SJ_4 <- dplyr::filter(SE_nb_SJ, nb_find == 4)
                                   # If none
                                   if(dim(SE_nb_SJ_4)[1] != 0){
                                     # Filter these SJ events
                                     SE_ls_SJ_4 <- apply(as.data.frame(SE_nb_SJ_4), 1, FUN = function(onecoord){
                                       dplyr::filter(SE_test,
                                                     exonStart_0base == as.numeric(onecoord[c("exonStart_0base")]) &
                                                       exonEnd == as.numeric(onecoord[c("exonEnd")]) &
                                                       upstreamES == as.numeric(onecoord[c("upstreamES")]) &
                                                       upstreamEE == as.numeric(onecoord[c("upstreamEE")]) &
                                                       downstreamES == as.numeric(onecoord[c("downstreamES")]) & 
                                                       downstreamEE == as.numeric(onecoord[c("downstreamEE")])) 
                                     })
                                     # As data.frame
                                     SE_df_SJ_4 <- bind_rows(SE_ls_SJ_4, .id = "SJ_ID")
                                   } else { 
                                     # As data.frame
                                     SE_df_SJ_4 <- data.frame()
                                   }
                                 } else { 
                                   # As data.frame
                                   SE_df_SJ_4 <- data.frame()
                                 }
                                 # RI 
                                 # RI: riExonStart_0base riExonEnd upstreamES upstreamEE downstreamES downstreamEE
                                 # The inclusion form includes (retains) the intron (upstreamEE, downstreamES)
                                 if("RI" %in% SJ_name_of_interest){
                                   # Filter
                                   RI_test <- dplyr::filter(df_of_interest, SJ_name == "RI")
                                   # Count the number of unique combination of coordinates
                                   # nb_find = 4 means that it was found 4 times (in the 4 conditions)
                                   RI_nb_SJ <- RI_test %>% 
                                     group_by(riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE) %>% 
                                     summarize(nb_find = n())
                                   # Filter nb_find == 4
                                   RI_nb_SJ_4 <- dplyr::filter(RI_nb_SJ, nb_find == 4)
                                   # If none
                                   if(dim(RI_nb_SJ_4)[1] != 0){
                                     # Filter these SJ events
                                     RI_ls_SJ_4 <- apply(as.data.frame(RI_nb_SJ_4), 1, FUN = function(onecoord){
                                       dplyr::filter(RI_test,
                                                     riExonStart_0base == as.numeric(onecoord[c("riExonStart_0base")]) &
                                                       riExonEnd == as.numeric(onecoord[c("riExonEnd")]) &
                                                       upstreamES == as.numeric(onecoord[c("upstreamES")]) &
                                                       upstreamEE == as.numeric(onecoord[c("upstreamEE")]) &
                                                       downstreamES == as.numeric(onecoord[c("downstreamES")]) & 
                                                       downstreamEE == as.numeric(onecoord[c("downstreamEE")])) 
                                     })
                                     # As data.frame
                                     RI_df_SJ_4 <- bind_rows(RI_ls_SJ_4, .id = "SJ_ID")
                                   } else { 
                                     # As data.frame
                                     RI_df_SJ_4 <- data.frame()
                                   }
                                 } else { 
                                   # As data.frame
                                   RI_df_SJ_4 <- data.frame()
                                 }
                                 # MXE
                                 # MXE: 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
                                 # If the strand is + then the inclusion form includes the 1st exon (1stExonStart_0base, 1stExonEnd) and skips the 2nd exon
                                 # If the strand is - then the inclusion form includes the 2nd exon (2ndExonStart_0base, 2ndExonEnd) and skips the 1st exon
                                 if("MXE" %in% SJ_name_of_interest){
                                   # Filter
                                   MXE_test <- dplyr::filter(df_of_interest, SJ_name == "MXE")
                                   # Count the number of unique combination of coordinates
                                   # nb_find = 4 means that it was found 4 times (in the 4 conditions)
                                   MXE_nb_SJ <- MXE_test %>% 
                                     group_by(X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE) %>% 
                                     summarize(nb_find = n())
                                   # Filter nb_find == 4
                                   MXE_nb_SJ_4 <- dplyr::filter(MXE_nb_SJ, nb_find == 4)
                                   # If none
                                   if(dim(MXE_nb_SJ_4)[1] != 0){
                                     # Filter these SJ events
                                     MXE_ls_SJ_4 <- apply(as.data.frame(MXE_nb_SJ_4), 1, FUN = function(onecoord){
                                       dplyr::filter(MXE_test,
                                                     X1stExonStart_0base == as.numeric(onecoord[c("X1stExonStart_0base")]) &
                                                       X1stExonEnd == as.numeric(onecoord[c("X1stExonEnd")]) &
                                                       X2ndExonStart_0base == as.numeric(onecoord[c("X2ndExonStart_0base")]) &
                                                       X2ndExonEnd == as.numeric(onecoord[c("X2ndExonEnd")]) &
                                                       upstreamES == as.numeric(onecoord[c("upstreamES")]) & 
                                                       upstreamEE == as.numeric(onecoord[c("upstreamEE")]) &
                                                       downstreamES == as.numeric(onecoord[c("downstreamES")]) & 
                                                       downstreamEE == as.numeric(onecoord[c("downstreamEE")])) 
                                     })
                                     # As data.frame
                                     MXE_df_SJ_4 <- bind_rows(MXE_ls_SJ_4, .id = "SJ_ID")
                                   } else { 
                                     # As data.frame
                                     MXE_df_SJ_4 <- data.frame()
                                   }
                                 } else { 
                                   # As data.frame
                                   MXE_df_SJ_4 <- data.frame()
                                 }
                                 # A3SS or A5SS
                                 # A3SS, A5SS: longExonStart_0base longExonEnd shortES shortEE flankingES flankingEE
                                 # The inclusion form includes the long exon (longExonStart_0base, longExonEnd) instead of the short exon (shortES shortEE)
                                 if("A3SS" %in% SJ_name_of_interest | "A5SS" %in% SJ_name_of_interest){
                                   # Filter
                                   ASS_test <- dplyr::filter(df_of_interest, SJ_name == "A3SS" | SJ_name == "A5SS")
                                   # Count the number of unique combination of coordinates
                                   # nb_find = 4 means that it was found 4 times (in the 4 conditions)
                                   ASS_nb_SJ <- ASS_test %>% 
                                     group_by(longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, SJ_name) %>% 
                                     summarize(nb_find = n())
                                   # Filter nb_find == 4
                                   ASS_nb_SJ_4 <- dplyr::filter(ASS_nb_SJ, nb_find == 4)
                                   # If none
                                   if(dim(ASS_nb_SJ_4)[1] != 0){
                                     # Filter these SJ events
                                     ASS_ls_SJ_4 <- apply(as.data.frame(ASS_nb_SJ_4), 1, FUN = function(onecoord){
                                       dplyr::filter(ASS_test,
                                                     longExonStart_0base == as.numeric(onecoord[c("longExonStart_0base")]) &
                                                       longExonEnd == as.numeric(onecoord[c("longExonEnd")]) &
                                                       shortES == as.numeric(onecoord[c("shortES")]) &
                                                       shortEE == as.numeric(onecoord[c("shortEE")]) &
                                                       flankingES == as.numeric(onecoord[c("flankingES")]) & 
                                                       flankingEE == as.numeric(onecoord[c("flankingEE")]) & 
                                                       SJ_name == (onecoord[c("SJ_name")])) 
                                     })
                                     # As data.frame
                                     ASS_df_SJ_4 <- bind_rows(ASS_ls_SJ_4, .id = "SJ_ID")
                                   } else { 
                                     # As data.frame
                                     ASS_df_SJ_4 <- data.frame()
                                   }
                                 } else { 
                                   # As data.frame
                                   ASS_df_SJ_4 <- data.frame()
                                 }
                                 # Combine
                                 bind_rows(SE_df_SJ_4,
                                           RI_df_SJ_4,
                                           MXE_df_SJ_4,
                                           ASS_df_SJ_4)
                               })
saveRDS(object = ls_sameevent_pergene, file = "rmats_directory/ls_sameevent_pergene.rds")

