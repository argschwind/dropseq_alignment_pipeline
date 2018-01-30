#' Extract STAR mapping categories
#' 
#' Extract percentage of reads falling into different mapping categories from
#' the STAR log file \code{Log.final.out}.
#' 
#' @param star_log_file Path to the \code{Log.final.out} STAR log file.
#' 
get_mapping_cats <- function(star_log_file){
 
  # read log file
  star_log <- readLines(star_log_file)
  
  # extract total number of reads
  total_reads <- star_log[6]
  
  # extract lines containing desired data
  perc_reads <- c("uniquely_mapped" = star_log[10],
                   "multi_mapped" = star_log[25],
                   "too_many_loci" = star_log[27],
                   "umapped_mismatch" = star_log[29],
                   "unmapped_short" = star_log[30],
                   "unmapped_other" = star_log[31]
                   )
  
  # only retain numbers from strings
  total_reads <- gsub(total_reads, pattern = "[^0-9]", replacement = "")
  perc_reads <- gsub(perc_reads, pattern = "[^0-9.]", replacement = "")
  
  # get mapping categories and convert to factor with specified levels
  categories <- factor(names(perc_reads), levels = names(perc_reads))
  
  # create data.frame with percentage of reads in mapping cats
  data.frame(category = categories,
             perc_reads = as.numeric(perc_reads)
             )
  
  }
  
#' Plot STAR mapping categories
#' 
#' Plot percentage of reads falling into different mapping categories as
#' bar chart.
#' 
#' @param mapping_cats data.frame containing the percentage of reads in mapping
#' categories produced by get_mapping_cats().
#' 
#' @param sample Sample name for main title.
#' 
#' @param total_reads (optional) Number of total reads. If provided the total
#' number of reads is added to the plot title.
#' 
plot_mapping_cats <- function(mapping_cats, sample, total_reads = NULL){
  
  # generate title
  if (!is.null(total_reads)){
    
    title <- paste0(sample, " mapping results (total reads: ", total_reads, ")")
    
  }else{
    
    title <- paste0(sample, " mapping results")
    
  }
  
  # create barplot showing percentage of reads in mapping categories
  ggplot2::ggplot(mapping_cats, ggplot2::aes(x = category, y = perc_reads)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(y = "percentage of reads", x = "STAR mapping category",
                  title = title) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       vjust = 1)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 20),
                                limits = c(0, 100))
  
}
