#' Adapter trimming histogram
#' 
#' Plot the number of bases clipped from Drop-seq genome reads by
#' TrimStartingSequence from Drop-seq tools.
#' 
#' @param adapter_trim_file Path to the TrimStartingSequence summary file.
#' 
#' @param sample Sample name for main title.
#' 
#' @param total_reads (optional) Number of total reads. If provided, the
#' percentage of reads that were trimmed will be included in the plot title.
#' 
adapter_trim_hist <- function(adapter_trim_file, sample, total_reads = NULL){
  
  # read data
  adapter_trim <- read.table(adapter_trim_file, header = TRUE, skip = 5)
  
  # generate title
  if (!is.null(total_reads)){
    
    # calculate precentage of trimmed reads
    perc_trim <- round(sum(adapter_trim$VALUE) / total_reads * 100, digits = 3)
    
    # make title
    title <- paste0(sample, " adapter trimmed reads (", perc_trim, "%)")
    
  }else{
    
    title <- paste(sample, "adapter trimmed reads")
    
  }
  
  # plot histogram
  ggplot2::ggplot(adapter_trim, ggplot2::aes(BIN, VALUE)) + 
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "Number of clipped bases", y = "Reads", title = title)
  
}

#' PolyA trimming histogram
#' 
#' Plot post trimming read lengths of Drop-seq reads that were trimmed by
#' PolyATrimmer from Drop-seq tools.
#' 
#' @param polyA_trim_file Path to the PolyATrimmer summary file.
#' 
#' @param sample Sample name for main title.
#' 
#' @param total_reads (optional) Number of total reads. If provided, the
#' percentage of reads that were trimmed will be included in the plot title.
#' 
polyA_trim_hist <- function(polyA_trim_file, sample, total_reads = NULL){
  
  # read data
  polyA_trim <- read.table(polyA_trim_file, header = TRUE, skip = 5)
  
  # generate title
  if (!is.null(total_reads)){
    
    # calculate precentage of trimmed reads
    perc_trim <- round(sum(polyA_trim$VALUE) / total_reads * 100, digits = 3)
    
    # make title
    title <- paste0(sample, " polyA trimmed reads (", perc_trim, "%)")
    
  }else{
    
    # make title
    title <- paste(sample, "polyA trimmed reads")
    
  }
  
  # plot histogram
  ggplot2::ggplot(polyA_trim, ggplot2::aes(BIN, VALUE)) + 
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "Read lengths after trimming (bp)", y = "Reads",
                  title = title)
  
}
