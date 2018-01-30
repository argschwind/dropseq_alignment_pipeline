#' Cumulative fraction of reads per cell barcode
#'
#' Calculate cumulative fraction of reads per detected cell barcode. The
#' cumulative fraction can be plotted to estimate the number of sequenced cells
#' by identifying the "knee" of the distribution. For further details on
#' selection of cells, see the Drop-seq alignment cookbook provided by the
#' \href{http://mccarrolllab.com/dropseq/}{McCarroll lab}.
#'
#' @param read_counts_file Path to output file from BAMTagHistogram containing
#' the number of reads per cell barcode. This file should contain two columns,
#' where the number of reads must be in the first and the barcode identity in
#' the second column.
#'
#' @return A data.frame containing the original data in read_counts with an
#' added columns containing the cumulative sum and fraction.
#' 
cumfrac_rpc <- function(read_counts_file){

  # read data
  read_counts <- utils::read.table(read_counts_file, stringsAsFactors = FALSE)

  # set colnames for read_counts
  colnames(read_counts)[1:2] <- c("reads", "cell_barcode")

  # make sure that read_counts is ordered correctly
  read_counts <- read_counts[order(read_counts[, 1], decreasing = TRUE), ]

  # calculate cumulative sum
  reads_cumsum <- cumsum(read_counts[, 1])

  # calculate cumulative fraction
  reads_cumfrac <- reads_cumsum / max(reads_cumsum)

  # create output data.frame
  cbind(read_counts, reads_cumsum, reads_cumfrac)

}

#' Plot cumulative fraction of reads per cell barcode
#'
#' Plot cumulative fraction of reads per detected cell barcode. The cumulative
#' fraction can be plotted to estimate the number of sequenced cells by
#' identifying the "knee" of the distribution. For further details on selection
#' of cells, see the Drop-seq alignment cookbook provided by the
#' \href{http://mccarrolllab.com/dropseq/}{McCarroll lab}.
#'
#' @param cumfrac A data.frame containing the cumulative fraction of reads per
#' cell barcode calculated by cumfrac_rpc()
#'
#' @param nbcs Number of cell barcodes to be plotted on the x-axis.
#'
#' @param title String containing plot title.
#' 
plot_cumfrac_rpc <- function(cumfrac, nbcs, title = NULL){

  # only retain nbcs number of barcodes
  cumfrac <- cumfrac[1:nbcs, ]

  # add index used for plotting
  cumfrac$index <- 1:nrow(cumfrac)

  # define axis
  xaxis <- list(title = "Cell barcodes sorted by number of reads [descending]")
  yaxis <- list(title = "Cumulative fraction of reads",
                 range = c(0, 1))
  
  # create interactive line plot using plotly
  p <- plotly::plot_ly(cumfrac, x = ~index, y = ~reads_cumfrac,
                       type = "scatter", mode = "lines",
                       text = ~paste("Cumulative reads: ", reads_cumsum)
                       )
  
  # format title and axis
  plotly::layout(p, title = title, xaxis = xaxis, yaxis = yaxis)

}
