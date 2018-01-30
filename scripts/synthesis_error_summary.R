#' Read bead synthesis error summary
#' 
#' Read summary data on bead synthesis errors produced by
#' DetectBeadSynthesisErrors from Drop-seq tools.
#' 
#' @param synthesis_error_file Path to DetectBeadSynthesisErrors summary file.
#' 
#' @return A list with two data.frames containing bead synthesis metrics
#' and histogram data.
read_synthesis_error <- function(synthesis_error_file){
  
  # read summary data
  bc_error_dat <- readLines(synthesis_error_file)
  
  # remove first and last empty line
  bc_error_dat <- bc_error_dat[-c(1,length(bc_error_dat))]
  
  # remove comment lines
  bc_error_dat <- grep(bc_error_dat, pattern = "^#+\\s.+", perl = TRUE,
                       value = TRUE, invert = TRUE)
  
  # find empty line separating metrics from histogram tables
  sep_line <- which(bc_error_dat == "")
  
  # split strings into vectors
  bc_error_dat <- strsplit(bc_error_dat, split = "\t")
  
  # transform into numbers
  bc_error_dat <- lapply(bc_error_dat, FUN = utils::type.convert, as.is = TRUE)

  # extract metrics table and convert to data.frame
  metrics <- bc_error_dat[1:(sep_line - 1)]
  metrics_df <- as.data.frame(do.call(rbind, metrics[-1]))
  colnames(metrics_df) <- metrics[[1]]

  # extract histogram table and convert to data.frame
  hist <- bc_error_dat[(sep_line + 1):length(bc_error_dat)]
  hist_df <- as.data.frame(do.call(rbind, hist[-1]))
  colnames(hist_df) <- hist[[1]]
  
  # return list containing both data.frames
  list(metrics = metrics_df, hist = hist_df)

}

#' Synthesis error histogram
#' 
#' Cell barcode synthesis errors lead to fixed Ts at the end of UMI sequences.
#' This plot shows the number of cells with detected fixed Ts and the base
#' positions where the errors are deteced.
#' 
#' @param synthesis_error_dat List containing synthesis error data produced by
#' read_synthesis_error().
#' 
#' @param sample Sample name for main title.
#' 
#' @param cell_bc_length,mol_bc_length Length (bp) of cell and molecule (UMI)
#' barcodes.
#' 
synthesis_error_hist <- function(synthesis_error_dat, sample,
                                 cell_bc_length = 12, mol_bc_length = 8){
  
  # extract histogram data
  hist <- synthesis_error_dat$hist
  
  # create data.frame with barcode structure and position where synthesis
  # errors are detected (fixed Ts in UMI)
  barcode <- data.frame(base = 1:(cell_bc_length + mol_bc_length),
                        Barcode = factor(c(rep("cell", cell_bc_length),
                                           rep("UMI", mol_bc_length)),
                                         levels = c("cell", "UMI")),
                        fixed_T = 0L)

  barcode[hist[, 1] + cell_bc_length, "fixed_T"] <- hist[, 2]
  
  # define colors for plot
  bc_cols <- c("dodgerblue2", "firebrick1")
  axis_cols <- c(rep(bc_cols[1], cell_bc_length),
                 rep(bc_cols[2], mol_bc_length))
  
  # plot histogram
  ggplot2::ggplot(barcode, ggplot2::aes(x = base, y = fixed_T,
                                       fill = Barcode)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = bc_cols) +
    ggplot2::labs(title = paste0(sample, ": Detected fixed Ts in UMIs"),
                  x = "Base position in barcode",
                  y = "Number of cells with fixed Ts") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = axis_cols)) +
    ggplot2::scale_x_continuous(breaks = c(1:(cell_bc_length + mol_bc_length)))
  
}
