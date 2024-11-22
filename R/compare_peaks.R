#' @title Compare two peaks shape
#' @description
#' Compare two peaks shape using correlation.
#'
#' @param chr1 A XChromatogram object.
#' @param chr2 A XChromatogram object.
#' @param align Whether to align the vertices of two peaks. apex is true and direct is false.
#' @param method Retention time align method.
#'
#' @return A score.
#'
#' @examples
#' chr1 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP24545")[1]
#' chr2 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP15906")[1]
#' compare_peaks(chr1 = chr1, chr2 = chr2)
.compare_peaks <- function(chr1, chr2, align = c("apex", "direct")[1], method = c("closest", "approx")[1]){
  if(align == "apex"){
    chr1@rtime <- chr1@rtime - (chr1@chromPeaks[1, "rt"] - chr2@chromPeaks[1, "rt"])
  }
  MSnbase::compareChromatograms(chr1, chr2, ALIGNFUNARGS = list(method = method))
}
