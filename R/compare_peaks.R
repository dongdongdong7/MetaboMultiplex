#' @title Compare two peaks shape
#' @description
#' Compare two peaks shape using correlation.
#'
#' @param chr1 A XChromatogram object.
#' @param chr2 A XChromatogram object.
#' @param align Whether to align the vertices of two peaks.
#' Apex ("apex") aligns the vertices of two peaks directly.
#' Direct ("direct") dose not change the peak rtime dimension.
#' Cross uses ("cross") cross-correlation to compare the two peaks.
#' @param method Retention time align method.
#'
#' @return A score.
#'
#' @examples
#' chr1 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP24545")[1]
#' chr2 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP15906")[1]
#' compare_peaks(chr1 = chr1, chr2 = chr2)
.compare_peaks <- function(chr1, chr2, align = c("apex", "direct", "cross")[1], method = c("closest", "approx")[1]){
  if(align == "apex"){
    apex1_rt <- chr1@rtime[which.max(chr1@intensity)]
    apex2_rt <- chr2@rtime[which.max(chr2@intensity)]
    chr1@rtime <- chr1@rtime - (apex1_rt - apex2_rt)
    return(MSnbase::compareChromatograms(chr1, chr2, ALIGNFUNARGS = list(method = method)))
  }else if(align == "direct"){
    return(MSnbase::compareChromatograms(chr1, chr2, ALIGNFUNARGS = list(method = method)))
  }else if(align == "cross"){
    chr12 <- MSnbase::alignRt(x = chr1, y = chr2, method = method)
    cc <- ccf(chr12@intensity, chr2@intensity, na.action = na.pass, plot = FALSE)
    return(max(as.numeric(cc$acf), na.rm = TRUE))
  }else stop("align para wrong!")
}
