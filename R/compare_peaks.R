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
.compare_peaks_old <- function(chr1, chr2, align = c("apex", "direct", "cross")[1], method = c("closest", "approx")[1]){
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

#' @title Align two chromatograms
#' @description
#' Align two chromatograms on retention time using approx function.
#'
#' @param chr1 A XChromatogram object.
#' @param chr2 A XChromatogram object.
#'
#' @return A list containing the aligned chr1 and chr2.
#'
#' @examples
#' chr1 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP24545")[1]
#' chr2 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP15906")[1]
#' .alignRt(chr1 = chr1, chr2 = chr2)
.alignRt <- function(chr1, chr2){
  tmp1 <- approx(x = chr1@rtime, y = chr1@intensity, xout = chr2@rtime)
  chr1_rtime <- c(chr1@rtime, tmp1$x)
  chr1_intensity <- c(chr1@intensity, tmp1$y)
  chr1_intensity <- chr1_intensity[!duplicated(chr1_rtime)]
  chr1_rtime <- chr1_rtime[!duplicated(chr1_rtime)]
  chr1_intensity <- chr1_intensity[order(chr1_rtime)]
  chr1_rtime <- chr1_rtime[order(chr1_rtime)]

  tmp2 <- approx(x = chr2@rtime, y = chr2@intensity, xout = chr1@rtime)
  chr2_rtime <- c(chr2@rtime, tmp2$x)
  chr2_intensity <- c(chr2@intensity, tmp2$y)
  chr2_intensity <- chr2_intensity[!duplicated(chr2_rtime)]
  chr2_rtime <- chr2_rtime[!duplicated(chr2_rtime)]
  chr2_intensity <- chr2_intensity[order(chr2_rtime)]
  chr2_rtime <- chr2_rtime[order(chr2_rtime)]

  chr1@rtime <- chr1_rtime
  chr1@intensity <- chr1_intensity
  chr2@rtime <- chr2_rtime
  chr2@intensity <- chr2_intensity

  return(list(chr1 = chr1, chr2 = chr2))
}

#' @title Compare two peaks' shape
#' @description
#' Compare two peaks' shape using pearson or cross-correlation function.
#'
#' @param chr1 A XChromatogram object.
#' @param chr2 A XChromatogram object.
#' @param method Approach to calculating peak shape similarity ("pearson" and "cross-correlation").
#' @param apexAlign Whether to align the peaks' apex.
#' @param plot Whether to plot two aligned peaks.
#'
#' @return Peak shape similarity score.
#'
#' @examples
#' chr1 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP24545")[1]
#' chr2 <- xcms::chromPeakChromatograms(data$rawData, peaks = "CP15906")[1]
#' .comparePeaks(chr1 = chr1, chr2 = chr2, method = "cross-correlation", apexAlign = FALSE, plot = TRUE)
.comparePeaks <- function(chr1, chr2, method = c("pearson", "cross-correlation")[1], apexAlign = FALSE, plot = FALSE){
  if(apexAlign){
    apex1_rt <- chr1@rtime[which.max(chr1@intensity)]
    apex2_rt <- chr2@rtime[which.max(chr2@intensity)]
    chr1@rtime <- chr1@rtime - (apex1_rt - apex2_rt)
  }
  tmp <- .alignRt(chr1 = chr1, chr2 = chr2)
  if(plot){
    ref_idx <- which.max(c(max(tmp$chr1@intensity, na.rm = TRUE), max(tmp$chr2@intensity, na.rm = TRUE)))
    other_idx <- c(1, 2)[-ref_idx]
    plot(x = tmp[[ref_idx]]@rtime, y = tmp[[ref_idx]]@intensity, type = "l")
    lines(x = tmp[[other_idx]]@rtime, y = tmp[[other_idx]]@intensity)
  }
  if(method == "pearson"){
    return(cor(x = tmp$chr1@intensity, y = tmp$chr2@intensity, use = "pairwise.complete.obs", method = "pearson"))
  }else if(method == "cross-correlation"){
    return(max(as.numeric(ccf(x = tmp$chr1@intensity, y = tmp$chr2@intensity, plot = FALSE, na.action = na.pass)$acf), na.rm = TRUE))
  }else stop("Method is wrong!")
}
